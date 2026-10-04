// Règles du jeu : XP, niveaux, séries de jours, badges.
// Ce module est pur (aucune dépendance React Native) pour rester simple à tester.
import type { TrackId } from '../data/types';
import { addDays, dayKey, daysBetween } from './dates';

export const XP = {
  correctAnswer: 10,
  perfectQuiz: 20,
  ficheRead: 10,
  flashcards: 5,
  dailyChallenge: 50,
  mockExam: 30,
  dailyGoalReached: 20,
};

export const DAILY_GOAL_OPTIONS = [30, 50, 100, 150];
export const MAX_FREEZES = 2;

export interface Profile {
  name: string;
  track: TrackId;
  examDate: string;
  dailyGoal: number;
}

export interface ProgressState {
  version: 1;
  profile: Profile | null;
  xp: number;
  streak: { current: number; best: number; lastDay: string | null; freezes: number };
  today: { day: string; xp: number; goalBonusGiven: boolean; challengeDone: boolean };
  /** chapterId → date de première lecture */
  fichesRead: Record<string, string>;
  /** chapterId → meilleur score en % */
  quizBest: Record<string, number>;
  quizCount: number;
  perfectCount: number;
  challengesDone: number;
  /** subjectId → meilleure note /20 à l'examen blanc */
  examBest: Record<string, number>;
  /** questionId → nombre d'erreurs (retiré quand l'élève la réussit en révision) */
  mistakes: Record<string, number>;
  subjectsTouched: string[];
  /** badgeId → date d'obtention */
  badges: Record<string, string>;
  /** jour → XP gagnés (30 derniers jours) */
  history: Record<string, number>;
}

export function initialState(): ProgressState {
  return {
    version: 1,
    profile: null,
    xp: 0,
    streak: { current: 0, best: 0, lastDay: null, freezes: 0 },
    today: { day: dayKey(), xp: 0, goalBonusGiven: false, challengeDone: false },
    fichesRead: {},
    quizBest: {},
    quizCount: 0,
    perfectCount: 0,
    challengesDone: 0,
    examBest: {},
    mistakes: {},
    subjectsTouched: [],
    badges: {},
    history: {},
  };
}

// ---------- Niveaux ----------

const LEVEL_TITLES = [
  'Débutant',
  'Apprenti',
  'Élève sérieux',
  'Bosseur',
  'Passionné',
  'Mention Assez-Bien',
  'Mention Bien',
  'Mention Très Bien',
  'Major de promo',
  'Lauréat du Concours général',
];

/** XP total nécessaire pour atteindre le niveau n (niveau 1 = 0 XP, 2 = 100, 3 = 300, 4 = 600…). */
export function levelThreshold(level: number): number {
  return 50 * (level - 1) * level;
}

export function levelInfo(xp: number) {
  let level = 1;
  while (xp >= levelThreshold(level + 1)) level++;
  const start = levelThreshold(level);
  const end = levelThreshold(level + 1);
  return {
    level,
    title: LEVEL_TITLES[Math.min(level, LEVEL_TITLES.length) - 1],
    progress: (xp - start) / (end - start),
    toNext: end - xp,
  };
}

// ---------- Série de jours ----------

/** Série réellement affichée aujourd'hui (0 si elle est perdue, en tenant compte des gels). */
export function effectiveStreak(state: ProgressState, today = dayKey()): number {
  const { lastDay, current, freezes } = state.streak;
  if (!lastDay) return 0;
  const gap = daysBetween(lastDay, today);
  if (gap <= 1) return current;
  return gap - 1 <= freezes ? current : 0;
}

/** Met à jour la série quand l'élève gagne de l'XP aujourd'hui. */
function touchStreak(s: ProgressState, today: string): string[] {
  const events: string[] = [];
  const { lastDay } = s.streak;
  if (lastDay === today) return events;
  const missed = lastDay ? daysBetween(lastDay, today) - 1 : 0;
  if (lastDay && missed === 0) {
    s.streak.current += 1;
  } else if (lastDay && missed > 0 && missed <= s.streak.freezes) {
    s.streak.freezes -= missed;
    s.streak.current += 1;
    events.push(`🧊 ${missed > 1 ? `${missed} gels utilisés` : 'Gel de série utilisé'} : ta série est sauvée !`);
  } else {
    s.streak.current = 1;
  }
  s.streak.lastDay = today;
  s.streak.best = Math.max(s.streak.best, s.streak.current);
  if (s.streak.current % 7 === 0 && s.streak.freezes < MAX_FREEZES) {
    s.streak.freezes += 1;
    events.push('🧊 7 jours d’affilée : tu gagnes un gel de série !');
  }
  return events;
}

/** Remet à zéro les compteurs du jour si on a changé de jour. */
export function rollDay(s: ProgressState, today = dayKey()): ProgressState {
  if (s.today.day === today) return s;
  return { ...s, today: { day: today, xp: 0, goalBonusGiven: false, challengeDone: false } };
}

// ---------- Badges ----------

export interface Badge {
  id: string;
  icon: string;
  name: string;
  description: string;
  earned: (s: ProgressState) => boolean;
}

export const BADGES: Badge[] = [
  { id: 'first-fiche', icon: '📄', name: 'Premier pas', description: 'Lire ta première fiche', earned: (s) => Object.keys(s.fichesRead).length >= 1 },
  { id: 'fiches-10', icon: '📚', name: 'Rat de bibliothèque', description: 'Lire 10 fiches', earned: (s) => Object.keys(s.fichesRead).length >= 10 },
  { id: 'fiches-30', icon: '🎓', name: 'Encyclopédie', description: 'Lire 30 fiches', earned: (s) => Object.keys(s.fichesRead).length >= 30 },
  { id: 'first-quiz', icon: '✅', name: 'Premier quiz', description: 'Terminer un quiz', earned: (s) => s.quizCount >= 1 },
  { id: 'quiz-25', icon: '🧠', name: 'Machine à quiz', description: 'Terminer 25 quiz', earned: (s) => s.quizCount >= 25 },
  { id: 'perfect', icon: '💯', name: 'Sans faute', description: 'Réussir un quiz à 100 %', earned: (s) => s.perfectCount >= 1 },
  { id: 'perfect-10', icon: '🏹', name: 'Tireur d’élite', description: '10 quiz sans faute', earned: (s) => s.perfectCount >= 10 },
  { id: 'streak-3', icon: '🔥', name: 'Ça chauffe', description: 'Série de 3 jours', earned: (s) => s.streak.best >= 3 },
  { id: 'streak-7', icon: '⚡', name: 'Semaine parfaite', description: 'Série de 7 jours', earned: (s) => s.streak.best >= 7 },
  { id: 'streak-30', icon: '🦁', name: 'Lion de la Teranga', description: 'Série de 30 jours', earned: (s) => s.streak.best >= 30 },
  { id: 'challenge-1', icon: '🎯', name: 'Défi relevé', description: 'Terminer un défi du jour', earned: (s) => s.challengesDone >= 1 },
  { id: 'challenge-10', icon: '🏆', name: 'Champion des défis', description: 'Terminer 10 défis du jour', earned: (s) => s.challengesDone >= 10 },
  { id: 'exam-pass', icon: '📝', name: 'Admis !', description: 'Avoir au moins 10/20 à un examen blanc', earned: (s) => Object.values(s.examBest).some((n) => n >= 10) },
  { id: 'exam-tb', icon: '🌟', name: 'Mention Très Bien', description: 'Avoir au moins 16/20 à un examen blanc', earned: (s) => Object.values(s.examBest).some((n) => n >= 16) },
  { id: 'polyvalent', icon: '🧭', name: 'Polyvalent', description: 'Réviser 5 matières différentes', earned: (s) => s.subjectsTouched.length >= 5 },
  { id: 'xp-1000', icon: '💎', name: 'Millionnaire… en XP', description: 'Cumuler 1 000 XP', earned: (s) => s.xp >= 1000 },
];

// ---------- Gains d'XP ----------

export interface Reward {
  xp: number;
  messages: string[];
  newBadges: Badge[];
  levelUp: number | null;
}

/**
 * Applique un gain d'XP « brut » et tous ses effets de bord :
 * série, objectif du jour, historique, badges, montée de niveau.
 * `mutate` modifie la copie de l'état (compteurs propres à l'activité) avant le calcul des badges.
 */
export function applyGain(
  prev: ProgressState,
  baseXp: number,
  baseMessages: string[],
  mutate: (s: ProgressState) => void,
  today = dayKey(),
): { state: ProgressState; reward: Reward } {
  const s: ProgressState = structuredCloneState(rollDay(prev, today));
  mutate(s);
  const messages = [...baseMessages];
  let xp = baseXp;

  if (xp > 0) {
    messages.push(...touchStreak(s, today));
    const goal = s.profile?.dailyGoal ?? 50;
    if (!s.today.goalBonusGiven && s.today.xp + xp >= goal) {
      s.today.goalBonusGiven = true;
      xp += XP.dailyGoalReached;
      messages.push(`🎉 Objectif du jour atteint : +${XP.dailyGoalReached} XP bonus`);
    }
  }

  const levelBefore = levelInfo(s.xp).level;
  s.xp += xp;
  s.today.xp += xp;
  s.history[today] = (s.history[today] ?? 0) + xp;
  pruneHistory(s, today);
  const levelAfter = levelInfo(s.xp).level;

  const newBadges = BADGES.filter((b) => !s.badges[b.id] && b.earned(s));
  for (const b of newBadges) s.badges[b.id] = today;

  return { state: s, reward: { xp, messages, newBadges, levelUp: levelAfter > levelBefore ? levelAfter : null } };
}

function pruneHistory(s: ProgressState, today: string) {
  const limit = addDays(today, -30);
  for (const k of Object.keys(s.history)) if (k < limit) delete s.history[k];
}

function structuredCloneState(s: ProgressState): ProgressState {
  return JSON.parse(JSON.stringify(s));
}

export type QuizMode = 'chapter' | 'daily' | 'exam' | 'review';

export interface AnswerResult {
  questionId: string;
  subjectId: string;
  correct: boolean;
}

export function quizGain(
  prev: ProgressState,
  mode: QuizMode,
  results: AnswerResult[],
  opts: { chapterId?: string; subjectId?: string } = {},
  today = dayKey(),
) {
  const correct = results.filter((r) => r.correct).length;
  const total = results.length;
  const percent = total ? Math.round((correct / total) * 100) : 0;
  const messages: string[] = [];
  let xp = correct * XP.correctAnswer;
  messages.push(`${correct} bonne${correct > 1 ? 's' : ''} réponse${correct > 1 ? 's' : ''} : +${xp} XP`);

  const perfect = total >= 5 && correct === total;
  if (perfect) {
    xp += XP.perfectQuiz;
    messages.push(`💯 Sans faute : +${XP.perfectQuiz} XP`);
  }
  const challengeAlreadyDone = rollDay(prev, today).today.challengeDone;
  if (mode === 'daily' && !challengeAlreadyDone) {
    xp += XP.dailyChallenge;
    messages.push(`🎯 Défi du jour terminé : +${XP.dailyChallenge} XP`);
  }
  if (mode === 'exam') {
    xp += XP.mockExam;
    messages.push(`📝 Examen blanc terminé : +${XP.mockExam} XP`);
  }

  return applyGain(
    prev,
    xp,
    messages,
    (s) => {
      s.quizCount += 1;
      if (perfect) s.perfectCount += 1;
      if (mode === 'daily' && !challengeAlreadyDone) {
        s.today.challengeDone = true;
        s.challengesDone += 1;
      }
      if (mode === 'chapter' && opts.chapterId) {
        s.quizBest[opts.chapterId] = Math.max(s.quizBest[opts.chapterId] ?? 0, percent);
      }
      if (mode === 'exam' && opts.subjectId) {
        const note = Math.round((correct / Math.max(total, 1)) * 20 * 2) / 2;
        s.examBest[opts.subjectId] = Math.max(s.examBest[opts.subjectId] ?? 0, note);
      }
      for (const r of results) {
        if (r.correct) {
          if (mode === 'review') delete s.mistakes[r.questionId];
        } else {
          s.mistakes[r.questionId] = (s.mistakes[r.questionId] ?? 0) + 1;
        }
        if (!s.subjectsTouched.includes(r.subjectId)) s.subjectsTouched.push(r.subjectId);
      }
    },
    today,
  );
}

/** Mention correspondant à une note sur 20 (barème des examens sénégalais). */
export function mention(note: number): string {
  if (note >= 16) return 'Très Bien';
  if (note >= 14) return 'Bien';
  if (note >= 12) return 'Assez-Bien';
  if (note >= 10) return 'Passable';
  return 'Insuffisant';
}

/** Compteurs du jour, remis à zéro si l'application est restée ouverte après minuit. */
export function todayStats(s: ProgressState, today = dayKey()) {
  return rollDay(s, today).today;
}
