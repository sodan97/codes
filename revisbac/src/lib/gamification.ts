// Règles du jeu : XP, niveaux, séries de jours, badges, quêtes et coffre du jour, suivi des flashcards.
// Ce module est pur (aucune dépendance React Native) pour rester simple à tester.
import { getSubjects } from '../data/catalog';
import { getTrack, isTrackId } from '../data/tracks';
import type { TrackId } from '../data/types';
import { addDays, dayKey, daysBetween, isDayKey } from './dates';
import {
  advanceQuests,
  allQuestsDone,
  generateQuests,
  QUEST_KINDS,
  QUEST_TIERS,
  questTitle,
  type ActivityEvents,
  type Quest,
} from './quests';
import { seededRandom } from './random';
import { capDue, scheduleCard, type CardBox, type CardEntry } from './srs';
import { chapterStars, STAR_QUIZ_PERCENT } from './stats';

export const XP = {
  correctAnswer: 10,
  perfectQuiz: 20,
  ficheRead: 10,
  flashcards: 5,
  dailyChallenge: 50,
  mockExam: 30,
  dailyGoalReached: 20,
  /** Par quête du jour accomplie. */
  quest: 10,
};

/** Coffre du jour (les 3 quêtes faites) : entre 20 et 40 XP, ou un gel de série. */
export const CHEST_XP_MIN = 20;
export const CHEST_XP_MAX = 40;
/** Chance que le coffre contienne un gel, quand l'élève en a moins de MAX_FREEZES. */
const CHEST_FREEZE_CHANCE = 0.3;

export const DAILY_GOAL_OPTIONS = [30, 50, 100, 150];
export const MAX_FREEZES = 2;
/** Version du schéma de sauvegarde (voir migrate). */
export const STATE_VERSION = 3;
/** Nombre de notes d'examen blanc gardées par matière. */
export const EXAM_HISTORY_SIZE = 5;

export interface ReminderSettings {
  enabled: boolean;
  hour: number;
  minute: number;
}

export interface Profile {
  name: string;
  track: TrackId;
  examDate: string;
  dailyGoal: number;
  /** Matières que l'élève ne passe pas (ex. LV2) : masquées partout, progression conservée. */
  hiddenSubjects: string[];
  /** Rappel quotidien local, null tant que l'élève n'a rien choisi. */
  reminder: ReminderSettings | null;
  /** Vibrations aux réponses et célébrations. */
  haptics: boolean;
}

/** Erreur suivie en répétition espacée (Leitner 3 boîtes : revue à J+1, puis J+3, puis J+7). */
export interface MistakeEntry {
  box: 1 | 2 | 3;
  /** Jour (AAAA-MM-JJ) à partir duquel la question est à revoir. */
  due: string;
  /** Nombre total d'erreurs sur cette question. */
  fails: number;
  /** Dernier jour où l'entrée a été créée ou a changé de boîte. */
  last?: string;
}

/** Score officiel du défi du jour (première tentative uniquement). */
export interface DailyResult {
  correct: number;
  total: number;
  /** Suite de ✅ et ❌ dans l'ordre des questions (identique pour tous les candidats). */
  grid: string;
}

export interface ExamNote {
  day: string;
  note: number;
}

export interface ProgressState {
  version: 3;
  profile: Profile | null;
  xp: number;
  streak: { current: number; best: number; lastDay: string | null; freezes: number };
  /**
   * Compteurs du jour. `done` liste les activités déjà terminées : 'daily', `quiz:{chapterId}`, `flash:{chapterId}`, `exam:{subjectId}`.
   * `quests` : quêtes du jour (vides sans profil), `chestOpened` : coffre du jour déjà ouvert.
   */
  today: { day: string; xp: number; goalBonusGiven: boolean; challengeDone: boolean; done: string[]; quests: Quest[]; chestOpened: boolean };
  /** chapterId → date de première lecture */
  fichesRead: Record<string, string>;
  /** chapterId → meilleur score en % */
  quizBest: Record<string, number>;
  quizCount: number;
  perfectCount: number;
  challengesDone: number;
  /** subjectId → meilleure note /20 à l'examen blanc */
  examBest: Record<string, number>;
  /** subjectId → dernières notes (les plus anciennes d'abord, 5 au plus) */
  examHistory: Record<string, ExamNote[]>;
  /** subjectId → ids des questions du dernier examen blanc */
  examLast: Record<string, string[]>;
  /** questionId → suivi de l'erreur (sort de la liste après 3 réussites espacées) */
  mistakes: Record<string, MistakeEntry>;
  /** jour → score officiel du défi du jour (30 derniers jours) */
  dailyResults: Record<string, DailyResult>;
  subjectsTouched: string[];
  /** badgeId → date d'obtention */
  badges: Record<string, string>;
  /** jour → XP gagnés (30 derniers jours) */
  history: Record<string, number>;
  /** clé de carte (voir srs.cardKey) → suivi de la flashcard (Leitner 5 boîtes) */
  cards: Record<string, CardEntry>;
  /** chapterId → dernier jour où le paquet de flashcards a été terminé */
  flashcardsDone: Record<string, string>;
  /** chapterId → jours (5 derniers) où le quiz de chapitre a été réussi à 80 % ou plus (étoiles, voir stats.ts) */
  quiz80Days: Record<string, string[]>;
}

/** Nombre de jours gardés dans quiz80Days pour un chapitre. */
const QUIZ80_DAYS_KEPT = 5;

export function initialState(today = dayKey()): ProgressState {
  return {
    version: STATE_VERSION,
    profile: null,
    xp: 0,
    streak: { current: 0, best: 0, lastDay: null, freezes: 0 },
    today: emptyDay(today),
    fichesRead: {},
    quizBest: {},
    quizCount: 0,
    perfectCount: 0,
    challengesDone: 0,
    examBest: {},
    examHistory: {},
    examLast: {},
    mistakes: {},
    dailyResults: {},
    subjectsTouched: [],
    badges: {},
    history: {},
    cards: {},
    flashcardsDone: {},
    quiz80Days: {},
  };
}

function emptyDay(day: string): ProgressState['today'] {
  return { day, xp: 0, goalBonusGiven: false, challengeDone: false, done: [], quests: [], chestOpened: false };
}

/** Profil complet avec les réglages par défaut (matières toutes visibles, pas de rappel, vibrations actives). */
export function createProfile(fields: Pick<Profile, 'name' | 'track' | 'examDate' | 'dailyGoal'> & Partial<Profile>): Profile {
  return { hiddenSubjects: [], reminder: null, haptics: true, ...fields };
}

// ---------- Lecture d'une sauvegarde ----------

type Json = Record<string, unknown>;

const isObject = (v: unknown): v is Json => typeof v === 'object' && v !== null && !Array.isArray(v);
const isFiniteNumber = (v: unknown): v is number => typeof v === 'number' && Number.isFinite(v);
const clamp = (n: number, min: number, max: number) => Math.min(max, Math.max(min, n));
/** Entier ≥ 0 (ou la valeur par défaut si ce n'est pas un nombre). */
const count = (v: unknown, fallback = 0) => (isFiniteNumber(v) ? Math.max(0, Math.floor(v)) : fallback);
const bool = (v: unknown, fallback: boolean) => (typeof v === 'boolean' ? v : fallback);
const stringList = (v: unknown): string[] => (Array.isArray(v) ? [...new Set(v.filter((x): x is string => typeof x === 'string'))] : []);

/** Recopie un dictionnaire en ne gardant que les entrées valides (transformées par `read`, undefined = ignorée). */
function record<T>(v: unknown, read: (value: unknown, key: string) => T | undefined): Record<string, T> {
  const out: Record<string, T> = {};
  if (!isObject(v)) return out;
  for (const [k, value] of Object.entries(v)) {
    const r = read(value, k);
    if (r !== undefined) out[k] = r;
  }
  return out;
}

function readMistake(v: unknown, today: string): MistakeEntry | undefined {
  if (!isObject(v)) return undefined;
  const box = v.box === 2 || v.box === 3 ? v.box : 1;
  return { box, due: isDayKey(v.due) ? v.due : today, fails: count(v.fails, 1), ...(isDayKey(v.last) ? { last: v.last } : {}) };
}

function readCard(v: unknown, today: string): CardEntry | undefined {
  if (!isObject(v)) return undefined;
  const box = (isFiniteNumber(v.box) ? clamp(Math.floor(v.box), 1, 5) : 1) as CardBox;
  return { box, due: isDayKey(v.due) ? v.due : today, last: isDayKey(v.last) ? v.last : today };
}

function readQuest(v: unknown): Quest | undefined {
  if (!isObject(v) || typeof v.id !== 'string') return undefined;
  const kind = QUEST_KINDS.find((k) => k === v.kind);
  const tier = QUEST_TIERS.find((t) => t === v.tier);
  if (!kind || !tier) return undefined;
  const target = Math.max(1, count(v.target, 1));
  const progress = Math.min(count(v.progress), target);
  return {
    id: v.id,
    kind,
    tier,
    target,
    progress,
    done: bool(v.done, false) || progress >= target,
    ...(typeof v.subjectId === 'string' ? { subjectId: v.subjectId } : {}),
  };
}

function readProfile(v: unknown): Profile | null {
  if (!isObject(v) || !isTrackId(v.track)) return null;
  const track = v.track;
  const name = typeof v.name === 'string' && v.name.trim() ? v.name.trim().slice(0, 30) : 'Champion';
  const dailyGoal = DAILY_GOAL_OPTIONS.includes(v.dailyGoal as number) ? (v.dailyGoal as number) : 50;
  const examDate = isDayKey(v.examDate) ? v.examDate : getTrack(track).defaultExamDate;
  const r = v.reminder;
  const reminder =
    isObject(r) && isFiniteNumber(r.hour) && isFiniteNumber(r.minute)
      ? { enabled: bool(r.enabled, false), hour: clamp(Math.floor(r.hour), 0, 23), minute: clamp(Math.floor(r.minute), 0, 59) }
      : null;
  return { name, track, examDate, dailyGoal, hiddenSubjects: stringList(v.hiddenSubjects), reminder, haptics: bool(v.haptics, true) };
}

/** v1 → v2 : `mistakes[id]` passe d'un nombre d'erreurs à une entrée Leitner en boîte 1, à revoir aujourd'hui. */
function migrateV1(data: Json, today: string): Json {
  const mistakes = record(data.mistakes, (n) => (isFiniteNumber(n) ? { box: 1, due: today, fails: Math.max(1, Math.floor(n)) } : n));
  return { ...data, mistakes, version: 2 };
}

/**
 * Transforme une sauvegarde (quelle que soit sa version) en état valide et à jour.
 * Chaque champ est recopié et vérifié en profondeur ; les valeurs aberrantes sont bornées.
 * Les identifiants inconnus du catalogue sont gardés (filtrés à l'affichage, voir selectors.ts).
 * Lève une erreur si `raw` n'est pas un objet.
 */
export function migrate(raw: unknown, today = dayKey()): ProgressState {
  if (!isObject(raw)) throw new Error('Sauvegarde illisible : un objet était attendu.');
  let data = raw;
  const version = isFiniteNumber(raw.version) ? raw.version : 1;
  switch (version) {
    case 1:
      data = migrateV1(data, today);
    // falls through
    case 2:
    // v2 → v3 : flashcards suivies, jours de quiz réussis, quêtes et coffre du jour (valeurs par défaut).
    // falls through
    case 3:
      break;
    default:
    // Version plus récente que l'application : on lit tout ce qui est compréhensible.
  }

  const s = initialState(today);
  s.profile = readProfile(data.profile);
  s.xp = count(data.xp);

  if (isObject(data.streak)) {
    const st = data.streak;
    s.streak.current = count(st.current);
    s.streak.best = Math.max(count(st.best), s.streak.current);
    s.streak.lastDay = isDayKey(st.lastDay) ? st.lastDay : null;
    s.streak.freezes = clamp(count(st.freezes), 0, MAX_FREEZES);
  }
  if (isObject(data.today) && isDayKey(data.today.day)) {
    const t = data.today;
    s.today = {
      day: t.day as string,
      xp: count(t.xp),
      goalBonusGiven: bool(t.goalBonusGiven, false),
      challengeDone: bool(t.challengeDone, false),
      done: stringList(t.done),
      quests: Array.isArray(t.quests) ? t.quests.map(readQuest).filter((q): q is Quest => !!q) : [],
      chestOpened: bool(t.chestOpened, false),
    };
  }

  s.fichesRead = record(data.fichesRead, (v) => (isDayKey(v) ? v : undefined));
  s.quizBest = record(data.quizBest, (v) => (isFiniteNumber(v) ? clamp(Math.round(v), 0, 100) : undefined));
  s.quizCount = count(data.quizCount);
  s.perfectCount = count(data.perfectCount);
  s.challengesDone = count(data.challengesDone);
  s.examBest = record(data.examBest, (v) => (isFiniteNumber(v) ? clamp(v, 0, 20) : undefined));
  s.examHistory = record(data.examHistory, (v) => {
    if (!Array.isArray(v)) return undefined;
    const notes = v
      .filter((e): e is Json => isObject(e) && isDayKey(e.day) && isFiniteNumber(e.note))
      .map((e) => ({ day: e.day as string, note: clamp(e.note as number, 0, 20) }));
    return notes.length ? notes.slice(-EXAM_HISTORY_SIZE) : undefined;
  });
  s.examLast = record(data.examLast, (v) => (Array.isArray(v) ? stringList(v) : undefined));
  s.mistakes = record(data.mistakes, (v) => readMistake(v, today));
  s.dailyResults = record(data.dailyResults, (v, day) => {
    if (!isDayKey(day) || !isObject(v)) return undefined;
    const total = count(v.total);
    return { correct: Math.min(count(v.correct), total), total, grid: typeof v.grid === 'string' ? v.grid : '' };
  });
  s.subjectsTouched = stringList(data.subjectsTouched);
  s.badges = record(data.badges, (v) => (isDayKey(v) ? v : undefined));
  s.history = record(data.history, (v, day) => (isDayKey(day) && isFiniteNumber(v) ? Math.max(0, v) : undefined));
  s.cards = record(data.cards, (v) => readCard(v, today));
  s.flashcardsDone = record(data.flashcardsDone, (v) => (isDayKey(v) ? v : undefined));
  // v2 : un paquet terminé le jour même n'était noté que dans today.done (avant que rollDay ne le vide).
  if (version < 3 && isObject(data.today) && isDayKey(data.today.day)) {
    for (const key of s.today.done) {
      if (key.startsWith('flash:')) s.flashcardsDone[key.slice(6)] ??= s.today.day;
    }
  }
  s.quiz80Days = record(data.quiz80Days, (v) => {
    if (!Array.isArray(v)) return undefined;
    const days = [...new Set(v.filter(isDayKey))].sort().slice(-QUIZ80_DAYS_KEPT);
    return days.length ? days : undefined;
  });

  return ensureQuests(rollDay(s, today), today);
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

/** Met à jour la série quand l'élève termine une activité aujourd'hui. */
function touchStreak(s: ProgressState, today: string): string[] {
  const events: string[] = [];
  const { lastDay } = s.streak;
  // Horloge qui recule (fuseau, date mal réglée) : on ne casse pas la série.
  if (lastDay !== null && today <= lastDay) return events;
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

/**
 * Remet à zéro les compteurs du jour si on a changé de jour, et tire les quêtes du nouveau jour
 * (celles de la veille disparaissent, faites ou non, sans message d'échec).
 */
export function rollDay(s: ProgressState, today = dayKey()): ProgressState {
  // On n'avance que vers un jour plus récent : une horloge qui recule ne rouvre pas la journée.
  if (today <= s.today.day) return s;
  const next = { ...s, today: emptyDay(today) };
  next.today.quests = generateQuests(next, today);
  return next;
}

/**
 * Tire les quêtes du jour si l'élève a un profil et qu'il n'en a pas encore (sauvegarde v2, profil créé en cours de journée).
 * Renvoie le même état s'il n'y a rien à faire.
 */
export function ensureQuests(s: ProgressState, today = dayKey()): ProgressState {
  if (!s.profile || s.today.quests.length > 0 || s.today.day !== today) return s;
  const quests = generateQuests(s, today);
  return quests.length ? { ...s, today: { ...s.today, quests } } : s;
}

/**
 * Nouveau profil (onboarding, réglages, changement d'examen).
 * • Examen changé, quêtes du jour pas encore entamées : elles sont tirées à nouveau (elles visent les matières de l'examen).
 * • Sinon, les quêtes du jour sont gardées (pas de double gain, un coffre prêt reste prêt) ; seule la quête de variété
 *   non faite qui vise une matière absente ou masquée est tirée à nouveau parmi les matières visibles.
 */
export function withProfile(prev: ProgressState, profile: Profile, today = dayKey()): ProgressState {
  const rolled = rollDay(prev, today);
  const trackChanged = rolled.profile?.track !== profile.track;
  const s: ProgressState = { ...rolled, profile };
  if (s.today.day !== today) return s;
  const started = s.today.chestOpened || s.today.quests.some((q) => q.done || q.progress > 0);
  if (trackChanged && !started) {
    const fresh = generateQuests(s, today);
    return fresh.length ? { ...s, today: { ...s.today, quests: fresh } } : s;
  }
  const visible = new Set(getSubjects(profile.track).filter((sub) => !profile.hiddenSubjects.includes(sub.id)).map((sub) => sub.id));
  const stale = s.today.quests.find((q) => q.kind === 'subject' && !q.done && (!q.subjectId || !visible.has(q.subjectId)));
  if (!stale) return ensureQuests(s, today);
  // Le tirage est fixé par jour et par examen : la quête de variété est recalculée parmi les matières visibles.
  const fresh = generateQuests(s, today).find((q) => q.tier === 'variete');
  const keep = fresh && !s.today.quests.some((o) => o !== stale && o.id === fresh.id);
  const quests = s.today.quests.flatMap((q) => (q === stale ? (keep ? [fresh] : []) : [q]));
  return { ...s, today: { ...s.today, quests } };
}

// ---------- Badges ----------

export interface BadgeProgress {
  value: number;
  target: number;
  /** Texte prêt à afficher sous un badge verrouillé : « 7/10 fiches ». */
  label: string;
}

export interface Badge {
  id: string;
  icon: string;
  name: string;
  description: string;
  earned: (s: ProgressState) => boolean;
  /** Avancement vers le badge (badges à compteur), pour l'afficher tant qu'il est verrouillé. */
  progress?: (s: ProgressState) => BadgeProgress;
}

const fmt = (n: number) => String(n).replace(/\B(?=(\d{3})+(?!\d))/g, '\u00a0');
/** Avancement « valeur/cible unité », la valeur étant plafonnée à la cible. */
const counter = (value: number, target: number, unit: string): BadgeProgress => {
  const v = Math.min(value, target);
  return { value: v, target, label: `${fmt(v)}/${fmt(target)} ${unit}` };
};
const fichesCount = (s: ProgressState) => Object.keys(s.fichesRead).length;

/** Meilleur nombre d'étoiles parmi les chapitres travaillés. */
function bestChapterStars(s: ProgressState): number {
  let best = 0;
  for (const id of Object.keys(s.fichesRead)) best = Math.max(best, chapterStars(id, s, s.today.day).stars);
  return best;
}

/** Matière de l'examen de l'élève la plus proche d'être entièrement à 3 étoiles. */
function bestSubjectGold(s: ProgressState): { gold: number; total: number } {
  let best = { gold: 0, total: 0 };
  if (!s.profile) return best;
  for (const subject of getSubjects(s.profile.track)) {
    const total = subject.chapters.length;
    if (!total) continue;
    const gold = subject.chapters.filter((c) => chapterStars(c.id, s, s.today.day).stars === 3).length;
    if (best.total === 0 || gold / total > best.gold / best.total) best = { gold, total };
  }
  return best;
}

export const BADGES: Badge[] = [
  { id: 'first-fiche', icon: '📄', name: 'Premier pas', description: 'Lire ta première fiche', earned: (s) => fichesCount(s) >= 1 },
  { id: 'fiches-10', icon: '📚', name: 'Rat de bibliothèque', description: 'Lire 10 fiches', earned: (s) => fichesCount(s) >= 10, progress: (s) => counter(fichesCount(s), 10, 'fiches') },
  { id: 'fiches-30', icon: '🎓', name: 'Encyclopédie', description: 'Lire 30 fiches', earned: (s) => fichesCount(s) >= 30, progress: (s) => counter(fichesCount(s), 30, 'fiches') },
  { id: 'first-quiz', icon: '✅', name: 'Premier quiz', description: 'Terminer un quiz', earned: (s) => s.quizCount >= 1 },
  { id: 'quiz-25', icon: '🧠', name: 'Machine à quiz', description: 'Terminer 25 quiz', earned: (s) => s.quizCount >= 25, progress: (s) => counter(s.quizCount, 25, 'quiz') },
  { id: 'perfect', icon: '💯', name: 'Sans faute', description: 'Réussir un quiz à 100 %', earned: (s) => s.perfectCount >= 1 },
  { id: 'perfect-10', icon: '🏹', name: 'Tireur d’élite', description: '10 quiz sans faute', earned: (s) => s.perfectCount >= 10, progress: (s) => counter(s.perfectCount, 10, 'sans faute') },
  { id: 'streak-3', icon: '🔥', name: 'Ça chauffe', description: 'Série de 3 jours', earned: (s) => s.streak.best >= 3, progress: (s) => counter(s.streak.best, 3, 'jours') },
  { id: 'streak-7', icon: '⚡', name: 'Semaine parfaite', description: 'Série de 7 jours', earned: (s) => s.streak.best >= 7, progress: (s) => counter(s.streak.best, 7, 'jours') },
  { id: 'streak-30', icon: '🦁', name: 'Lion de la Teranga', description: 'Série de 30 jours', earned: (s) => s.streak.best >= 30, progress: (s) => counter(s.streak.best, 30, 'jours') },
  { id: 'challenge-1', icon: '🎯', name: 'Défi relevé', description: 'Terminer un défi du jour', earned: (s) => s.challengesDone >= 1 },
  { id: 'challenge-10', icon: '🏆', name: 'Champion des défis', description: 'Terminer 10 défis du jour', earned: (s) => s.challengesDone >= 10, progress: (s) => counter(s.challengesDone, 10, 'défis') },
  { id: 'exam-pass', icon: '📝', name: 'Admis !', description: 'Avoir au moins 10/20 à un examen blanc', earned: (s) => Object.values(s.examBest).some((n) => n >= 10) },
  { id: 'exam-tb', icon: '🌟', name: 'Mention Très Bien', description: 'Avoir au moins 16/20 à un examen blanc', earned: (s) => Object.values(s.examBest).some((n) => n >= 16) },
  { id: 'polyvalent', icon: '🧭', name: 'Polyvalent', description: 'Réviser 5 matières différentes', earned: (s) => s.subjectsTouched.length >= 5, progress: (s) => counter(s.subjectsTouched.length, 5, 'matières') },
  { id: 'xp-1000', icon: '💎', name: 'Millionnaire… en XP', description: 'Cumuler 1 000 XP', earned: (s) => s.xp >= 1000, progress: (s) => counter(s.xp, 1000, 'XP') },
  {
    id: 'chapter-gold',
    icon: '🥇',
    name: 'Chapitre en or',
    description: 'Obtenir 3 étoiles sur un chapitre',
    earned: (s) => bestChapterStars(s) >= 3,
    progress: (s) => counter(bestChapterStars(s), 3, 'étoiles'),
  },
  {
    id: 'subject-master',
    icon: '👑',
    name: 'Matière maîtrisée',
    description: 'Obtenir 3 étoiles sur tous les chapitres d’une matière',
    earned: (s) => {
      const { gold, total } = bestSubjectGold(s);
      return total > 0 && gold === total;
    },
    progress: (s) => {
      const { gold, total } = bestSubjectGold(s);
      return counter(gold, Math.max(total, 1), 'chapitres en or');
    },
  },
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
 * `mutate` modifie la copie de l'état (compteurs propres à l'activité) avant le calcul des badges ;
 * il peut renvoyer des messages à ajouter à la récompense.
 * `opts.activity` : une activité a été terminée, la série est prolongée même sans XP
 * (par défaut : seulement s'il y a des XP). Un jour présent dans `history` est un jour d'activité, même à 0 XP.
 * `opts.events` : ce qui s'est passé pendant l'activité, pour faire avancer les quêtes du jour
 * (+10 XP par quête accomplie). `mutate` peut le compléter (ex. erreurs corrigées).
 */
export function applyGain(
  prev: ProgressState,
  baseXp: number,
  baseMessages: string[],
  mutate: (s: ProgressState) => string[] | void,
  today = dayKey(),
  opts: { activity?: boolean; events?: ActivityEvents } = {},
): { state: ProgressState; reward: Reward } {
  const s: ProgressState = structuredCloneState(ensureQuests(rollDay(prev, today), today));
  const messages = [...baseMessages, ...(mutate(s) ?? [])];
  let xp = baseXp;
  const active = opts.activity ?? xp > 0;

  if (opts.events) {
    const hadAll = allQuestsDone(s.today.quests);
    for (const q of advanceQuests(s.today.quests, opts.events)) {
      xp += XP.quest;
      messages.push(`🗺️ Quête accomplie : ${questTitle(q)} : +${XP.quest} XP`);
    }
    if (!hadAll && allQuestsDone(s.today.quests) && !s.today.chestOpened) messages.push('🎁 Les 3 quêtes du jour sont faites : ouvre ton coffre !');
  }

  if (active) messages.push(...touchStreak(s, today));
  if (xp > 0) {
    // L'objectif du jour reste calculé sur l'XP.
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
  if (active) s.history[today] = (s.history[today] ?? 0) + xp;
  prune(s.history, today);
  prune(s.dailyResults, today);
  const levelAfter = levelInfo(s.xp).level;

  const newBadges = BADGES.filter((b) => !s.badges[b.id] && b.earned(s));
  for (const b of newBadges) s.badges[b.id] = today;

  return { state: s, reward: { xp, messages, newBadges, levelUp: levelAfter > levelBefore ? levelAfter : null } };
}

/** Garde seulement les 30 derniers jours d'un dictionnaire indexé par jour. */
function prune(byDay: Record<string, unknown>, today: string) {
  const limit = addDays(today, -30);
  for (const k of Object.keys(byDay)) if (k < limit) delete byDay[k];
}

function structuredCloneState(s: ProgressState): ProgressState {
  return JSON.parse(JSON.stringify(s));
}

function markDone(s: ProgressState, key: string) {
  if (!s.today.done.includes(key)) s.today.done.push(key);
}

function touchSubject(s: ProgressState, subjectId: string) {
  if (!s.subjectsTouched.includes(subjectId)) s.subjectsTouched.push(subjectId);
}

/** Lecture d'une fiche : +10 XP la première fois seulement (null si déjà lue). */
export function ficheGain(prev: ProgressState, chapterId: string, subjectId: string, today = dayKey()) {
  if (prev.fichesRead[chapterId]) return null;
  return applyGain(
    prev,
    XP.ficheRead,
    [`📄 Fiche lue : +${XP.ficheRead} XP`],
    (s) => {
      s.fichesRead[chapterId] = today;
      touchSubject(s, subjectId);
    },
    today,
    { events: { kind: 'fiche', chapterId, subjectId } },
  );
}

/**
 * Paquet de flashcards terminé : +5 XP une fois par jour et par chapitre, la série compte toujours.
 * Les cartes elles-mêmes ne rapportent pas d'XP (voir cardReview).
 */
export function flashcardsGain(prev: ProgressState, subjectId: string, chapterId: string, today = dayKey()) {
  const key = `flash:${chapterId}`;
  const already = rollDay(prev, today).today.done.includes(key);
  const xp = already ? 0 : XP.flashcards;
  const message = already ? 'Paquet déjà terminé aujourd’hui : pas d’XP, mais bravo pour la révision !' : `🃏 Flashcards terminées : +${XP.flashcards} XP`;
  return applyGain(
    prev,
    xp,
    [message],
    (s) => {
      markDone(s, key);
      touchSubject(s, subjectId);
      s.flashcardsDone[chapterId] = today;
    },
    today,
    { activity: true, events: { kind: 'flashcards', chapterId, subjectId } },
  );
}

/** Réponse à une flashcard (« Je savais » / « À revoir ») : met à jour son suivi, sans XP. */
export function cardReview(prev: ProgressState, key: string, knew: boolean, today = dayKey()): ProgressState {
  const s = rollDay(prev, today);
  const entry = scheduleCard(s.cards[key], knew, today, s.profile?.examDate);
  if (entry === s.cards[key]) return s;
  return { ...s, cards: { ...s.cards, [key]: entry } };
}

export type QuizMode = 'chapter' | 'daily' | 'exam' | 'review' | 'express';

export interface AnswerResult {
  questionId: string;
  subjectId: string;
  correct: boolean;
  /** Question non traitée (fin du temps en examen) : compte 0 mais n'est jamais ajoutée aux erreurs. */
  skipped?: boolean;
}

/** Note sur 20 arrondie au demi-point. Les questions non traitées comptent 0. */
export function examNote(correct: number, total: number): number {
  return Math.round((correct / Math.max(total, 1)) * 20 * 2) / 2;
}

/** Note au format français, sans « /20 » : 13.5 → « 13,5 », 12 → « 12 ». */
export function formatNote(n: number): string {
  return String(Math.round(n * 2) / 2).replace('.', ',');
}

/** Évolution des dernières notes, « 8 → 11,5 → 13 » (null s'il y en a moins de 2). */
export function noteTrend(history: ExamNote[] | undefined): string | null {
  if (!history || history.length < 2) return null;
  return history.map((e) => formatNote(e.note)).join(' → ');
}

/** Intervalle (en jours) avant la prochaine révision, selon la boîte atteinte. */
const BOX_DAYS: Record<1 | 2 | 3, number> = { 1: 1, 2: 3, 3: 7 };

/** Échéance plafonnée pour que tout soit revu avant l'épreuve : au plus tard examDate - 2, jamais avant demain. */
function cappedDue(s: ProgressState, due: string, today: string): string {
  return capDue(due, today, s.profile?.examDate);
}

export type MistakeOutcome = 'none' | 'added' | 'promoted' | 'mastered';

/**
 * Répétition espacée des erreurs (Leitner 3 boîtes), appliquée à chaque réponse quel que soit le mode.
 * Modifie `s.mistakes` (s est la copie de travail de quizGain) et renvoie ce qui s'est passé.
 * • non traitée (skipped) : rien ;
 * • mauvaise réponse : boîte 1, à revoir demain, fails + 1 (entrée créée si besoin) ;
 * • bonne réponse sur une erreur suivie : au plus une avance par jour, avant l'échéance si elle n'a pas
 *   changé aujourd'hui ; boîte 1 → 2 (J+3), 2 → 3 (J+7), 3 → sortie de la liste,
 *   sauf en examen blanc où elle reste en boîte 3 (J+7).
 */
export function applyAnswerToMistakes(s: ProgressState, r: AnswerResult, mode: QuizMode, today: string): MistakeOutcome {
  if (r.skipped) return 'none';
  const entry = s.mistakes[r.questionId];
  if (!r.correct) {
    s.mistakes[r.questionId] = { box: 1, due: cappedDue(s, addDays(today, BOX_DAYS[1]), today), fails: (entry?.fails ?? 0) + 1, last: today };
    return 'added';
  }
  if (!entry) return 'none';
  // Déjà créée ou promue aujourd'hui : on attend un autre jour (pas de sortie en rejouant la même série).
  if (entry.last === today) return 'none';
  if (entry.box === 3 && mode !== 'exam') {
    delete s.mistakes[r.questionId];
    return 'mastered';
  }
  const box = entry.box === 1 ? 2 : 3;
  s.mistakes[r.questionId] = { ...entry, box, due: cappedDue(s, addDays(today, BOX_DAYS[box]), today), last: today };
  return 'promoted';
}

const plural = (n: number, word: string) => `${n} ${word}${n > 1 ? 's' : ''}`;

/**
 * Fin d'un quiz, quel que soit le mode. Règles anti-rejeu (clés de today.done) :
 * • défi du jour rejoué, ou commencé un autre jour (`opts.day`, ex. avant minuit) : 0 XP, compteurs et score officiel inchangés ;
 * • quiz de chapitre refait le même jour : XP des bonnes réponses ÷ 2, sans bonus ni compteurs ;
 * • examen blanc : bonus et compteur une fois par jour et par matière ; sans aucune réponse, ni note ni bonus.
 * Toute session avec au moins une réponse prolonge la série, même à 0 XP.
 * `opts.maxCombo` : plus longue suite de bonnes réponses du quiz (quête « Enchaîne 5 bonnes réponses ») ;
 * à défaut, elle est estimée sur l'ordre des résultats.
 */
export function quizGain(
  prev: ProgressState,
  mode: QuizMode,
  results: AnswerResult[],
  opts: { chapterId?: string; subjectId?: string; day?: string; maxCombo?: number } = {},
  today = dayKey(),
) {
  const done = rollDay(prev, today).today.done;
  const correct = results.filter((r) => r.correct).length;
  const total = results.length;
  const percent = total ? Math.round((correct / total) * 100) : 0;
  const answered = results.some((r) => !r.skipped);

  const challengeAlreadyDone = rollDay(prev, today).today.challengeDone;
  // Défi commencé un autre jour (ex. avant minuit) : ce sont les questions de ce jour-là, il ne compte pas.
  const dailyStale = mode === 'daily' && !!opts.day && opts.day !== today;
  const dailyReplay = mode === 'daily' && (challengeAlreadyDone || dailyStale);
  const chapterReplay = mode === 'chapter' && !!opts.chapterId && done.includes(`quiz:${opts.chapterId}`);
  const examRepeat = mode === 'exam' && !!opts.subjectId && done.includes(`exam:${opts.subjectId}`);
  // Examen où aucune question n'a été traitée (chrono expiré) : ni note, ni bonus, ni série.
  const examBlank = mode === 'exam' && !answered;
  const counted = !dailyReplay && !chapterReplay && !examRepeat && !examBlank;

  const messages: string[] = [];
  let xp = 0;
  let perfect = false;
  if (dailyReplay) {
    messages.push(dailyStale ? '🔁 Entraînement : ce défi date d’hier, pas d’XP' : '🔁 Entraînement : pas d’XP pour un défi déjà relevé');
  } else if (chapterReplay) {
    xp = Math.floor((correct * XP.correctAnswer) / 2);
    messages.push(`${plural(correct, 'bonne')} ${correct > 1 ? 'réponses' : 'réponse'} : +${xp} XP`);
    messages.push('XP réduits de moitié : tu as déjà fait ce quiz aujourd’hui');
  } else {
    xp = correct * XP.correctAnswer;
    messages.push(`${plural(correct, 'bonne')} ${correct > 1 ? 'réponses' : 'réponse'} : +${xp} XP`);
    perfect = total >= 5 && correct === total;
    if (perfect) {
      xp += XP.perfectQuiz;
      messages.push(`💯 Sans faute : +${XP.perfectQuiz} XP`);
    }
  }
  if (mode === 'daily' && !dailyReplay) {
    xp += XP.dailyChallenge;
    messages.push(`🎯 Défi du jour terminé : +${XP.dailyChallenge} XP`);
  }
  if (mode === 'exam') {
    if (examBlank) {
      messages.push('📝 Examen non traité : pas de note ni de bonus');
    } else if (examRepeat) {
      messages.push('📝 Examen blanc terminé : bonus déjà reçu aujourd’hui pour cette matière');
    } else {
      xp += XP.mockExam;
      messages.push(`📝 Examen blanc terminé : +${XP.mockExam} XP`);
    }
  }

  const events: ActivityEvents = {
    kind: 'quiz',
    mode,
    chapterId: opts.chapterId,
    subjectId: opts.subjectId ?? (mode === 'chapter' ? results[0]?.subjectId : undefined),
    correct: results.filter((r) => r.correct && !r.skipped).length,
    percent,
    maxCombo: opts.maxCombo ?? longestRun(results),
    answered,
  };

  return applyGain(
    prev,
    xp,
    messages,
    (s) => {
      if (counted) {
        s.quizCount += 1;
        if (perfect) s.perfectCount += 1;
      }
      if (mode === 'daily' && !dailyReplay) {
        s.today.challengeDone = true;
        s.challengesDone += 1;
        s.dailyResults[today] = { correct, total, grid: results.map((r) => (r.correct ? '✅' : '❌')).join('') };
        markDone(s, 'daily');
      }
      if (mode === 'chapter' && opts.chapterId) {
        s.quizBest[opts.chapterId] = Math.max(s.quizBest[opts.chapterId] ?? 0, percent);
        markDone(s, `quiz:${opts.chapterId}`);
        if (answered && percent >= STAR_QUIZ_PERCENT) {
          const days = s.quiz80Days[opts.chapterId] ?? [];
          if (!days.includes(today)) s.quiz80Days[opts.chapterId] = [...days, today].sort().slice(-QUIZ80_DAYS_KEPT);
        }
      }
      if (mode === 'exam' && opts.subjectId && !examBlank) {
        const note = examNote(correct, total);
        s.examBest[opts.subjectId] = Math.max(s.examBest[opts.subjectId] ?? 0, note);
        s.examHistory[opts.subjectId] = [...(s.examHistory[opts.subjectId] ?? []), { day: today, note }].slice(-EXAM_HISTORY_SIZE);
        s.examLast[opts.subjectId] = results.map((r) => r.questionId);
        markDone(s, `exam:${opts.subjectId}`);
      }
      let mastered = 0;
      let fixed = 0;
      for (const r of results) {
        const outcome = applyAnswerToMistakes(s, r, mode, today);
        if (outcome === 'mastered') mastered++;
        if (outcome === 'mastered' || outcome === 'promoted') fixed++;
        if (!r.skipped) touchSubject(s, r.subjectId);
      }
      events.mistakesFixed = fixed;
      return mastered > 0 ? [`✅ ${plural(mastered, 'question')} ${mastered > 1 ? 'maîtrisées' : 'maîtrisée'} : ${mastered > 1 ? 'elles sortent' : 'elle sort'} de ta liste`] : [];
    },
    today,
    { activity: answered, events },
  );
}

/** Plus longue suite de bonnes réponses dans une liste de résultats. */
export function longestRun(results: AnswerResult[]): number {
  let best = 0;
  let run = 0;
  for (const r of results) {
    run = r.correct ? run + 1 : 0;
    best = Math.max(best, run);
  }
  return best;
}

export type ChestState = 'locked' | 'ready' | 'opened';

/** Coffre du jour : fermé tant que les 3 quêtes ne sont pas faites, prêt, puis ouvert. */
export function chestState(s: ProgressState, today = dayKey()): ChestState {
  const t = rollDay(s, today).today;
  if (t.chestOpened) return 'opened';
  return allQuestsDone(t.quests) ? 'ready' : 'locked';
}

/**
 * Ouverture du coffre du jour (null s'il n'est pas prêt). Contenu tiré d'après le jour et l'examen :
 * un gel de série (si l'élève en a moins de MAX_FREEZES), sinon de 20 à 40 XP.
 */
export function openChest(prev: ProgressState, today = dayKey()) {
  if (!prev.profile || chestState(prev, today) !== 'ready') return null;
  const rand = seededRandom(`${today}:chest:${prev.profile.track}`);
  const freeze = rand() < CHEST_FREEZE_CHANCE && prev.streak.freezes < MAX_FREEZES;
  const xp = freeze ? 0 : CHEST_XP_MIN + Math.floor(rand() * (CHEST_XP_MAX - CHEST_XP_MIN + 1));
  return applyGain(
    prev,
    xp,
    [freeze ? '🎁 Coffre du jour : un gel de série 🧊 ! Ta série est protégée un jour de plus.' : `🎁 Coffre du jour : +${xp} XP`],
    (s) => {
      s.today.chestOpened = true;
      if (freeze) s.streak.freezes += 1;
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
