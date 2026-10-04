// Quêtes du jour : 3 missions variées (une facile, une moyenne, une de variété), renouvelées à minuit.
// Module pur et testé. Les quêtes non faites disparaissent avec la journée, sans message d'échec.
import { getSubject } from '../data/catalog';
import { daysBetween } from './dates';
import type { ProgressState, QuizMode } from './gamification';
import { nextStep } from './plan';
import { seededRandom, shuffle } from './random';
import { activeMistakes, visibleSubjects } from './selectors';
import { STAR_QUIZ_PERCENT, subjectStars } from './stats';

export type QuestKind = 'correct' | 'fiche' | 'flashcards' | 'quiz80' | 'fix' | 'combo' | 'subject' | 'exam';
export type QuestTier = 'facile' | 'moyenne' | 'variete';

export interface Quest {
  /** `${jour}:${kind}` */
  id: string;
  kind: QuestKind;
  tier: QuestTier;
  target: number;
  progress: number;
  done: boolean;
  /** Matière visée (quête « Fais un quiz en … »). */
  subjectId?: string;
}

export const QUEST_KINDS: QuestKind[] = ['correct', 'fiche', 'flashcards', 'quiz80', 'fix', 'combo', 'subject', 'exam'];
export const QUEST_TIERS: QuestTier[] = ['facile', 'moyenne', 'variete'];

/** Objectif de chaque quête. */
export const QUEST_TARGETS: Record<QuestKind, number> = {
  correct: 10,
  fiche: 1,
  flashcards: 1,
  quiz80: 1,
  fix: 3,
  combo: 5,
  subject: 1,
  exam: 1,
};

/** Examen blanc proposé en quête à partir de J-60. */
export const QUEST_EXAM_DAYS = 60;

/** Ce qui s'est passé pendant une activité terminée : sert à faire avancer les quêtes. */
export interface ActivityEvents {
  kind: 'quiz' | 'fiche' | 'flashcards';
  mode?: QuizMode;
  chapterId?: string;
  subjectId?: string;
  /** Bonnes réponses (questions traitées). */
  correct?: number;
  /** Score du quiz en %. */
  percent?: number;
  /** Plus longue suite de bonnes réponses. */
  maxCombo?: number;
  /** Erreurs suivies corrigées (bonne réponse qui fait avancer ou sortir une erreur). */
  mistakesFixed?: number;
  /** Au moins une question traitée. */
  answered?: boolean;
}

/**
 * Quêtes du jour, tirées avec seededRandom(`${day}:quests:${track}`) parmi celles qui ont du sens pour l'élève :
 * • facile : 10 bonnes réponses, un paquet de flashcards, une nouvelle fiche (s'il en reste) ;
 * • moyenne : ≥ 80 % à un quiz de chapitre, 5 bonnes réponses d'affilée, 3 erreurs corrigées (si au moins 3 sont dues),
 *   un examen blanc (à partir de J-60) ;
 * • variété : un quiz dans la matière la moins travaillée (à défaut, une autre quête facile).
 * Vide sans profil.
 */
export function generateQuests(s: ProgressState, day: string): Quest[] {
  if (!s.profile) return [];
  const rand = seededRandom(`${day}:quests:${s.profile.track}`);
  const subjects = visibleSubjects(s.profile).filter((sub) => sub.chapters.length > 0);
  const remaining = subjects.some((sub) => sub.chapters.some((c) => !s.fichesRead[c.id]));
  const due = activeMistakes(s, day).due.length;
  const daysLeft = daysBetween(day, s.profile.examDate);

  const easy: QuestKind[] = ['correct', 'flashcards', ...(remaining ? (['fiche'] as const) : [])];
  const medium: QuestKind[] = [
    'quiz80',
    'combo',
    ...(due >= QUEST_TARGETS.fix ? (['fix'] as const) : []),
    ...(daysLeft >= 0 && daysLeft <= QUEST_EXAM_DAYS ? (['exam'] as const) : []),
  ];
  const pick = (pool: QuestKind[]) => pool[Math.floor(rand() * pool.length)];
  const quest = (kind: QuestKind, tier: QuestTier, subjectId?: string): Quest => ({
    id: `${day}:${kind}`,
    kind,
    tier,
    target: QUEST_TARGETS[kind],
    progress: 0,
    done: false,
    ...(subjectId ? { subjectId } : {}),
  });

  const first = pick(easy);
  const second = pick(medium);
  const quests = [quest(first, 'facile'), quest(second, 'moyenne')];

  // Matière la moins travaillée (part d'étoiles la plus basse) ; à égalité, le tirage du jour départage.
  let weakest: { id: string; ratio: number } | null = null;
  for (const sub of shuffle(subjects, rand)) {
    const { stars, max } = subjectStars(sub, s);
    const ratio = stars / max;
    if (!weakest || ratio < weakest.ratio) weakest = { id: sub.id, ratio };
  }
  if (weakest) quests.push(quest('subject', 'variete', weakest.id));
  else {
    const others = easy.filter((k) => k !== first);
    if (others.length) quests.push(quest(pick(others), 'variete'));
  }
  return quests;
}

/**
 * Fait avancer les quêtes d'après une activité (modifie `quests`) et renvoie celles qui viennent d'être accomplies.
 */
export function advanceQuests(quests: Quest[], ev: ActivityEvents): Quest[] {
  const completed: Quest[] = [];
  const quiz = ev.kind === 'quiz';
  for (const q of quests) {
    if (q.done) continue;
    let progress = q.progress;
    switch (q.kind) {
      case 'correct':
        progress += quiz ? (ev.correct ?? 0) : 0;
        break;
      case 'fiche':
        if (ev.kind === 'fiche') progress += 1;
        break;
      case 'flashcards':
        if (ev.kind === 'flashcards') progress += 1;
        break;
      case 'quiz80':
        if (quiz && ev.mode === 'chapter' && (ev.percent ?? 0) >= STAR_QUIZ_PERCENT) progress += 1;
        break;
      case 'fix':
        progress += quiz ? (ev.mistakesFixed ?? 0) : 0;
        break;
      case 'combo':
        if (quiz) progress = Math.max(progress, ev.maxCombo ?? 0);
        break;
      case 'subject':
        if (quiz && ev.answered && (ev.mode === 'chapter' || ev.mode === 'exam') && ev.subjectId === q.subjectId) progress += 1;
        break;
      case 'exam':
        if (quiz && ev.answered && ev.mode === 'exam') progress += 1;
        break;
    }
    q.progress = Math.min(q.target, progress);
    if (q.progress >= q.target) {
      q.done = true;
      completed.push(q);
    }
  }
  return completed;
}

/** Les quêtes du jour sont toutes faites (faux s'il n'y en a pas). */
export function allQuestsDone(quests: Quest[]): boolean {
  return quests.length > 0 && quests.every((q) => q.done);
}

/** Intitulé d'une quête : « Donne 10 bonnes réponses ». */
export function questTitle(q: Quest): string {
  switch (q.kind) {
    case 'correct':
      return `Donne ${q.target} bonnes réponses`;
    case 'fiche':
      return 'Lis une nouvelle fiche';
    case 'flashcards':
      return 'Termine un paquet de flashcards';
    case 'quiz80':
      return `Obtiens au moins ${STAR_QUIZ_PERCENT} % à un quiz de chapitre`;
    case 'fix':
      return `Corrige ${q.target} erreurs`;
    case 'combo':
      return `Enchaîne ${q.target} bonnes réponses`;
    case 'subject': {
      const subject = q.subjectId ? getSubject(q.subjectId) : undefined;
      return subject ? `Fais un quiz en ${subject.name}` : 'Fais un quiz de chapitre';
    }
    case 'exam':
      return 'Passe un examen blanc';
  }
}

/** Avancement affiché : « 3/10 », « ✓ » une fois faite. */
export function questProgressLabel(q: Quest): string {
  return q.done ? '✓' : `${q.progress}/${q.target}`;
}

/** Libellé de difficulté : « Facile », « Moyenne », « Variété ». */
export function questTierLabel(tier: QuestTier): string {
  return tier === 'facile' ? 'Facile' : tier === 'moyenne' ? 'Moyenne' : 'Variété';
}

/** Écran où mène une ligne de quête (à passer tel quel à router.push). */
export type QuestLink =
  | { pathname: '/quiz'; params: { mode: QuizMode; id?: string } }
  | { pathname: '/fiche/[id]'; params: { id: string } }
  | { pathname: '/flashcards/[id]'; params: { id: string } }
  | { pathname: '/matiere/[id]'; params: { id: string } }
  | { pathname: '/matieres' }
  | { pathname: '/defis' };

/** Écran concerné par une quête, d'après l'état de l'élève (prochaine fiche, dernier chapitre lu…). */
export function questLink(q: Quest, s: ProgressState, today: string): QuestLink {
  const subjects = visibleSubjects(s.profile);
  const step = nextStep(s, subjects, today);
  // Chapitre lu le plus récemment, pour les flashcards.
  const lastRead = subjects
    .flatMap((sub) => sub.chapters.filter((c) => s.fichesRead[c.id]).map((c) => ({ id: c.id, day: s.fichesRead[c.id] })))
    .sort((a, b) => b.day.localeCompare(a.day))[0];
  switch (q.kind) {
    case 'correct':
    case 'combo':
      return { pathname: '/quiz', params: { mode: 'express' } };
    case 'fiche':
      return step?.kind === 'fiche' ? { pathname: '/fiche/[id]', params: { id: step.chapterId } } : { pathname: '/matieres' };
    case 'flashcards':
      return lastRead ? { pathname: '/flashcards/[id]', params: { id: lastRead.id } } : { pathname: '/matieres' };
    case 'quiz80':
      return step?.kind === 'quiz' ? { pathname: '/quiz', params: { mode: 'chapter', id: step.chapterId } } : { pathname: '/matieres' };
    case 'fix':
      return { pathname: '/quiz', params: { mode: 'review' } };
    case 'subject':
      return q.subjectId ? { pathname: '/matiere/[id]', params: { id: q.subjectId } } : { pathname: '/matieres' };
    case 'exam': {
      // Examen blanc dans une matière jamais passée, sinon dans celle dont la meilleure note est la plus basse.
      const withQuiz = subjects.filter((sub) => sub.chapters.length > 0);
      const target =
        withQuiz.find((sub) => s.examBest[sub.id] === undefined) ??
        [...withQuiz].sort((a, b) => (s.examBest[a.id] ?? 0) - (s.examBest[b.id] ?? 0))[0];
      return target ? { pathname: '/quiz', params: { mode: 'exam', id: target.id } } : { pathname: '/defis' };
    }
  }
}
