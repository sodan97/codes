// Plan de révision : prochaine étape conseillée et rythme jusqu'au jour J.
// Module pur (testé) : les matières sont passées en paramètre (matières visibles de l'élève).
import type { Chapter, Subject } from '../data/types';
import { daysBetween } from './dates';
import type { ProgressState } from './gamification';

type PlanState = Pick<ProgressState, 'fichesRead' | 'quizBest'>;

export interface NextStep {
  kind: 'quiz' | 'fiche';
  chapterId: string;
  subjectId: string;
  title: string;
  reason: string;
}

/** Sous ce score, une fiche lue mérite d'abord un quiz. */
const QUIZ_AFTER_READ_BELOW = 50;

/**
 * Prochaine étape conseillée :
 * a) la fiche lue le plus récemment dont le quiz n'est pas fait ou sous 50 % → quiz ;
 * b) sinon le premier chapitre non lu (ordre du programme) de la matière la moins avancée → fiche ;
 * c) sinon (tout est lu) le chapitre au meilleur score le plus bas, s'il est sous 100 % → quiz ;
 * d) sinon rien.
 * `today` est gardé dans la signature pour les évolutions (séance du jour), il n'influe pas encore.
 */
export function nextStep(state: PlanState, subjects: Subject[], today: string): NextStep | null {
  const chapters = subjects.flatMap((subject) => subject.chapters.map((chapter) => ({ subject, chapter })));
  const step = (kind: NextStep['kind'], subject: Subject, chapter: Chapter, reason: string): NextStep => ({
    kind,
    chapterId: chapter.id,
    subjectId: subject.id,
    title: chapter.title,
    reason,
  });

  // a) Fiche lue à tester. À date égale, l'ordre du catalogue départage (tri stable).
  const toTest = chapters
    .filter(({ chapter }) => {
      const best = state.quizBest[chapter.id];
      return !!state.fichesRead[chapter.id] && (best === undefined || best < QUIZ_AFTER_READ_BELOW);
    })
    .sort((x, y) => state.fichesRead[y.chapter.id].localeCompare(state.fichesRead[x.chapter.id]));
  if (toTest.length) return step('quiz', toTest[0].subject, toTest[0].chapter, 'Fiche lue : teste-toi maintenant');

  // b) Matière la moins avancée (fiches lues / total), à égalité l'ordre du catalogue.
  let behind: { subject: Subject; ratio: number } | null = null;
  for (const subject of subjects) {
    const total = subject.chapters.length;
    const read = subject.chapters.filter((c) => state.fichesRead[c.id]).length;
    if (total === 0 || read === total) continue;
    const ratio = read / total;
    if (!behind || ratio < behind.ratio) behind = { subject, ratio };
  }
  if (behind) {
    const index = behind.subject.chapters.findIndex((c) => !state.fichesRead[c.id]);
    const chapter = behind.subject.chapters[index];
    return step('fiche', behind.subject, chapter, `Chapitre ${index + 1} de ${behind.subject.name}`);
  }

  // c) Tout est lu : le chapitre le plus fragile.
  let weakest: (typeof chapters)[number] | null = null;
  for (const c of chapters) {
    if (!weakest || (state.quizBest[c.chapter.id] ?? 0) < (state.quizBest[weakest.chapter.id] ?? 0)) weakest = c;
  }
  if (weakest && (state.quizBest[weakest.chapter.id] ?? 0) < 100) {
    return step('quiz', weakest.subject, weakest.chapter, 'Ton chapitre le plus fragile');
  }
  return null;
}

export type Phase = 'Découverte' | 'Consolidation' | 'Sprint final' | 'Veille' | 'Jour J' | 'Passé';

export function phaseOf(daysLeft: number): Phase {
  if (daysLeft > 90) return 'Découverte';
  if (daysLeft >= 30) return 'Consolidation';
  if (daysLeft >= 2) return 'Sprint final';
  if (daysLeft === 1) return 'Veille';
  if (daysLeft === 0) return 'Jour J';
  return 'Passé';
}

export interface Pace {
  daysLeft: number;
  /** Fiches pas encore lues dans les matières données. */
  remaining: number;
  weeksLeft: number;
  /** Fiches à lire par semaine pour avoir tout vu le jour J. */
  perWeek: number;
  phase: Phase;
}

/** Rythme de lecture des fiches jusqu'à l'examen. */
export function pace(state: PlanState, subjects: Subject[], examDate: string, today: string): Pace {
  const daysLeft = daysBetween(today, examDate);
  const remaining = subjects.reduce((n, s) => n + s.chapters.filter((c) => !state.fichesRead[c.id]).length, 0);
  const weeksLeft = Math.max(1, Math.ceil(daysLeft / 7));
  return { daysLeft, remaining, weeksLeft, perWeek: Math.ceil(remaining / weeksLeft), phase: phaseOf(daysLeft) };
}
