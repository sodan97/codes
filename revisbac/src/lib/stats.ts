import type { Subject } from '../data/types';
import type { ProgressState } from './gamification';

/** Avancement d'une matière : fiches lues et maîtrise moyenne aux quiz (0 → 1). */
export function subjectProgress(subject: Subject, state: ProgressState) {
  const total = subject.chapters.length;
  const read = subject.chapters.filter((c) => state.fichesRead[c.id]).length;
  const mastery = total ? subject.chapters.reduce((sum, c) => sum + (state.quizBest[c.id] ?? 0), 0) / (total * 100) : 0;
  return { total, read, mastery };
}
