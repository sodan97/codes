// Lectures de l'état filtrées par le catalogue actuel.
// Les identifiants inconnus (contenu retiré ou renommé pendant la relecture, autre examen)
// restent dans la sauvegarde mais ne sont ni comptés ni proposés.
import { getQuestion, getSubjects, type QuestionRef } from '../data/catalog';
import type { Subject } from '../data/types';
import type { MistakeEntry, Profile, ProgressState } from './gamification';
import type { QuizContext } from './quizBuilder';
import { chapterDueCount } from './srs';

/** Matières de l'examen de l'élève, sans celles qu'il a masquées. */
export function visibleSubjects(profile: Profile | null): Subject[] {
  if (!profile) return [];
  return getSubjects(profile.track).filter((s) => !profile.hiddenSubjects.includes(s.id));
}

export interface ActiveMistakes {
  /** Erreurs à revoir aujourd'hui (échéance passée), les plus anciennes d'abord puis les plus ratées. */
  due: QuestionRef[];
  /** Erreurs suivies dont l'échéance n'est pas encore arrivée. */
  waiting: number;
}

/** Erreurs suivies, limitées aux questions existantes des matières données. */
export function dueMistakes(mistakes: Record<string, MistakeEntry>, subjects: Subject[], today: string): ActiveMistakes {
  const allowed = new Set(subjects.map((s) => s.id));
  const due: { ref: QuestionRef; entry: MistakeEntry }[] = [];
  let waiting = 0;
  for (const [id, entry] of Object.entries(mistakes)) {
    const ref = getQuestion(id);
    if (!ref || !allowed.has(ref.subject.id)) continue;
    if (entry.due <= today) due.push({ ref, entry });
    else waiting++;
  }
  due.sort((a, b) => (a.entry.due === b.entry.due ? b.entry.fails - a.entry.fails : a.entry.due < b.entry.due ? -1 : 1));
  return { due: due.map((d) => d.ref), waiting };
}

/** Erreurs de l'élève sur son examen et ses matières visibles. */
export function activeMistakes(state: ProgressState, today: string): ActiveMistakes {
  return dueMistakes(state.mistakes, visibleSubjects(state.profile), today);
}

/** Nombre de fiches lues parmi les chapitres existants de l'examen de l'élève. */
export function readCount(state: ProgressState, profile: Profile | null = state.profile): number {
  if (!profile) return 0;
  return getSubjects(profile.track).reduce((n, s) => n + s.chapters.filter((c) => state.fichesRead[c.id]).length, 0);
}

/** Contexte de construction d'un quiz (null sans profil). */
export function quizContext(state: ProgressState, today: string): QuizContext | null {
  if (!state.profile) return null;
  return {
    track: state.profile.track,
    subjects: visibleSubjects(state.profile),
    mistakes: state.mistakes,
    fichesRead: state.fichesRead,
    quizBest: state.quizBest,
    examLast: state.examLast,
    today,
  };
}

export interface DueCards {
  /** Flashcards à revoir aujourd'hui, dans les matières visibles. */
  total: number;
  /** Détail par chapitre, les plus chargés d'abord (pour ouvrir directement le bon paquet). */
  byChapter: { chapterId: string; subjectId: string; count: number }[];
}

/** Flashcards déjà vues dont l'échéance est arrivée (cartes du contenu actuel, matières visibles). */
export function dueCards(state: ProgressState, today: string): DueCards {
  const byChapter: DueCards['byChapter'] = [];
  for (const subject of visibleSubjects(state.profile)) {
    for (const chapter of subject.chapters) {
      const count = chapterDueCount(chapter.id, chapter.flashcards, state.cards, today);
      if (count > 0) byChapter.push({ chapterId: chapter.id, subjectId: subject.id, count });
    }
  }
  byChapter.sort((a, b) => b.count - a.count);
  return { total: byChapter.reduce((n, c) => n + c.count, 0), byChapter };
}

/** « 🧠 3 cartes à revoir » (null s'il n'y en a pas). */
export function dueCardsText(n: number): string | null {
  if (n <= 0) return null;
  return `🧠 ${n} carte${n > 1 ? 's' : ''} à revoir`;
}
