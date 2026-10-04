import { getChapter, getQuestion, getSubject, getSubjects, questionsOfSubject, type QuestionRef } from '../data/catalog';
import type { TrackId } from '../data/types';
import type { QuizMode } from './gamification';
import { seededRandom, shuffle } from './random';

export const DAILY_SIZE = 10;
export const EXAM_SIZE = 20;
export const EXAM_SECONDS_PER_QUESTION = 45;
export const REVIEW_SIZE = 15;

export interface QuizSession {
  title: string;
  color?: string;
  questions: QuestionRef[];
  /** Durée limite en secondes (examen blanc). */
  timeLimit?: number;
}

/** Construit la liste de questions d'une session selon le mode. */
export function buildQuiz(mode: QuizMode, id: string | undefined, track: TrackId, mistakes: Record<string, number>, today: string): QuizSession {
  switch (mode) {
    case 'chapter': {
      const ref = getChapter(id ?? '');
      if (!ref) return { title: 'Quiz', questions: [] };
      return {
        title: ref.chapter.title,
        color: ref.subject.color,
        questions: shuffle(ref.chapter.quiz).map((question) => ({ question, ...ref })),
      };
    }
    case 'exam': {
      const subject = getSubject(id ?? '');
      if (!subject) return { title: 'Examen blanc', questions: [] };
      const questions = shuffle(questionsOfSubject(subject)).slice(0, EXAM_SIZE);
      return {
        title: `Examen blanc · ${subject.name}`,
        color: subject.color,
        questions,
        timeLimit: questions.length * EXAM_SECONDS_PER_QUESTION,
      };
    }
    case 'review': {
      const refs = Object.keys(mistakes)
        .map(getQuestion)
        .filter((r): r is QuestionRef => !!r);
      return { title: 'Revoir mes erreurs', questions: shuffle(refs).slice(0, REVIEW_SIZE) };
    }
    case 'daily': {
      // Même défi pour tous les élèves d'un même examen, le même jour.
      const rand = seededRandom(`${today}:${track}`);
      const pools = shuffle(getSubjects(track), rand).map((s) => shuffle(questionsOfSubject(s), rand));
      const picked: QuestionRef[] = [];
      for (let round = 0; picked.length < DAILY_SIZE && pools.some((p) => p.length > round); round++) {
        for (const pool of pools) {
          if (pool[round] && picked.length < DAILY_SIZE) picked.push(pool[round]);
        }
      }
      return { title: 'Défi du jour', questions: shuffle(picked, rand) };
    }
  }
}
