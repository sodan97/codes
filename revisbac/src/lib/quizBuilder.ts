import { getChapter, getSubject, getSubjects, isOptional, questionsOfSubject, type QuestionRef } from '../data/catalog';
import type { Chapter, Subject, TrackId } from '../data/types';
import type { MistakeEntry, QuizMode } from './gamification';
import { seededRandom, shuffle } from './random';
import { dueMistakes } from './selectors';

export const DAILY_SIZE = 10;
export const EXAM_SIZE = 20;
export const EXAM_SECONDS_PER_QUESTION = 45;
export const REVIEW_SIZE = 15;
export const EXPRESS_SIZE = 5;
/** Révision express : au plus 2 questions par chapitre. */
const EXPRESS_PER_CHAPTER = 2;
/** Révision express : au plus 2 erreurs dues. */
const EXPRESS_MISTAKES = 2;

export interface QuizSession {
  title: string;
  color?: string;
  questions: QuestionRef[];
  /** Durée limite en secondes (examen blanc). */
  timeLimit?: number;
}

/** Ce dont buildQuiz a besoin de l'état de l'élève (voir selectors.quizContext). */
export interface QuizContext {
  track: TrackId;
  /** Matières visibles (examen de l'élève, sans les matières masquées). */
  subjects: Subject[];
  mistakes: Record<string, MistakeEntry>;
  fichesRead: Record<string, string>;
  quizBest: Record<string, number>;
  examLast: Record<string, string[]>;
  today: string;
}

/** Construit la liste de questions d'une session selon le mode. */
export function buildQuiz(mode: QuizMode, id: string | undefined, ctx: QuizContext): QuizSession {
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
      const questions = shuffle(examQuestions(subject, ctx.examLast[subject.id] ?? []));
      return {
        title: `Examen blanc · ${subject.name}`,
        color: subject.color,
        questions,
        timeLimit: questions.length * EXAM_SECONDS_PER_QUESTION,
      };
    }
    case 'review': {
      // Seulement les erreurs arrivées à échéance, les plus urgentes d'abord, puis mélangées.
      const { due } = dueMistakes(ctx.mistakes, ctx.subjects, ctx.today);
      return { title: 'Revoir mes erreurs', questions: shuffle(due.slice(0, REVIEW_SIZE)) };
    }
    case 'express':
      return { title: 'Révision express', questions: buildExpress(ctx) };
    case 'daily': {
      // Même défi pour tous les élèves d'un même examen, le même jour : tiré sur le tronc commun
      // (sans LV2), quelles que soient les matières masquées par chacun.
      const rand = seededRandom(`${ctx.today}:${ctx.track}`);
      const common = getSubjects(ctx.track).filter((s) => !isOptional(s.id));
      const pools = shuffle(common, rand).map((s) => shuffle(questionsOfSubject(s), rand));
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

/**
 * Nombre de questions par chapitre pour un examen blanc, proportionnel à la taille du chapitre :
 * max(1, round(size × n / total)), réajusté pour totaliser min(size, total).
 * S'il y a plus de chapitres que de questions à tirer, les plus petits chapitres sont laissés de côté.
 */
export function examQuotas(counts: number[], size = EXAM_SIZE): number[] {
  const total = counts.reduce((a, b) => a + b, 0);
  if (total === 0) return counts.map(() => 0);
  const target = Math.min(size, total);
  const exact = counts.map((n) => (size * n) / total);
  const quotas = counts.map((n, i) => (n === 0 ? 0 : Math.min(n, Math.max(1, Math.round(exact[i])))));
  let sum = quotas.reduce((a, b) => a + b, 0);
  while (sum > target) {
    // On retire au chapitre le plus servi par rapport à sa part exacte, en lui laissant une question…
    let i = pickIndex(quotas, (j) => quotas[j] > 1, (j) => quotas[j] - exact[j]);
    // … sauf s'il y a plus de chapitres que de questions.
    if (i < 0) i = pickIndex(quotas, (j) => quotas[j] > 0, (j) => -exact[j]);
    quotas[i]--;
    sum--;
  }
  while (sum < target) {
    const i = pickIndex(quotas, (j) => quotas[j] < counts[j], (j) => exact[j] - quotas[j]);
    quotas[i]++;
    sum++;
  }
  return quotas;
}

/** Indice qui maximise `score` parmi ceux qui vérifient `ok` (-1 si aucun). */
function pickIndex(items: unknown[], ok: (i: number) => boolean, score: (i: number) => number): number {
  let best = -1;
  for (let i = 0; i < items.length; i++) if (ok(i) && (best < 0 || score(i) > score(best))) best = i;
  return best;
}

/**
 * Examen blanc : tirage stratifié par chapitre (voir examQuotas). Dans chaque chapitre,
 * les questions absentes du dernier examen de la matière passent en premier.
 */
function examQuestions(subject: Subject, lastIds: string[]): QuestionRef[] {
  const last = new Set(lastIds);
  const quotas = examQuotas(subject.chapters.map((c) => c.quiz.length));
  return subject.chapters.flatMap((chapter, i) => {
    const fresh = shuffle(chapter.quiz.filter((q) => !last.has(q.id)));
    const seen = shuffle(chapter.quiz.filter((q) => last.has(q.id)));
    return [...fresh, ...seen].slice(0, quotas[i]).map((question) => ({ question, chapter, subject }));
  });
}

/**
 * Révision express : 5 questions, au plus 2 par chapitre, sans doublon.
 * Priorités : erreurs dues, chapitres lus encore fragiles (quiz < 80 %), autres chapitres lus,
 * puis à défaut le premier chapitre de chaque matière visible.
 */
function buildExpress(ctx: QuizContext): QuestionRef[] {
  const picked: QuestionRef[] = [];
  const perChapter = new Map<string, number>();
  const take = (refs: QuestionRef[], max = EXPRESS_SIZE) => {
    let n = 0;
    for (const ref of refs) {
      if (picked.length >= EXPRESS_SIZE || n >= max) return;
      const used = perChapter.get(ref.chapter.id) ?? 0;
      if (used >= EXPRESS_PER_CHAPTER || picked.some((p) => p.question.id === ref.question.id)) continue;
      perChapter.set(ref.chapter.id, used + 1);
      picked.push(ref);
      n++;
    }
  };
  const questionsOf = (subject: Subject, chapter: Chapter): QuestionRef[] => chapter.quiz.map((question) => ({ question, chapter, subject }));

  take(dueMistakes(ctx.mistakes, ctx.subjects, ctx.today).due, EXPRESS_MISTAKES);
  const read = ctx.subjects.flatMap((s) => s.chapters.filter((c) => ctx.fichesRead[c.id]).map((c) => ({ s, c })));
  const fragile = read.filter(({ c }) => (ctx.quizBest[c.id] ?? 0) < 80);
  take(shuffle(fragile.flatMap(({ s, c }) => questionsOf(s, c))));
  take(shuffle(read.flatMap(({ s, c }) => questionsOf(s, c))));
  take(shuffle(ctx.subjects.flatMap((s) => (s.chapters[0] ? questionsOf(s, s.chapters[0]) : []))));
  // Dernier recours (une seule matière visible, par exemple) : n'importe quel chapitre visible.
  take(shuffle(ctx.subjects.flatMap((s) => s.chapters.flatMap((c) => questionsOf(s, c)))));
  return shuffle(picked);
}
