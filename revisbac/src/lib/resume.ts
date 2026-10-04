// Reprise d'un quiz ou d'un examen interrompu (appli fermée par Android, appel, car rapide).
// Module pur et testé : la lecture et l'écriture sur le téléphone sont dans state/session.ts.
import { dayKey, isDayKey } from './dates';
import type { AnswerResult, ProgressState, QuizMode } from './gamification';
import { rebuildSession, type QuizSession } from './quizBuilder';

/** Une session de plus de 24 h n'est plus proposée. */
export const SESSION_MAX_AGE_MS = 24 * 60 * 60 * 1000;
/** Examen dont le temps est écoulé depuis plus de 30 min : noté sans proposition de reprise. */
export const EXAM_GRACE_MS = 30 * 60 * 1000;
/** Examen terminé appli fermée : noté à la réouverture pendant 7 jours, ensuite effacé. */
export const EXAM_GRADE_MAX_AGE_MS = 7 * 24 * 60 * 60 * 1000;
/** Une fiche ouverte il y a plus de 7 jours n'est plus proposée. */
export const LAST_FICHE_MAX_AGE_MS = 7 * 24 * 60 * 60 * 1000;

const MODES: QuizMode[] = ['chapter', 'daily', 'exam', 'review', 'express'];

/** Session de quiz en cours, enregistrée après chaque réponse validée (clé 'revisbac/session/v1'). */
export interface SavedSession {
  mode: QuizMode;
  /** Chapitre (quiz de chapitre) ou matière (examen blanc), null sinon. */
  id: string | null;
  /**
   * Jour du début (AAAA-MM-JJ) : un défi d'un autre jour n'est pas repris (le défi du jour l'emporte) ;
   * commencé avant minuit et fini après sans fermer l'appli, il reste celui de ce jour-là.
   */
  day: string;
  /** Questions de la session, dans l'ordre de session.questions. */
  questionIds: string[];
  /** File des questions restantes (indices dans questionIds). */
  order: number[];
  /** Réponse donnée à chaque question (true/false), null tant qu'elle n'est pas traitée. */
  answers: (boolean | null)[];
  /** Début de la session (ms depuis 1970). */
  startedAt: number;
  /** Examen blanc : fin du temps imparti (ms depuis 1970). */
  deadline?: number;
  /** Suite de bonnes réponses en cours et record de la session. */
  combo: number;
  maxCombo: number;
}

const isObject = (v: unknown): v is Record<string, unknown> => typeof v === 'object' && v !== null && !Array.isArray(v);
const isCount = (v: unknown): v is number => typeof v === 'number' && Number.isInteger(v) && v >= 0;
const isTime = (v: unknown): v is number => typeof v === 'number' && Number.isFinite(v) && v > 0;

/** Vérifie une session relue (null si elle est incohérente ou abîmée). */
export function parseSavedSession(raw: unknown): SavedSession | null {
  if (!isObject(raw)) return null;
  const mode = MODES.find((m) => m === raw.mode);
  if (!mode || !isDayKey(raw.day) || !isTime(raw.startedAt)) return null;
  const id = typeof raw.id === 'string' && raw.id ? raw.id : null;
  if ((mode === 'chapter' || mode === 'exam') && !id) return null;
  const { questionIds, order, answers } = raw;
  if (!Array.isArray(questionIds) || questionIds.length === 0 || !questionIds.every((q) => typeof q === 'string')) return null;
  const n = questionIds.length;
  if (!Array.isArray(answers) || answers.length !== n || !answers.every((a) => a === null || typeof a === 'boolean')) return null;
  if (!Array.isArray(order) || !order.every((i) => isCount(i) && i < n) || new Set(order).size !== order.length) return null;
  // Les questions restantes sont exactement celles sans réponse.
  if (order.some((i) => answers[i] !== null) || answers.filter((a) => a === null).length !== order.length) return null;
  if (raw.deadline !== undefined && !isTime(raw.deadline)) return null;
  if (mode === 'exam' && raw.deadline === undefined) return null;
  return {
    mode,
    id,
    day: raw.day,
    questionIds: [...questionIds] as string[],
    order: [...order] as number[],
    answers: [...answers] as (boolean | null)[],
    startedAt: raw.startedAt,
    ...(isTime(raw.deadline) ? { deadline: raw.deadline } : {}),
    combo: isCount(raw.combo) ? raw.combo : 0,
    maxCombo: isCount(raw.maxCombo) ? raw.maxCombo : 0,
  };
}

/** La session enregistrée correspond-elle au quiz ouvert (même mode, même chapitre ou matière) ? */
export function matchesSession(saved: SavedSession, mode: QuizMode, id: string | undefined): boolean {
  return saved.mode === mode && saved.id === (id ?? null);
}

/** Nombre de questions déjà traitées. */
export function answeredCount(saved: SavedSession): number {
  return saved.answers.filter((a) => a !== null).length;
}

/**
 * Que faire d'une session retrouvée à l'ouverture ?
 * • 'resume' : proposer « Reprendre (k/n) » ou « Recommencer » (moins de 24 h, temps d'examen non écoulé) ;
 * • 'grade' : examen dont le temps s'est écoulé pendant que l'appli était fermée (moins de 30 min),
 *   ou session dont toutes les questions sont traitées : noter tout de suite et montrer le résultat ;
 * • 'gradeSilently' : examen terminé depuis plus de 30 min (commencé il y a moins de 7 jours) :
 *   le noter sans proposer de reprise ;
 * • 'discard' : trop ancienne (hors examen terminé), défi du jour d'un autre jour, ou rien à noter : l'effacer.
 * Les questions non traitées d'un examen noté comptent comme non traitées (skipped).
 */
export type ResumeDecision = 'resume' | 'grade' | 'gradeSilently' | 'discard';

export function resumeDecision(saved: SavedSession, now: number, today: string = dayKey(new Date(now))): ResumeDecision {
  // Défi d'un autre jour : il ne rapporterait plus rien, le défi d'aujourd'hui l'emporte.
  if (saved.mode === 'daily' && saved.day !== today) return 'discard';
  const tooOld = now - saved.startedAt > SESSION_MAX_AGE_MS || now < saved.startedAt - SESSION_MAX_AGE_MS;
  const answered = answeredCount(saved);
  if (saved.mode === 'exam' && saved.deadline !== undefined && now >= saved.deadline) {
    // Rien de traité : ni note ni bonus, inutile de le noter.
    // Examen trop ancien (plus de 7 jours) : effacé.
    if (answered === 0 || now - saved.startedAt > EXAM_GRADE_MAX_AGE_MS) return 'discard';
    return now - saved.deadline <= EXAM_GRACE_MS ? 'grade' : 'gradeSilently';
  }
  if (tooOld) return 'discard';
  if (saved.order.length === 0) return 'grade';
  if (answered === 0 && saved.mode !== 'exam') return 'discard';
  return 'resume';
}

/**
 * Examen blanc en cours (avec des réponses) alors que l'élève ouvre un autre quiz : une seule session est
 * enregistrée, la première réponse de l'autre quiz l'écraserait. L'écran demande d'abord de le reprendre ou de le noter.
 */
export function isPendingOtherExam(saved: SavedSession, mode: QuizMode, id: string | undefined, now: number): boolean {
  if (saved.mode !== 'exam' || matchesSession(saved, mode, id)) return false;
  return answeredCount(saved) > 0 && resumeDecision(saved, now) !== 'discard';
}

/** Temps d'examen restant en ms (0 si écoulé, null hors examen). */
export function remainingMs(saved: SavedSession, now: number): number | null {
  return saved.deadline === undefined ? null : Math.max(0, saved.deadline - now);
}

/** Session reconstruite, prête à être reprise par l'écran de quiz. */
export interface RestoredSession {
  saved: SavedSession;
  session: QuizSession;
  /** File des questions restantes (indices dans session.questions). */
  order: number[];
  /** Réponses alignées sur session.questions (null = pas encore traitée). */
  answers: (boolean | null)[];
  answered: number;
  total: number;
}

/**
 * Reconstruit la session (rebuildSession) et recale file et réponses sur les questions encore présentes
 * dans le contenu. Null s'il ne reste aucune question.
 */
export function restoreSession(saved: SavedSession): RestoredSession | null {
  const session = rebuildSession(saved.mode, saved.id ?? undefined, saved.questionIds);
  if (session.questions.length === 0) return null;
  const newIndex = new Map(session.questions.map((q, i) => [q.question.id, i]));
  const answers = session.questions.map((q) => saved.answers[saved.questionIds.indexOf(q.question.id)] ?? null);
  const order = saved.order.map((i) => newIndex.get(saved.questionIds[i])).filter((i): i is number => i !== undefined);
  const answered = answers.filter((a) => a !== null).length;
  return { saved, session, order, answers, answered, total: session.questions.length };
}

/** Résultats d'une session reprise, à passer à finishQuiz : les questions sans réponse sont non traitées (skipped). */
export function sessionResults(restored: RestoredSession): AnswerResult[] {
  return restored.session.questions.map((q, i) => {
    const answer = restored.answers[i];
    const base = { questionId: q.question.id, subjectId: q.subject.id };
    return answer === null ? { ...base, correct: false, skipped: true } : { ...base, correct: answer };
  });
}

/** Tout ce qu'il faut pour noter une session sans l'écran de quiz : finishQuiz(g.mode, g.results, g.opts). */
export function gradeSaved(restored: RestoredSession): {
  mode: QuizMode;
  results: AnswerResult[];
  opts: { chapterId?: string; subjectId?: string; day: string; maxCombo: number };
} {
  const { saved } = restored;
  return {
    mode: saved.mode,
    results: sessionResults(restored),
    opts: {
      ...(saved.mode === 'chapter' && saved.id ? { chapterId: saved.id } : {}),
      ...(saved.mode === 'exam' && saved.id ? { subjectId: saved.id } : {}),
      day: saved.day,
      maxCombo: saved.maxCombo,
    },
  };
}

/** Carte de l'accueil : « ▶ Reprendre ton quiz : Les suites, 6/10 ». */
export function resumeCardText(restored: RestoredSession): string {
  return `▶ Reprendre ton quiz : ${restored.session.title}, ${restored.answered}/${restored.total}`;
}

/** Bouton de l'écran de quiz : « Reprendre (6/10) ». */
export function resumeButtonLabel(restored: RestoredSession): string {
  return `Reprendre (${restored.answered}/${restored.total})`;
}

/** Message de sortie d'un quiz, maintenant que la progression est gardée. */
export const QUIT_KEPT_MESSAGE = 'Ta progression est gardée, tu pourras reprendre.';

/** Dernière fiche ouverte (clé 'revisbac/last-fiche/v1'). */
export interface LastFiche {
  chapterId: string;
  /** Ouverture (ms depuis 1970). */
  openedAt: number;
}

export function parseLastFiche(raw: unknown): LastFiche | null {
  if (!isObject(raw) || typeof raw.chapterId !== 'string' || !raw.chapterId || !isTime(raw.openedAt)) return null;
  return { chapterId: raw.chapterId, openedAt: raw.openedAt };
}

/**
 * Fiche à proposer « Continuer ta fiche » : la dernière ouverte, si elle n'est pas encore lue jusqu'au bout,
 * appartient à une matière visible et date de moins de 7 jours. `chapterIds` : chapitres visibles de l'élève.
 */
export function ficheToContinue(
  last: LastFiche | null,
  state: Pick<ProgressState, 'fichesRead'>,
  chapterIds: ReadonlySet<string>,
  now: number,
): string | null {
  if (!last || state.fichesRead[last.chapterId] || !chapterIds.has(last.chapterId)) return null;
  return now - last.openedAt <= LAST_FICHE_MAX_AGE_MS ? last.chapterId : null;
}
