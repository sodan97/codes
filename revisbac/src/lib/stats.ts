// Avancement et étoiles de maîtrise. Module pur : les matières sont passées en paramètre.
import type { Subject } from '../data/types';
import { dayKey, daysBetween } from './dates';
import type { ProgressState } from './gamification';

/** Ce dont le calcul des étoiles a besoin de l'état. */
export type StarState = Pick<ProgressState, 'fichesRead' | 'quizBest' | 'quiz80Days' | 'flashcardsDone'>;

/** Score minimal d'un quiz de chapitre pour les étoiles 2 et 3. */
export const STAR_QUIZ_PERCENT = 80;
/** Au-delà, un chapitre à 3 étoiles non pratiqué mérite un petit rappel (étoiles estompées, jamais retirées). */
export const STAR_FADE_DAYS = 21;

export type StarCount = 0 | 1 | 2 | 3;

export interface ChapterStars {
  stars: StarCount;
  /** 3 étoiles, mais rien pratiqué depuis plus de 21 jours : afficher les étoiles estompées. */
  faded: boolean;
  /** « Prochaine étoile : … » (null à 3 étoiles). */
  next: string | null;
  /** « 🔄 petit rappel conseillé » quand les étoiles sont estompées, sinon null. */
  reminder: string | null;
  /** Dernier jour de pratique connu (quiz ≥ 80 % ou paquet de flashcards terminé). */
  lastPractice: string | null;
}

/**
 * Étoiles d'un chapitre, chacune suppose la précédente :
 * ⭐ fiche lue ; ⭐⭐ quiz de chapitre ≥ 80 % ; ⭐⭐⭐ quiz ≥ 80 % sur 2 jours différents et paquet de flashcards terminé.
 */
export function chapterStars(chapterId: string, state: StarState, today: string = dayKey()): ChapterStars {
  const read = !!state.fichesRead[chapterId];
  const quiz80 = (state.quizBest[chapterId] ?? 0) >= STAR_QUIZ_PERCENT;
  const days = state.quiz80Days[chapterId] ?? [];
  const flash = state.flashcardsDone[chapterId];
  const missingDays = Math.max(0, 2 - days.length);

  let stars: StarCount = 0;
  if (read) stars = 1;
  if (read && quiz80) stars = 2;
  if (read && quiz80 && missingDays === 0 && flash) stars = 3;

  let next: string | null = null;
  if (stars === 0) next = 'Prochaine étoile : lis la fiche';
  else if (stars === 1) next = `Prochaine étoile : réussis le quiz à ${STAR_QUIZ_PERCENT} % ou plus`;
  else if (stars === 2) {
    const quiz =
      missingDays === 2
        ? `réussis le quiz à ${STAR_QUIZ_PERCENT} % sur 2 jours différents`
        : missingDays === 1
          ? `refais le quiz à ${STAR_QUIZ_PERCENT} % un autre jour`
          : null;
    const deck = flash ? null : 'termine le paquet de flashcards';
    next = `Prochaine étoile : ${[deck, quiz].filter(Boolean).join(' et ')}`;
  }

  const practice = [days[days.length - 1], flash].filter((d): d is string => !!d).sort();
  const lastPractice = practice.length ? practice[practice.length - 1] : null;
  const faded = stars === 3 && !!lastPractice && daysBetween(lastPractice, today) > STAR_FADE_DAYS;
  return { stars, faded, next, reminder: faded ? '🔄 petit rappel conseillé' : null, lastPractice };
}

/** Étoiles d'une matière (3 par chapitre au plus). */
export function subjectStars(subject: Subject, state: StarState): { stars: number; max: number } {
  const stars = subject.chapters.reduce((n, c) => n + chapterStars(c.id, state).stars, 0);
  return { stars, max: 3 * subject.chapters.length };
}

/**
 * Avancement d'une matière : fiches lues, étoiles et maîtrise (0 → 1).
 * La maîtrise vaut étoiles / (3 × chapitres) : elle progresse avec la régularité, pas seulement avec un bon score.
 */
export function subjectProgress(subject: Subject, state: StarState) {
  const total = subject.chapters.length;
  const read = subject.chapters.filter((c) => state.fichesRead[c.id]).length;
  const { stars, max } = subjectStars(subject, state);
  return { total, read, stars, maxStars: max, mastery: max ? stars / max : 0 };
}

/** « ⭐ 14/24 » */
export function starsLabel(stars: number, max: number): string {
  return `⭐ ${stars}/${max}`;
}

/** Libellé pour le lecteur d'écran : « 2 étoiles sur 3 ». */
export function starsAccessibilityLabel(stars: number, max = 3): string {
  return `${stars} étoile${stars > 1 ? 's' : ''} sur ${max}`;
}
