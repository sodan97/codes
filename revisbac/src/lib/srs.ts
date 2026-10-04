// Répétition espacée des flashcards (Leitner 5 boîtes) : module pur et testé.
import type { Flashcard } from '../data/types';
import { addDays } from './dates';
import { fnv1a } from './hash';
import { shuffle } from './random';

export type CardBox = 1 | 2 | 3 | 4 | 5;

/** Suivi d'une carte : boîte atteinte et jour à partir duquel elle est à revoir. */
export interface CardEntry {
  box: CardBox;
  /** Jour (AAAA-MM-JJ) à partir duquel la carte est à revoir. */
  due: string;
  /** Dernier jour où la carte a été revue. */
  last: string;
}

/** Intervalle (en jours) avant la prochaine révision, selon la boîte atteinte. */
export const CARD_BOX_DAYS: Record<CardBox, number> = { 1: 1, 2: 3, 3: 7, 4: 14, 5: 30 };

const ACCENTS: Record<string, string> = {
  à: 'a', á: 'a', â: 'a', ã: 'a', ä: 'a', å: 'a',
  ç: 'c',
  è: 'e', é: 'e', ê: 'e', ë: 'e',
  ì: 'i', í: 'i', î: 'i', ï: 'i',
  ñ: 'n',
  ò: 'o', ó: 'o', ô: 'o', õ: 'o', ö: 'o',
  ù: 'u', ú: 'u', û: 'u', ü: 'u',
  ý: 'y', ÿ: 'y',
  œ: 'oe', æ: 'ae',
};

/** Ponctuation courante, sans importance pour distinguer deux cartes. */
const PUNCTUATION = /[\s.,;:!?'’"«»()[\]\-–—…/]/;

/**
 * Identifiant lisible d'un texte : minuscules sans accents, mots séparés par des tirets (60 caractères au plus).
 * Sans String.normalize, pour donner le même résultat sous Hermes, Node et dans le navigateur.
 * Les symboles (formules, exposants, flèches) sont perdus dans le slug : une empreinte courte du texte
 * est alors ajoutée, pour que « (a − b)² = ? » et « (a + b)² = ? » restent distincts.
 */
export function slug(text: string): string {
  let out = '';
  let lossy = false;
  for (const ch of text.toLowerCase()) {
    if (/[a-z0-9]/.test(ch)) out += ch;
    else if (ACCENTS[ch]) out += ACCENTS[ch];
    else {
      if (!PUNCTUATION.test(ch)) lossy = true;
      out += '-';
    }
  }
  const base = out.replace(/-+/g, '-').replace(/^-|-$/g, '').slice(0, 60).replace(/-$/, '');
  if (!lossy && base) return base;
  const hash = fnv1a(text).slice(0, 6);
  return base ? `${base}-${hash}` : hash;
}

/** Clé de suivi d'une carte : son id s'il existe, sinon `${chapterId}#${slug(recto)}` (l'ordre des cartes peut changer). */
export function cardKey(chapterId: string, card: Flashcard): string {
  return card.id ? card.id : `${chapterId}#${slug(card.front)}`;
}

/** Échéance plafonnée pour que tout soit revu avant l'épreuve : au plus tard examDate - 2, jamais avant demain. */
export function capDue(due: string, today: string, examDate?: string): string {
  if (!examDate || examDate <= today) return due;
  const limit = addDays(examDate, -2);
  const tomorrow = addDays(today, 1);
  const capped = due < limit ? due : limit;
  return capped < tomorrow ? tomorrow : capped;
}

/**
 * Réponse à une carte :
 * • « À revoir » : boîte 1, à revoir demain ;
 * • « Je savais » : boîte suivante (une carte jamais vue part de la boîte 1), au plus une avance par jour :
 *   une carte déjà revue aujourd'hui (ratée puis retrouvée en fin de paquet, paquet refait) ne bouge plus.
 */
export function scheduleCard(entry: CardEntry | undefined, knew: boolean, today: string, examDate?: string): CardEntry {
  if (!knew) return { box: 1, due: capDue(addDays(today, CARD_BOX_DAYS[1]), today, examDate), last: today };
  if (entry && entry.last === today) return entry;
  const box = Math.min(5, (entry?.box ?? 1) + 1) as CardBox;
  return { box, due: capDue(addDays(today, CARD_BOX_DAYS[box]), today, examDate), last: today };
}

/** Carte déjà vue dont l'échéance est arrivée. */
export function isCardDue(entry: CardEntry | undefined, today: string): boolean {
  return !!entry && entry.due <= today;
}

/**
 * Ordre du paquet d'un chapitre : d'abord les cartes à revoir (les plus en retard, puis les plus fragiles),
 * puis les cartes jamais vues, puis les autres (les plus proches de leur échéance d'abord).
 * Toutes les cartes sont proposées : le paquet reste complet.
 */
export function deckOrder<T extends Flashcard>(
  chapterId: string,
  cards: readonly T[],
  entries: Record<string, CardEntry>,
  today: string,
  rand: () => number = Math.random,
): T[] {
  const withEntry = shuffle(cards, rand).map((card) => ({ card, entry: entries[cardKey(chapterId, card)] }));
  const due = withEntry.filter((c) => isCardDue(c.entry, today));
  const fresh = withEntry.filter((c) => !c.entry);
  const later = withEntry.filter((c) => c.entry && !isCardDue(c.entry, today));
  due.sort((a, b) => (a.entry!.due === b.entry!.due ? a.entry!.box - b.entry!.box : a.entry!.due < b.entry!.due ? -1 : 1));
  later.sort((a, b) => (a.entry!.due < b.entry!.due ? -1 : a.entry!.due > b.entry!.due ? 1 : 0));
  return [...due, ...fresh, ...later].map((c) => c.card);
}

/** Nombre de cartes à revoir aujourd'hui dans un paquet. */
export function chapterDueCount(chapterId: string, cards: readonly Flashcard[], entries: Record<string, CardEntry>, today: string): number {
  return cards.filter((card) => isCardDue(entries[cardKey(chapterId, card)], today)).length;
}
