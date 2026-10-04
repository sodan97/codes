// Modèle de contenu de l'application.
// Tout le contenu pédagogique (fiches, flashcards, quiz) est décrit avec ces types
// pour pouvoir être enrichi plus tard sans toucher au code des écrans.

/** Examen préparé par l'élève. */
export type TrackId = 'bfm' | 'bac-s' | 'bac-l';

/** Un bloc d'une fiche de révision. */
export type FicheBlock =
  /** Paragraphe court. */
  | { kind: 'text'; text: string }
  /** Liste à puces, avec un titre optionnel. */
  | { kind: 'list'; title?: string; items: string[] }
  /** Formule encadrée (notation Unicode : x², √, ∫, ≤, →, etc.). */
  | { kind: 'formula'; label?: string; formula: string; note?: string }
  /** Définition à connaître par cœur. */
  | { kind: 'definition'; term: string; definition: string }
  /** Date clé (histoire) ou repère chronologique. */
  | { kind: 'date'; date: string; event: string }
  /** Méthode / astuce pour l'examen. */
  | { kind: 'tip'; text: string }
  /** Piège fréquent / erreur à éviter. */
  | { kind: 'warning'; text: string }
  /** Exemple ou application. */
  | { kind: 'example'; title?: string; text: string };

export interface FicheSection {
  title: string;
  blocks: FicheBlock[];
}

interface QuestionBase {
  /** Identifiant unique dans toute l'application, ex. "maths-s-suites-q1". */
  id: string;
  /** Explication affichée après la réponse (toujours renseignée : on apprend de ses erreurs). */
  explanation: string;
}

/** Question à choix multiple : une seule bonne réponse. */
export interface QcmQuestion extends QuestionBase {
  type: 'qcm';
  prompt: string;
  choices: string[];
  /** Index de la bonne réponse dans `choices`. */
  answer: number;
}

/** Affirmation à juger vraie ou fausse. */
export interface TrueFalseQuestion extends QuestionBase {
  type: 'vrai-faux';
  prompt: string;
  answer: boolean;
}

/**
 * Texte à trous. Chaque trou est noté "___" (trois tirets bas) dans `prompt`.
 * `answers[i]` est la bonne réponse du i-ème trou.
 * `bank` contient toutes les bonnes réponses + des distracteurs (l'élève tape sur les mots).
 */
export interface FillBlankQuestion extends QuestionBase {
  type: 'trous';
  prompt: string;
  answers: string[];
  bank: string[];
}

export type Question = QcmQuestion | TrueFalseQuestion | FillBlankQuestion;

export interface Flashcard {
  /**
   * Identifiant facultatif, unique dans toute l'application et préfixé par l'id du chapitre
   * (ex. "maths-s-suites-c1"). Sans id, la carte est suivie par son recto (voir lib/srs.ts cardKey) :
   * corriger le recto remet alors son suivi à zéro.
   */
  id?: string;
  front: string;
  back: string;
}

export interface Chapter {
  /** Identifiant unique, ex. "maths-s-suites". */
  id: string;
  title: string;
  /** Une phrase qui dit de quoi parle le chapitre. */
  summary: string;
  /** « L'essentiel en 30 secondes » : 3 à 5 idées à retenir absolument. */
  essentials: string[];
  sections: FicheSection[];
  flashcards: Flashcard[];
  quiz: Question[];
}

export interface Subject {
  /** Identifiant unique, ex. "maths-s". */
  id: string;
  name: string;
  /** Emoji affiché comme icône. */
  icon: string;
  /** Couleur d'accent (hex). */
  color: string;
  /** Examens concernés par cette matière. */
  tracks: TrackId[];
  chapters: Chapter[];
}
