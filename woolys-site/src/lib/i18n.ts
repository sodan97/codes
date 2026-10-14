/**
 * Préparation de l'internationalisation.
 *
 * Le site est publié en français. Pour ajouter l'anglais :
 *  1. ajouter "en" dans `locales` ;
 *  2. créer `src/lib/dictionaries/en.ts` sur le modèle de `fr.ts` ;
 *  3. déplacer les routes sous `src/app/[locale]/` et traduire les fichiers de `src/data/`
 *     (chaque fichier exporte des contenus typés, faciles à dupliquer par langue) ;
 *  4. compléter `alternates.languages` dans `src/lib/seo.ts`.
 */
import { fr } from "./dictionaries/fr";

export const locales = ["fr"] as const;
export type Locale = (typeof locales)[number];
export const defaultLocale: Locale = "fr";

const dictionaries = { fr } as const;

export type Dictionary = typeof fr;

export function getDictionary(locale: Locale = defaultLocale): Dictionary {
  return dictionaries[locale];
}

/** Dictionnaire de la langue par défaut, utilisé par les composants. */
export const t = getDictionary();
