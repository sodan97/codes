import AsyncStorage from '@react-native-async-storage/async-storage';
import { createContext, useCallback, useContext, useEffect, useMemo, useState, type ReactNode } from 'react';
import { Appearance, Platform, useColorScheme } from 'react-native';

import { palettes, type Colors, type Scheme } from '../theme';

/** Choix de l'élève : suivre le téléphone (par défaut), ou forcer le clair ou le sombre. */
export type ThemePreference = 'auto' | 'light' | 'dark';

// Rangée à part de la progression : une sauvegarde abîmée ne touche pas l'apparence, et inversement.
const STORAGE_KEY = 'revisbac/theme';
const PREFERENCES: ThemePreference[] = ['auto', 'light', 'dark'];

interface ThemeContextValue {
  colors: Colors;
  scheme: Scheme;
  preference: ThemePreference;
  setPreference: (preference: ThemePreference) => void;
}

const ThemeContext = createContext<ThemeContextValue | null>(null);

/** Préférence appliquée aux éléments natifs (clavier, alertes, interrupteurs). */
function applyNative(preference: ThemePreference) {
  if (Platform.OS === 'web') return;
  try {
    Appearance.setColorScheme(preference === 'auto' ? 'unspecified' : preference);
  } catch {
    // Pas de forçage possible : seule l'application change de thème.
  }
}

export function ThemeProvider({ children }: { children: ReactNode }) {
  const system = useColorScheme();
  // « Automatique » en attendant la lecture : c'est le cas de la grande majorité des élèves.
  const [preference, setPreferenceState] = useState<ThemePreference>('auto');
  // Préférence relue : rien n'est affiché avant, pour ne pas montrer un instant le mauvais thème.
  const [ready, setReady] = useState(false);

  useEffect(() => {
    let alive = true;
    AsyncStorage.getItem(STORAGE_KEY)
      .then((raw) => {
        if (!alive || !PREFERENCES.includes(raw as ThemePreference)) return;
        setPreferenceState(raw as ThemePreference);
        applyNative(raw as ThemePreference);
      })
      .catch(() => {})
      // Même en cas d'erreur de lecture : l'application démarre en « Automatique ».
      .finally(() => {
        if (alive) setReady(true);
      });
    return () => {
      alive = false;
    };
  }, []);

  const setPreference = useCallback((next: ThemePreference) => {
    setPreferenceState(next);
    applyNative(next);
    AsyncStorage.setItem(STORAGE_KEY, next).catch(() => {});
  }, []);

  // Un forçage natif fait aussi renvoyer la préférence par useColorScheme : la préférence passe d'abord.
  const scheme: Scheme = preference === 'auto' ? (system === 'dark' ? 'dark' : 'light') : preference;
  const value = useMemo(() => ({ colors: palettes[scheme], scheme, preference, setPreference }), [scheme, preference, setPreference]);
  // Lecture très courte : l'écran de démarrage reste affiché (la progression et SplashGate ne sont pas encore montés).
  if (!ready) return null;
  return <ThemeContext.Provider value={value}>{children}</ThemeContext.Provider>;
}

const fallbacks: Record<Scheme, ThemeContextValue> = {
  light: { colors: palettes.light, scheme: 'light', preference: 'auto', setPreference: () => {} },
  dark: { colors: palettes.dark, scheme: 'dark', preference: 'auto', setPreference: () => {} },
};

/** Couleurs du thème courant. Hors ThemeProvider (filet d'erreur racine) : thème du téléphone. */
export function useTheme(): ThemeContextValue {
  return useContext(ThemeContext) ?? fallbacks[Appearance.getColorScheme() === 'dark' ? 'dark' : 'light'];
}

// Une feuille de styles par fabrique et par thème, créée à la première demande puis réutilisée.
const styleCache = new WeakMap<object, Partial<Record<Scheme, unknown>>>();

/**
 * Styles construits à partir du thème. `makeStyles` est déclarée au niveau du module :
 *   const makeStyles = (colors: Colors) => StyleSheet.create({ title: { color: colors.text } });
 *   const styles = useStyles(makeStyles);
 */
export function useStyles<T>(makeStyles: (colors: Colors) => T): T {
  const { colors, scheme } = useTheme();
  let entry = styleCache.get(makeStyles);
  if (!entry) {
    entry = {};
    styleCache.set(makeStyles, entry);
  }
  if (!(scheme in entry)) entry[scheme] = makeStyles(colors);
  return entry[scheme] as T;
}
