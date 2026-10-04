// Session de quiz en cours et dernière fiche ouverte, gardées sur le téléphone pour pouvoir reprendre
// après une interruption. Les règles (validation, expiration, notation automatique) sont dans lib/resume.ts.
import AsyncStorage from '@react-native-async-storage/async-storage';
import { useFocusEffect } from 'expo-router';
import { useCallback, useState } from 'react';

import { parseLastFiche, parseSavedSession, type LastFiche, type SavedSession } from '../lib/resume';

const SESSION_KEY = 'revisbac/session/v1';
const LAST_FICHE_KEY = 'revisbac/last-fiche/v1';

// Les écritures et effacements passent l'un après l'autre, dans l'ordre des appels :
// un effacement juste après la dernière réponse n'est jamais écrasé par elle.
let queue: Promise<void> = Promise.resolve();

function enqueue(task: () => Promise<void>) {
  queue = queue.then(task).catch((e) => console.warn('Session de quiz : écriture impossible', e));
}

async function read<T>(key: string, parse: (raw: unknown) => T | null): Promise<T | null> {
  await queue;
  try {
    const raw = await AsyncStorage.getItem(key);
    return raw ? parse(JSON.parse(raw)) : null;
  } catch (e) {
    console.warn('Session de quiz : lecture impossible', e);
    return null;
  }
}

/** Enregistre la session en cours (sans attendre la fin de l'écriture). À appeler après chaque réponse validée. */
export function saveSession(session: SavedSession): void {
  const json = JSON.stringify(session);
  enqueue(() => AsyncStorage.setItem(SESSION_KEY, json));
}

/** Session enregistrée, vérifiée (null s'il n'y en a pas ou si elle est abîmée). */
export function loadSession(): Promise<SavedSession | null> {
  return read(SESSION_KEY, parseSavedSession);
}

/** Efface la session (quiz terminé, « Quitter » confirmé, « Recommencer », progression remplacée). */
export function clearSession(): void {
  enqueue(() => AsyncStorage.removeItem(SESSION_KEY));
}

/** Note la fiche ouverte (« Continuer ta fiche » sur l'accueil). */
export function saveLastFiche(chapterId: string): void {
  const value: LastFiche = { chapterId, openedAt: Date.now() };
  const json = JSON.stringify(value);
  enqueue(() => AsyncStorage.setItem(LAST_FICHE_KEY, json));
}

export function loadLastFiche(): Promise<LastFiche | null> {
  return read(LAST_FICHE_KEY, parseLastFiche);
}

export function clearLastFiche(): void {
  enqueue(() => AsyncStorage.removeItem(LAST_FICHE_KEY));
}

/**
 * Session enregistrée et dernière fiche ouverte, relues à chaque fois que l'écran reprend le focus.
 * `reload` relit après un changement (session notée ou effacée).
 */
export function useSavedSession(): { session: SavedSession | null; lastFiche: LastFiche | null; reload: () => void } {
  const [session, setSession] = useState<SavedSession | null>(null);
  const [lastFiche, setLastFiche] = useState<LastFiche | null>(null);
  const [token, setToken] = useState(0);
  useFocusEffect(
    useCallback(() => {
      let active = true;
      void Promise.all([loadSession(), loadLastFiche()]).then(([s, f]) => {
        if (!active) return;
        setSession(s);
        setLastFiche(f);
      });
      return () => {
        active = false;
      };
      // `token` force une nouvelle lecture (reload).
      // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [token]),
  );
  return { session, lastFiche, reload: useCallback(() => setToken((t) => t + 1), []) };
}
