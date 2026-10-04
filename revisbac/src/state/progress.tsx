import AsyncStorage from '@react-native-async-storage/async-storage';
import { createContext, useCallback, useContext, useEffect, useMemo, useRef, useState, type ReactNode } from 'react';
import { AppState } from 'react-native';

import { dayKey } from '../lib/dates';
import {
  cardReview,
  ensureQuests,
  ficheGain,
  flashcardsGain,
  initialState,
  migrate,
  openChest,
  quizGain,
  rollDay,
  STATE_VERSION,
  withProfile,
  type AnswerResult,
  type Profile,
  type ProgressState,
  type QuizMode,
  type Reward,
} from '../lib/gamification';
import { rescheduleReminders } from '../lib/reminders';
import { clearLastFiche, clearSession } from './session';

// La clé reste « v1 » pour ne rien perdre : la version du schéma est stockée dans l'état lui-même.
const STORAGE_KEY = 'revisbac/progress/v1';
const BACKUP_V1_KEY = 'revisbac/progress/backup-v1';
const corruptKey = (day: string) => `revisbac/progress/corrupt-${day}`;
/** Délai avant de replanifier les rappels : plusieurs commits rapprochés ne donnent qu'une replanification. */
const REMINDER_DEBOUNCE_MS = 2000;

interface ProgressContextValue {
  state: ProgressState;
  loaded: boolean;
  /** La sauvegarde n'a pas pu être lue : rien n'est écrit tant que ce n'est pas résolu. */
  loadError: boolean;
  retryLoad: () => void;
  setProfile: (profile: Profile) => void;
  /** Modifie quelques champs du profil (réglages). Sans effet si l'élève n'a pas encore de profil. */
  updateProfile: (patch: Partial<Profile>) => void;
  markFicheRead: (chapterId: string, subjectId: string) => Reward | null;
  finishFlashcards: (subjectId: string, chapterId: string) => Reward;
  /**
   * Fin d'un quiz (efface aussi la session enregistrée pour la reprise).
   * `opts.day` : jour où la session a commencé (défi du jour commencé avant minuit = entraînement).
   * `opts.maxCombo` : plus longue suite de bonnes réponses (quête « Enchaîne 5 bonnes réponses »).
   */
  finishQuiz: (mode: QuizMode, results: AnswerResult[], opts?: { chapterId?: string; subjectId?: string; day?: string; maxCombo?: number }) => Reward;
  /** Réponse à une flashcard (« Je savais » / « À revoir ») : suivi espacé de la carte, sans XP. `cardKey` : srs.cardKey(chapterId, card). */
  reviewCard: (cardKey: string, knew: boolean) => void;
  /** Ouvre le coffre du jour quand les 3 quêtes sont faites (null sinon). */
  openChest: () => Reward | null;
  /** Remplace toute la progression par une sauvegarde (backup.importCode) ; les rappels sont replanifiés. */
  restore: (state: ProgressState) => void;
  reset: () => void;
}

const ProgressContext = createContext<ProgressContextValue | null>(null);

export function ProgressProvider({ children }: { children: ReactNode }) {
  const [state, setState] = useState<ProgressState>(initialState);
  const [loaded, setLoaded] = useState(false);
  const [loadError, setLoadError] = useState(false);
  // Référence toujours à jour : les actions calculent le nouvel état de façon synchrone
  // pour pouvoir renvoyer la récompense à l'écran appelant.
  const ref = useRef(state);
  // Tant que la lecture a échoué, on n'écrit rien : la sauvegarde existe peut-être encore.
  const writeBlocked = useRef(true);
  // Sauvegarde v1 à recopier avant la première écriture au nouveau format.
  const pendingBackup = useRef<string | null>(null);
  // Les écritures passent l'une après l'autre, dans l'ordre des commits.
  const writes = useRef<Promise<void>>(Promise.resolve());
  const reminderTimer = useRef<ReturnType<typeof setTimeout> | null>(null);

  // Rappels du jour et des 6 suivants, recalculés d'après l'état le plus récent.
  const scheduleReminders = useCallback(() => {
    if (reminderTimer.current) clearTimeout(reminderTimer.current);
    reminderTimer.current = setTimeout(() => {
      reminderTimer.current = null;
      if (!writeBlocked.current) void rescheduleReminders(ref.current);
    }, REMINDER_DEBOUNCE_MS);
  }, []);

  useEffect(
    () => () => {
      if (reminderTimer.current) clearTimeout(reminderTimer.current);
    },
    [],
  );

  const fail = useCallback((e: unknown) => {
    console.warn('Lecture de la progression impossible', e);
    writeBlocked.current = true;
    setLoadError(true);
    setLoaded(true);
  }, []);

  const hydrate = useCallback(
    async (raw: string | null) => {
      let next = initialState();
      if (raw) {
        try {
          const parsed: unknown = JSON.parse(raw);
          next = migrate(parsed);
          const version = (parsed as { version?: unknown }).version;
          if (typeof version !== 'number' || version < STATE_VERSION) pendingBackup.current = raw;
        } catch (e) {
          // Sauvegarde abîmée : on la met de côté avant de repartir d'un état vierge.
          console.warn('Progression illisible, copie mise de côté', e);
          try {
            await AsyncStorage.setItem(corruptKey(dayKey()), raw);
          } catch (err) {
            // Sans copie de côté, on n'écrase surtout pas l'original.
            fail(err);
            return;
          }
        }
      }
      ref.current = next;
      setState(next);
      writeBlocked.current = false;
      setLoadError(false);
      setLoaded(true);
      scheduleReminders();
    },
    [fail, scheduleReminders],
  );

  const load = useCallback(() => {
    AsyncStorage.getItem(STORAGE_KEY).then(hydrate, fail);
  }, [hydrate, fail]);

  useEffect(() => {
    load();
  }, [load]);

  // Retour dans l'application après minuit : défi, objectif et série du nouveau jour.
  useEffect(() => {
    const sub = AppState.addEventListener('change', (status) => {
      if (status !== 'active') return;
      const rolled = rollDay(ref.current);
      if (rolled !== ref.current) {
        ref.current = rolled;
        setState(rolled);
      }
      scheduleReminders();
    });
    return () => sub.remove();
  }, [scheduleReminders]);

  const commit = useCallback((next: ProgressState) => {
    if (writeBlocked.current) return;
    ref.current = next;
    setState(next);
    const json = JSON.stringify(next);
    writes.current = writes.current
      .then(async () => {
        if (pendingBackup.current) {
          await AsyncStorage.setItem(BACKUP_V1_KEY, pendingBackup.current);
          pendingBackup.current = null;
        }
        await AsyncStorage.setItem(STORAGE_KEY, json);
      })
      .catch((e) => console.warn('Sauvegarde de la progression impossible', e));
    scheduleReminders();
  }, [scheduleReminders]);

  const value = useMemo<ProgressContextValue>(
    () => ({
      state,
      loaded,
      loadError,
      retryLoad: () => {
        load();
      },
      setProfile: (profile) => commit(withProfile(ref.current, profile, dayKey())),
      updateProfile: (patch) => {
        const profile = ref.current.profile;
        if (profile) commit(withProfile(ref.current, { ...profile, ...patch }, dayKey()));
      },
      markFicheRead: (chapterId, subjectId) => {
        const gain = ficheGain(ref.current, chapterId, subjectId, dayKey());
        if (!gain) return null;
        commit(gain.state);
        return gain.reward;
      },
      finishFlashcards: (subjectId, chapterId) => {
        const { state: next, reward } = flashcardsGain(ref.current, subjectId, chapterId);
        commit(next);
        return reward;
      },
      finishQuiz: (mode, results, opts) => {
        const { state: next, reward } = quizGain(ref.current, mode, results, opts);
        commit(next);
        clearSession();
        return reward;
      },
      reviewCard: (cardKey, knew) => {
        const next = cardReview(ref.current, cardKey, knew, dayKey());
        if (next !== ref.current) commit(next);
      },
      openChest: () => {
        const gain = openChest(ref.current, dayKey());
        if (!gain) return null;
        commit(gain.state);
        return gain.reward;
      },
      restore: (next) => {
        const today = dayKey();
        commit(ensureQuests(rollDay(next, today), today));
        clearSession();
        clearLastFiche();
      },
      reset: () => {
        commit(initialState());
        clearSession();
        clearLastFiche();
      },
    }),
    [state, loaded, loadError, commit, load],
  );

  return <ProgressContext.Provider value={value}>{children}</ProgressContext.Provider>;
}

export function useProgress(): ProgressContextValue {
  const ctx = useContext(ProgressContext);
  if (!ctx) throw new Error('useProgress doit être utilisé dans <ProgressProvider>');
  return ctx;
}
