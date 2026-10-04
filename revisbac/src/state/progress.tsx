import AsyncStorage from '@react-native-async-storage/async-storage';
import { createContext, useCallback, useContext, useEffect, useMemo, useRef, useState, type ReactNode } from 'react';

import {
  XP,
  applyGain,
  initialState,
  quizGain,
  rollDay,
  type AnswerResult,
  type Profile,
  type ProgressState,
  type QuizMode,
  type Reward,
} from '../lib/gamification';

const STORAGE_KEY = 'revisbac/progress/v1';

interface ProgressContextValue {
  state: ProgressState;
  loaded: boolean;
  setProfile: (profile: Profile) => void;
  markFicheRead: (chapterId: string, subjectId: string) => Reward | null;
  finishFlashcards: (subjectId: string) => Reward;
  finishQuiz: (mode: QuizMode, results: AnswerResult[], opts?: { chapterId?: string; subjectId?: string }) => Reward;
  reset: () => void;
}

const ProgressContext = createContext<ProgressContextValue | null>(null);

export function ProgressProvider({ children }: { children: ReactNode }) {
  const [state, setState] = useState<ProgressState>(initialState);
  const [loaded, setLoaded] = useState(false);
  // Référence toujours à jour : les actions calculent le nouvel état de façon synchrone
  // pour pouvoir renvoyer la récompense à l'écran appelant.
  const ref = useRef(state);

  useEffect(() => {
    AsyncStorage.getItem(STORAGE_KEY)
      .then((raw) => {
        if (raw) {
          const saved = JSON.parse(raw) as ProgressState;
          const next = rollDay({ ...initialState(), ...saved });
          ref.current = next;
          setState(next);
        }
      })
      .catch(() => {})
      .finally(() => setLoaded(true));
  }, []);

  const commit = useCallback((next: ProgressState) => {
    ref.current = next;
    setState(next);
    AsyncStorage.setItem(STORAGE_KEY, JSON.stringify(next)).catch(() => {});
  }, []);

  const value = useMemo<ProgressContextValue>(
    () => ({
      state,
      loaded,
      setProfile: (profile) => commit({ ...ref.current, profile }),
      markFicheRead: (chapterId, subjectId) => {
        if (ref.current.fichesRead[chapterId]) return null;
        const { state: next, reward } = applyGain(ref.current, XP.ficheRead, [`📄 Fiche lue : +${XP.ficheRead} XP`], (s) => {
          s.fichesRead[chapterId] = s.today.day;
          if (!s.subjectsTouched.includes(subjectId)) s.subjectsTouched.push(subjectId);
        });
        commit(next);
        return reward;
      },
      finishFlashcards: (subjectId) => {
        const { state: next, reward } = applyGain(ref.current, XP.flashcards, [`🃏 Flashcards terminées : +${XP.flashcards} XP`], (s) => {
          if (!s.subjectsTouched.includes(subjectId)) s.subjectsTouched.push(subjectId);
        });
        commit(next);
        return reward;
      },
      finishQuiz: (mode, results, opts) => {
        const { state: next, reward } = quizGain(ref.current, mode, results, opts);
        commit(next);
        return reward;
      },
      reset: () => commit({ ...initialState() }),
    }),
    [state, loaded, commit],
  );

  return <ProgressContext.Provider value={value}>{children}</ProgressContext.Provider>;
}

export function useProgress(): ProgressContextValue {
  const ctx = useContext(ProgressContext);
  if (!ctx) throw new Error('useProgress doit être utilisé dans <ProgressProvider>');
  return ctx;
}
