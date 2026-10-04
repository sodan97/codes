import { memo, useEffect, useRef, useState } from 'react';
import { AppState, View } from 'react-native';

import { useTheme } from '../state/theme';
import { Pill } from './ui';

const secondsLeft = (deadline: number) => Math.max(0, Math.ceil((deadline - Date.now()) / 1000));

function formatTime(s: number) {
  return `${Math.floor(s / 60)}:${String(s % 60).padStart(2, '0')}`;
}

/**
 * Chrono de l'examen blanc, calculé sur une échéance (Date.now() au moment du départ + durée) :
 * il reste juste en arrière-plan ou écran verrouillé. Seul ce composant se redessine chaque seconde.
 * `onExpire` est appelé une seule fois, à zéro.
 */
export const Countdown = memo(function Countdown({ deadline, onExpire }: { deadline: number; onExpire: () => void }) {
  const { colors } = useTheme();
  const [remaining, setRemaining] = useState(() => secondsLeft(deadline));
  const expired = useRef(false);
  const onExpireRef = useRef(onExpire);
  useEffect(() => {
    onExpireRef.current = onExpire;
  }, [onExpire]);

  useEffect(() => {
    const tick = () => {
      const r = secondsLeft(deadline);
      setRemaining(r);
      if (r === 0 && !expired.current) {
        expired.current = true;
        onExpireRef.current();
      }
    };
    const timer = setInterval(tick, 1000);
    // Retour dans l'application : on recalcule tout de suite.
    const sub = AppState.addEventListener('change', (status) => {
      if (status === 'active') tick();
    });
    return () => {
      clearInterval(timer);
      sub.remove();
    };
  }, [deadline]);

  return (
    <View accessible accessibilityLabel={`Temps restant : ${Math.floor(remaining / 60)} min ${remaining % 60} s`}>
      <Pill label={`⏱️ ${formatTime(remaining)}`} color={remaining < 60 ? colors.red : colors.text} />
    </View>
  );
});
