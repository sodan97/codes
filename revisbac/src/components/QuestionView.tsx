import { useEffect, useState } from 'react';
import { Animated, Pressable, StyleSheet, Text, View } from 'react-native';

import type { FillBlankQuestion, QcmQuestion, Question, TrueFalseQuestion } from '../data/types';
import { bad, good } from '../lib/feedback';
import { shuffle } from '../lib/random';
import { colors, nativeDriver, radius } from '../theme';
import { Button, useReduceMotion } from './ui';

const TYPE_LABEL: Record<Question['type'], string> = {
  qcm: 'QCM · une seule bonne réponse',
  'vrai-faux': 'Vrai ou faux ?',
  trous: 'Complète la phrase',
};

/**
 * Une question interactive. L'élève répond, valide, lit la correction puis passe à la suite.
 * À utiliser avec `key={question.id}` pour repartir d'un état vierge à chaque question.
 * • `deferFeedback` (examen blanc) : « Valider » passe directement à la suite, sans correction.
 * • `onChecked` : appelé quand la correction s'affiche (pour faire défiler jusqu'à « Continuer »).
 * • `chapterLabel` : chapitre révélé dans la correction seulement (pas d'indice avant la réponse).
 * • `haptics` : réglage « Vibrations » du profil.
 */
export function QuestionView({
  question,
  color,
  onNext,
  deferFeedback = false,
  onChecked,
  chapterLabel,
  haptics = true,
}: {
  question: Question;
  color: string;
  onNext: (correct: boolean) => void;
  deferFeedback?: boolean;
  onChecked?: () => void;
  chapterLabel?: string;
  haptics?: boolean;
}) {
  const [checked, setChecked] = useState<boolean | null>(null);
  const reduceMotion = useReduceMotion();
  const [scale] = useState(() => new Animated.Value(1));
  const [shift] = useState(() => new Animated.Value(0));

  const check = (correct: boolean) => {
    if (deferFeedback) {
      onNext(correct);
      return;
    }
    setChecked(correct);
    if (correct) good(haptics);
    else bad(haptics);
    if (reduceMotion) return;
    if (correct) {
      // Petit rebond 1 → 1,04 → 1 du bloc de correction.
      Animated.sequence([
        Animated.timing(scale, { toValue: 1.04, duration: 120, useNativeDriver: nativeDriver }),
        Animated.spring(scale, { toValue: 1, friction: 4, useNativeDriver: nativeDriver }),
      ]).start();
    } else {
      // Léger tremblement horizontal (±6 px sur 250 ms).
      Animated.sequence(
        [6, -6, 6, -6, 0].map((toValue) => Animated.timing(shift, { toValue, duration: 50, useNativeDriver: nativeDriver })),
      ).start();
    }
  };

  // Après l'affichage de la correction (le bloc est alors mis en page).
  useEffect(() => {
    if (checked !== null) onChecked?.();
    // onChecked n'est appelé qu'une fois, au moment de la correction.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [checked]);

  return (
    <View style={{ gap: 14 }}>
      <Text style={[styles.typeLabel, { color }]}>{TYPE_LABEL[question.type]}</Text>
      {question.type === 'qcm' && <Qcm q={question} color={color} checked={checked} onCheck={check} />}
      {question.type === 'vrai-faux' && <TrueFalse q={question} color={color} checked={checked} onCheck={check} />}
      {question.type === 'trous' && <FillBlank q={question} color={color} checked={checked} onCheck={check} />}

      {checked !== null && (
        <Animated.View
          accessibilityLiveRegion="polite"
          style={[
            styles.feedback,
            { backgroundColor: checked ? colors.primarySoft : colors.redSoft, transform: [{ scale }, { translateX: shift }] },
          ]}
        >
          <Text style={[styles.feedbackTitle, { color: checked ? colors.primaryDark : colors.red }]}>
            {checked ? '✅ Bonne réponse !' : '❌ Pas tout à fait…'}
          </Text>
          <Text style={styles.body}>{question.explanation}</Text>
          {chapterLabel ? <Text style={styles.chapter}>Chapitre : {chapterLabel}</Text> : null}
        </Animated.View>
      )}
      {checked !== null && <Button testID="continue" label="Continuer" onPress={() => onNext(checked)} color={checked ? colors.primary : colors.text} />}
    </View>
  );
}

interface Props<Q> {
  q: Q;
  color: string;
  checked: boolean | null;
  onCheck: (correct: boolean) => void;
}

function Qcm({ q, color, checked, onCheck }: Props<QcmQuestion>) {
  // Les choix sont mélangés pour qu'on ne retienne pas la position de la réponse.
  const [order] = useState(() => shuffle(q.choices.map((_, i) => i)));
  const [selected, setSelected] = useState<number | null>(null);
  const done = checked !== null;

  return (
    <View style={{ gap: 10 }}>
      <Text style={styles.prompt}>{q.prompt}</Text>
      {order.map((i) => {
        const isAnswer = i === q.answer;
        const isSelected = i === selected;
        const state = done ? (isAnswer ? 'good' : isSelected ? 'bad' : 'idle') : isSelected ? 'selected' : 'idle';
        return (
          <Choice key={i} label={q.choices[i]} state={state} color={color} disabled={done} onPress={() => setSelected(i)} />
        );
      })}
      {!done && <Button testID="validate" label="Valider" disabled={selected === null} onPress={() => onCheck(selected === q.answer)} color={color} />}
    </View>
  );
}

function TrueFalse({ q, color, checked, onCheck }: Props<TrueFalseQuestion>) {
  const [selected, setSelected] = useState<boolean | null>(null);
  const done = checked !== null;
  const stateOf = (v: boolean) => (done ? (v === q.answer ? 'good' : v === selected ? 'bad' : 'idle') : v === selected ? 'selected' : 'idle');

  return (
    <View style={{ gap: 10 }}>
      <Text style={styles.prompt}>{q.prompt}</Text>
      <View style={{ flexDirection: 'row', gap: 10 }}>
        <View style={{ flex: 1 }}>
          <Choice label="👍 Vrai" state={stateOf(true)} color={color} disabled={done} onPress={() => setSelected(true)} center />
        </View>
        <View style={{ flex: 1 }}>
          <Choice label="👎 Faux" state={stateOf(false)} color={color} disabled={done} onPress={() => setSelected(false)} center />
        </View>
      </View>
      {!done && <Button testID="validate" label="Valider" disabled={selected === null} onPress={() => onCheck(selected === q.answer)} color={color} />}
    </View>
  );
}

function FillBlank({ q, color, checked, onCheck }: Props<FillBlankQuestion>) {
  const [bank] = useState(() => shuffle(q.bank));
  // fills[i] = index dans `bank` du mot placé dans le trou i
  const [fills, setFills] = useState<(number | null)[]>(() => q.answers.map(() => null));
  const done = checked !== null;
  const parts = q.prompt.split('___');
  // Prochain trou rempli par un tap sur un mot (surligné).
  const nextSlot = fills.indexOf(null);
  const isRight = (slot: number) => fills[slot] !== null && bank[fills[slot]!] === q.answers[slot];

  const place = (bankIndex: number) => {
    if (nextSlot === -1) return;
    const next = [...fills];
    next[nextSlot] = bankIndex;
    setFills(next);
  };
  const clear = (slot: number) => {
    if (done) return;
    const next = [...fills];
    next[slot] = null;
    setFills(next);
  };
  const validate = () => onCheck(fills.every((_, i) => isRight(i)));

  return (
    <View style={{ gap: 14 }}>
      <Text style={styles.prompt}>
        {parts.map((part, i) => (
          <Text key={i}>
            {part}
            {i < parts.length - 1 && (
              <Text
                style={[
                  styles.blank,
                  { color },
                  !done && i === nextSlot && { backgroundColor: color + '1A' },
                  done && { color: isRight(i) ? colors.primary : colors.red },
                ]}
              >
                {fills[i] !== null ? ` ${bank[fills[i]!]} ` : parts.length > 2 ? ` (${i + 1}) ______ ` : ' ______ '}
              </Text>
            )}
          </Text>
        ))}
      </Text>
      {/* Les trous en grandes cases, faciles à toucher et lisibles par TalkBack. */}
      <View style={styles.slots}>
        {fills.map((b, i) => {
          const active = !done && i === nextSlot;
          const border = done ? (isRight(i) ? colors.primary : colors.red) : active ? color : colors.border;
          return (
            <Pressable
              key={i}
              testID="slot"
              accessibilityRole="button"
              accessibilityLabel={`Trou ${i + 1} : ${b !== null ? bank[b] : 'vide'}`}
              accessibilityHint={!done && b !== null ? 'Touche pour vider ce trou' : undefined}
              accessibilityState={{ disabled: done || b === null, selected: active }}
              disabled={done || b === null}
              onPress={() => clear(i)}
              style={[styles.slot, { borderColor: border }, active && { backgroundColor: color + '14' }]}
            >
              <Text style={styles.slotText}>
                <Text style={{ color: colors.muted }}>Trou {i + 1} : </Text>
                {b !== null ? bank[b] : '…'}
              </Text>
            </Pressable>
          );
        })}
      </View>
      {done && checked === false && (
        <Text style={styles.body}>
          <Text style={{ fontWeight: '800' }}>Bonne(s) réponse(s) : </Text>
          {q.answers.join(' · ')}
        </Text>
      )}
      {!done && (
        <>
          <Text style={styles.hint}>Touche un mot pour le placer dans le trou surligné, touche une case pour la vider.</Text>
          <View style={styles.bank}>
            {bank.map((word, i) => {
              const used = fills.includes(i);
              return (
                <Pressable
                  testID="chip"
                  key={i}
                  accessibilityRole="button"
                  accessibilityState={{ disabled: used }}
                  disabled={used}
                  onPress={() => place(i)}
                  style={({ pressed }) => [styles.chip, { borderColor: color }, used && styles.chipUsed, pressed && { opacity: 0.7 }]}
                >
                  <Text style={[styles.chipText, used && { color: colors.border }]}>{word}</Text>
                </Pressable>
              );
            })}
          </View>
          <Button testID="validate" label="Valider" disabled={fills.includes(null)} onPress={validate} color={color} />
        </>
      )}
    </View>
  );
}

type ChoiceState = 'idle' | 'selected' | 'good' | 'bad';

function Choice({
  label,
  state,
  color,
  disabled,
  onPress,
  center,
}: {
  label: string;
  state: ChoiceState;
  color: string;
  disabled: boolean;
  onPress: () => void;
  center?: boolean;
}) {
  const border = state === 'good' ? colors.primary : state === 'bad' ? colors.red : state === 'selected' ? color : colors.border;
  const bg = state === 'good' ? colors.primarySoft : state === 'bad' ? colors.redSoft : state === 'selected' ? color + '14' : colors.card;
  return (
    <Pressable
      testID="choice"
      accessibilityRole="button"
      accessibilityState={{ selected: state === 'selected' }}
      disabled={disabled}
      onPress={onPress}
      style={({ pressed }) => [styles.choice, { borderColor: border, backgroundColor: bg }, center && { alignItems: 'center' }, pressed && { opacity: 0.8 }]}
    >
      <Text style={styles.choiceText}>
        {state === 'good' ? '✓ ' : state === 'bad' ? '✗ ' : ''}
        {label}
      </Text>
    </Pressable>
  );
}

const styles = StyleSheet.create({
  typeLabel: { fontSize: 12, fontWeight: '800', textTransform: 'uppercase', letterSpacing: 0.6 },
  prompt: { fontSize: 19, lineHeight: 30, fontWeight: '600', color: colors.text },
  body: { fontSize: 15, lineHeight: 22, color: colors.text },
  hint: { fontSize: 13, color: colors.muted },
  choice: { borderWidth: 2, borderRadius: radius.md, paddingVertical: 14, paddingHorizontal: 14 },
  choiceText: { fontSize: 16, color: colors.text, fontWeight: '500' },
  feedback: { borderRadius: radius.md, padding: 14, gap: 6 },
  feedbackTitle: { fontSize: 16, fontWeight: '800' },
  chapter: { fontSize: 13, fontWeight: '700', color: colors.muted },
  blank: { fontWeight: '800', textDecorationLine: 'underline', borderWidth: 0 },
  slots: { gap: 8 },
  slot: { minHeight: 44, borderWidth: 2, borderRadius: radius.sm, paddingHorizontal: 12, paddingVertical: 8, justifyContent: 'center', backgroundColor: colors.card },
  slotText: { fontSize: 15, fontWeight: '600', color: colors.text },
  bank: { flexDirection: 'row', flexWrap: 'wrap', gap: 8 },
  chip: { borderWidth: 2, borderRadius: 999, paddingVertical: 8, paddingHorizontal: 14, backgroundColor: colors.card },
  chipUsed: { borderColor: colors.border, backgroundColor: colors.bg },
  chipText: { fontSize: 15, fontWeight: '600', color: colors.text },
});
