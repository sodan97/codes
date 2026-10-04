import { useState } from 'react';
import { Pressable, StyleSheet, Text, View } from 'react-native';

import type { FillBlankQuestion, QcmQuestion, Question, TrueFalseQuestion } from '../data/types';
import { shuffle } from '../lib/random';
import { colors, radius } from '../theme';
import { Button } from './ui';

const TYPE_LABEL: Record<Question['type'], string> = {
  qcm: 'QCM · une seule bonne réponse',
  'vrai-faux': 'Vrai ou faux ?',
  trous: 'Complète la phrase',
};

/**
 * Une question interactive. L'élève répond, valide, lit la correction puis passe à la suite.
 * À utiliser avec `key={question.id}` pour repartir d'un état vierge à chaque question.
 */
export function QuestionView({ question, color, onNext }: { question: Question; color: string; onNext: (correct: boolean) => void }) {
  const [checked, setChecked] = useState<boolean | null>(null);

  return (
    <View style={{ gap: 14 }}>
      <Text style={[styles.typeLabel, { color }]}>{TYPE_LABEL[question.type]}</Text>
      {question.type === 'qcm' && <Qcm q={question} color={color} checked={checked} onCheck={setChecked} />}
      {question.type === 'vrai-faux' && <TrueFalse q={question} color={color} checked={checked} onCheck={setChecked} />}
      {question.type === 'trous' && <FillBlank q={question} color={color} checked={checked} onCheck={setChecked} />}

      {checked !== null && (
        <View style={[styles.feedback, { backgroundColor: checked ? colors.primarySoft : colors.redSoft }]}>
          <Text style={[styles.feedbackTitle, { color: checked ? colors.primaryDark : colors.red }]}>
            {checked ? '✅ Bonne réponse !' : '❌ Pas tout à fait…'}
          </Text>
          <Text style={styles.body}>{question.explanation}</Text>
        </View>
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

  const place = (bankIndex: number) => {
    const slot = fills.indexOf(null);
    if (slot === -1) return;
    const next = [...fills];
    next[slot] = bankIndex;
    setFills(next);
  };
  const clear = (slot: number) => {
    if (done) return;
    const next = [...fills];
    next[slot] = null;
    setFills(next);
  };
  const validate = () => onCheck(fills.every((b, i) => b !== null && bank[b] === q.answers[i]));

  return (
    <View style={{ gap: 14 }}>
      <Text style={styles.prompt}>
        {parts.map((part, i) => (
          <Text key={i}>
            {part}
            {i < parts.length - 1 && (
              <Text
                onPress={() => clear(i)}
                style={[
                  styles.blank,
                  { borderColor: color, color },
                  done && { color: fills[i] !== null && bank[fills[i]!] === q.answers[i] ? colors.primary : colors.red },
                ]}
              >
                {fills[i] !== null ? ` ${bank[fills[i]!]} ` : ' ______ '}
              </Text>
            )}
          </Text>
        ))}
      </Text>
      {done && checked === false && (
        <Text style={styles.body}>
          <Text style={{ fontWeight: '800' }}>Bonne(s) réponse(s) : </Text>
          {q.answers.join(' · ')}
        </Text>
      )}
      {!done && (
        <>
          <Text style={styles.hint}>Touche un mot pour le placer, touche un trou pour le vider.</Text>
          <View style={styles.bank}>
            {bank.map((word, i) => {
              const used = fills.includes(i);
              return (
                <Pressable
                  testID="chip"
                  key={word}
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
  blank: { fontWeight: '800', textDecorationLine: 'underline', borderWidth: 0 },
  bank: { flexDirection: 'row', flexWrap: 'wrap', gap: 8 },
  chip: { borderWidth: 2, borderRadius: 999, paddingVertical: 8, paddingHorizontal: 14, backgroundColor: colors.card },
  chipUsed: { borderColor: colors.border, backgroundColor: colors.bg },
  chipText: { fontSize: 15, fontWeight: '600', color: colors.text },
});
