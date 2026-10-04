import { router, Stack, useLocalSearchParams } from 'expo-router';
import { useCallback, useEffect, useRef, useState } from 'react';
import { StyleSheet, Text, View } from 'react-native';

import { QuestionView } from '../components/QuestionView';
import { Button, Card, Pill, ProgressBar, Screen, styles as ui } from '../components/ui';
import type { QuestionRef } from '../data/catalog';
import { dayKey } from '../lib/dates';
import { mention, type AnswerResult, type QuizMode, type Reward } from '../lib/gamification';
import { buildQuiz, type QuizSession } from '../lib/quizBuilder';
import { useProgress } from '../state/progress';
import { colors } from '../theme';

const MODES: QuizMode[] = ['chapter', 'daily', 'exam', 'review'];

export default function QuizScreen() {
  const params = useLocalSearchParams<{ mode?: string; id?: string }>();
  const mode: QuizMode = MODES.includes(params.mode as QuizMode) ? (params.mode as QuizMode) : 'chapter';
  const { state } = useProgress();
  // La session est figée au montage : elle ne doit pas changer pendant que l'élève répond.
  const [session] = useState(() => buildQuiz(mode, params.id, state.profile!.track, state.mistakes, dayKey()));

  if (session.questions.length === 0) {
    return (
      <Screen edges={['bottom']}>
        <Stack.Screen options={{ title: session.title }} />
        <Text style={{ fontSize: 50, textAlign: 'center', marginTop: 40 }}>🎉</Text>
        <Text style={[ui.h2, { textAlign: 'center' }]}>
          {mode === 'review' ? 'Aucune erreur à revoir, bravo !' : 'Pas encore de questions ici.'}
        </Text>
        <Button label="Retour" onPress={() => router.back()} />
      </Screen>
    );
  }
  return <QuizRunner mode={mode} session={session} subjectId={mode === 'exam' ? params.id : undefined} chapterId={mode === 'chapter' ? params.id : undefined} />;
}

function QuizRunner({ mode, session, chapterId, subjectId }: { mode: QuizMode; session: QuizSession; chapterId?: string; subjectId?: string }) {
  const { finishQuiz } = useProgress();
  const [index, setIndex] = useState(0);
  const [results, setResults] = useState<AnswerResult[]>([]);
  const [combo, setCombo] = useState(0);
  const [reward, setReward] = useState<Reward | null>(null);
  const [remaining, setRemaining] = useState(session.timeLimit ?? 0);
  const finished = useRef(false);
  const total = session.questions.length;
  const current: QuestionRef | undefined = session.questions[index];
  const color = session.color ?? current?.subject.color ?? colors.primary;

  const finish = useCallback(
    (all: AnswerResult[]) => {
      if (finished.current) return;
      finished.current = true;
      setResults(all);
      setReward(finishQuiz(mode, all, { chapterId, subjectId }));
    },
    [finishQuiz, mode, chapterId, subjectId],
  );

  // Chronomètre de l'examen blanc : à zéro, les questions restantes comptent comme fausses.
  useEffect(() => {
    if (!session.timeLimit || reward) return;
    const timer = setInterval(() => setRemaining((r) => Math.max(0, r - 1)), 1000);
    return () => clearInterval(timer);
  }, [session, reward]);
  useEffect(() => {
    if (!session.timeLimit || remaining > 0) return;
    const missing = session.questions.slice(results.length).map((q) => ({ questionId: q.question.id, subjectId: q.subject.id, correct: false }));
    finish([...results, ...missing]);
  }, [remaining, session, results, finish]);

  const next = (correct: boolean) => {
    const all = [...results, { questionId: current!.question.id, subjectId: current!.subject.id, correct }];
    setCombo(correct ? combo + 1 : 0);
    if (index + 1 >= total) finish(all);
    else {
      setResults(all);
      setIndex(index + 1);
    }
  };

  if (reward) return <Results mode={mode} session={session} results={results} reward={reward} color={color} />;

  return (
    <Screen edges={['bottom']}>
      <Stack.Screen options={{ title: session.title }} />
      <View style={[ui.row, { justifyContent: 'space-between' }]}>
        <Text style={styles.counter}>
          Question {index + 1}/{total}
        </Text>
        <View style={ui.row}>
          {combo >= 2 && <Pill label={`🔥 Combo x${combo}`} color={colors.gold} bg={colors.goldSoft} />}
          {session.timeLimit ? <Pill label={`⏱️ ${formatTime(remaining)}`} color={remaining < 60 ? colors.red : colors.text} /> : null}
        </View>
      </View>
      <ProgressBar value={index / total} color={color} height={10} />
      {mode !== 'chapter' && current && (
        <Text style={[ui.muted, { color: current.subject.color, fontWeight: '700' }]}>
          {current.subject.icon} {current.subject.name} · {current.chapter.title}
        </Text>
      )}
      <Card>{current && <QuestionView key={current.question.id + index} question={current.question} color={color} onNext={next} />}</Card>
    </Screen>
  );
}

function Results({ mode, session, results, reward, color }: { mode: QuizMode; session: QuizSession; results: AnswerResult[]; reward: Reward; color: string }) {
  const correct = results.filter((r) => r.correct).length;
  const total = results.length;
  const percent = Math.round((correct / total) * 100);
  const note = Math.round((correct / total) * 20 * 2) / 2;
  const wrong = results.filter((r) => !r.correct).map((r) => session.questions.find((q) => q.question.id === r.questionId)!);
  const emoji = percent === 100 ? '🏆' : percent >= 80 ? '🌟' : percent >= 50 ? '👍' : '💪';
  const message =
    percent === 100 ? 'Parfait, aucune erreur !' : percent >= 80 ? 'Excellent travail !' : percent >= 50 ? 'Pas mal, continue comme ça !' : 'Relis la fiche et réessaie, tu vas y arriver !';

  return (
    <Screen edges={['bottom']}>
      <Stack.Screen options={{ title: 'Résultats' }} />
      <View style={{ alignItems: 'center', gap: 6, marginTop: 8 }}>
        <Text style={{ fontSize: 64 }}>{emoji}</Text>
        {mode === 'exam' ? (
          <>
            <Text style={[styles.score, { color }]}>{note}/20</Text>
            <Pill label={`Mention : ${mention(note)}`} color={note >= 10 ? colors.primary : colors.red} />
          </>
        ) : (
          <Text style={[styles.score, { color }]}>
            {correct}/{total}
          </Text>
        )}
        <Text style={[ui.h2, { textAlign: 'center' }]}>{message}</Text>
      </View>

      <Card style={{ gap: 6, backgroundColor: colors.goldSoft }}>
        <Text style={styles.xp}>+{reward.xp} XP</Text>
        {reward.levelUp != null && <Text style={styles.line}>⬆️ Tu passes au niveau {reward.levelUp} !</Text>}
        {reward.messages.map((m) => (
          <Text key={m} style={styles.line}>
            {m}
          </Text>
        ))}
        {reward.newBadges.map((b) => (
          <Text key={b.id} style={[styles.line, { fontWeight: '800' }]}>
            {b.icon} Nouveau badge : {b.name}
          </Text>
        ))}
      </Card>

      {wrong.length > 0 && (
        <>
          <Text style={[ui.h2, { marginTop: 4 }]}>📌 À retenir</Text>
          {wrong.map((q, i) => (
            <Card key={q.question.id + i} style={{ gap: 6, borderLeftWidth: 4, borderLeftColor: colors.red }}>
              <Text style={{ fontWeight: '700', color: colors.text }}>{q.question.prompt}</Text>
              <Text style={{ color: colors.primaryDark, fontWeight: '700' }}>→ {correctAnswerText(q)}</Text>
              <Text style={ui.muted}>{q.question.explanation}</Text>
            </Card>
          ))}
          <Text style={ui.muted}>Ces questions sont ajoutées à « Revoir mes erreurs » (onglet Défis).</Text>
        </>
      )}

      <Button label="Terminer" color={color} onPress={() => router.back()} />
    </Screen>
  );
}

function correctAnswerText({ question: q }: QuestionRef): string {
  if (q.type === 'qcm') return q.choices[q.answer];
  if (q.type === 'vrai-faux') return q.answer ? 'Vrai' : 'Faux';
  return q.answers.join(' · ');
}

function formatTime(s: number) {
  return `${Math.floor(s / 60)}:${String(s % 60).padStart(2, '0')}`;
}

const styles = StyleSheet.create({
  counter: { fontSize: 15, fontWeight: '800', color: colors.text },
  score: { fontSize: 48, fontWeight: '900' },
  xp: { fontSize: 26, fontWeight: '900', color: colors.primary },
  line: { fontSize: 15, color: colors.text },
});
