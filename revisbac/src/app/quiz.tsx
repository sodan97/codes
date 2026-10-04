import { Redirect, router, Stack, useLocalSearchParams, useNavigation } from 'expo-router';
import { usePreventRemove } from 'expo-router/react-navigation';
import { useCallback, useEffect, useRef, useState } from 'react';
import { ActivityIndicator, Alert, Animated, Platform, Pressable, ScrollView, StyleSheet, Text, View } from 'react-native';

import { Countdown } from '../components/Countdown';
import { QuestionView } from '../components/QuestionView';
import { Button, Card, Pill, ProgressBar, Screen, styles as ui, useReduceMotion } from '../components/ui';
import { getSubject, getTrack, type QuestionRef } from '../data/catalog';
import type { Chapter } from '../data/types';
import { dayKey } from '../lib/dates';
import { celebrate } from '../lib/feedback';
import {
  effectiveStreak,
  examNote,
  formatNote,
  mention,
  noteTrend,
  todayStats,
  type AnswerResult,
  type QuizMode,
  type Reward,
} from '../lib/gamification';
import { buildQuiz, EXPRESS_SIZE, type QuizSession } from '../lib/quizBuilder';
import { quizContext } from '../lib/selectors';
import { dailyShareText, examShareText, share } from '../lib/share';
import { useProgress } from '../state/progress';
import { colors, nativeDriver, radius } from '../theme';

const MODES: QuizMode[] = ['chapter', 'daily', 'exam', 'review', 'express'];

/** Quitte le quiz ; ouvert par lien profond, il est seul dans la pile : on revient à l'accueil. */
function exitQuiz() {
  if (router.canGoBack()) router.back();
  else router.replace('/');
}

export default function QuizScreen() {
  // `n` : jeton qui force une nouvelle session (« ⚡ Encore 5 ? »).
  const params = useLocalSearchParams<{ mode?: string; id?: string; n?: string }>();
  const mode: QuizMode = MODES.includes(params.mode as QuizMode) ? (params.mode as QuizMode) : 'chapter';
  const { state, loaded, loadError } = useProgress();
  if (!loaded) return <ActivityIndicator color={colors.primary} style={{ flex: 1 }} />;
  // Lecture impossible : l'écran d'accueil propose de réessayer, surtout pas d'onboarding.
  if (loadError) return <Redirect href="/" />;
  // Lien profond ouvert avant l'onboarding.
  if (!state.profile) return <Redirect href="/onboarding" />;
  return <QuizSessionScreen key={`${mode}:${params.id ?? ''}:${params.n ?? ''}`} mode={mode} id={params.id} />;
}

function QuizSessionScreen({ mode, id }: { mode: QuizMode; id?: string }) {
  const { state } = useProgress();
  // Jour du début de la session : un défi commencé avant minuit et fini après reste celui de ce jour-là.
  const [day] = useState(() => dayKey());
  // La session est figée au montage : elle ne doit pas changer pendant que l'élève répond.
  const [session] = useState(() => buildQuiz(mode, id, quizContext(state, day)!));
  // Défi déjà relevé aujourd'hui : simple entraînement, sans XP.
  const [practice] = useState(() => mode === 'daily' && todayStats(state, day).challengeDone);

  if (session.questions.length === 0) {
    return (
      <Screen edges={['bottom']}>
        <Stack.Screen options={{ title: session.title }} />
        <Text style={{ fontSize: 50, textAlign: 'center', marginTop: 40 }}>🎉</Text>
        <Text style={[ui.h2, { textAlign: 'center' }]}>
          {mode === 'review' ? 'Aucune erreur à revoir, bravo !' : 'Pas encore de questions ici.'}
        </Text>
        <Button label="Retour" onPress={exitQuiz} />
      </Screen>
    );
  }
  return (
    <QuizRunner
      mode={mode}
      session={session}
      practice={practice}
      day={day}
      subjectId={mode === 'exam' ? id : undefined}
      chapterId={mode === 'chapter' ? id : undefined}
    />
  );
}

function QuizRunner({
  mode,
  session,
  practice,
  day,
  chapterId,
  subjectId,
}: {
  mode: QuizMode;
  session: QuizSession;
  practice: boolean;
  day: string;
  chapterId?: string;
  subjectId?: string;
}) {
  const { state, finishQuiz } = useProgress();
  const navigation = useNavigation();
  const haptics = state.profile?.haptics ?? true;
  const isExam = mode === 'exam';
  const total = session.questions.length;
  // File des questions restantes (indices dans session.questions) : « Passer » renvoie la question au bout.
  const [order, setOrder] = useState(() => session.questions.map((_, i) => i));
  // Réponse donnée à chaque question (null = pas encore répondu).
  const answers = useRef<(boolean | null)[]>(session.questions.map(() => null));
  const [combo, setCombo] = useState(0);
  // Une réponse a été validée : quitter fait perdre la série en cours.
  const [touched, setTouched] = useState(false);
  const [results, setResults] = useState<AnswerResult[]>([]);
  const [reward, setReward] = useState<Reward | null>(null);
  // Examen blanc : échéance fixée au tap sur « Commencer » (null tant que l'introduction est affichée).
  const [deadline, setDeadline] = useState<number | null>(null);
  const finished = useRef(false);
  const scrollRef = useRef<ScrollView>(null);
  const answered = total - order.length;
  const current: QuestionRef | undefined = session.questions[order[0]];
  const color = session.color ?? current?.subject.color ?? colors.primary;

  const finish = useCallback(() => {
    if (finished.current) return;
    finished.current = true;
    // Questions sans réponse à la fin du temps : elles comptent 0 mais ne vont jamais dans les erreurs.
    const all: AnswerResult[] = session.questions.map((q, i) => {
      const answer = answers.current[i];
      const base = { questionId: q.question.id, subjectId: q.subject.id };
      return answer === null ? { ...base, correct: false, skipped: true } : { ...base, correct: answer };
    });
    setResults(all);
    setReward(finishQuiz(mode, all, { chapterId, subjectId, day }));
  }, [session, finishQuiz, mode, chapterId, subjectId, day]);

  // Alerte « Quitter ? » affichée : une expiration du chrono est mise en attente.
  const quitPrompt = useRef(false);
  const expiredDuringPrompt = useRef(false);
  const onExpire = useCallback(() => {
    if (quitPrompt.current) expiredDuringPrompt.current = true;
    else finish();
  }, [finish]);

  // Sortie protégée (flèche de l'en-tête, retour Android) dès que l'élève a commencé, jusqu'aux résultats.
  const started = isExam ? deadline !== null : touched || answered > 0;
  usePreventRemove(started && !reward, ({ data }) => {
    const leave = () => navigation.dispatch(data.action);
    const [title, message, stay] = isExam
      ? ['Quitter l’examen ?', 'Il ne sera pas noté et tes réponses seront perdues.', 'Continuer l’examen']
      : ['Quitter le quiz ?', 'Tes réponses de cette série ne seront pas comptées.', 'Continuer le quiz'];
    if (Platform.OS === 'web') {
      // Alert n'affiche rien sur le web.
      if (globalThis.confirm?.(`${title}\n${message}`)) leave();
      return;
    }
    quitPrompt.current = true;
    const stayHere = () => {
      quitPrompt.current = false;
      if (expiredDuringPrompt.current) finish();
    };
    Alert.alert(
      title,
      message,
      [
        { text: stay, style: 'cancel', onPress: stayHere },
        {
          text: 'Quitter',
          style: 'destructive',
          onPress: () => {
            quitPrompt.current = false;
            leave();
          },
        },
      ],
      { cancelable: true, onDismiss: stayHere },
    );
  });

  // Correction affichée : on fait défiler pour que « Continuer » soit visible.
  const onChecked = useCallback(() => {
    setTouched(true);
    requestAnimationFrame(() => scrollRef.current?.scrollToEnd({ animated: true }));
  }, []);

  const next = (correct: boolean) => {
    answers.current[order[0]] = correct;
    setCombo(correct ? combo + 1 : 0);
    const rest = order.slice(1);
    if (rest.length === 0) finish();
    else {
      setOrder(rest);
      scrollRef.current?.scrollTo({ y: 0, animated: false });
    }
  };
  const skip = () => {
    setOrder([...order.slice(1), order[0]]);
    scrollRef.current?.scrollTo({ y: 0, animated: false });
  };

  if (reward) {
    return <Results mode={mode} session={session} results={results} reward={reward} color={color} practice={practice} day={day} subjectId={subjectId} />;
  }
  if (isExam && deadline === null) {
    return <ExamIntro session={session} color={color} onStart={() => setDeadline(Date.now() + (session.timeLimit ?? 0) * 1000)} />;
  }

  return (
    <Screen edges={['bottom']} scrollRef={scrollRef}>
      <Stack.Screen options={{ title: session.title }} />
      <View style={[ui.row, { justifyContent: 'space-between', flexWrap: 'wrap' }]}>
        <Text style={styles.counter}>{isExam ? `Répondu ${answered}/${total}` : `Question ${answered + 1}/${total}`}</Text>
        <View style={[ui.row, { flexShrink: 1, flexWrap: 'wrap', justifyContent: 'flex-end' }]}>
          {practice && <Pill label="Entraînement · sans XP" color={colors.muted} />}
          {!isExam && <ComboPill combo={combo} />}
          {isExam && deadline !== null && <Countdown deadline={deadline} onExpire={onExpire} />}
        </View>
      </View>
      <ProgressBar value={answered / total} color={color} height={10} />
      {/* Matière seulement : le chapitre serait un indice, il est révélé dans la correction. */}
      {mode !== 'chapter' && !isExam && current && (
        <Text style={[ui.muted, { color: current.subject.color, fontWeight: '700' }]}>
          {current.subject.icon} {current.subject.name}
        </Text>
      )}
      <Card>
        {current && (
          <QuestionView
            key={current.question.id}
            question={current.question}
            color={color}
            onNext={next}
            deferFeedback={isExam}
            onChecked={onChecked}
            chapterLabel={mode === 'chapter' ? undefined : current.chapter.title}
            haptics={haptics}
          />
        )}
      </Card>
      {isExam && order.length > 1 && <Button label="Passer (j’y reviendrai)" variant="ghost" color={colors.muted} onPress={skip} />}
    </Screen>
  );
}

function ExamIntro({ session, color, onStart }: { session: QuizSession; color: string; onStart: () => void }) {
  const minutes = Math.ceil((session.timeLimit ?? 0) / 60);
  return (
    <Screen edges={['bottom']}>
      <Stack.Screen options={{ title: 'Examen blanc' }} />
      <Card style={{ gap: 10 }}>
        <Text style={ui.h2}>{session.title}</Text>
        <Text style={[ui.body, { fontWeight: '800' }]}>
          {session.questions.length} questions · {minutes} min · note sur 20
        </Text>
        <Text style={ui.body}>Pas de correction pendant l’épreuve : tu verras tout à la fin.</Text>
        <Text style={ui.body}>Tu peux passer une question et y revenir plus tard.</Text>
        <Text style={ui.muted}>Questions de révision RéviBac tirées de tes chapitres (pas un sujet officiel).</Text>
        <Button label="Commencer" color={color} onPress={onStart} />
      </Card>
    </Screen>
  );
}

const COMBO_WORDS: Record<number, string> = { 3: 'Bien joué !', 5: 'En feu !', 10: 'Inarrêtable !' };

/** Pastille de combo : petit « pop » aux paliers 3, 5 et 10 (rien quand le combo est perdu). */
function ComboPill({ combo }: { combo: number }) {
  const reduceMotion = useReduceMotion();
  const [scale] = useState(() => new Animated.Value(1));
  const word = COMBO_WORDS[combo];
  useEffect(() => {
    if (!word || reduceMotion) return;
    const animation = Animated.sequence([
      Animated.timing(scale, { toValue: 1.25, duration: 120, useNativeDriver: nativeDriver }),
      Animated.spring(scale, { toValue: 1, friction: 4, useNativeDriver: nativeDriver }),
    ]);
    animation.start();
    return () => animation.stop();
  }, [combo, word, reduceMotion, scale]);
  if (combo < 2) return null;
  return (
    <Animated.View style={{ transform: [{ scale }] }}>
      <Pill label={`🔥 Combo x${combo}${word ? ` · ${word}` : ''}`} color={colors.goldText} bg={colors.goldSoft} />
    </Animated.View>
  );
}

function Results({
  mode,
  session,
  results,
  reward,
  color,
  practice,
  day,
  subjectId,
}: {
  mode: QuizMode;
  session: QuizSession;
  results: AnswerResult[];
  reward: Reward;
  color: string;
  practice: boolean;
  /** Jour où la session a commencé (celui du défi). */
  day: string;
  subjectId?: string;
}) {
  const { state } = useProgress();
  // Défi commencé un autre jour (avant minuit) : entraînement, comme un défi rejoué.
  const [stale] = useState(() => day !== dayKey());
  const training = practice || stale;
  const [showRight, setShowRight] = useState(false);
  const haptics = state.profile?.haptics ?? true;
  const correct = results.filter((r) => r.correct).length;
  const total = results.length;
  const percent = Math.round((correct / total) * 100);
  const note = examNote(correct, total);
  const refOf = (r: AnswerResult) => session.questions.find((q) => q.question.id === r.questionId)!;
  const wrong = results.filter((r) => !r.correct && !r.skipped).map(refOf);
  const skipped = results.filter((r) => r.skipped).map(refOf);
  const right = results.filter((r) => r.correct).map(refOf);
  const official = state.dailyResults[day];
  const emoji = percent === 100 ? '🏆' : percent >= 80 ? '🌟' : percent >= 50 ? '👍' : '💪';
  const message =
    percent === 100 ? 'Parfait, aucune erreur !' : percent >= 80 ? 'Excellent travail !' : percent >= 50 ? 'Pas mal, continue comme ça !' : 'Relis la fiche et réessaie, tu vas y arriver !';

  useEffect(() => {
    if (reward.levelUp != null) celebrate(haptics);
  }, [reward.levelUp, haptics]);

  const shareDaily = () => {
    if (!state.profile) return;
    void share(
      dailyShareText({
        day,
        trackLabel: getTrack(state.profile.track).label,
        correct,
        total,
        grid: results.map((r) => (r.correct ? '✅' : '❌')).join(''),
        streak: effectiveStreak(state, day),
      }),
    );
  };
  const shareExam = () => {
    const subjectName = getSubject(subjectId ?? '')?.name ?? session.questions[0].subject.name;
    void share(examShareText({ subjectName, note }));
  };

  return (
    <Screen edges={['bottom']}>
      <Stack.Screen options={{ title: 'Résultats' }} />
      {mode === 'daily' && training && (stale || official) && (
        <Card style={{ backgroundColor: colors.bg, borderWidth: 1, borderColor: colors.border }}>
          <Text style={[ui.body, { fontWeight: '700', textAlign: 'center' }]}>
            {stale
              ? 'Entraînement · ce défi date d’hier, il ne compte pas pour aujourd’hui'
              : `Entraînement · ton score officiel du jour reste ${official!.correct}/${official!.total}`}
          </Text>
        </Card>
      )}
      <View style={{ alignItems: 'center', gap: 6, marginTop: 8 }}>
        <Text style={{ fontSize: 64 }}>{emoji}</Text>
        {mode === 'exam' ? (
          <>
            <Text style={[styles.score, { color }]}>{formatNote(note)}/20</Text>
            <Pill label={`Mention : ${mention(note)}`} color={note >= 10 ? colors.primary : colors.red} />
            <Text style={[ui.muted, { textAlign: 'center' }]}>
              Note indicative sur des questions de cours : l’épreuve réelle comporte aussi des exercices rédigés.
            </Text>
            {noteTrend(state.examHistory[subjectId ?? '']) && (
              <Text style={[ui.body, { fontWeight: '700' }]}>Tes dernières notes : {noteTrend(state.examHistory[subjectId ?? ''])}</Text>
            )}
          </>
        ) : (
          <Text style={[styles.score, { color }]}>
            {correct}/{total}
          </Text>
        )}
        <Text style={[ui.h2, { textAlign: 'center' }]}>{message}</Text>
      </View>

      {mode === 'daily' && !training && <Button label="📤 Partager mon score" variant="secondary" color={color} onPress={shareDaily} />}
      {mode === 'exam' && <Button label="📤 Partager ma note" variant="secondary" color={color} onPress={shareExam} />}

      <Card style={{ gap: 6, backgroundColor: colors.goldSoft }}>
        <Text style={styles.xp}>+{reward.xp} XP</Text>
        {reward.levelUp != null && <Text style={styles.line}>🎉 Tu passes au niveau {reward.levelUp} !</Text>}
        {reward.messages.map((m, i) => (
          <Text key={i} style={styles.line}>
            {m}
          </Text>
        ))}
        {reward.newBadges.map((b) => (
          <Text key={b.id} style={[styles.line, { fontWeight: '800' }]}>
            {b.icon} Nouveau badge : {b.name}
          </Text>
        ))}
      </Card>

      {mode === 'exam' && <Diagnostic session={session} results={results} />}

      {wrong.length > 0 && (
        <>
          <Text style={[ui.h2, { marginTop: 4 }]}>📌 À retenir</Text>
          {wrong.map((q) => (
            <AnswerCard key={q.question.id} refQ={q} border={colors.red} />
          ))}
          <Text style={ui.muted}>
            {mode === 'review'
              ? 'Ces questions restent dans ta liste : elles reviendront demain.'
              : 'Ces questions reviendront demain dans « À revoir ».'}
          </Text>
        </>
      )}

      {skipped.length > 0 && (
        <>
          <Text style={[ui.h2, { marginTop: 4 }]}>⏱️ Non traitées ({skipped.length})</Text>
          {skipped.map((q) => (
            <AnswerCard key={q.question.id} refQ={q} border={colors.muted} />
          ))}
          <Text style={ui.muted}>Le temps était écoulé : elles comptent 0, mais ne sont pas ajoutées à tes erreurs.</Text>
        </>
      )}

      {right.length > 0 && (
        <>
          <Button
            label={showRight ? 'Masquer les questions réussies' : `Voir aussi ${right.length > 1 ? `les ${right.length} questions réussies` : 'la question réussie'}`}
            variant="ghost"
            color={colors.primaryDark}
            onPress={() => setShowRight(!showRight)}
          />
          {showRight && right.map((q) => <AnswerCard key={q.question.id} refQ={q} border={colors.primary} />)}
        </>
      )}

      {mode === 'express' && (
        <Button
          label={`⚡ Encore ${EXPRESS_SIZE} ?`}
          color={color}
          onPress={() => router.replace({ pathname: '/quiz', params: { mode: 'express', n: String(Date.now()) } })}
        />
      )}
      <Button label="Terminer" variant={mode === 'express' ? 'secondary' : 'primary'} color={color} onPress={exitQuiz} />
    </Screen>
  );
}

/** Examen blanc : réussite par chapitre, avec un lien vers la fiche des chapitres sous 50 %. */
function Diagnostic({ session, results }: { session: QuizSession; results: AnswerResult[] }) {
  const subject = session.questions[0]?.subject;
  if (!subject) return null;
  const stats = new Map<string, { ok: number; n: number }>();
  for (const r of results) {
    const chapterId = session.questions.find((q) => q.question.id === r.questionId)?.chapter.id;
    if (!chapterId) continue;
    const s = stats.get(chapterId) ?? { ok: 0, n: 0 };
    s.n++;
    if (r.correct) s.ok++;
    stats.set(chapterId, s);
  }
  const rows: { chapter: Chapter; ok: number; n: number }[] = subject.chapters
    .filter((c) => stats.has(c.id))
    .map((c) => ({ chapter: c, ...stats.get(c.id)! }));

  return (
    <>
      <Text style={[ui.h2, { marginTop: 4 }]}>🔎 Diagnostic par chapitre</Text>
      <Card style={{ gap: 14 }}>
        {rows.map(({ chapter, ok, n }) => {
          const weak = ok / n < 0.5;
          return (
            <View key={chapter.id} style={{ gap: 6 }}>
              <Text style={{ fontWeight: '700', color: colors.text }}>
                {chapter.title} · {ok}/{n}
              </Text>
              <ProgressBar value={ok / n} color={weak ? colors.red : ok / n >= 0.8 ? colors.primary : colors.gold} />
              {weak && (
                <Button
                  label="📄 Revoir la fiche"
                  variant="secondary"
                  color={subject.color}
                  style={{ paddingVertical: 8 }}
                  onPress={() => router.push({ pathname: '/fiche/[id]', params: { id: chapter.id } })}
                />
              )}
            </View>
          );
        })}
      </Card>
    </>
  );
}

/** Question corrigée : énoncé, bonne réponse, explication et lien vers la fiche du chapitre. */
function AnswerCard({ refQ, border }: { refQ: QuestionRef; border: string }) {
  return (
    <Card style={{ gap: 6, borderLeftWidth: 4, borderLeftColor: border }}>
      <Text style={{ fontWeight: '700', color: colors.text }}>{refQ.question.prompt}</Text>
      <Text style={{ color: colors.primaryDark, fontWeight: '700' }}>→ {correctAnswerText(refQ)}</Text>
      <Text style={ui.muted}>{refQ.question.explanation}</Text>
      <Pressable
        accessibilityRole="link"
        onPress={() => router.push({ pathname: '/fiche/[id]', params: { id: refQ.chapter.id } })}
        style={({ pressed }) => [styles.link, pressed && { opacity: 0.6 }]}
      >
        <Text style={styles.linkText}>📄 Revoir la fiche : {refQ.chapter.title}</Text>
      </Pressable>
    </Card>
  );
}

function correctAnswerText({ question: q }: QuestionRef): string {
  if (q.type === 'qcm') return q.choices[q.answer];
  if (q.type === 'vrai-faux') return q.answer ? 'Vrai' : 'Faux';
  return q.answers.join(' · ');
}

const styles = StyleSheet.create({
  counter: { fontSize: 15, fontWeight: '800', color: colors.text },
  score: { fontSize: 48, fontWeight: '900' },
  xp: { fontSize: 26, fontWeight: '900', color: colors.primary },
  line: { fontSize: 15, color: colors.text },
  link: { minHeight: 44, justifyContent: 'center', borderRadius: radius.sm },
  linkText: { fontWeight: '700', color: colors.primaryDark },
});
