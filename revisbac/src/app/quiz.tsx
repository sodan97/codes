import { Redirect, router, Stack, useLocalSearchParams, useNavigation } from 'expo-router';
import { usePreventRemove } from 'expo-router/react-navigation';
import { useCallback, useEffect, useRef, useState } from 'react';
import { ActivityIndicator, Alert, Animated, Platform, Pressable, ScrollView, StyleSheet, Text, View } from 'react-native';

import { Countdown } from '../components/Countdown';
import { QuestionView } from '../components/QuestionView';
import { RewardModal } from '../components/RewardModal';
import { ShareButton } from '../components/ShareButton';
import { Button, Card, MAX_FONT_SCALE, Pill, ProgressBar, Screen, useReduceMotion, useUi } from '../components/ui';
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
import {
  gradeSaved,
  isPendingOtherExam,
  matchesSession,
  QUIT_KEPT_MESSAGE,
  remainingMs,
  restoreSession,
  resumeButtonLabel,
  resumeDecision,
  type RestoredSession,
  type SavedSession,
} from '../lib/resume';
import { quizContext } from '../lib/selectors';
import { dailyShareText, examShareText } from '../lib/share';
import { useProgress } from '../state/progress';
import { clearSession, loadSession, saveSession } from '../state/session';
import { useStyles, useTheme } from '../state/theme';
import { accentText, nativeDriver, radius, type Colors } from '../theme';

const MODES: QuizMode[] = ['chapter', 'daily', 'exam', 'review', 'express'];

/** Quitte le quiz ; ouvert par lien profond, il est seul dans la pile : on revient à l'accueil. */
function exitQuiz() {
  if (router.canGoBack()) router.back();
  else router.replace('/');
}

export default function QuizScreen() {
  // `n` : jeton qui force une nouvelle session (« ⚡ Encore 5 ? »).
  // `resume=1` : reprise demandée depuis l'accueil, sans redemander « Reprendre ou Recommencer ».
  const params = useLocalSearchParams<{ mode?: string; id?: string; n?: string; resume?: string }>();
  const mode: QuizMode = MODES.includes(params.mode as QuizMode) ? (params.mode as QuizMode) : 'chapter';
  const { state, loaded, loadError } = useProgress();
  const { colors } = useTheme();
  if (!loaded) return <ActivityIndicator color={colors.primary} style={{ flex: 1, backgroundColor: colors.bg }} />;
  // Lecture impossible : l'écran d'accueil propose de réessayer, surtout pas d'onboarding.
  if (loadError) return <Redirect href="/" />;
  // Lien profond ouvert avant l'onboarding.
  if (!state.profile) return <Redirect href="/onboarding" />;
  return <QuizSessionScreen key={`${mode}:${params.id ?? ''}:${params.n ?? ''}`} mode={mode} id={params.id} resume={params.resume === '1'} />;
}

/** Ce que l'écran propose à l'ouverture, selon la session enregistrée. */
type Start =
  | { kind: 'fresh' }
  /** Choix « Reprendre » / « Recommencer » affiché dans l'écran. */
  | { kind: 'ask'; restored: RestoredSession }
  | { kind: 'resume'; restored: RestoredSession }
  /** Temps d'examen écoulé ou toutes les questions traitées : noter tout de suite. */
  | { kind: 'grade'; restored: RestoredSession }
  /** Un examen blanc est en cours ailleurs : le reprendre, ou le noter avant de commencer ce quiz. `done` : temps écoulé ou tout répondu. */
  | { kind: 'other'; restored: RestoredSession; done: boolean };

/** Décide de la reprise ; efface la session enregistrée quand elle ne peut plus servir. */
function startFor(saved: SavedSession | null, mode: QuizMode, id: string | undefined, resume: boolean): Start {
  // Examen blanc en cours : il ne doit pas être écrasé par la première réponse de ce quiz sans que l'élève ait choisi.
  if (saved && isPendingOtherExam(saved, mode, id, Date.now())) {
    const restored = restoreSession(saved);
    if (restored && restored.answered > 0) {
      const done = resumeDecision(saved, Date.now()) !== 'resume' || restored.order.length === 0;
      return { kind: 'other', restored, done };
    }
  }
  // Autre session (quiz) : elle reste proposée sur l'accueil jusqu'à la première réponse de celui-ci.
  if (!saved || !matchesSession(saved, mode, id)) return { kind: 'fresh' };
  const decision = resumeDecision(saved, Date.now());
  const restored = decision === 'discard' ? null : restoreSession(saved);
  if (!restored) {
    clearSession();
    return { kind: 'fresh' };
  }
  if (decision !== 'resume' || restored.order.length === 0) return { kind: 'grade', restored };
  return { kind: resume ? 'resume' : 'ask', restored };
}

function QuizSessionScreen({ mode, id, resume }: { mode: QuizMode; id?: string; resume: boolean }) {
  const { colors } = useTheme();
  const { state, finishQuiz } = useProgress();
  const [start, setStart] = useState<Start | null>(null);
  // Note de l'examen blanc en cours, noté avant de commencer ce quiz.
  const [graded, setGraded] = useState<Reward | null>(null);

  useEffect(() => {
    let active = true;
    void loadSession().then((saved) => {
      if (active) setStart(startFor(saved, mode, id, resume));
    });
    return () => {
      active = false;
    };
  }, [mode, id, resume]);

  if (!start) return <ActivityIndicator color={colors.primary} style={{ flex: 1, backgroundColor: colors.bg }} />;
  if (start.kind === 'fresh') {
    return (
      <>
        <FreshQuiz mode={mode} id={id} />
        <RewardModal
          reward={graded}
          onClose={() => setGraded(null)}
          haptics={state.profile?.haptics ?? true}
          icon="📝"
          title="Ton examen blanc est noté"
        />
      </>
    );
  }
  if (start.kind === 'other') {
    const { restored } = start;
    return (
      <OtherExamChoice
        restored={restored}
        done={start.done}
        onGrade={() => {
          // finishQuiz efface aussi la session enregistrée : ce quiz peut commencer.
          const g = gradeSaved(restored);
          const note = examNote(g.results.filter((r) => r.correct).length, g.results.length);
          const r = finishQuiz(g.mode, g.results, g.opts);
          setGraded({ ...r, messages: [`${restored.session.title} : ${formatNote(note)}/20 · ${mention(note)}`, ...r.messages] });
          setStart({ kind: 'fresh' });
        }}
      />
    );
  }
  if (start.kind === 'ask') {
    return (
      <ResumeChoice
        restored={start.restored}
        onResume={() => setStart({ kind: 'resume', restored: start.restored })}
        onRestart={() => {
          clearSession();
          setStart({ kind: 'fresh' });
        }}
      />
    );
  }
  return <RestoredQuiz restored={start.restored} grade={start.kind === 'grade'} />;
}

/** Nouvelle session, tirée au montage. */
function FreshQuiz({ mode, id }: { mode: QuizMode; id?: string }) {
  const ui = useUi();
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
  return <QuizRunner mode={mode} id={id} session={session} practice={practice} day={day} />;
}

/** Session reprise (ou notée tout de suite), reconstruite depuis la sauvegarde. */
function RestoredQuiz({ restored, grade }: { restored: RestoredSession; grade: boolean }) {
  const { state } = useProgress();
  const { saved } = restored;
  // Figé au montage : le défi devient « relevé » dès la fin de cette session.
  const [practice] = useState(() => saved.mode === 'daily' && todayStats(state, saved.day).challengeDone);
  return (
    <QuizRunner
      mode={saved.mode}
      id={saved.id ?? undefined}
      session={restored.session}
      practice={practice}
      day={saved.day}
      restored={restored}
      autoFinish={grade}
    />
  );
}

/** Durée restante lisible : « 7 min », « moins d'une minute ». */
function minutesText(ms: number): string {
  const min = Math.floor(ms / 60_000);
  return min < 1 ? 'moins d’une minute' : `${min} min`;
}

/** Choix affiché à l'ouverture d'un quiz interrompu : reprendre là où on en était, ou recommencer. */
function ResumeChoice({ restored, onResume, onRestart }: { restored: RestoredSession; onResume: () => void; onRestart: () => void }) {
  const { colors } = useTheme();
  const ui = useUi();
  const { saved, session, answered, total } = restored;
  const color = session.color ?? colors.primary;
  const isExam = saved.mode === 'exam';
  const [left] = useState(() => remainingMs(saved, Date.now()));
  return (
    <Screen edges={['bottom']}>
      <Stack.Screen options={{ title: session.title }} />
      <Card style={{ gap: 10 }}>
        <Text style={ui.h2} accessibilityRole="header">
          On reprend là où tu en étais ?
        </Text>
        <Text style={ui.body}>
          {session.title} : {answered}/{total} {isExam ? 'questions répondues' : 'questions faites'}.
        </Text>
        <ProgressBar value={answered / total} color={color} height={10} accessibilityLabel="Questions déjà faites" />
        {left !== null && (
          <Text style={[ui.body, { fontWeight: '700' }]}>
            ⏱️ Le chrono a continué : il te reste {minutesText(left)}.
          </Text>
        )}
        <Button label={`▶ ${resumeButtonLabel(restored)}`} color={color} onPress={onResume} />
        <Button label="Recommencer" variant="secondary" color={color} onPress={onRestart} />
        <Text style={ui.muted}>
          {isExam
            ? 'Recommencer efface tes réponses et lance un nouvel examen blanc.'
            : 'Recommencer efface tes réponses de cette série et repart de zéro.'}
        </Text>
      </Card>
    </Screen>
  );
}

/** Ouverture d'un autre quiz pendant un examen blanc en cours : le reprendre, ou le noter avant de commencer. */
function OtherExamChoice({ restored, done, onGrade }: { restored: RestoredSession; done: boolean; onGrade: () => void }) {
  const { colors } = useTheme();
  const ui = useUi();
  const { saved, session, answered, total } = restored;
  const color = session.color ?? colors.primary;
  const [left] = useState(() => remainingMs(saved, Date.now()));
  // Temps écoulé : l'écran de l'examen le note tout de suite ; sinon il reprend sans redemander.
  const openExam = () =>
    router.replace({ pathname: '/quiz', params: { mode: 'exam', ...(saved.id ? { id: saved.id } : {}), ...(done ? {} : { resume: '1' }) } });
  return (
    <Screen edges={['bottom']}>
      <Stack.Screen options={{ title: 'Examen blanc en cours' }} />
      <Card style={{ gap: 10 }}>
        <Text style={ui.h2} accessibilityRole="header">
          Tu as un examen blanc en cours
        </Text>
        <Text style={ui.body}>
          {session.title} : {answered}/{total} questions répondues.
        </Text>
        <ProgressBar value={answered / total} color={color} height={10} accessibilityLabel="Questions déjà répondues" />
        <Text style={[ui.body, { fontWeight: '700' }]}>
          {done || left === null ? '⏱️ Le temps est écoulé : il ne reste plus qu’à le noter.' : `⏱️ Le chrono continue : il te reste ${minutesText(left)}.`}
        </Text>
        <Button label={done ? '📝 Voir ma note' : '▶ Reprendre l’examen'} color={color} onPress={openExam} />
        <Button label="Le noter et commencer ce quiz" variant="secondary" color={color} onPress={onGrade} />
        <Text style={ui.muted}>Les questions sans réponse comptent comme non traitées.</Text>
      </Card>
    </Screen>
  );
}

function QuizRunner({
  mode,
  id,
  session,
  practice,
  day,
  restored,
  autoFinish = false,
}: {
  mode: QuizMode;
  /** Chapitre (quiz de chapitre) ou matière (examen blanc). */
  id?: string;
  session: QuizSession;
  practice: boolean;
  day: string;
  /** Session reprise : file, réponses, combo et échéance de la sauvegarde. */
  restored?: RestoredSession;
  /** Noter tout de suite la session reprise (temps écoulé ou toutes les questions traitées). */
  autoFinish?: boolean;
}) {
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
  const { state, finishQuiz } = useProgress();
  const navigation = useNavigation();
  const haptics = state.profile?.haptics ?? true;
  const isExam = mode === 'exam';
  const chapterId = mode === 'chapter' ? id : undefined;
  const subjectId = isExam ? id : undefined;
  const total = session.questions.length;
  // File des questions restantes (indices dans session.questions) : « Passer » renvoie la question au bout.
  const [order, setOrder] = useState(() => restored?.order ?? session.questions.map((_, i) => i));
  // Réponse donnée à chaque question (null = pas encore répondu).
  const answers = useRef<(boolean | null)[]>(restored ? [...restored.answers] : session.questions.map(() => null));
  // Suite de bonnes réponses en cours et record de la session (quête « Enchaîne 5 bonnes réponses »).
  const comboRef = useRef(restored?.saved.combo ?? 0);
  const maxCombo = useRef(restored?.saved.maxCombo ?? 0);
  const [combo, setCombo] = useState(restored?.saved.combo ?? 0);
  // Début de la session (fixé à la première sauvegarde, ou au départ du chrono).
  const startedAt = useRef<number | null>(restored?.saved.startedAt ?? null);
  // Une réponse a été validée : la sortie est confirmée (la progression reste gardée).
  const [touched, setTouched] = useState((restored?.answered ?? 0) > 0);
  const [results, setResults] = useState<AnswerResult[]>([]);
  const [reward, setReward] = useState<Reward | null>(null);
  // Examen blanc : échéance fixée au tap sur « Commencer » (null tant que l'introduction est affichée).
  const [deadline, setDeadline] = useState<number | null>(restored?.saved.deadline ?? null);
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
    // finishQuiz efface aussi la session enregistrée.
    setReward(finishQuiz(mode, all, { chapterId, subjectId, day, maxCombo: maxCombo.current }));
  }, [session, finishQuiz, mode, chapterId, subjectId, day]);

  // Reprise d'une session à noter tout de suite (temps d'examen écoulé pendant que l'appli était fermée…).
  useEffect(() => {
    if (autoFinish) finish();
  }, [autoFinish, finish]);

  /** Enregistre la session pour pouvoir la reprendre après une interruption (sans attendre l'écriture). */
  const persist = (queue: number[], end: number | null = deadline) => {
    if (finished.current) return;
    startedAt.current ??= Date.now();
    saveSession({
      mode,
      id: id ?? null,
      day,
      questionIds: session.questions.map((q) => q.question.id),
      order: queue,
      answers: [...answers.current],
      startedAt: startedAt.current,
      ...(end !== null ? { deadline: end } : {}),
      combo: comboRef.current,
      maxCombo: maxCombo.current,
    });
  };

  /** Réponse validée à la question en tête de file : comptée et enregistrée une seule fois. */
  const record = (correct: boolean) => {
    const index = order[0];
    if (index === undefined || answers.current[index] !== null) return;
    answers.current[index] = correct;
    comboRef.current = correct ? comboRef.current + 1 : 0;
    maxCombo.current = Math.max(maxCombo.current, comboRef.current);
    persist(order.slice(1));
  };

  // Alerte « Quitter ? » affichée : une expiration du chrono est mise en attente.
  const quitPrompt = useRef(false);
  const expiredDuringPrompt = useRef(false);
  const onExpire = useCallback(() => {
    if (quitPrompt.current) expiredDuringPrompt.current = true;
    else finish();
  }, [finish]);

  // Sortie confirmée (flèche de l'en-tête, retour Android) dès que l'élève a commencé, jusqu'aux résultats.
  // La session reste enregistrée : l'élève pourra reprendre.
  const started = isExam ? deadline !== null : touched || answered > 0;
  usePreventRemove(started && !reward, ({ data }) => {
    const leave = () => navigation.dispatch(data.action);
    const [title, message, stay] = isExam
      ? ['Quitter l’examen ?', `${QUIT_KEPT_MESSAGE} Le chrono continue de tourner pendant ce temps.`, 'Continuer l’examen']
      : ['Quitter le quiz ?', QUIT_KEPT_MESSAGE, 'Continuer le quiz'];
    if (Platform.OS === 'web') {
      // Pas de boîte de confirmation fiable sur le web (Alert n'y affiche rien, confirm() peut être bloqué
      // dans un cadre) : la session étant gardée, l'élève sort directement et pourra reprendre.
      leave();
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
          onPress: () => {
            quitPrompt.current = false;
            leave();
          },
        },
      ],
      { cancelable: true, onDismiss: stayHere },
    );
  });

  // Correction affichée : la réponse est validée (on l'enregistre) et on fait défiler jusqu'à « Continuer ».
  const onChecked = (correct: boolean) => {
    setTouched(true);
    record(correct);
    requestAnimationFrame(() => scrollRef.current?.scrollToEnd({ animated: true }));
  };

  const next = (correct: boolean) => {
    // Examen blanc : pas de correction, la réponse est validée ici.
    record(correct);
    setCombo(comboRef.current);
    const rest = order.slice(1);
    if (rest.length === 0) finish();
    else {
      setOrder(rest);
      scrollRef.current?.scrollTo({ y: 0, animated: false });
    }
  };
  const skip = () => {
    const rest = [...order.slice(1), order[0]];
    setOrder(rest);
    persist(rest);
    scrollRef.current?.scrollTo({ y: 0, animated: false });
  };
  const startExam = () => {
    startedAt.current = Date.now();
    const end = startedAt.current + (session.timeLimit ?? 0) * 1000;
    setDeadline(end);
    // Examen commencé : il peut être repris (ou noté) même sans aucune réponse.
    persist(order, end);
  };

  if (reward) {
    return <Results mode={mode} session={session} results={results} reward={reward} color={color} practice={practice} day={day} subjectId={subjectId} />;
  }
  if (autoFinish) return <ActivityIndicator color={colors.primary} style={{ flex: 1, backgroundColor: colors.bg }} />;
  if (isExam && deadline === null) {
    return <ExamIntro session={session} color={color} onStart={startExam} />;
  }

  return (
    <Screen edges={['bottom']} scrollRef={scrollRef}>
      <Stack.Screen options={{ title: session.title }} />
      <View style={[ui.row, { justifyContent: 'space-between', flexWrap: 'wrap' }]}>
        <Text style={styles.counter} maxFontSizeMultiplier={MAX_FONT_SCALE}>
          {isExam ? `Répondu ${answered}/${total}` : `Question ${answered + 1}/${total}`}
        </Text>
        <View style={[ui.row, { flexShrink: 1, flexWrap: 'wrap', justifyContent: 'flex-end' }]}>
          {practice && <Pill label="Entraînement · sans XP" color={colors.muted} />}
          {!isExam && <ComboPill combo={combo} />}
          {isExam && deadline !== null && <Countdown deadline={deadline} onExpire={onExpire} />}
        </View>
      </View>
      <ProgressBar value={answered / total} color={color} height={10} accessibilityLabel="Avancement du quiz" />
      {/* Matière seulement : le chapitre serait un indice, il est révélé dans la correction. */}
      {mode !== 'chapter' && !isExam && current && (
        <Text style={[ui.muted, { color: accentText(current.subject.color, colors), fontWeight: '700' }]}>
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
  const ui = useUi();
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
  const { colors } = useTheme();
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
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
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

  const dailyMessage = () => {
    if (!state.profile) return null;
    return dailyShareText({
      day,
      trackLabel: getTrack(state.profile.track).label,
      correct,
      total,
      grid: results.map((r) => (r.correct ? '✅' : '❌')).join(''),
      streak: effectiveStreak(state, day),
    });
  };
  const examMessage = () => {
    const subjectName = getSubject(subjectId ?? '')?.name ?? session.questions[0].subject.name;
    return examShareText({ subjectName, note });
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
            <Text style={[styles.score, { color: accentText(color, colors) }]} maxFontSizeMultiplier={MAX_FONT_SCALE}>
              {formatNote(note)}/20
            </Text>
            <Pill label={`Mention : ${mention(note)}`} color={note >= 10 ? colors.primary : colors.red} />
            <Text style={[ui.muted, { textAlign: 'center' }]}>
              Note indicative sur des questions de cours : l’épreuve réelle comporte aussi des exercices rédigés.
            </Text>
            {noteTrend(state.examHistory[subjectId ?? '']) && (
              <Text style={[ui.body, { fontWeight: '700' }]}>Tes dernières notes : {noteTrend(state.examHistory[subjectId ?? ''])}</Text>
            )}
          </>
        ) : (
          <Text style={[styles.score, { color: accentText(color, colors) }]} maxFontSizeMultiplier={MAX_FONT_SCALE}>
            {correct}/{total}
          </Text>
        )}
        <Text style={[ui.h2, { textAlign: 'center' }]}>{message}</Text>
      </View>

      {mode === 'daily' && !training && <ShareButton label="📤 Partager mon score" color={color} message={dailyMessage} />}
      {mode === 'exam' && <ShareButton label="📤 Partager ma note" color={color} message={examMessage} />}

      <Card style={{ gap: 6, backgroundColor: colors.goldSoft }}>
        <Text style={styles.xp} maxFontSizeMultiplier={MAX_FONT_SCALE}>
          +{reward.xp} XP
        </Text>
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
  const { colors } = useTheme();
  const ui = useUi();
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
              <ProgressBar
                value={ok / n}
                color={weak ? colors.red : ok / n >= 0.8 ? colors.primary : colors.gold}
                accessibilityLabel={`${chapter.title} : ${ok} sur ${n}`}
              />
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
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
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

const makeStyles = (colors: Colors) =>
  StyleSheet.create({
    counter: { fontSize: 15, fontWeight: '800', color: colors.text },
    score: { fontSize: 48, fontWeight: '900' },
    xp: { fontSize: 26, fontWeight: '900', color: colors.primary },
    line: { fontSize: 15, color: colors.text },
    link: { minHeight: 44, justifyContent: 'center', borderRadius: radius.sm },
    linkText: { fontWeight: '700', color: colors.primaryDark },
  });
