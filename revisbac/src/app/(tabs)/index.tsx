import { router, useFocusEffect } from 'expo-router';
import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { Pressable, StyleSheet, Text, View } from 'react-native';

import { QuestsCard } from '../../components/QuestsCard';
import { RewardModal } from '../../components/RewardModal';
import { ShareButton } from '../../components/ShareButton';
import { Button, Card, MAX_FONT_SCALE, Pill, ProgressBar, Screen, SectionTitle, useUi } from '../../components/ui';
import { getChapter, getTrack } from '../../data/catalog';
import { addDays, dayKey, daysBetween, formatDay, weekdayLetter } from '../../lib/dates';
import { effectiveStreak, examNote, formatNote, levelInfo, mention, todayStats, XP, type ProgressState, type Reward } from '../../lib/gamification';
import { nextStep, pace } from '../../lib/plan';
import { DAILY_SIZE, EXPRESS_SIZE } from '../../lib/quizBuilder';
import {
  ficheToContinue,
  gradeSaved,
  remainingMs,
  restoreSession,
  resumeCardText,
  resumeDecision,
  type RestoredSession,
  type SavedSession,
} from '../../lib/resume';
import { activeMistakes, dueCards, dueCardsText, visibleSubjects, type DueCards } from '../../lib/selectors';
import { dailyShareText } from '../../lib/share';
import { useProgress } from '../../state/progress';
import { clearSession, useSavedSession } from '../../state/session';
import { useStyles, useTheme } from '../../state/theme';
import { accentText, lightColors, radius, textOn, type Colors } from '../../theme';

/** Révision express : au plus 2 erreurs dues parmi les 5 questions (même règle que quizBuilder). */
const EXPRESS_MISTAKES = 2;
/** Chapitres listés dans la carte « cartes à revoir » ; les autres sont dans l'onglet Réviser. */
const DUE_CHAPTERS_SHOWN = 3;

const TIPS = [
  'Révise un peu chaque jour plutôt que tout la veille : ta mémoire retient mieux en plusieurs fois.',
  'Après une fiche, fais tout de suite le quiz : se tester est la meilleure façon de mémoriser.',
  'Explique une notion à voix haute comme si tu l’enseignais à un camarade : si tu bloques, relis la fiche.',
  'Tes erreurs sont précieuses : passe régulièrement par « Revoir mes erreurs » dans l’onglet Défis.',
  'Le jour de l’examen, lis tout le sujet avant de commencer et commence par ce que tu maîtrises le mieux.',
  'Dors bien la veille d’un examen : le sommeil consolide ce que tu as appris.',
  'Fais des pauses de 5 minutes toutes les 25 minutes de révision (méthode Pomodoro).',
];

export default function Home() {
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
  const { state } = useProgress();
  const profile = state.profile!;
  const today = dayKey();
  const track = getTrack(profile.track);
  const subjects = visibleSubjects(profile);
  const lvl = levelInfo(state.xp);
  const streak = effectiveStreak(state, today);
  const activeToday = state.streak.lastDay === today;
  const todayData = todayStats(state, today);
  const dailyResult = state.dailyResults[today];
  const mistakes = activeMistakes(state, today);
  const dueCount = mistakes.due.length;
  const plan = pace(state, subjects, profile.examDate, today);
  const step = nextStep(state, subjects, today);
  const stepSubject = step ? subjects.find((s) => s.id === step.subjectId) : undefined;
  const cards = dueCards(state, today);
  const { session: saved, lastFiche, reload } = useSavedSession();
  const now = useFocusNow();
  // Fiche ouverte mais pas lue jusqu'au bout (moins de 7 jours, matière visible).
  const ficheId = ficheToContinue(lastFiche, state, new Set(subjects.flatMap((s) => s.chapters.map((c) => c.id))), now);
  const fiche = ficheId ? getChapter(ficheId) : undefined;
  // La prochaine étape ne répète pas « Continuer ta fiche ».
  const showStep = step && stepSubject && !(step.kind === 'fiche' && step.chapterId === ficheId);

  // Origine des questions de la révision express (voir quizBuilder).
  const anyRead = subjects.some((s) => s.chapters.some((c) => state.fichesRead[c.id]));
  const expressMistakes = Math.min(EXPRESS_MISTAKES, dueCount);
  const expressRest = EXPRESS_SIZE - expressMistakes;
  const expressOrigin =
    expressMistakes === 0
      ? anyRead
        ? `${EXPRESS_SIZE} questions de tes chapitres`
        : `${EXPRESS_SIZE} questions pour découvrir tes matières`
      : `${expressMistakes} erreur${expressMistakes > 1 ? 's' : ''} à revoir + ${expressRest} question${expressRest > 1 ? 's' : ''} ${anyRead ? 'de tes chapitres' : 'pour découvrir tes matières'}`;

  const tip = TIPS[Math.abs(daysBetween('2026-01-01', today)) % TIPS.length];
  const week = Array.from({ length: 7 }, (_, i) => addDays(today, i - 6));

  const dailyMessage = () => (dailyResult ? dailyShareText({ day: today, trackLabel: track.label, ...dailyResult, streak }) : null);

  return (
    <Screen>
      <View style={[ui.row, { justifyContent: 'space-between' }]}>
        <View style={{ flex: 1, gap: 4 }}>
          <Text style={ui.h1}>Salut {profile.name} 👋</Text>
          <Pill label={`${track.emoji} ${track.label}`} />
        </View>
        <View style={styles.streak} accessible accessibilityLabel={`Série : ${streak} jour${streak > 1 ? 's' : ''}${activeToday ? ', activité du jour faite' : ''}`}>
          <Text style={{ fontSize: 26, opacity: streak > 0 ? 1 : 0.35 }} maxFontSizeMultiplier={MAX_FONT_SCALE}>
            🔥
          </Text>
          <Text style={styles.streakNumber} maxFontSizeMultiplier={MAX_FONT_SCALE}>
            {streak}
            {activeToday ? ' ✓' : ''}
          </Text>
        </View>
      </View>
      {!activeToday && streak > 0 && <Text style={[ui.muted, { textAlign: 'right', marginTop: -8 }]}>Une activité aujourd’hui prolonge ta série</Text>}

      <ResumeQuiz saved={saved} now={now} reload={reload} />

      {fiche && (
        <Card
          onPress={() => router.push({ pathname: '/fiche/[id]', params: { id: fiche.chapter.id } })}
          accessibilityLabel={`Continuer ta fiche : ${fiche.chapter.title}, ${fiche.subject.name}`}
        >
          <View style={ui.row}>
            <Text style={{ fontSize: 28 }}>{fiche.subject.icon}</Text>
            <View style={{ flex: 1, gap: 2 }}>
              <Text style={styles.cardTitle}>📖 Continuer ta fiche : {fiche.chapter.title}</Text>
              <Text style={[ui.muted, { color: accentText(fiche.subject.color, colors), fontWeight: '700' }]}>{fiche.subject.name}</Text>
            </View>
            <Text style={styles.chevron}>›</Text>
          </View>
        </Card>
      )}

      <Card
        style={styles.express}
        onPress={() => router.push({ pathname: '/quiz', params: { mode: 'express' } })}
        accessibilityLabel={`Révision express, ${EXPRESS_SIZE} questions, environ 3 minutes : ${expressOrigin}`}
      >
        <View style={ui.row}>
          <View style={{ flex: 1, gap: 2 }}>
            <Text style={styles.expressTitle}>⚡ Révision express · {EXPRESS_SIZE} questions · ~3 min</Text>
            <Text style={styles.expressText}>{expressOrigin}</Text>
          </View>
          <Text style={[styles.chevron, { color: textOn(colors.primary) }]}>›</Text>
        </View>
      </Card>

      <QuestsCard />

      {/* Niveau et objectif du jour réunis : une seule carte de progression. */}
      <Card>
        <View style={[ui.row, { justifyContent: 'space-between' }]}>
          <Text style={styles.levelTitle}>
            Niveau {lvl.level} · {lvl.title}
          </Text>
          <Text style={styles.xp}>{state.xp} XP</Text>
        </View>
        <View style={{ marginVertical: 8 }}>
          <ProgressBar value={lvl.progress} color={colors.gold} height={10} accessibilityLabel={`Progression vers le niveau ${lvl.level + 1}`} />
        </View>
        <Text style={ui.muted}>Encore {lvl.toNext} XP pour le niveau {lvl.level + 1}</Text>
        <View style={styles.divider} />
        <View style={[ui.row, { justifyContent: 'space-between' }]}>
          <Text style={styles.cardTitle}>Objectif du jour</Text>
          <Text style={styles.goalText}>
            {Math.min(todayData.xp, profile.dailyGoal)} / {profile.dailyGoal} XP {todayData.goalBonusGiven ? '✅' : ''}
          </Text>
        </View>
        <View style={{ marginVertical: 8 }}>
          <ProgressBar value={todayData.xp / profile.dailyGoal} height={10} accessibilityLabel="Objectif du jour" />
        </View>
        {/* Une seule phrase pour TalkBack au lieu de 7 pastilles muettes. */}
        <View style={styles.week} accessible accessibilityLabel={`Actif ${week.filter((d) => isActive(state, d)).length} jours sur 7`}>
          {week.map((d) => {
            const active = isActive(state, d);
            const frozen = !active && savedByFreeze(state, d);
            return (
              <View key={d} style={{ alignItems: 'center', gap: 4 }}>
                <View
                  style={[
                    styles.weekDot,
                    active && { backgroundColor: colors.primary },
                    frozen && styles.frozenDot,
                    d === today && { borderColor: colors.gold, borderWidth: 2 },
                  ]}
                >
                  <Text style={{ fontSize: 12 }} maxFontSizeMultiplier={MAX_FONT_SCALE}>
                    {active ? '🔥' : frozen ? '🧊' : ''}
                  </Text>
                </View>
                <Text style={ui.muted} maxFontSizeMultiplier={MAX_FONT_SCALE}>
                  {weekdayLetter(d)}
                </Text>
              </View>
            );
          })}
        </View>
        {state.streak.freezes > 0 && (
          <Text style={[ui.muted, { marginTop: 8 }]}>
            🧊 {state.streak.freezes} gel{state.streak.freezes > 1 ? 's' : ''} de série en réserve (protège ta série si tu rates un jour)
          </Text>
        )}
      </Card>

      <Card style={[styles.challenge, todayData.challengeDone && { backgroundColor: colors.primarySoft }]}>
        <Text style={styles.cardTitle}>🎯 Défi du jour</Text>
        {todayData.challengeDone ? (
          dailyResult ? (
            <>
              <Text style={[ui.body, { fontWeight: '800', color: colors.primary }]}>
                ✓ Défi relevé : {dailyResult.correct}/{dailyResult.total}
              </Text>
              <Text style={ui.muted}>Reviens demain pour un nouveau défi.</Text>
              <ShareButton label="📤 Partager mon score" message={dailyMessage} />
            </>
          ) : (
            <Text style={ui.body}>Bravo, défi relevé ! Reviens demain pour un nouveau défi.</Text>
          )
        ) : (
          <>
            <Text style={ui.body}>
              {DAILY_SIZE} questions surprises sur tes matières du tronc commun. +{XP.dailyChallenge} XP bonus !
            </Text>
            <Button label="Relever le défi" variant="gold" onPress={() => router.push({ pathname: '/quiz', params: { mode: 'daily' } })} />
          </>
        )}
      </Card>

      <Card style={styles.countdown}>
        {plan.daysLeft >= 0 ? (
          <View style={[ui.row, { gap: 14 }]}>
            <Text style={styles.countdownNumber} maxFontSizeMultiplier={MAX_FONT_SCALE}>
              {plan.daysLeft === 0 ? 'Jour J' : `J-${plan.daysLeft}`}
            </Text>
            <View style={{ flex: 1, gap: 2 }}>
              <Text style={[styles.cardTitle, styles.onHero]}>
                {plan.daysLeft === 0 ? `C’est le ${track.label} aujourd’hui` : `avant le ${track.label}`} · Phase : {plan.phase}
              </Text>
              <Text style={styles.countdownText}>{rhythmText(plan, profile.name)}</Text>
              <Text style={styles.countdownDate}>
                {profile.examDate === track.defaultExamDate ? 'Date indicative' : 'Date'} : {formatDay(profile.examDate)} (modifiable dans Profil)
              </Text>
            </View>
          </View>
        ) : (
          <View style={{ gap: 10 }}>
            <Text style={[styles.cardTitle, styles.onHero]}>Ton examen est passé ? Règle la date de ta prochaine session.</Text>
            <Button label="Régler la date dans Profil" variant="secondary" color={colors.primaryDark} onPress={() => router.navigate('/profil')} />
          </View>
        )}
      </Card>

      {showStep && (
        <>
          <SectionTitle>👉 Ta prochaine étape</SectionTitle>
          <Card style={{ gap: 12 }}>
            <View style={ui.row}>
              <Text style={{ fontSize: 30 }}>{stepSubject.icon}</Text>
              <View style={{ flex: 1 }}>
                <Text style={[ui.muted, { color: accentText(stepSubject.color, colors), fontWeight: '700' }]}>{stepSubject.name}</Text>
                <Text style={styles.cardTitle}>{step.title}</Text>
                <Text style={ui.muted}>{step.reason}</Text>
              </View>
            </View>
            {step.kind === 'fiche' ? (
              <Button label="Ouvrir la fiche" onPress={() => router.push({ pathname: '/fiche/[id]', params: { id: step.chapterId } })} />
            ) : (
              <Button label="Lancer le quiz" onPress={() => router.push({ pathname: '/quiz', params: { mode: 'chapter', id: step.chapterId } })} />
            )}
          </Card>
        </>
      )}

      {dueCount > 0 && (
        <Card
          onPress={() => router.push({ pathname: '/quiz', params: { mode: 'review' } })}
          style={{ backgroundColor: colors.redSoft }}
          accessibilityLabel={`À revoir aujourd’hui : ${dueCount} question${dueCount > 1 ? 's' : ''}`}
        >
          <View style={ui.row}>
            <Text style={{ fontSize: 28 }}>🔁</Text>
            <View style={{ flex: 1 }}>
              <Text style={styles.cardTitle}>
                À revoir aujourd’hui : {dueCount} question{dueCount > 1 ? 's' : ''}
              </Text>
              {mistakes.waiting > 0 && (
                <Text style={ui.muted}>
                  {mistakes.waiting} autre{mistakes.waiting > 1 ? 's' : ''} en attente
                </Text>
              )}
            </View>
            <Text style={styles.chevron}>›</Text>
          </View>
        </Card>
      )}

      <DueCardsCard due={cards} />

      <Card style={{ backgroundColor: colors.goldSoft }}>
        <Text style={styles.cardTitle}>💡 Conseil du jour</Text>
        <Text style={[ui.body, { marginTop: 4 }]}>{tip}</Text>
      </Card>
    </Screen>
  );
}

/** Heure relevée à chaque retour sur l'écran : sert à l'expiration des sessions et de la fiche à continuer. */
function useFocusNow(): number {
  const [now, setNow] = useState(() => Date.now());
  useFocusEffect(useCallback(() => setNow(Date.now()), []));
  return now;
}

type ResumeView =
  | { kind: 'none' }
  | { kind: 'drop' }
  | { kind: 'gradeSilently'; restored: RestoredSession }
  | { kind: 'resume' | 'grade'; restored: RestoredSession; remaining: number | null };

/** Que montrer pour la session enregistrée (voir resume.resumeDecision). */
function resumeView(saved: SavedSession | null, now: number): ResumeView {
  if (!saved) return { kind: 'none' };
  const decision = resumeDecision(saved, now);
  if (decision === 'discard') return { kind: 'drop' };
  const restored = restoreSession(saved);
  if (!restored) return { kind: 'drop' };
  if (decision === 'gradeSilently') return { kind: 'gradeSilently', restored };
  return { kind: decision, restored, remaining: remainingMs(saved, now) };
}

/**
 * Quiz ou examen interrompu : « ▶ Reprendre ton quiz : {titre}, {k}/{n} ». La carte reprend le quiz directement (resume=1).
 * Un examen dont le temps est écoulé depuis plus de 30 min est noté ici, sans proposition de reprise.
 */
function ResumeQuiz({ saved, now, reload }: { saved: SavedSession | null; now: number; reload: () => void }) {
  const styles = useStyles(makeStyles);
  const ui = useUi();
  const { state, finishQuiz } = useProgress();
  const [reward, setReward] = useState<Reward | null>(null);
  // Une session n'est notée qu'une fois, même si l'effet repasse avant la relecture.
  const graded = useRef<number | null>(null);
  // finishQuiz change à chaque modification de l'état : la référence évite de relancer l'effet avec une session périmée.
  const finish = useRef(finishQuiz);
  useEffect(() => {
    finish.current = finishQuiz;
  }, [finishQuiz]);
  const view = useMemo(() => resumeView(saved, now), [saved, now]);

  // Effacement et notation seulement quand la session vient d'être relue (jamais sur une copie ancienne).
  useEffect(() => {
    const current = resumeView(saved, Date.now());
    if (current.kind === 'drop') {
      clearSession();
      reload();
    } else if (current.kind === 'gradeSilently' && graded.current !== current.restored.saved.startedAt) {
      graded.current = current.restored.saved.startedAt;
      const g = gradeSaved(current.restored);
      const correct = g.results.filter((r) => r.correct).length;
      const note = examNote(correct, g.results.length);
      const r = finish.current(g.mode, g.results, g.opts);
      setReward({ ...r, messages: [`${current.restored.session.title} : ${formatNote(note)}/20 · ${mention(note)}`, ...r.messages] });
      reload();
    }
  }, [saved, reload]);

  const modal = (
    <RewardModal
      reward={reward}
      onClose={() => setReward(null)}
      haptics={state.profile?.haptics ?? true}
      icon="📝"
      title="Ton examen blanc est noté"
    />
  );
  if (view.kind !== 'resume' && view.kind !== 'grade') return modal;

  const { restored, remaining } = view;
  const { mode, id } = restored.saved;
  // Reprise demandée ici : l'écran de quiz ne redemande pas « Reprendre ou Recommencer ».
  const open = () =>
    router.push({ pathname: '/quiz', params: { mode, ...(id ? { id } : {}), ...(view.kind === 'resume' ? { resume: '1' } : {}) } });
  const minutes = remaining === null ? null : Math.max(1, Math.ceil(remaining / 60000));
  const title = view.kind === 'resume' ? resumeCardText(restored) : `📝 ${restored.session.title} : ${mode === 'exam' ? 'le temps est écoulé' : 'tout est fait'}`;
  const detail =
    view.kind === 'grade'
      ? mode === 'exam'
        ? 'Touche pour voir ta note.'
        : 'Touche pour voir ton résultat.'
      : minutes !== null
        ? `⏱ Il te reste environ ${minutes} min : le chrono de l’examen continue.`
        : 'Tu reprends là où tu t’étais arrêté.';

  return (
    <>
      <Card onPress={open} style={styles.resume} accessibilityLabel={`${title.replace('▶ ', '')}. ${detail}`}>
        <View style={ui.row}>
          <View style={{ flex: 1, gap: 2 }}>
            <Text style={styles.cardTitle}>{title}</Text>
            <Text style={ui.muted}>{detail}</Text>
          </View>
          <Text style={styles.chevron}>›</Text>
        </View>
      </Card>
      {modal}
    </>
  );
}

/** « 🧠 {n} cartes à revoir » : un chapitre → ses flashcards ; plusieurs → un lien par chapitre (les plus chargés d'abord). */
function DueCardsCard({ due }: { due: DueCards }) {
  const styles = useStyles(makeStyles);
  const ui = useUi();
  const text = dueCardsText(due.total);
  if (!text) return null;
  const openDeck = (chapterId: string) => router.push({ pathname: '/flashcards/[id]', params: { id: chapterId } });

  if (due.byChapter.length === 1) {
    const ref = getChapter(due.byChapter[0].chapterId);
    return (
      <Card onPress={() => openDeck(due.byChapter[0].chapterId)} accessibilityLabel={`${text.replace('🧠 ', '')}${ref ? ` : ${ref.chapter.title}` : ''}`}>
        <View style={ui.row}>
          <View style={{ flex: 1, gap: 2 }}>
            <Text style={styles.cardTitle}>{text}</Text>
            {ref && <Text style={ui.muted}>Flashcards : {ref.chapter.title}</Text>}
          </View>
          <Text style={styles.chevron}>›</Text>
        </View>
      </Card>
    );
  }

  const shown = due.byChapter.slice(0, DUE_CHAPTERS_SHOWN);
  const more = due.byChapter.length - shown.length;
  return (
    <Card style={{ gap: 8 }}>
      <Text style={styles.cardTitle} accessibilityRole="header">
        {text}
      </Text>
      <Text style={ui.muted}>Quelques minutes suffisent : commence par le chapitre le plus chargé.</Text>
      {shown.map((c) => {
        const ref = getChapter(c.chapterId);
        const name = ref?.chapter.title ?? c.chapterId;
        return (
          <Pressable
            key={c.chapterId}
            onPress={() => openDeck(c.chapterId)}
            accessibilityRole="button"
            accessibilityLabel={`${name}${ref ? `, ${ref.subject.name}` : ''} : ${c.count} carte${c.count > 1 ? 's' : ''} à revoir`}
            style={({ pressed }) => [styles.dueRow, pressed && { opacity: 0.8 }]}
          >
            <Text style={{ fontSize: 22 }}>{ref?.subject.icon ?? '🧠'}</Text>
            <Text style={[ui.body, { flex: 1, fontWeight: '700' }]}>{name}</Text>
            <Text style={styles.dueCount} maxFontSizeMultiplier={MAX_FONT_SCALE}>
              {c.count}
            </Text>
            <Text style={styles.chevron}>›</Text>
          </Pressable>
        );
      })}
      {more > 0 && (
        <Button
          label={`Et ${more} autre${more > 1 ? 's' : ''} chapitre${more > 1 ? 's' : ''} : voir dans Réviser`}
          variant="ghost"
          onPress={() => router.navigate('/matieres')}
        />
      )}
    </Card>
  );
}

/** Rythme conseillé selon la phase. */
function rhythmText(plan: ReturnType<typeof pace>, name: string): string {
  if (plan.phase === 'Jour J') return `Bonne chance, ${name} ! 🍀`;
  if (plan.phase === 'Veille') return 'Relis tes essentiels et dors tôt 💤';
  if (plan.remaining === 0) return 'Toutes les fiches sont lues : place aux examens blancs et aux révisions !';
  const fiches = plan.remaining > 1 ? `${plan.remaining} fiches restantes` : '1 fiche restante';
  return `${fiches} : ${plan.perWeek} par semaine ${plan.perWeek > 1 ? 'suffisent' : 'suffit'}.`;
}

/** Toute clé de `history` est un jour d'activité, même à 0 XP (voir applyGain). */
function isActive(state: ProgressState, day: string): boolean {
  return day in state.history || state.streak.lastDay === day;
}

/**
 * Jour sans activité couvert par la série qui se termine à streak.lastDay : il a été sauvé par un gel.
 * La série compte les jours actifs ; s'il y en a moins que `current` après ce jour, elle remonte plus loin.
 */
function savedByFreeze(state: ProgressState, day: string): boolean {
  const { lastDay, current } = state.streak;
  if (!lastDay || day >= lastDay || current <= 1) return false;
  const activeBetween = Object.keys(state.history).filter((d) => d > day && d < lastDay).length;
  return activeBetween + 1 < current;
}

const makeStyles = (colors: Colors) =>
  StyleSheet.create({
    streak: { alignItems: 'center', backgroundColor: colors.card, borderRadius: 16, paddingHorizontal: 14, paddingVertical: 6 },
    streakNumber: { fontSize: 18, fontWeight: '900', color: colors.text },
    levelTitle: { fontSize: 16, fontWeight: '800', color: colors.text },
    xp: { fontSize: 16, fontWeight: '900', color: colors.goldText },
    cardTitle: { fontSize: 16, fontWeight: '800', color: colors.text },
    goalText: { fontSize: 14, fontWeight: '700', color: colors.primary },
    week: { flexDirection: 'row', justifyContent: 'space-between', marginTop: 4 },
    weekDot: { minWidth: 30, minHeight: 30, borderRadius: 15, backgroundColor: colors.border, alignItems: 'center', justifyContent: 'center' },
    challenge: { gap: 8, borderWidth: 2, borderColor: colors.gold },
    // Bandeau vert foncé dans les deux thèmes : texte clair.
    countdown: { backgroundColor: colors.hero },
    countdownNumber: { fontSize: 30, fontWeight: '900', color: textOn(colors.hero) },
    countdownText: { color: textOn(colors.hero), fontSize: 14, fontWeight: '600' },
    countdownDate: { color: lightColors.goldSoft, fontSize: 12 },
    onHero: { color: textOn(colors.hero) },
    express: { backgroundColor: colors.primary },
    expressTitle: { fontSize: 16, fontWeight: '800', color: textOn(colors.primary) },
    expressText: { fontSize: 13, color: textOn(colors.primary) },
    frozenDot: { backgroundColor: colors.card, borderWidth: 1, borderColor: colors.border },
    chevron: { fontSize: 28, color: colors.muted },
    divider: { height: 1, backgroundColor: colors.border, marginVertical: 12 },
    resume: { borderWidth: 2, borderColor: colors.primary },
    dueRow: {
      flexDirection: 'row',
      alignItems: 'center',
      gap: 10,
      paddingVertical: 8,
      paddingHorizontal: 10,
      borderRadius: radius.sm,
      borderWidth: 1,
      borderColor: colors.border,
      minHeight: 48,
    },
    dueCount: { fontSize: 15, fontWeight: '900', color: colors.primaryDark, fontVariant: ['tabular-nums'] },
  });
