import { router } from 'expo-router';
import { StyleSheet, Text, View } from 'react-native';

import { Button, Card, Pill, ProgressBar, Screen, SectionTitle, styles as ui } from '../../components/ui';
import { getTrack } from '../../data/catalog';
import { addDays, dayKey, daysBetween, formatDay, weekdayLetter } from '../../lib/dates';
import { effectiveStreak, levelInfo, todayStats, XP, type ProgressState } from '../../lib/gamification';
import { nextStep, pace } from '../../lib/plan';
import { DAILY_SIZE, EXPRESS_SIZE } from '../../lib/quizBuilder';
import { activeMistakes, visibleSubjects } from '../../lib/selectors';
import { dailyShareText, share } from '../../lib/share';
import { useProgress } from '../../state/progress';
import { colors } from '../../theme';

/** Révision express : au plus 2 erreurs dues parmi les 5 questions (même règle que quizBuilder). */
const EXPRESS_MISTAKES = 2;

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

  const shareDaily = () => {
    if (!dailyResult) return;
    void share(dailyShareText({ day: today, trackLabel: track.label, ...dailyResult, streak }));
  };

  return (
    <Screen>
      <View style={[ui.row, { justifyContent: 'space-between' }]}>
        <View style={{ flex: 1, gap: 4 }}>
          <Text style={ui.h1}>Salut {profile.name} 👋</Text>
          <Pill label={`${track.emoji} ${track.label}`} />
        </View>
        <View style={styles.streak}>
          <Text style={{ fontSize: 26, opacity: streak > 0 ? 1 : 0.35 }}>🔥</Text>
          <Text style={styles.streakNumber}>
            {streak}
            {activeToday ? ' ✓' : ''}
          </Text>
        </View>
      </View>
      {!activeToday && streak > 0 && <Text style={[ui.muted, { textAlign: 'right', marginTop: -8 }]}>Une activité aujourd’hui prolonge ta série</Text>}

      <Card style={styles.express} onPress={() => router.push({ pathname: '/quiz', params: { mode: 'express' } })}>
        <View style={ui.row}>
          <View style={{ flex: 1, gap: 2 }}>
            <Text style={styles.expressTitle}>⚡ Révision express · {EXPRESS_SIZE} questions · ~3 min</Text>
            <Text style={styles.expressText}>{expressOrigin}</Text>
          </View>
          <Text style={[styles.chevron, { color: '#fff' }]}>›</Text>
        </View>
      </Card>

      <Card>
        <View style={[ui.row, { justifyContent: 'space-between' }]}>
          <Text style={styles.levelTitle}>
            Niveau {lvl.level} · {lvl.title}
          </Text>
          <Text style={styles.xp}>{state.xp} XP</Text>
        </View>
        <View style={{ marginVertical: 8 }}>
          <ProgressBar value={lvl.progress} color={colors.gold} height={10} />
        </View>
        <Text style={ui.muted}>Encore {lvl.toNext} XP pour le niveau {lvl.level + 1}</Text>
      </Card>

      <Card>
        <View style={[ui.row, { justifyContent: 'space-between' }]}>
          <Text style={styles.cardTitle}>Objectif du jour</Text>
          <Text style={styles.goalText}>
            {Math.min(todayData.xp, profile.dailyGoal)} / {profile.dailyGoal} XP {todayData.goalBonusGiven ? '✅' : ''}
          </Text>
        </View>
        <View style={{ marginVertical: 8 }}>
          <ProgressBar value={todayData.xp / profile.dailyGoal} height={10} />
        </View>
        <View style={styles.week}>
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
                  <Text style={{ fontSize: 12 }}>{active ? '🔥' : frozen ? '🧊' : ''}</Text>
                </View>
                <Text style={ui.muted}>{weekdayLetter(d)}</Text>
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
              <Button label="📤 Partager mon score" variant="secondary" onPress={shareDaily} />
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
            <Text style={styles.countdownNumber}>{plan.daysLeft === 0 ? 'Jour J' : `J-${plan.daysLeft}`}</Text>
            <View style={{ flex: 1, gap: 2 }}>
              <Text style={[styles.cardTitle, { color: '#fff' }]}>
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
            <Text style={[styles.cardTitle, { color: '#fff' }]}>Ton examen est passé ? Règle la date de ta prochaine session.</Text>
            <Button label="Régler la date dans Profil" variant="secondary" color={colors.primaryDark} onPress={() => router.navigate('/profil')} />
          </View>
        )}
      </Card>

      {step && stepSubject && (
        <>
          <SectionTitle>👉 Ta prochaine étape</SectionTitle>
          <Card style={{ gap: 12 }}>
            <View style={ui.row}>
              <Text style={{ fontSize: 30 }}>{stepSubject.icon}</Text>
              <View style={{ flex: 1 }}>
                <Text style={[ui.muted, { color: stepSubject.color, fontWeight: '700' }]}>{stepSubject.name}</Text>
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
        <Card onPress={() => router.push({ pathname: '/quiz', params: { mode: 'review' } })} style={{ backgroundColor: colors.redSoft }}>
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

      <Card style={{ backgroundColor: colors.goldSoft }}>
        <Text style={styles.cardTitle}>💡 Conseil du jour</Text>
        <Text style={[ui.body, { marginTop: 4 }]}>{tip}</Text>
      </Card>
    </Screen>
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

const styles = StyleSheet.create({
  streak: { alignItems: 'center', backgroundColor: colors.card, borderRadius: 16, paddingHorizontal: 14, paddingVertical: 6 },
  streakNumber: { fontSize: 18, fontWeight: '900', color: colors.text },
  levelTitle: { fontSize: 16, fontWeight: '800', color: colors.text },
  xp: { fontSize: 16, fontWeight: '900', color: colors.goldText },
  cardTitle: { fontSize: 16, fontWeight: '800', color: colors.text },
  goalText: { fontSize: 14, fontWeight: '700', color: colors.primary },
  week: { flexDirection: 'row', justifyContent: 'space-between', marginTop: 4 },
  weekDot: { width: 30, height: 30, borderRadius: 15, backgroundColor: colors.border, alignItems: 'center', justifyContent: 'center' },
  challenge: { gap: 8, borderWidth: 2, borderColor: colors.gold },
  countdown: { backgroundColor: colors.primaryDark },
  countdownNumber: { fontSize: 30, fontWeight: '900', color: '#fff' },
  countdownText: { color: '#fff', fontSize: 14, fontWeight: '600' },
  countdownDate: { color: colors.goldSoft, fontSize: 12 },
  express: { backgroundColor: colors.primary },
  expressTitle: { fontSize: 16, fontWeight: '800', color: '#fff' },
  expressText: { fontSize: 13, color: '#fff' },
  frozenDot: { backgroundColor: colors.card, borderWidth: 1, borderColor: colors.border },
  chevron: { fontSize: 28, color: colors.muted },
});
