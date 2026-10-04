import { router } from 'expo-router';
import { StyleSheet, Text, View } from 'react-native';

import { Button, Card, Pill, ProgressBar, Screen, SectionTitle, styles as ui } from '../../components/ui';
import { getSubjects, getTrack } from '../../data/catalog';
import { addDays, dayKey, daysBetween, formatDay, weekdayLetter } from '../../lib/dates';
import { effectiveStreak, levelInfo, todayStats } from '../../lib/gamification';
import { useProgress } from '../../state/progress';
import { colors } from '../../theme';

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
  const lvl = levelInfo(state.xp);
  const streak = effectiveStreak(state, today);
  const todayData = todayStats(state, today);
  const daysLeft = daysBetween(today, profile.examDate);
  const mistakes = Object.keys(state.mistakes).length;

  // Prochaine fiche suggérée : une fiche pas encore lue, qui change chaque jour.
  const unread = getSubjects(profile.track).flatMap((s) => s.chapters.filter((c) => !state.fichesRead[c.id]).map((c) => ({ c, s })));
  const suggestion = unread.length ? unread[daysBetween('2026-01-01', today) % unread.length] : null;
  const tip = TIPS[Math.abs(daysBetween('2026-01-01', today)) % TIPS.length];
  const week = Array.from({ length: 7 }, (_, i) => addDays(today, i - 6));

  return (
    <Screen>
      <View style={[ui.row, { justifyContent: 'space-between' }]}>
        <View style={{ flex: 1 }}>
          <Text style={ui.h1}>Salut {profile.name} 👋</Text>
          <Pill label={`${track.emoji} ${track.label}`} />
        </View>
        <View style={styles.streak}>
          <Text style={{ fontSize: 26 }}>{streak > 0 ? '🔥' : '🩶'}</Text>
          <Text style={styles.streakNumber}>{streak}</Text>
        </View>
      </View>

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
            const active = (state.history[d] ?? 0) > 0;
            return (
              <View key={d} style={{ alignItems: 'center', gap: 4 }}>
                <View style={[styles.weekDot, active && { backgroundColor: colors.primary }, d === today && { borderColor: colors.gold, borderWidth: 2 }]}>
                  <Text style={{ fontSize: 12 }}>{active ? '🔥' : ''}</Text>
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
          <Text style={ui.body}>Bravo, défi relevé ! Reviens demain pour un nouveau défi.</Text>
        ) : (
          <>
            <Text style={ui.body}>10 questions surprises sur toutes tes matières. +50 XP bonus !</Text>
            <Button label="Relever le défi" variant="gold" onPress={() => router.push({ pathname: '/quiz', params: { mode: 'daily' } })} />
          </>
        )}
      </Card>

      {daysLeft >= 0 && (
        <Card style={styles.countdown}>
          <Text style={styles.countdownNumber}>J-{daysLeft}</Text>
          <View style={{ flex: 1 }}>
            <Text style={[styles.cardTitle, { color: '#fff' }]}>avant le {track.label}</Text>
            <Text style={{ color: '#ffffffcc', fontSize: 13 }}>Date indicative : {formatDay(profile.examDate)} (modifiable dans Profil)</Text>
          </View>
        </Card>
      )}

      {suggestion && (
        <>
          <SectionTitle>Suggestion du jour</SectionTitle>
          <Card onPress={() => router.push({ pathname: '/fiche/[id]', params: { id: suggestion.c.id } })}>
            <View style={ui.row}>
              <Text style={{ fontSize: 30 }}>{suggestion.s.icon}</Text>
              <View style={{ flex: 1 }}>
                <Text style={[ui.muted, { color: suggestion.s.color, fontWeight: '700' }]}>{suggestion.s.name}</Text>
                <Text style={styles.cardTitle}>{suggestion.c.title}</Text>
                <Text style={ui.muted}>{suggestion.c.summary}</Text>
              </View>
              <Text style={styles.chevron}>›</Text>
            </View>
          </Card>
        </>
      )}

      {mistakes > 0 && (
        <Card onPress={() => router.push({ pathname: '/quiz', params: { mode: 'review' } })} style={{ backgroundColor: colors.redSoft }}>
          <View style={ui.row}>
            <Text style={{ fontSize: 28 }}>🔁</Text>
            <View style={{ flex: 1 }}>
              <Text style={styles.cardTitle}>Revoir mes erreurs</Text>
              <Text style={ui.muted}>
                {mistakes} question{mistakes > 1 ? 's' : ''} à retravailler
              </Text>
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

const styles = StyleSheet.create({
  streak: { alignItems: 'center', backgroundColor: colors.card, borderRadius: 16, paddingHorizontal: 14, paddingVertical: 6 },
  streakNumber: { fontSize: 18, fontWeight: '900', color: colors.text },
  levelTitle: { fontSize: 16, fontWeight: '800', color: colors.text },
  xp: { fontSize: 16, fontWeight: '900', color: colors.gold },
  cardTitle: { fontSize: 16, fontWeight: '800', color: colors.text },
  goalText: { fontSize: 14, fontWeight: '700', color: colors.primary },
  week: { flexDirection: 'row', justifyContent: 'space-between', marginTop: 4 },
  weekDot: { width: 30, height: 30, borderRadius: 15, backgroundColor: colors.border, alignItems: 'center', justifyContent: 'center' },
  challenge: { gap: 8, borderWidth: 2, borderColor: colors.gold },
  countdown: { flexDirection: 'row', alignItems: 'center', gap: 14, backgroundColor: colors.primary },
  countdownNumber: { fontSize: 30, fontWeight: '900', color: colors.gold },
  chevron: { fontSize: 28, color: colors.muted },
});
