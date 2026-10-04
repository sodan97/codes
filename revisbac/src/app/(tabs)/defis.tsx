import { router } from 'expo-router';
import { StyleSheet, Text, View } from 'react-native';

import { Button, Card, Screen, SectionTitle, styles as ui } from '../../components/ui';
import { getSubjects } from '../../data/catalog';
import { addDays, dayKey, weekdayLetter } from '../../lib/dates';
import { mention, todayStats, XP } from '../../lib/gamification';
import { DAILY_SIZE, EXAM_SIZE } from '../../lib/quizBuilder';
import { useProgress } from '../../state/progress';
import { colors } from '../../theme';

export default function Challenges() {
  const { state } = useProgress();
  const today = dayKey();
  const todayData = todayStats(state, today);
  const mistakes = Object.keys(state.mistakes).length;
  const subjects = getSubjects(state.profile!.track);
  const week = Array.from({ length: 7 }, (_, i) => addDays(today, i - 6));
  const weekXp = week.map((d) => state.history[d] ?? 0);
  const maxXp = Math.max(...weekXp, state.profile!.dailyGoal);

  return (
    <Screen>
      <Text style={ui.h1}>Défis</Text>

      <Card style={[styles.daily, todayData.challengeDone && { backgroundColor: colors.primarySoft }]}>
        <Text style={styles.title}>🎯 Défi du jour</Text>
        <Text style={ui.body}>
          {DAILY_SIZE} questions mélangées sur toutes tes matières. Le même défi pour tous les candidats aujourd’hui : compare ton score avec tes amis !
        </Text>
        {todayData.challengeDone ? (
          <>
            <Text style={[ui.body, { fontWeight: '800', color: colors.primary }]}>✓ Défi relevé ! Reviens demain.</Text>
            <Button label="Rejouer pour m’entraîner" variant="secondary" onPress={() => router.push({ pathname: '/quiz', params: { mode: 'daily' } })} />
          </>
        ) : (
          <Button label={`Relever le défi (+${XP.dailyChallenge} XP)`} variant="gold" onPress={() => router.push({ pathname: '/quiz', params: { mode: 'daily' } })} />
        )}
      </Card>

      <Card style={{ gap: 8 }}>
        <Text style={styles.title}>🔁 Revoir mes erreurs</Text>
        <Text style={ui.body}>
          {mistakes > 0
            ? `${mistakes} question${mistakes > 1 ? 's' : ''} à retravailler. Une question réussie ici sort de la liste.`
            : 'Aucune erreur en attente. Les questions que tu rates aux quiz apparaîtront ici.'}
        </Text>
        <Button label="Retravailler mes erreurs" disabled={mistakes === 0} color={colors.red} onPress={() => router.push({ pathname: '/quiz', params: { mode: 'review' } })} />
      </Card>

      <Card style={{ gap: 10 }}>
        <Text style={styles.title}>📊 Ma semaine</Text>
        <View style={styles.chart}>
          {week.map((d, i) => (
            <View key={d} style={{ alignItems: 'center', flex: 1, gap: 4 }}>
              <Text style={styles.barValue}>{weekXp[i] || ''}</Text>
              <View style={styles.barTrack}>
                <View style={[styles.bar, { height: `${(weekXp[i] / maxXp) * 100}%`, backgroundColor: weekXp[i] >= state.profile!.dailyGoal ? colors.primary : colors.gold }]} />
              </View>
              <Text style={[ui.muted, d === today && { fontWeight: '900', color: colors.text }]}>{weekdayLetter(d)}</Text>
            </View>
          ))}
        </View>
        <Text style={ui.muted}>Total : {weekXp.reduce((a, b) => a + b, 0)} XP sur 7 jours · vert = objectif atteint</Text>
      </Card>

      <SectionTitle>📝 Examens blancs</SectionTitle>
      <Text style={ui.muted}>{EXAM_SIZE} questions chronométrées par matière, notées sur 20 avec mention.</Text>
      {subjects.map((s) => {
        const best = state.examBest[s.id];
        return (
          <Card key={s.id} onPress={() => router.push({ pathname: '/quiz', params: { mode: 'exam', id: s.id } })}>
            <View style={ui.row}>
              <Text style={{ fontSize: 26 }}>{s.icon}</Text>
              <Text style={[styles.title, { flex: 1 }]}>{s.name}</Text>
              <Text style={{ fontWeight: '800', color: best === undefined ? colors.muted : best >= 10 ? colors.primary : colors.red }}>
                {best === undefined ? 'Pas encore passé' : `${best}/20 · ${mention(best)}`}
              </Text>
            </View>
          </Card>
        );
      })}
    </Screen>
  );
}

const styles = StyleSheet.create({
  title: { fontSize: 16, fontWeight: '800', color: colors.text },
  daily: { gap: 8, borderWidth: 2, borderColor: colors.gold },
  chart: { flexDirection: 'row', height: 130, alignItems: 'flex-end', gap: 6 },
  barTrack: { width: '70%', flex: 1, justifyContent: 'flex-end' },
  bar: { width: '100%', borderRadius: 6, minHeight: 2 },
  barValue: { fontSize: 11, color: colors.muted, fontWeight: '700' },
});
