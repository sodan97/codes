import { router } from 'expo-router';
import { StyleSheet, Text, View } from 'react-native';

import { QuestsCard } from '../../components/QuestsCard';
import { ShareButton } from '../../components/ShareButton';
import { Button, Card, MAX_FONT_SCALE, Screen, SectionTitle, useUi } from '../../components/ui';
import { getSubjects, getTrack, isOptional } from '../../data/catalog';
import { addDays, dayKey, weekdayLetter } from '../../lib/dates';
import { effectiveStreak, formatNote, mention, noteTrend, todayStats, XP, type DailyResult } from '../../lib/gamification';
import { DAILY_SIZE, EXAM_SIZE } from '../../lib/quizBuilder';
import { activeMistakes, visibleSubjects } from '../../lib/selectors';
import { dailyShareText } from '../../lib/share';
import { useProgress } from '../../state/progress';
import { useStyles, useTheme } from '../../state/theme';
import type { Colors } from '../../theme';

export default function Challenges() {
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
  const { state } = useProgress();
  const today = dayKey();
  const todayData = todayStats(state, today);
  const { due, waiting } = activeMistakes(state, today);
  const profile = state.profile!;
  const track = getTrack(profile.track);
  const subjects = visibleSubjects(profile);
  // Le défi du jour est tiré sur le tronc commun : on le précise seulement si l'examen a une LV2.
  const hasOptional = getSubjects(profile.track).some((s) => isOptional(s.id));
  const week = Array.from({ length: 7 }, (_, i) => addDays(today, i - 6));
  const weekXp = week.map((d) => state.history[d] ?? 0);
  const maxXp = Math.max(...weekXp, profile.dailyGoal);
  const official = state.dailyResults[today];
  const record = bestResult(Object.values(state.dailyResults));

  const dailyMessage = () =>
    official ? dailyShareText({ day: today, trackLabel: track.label, ...official, streak: effectiveStreak(state, today) }) : null;

  return (
    <Screen>
      <Text style={ui.h1}>Défis</Text>

      <Card style={[styles.daily, todayData.challengeDone && { backgroundColor: colors.primarySoft }]}>
        <Text style={styles.title}>🎯 Défi du jour</Text>
        <Text style={ui.body}>
          {DAILY_SIZE} questions mélangées sur {hasOptional ? 'tes matières du tronc commun (sans LV2)' : 'tes matières'}. Le même défi pour tous les
          candidats {track.label} aujourd’hui : partage ton score à tes amis !
        </Text>
        {todayData.challengeDone ? (
          <>
            <Text style={[ui.body, { fontWeight: '800', color: colors.primary }]}>
              {official ? `Ton score du jour : ${official.correct}/${official.total}` : '✓ Défi relevé ! Reviens demain.'}
            </Text>
            {official && <ShareButton label="📤 Partager mon score" message={dailyMessage} />}
            <View
              style={styles.days}
              accessible
              accessibilityLabel={`Défis des 7 derniers jours : ${week.filter((d) => state.dailyResults[d]).length} relevés sur 7`}
            >
              {week.map((d) => {
                const r = state.dailyResults[d];
                return (
                  <View key={d} style={styles.day}>
                    <Text style={[ui.muted, d === today && { fontWeight: '900', color: colors.text }]}>{weekdayLetter(d)}</Text>
                    <View style={[styles.dayScore, r && { backgroundColor: colors.card, borderColor: colors.primary }]}>
                      <Text style={[styles.dayScoreText, !r && { color: colors.muted }]} maxFontSizeMultiplier={MAX_FONT_SCALE}>{r ? `${r.correct}/${r.total}` : '–'}</Text>
                    </View>
                  </View>
                );
              })}
            </View>
            {record && (
              <Text style={ui.muted}>
                Record : {record.correct}/{record.total}
              </Text>
            )}
            <Button label="Rejouer pour m’entraîner (sans XP)" variant="secondary" onPress={() => router.push({ pathname: '/quiz', params: { mode: 'daily' } })} />
          </>
        ) : (
          <Button label={`Relever le défi (+${XP.dailyChallenge} XP)`} variant="gold" onPress={() => router.push({ pathname: '/quiz', params: { mode: 'daily' } })} />
        )}
      </Card>

      <QuestsCard />

      <Card style={{ gap: 8 }}>
        <Text style={styles.title}>🔁 Revoir mes erreurs</Text>
        <Text style={ui.body}>
          {due.length} à revoir aujourd’hui · {waiting} en attente. Une question ratée revient le lendemain, puis 3 et 7 jours après si tu la réussis.
        </Text>
        <Button
          label={due.length > 0 ? 'Retravailler mes erreurs' : 'Rien à revoir aujourd’hui 👍'}
          disabled={due.length === 0}
          color={colors.red}
          onPress={() => router.push({ pathname: '/quiz', params: { mode: 'review' } })}
        />
      </Card>

      <Card style={{ gap: 10 }}>
        <Text style={styles.title}>📊 Ma semaine</Text>
        <View
          style={styles.chart}
          accessible
          accessibilityLabel={`XP des 7 derniers jours : ${weekXp.join(', ')}. Objectif atteint ${weekXp.filter((x) => x >= profile.dailyGoal).length} jours sur 7`}
        >
          {week.map((d, i) => (
            <View key={d} style={{ alignItems: 'center', flex: 1, gap: 4 }}>
              <Text style={styles.barValue} maxFontSizeMultiplier={MAX_FONT_SCALE}>
                {weekXp[i] || ''}
              </Text>
              <View style={styles.barTrack}>
                <View style={[styles.bar, { height: `${(weekXp[i] / maxXp) * 100}%`, backgroundColor: weekXp[i] >= profile.dailyGoal ? colors.primary : colors.gold }]} />
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
        const trend = noteTrend(state.examHistory[s.id]);
        return (
          <Card
            key={s.id}
            style={{ gap: 4 }}
            onPress={() => router.push({ pathname: '/quiz', params: { mode: 'exam', id: s.id } })}
            accessibilityLabel={`Examen blanc de ${s.name}, ${best === undefined ? 'pas encore passé' : `meilleure note ${formatNote(best)} sur 20, ${mention(best)}`}`}
          >
            <View style={ui.row}>
              <Text style={{ fontSize: 26 }}>{s.icon}</Text>
              <Text style={[styles.title, { flex: 1 }]}>{s.name}</Text>
              <Text style={{ fontWeight: '800', color: best === undefined ? colors.muted : best >= 10 ? colors.primary : colors.red }}>
                {best === undefined ? 'Pas encore passé' : `${formatNote(best)}/20 · ${mention(best)}`}
              </Text>
            </View>
            {trend && <Text style={[ui.muted, { textAlign: 'right' }]}>Tes dernières notes : {trend}</Text>}
          </Card>
        );
      })}
    </Screen>
  );
}

/** Meilleur score du défi (taux de réussite, puis nombre de bonnes réponses). */
function bestResult(results: DailyResult[]): DailyResult | null {
  let best: DailyResult | null = null;
  for (const r of results) {
    if (r.total === 0) continue;
    if (!best || r.correct / r.total > best.correct / best.total || (r.correct / r.total === best.correct / best.total && r.correct > best.correct)) best = r;
  }
  return best;
}

const makeStyles = (colors: Colors) =>
  StyleSheet.create({
    title: { fontSize: 16, fontWeight: '800', color: colors.text },
    daily: { gap: 8, borderWidth: 2, borderColor: colors.gold },
    // Colonnes étirées sur toute la hauteur : la barre (flex: 1) prend la place restante.
    chart: { flexDirection: 'row', height: 130, gap: 6 },
    barTrack: { width: '70%', flex: 1, justifyContent: 'flex-end' },
    bar: { width: '100%', borderRadius: 6, minHeight: 2 },
    barValue: { fontSize: 11, color: colors.muted, fontWeight: '700' },
    days: { flexDirection: 'row', gap: 4 },
    day: { flex: 1, alignItems: 'center', gap: 4 },
    dayScore: { alignSelf: 'stretch', alignItems: 'center', paddingVertical: 4, borderRadius: 8, borderWidth: 1, borderColor: colors.border },
    dayScoreText: { fontSize: 12, fontWeight: '800', color: colors.text },
  });
