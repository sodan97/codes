import { router } from 'expo-router';
import { useState } from 'react';
import { Pressable, StyleSheet, Text, View } from 'react-native';

import { Button, Card, ProgressBar, Screen, SectionTitle, styles as ui } from '../../components/ui';
import { getTrack } from '../../data/catalog';
import { addDays, formatDay } from '../../lib/dates';
import { BADGES, DAILY_GOAL_OPTIONS, effectiveStreak, levelInfo } from '../../lib/gamification';
import { useProgress } from '../../state/progress';
import { colors, radius } from '../../theme';

export default function ProfileScreen() {
  const { state, setProfile, reset } = useProgress();
  const [confirmReset, setConfirmReset] = useState(false);
  const profile = state.profile!;
  const track = getTrack(profile.track);
  const lvl = levelInfo(state.xp);
  const earned = BADGES.filter((b) => state.badges[b.id]).length;

  const stats = [
    { label: 'XP total', value: state.xp, icon: '⭐' },
    { label: 'Série actuelle', value: effectiveStreak(state), icon: '🔥' },
    { label: 'Meilleure série', value: state.streak.best, icon: '🏅' },
    { label: 'Fiches lues', value: Object.keys(state.fichesRead).length, icon: '📄' },
    { label: 'Quiz terminés', value: state.quizCount, icon: '✅' },
    { label: 'Sans faute', value: state.perfectCount, icon: '💯' },
    { label: 'Défis relevés', value: state.challengesDone, icon: '🎯' },
    { label: 'Gels de série', value: state.streak.freezes, icon: '🧊' },
  ];

  return (
    <Screen>
      <Card style={{ alignItems: 'center', gap: 6 }}>
        <View style={styles.avatar}>
          <Text style={styles.avatarText}>{profile.name.slice(0, 1).toUpperCase()}</Text>
        </View>
        <Text style={ui.h2}>{profile.name}</Text>
        <Text style={ui.muted}>
          {track.emoji} {track.label} · Niveau {lvl.level} · {lvl.title}
        </Text>
        <View style={{ alignSelf: 'stretch', marginTop: 6 }}>
          <ProgressBar value={lvl.progress} color={colors.gold} height={10} />
        </View>
      </Card>

      <View style={styles.grid}>
        {stats.map((s) => (
          <View key={s.label} style={styles.stat}>
            <Text style={{ fontSize: 22 }}>{s.icon}</Text>
            <Text style={styles.statValue}>{s.value}</Text>
            <Text style={styles.statLabel}>{s.label}</Text>
          </View>
        ))}
      </View>

      <SectionTitle right={<Text style={ui.muted}>{earned}/{BADGES.length}</Text>}>Badges</SectionTitle>
      <View style={styles.grid}>
        {BADGES.map((b) => {
          const got = !!state.badges[b.id];
          return (
            <View key={b.id} style={[styles.badge, !got && { opacity: 0.45 }]}>
              <Text style={{ fontSize: 30 }}>{got ? b.icon : '🔒'}</Text>
              <Text style={styles.badgeName}>{b.name}</Text>
              <Text style={styles.badgeDesc}>{b.description}</Text>
            </View>
          );
        })}
      </View>

      <SectionTitle>Réglages</SectionTitle>
      <Card style={{ gap: 10 }}>
        <Text style={styles.settingTitle}>Objectif quotidien</Text>
        <View style={{ flexDirection: 'row', gap: 8 }}>
          {DAILY_GOAL_OPTIONS.map((g) => (
            <Pressable
              key={g}
              onPress={() => setProfile({ ...profile, dailyGoal: g })}
              style={[styles.choice, profile.dailyGoal === g && { backgroundColor: colors.primary, borderColor: colors.primary }]}
            >
              <Text style={{ fontWeight: '800', color: profile.dailyGoal === g ? '#fff' : colors.text }}>{g} XP</Text>
            </Pressable>
          ))}
        </View>

        <Text style={[styles.settingTitle, { marginTop: 6 }]}>Date de l’examen</Text>
        <Text style={ui.body}>{formatDay(profile.examDate)}</Text>
        <Text style={ui.muted}>Ajuste-la quand le calendrier officiel est publié.</Text>
        <View style={{ flexDirection: 'row', gap: 8 }}>
          {[-7, -1, 1, 7].map((n) => (
            <Pressable key={n} onPress={() => setProfile({ ...profile, examDate: addDays(profile.examDate, n) })} style={styles.choice}>
              <Text style={{ fontWeight: '800', color: colors.text }}>
                {n > 0 ? '+' : '−'}
                {Math.abs(n)} j
              </Text>
            </Pressable>
          ))}
        </View>

        <Button label="Changer de prénom ou d’examen" variant="secondary" onPress={() => router.push('/onboarding')} style={{ marginTop: 6 }} />
        {confirmReset ? (
          <View style={{ gap: 8 }}>
            <Text style={[ui.body, { color: colors.red, fontWeight: '700' }]}>Tout ton XP, tes badges et ta série seront effacés. Sûr ?</Text>
            <View style={{ flexDirection: 'row', gap: 8 }}>
              <Button label="Annuler" variant="secondary" style={{ flex: 1 }} onPress={() => setConfirmReset(false)} />
              <Button label="Oui, effacer" color={colors.red} style={{ flex: 1 }} onPress={() => reset()} />
            </View>
          </View>
        ) : (
          <Button label="Réinitialiser ma progression" variant="ghost" color={colors.red} onPress={() => setConfirmReset(true)} />
        )}
      </Card>

      <Text style={[ui.muted, { textAlign: 'center' }]}>
        RéviBac · contenu de révision conforme aux grandes lignes du programme sénégalais, à compléter et valider avec des enseignants.
      </Text>
    </Screen>
  );
}

const styles = StyleSheet.create({
  avatar: { width: 72, height: 72, borderRadius: 36, backgroundColor: colors.primary, alignItems: 'center', justifyContent: 'center' },
  avatarText: { fontSize: 32, fontWeight: '900', color: '#fff' },
  grid: { flexDirection: 'row', flexWrap: 'wrap', gap: 10 },
  stat: { width: '23%', flexGrow: 1, minWidth: 76, backgroundColor: colors.card, borderRadius: radius.md, padding: 10, alignItems: 'center' },
  statValue: { fontSize: 20, fontWeight: '900', color: colors.text },
  statLabel: { fontSize: 11, color: colors.muted, textAlign: 'center' },
  badge: { width: '30%', flexGrow: 1, minWidth: 100, backgroundColor: colors.card, borderRadius: radius.md, padding: 10, alignItems: 'center', gap: 2 },
  badgeName: { fontSize: 13, fontWeight: '800', color: colors.text, textAlign: 'center' },
  badgeDesc: { fontSize: 11, color: colors.muted, textAlign: 'center' },
  settingTitle: { fontSize: 15, fontWeight: '800', color: colors.text },
  choice: { flex: 1, borderWidth: 2, borderColor: colors.border, borderRadius: radius.sm, paddingVertical: 10, alignItems: 'center', backgroundColor: colors.card },
});
