import { router } from 'expo-router';
import { useEffect, useState } from 'react';
import { AppState, Linking, Platform, Pressable, StyleSheet, Switch, Text, View } from 'react-native';

import { SubjectPicker } from '../../components/SubjectPicker';
import { Button, Card, ProgressBar, Screen, SectionTitle, styles as ui } from '../../components/ui';
import { getTrack } from '../../data/catalog';
import { addDays, addMonths, dayKey, formatDay } from '../../lib/dates';
import { BADGES, DAILY_GOAL_OPTIONS, effectiveStreak, levelInfo, type ReminderSettings } from '../../lib/gamification';
import {
  DEFAULT_REMINDER_HOUR,
  formatHour,
  REMINDER_HOURS,
  reminderPermission,
  requestReminderPermission,
  type ReminderPermission,
} from '../../lib/reminders';
import { readCount } from '../../lib/selectors';
import { useProgress } from '../../state/progress';
import { colors, radius } from '../../theme';

export default function ProfileScreen() {
  const { state, updateProfile, reset } = useProgress();
  const [confirmReset, setConfirmReset] = useState(false);
  const profile = state.profile!;
  const track = getTrack(profile.track);
  const lvl = levelInfo(state.xp);
  const earned = BADGES.filter((b) => state.badges[b.id]).length;
  const today = dayKey();
  const dateMoves = [
    { label: '−1 mois', date: addMonths(profile.examDate, -1) },
    { label: '−7 j', date: addDays(profile.examDate, -7) },
    { label: '+7 j', date: addDays(profile.examDate, 7) },
    { label: '+1 mois', date: addMonths(profile.examDate, 1) },
  ];

  const stats = [
    { label: 'XP total', value: state.xp, icon: '⭐' },
    { label: 'Série actuelle', value: effectiveStreak(state), icon: '🔥' },
    { label: 'Meilleure série', value: state.streak.best, icon: '🏅' },
    { label: 'Fiches lues', value: readCount(state, profile), icon: '📄' },
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
              onPress={() => updateProfile({ dailyGoal: g })}
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
          {dateMoves.map((m) => {
            // La date ne descend pas sous aujourd'hui.
            const disabled = m.date < profile.examDate && m.date < today;
            return (
              <Pressable
                key={m.label}
                disabled={disabled}
                onPress={() => updateProfile({ examDate: m.date })}
                style={[styles.choice, disabled && { opacity: 0.4 }]}
              >
                <Text style={{ fontWeight: '800', color: colors.text }}>{m.label}</Text>
              </Pressable>
            );
          })}
        </View>
        {profile.examDate !== track.defaultExamDate && (
          <Text style={styles.link} onPress={() => updateProfile({ examDate: track.defaultExamDate })}>
            Revenir à la date indicative ({formatDay(track.defaultExamDate)})
          </Text>
        )}
      </Card>

      <Card style={{ gap: 10 }}>
        <Text style={styles.settingTitle}>Mes matières</Text>
        <SubjectPicker track={profile.track} hidden={profile.hiddenSubjects} onChange={(hiddenSubjects) => updateProfile({ hiddenSubjects })} />
      </Card>

      {Platform.OS !== 'web' && (
        <>
          <ReminderCard reminder={profile.reminder} onChange={(reminder) => updateProfile({ reminder })} />
          <Card style={[ui.row, { justifyContent: 'space-between' }]}>
            <View style={{ flex: 1 }}>
              <Text style={styles.settingTitle}>Vibrations</Text>
              <Text style={ui.muted}>Aux réponses et aux récompenses</Text>
            </View>
            <Switch
              value={profile.haptics}
              onValueChange={(haptics) => updateProfile({ haptics })}
              trackColor={{ true: colors.primary, false: colors.border }}
              thumbColor="#fff"
            />
          </Card>
        </>
      )}

      <Card style={{ gap: 10 }}>
        <Button label="Modifier mon prénom, mon examen ou mes matières" variant="secondary" onPress={() => router.push('/onboarding')} />
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

/** Rappel quotidien : interrupteur et heure. La permission n'est demandée qu'à l'activation. */
function ReminderCard({ reminder, onChange }: { reminder: ReminderSettings | null; onChange: (reminder: ReminderSettings) => void }) {
  const [permission, setPermission] = useState<ReminderPermission | null>(null);
  const [busy, setBusy] = useState(false);
  const enabled = !!reminder?.enabled;
  const hour = reminder?.hour ?? DEFAULT_REMINDER_HOUR;

  // Relue au retour dans l'appli : l'élève a pu l'autoriser dans les réglages du téléphone.
  useEffect(() => {
    let alive = true;
    const refresh = () => {
      reminderPermission().then((p) => alive && setPermission(p));
    };
    refresh();
    const sub = AppState.addEventListener('change', (status) => status === 'active' && refresh());
    return () => {
      alive = false;
      sub.remove();
    };
  }, []);

  const apply = async (on: boolean, h: number) => {
    const minute = h === reminder?.hour ? reminder.minute : 0;
    if (!on) {
      onChange({ enabled: false, hour: h, minute });
      return;
    }
    setBusy(true);
    const granted = await requestReminderPermission();
    setBusy(false);
    setPermission(granted ? 'granted' : 'denied');
    onChange({ enabled: granted, hour: h, minute });
  };

  return (
    <Card style={{ gap: 10 }}>
      <View style={[ui.row, { justifyContent: 'space-between' }]}>
        <View style={{ flex: 1 }}>
          <Text style={styles.settingTitle}>Rappel quotidien</Text>
          <Text style={ui.muted}>{enabled ? `Chaque jour à ${formatHour(hour, reminder?.minute)}, au plus une fois` : 'Désactivé'}</Text>
        </View>
        <Switch
          value={enabled}
          disabled={busy}
          onValueChange={(on) => void apply(on, hour)}
          trackColor={{ true: colors.primary, false: colors.border }}
          thumbColor="#fff"
        />
      </View>
      <View style={styles.hours}>
        {REMINDER_HOURS.map((h) => {
          const on = enabled && hour === h;
          return (
            <Pressable key={h} disabled={busy} onPress={() => void apply(true, h)} style={[styles.hour, on && styles.choiceOn]}>
              <Text style={{ fontWeight: '800', color: on ? '#fff' : colors.text }}>{formatHour(h)}</Text>
            </Pressable>
          );
        })}
      </View>
      {permission === 'denied' && (
        <>
          <Text style={[ui.body, { color: colors.red }]}>Les notifications sont bloquées pour RéviBac : autorise-les dans les réglages du téléphone.</Text>
          <Button label="Ouvrir les réglages du téléphone" variant="ghost" onPress={() => void Linking.openSettings().catch(() => {})} />
        </>
      )}
    </Card>
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
  choiceOn: { backgroundColor: colors.primary, borderColor: colors.primary },
  hours: { flexDirection: 'row', flexWrap: 'wrap', gap: 8 },
  hour: { borderWidth: 2, borderColor: colors.border, borderRadius: radius.sm, paddingVertical: 10, paddingHorizontal: 12, backgroundColor: colors.card },
  link: { fontSize: 14, fontWeight: '700', color: colors.primary, textDecorationLine: 'underline' },
});
