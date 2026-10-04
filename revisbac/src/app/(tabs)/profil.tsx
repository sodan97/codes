import { router } from 'expo-router';
import { useEffect, useState } from 'react';
import { AppState, Linking, Platform, Pressable, Share, StyleSheet, Switch, Text, TextInput, View } from 'react-native';

import { RestoreBackup } from '../../components/RestoreBackup';
import { SubjectPicker } from '../../components/SubjectPicker';
import { Button, Card, MAX_FONT_SCALE, ProgressBar, Screen, SectionTitle, useUi } from '../../components/ui';
import { getTrack } from '../../data/catalog';
import { BACKUP_SAFE_LENGTH, backupShareMessage, exportCode } from '../../lib/backup';
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
import { useStyles, useTheme, type ThemePreference } from '../../state/theme';
import { lightColors, radius, textOn, type Colors } from '../../theme';

export default function ProfileScreen() {
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
  const { state, updateProfile, reset } = useProgress();
  const [confirmReset, setConfirmReset] = useState(false);
  const [restoring, setRestoring] = useState(false);
  const [restored, setRestored] = useState<string | null>(null);
  const [backup, setBackup] = useState<BackupState | null>(null);
  const profile = state.profile!;
  const track = getTrack(profile.track);
  const lvl = levelInfo(state.xp);
  const earned = BADGES.filter((b) => state.badges[b.id]).length;
  const today = dayKey();
  const dateMoves = [
    { label: '−1 mois', a11y: 'Un mois plus tôt', date: addMonths(profile.examDate, -1) },
    { label: '−7 j', a11y: 'Une semaine plus tôt', date: addDays(profile.examDate, -7) },
    { label: '+7 j', a11y: 'Une semaine plus tard', date: addDays(profile.examDate, 7) },
    { label: '+1 mois', a11y: 'Un mois plus tard', date: addMonths(profile.examDate, 1) },
  ];

  // Partage du code de sauvegarde. Le résultat du partage n'est fiable que sur iOS : sur Android, Share.share
  // résout toujours sharedAction, même si la feuille est fermée sans rien envoyer ; sur le web, le partage du
  // navigateur est absent ou peu fiable. Ailleurs qu'un envoi confirmé sur iOS, le code s'affiche donc pour être copié.
  // Code très long : il est toujours affiché, un message risquerait de le couper.
  const saveBackup = async (place: BackupState['place']) => {
    const code = exportCode(state);
    let shared = false;
    let opened = false;
    // Sur le web, le code est simplement affiché pour être copié.
    if (Platform.OS !== 'web') {
      try {
        const result = await Share.share({ message: backupShareMessage(code) });
        opened = true;
        shared = Platform.OS === 'ios' && result.action === Share.sharedAction;
      } catch (e) {
        console.warn('Partage de la sauvegarde impossible', e);
      }
    }
    const long = code.length > BACKUP_SAFE_LENGTH;
    setBackup({ place, code, showCode: !shared || long, shareOpened: opened, long });
  };

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
        <View style={styles.avatar} importantForAccessibility="no-hide-descendants" accessibilityElementsHidden>
          <Text style={styles.avatarText} maxFontSizeMultiplier={MAX_FONT_SCALE}>{profile.name.slice(0, 1).toUpperCase()}</Text>
        </View>
        <Text style={ui.h2}>{profile.name}</Text>
        <Text style={ui.muted}>
          {track.emoji} {track.label} · Niveau {lvl.level} · {lvl.title}
        </Text>
        <View style={{ alignSelf: 'stretch', marginTop: 6 }}>
          <ProgressBar value={lvl.progress} color={colors.gold} height={10} accessibilityLabel={`Progression vers le niveau ${lvl.level + 1}`} />
        </View>
      </Card>

      <View style={styles.grid}>
        {stats.map((s) => (
          <View key={s.label} style={styles.stat} accessible accessibilityLabel={`${s.label} : ${s.value}`}>
            <Text style={{ fontSize: 22 }} maxFontSizeMultiplier={MAX_FONT_SCALE}>
              {s.icon}
            </Text>
            <Text style={styles.statValue} maxFontSizeMultiplier={MAX_FONT_SCALE}>
              {s.value}
            </Text>
            <Text style={styles.statLabel} maxFontSizeMultiplier={MAX_FONT_SCALE}>
              {s.label}
            </Text>
          </View>
        ))}
      </View>

      <SectionTitle right={<Text style={ui.muted}>{earned}/{BADGES.length}</Text>}>Badges</SectionTitle>
      <View style={styles.grid}>
        {BADGES.map((b) => {
          const got = !!state.badges[b.id];
          // Badge verrouillé à compteur : on montre le chemin déjà fait (« 7/10 fiches »).
          const progress = !got && b.progress ? b.progress(state) : null;
          return (
            <View
              key={b.id}
              style={styles.badge}
              accessible
              accessibilityLabel={`${b.name}, ${got ? 'obtenu' : 'à débloquer'} : ${b.description}${progress ? `, ${progress.label}` : ''}`}
            >
              <Text style={[{ fontSize: 30 }, !got && { opacity: 0.45 }]} maxFontSizeMultiplier={MAX_FONT_SCALE}>
                {got ? b.icon : '🔒'}
              </Text>
              <Text style={[styles.badgeName, !got && { color: colors.muted }]}>{b.name}</Text>
              <Text style={styles.badgeDesc}>{b.description}</Text>
              {progress && (
                <View style={styles.badgeProgress}>
                  <ProgressBar value={progress.target ? progress.value / progress.target : 0} color={colors.gold} height={5} />
                  <Text style={styles.badgeProgressText} maxFontSizeMultiplier={MAX_FONT_SCALE}>
                    {progress.label}
                  </Text>
                </View>
              )}
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
              accessibilityRole="button"
              accessibilityState={{ selected: profile.dailyGoal === g }}
              style={[styles.choice, profile.dailyGoal === g && styles.choiceOn]}
            >
              <Text style={[styles.choiceText, profile.dailyGoal === g && styles.choiceTextOn]}>{g} XP</Text>
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
                accessibilityRole="button"
                accessibilityLabel={m.a11y}
                accessibilityState={{ disabled }}
                style={[styles.choice, disabled && { opacity: 0.4 }]}
              >
                <Text style={styles.choiceText}>{m.label}</Text>
              </Pressable>
            );
          })}
        </View>
        {profile.examDate !== track.defaultExamDate && (
          <Text style={styles.link} accessibilityRole="link" onPress={() => updateProfile({ examDate: track.defaultExamDate })}>
            Revenir à la date indicative ({formatDay(track.defaultExamDate)})
          </Text>
        )}
      </Card>

      <Card style={{ gap: 10 }}>
        <Text style={styles.settingTitle}>Mes matières</Text>
        <SubjectPicker track={profile.track} hidden={profile.hiddenSubjects} onChange={(hiddenSubjects) => updateProfile({ hiddenSubjects })} />
      </Card>

      <AppearanceCard />

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
              thumbColor={lightColors.card}
            />
          </Card>
        </>
      )}

      <Card style={{ gap: 10 }}>
        <Text style={styles.settingTitle}>Sauvegarde</Text>
        <Text style={ui.muted}>
          Ta progression reste sur ce téléphone. Garde un code de sauvegarde pour la retrouver après une réinstallation ou sur un autre téléphone,
          sans compte ni Internet.
        </Text>
        {restored ? (
          <Text style={[ui.body, { color: colors.primaryDark, fontWeight: '700' }]} accessibilityLiveRegion="polite">
            ✓ Progression restaurée. Content de te revoir, {restored} !
          </Text>
        ) : null}
        {restoring ? (
          <RestoreBackup
            onRestored={(name) => {
              setRestoring(false);
              setRestored(name);
              setBackup(null);
            }}
            onCancel={() => setRestoring(false)}
          />
        ) : (
          <>
            <Button label="💾 Sauvegarder ma progression" onPress={() => void saveBackup('settings')} />
            {backup?.place === 'settings' && <BackupResult backup={backup} />}
            <Button
              label="Restaurer une sauvegarde"
              variant="secondary"
              onPress={() => {
                setRestoring(true);
                setRestored(null);
              }}
            />
          </>
        )}
      </Card>

      <Card style={{ gap: 10 }}>
        <Button label="Modifier mon prénom, mon examen ou mes matières" variant="secondary" onPress={() => router.push('/onboarding')} />
        {confirmReset ? (
          <View style={{ gap: 8 }}>
            <Text style={[ui.body, { color: colors.red, fontWeight: '700' }]}>Tout ton XP, tes badges et ta série seront effacés. Sûr ?</Text>
            <Text style={ui.muted}>Garde d’abord un code de sauvegarde : tu pourras tout retrouver si tu changes d’avis.</Text>
            <Button label="💾 Sauvegarder avant d’effacer" variant="secondary" onPress={() => void saveBackup('reset')} />
            {backup?.place === 'reset' && <BackupResult backup={backup} />}
            <View style={{ flexDirection: 'row', gap: 8 }}>
              <Button label="Annuler" variant="secondary" style={{ flex: 1 }} onPress={() => setConfirmReset(false)} />
              <Button
                label="Oui, effacer"
                color={colors.red}
                style={{ flex: 1 }}
                onPress={() => {
                  setBackup(null);
                  reset();
                }}
              />
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

interface BackupState {
  /** Bouton d'origine : réglages ou réinitialisation. */
  place: 'settings' | 'reset';
  code: string;
  /** Envoi non confirmé (Android, web, partage annulé ou impossible) ou code très long : afficher le code à copier. */
  showCode: boolean;
  /** La feuille de partage s'est ouverte (téléphone) : l'élève a peut-être déjà envoyé le message. */
  shareOpened: boolean;
  /** Code plus long que BACKUP_SAFE_LENGTH. */
  long: boolean;
}

/** Après « Sauvegarder » : confirmation du partage, ou code à copier quand le partage n'a pas abouti. */
function BackupResult({ backup }: { backup: BackupState }) {
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
  if (!backup.showCode) {
    return (
      <Text style={[ui.body, { color: colors.primaryDark, fontWeight: '700' }]} accessibilityLiveRegion="polite">
        ✓ Sauvegarde envoyée. Garde bien ce message : il te suffira de le coller dans « Restaurer une sauvegarde ».
      </Text>
    );
  }
  return (
    <View style={{ gap: 6 }} accessibilityLiveRegion="polite">
      <Text style={styles.settingTitle}>Copie ce code</Text>
      {backup.shareOpened && (
        <Text style={ui.body}>Si tu as bien envoyé le message, c’est bon. Sinon, copie ce code pour le garder.</Text>
      )}
      <Text style={ui.muted}>
        {backup.long
          ? 'Ta progression est bien remplie, bravo ! Le code est long : colle-le de préférence dans une note, un message risquerait de le couper.'
          : 'Sélectionne-le en entier, copie-le et colle-le dans une note ou un message à toi-même.'}
      </Text>
      <TextInput
        value={backup.code}
        // Champ en lecture seule : le texte reste sélectionnable et copiable (sur Android, un champ
        // non modifiable ne se sélectionne plus, d'où la valeur figée plutôt que readOnly).
        onChangeText={() => {}}
        readOnly={Platform.OS === 'web'}
        multiline
        selectTextOnFocus
        showSoftInputOnFocus={false}
        autoCorrect={false}
        spellCheck={false}
        accessibilityLabel="Code de sauvegarde"
        style={styles.code}
      />
    </View>
  );
}

const APPEARANCES: { value: ThemePreference; label: string }[] = [
  { value: 'auto', label: 'Automatique' },
  { value: 'light', label: 'Clair' },
  { value: 'dark', label: 'Sombre' },
];

/** Apparence : suit le téléphone par défaut, ou reste en clair ou en sombre. Rangée à part de la progression. */
function AppearanceCard() {
  const { preference, setPreference } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
  return (
    <Card style={{ gap: 10 }}>
      <Text style={styles.settingTitle}>Apparence</Text>
      <View style={{ flexDirection: 'row', gap: 8 }} accessibilityRole="radiogroup">
        {APPEARANCES.map((a) => {
          const on = preference === a.value;
          return (
            <Pressable
              key={a.value}
              onPress={() => setPreference(a.value)}
              accessibilityRole="radio"
              accessibilityState={{ checked: on, selected: on }}
              style={[styles.choice, on && styles.choiceOn]}
            >
              <Text style={[styles.choiceText, on && styles.choiceTextOn]}>{a.label}</Text>
            </Pressable>
          );
        })}
      </View>
      <Text style={ui.muted}>« Automatique » suit le réglage de ton téléphone. Le mode sombre est plus reposant le soir.</Text>
    </Card>
  );
}

/** Rappel quotidien : interrupteur et heure. La permission n'est demandée qu'à l'activation. */
function ReminderCard({ reminder, onChange }: { reminder: ReminderSettings | null; onChange: (reminder: ReminderSettings) => void }) {
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
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
          thumbColor={lightColors.card}
        />
      </View>
      <View style={styles.hours}>
        {REMINDER_HOURS.map((h) => {
          const on = enabled && hour === h;
          return (
            <Pressable
              key={h}
              disabled={busy}
              onPress={() => void apply(true, h)}
              accessibilityRole="button"
              accessibilityState={{ selected: on, disabled: busy }}
              style={[styles.hour, on && styles.choiceOn]}
            >
              <Text style={[styles.choiceText, on && styles.choiceTextOn]}>{formatHour(h)}</Text>
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

const makeStyles = (colors: Colors) =>
  StyleSheet.create({
    avatar: { minWidth: 72, minHeight: 72, borderRadius: 36, backgroundColor: colors.primary, alignItems: 'center', justifyContent: 'center' },
    avatarText: { fontSize: 32, fontWeight: '900', color: textOn(colors.primary) },
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
    choiceText: { fontWeight: '800', color: colors.text },
    choiceTextOn: { color: textOn(colors.primary) },
    hours: { flexDirection: 'row', flexWrap: 'wrap', gap: 8 },
    hour: { borderWidth: 2, borderColor: colors.border, borderRadius: radius.sm, paddingVertical: 10, paddingHorizontal: 12, backgroundColor: colors.card },
    badgeProgress: { alignSelf: 'stretch', gap: 3, marginTop: 4 },
    badgeProgressText: { fontSize: 11, fontWeight: '700', color: colors.goldText, textAlign: 'center' },
    code: {
      maxHeight: 140,
      borderWidth: 2,
      borderColor: colors.border,
      borderRadius: radius.sm,
      padding: 10,
      fontSize: 12,
      fontFamily: Platform.OS === 'ios' ? 'Menlo' : 'monospace',
      textAlignVertical: 'top',
      backgroundColor: colors.bg,
      color: colors.text,
    },
    link: { fontSize: 14, fontWeight: '700', color: colors.primaryDark, textDecorationLine: 'underline' },
  });
