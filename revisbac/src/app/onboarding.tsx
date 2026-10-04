import { router } from 'expo-router';
import { useEffect, useState } from 'react';
import { ActivityIndicator, BackHandler, KeyboardAvoidingView, Platform, Pressable, StyleSheet, Text, TextInput, View } from 'react-native';

import { SubjectPicker } from '../components/SubjectPicker';
import { Button, Screen, styles as ui } from '../components/ui';
import { tracks } from '../data/catalog';
import type { TrackId } from '../data/types';
import { createProfile, DAILY_GOAL_OPTIONS, type ReminderSettings } from '../lib/gamification';
import { DEFAULT_REMINDER_HOUR, formatHour, REMINDER_HOURS, requestReminderPermission } from '../lib/reminders';
import { useProgress } from '../state/progress';
import { colors, radius } from '../theme';

const GOAL_LABELS: Record<number, string> = { 30: 'Tranquille', 50: 'Régulier', 100: 'Sérieux', 150: 'Intense' };
// Pas de notifications sur le web : l'étape du rappel n'y est pas proposée.
const WITH_REMINDER = Platform.OS !== 'web';
const STEPS = WITH_REMINDER ? 5 : 4;

/** Quitte l'onboarding sans empiler une 2e instance des onglets quand il a été ouvert depuis le Profil. */
function leave() {
  if (router.canGoBack()) router.back();
  else router.replace('/');
}

export default function Onboarding() {
  const { loaded } = useProgress();
  // Le formulaire s'initialise depuis le profil : on attend qu'il soit chargé (rechargement web, lien profond).
  if (!loaded) return <ActivityIndicator color={colors.primary} style={{ flex: 1, backgroundColor: colors.bg }} />;
  return <OnboardingForm />;
}

function OnboardingForm() {
  const { state, setProfile } = useProgress();
  const editing = !!state.profile;
  const [step, setStep] = useState(0);
  const [name, setName] = useState(state.profile?.name ?? '');
  const [track, setTrack] = useState<TrackId | null>(state.profile?.track ?? null);
  const [hidden, setHidden] = useState<string[]>(state.profile?.hiddenSubjects ?? []);
  const [goal, setGoal] = useState(state.profile?.dailyGoal ?? 50);
  const previousReminder = state.profile?.reminder ?? null;
  // Heure du rappel, null = pas de rappel. 19 h est proposé par défaut.
  const [reminderHour, setReminderHour] = useState<number | null>(
    previousReminder ? (previousReminder.enabled ? previousReminder.hour : null) : DEFAULT_REMINDER_HOUR,
  );
  const [saving, setSaving] = useState(false);

  // Bouton retour d'Android : revient à l'étape précédente au lieu de quitter l'onboarding.
  useEffect(() => {
    if (step === 0) return;
    const sub = BackHandler.addEventListener('hardwareBackPress', () => {
      setStep((s) => Math.max(0, s - 1));
      return true;
    });
    return () => sub.remove();
  }, [step]);

  const chooseTrack = (id: TrackId) => {
    setTrack(id);
    // Les matières masquées dépendent de l'examen : remises à zéro si l'élève en change.
    setHidden(id === state.profile?.track ? state.profile.hiddenSubjects : []);
  };

  const finish = async () => {
    if (!track || saving) return;
    setSaving(true);
    const info = tracks.find((t) => t.id === track)!;
    const sameTrack = state.profile?.track === track;
    let reminder: ReminderSettings | null = previousReminder;
    if (WITH_REMINDER) {
      if (reminderHour === null) {
        reminder = { enabled: false, hour: previousReminder?.hour ?? DEFAULT_REMINDER_HOUR, minute: previousReminder?.minute ?? 0 };
      } else {
        // La permission n'est demandée qu'une fois l'heure choisie ; un refus désactive simplement le rappel.
        const granted = await requestReminderPermission();
        reminder = { enabled: granted, hour: reminderHour, minute: 0 };
      }
    }
    // Les autres réglages existants (vibrations) sont gardés.
    setProfile(
      createProfile({
        ...state.profile,
        name: name.trim() || 'Champion',
        track,
        dailyGoal: goal,
        examDate: sameTrack ? state.profile!.examDate : info.defaultExamDate,
        hiddenSubjects: hidden,
        reminder,
      }),
    );
    leave();
  };

  return (
    <KeyboardAvoidingView style={styles.flex} behavior={Platform.OS === 'ios' ? 'padding' : 'height'}>
      <Screen edges={['top', 'bottom']}>
        <View style={styles.dots}>
          {Array.from({ length: STEPS }, (_, i) => (
            <View key={i} style={[styles.dot, i <= step && { backgroundColor: colors.primary }]} />
          ))}
        </View>

        {step === 0 && (
          <View style={{ gap: 16 }}>
            <Text style={styles.hero}>🇸🇳</Text>
            <Text style={[ui.h1, { textAlign: 'center' }]}>Bienvenue sur RéviBac</Text>
            <Text style={[ui.body, { textAlign: 'center', color: colors.muted }]}>
              Des fiches courtes, des quiz et des défis quotidiens pour réussir ton BFM ou ton Bac. Quelques minutes par jour suffisent !
            </Text>
            <Text style={styles.label}>Comment t’appelles-tu ?</Text>
            <TextInput
              value={name}
              onChangeText={setName}
              placeholder="Ton prénom"
              placeholderTextColor={colors.muted}
              style={styles.input}
              maxLength={30}
              returnKeyType="next"
              onSubmitEditing={() => setStep(1)}
            />
            <Button label="Continuer" onPress={() => setStep(1)} />
            {editing && <Button label="Annuler" variant="ghost" onPress={leave} />}
          </View>
        )}

        {step === 1 && (
          <View style={{ gap: 12 }}>
            <Text style={ui.h1}>Quel examen prépares-tu ?</Text>
            {tracks.map((t) => (
              <Pressable
                key={t.id}
                onPress={() => chooseTrack(t.id)}
                style={[styles.option, track === t.id && { borderColor: colors.primary, backgroundColor: colors.primarySoft }]}
              >
                <Text style={{ fontSize: 30 }}>{t.emoji}</Text>
                <View style={{ flex: 1 }}>
                  <Text style={styles.optionTitle}>{t.label}</Text>
                  <Text style={ui.muted}>{t.description}</Text>
                </View>
              </Pressable>
            ))}
            <Button label="Continuer" disabled={!track} onPress={() => setStep(2)} />
            <Button label="Retour" variant="ghost" onPress={() => setStep(0)} />
          </View>
        )}

        {step === 2 && track && (
          <View style={{ gap: 12 }}>
            <Text style={ui.h1}>Tes matières</Text>
            <SubjectPicker track={track} hidden={hidden} onChange={setHidden} />
            <Button label="Continuer" onPress={() => setStep(3)} />
            <Button label="Retour" variant="ghost" onPress={() => setStep(1)} />
          </View>
        )}

        {step === 3 && (
          <View style={{ gap: 12 }}>
            <Text style={ui.h1}>Ton objectif quotidien</Text>
            <Text style={[ui.body, { color: colors.muted }]}>
              Atteins-le chaque jour pour gagner un bonus et faire grandir ta série 🔥. Tu pourras le changer plus tard.
            </Text>
            {DAILY_GOAL_OPTIONS.map((g) => (
              <Pressable
                key={g}
                onPress={() => setGoal(g)}
                style={[styles.option, goal === g && { borderColor: colors.primary, backgroundColor: colors.primarySoft }]}
              >
                <View style={{ flex: 1 }}>
                  <Text style={styles.optionTitle}>{GOAL_LABELS[g]}</Text>
                  <Text style={ui.muted}>≈ {Math.round(g / 10)} bonnes réponses par jour</Text>
                </View>
                <Text style={styles.goalXp}>{g} XP</Text>
              </Pressable>
            ))}
            {WITH_REMINDER ? (
              <Button label="Continuer" onPress={() => setStep(4)} />
            ) : (
              <Button label="C’est parti ! 🚀" disabled={saving} onPress={finish} />
            )}
            <Button label="Retour" variant="ghost" onPress={() => setStep(2)} />
          </View>
        )}

        {step === 4 && (
          <View style={{ gap: 12 }}>
            <Text style={ui.h1}>Un petit rappel chaque jour ?</Text>
            <Text style={[ui.body, { color: colors.muted }]}>
              Une seule notification par jour, à l’heure de ton choix. Tu pourras la changer dans Profil.
            </Text>
            <View style={styles.hours}>
              {REMINDER_HOURS.map((h) => (
                <Pressable key={h} onPress={() => setReminderHour(h)} style={[styles.hour, reminderHour === h && styles.hourOn]}>
                  <Text style={[styles.hourText, reminderHour === h && { color: '#fff' }]}>{formatHour(h)}</Text>
                </Pressable>
              ))}
              <Pressable onPress={() => setReminderHour(null)} style={[styles.hour, reminderHour === null && styles.hourOn]}>
                <Text style={[styles.hourText, reminderHour === null && { color: '#fff' }]}>Pas de rappel</Text>
              </Pressable>
            </View>
            <Button label="C’est parti ! 🚀" disabled={saving} onPress={finish} />
            <Button label="Retour" variant="ghost" onPress={() => setStep(3)} />
          </View>
        )}
      </Screen>
    </KeyboardAvoidingView>
  );
}

const styles = StyleSheet.create({
  flex: { flex: 1, backgroundColor: colors.bg },
  dots: { flexDirection: 'row', gap: 6, justifyContent: 'center', marginVertical: 8 },
  dot: { width: 28, height: 6, borderRadius: 3, backgroundColor: colors.border },
  hero: { fontSize: 64, textAlign: 'center', marginTop: 24 },
  label: { fontSize: 15, fontWeight: '700', color: colors.text, marginTop: 8 },
  input: {
    borderWidth: 2,
    borderColor: colors.border,
    borderRadius: radius.md,
    paddingHorizontal: 14,
    paddingVertical: 12,
    fontSize: 17,
    backgroundColor: colors.card,
    color: colors.text,
  },
  option: {
    flexDirection: 'row',
    alignItems: 'center',
    gap: 14,
    borderWidth: 2,
    borderColor: colors.border,
    borderRadius: radius.md,
    padding: 14,
    backgroundColor: colors.card,
  },
  optionTitle: { fontSize: 17, fontWeight: '800', color: colors.text },
  goalXp: { fontSize: 16, fontWeight: '800', color: colors.primary },
  hours: { flexDirection: 'row', flexWrap: 'wrap', gap: 8 },
  hour: {
    borderWidth: 2,
    borderColor: colors.border,
    borderRadius: radius.sm,
    paddingVertical: 10,
    paddingHorizontal: 14,
    backgroundColor: colors.card,
  },
  hourOn: { backgroundColor: colors.primary, borderColor: colors.primary },
  hourText: { fontWeight: '800', color: colors.text },
});
