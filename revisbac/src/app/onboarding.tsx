import { router } from 'expo-router';
import { useState } from 'react';
import { Pressable, StyleSheet, Text, TextInput, View } from 'react-native';

import { Button, Screen, styles as ui } from '../components/ui';
import { tracks } from '../data/catalog';
import type { TrackId } from '../data/types';
import { DAILY_GOAL_OPTIONS } from '../lib/gamification';
import { useProgress } from '../state/progress';
import { colors, radius } from '../theme';

const GOAL_LABELS: Record<number, string> = { 30: 'Tranquille', 50: 'Régulier', 100: 'Sérieux', 150: 'Intense' };

export default function Onboarding() {
  const { state, setProfile } = useProgress();
  const [step, setStep] = useState(0);
  const [name, setName] = useState(state.profile?.name ?? '');
  const [track, setTrack] = useState<TrackId | null>(state.profile?.track ?? null);
  const [goal, setGoal] = useState(state.profile?.dailyGoal ?? 50);

  const finish = () => {
    if (!track) return;
    const info = tracks.find((t) => t.id === track)!;
    setProfile({
      name: name.trim() || 'Champion',
      track,
      dailyGoal: goal,
      examDate: state.profile?.track === track ? state.profile.examDate : info.defaultExamDate,
    });
    router.replace('/');
  };

  return (
    <Screen edges={['top', 'bottom']}>
      <View style={styles.dots}>
        {[0, 1, 2].map((i) => (
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
        </View>
      )}

      {step === 1 && (
        <View style={{ gap: 12 }}>
          <Text style={ui.h1}>Quel examen prépares-tu ?</Text>
          {tracks.map((t) => (
            <Pressable
              key={t.id}
              onPress={() => setTrack(t.id)}
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

      {step === 2 && (
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
          <Button label="C’est parti ! 🚀" onPress={finish} />
          <Button label="Retour" variant="ghost" onPress={() => setStep(1)} />
        </View>
      )}
    </Screen>
  );
}

const styles = StyleSheet.create({
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
});
