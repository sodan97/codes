import { router, Stack, useLocalSearchParams } from 'expo-router';
import { useEffect, useState } from 'react';
import { StyleSheet, Text, View } from 'react-native';

import { FicheBlockView } from '../../components/FicheBlockView';
import { RewardModal } from '../../components/RewardModal';
import { Button, Card, Screen, useUi } from '../../components/ui';
import { getChapter } from '../../data/catalog';
import { XP, type Reward } from '../../lib/gamification';
import { useProgress } from '../../state/progress';
import { saveLastFiche } from '../../state/session';
import { useStyles, useTheme } from '../../state/theme';
import { accentText, radius, type Colors } from '../../theme';

export default function FicheScreen() {
  const { id } = useLocalSearchParams<{ id: string }>();
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
  const { state, markFicheRead } = useProgress();
  const [reward, setReward] = useState<Reward | null>(null);
  const ref = getChapter(id);
  const chapterId = ref?.chapter.id;

  // Dernière fiche ouverte : l'accueil propose « Continuer ta fiche » tant qu'elle n'est pas lue.
  useEffect(() => {
    if (chapterId) saveLastFiche(chapterId);
  }, [chapterId]);

  if (!ref) return <Text style={[ui.body, { padding: 20 }]}>Fiche introuvable.</Text>;
  const { chapter, subject } = ref;
  // Couleur de la matière pour le texte (lisible dans les deux thèmes).
  const accent = accentText(subject.color, colors);
  const read = !!state.fichesRead[chapter.id];
  const index = subject.chapters.findIndex((c) => c.id === chapter.id);
  const next = subject.chapters[index + 1];

  return (
    <Screen edges={['bottom']}>
      <Stack.Screen options={{ title: subject.name }} />
      <View style={{ gap: 4 }}>
        <Text style={[styles.kicker, { color: accent }]}>
          {subject.icon} Chapitre {index + 1}
        </Text>
        <Text style={ui.h1}>{chapter.title}</Text>
        <Text style={ui.muted}>{chapter.summary}</Text>
      </View>

      <View style={[styles.essentials, { borderColor: subject.color }]}>
        <Text style={[styles.essentialsTitle, { color: accent }]}>⏱️ L’essentiel en 30 secondes</Text>
        {chapter.essentials.map((e, i) => (
          <View key={i} style={{ flexDirection: 'row', gap: 8 }}>
            <Text style={{ color: accent, fontWeight: '900' }}>{i + 1}.</Text>
            <Text style={[ui.body, { flex: 1, fontWeight: '600' }]}>{e}</Text>
          </View>
        ))}
      </View>

      {chapter.sections.map((section, i) => (
        <Card key={i} style={{ gap: 12 }}>
          <Text style={styles.sectionTitle}>{section.title}</Text>
          {section.blocks.map((block, j) => (
            <FicheBlockView key={j} block={block} color={subject.color} />
          ))}
        </Card>
      ))}

      <Card style={{ gap: 10 }}>
        {read ? (
          <Text style={[ui.body, { textAlign: 'center', color: colors.primaryDark, fontWeight: '700' }]}>✓ Fiche lue — teste-toi maintenant !</Text>
        ) : (
          <Button label={`J’ai lu cette fiche ✓  (+${XP.ficheRead} XP)`} onPress={() => setReward(markFicheRead(chapter.id, subject.id))} />
        )}
        <View style={{ flexDirection: 'row', gap: 10 }}>
          <Button
            label="🃏 Flashcards"
            variant="secondary"
            color={subject.color}
            style={{ flex: 1 }}
            onPress={() => router.push({ pathname: '/flashcards/[id]', params: { id: chapter.id } })}
          />
          <Button
            label="✅ Quiz"
            color={subject.color}
            style={{ flex: 1 }}
            onPress={() => router.push({ pathname: '/quiz', params: { mode: 'chapter', id: chapter.id } })}
          />
        </View>
        {next && (
          <Button label={`Chapitre suivant : ${next.title} ›`} variant="ghost" color={colors.muted} onPress={() => router.replace({ pathname: '/fiche/[id]', params: { id: next.id } })} />
        )}
      </Card>

      <RewardModal reward={reward} onClose={() => setReward(null)} haptics={state.profile?.haptics ?? true} />
    </Screen>
  );
}

const makeStyles = (colors: Colors) =>
  StyleSheet.create({
    kicker: { fontSize: 13, fontWeight: '800', textTransform: 'uppercase', letterSpacing: 0.6 },
    essentials: { borderWidth: 2, borderRadius: radius.md, padding: 14, gap: 8, backgroundColor: colors.card },
    essentialsTitle: { fontSize: 15, fontWeight: '900' },
    sectionTitle: { fontSize: 18, fontWeight: '800', color: colors.text },
  });
