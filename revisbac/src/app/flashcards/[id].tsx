import { router, Stack, useLocalSearchParams } from 'expo-router';
import { useState } from 'react';
import { Pressable, StyleSheet, Text, View } from 'react-native';

import { RewardModal } from '../../components/RewardModal';
import { Button, ProgressBar, Screen, styles as ui } from '../../components/ui';
import { getChapter } from '../../data/catalog';
import type { Flashcard } from '../../data/types';
import type { Reward } from '../../lib/gamification';
import { shuffle } from '../../lib/random';
import { useProgress } from '../../state/progress';
import { colors, radius, shadow } from '../../theme';

export default function FlashcardsScreen() {
  const { id } = useLocalSearchParams<{ id: string }>();
  const ref = getChapter(id);
  if (!ref) return <Text style={{ padding: 20 }}>Chapitre introuvable.</Text>;
  return <Deck key={id} chapterId={id} cards={ref.chapter.flashcards} color={ref.subject.color} subjectId={ref.subject.id} title={ref.chapter.title} />;
}

/** Paquet de cartes : on retourne la carte, on dit si on savait ; les cartes ratées reviennent à la fin. */
function Deck({ chapterId, cards, color, subjectId, title }: { chapterId: string; cards: Flashcard[]; color: string; subjectId: string; title: string }) {
  const { state, finishFlashcards } = useProgress();
  const [queue, setQueue] = useState(() => shuffle(cards));
  const [flipped, setFlipped] = useState(false);
  const [known, setKnown] = useState(0);
  const [reward, setReward] = useState<Reward | null>(null);
  const done = queue.length === 0;

  const answer = (knew: boolean) => {
    const [card, ...rest] = queue;
    setFlipped(false);
    if (knew) {
      setKnown((k) => k + 1);
      setQueue(rest);
      if (rest.length === 0) setReward(finishFlashcards(subjectId, chapterId));
    } else {
      setQueue([...rest, card]);
    }
  };

  return (
    <Screen edges={['bottom']}>
      <Stack.Screen options={{ title }} />
      <ProgressBar value={known / cards.length} color={color} height={10} />
      <Text style={ui.muted}>
        {known}/{cards.length} cartes maîtrisées · les cartes ratées reviennent à la fin
      </Text>

      {!done ? (
        <>
          <Pressable onPress={() => setFlipped((f) => !f)} style={[styles.card, { borderColor: color }, flipped && { backgroundColor: color }]}>
            <Text style={[styles.side, { color: flipped ? '#ffffffb0' : color }]}>{flipped ? 'RÉPONSE' : 'QUESTION'}</Text>
            <Text style={[styles.cardText, flipped && { color: '#fff' }]}>{flipped ? queue[0].back : queue[0].front}</Text>
            {!flipped && <Text style={ui.muted}>Touche la carte pour la retourner</Text>}
          </Pressable>
          {flipped ? (
            <View style={{ flexDirection: 'row', gap: 10 }}>
              <Button label="✗ À revoir" variant="secondary" color={colors.red} style={{ flex: 1 }} onPress={() => answer(false)} />
              <Button label="✓ Je savais" style={{ flex: 1 }} onPress={() => answer(true)} />
            </View>
          ) : (
            <Button label="Voir la réponse" color={color} onPress={() => setFlipped(true)} />
          )}
        </>
      ) : (
        <View style={{ paddingVertical: 32, gap: 14 }}>
          <Text style={{ fontSize: 60, textAlign: 'center' }}>🎉</Text>
          <Text style={[ui.h2, { textAlign: 'center' }]}>Toutes les cartes sont maîtrisées !</Text>
          <Button label="✅ Passer au quiz" color={color} onPress={() => router.replace({ pathname: '/quiz', params: { mode: 'chapter', id: chapterId } })} />
          <Button label="Recommencer" variant="secondary" color={color} onPress={() => { setQueue(shuffle(cards)); setKnown(0); }} />
        </View>
      )}
      <RewardModal reward={reward} onClose={() => setReward(null)} haptics={state.profile?.haptics ?? true} />
    </Screen>
  );
}

const styles = StyleSheet.create({
  card: {
    // Hauteur minimale de carte ; un texte long l'agrandit et l'écran défile.
    minHeight: 300,
    borderWidth: 3,
    borderRadius: radius.lg,
    backgroundColor: colors.card,
    padding: 24,
    alignItems: 'center',
    justifyContent: 'center',
    gap: 16,
    ...shadow,
  },
  side: { fontSize: 12, fontWeight: '900', letterSpacing: 1 },
  cardText: { fontSize: 22, fontWeight: '700', color: colors.text, textAlign: 'center', lineHeight: 32 },
});
