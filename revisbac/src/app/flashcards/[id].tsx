import { router, Stack, useLocalSearchParams } from 'expo-router';
import { useState } from 'react';
import { Pressable, StyleSheet, Text, View } from 'react-native';

import { RewardModal } from '../../components/RewardModal';
import { Button, Pill, ProgressBar, Screen, useUi } from '../../components/ui';
import { getChapter } from '../../data/catalog';
import type { Flashcard } from '../../data/types';
import { dayKey } from '../../lib/dates';
import type { Reward } from '../../lib/gamification';
import { dueCardsText } from '../../lib/selectors';
import { CARD_BOX_DAYS, cardKey, chapterDueCount, deckOrder, isCardDue, type CardEntry } from '../../lib/srs';
import { useProgress } from '../../state/progress';
import { useStyles, useTheme } from '../../state/theme';
import { accentText, radius, shadow, textOn, type Colors } from '../../theme';

export default function FlashcardsScreen() {
  const { id } = useLocalSearchParams<{ id: string }>();
  const ui = useUi();
  const ref = getChapter(id);
  if (!ref) return <Text style={[ui.body, { padding: 20 }]}>Chapitre introuvable.</Text>;
  return <Deck key={id} chapterId={id} cards={ref.chapter.flashcards} color={ref.subject.color} subjectId={ref.subject.id} title={ref.chapter.title} />;
}

/** Où en est une carte (répétition espacée), en clair pour l'élève. */
function cardStatus(entry: CardEntry | undefined, today: string): string {
  if (!entry) return '✨ Nouvelle carte';
  if (isCardDue(entry, today)) return `🔄 À revoir · boîte ${entry.box}/5`;
  return `Boîte ${entry.box}/5`;
}

/** Intervalles de révision lisibles : « 1, 3, 7, 14 puis 30 jours ». */
const BOX_DAYS_TEXT = (() => {
  const days = Object.values(CARD_BOX_DAYS);
  return `${days.slice(0, -1).join(', ')} puis ${days[days.length - 1]} jours`;
})();

/**
 * Paquet de cartes : on retourne la carte, on dit si on savait ; les cartes ratées reviennent à la fin.
 * Répétition espacée : les cartes à revoir passent d'abord, chaque réponse fait monter ou redescendre la carte de boîte.
 */
function Deck({ chapterId, cards, color, subjectId, title }: { chapterId: string; cards: Flashcard[]; color: string; subjectId: string; title: string }) {
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
  const { state, finishFlashcards, reviewCard } = useProgress();
  // Jour du début du paquet : l'ordre et le nombre de cartes à revoir sont figés à l'ouverture.
  const [today] = useState(() => dayKey());
  const [queue, setQueue] = useState(() => deckOrder(chapterId, cards, state.cards, today));
  const [dueCount, setDueCount] = useState(() => chapterDueCount(chapterId, cards, state.cards, today));
  const [flipped, setFlipped] = useState(false);
  const [known, setKnown] = useState(0);
  const [reward, setReward] = useState<Reward | null>(null);
  const done = queue.length === 0;
  const status = done ? '' : cardStatus(state.cards[cardKey(chapterId, queue[0])], today);
  const due = dueCardsText(dueCount);

  const answer = (knew: boolean) => {
    const [card, ...rest] = queue;
    setFlipped(false);
    // Boîte suivante si l'élève savait, retour en boîte 1 sinon (pas d'XP par carte).
    reviewCard(cardKey(chapterId, card), knew);
    if (knew) {
      setKnown((k) => k + 1);
      setQueue(rest);
      if (rest.length === 0) setReward(finishFlashcards(subjectId, chapterId));
    } else {
      setQueue([...rest, card]);
    }
  };

  const restart = () => {
    setQueue(deckOrder(chapterId, cards, state.cards, today));
    setDueCount(chapterDueCount(chapterId, cards, state.cards, today));
    setKnown(0);
  };

  return (
    <Screen edges={['bottom']}>
      <Stack.Screen options={{ title }} />
      <ProgressBar value={known / cards.length} color={color} height={10} accessibilityLabel="Cartes maîtrisées" />
      <Text style={ui.muted}>
        {known}/{cards.length} cartes maîtrisées · les cartes ratées reviennent à la fin
      </Text>
      {due && !done && <Text style={[ui.body, { fontWeight: '700' }]}>{due} : elles passent en premier.</Text>}

      {!done ? (
        <>
          <Pressable
            onPress={() => setFlipped((f) => !f)}
            accessibilityRole="button"
            accessibilityLabel={`${status}. ${flipped ? 'Réponse' : 'Question'} : ${flipped ? queue[0].back : queue[0].front}`}
            accessibilityHint="Retourne la carte"
            style={[styles.card, { borderColor: color }, flipped && { backgroundColor: color }]}
          >
            {!flipped && <Pill label={status} color={colors.muted} />}
            <Text style={[styles.side, { color: flipped ? textOn(color) : accentText(color, colors) }]}>{flipped ? 'RÉPONSE' : 'QUESTION'}</Text>
            <Text style={[styles.cardText, flipped && { color: textOn(color) }]}>{flipped ? queue[0].back : queue[0].front}</Text>
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
          <Text style={[ui.muted, { textAlign: 'center' }]}>
            Tes cartes reviendront de plus en plus espacées (après {BOX_DAYS_TEXT}) tant que tu les sais : c’est comme ça qu’on retient pour de bon.
          </Text>
          <Button label="✅ Passer au quiz" color={color} onPress={() => router.replace({ pathname: '/quiz', params: { mode: 'chapter', id: chapterId } })} />
          <Button label="Recommencer" variant="secondary" color={color} onPress={restart} />
        </View>
      )}
      <RewardModal reward={reward} onClose={() => setReward(null)} haptics={state.profile?.haptics ?? true} />
    </Screen>
  );
}

const makeStyles = (colors: Colors) =>
  StyleSheet.create({
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
