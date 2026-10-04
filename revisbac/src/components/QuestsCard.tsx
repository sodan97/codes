import { router } from 'expo-router';
import { useState } from 'react';
import { Pressable, StyleSheet, Text, View } from 'react-native';

import { dayKey } from '../lib/dates';
import { chestState, todayStats, XP, type Reward } from '../lib/gamification';
import { questLink, questProgressLabel, questTierLabel, questTitle, type Quest } from '../lib/quests';
import { useProgress } from '../state/progress';
import { useStyles, useTheme } from '../state/theme';
import { radius, textOn, type Colors } from '../theme';
import { RewardModal } from './RewardModal';
import { Button, Card, MAX_FONT_SCALE, ProgressBar, useUi } from './ui';

/**
 * Carte « Quêtes du jour » (accueil et onglet Défis) : 3 missions, chacune mène à l'écran concerné.
 * Quand les 3 sont faites, « Ouvrir le coffre » donne la récompense du jour (RewardModal).
 * Les quêtes non faites disparaissent à minuit : rien n'est signalé comme raté.
 */
export function QuestsCard() {
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
  const { state, openChest } = useProgress();
  const [reward, setReward] = useState<Reward | null>(null);
  const today = dayKey();
  const quests = todayStats(state, today).quests;
  if (quests.length === 0) return null;
  const done = quests.filter((q) => q.done).length;
  const chest = chestState(state, today);

  const open = () => {
    const r = openChest();
    if (r) setReward(r);
  };

  return (
    <Card style={{ gap: 10 }}>
      <View style={[ui.row, { justifyContent: 'space-between' }]}>
        <Text style={styles.title} accessibilityRole="header">
          🗺️ Quêtes du jour
        </Text>
        <Text style={styles.count} maxFontSizeMultiplier={MAX_FONT_SCALE} accessibilityLabel={`${done} quête${done > 1 ? 's' : ''} faite${done > 1 ? 's' : ''} sur ${quests.length}`}>
          {done}/{quests.length}
        </Text>
      </View>
      <Text style={ui.muted}>+{XP.quest} XP par quête, et un coffre quand les {quests.length} sont faites.</Text>

      {quests.map((q) => (
        <QuestRow key={q.id} quest={q} onPress={() => router.push(questLink(q, state, today))} />
      ))}

      {chest === 'ready' ? (
        <Button label="🎁 Ouvrir le coffre" variant="gold" onPress={open} />
      ) : (
        <View style={[styles.chest, chest === 'opened' && { backgroundColor: colors.goldSoft }]}>
          <Text style={[styles.chestIcon, chest === 'locked' && { opacity: 0.45 }]} maxFontSizeMultiplier={MAX_FONT_SCALE}>
            🎁
          </Text>
          <Text style={[ui.muted, { flex: 1 }]}>
            {chest === 'opened' ? 'Coffre ouvert, bravo ! De nouvelles quêtes t’attendent demain.' : `Termine les ${quests.length} quêtes pour ouvrir le coffre du jour.`}
          </Text>
        </View>
      )}

      <RewardModal reward={reward} onClose={() => setReward(null)} haptics={state.profile?.haptics ?? true} title="Coffre du jour" icon="🎁" />
    </Card>
  );
}

/** Une quête : difficulté, intitulé, avancement. Une quête faite n'est plus un lien. */
function QuestRow({ quest, onPress }: { quest: Quest; onPress: () => void }) {
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
  const title = questTitle(quest);
  const tier = questTierLabel(quest.tier);
  const progress = questProgressLabel(quest);
  const content = (
    <>
      <View style={[styles.check, quest.done && { backgroundColor: colors.primary, borderColor: colors.primary }]}>
        <Text style={[styles.checkText, quest.done && styles.checkTextDone]} maxFontSizeMultiplier={MAX_FONT_SCALE}>
          {quest.done ? '✓' : ''}
        </Text>
      </View>
      <View style={{ flex: 1, gap: 4 }}>
        <Text style={[ui.muted, styles.tier]}>{tier}</Text>
        <Text style={[styles.questTitle, quest.done && styles.questDone]}>{title}</Text>
        {!quest.done && quest.target > 1 && <ProgressBar value={quest.progress / quest.target} height={6} />}
      </View>
      <Text style={[styles.progress, quest.done && { color: colors.primaryDark }]} maxFontSizeMultiplier={MAX_FONT_SCALE}>
        {progress}
      </Text>
      {!quest.done && <Text style={styles.chevron}>›</Text>}
    </>
  );

  if (quest.done) {
    return (
      <View style={[styles.row, styles.rowDone]} accessible accessibilityLabel={`Quête ${tier.toLowerCase()} accomplie : ${title}`}>
        {content}
      </View>
    );
  }
  return (
    <Pressable
      onPress={onPress}
      accessibilityRole="button"
      accessibilityLabel={`Quête ${tier.toLowerCase()} : ${title}, ${quest.progress} sur ${quest.target}`}
      accessibilityHint="Ouvre l’écran pour avancer cette quête"
      style={({ pressed }) => [styles.row, pressed && { opacity: 0.8 }]}
    >
      {content}
    </Pressable>
  );
}

const makeStyles = (colors: Colors) =>
  StyleSheet.create({
    title: { fontSize: 16, fontWeight: '800', color: colors.text },
    count: { fontSize: 15, fontWeight: '900', color: colors.primaryDark },
    row: {
      flexDirection: 'row',
      alignItems: 'center',
      gap: 10,
      padding: 10,
      borderRadius: radius.sm,
      borderWidth: 1,
      borderColor: colors.border,
      minHeight: 48,
    },
    rowDone: { backgroundColor: colors.primarySoft, borderColor: colors.primarySoft },
    check: { minWidth: 26, minHeight: 26, borderRadius: 13, borderWidth: 2, borderColor: colors.border, alignItems: 'center', justifyContent: 'center' },
    checkText: { fontSize: 14, fontWeight: '900', color: colors.text },
    // Coche blanche ou sombre selon le vert du thème.
    checkTextDone: { color: textOn(colors.primary) },
    tier: { fontSize: 11, fontWeight: '700', textTransform: 'uppercase' },
    questTitle: { fontSize: 15, fontWeight: '700', color: colors.text },
    // Quête faite : posée sur primarySoft, le gris atténué y passerait sous 4,5:1.
    questDone: { color: colors.primaryDark },
    progress: { fontSize: 14, fontWeight: '800', color: colors.text, fontVariant: ['tabular-nums'] },
    chevron: { fontSize: 24, color: colors.muted },
    chest: { flexDirection: 'row', alignItems: 'center', gap: 10, padding: 10, borderRadius: radius.sm, borderWidth: 1, borderColor: colors.border, borderStyle: 'dashed' },
    chestIcon: { fontSize: 24 },
  });
