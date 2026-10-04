import { useEffect, useState } from 'react';
import { Animated, Modal, StyleSheet, Text, View } from 'react-native';

import { celebrate } from '../lib/feedback';
import { levelInfo, levelThreshold, type Badge, type Reward } from '../lib/gamification';
import { useStyles } from '../state/theme';
import { nativeDriver, radius, type Colors } from '../theme';
import { Button, MAX_FONT_SCALE, useReduceMotion } from './ui';

/**
 * Fenêtre de récompense affichée après une activité : XP, bonus, quêtes, badges, niveau.
 * `haptics` : réglage « Vibrations » du profil (vibration de célébration à la montée de niveau).
 * `icon` et `title` : en-tête facultatif (ex. 🎁 « Coffre du jour »).
 */
export function RewardModal({
  reward,
  onClose,
  haptics = true,
  icon,
  title,
}: {
  reward: Reward | null;
  onClose: () => void;
  haptics?: boolean;
  icon?: string;
  title?: string;
}) {
  const styles = useStyles(makeStyles);
  return (
    <Modal visible={!!reward} transparent animationType="fade" onRequestClose={onClose}>
      <View style={styles.backdrop}>{reward && <RewardContent reward={reward} onClose={onClose} haptics={haptics} icon={icon} title={title} />}</View>
    </Modal>
  );
}

/** Messages des quêtes et du coffre (voir gamification.applyGain), mis en avant dans un encadré. */
const isQuestMessage = (m: string) => m.startsWith('🗺️') || m.startsWith('🎁');

/** Contenu monté à chaque nouvelle récompense : les animations repartent de zéro. */
function RewardContent({
  reward,
  onClose,
  haptics,
  icon,
  title,
}: {
  reward: Reward;
  onClose: () => void;
  haptics: boolean;
  icon?: string;
  title?: string;
}) {
  const styles = useStyles(makeStyles);
  const reduceMotion = useReduceMotion();
  const [animatedXp, setAnimatedXp] = useState(0);
  const animateXp = !reduceMotion && reward.xp > 0;
  const shownXp = animateXp ? animatedXp : reward.xp;
  // « +0 XP » n'apporte rien quand un message dit ce qui a été gagné (ex. un gel de série dans le coffre).
  const showXp = reward.xp > 0 || reward.messages.length === 0;
  const quests = reward.messages.filter(isQuestMessage);
  const others = reward.messages.filter((m) => !isQuestMessage(m));

  // Le compteur « +XP » défile de 0 à xp en 800 ms.
  useEffect(() => {
    if (!animateXp) return;
    const value = new Animated.Value(0);
    const id = value.addListener(({ value: v }) => setAnimatedXp(Math.round(v)));
    const animation = Animated.timing(value, { toValue: reward.xp, duration: 800, useNativeDriver: false });
    animation.start();
    return () => {
      animation.stop();
      value.removeListener(id);
    };
  }, [reward.xp, animateXp]);

  useEffect(() => {
    if (reward.levelUp != null) celebrate(haptics);
  }, [reward.levelUp, haptics]);

  return (
    <View style={styles.box}>
      {reward.levelUp != null ? <Text style={{ fontSize: 56 }}>🎉</Text> : icon ? <Text style={{ fontSize: 56 }}>{icon}</Text> : null}
      {title && (
        <Text style={styles.title} accessibilityRole="header">
          {title}
        </Text>
      )}
      {showXp && (
        <Text style={styles.big} maxFontSizeMultiplier={MAX_FONT_SCALE} accessibilityLabel={`Plus ${reward.xp} XP`}>
          +{shownXp} XP
        </Text>
      )}
      {reward.levelUp != null && (
        <View style={styles.level}>
          <Text style={styles.levelText}>
            ⬆️ Niveau {reward.levelUp} : {levelInfo(levelThreshold(reward.levelUp)).title}
          </Text>
        </View>
      )}
      {others.map((m, i) => (
        <Text key={i} style={styles.msg}>
          {m}
        </Text>
      ))}
      {quests.length > 0 && (
        <View style={styles.quests}>
          {quests.map((m, i) => (
            <Text key={i} style={styles.questMsg}>
              {m}
            </Text>
          ))}
        </View>
      )}
      {reward.newBadges.length > 0 && (
        <View style={{ gap: 8, marginTop: 8, alignSelf: 'stretch' }}>
          <Text style={styles.badgeTitle}>Nouveau badge !</Text>
          {reward.newBadges.map((b, i) => (
            <BadgeRow key={b.id} badge={b} delay={i * 150} reduceMotion={reduceMotion} />
          ))}
        </View>
      )}
      <Button label="Super !" onPress={onClose} style={{ alignSelf: 'stretch', marginTop: 12 }} />
    </View>
  );
}

/** Badge qui apparaît avec un petit ressort (décalé de 150 ms d'un badge à l'autre). */
function BadgeRow({ badge, delay, reduceMotion }: { badge: Badge; delay: number; reduceMotion: boolean }) {
  const styles = useStyles(makeStyles);
  const [scale] = useState(() => new Animated.Value(reduceMotion ? 1 : 0));
  useEffect(() => {
    if (reduceMotion) {
      scale.setValue(1);
      return;
    }
    const animation = Animated.spring(scale, { toValue: 1, delay, friction: 5, useNativeDriver: nativeDriver });
    animation.start();
    return () => animation.stop();
  }, [scale, delay, reduceMotion]);

  return (
    <Animated.View style={[styles.badge, { transform: [{ scale }] }]}>
      <Text style={{ fontSize: 28 }}>{badge.icon}</Text>
      <View style={{ flex: 1 }}>
        <Text style={styles.badgeName}>{badge.name}</Text>
        <Text style={styles.badgeDesc}>{badge.description}</Text>
      </View>
    </Animated.View>
  );
}

const makeStyles = (colors: Colors) =>
  StyleSheet.create({
    backdrop: { flex: 1, backgroundColor: colors.backdrop, alignItems: 'center', justifyContent: 'center', padding: 24 },
    box: { backgroundColor: colors.card, borderRadius: radius.lg, padding: 22, alignItems: 'center', gap: 6, width: '100%', maxWidth: 420 },
    big: { fontSize: 40, fontWeight: '900', color: colors.primary, fontVariant: ['tabular-nums'] },
    title: { fontSize: 20, fontWeight: '800', color: colors.text, textAlign: 'center' },
    msg: { fontSize: 15, color: colors.text, textAlign: 'center' },
    quests: { alignSelf: 'stretch', gap: 6, marginTop: 6, backgroundColor: colors.primarySoft, borderRadius: radius.sm, padding: 10 },
    questMsg: { fontSize: 14, fontWeight: '700', color: colors.text, textAlign: 'center' },
    level: { backgroundColor: colors.goldSoft, borderRadius: radius.sm, paddingVertical: 6, paddingHorizontal: 12, marginBottom: 4 },
    levelText: { fontWeight: '800', color: colors.text },
    badgeTitle: { fontWeight: '800', color: colors.goldText, textAlign: 'center', fontSize: 16 },
    badgeName: { fontWeight: '800', color: colors.text },
    badgeDesc: { color: colors.muted, fontSize: 13 },
    badge: { flexDirection: 'row', gap: 12, alignItems: 'center', backgroundColor: colors.goldSoft, borderRadius: radius.sm, padding: 10 },
  });
