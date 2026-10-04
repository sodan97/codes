import { Modal, StyleSheet, Text, View } from 'react-native';

import { levelInfo, levelThreshold, type Reward } from '../lib/gamification';
import { colors, radius } from '../theme';
import { Button } from './ui';

/** Fenêtre de récompense affichée après une activité : XP, bonus, badges, niveau. */
export function RewardModal({ reward, onClose }: { reward: Reward | null; onClose: () => void }) {
  return (
    <Modal visible={!!reward} transparent animationType="fade" onRequestClose={onClose}>
      <View style={styles.backdrop}>
        {reward && (
          <View style={styles.box}>
            <Text style={styles.big}>+{reward.xp} XP</Text>
            {reward.levelUp != null && (
              <View style={styles.level}>
                <Text style={styles.levelText}>⬆️ Niveau {reward.levelUp} : {levelInfo(levelThreshold(reward.levelUp)).title}</Text>
              </View>
            )}
            {reward.messages.map((m) => (
              <Text key={m} style={styles.msg}>
                {m}
              </Text>
            ))}
            {reward.newBadges.length > 0 && (
              <View style={{ gap: 8, marginTop: 8, alignSelf: 'stretch' }}>
                <Text style={styles.badgeTitle}>Nouveau badge !</Text>
                {reward.newBadges.map((b) => (
                  <View key={b.id} style={styles.badge}>
                    <Text style={{ fontSize: 28 }}>{b.icon}</Text>
                    <View style={{ flex: 1 }}>
                      <Text style={{ fontWeight: '800', color: colors.text }}>{b.name}</Text>
                      <Text style={{ color: colors.muted, fontSize: 13 }}>{b.description}</Text>
                    </View>
                  </View>
                ))}
              </View>
            )}
            <Button label="Super !" onPress={onClose} style={{ alignSelf: 'stretch', marginTop: 12 }} />
          </View>
        )}
      </View>
    </Modal>
  );
}

const styles = StyleSheet.create({
  backdrop: { flex: 1, backgroundColor: 'rgba(0,0,0,0.45)', alignItems: 'center', justifyContent: 'center', padding: 24 },
  box: { backgroundColor: colors.card, borderRadius: radius.lg, padding: 22, alignItems: 'center', gap: 6, width: '100%', maxWidth: 420 },
  big: { fontSize: 40, fontWeight: '900', color: colors.primary },
  msg: { fontSize: 15, color: colors.text, textAlign: 'center' },
  level: { backgroundColor: colors.goldSoft, borderRadius: radius.sm, paddingVertical: 6, paddingHorizontal: 12, marginBottom: 4 },
  levelText: { fontWeight: '800', color: colors.text },
  badgeTitle: { fontWeight: '800', color: colors.gold, textAlign: 'center', fontSize: 16 },
  badge: { flexDirection: 'row', gap: 12, alignItems: 'center', backgroundColor: colors.goldSoft, borderRadius: radius.sm, padding: 10 },
});
