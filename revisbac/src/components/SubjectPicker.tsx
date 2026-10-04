import { Pressable, StyleSheet, Text, View } from 'react-native';

import { getSubjects } from '../data/catalog';
import type { TrackId } from '../data/types';
import { colors, radius } from '../theme';
import { styles as ui } from './ui';

/**
 * Choix des matières passées par l'élève : une ligne par matière de l'examen, cochée si elle est visible.
 * Au moins une matière reste toujours cochée.
 */
export function SubjectPicker({ track, hidden, onChange }: { track: TrackId; hidden: string[]; onChange: (hidden: string[]) => void }) {
  const subjects = getSubjects(track);
  const visibleCount = subjects.filter((s) => !hidden.includes(s.id)).length;

  const toggle = (id: string) => {
    if (hidden.includes(id)) onChange(hidden.filter((h) => h !== id));
    else if (visibleCount > 1) onChange([...hidden, id]);
  };

  return (
    <View style={{ gap: 8 }}>
      <Text style={[ui.body, { color: colors.muted }]}>
        Décoche les matières que tu ne passes pas (par exemple l’espagnol si tu fais arabe, allemand ou portugais).
      </Text>
      {subjects.map((s) => {
        const checked = !hidden.includes(s.id);
        const locked = checked && visibleCount === 1;
        return (
          <Pressable
            key={s.id}
            onPress={() => toggle(s.id)}
            disabled={locked}
            accessibilityRole="checkbox"
            accessibilityState={{ checked, disabled: locked }}
            style={[styles.row, checked && { borderColor: colors.primary }]}
          >
            <Text style={{ fontSize: 24 }}>{s.icon}</Text>
            <Text style={[styles.name, !checked && { color: colors.muted }]}>{s.name}</Text>
            <View style={[styles.box, checked && { backgroundColor: colors.primary, borderColor: colors.primary }]}>
              {checked && <Text style={styles.check}>✓</Text>}
            </View>
          </Pressable>
        );
      })}
    </View>
  );
}

const styles = StyleSheet.create({
  row: {
    flexDirection: 'row',
    alignItems: 'center',
    gap: 12,
    borderWidth: 2,
    borderColor: colors.border,
    borderRadius: radius.md,
    paddingHorizontal: 14,
    paddingVertical: 10,
    backgroundColor: colors.card,
  },
  name: { flex: 1, fontSize: 16, fontWeight: '700', color: colors.text },
  box: { width: 26, height: 26, borderRadius: 6, borderWidth: 2, borderColor: colors.border, alignItems: 'center', justifyContent: 'center' },
  check: { color: '#fff', fontWeight: '900', fontSize: 15 },
});
