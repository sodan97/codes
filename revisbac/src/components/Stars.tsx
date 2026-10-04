import { StyleSheet, Text, View } from 'react-native';

import { starsAccessibilityLabel } from '../lib/stats';
import { useStyles } from '../state/theme';
import type { Colors } from '../theme';
import { MAX_FONT_SCALE } from './ui';

/**
 * Étoiles de maîtrise d'un chapitre (★ gagnées, ☆ à gagner).
 * `faded` : chapitre à 3 étoiles pas pratiqué depuis longtemps, étoiles estompées mais jamais retirées.
 */
export function Stars({ stars, max = 3, faded = false, size = 18 }: { stars: number; max?: number; faded?: boolean; size?: number }) {
  const styles = useStyles(makeStyles);
  const label = starsAccessibilityLabel(stars, max) + (faded ? ', petit rappel conseillé' : '');
  return (
    <View style={styles.row} accessible accessibilityRole="image" accessibilityLabel={label}>
      {Array.from({ length: max }, (_, i) => {
        const on = i < stars;
        return (
          <Text
            key={i}
            style={[styles.star, { fontSize: size, lineHeight: Math.round(size * 1.25) }, on ? styles.on : styles.off, on && faded && styles.faded]}
            maxFontSizeMultiplier={MAX_FONT_SCALE}
          >
            {on ? '★' : '☆'}
          </Text>
        );
      })}
    </View>
  );
}

const makeStyles = (colors: Colors) =>
  StyleSheet.create({
    row: { flexDirection: 'row', alignItems: 'center', gap: 2 },
    star: { fontWeight: '700' },
    on: { color: colors.goldText },
    off: { color: colors.muted },
    faded: { opacity: 0.45 },
  });
