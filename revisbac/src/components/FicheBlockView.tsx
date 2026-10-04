import { StyleSheet, Text, View } from 'react-native';

import type { FicheBlock } from '../data/types';
import { colors, radius } from '../theme';

/** Affiche un bloc de fiche avec un style propre à sa nature (formule, date, piège…). */
export function FicheBlockView({ block, color }: { block: FicheBlock; color: string }) {
  switch (block.kind) {
    case 'text':
      return <Text style={styles.text}>{block.text}</Text>;
    case 'list':
      return (
        <View style={{ gap: 4 }}>
          {block.title && <Text style={styles.listTitle}>{block.title}</Text>}
          {block.items.map((item, i) => (
            <View key={i} style={styles.bulletRow}>
              <Text style={[styles.bullet, { color }]}>•</Text>
              <Text style={[styles.text, { flex: 1 }]}>{item}</Text>
            </View>
          ))}
        </View>
      );
    case 'formula':
      return (
        <View style={[styles.formula, { borderColor: color, backgroundColor: color + '0F' }]}>
          {block.label && <Text style={[styles.formulaLabel, { color }]}>{block.label}</Text>}
          <Text style={styles.formulaText} selectable>
            {block.formula}
          </Text>
          {block.note && <Text style={styles.note}>{block.note}</Text>}
        </View>
      );
    case 'definition':
      return (
        <View style={[styles.definition, { borderLeftColor: color }]}>
          <Text style={styles.term}>{block.term}</Text>
          <Text style={styles.text}>{block.definition}</Text>
        </View>
      );
    case 'date':
      return (
        <View style={styles.dateRow}>
          <View style={[styles.dateBadge, { backgroundColor: color }]}>
            <Text style={styles.dateText}>{block.date}</Text>
          </View>
          <Text style={[styles.text, { flex: 1 }]}>{block.event}</Text>
        </View>
      );
    case 'tip':
      return (
        <View style={[styles.callout, { backgroundColor: colors.primarySoft }]}>
          <Text style={styles.calloutTitle}>💡 Méthode</Text>
          <Text style={styles.text}>{block.text}</Text>
        </View>
      );
    case 'warning':
      return (
        <View style={[styles.callout, { backgroundColor: colors.redSoft }]}>
          <Text style={[styles.calloutTitle, { color: colors.red }]}>⚠️ Piège à éviter</Text>
          <Text style={styles.text}>{block.text}</Text>
        </View>
      );
    case 'example':
      return (
        <View style={[styles.callout, { backgroundColor: colors.goldSoft }]}>
          <Text style={styles.calloutTitle}>✏️ {block.title ?? 'Exemple'}</Text>
          <Text style={styles.text}>{block.text}</Text>
        </View>
      );
    default:
      // Bloc d'un type inconnu (contenu plus récent ou mal formé) : ignoré plutôt que de faire planter la fiche.
      if (__DEV__) console.warn('Bloc de fiche inconnu', (block as { kind?: unknown }).kind);
      return null;
  }
}

const styles = StyleSheet.create({
  text: { fontSize: 15, lineHeight: 22, color: colors.text },
  listTitle: { fontSize: 15, fontWeight: '700', color: colors.text },
  bulletRow: { flexDirection: 'row', gap: 8 },
  bullet: { fontSize: 18, lineHeight: 22, fontWeight: '900' },
  formula: { borderWidth: 1.5, borderRadius: radius.sm, padding: 12, gap: 4 },
  formulaLabel: { fontSize: 12, fontWeight: '800', textTransform: 'uppercase', letterSpacing: 0.5 },
  formulaText: { fontSize: 18, fontWeight: '600', color: colors.text, lineHeight: 26 },
  note: { fontSize: 13, color: colors.muted, lineHeight: 18 },
  definition: { borderLeftWidth: 4, paddingLeft: 12, gap: 2 },
  term: { fontSize: 15, fontWeight: '800', color: colors.text },
  dateRow: { flexDirection: 'row', gap: 10, alignItems: 'flex-start' },
  dateBadge: { borderRadius: 6, paddingHorizontal: 8, paddingVertical: 3, minWidth: 64, alignItems: 'center' },
  dateText: { color: '#fff', fontWeight: '800', fontSize: 13 },
  callout: { borderRadius: radius.sm, padding: 12, gap: 4 },
  calloutTitle: { fontSize: 13, fontWeight: '800', color: colors.primaryDark },
});
