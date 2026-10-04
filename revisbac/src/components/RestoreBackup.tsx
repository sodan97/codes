import { useMemo, useState } from 'react';
import { StyleSheet, Text, TextInput, View } from 'react-native';

import { importCode, previewText } from '../lib/backup';
import { dayKey } from '../lib/dates';
import { useProgress } from '../state/progress';
import { useStyles, useTheme } from '../state/theme';
import { radius, type Colors } from '../theme';
import { Button, useUi } from './ui';

/**
 * Restauration d'une sauvegarde : l'élève colle le message reçu (ou le code seul), voit à qui il appartient,
 * puis confirme. Rien n'est remplacé avant la confirmation.
 */
export function RestoreBackup({ onRestored, onCancel }: { onRestored: (name: string) => void; onCancel: () => void }) {
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
  const { state, restore } = useProgress();
  const [text, setText] = useState('');
  // Le code est vérifié dès qu'il est collé.
  const result = useMemo(() => (text.trim() ? importCode(text, dayKey()) : null), [text]);
  const hasProgress = !!state.profile;

  return (
    <View style={{ gap: 10 }}>
      <Text style={styles.label}>Colle ici le message de sauvegarde</Text>
      <TextInput
        value={text}
        onChangeText={setText}
        multiline
        placeholder="Le message contient un code qui commence par « RB »"
        placeholderTextColor={colors.muted}
        autoCapitalize="none"
        autoCorrect={false}
        spellCheck={false}
        accessibilityLabel="Message de sauvegarde"
        style={styles.input}
      />
      <Text style={ui.muted}>Appuie longuement dans le champ puis choisis « Coller ». Le message entier convient, pas besoin de l’arranger.</Text>

      {result && !result.ok && (
        <Text style={[ui.body, { color: colors.red }]} accessibilityLiveRegion="polite">
          {result.error}
        </Text>
      )}

      {result?.ok && (
        <View style={styles.preview} accessibilityLiveRegion="polite">
          <Text style={styles.previewTitle}>✓ Sauvegarde trouvée</Text>
          <Text style={ui.body}>{previewText(result.preview)}</Text>
          {hasProgress && <Text style={ui.muted}>Ta progression actuelle sur ce téléphone sera remplacée par celle-ci.</Text>}
          <Button
            label={hasProgress ? 'Remplacer ma progression actuelle' : 'Retrouver ma progression'}
            color={hasProgress ? colors.red : undefined}
            onPress={() => {
              restore(result.state);
              onRestored(result.preview.name);
            }}
          />
        </View>
      )}

      <Button label="Annuler" variant="ghost" color={colors.muted} onPress={onCancel} />
    </View>
  );
}

const makeStyles = (colors: Colors) =>
  StyleSheet.create({
    label: { fontSize: 15, fontWeight: '700', color: colors.text },
    input: {
      minHeight: 110,
      borderWidth: 2,
      borderColor: colors.border,
      borderRadius: radius.md,
      paddingHorizontal: 12,
      paddingVertical: 10,
      fontSize: 14,
      textAlignVertical: 'top',
      backgroundColor: colors.card,
      color: colors.text,
    },
    preview: { gap: 8, borderWidth: 2, borderColor: colors.primary, borderRadius: radius.md, padding: 12, backgroundColor: colors.primarySoft },
    previewTitle: { fontSize: 16, fontWeight: '800', color: colors.primaryDark },
  });
