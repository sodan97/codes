import { useEffect, useState } from 'react';
import { Platform, Text, View } from 'react-native';

import { share } from '../lib/share';
import { Button, useUi } from './ui';

/**
 * Bouton « Partager » : feuille de partage sur le téléphone. Sur le web, si le partage n'est pas
 * disponible, le texte est copié et un petit message le dit.
 */
export function ShareButton({ label, message, color }: { label: string; message: () => string | null; color?: string }) {
  const ui = useUi();
  const [note, setNote] = useState<string | null>(null);

  useEffect(() => {
    if (!note) return;
    const id = setTimeout(() => setNote(null), 5000);
    return () => clearTimeout(id);
  }, [note]);

  const onPress = async () => {
    const text = message();
    if (!text) return;
    const outcome = await share(text);
    if (outcome === 'copied') setNote('✓ Copié ! Colle-le dans WhatsApp ou un SMS.');
    else if (outcome === 'failed' && Platform.OS === 'web') setNote('Partage impossible ici : fais une capture d’écran de ton score.');
  };

  return (
    <View style={{ gap: 6 }}>
      <Button label={label} variant="secondary" color={color} onPress={() => void onPress()} />
      {note && (
        <Text style={[ui.muted, { textAlign: 'center' }]} accessibilityLiveRegion="polite">
          {note}
        </Text>
      )}
    </View>
  );
}
