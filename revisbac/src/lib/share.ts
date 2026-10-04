// Partage d'un texte (WhatsApp, SMS…) par la feuille de partage du téléphone.
import { Platform, Share } from 'react-native';

export { dailyShareText, examShareText } from './shareText';

/** shared : envoyé (ou feuille de partage ouverte) ; copied : copié dans le presse-papiers (web) ; failed : rien. */
export type ShareOutcome = 'shared' | 'copied' | 'failed';

/** Ouvre la feuille de partage ; sur le web sans partage disponible, copie le texte. */
export async function share(message: string): Promise<ShareOutcome> {
  if (Platform.OS === 'web') return shareOnWeb(message);
  try {
    const result = await Share.share({ message });
    return result.action === Share.sharedAction ? 'shared' : 'failed';
  } catch (e) {
    console.warn('Partage impossible', e);
    return 'failed';
  }
}

async function shareOnWeb(message: string): Promise<ShareOutcome> {
  const nav = typeof navigator === 'undefined' ? undefined : navigator;
  if (nav?.share) {
    try {
      await nav.share({ text: message });
      return 'shared';
    } catch {
      // Partage refusé ou fermé : on copie le texte à la place.
    }
  }
  try {
    if (!nav?.clipboard) return 'failed';
    await nav.clipboard.writeText(message);
    return 'copied';
  } catch {
    return 'failed';
  }
}
