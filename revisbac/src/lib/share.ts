// Partage d'un texte (WhatsApp, SMS…) par la feuille de partage du téléphone.
import { Share } from 'react-native';

export { dailyShareText, examShareText } from './shareText';

/** Ouvre la feuille de partage. Renvoie true si le texte a été partagé, false si annulé ou en cas d'échec. */
export async function share(message: string): Promise<boolean> {
  try {
    const result = await Share.share({ message });
    return result.action === Share.sharedAction;
  } catch (e) {
    console.warn('Partage impossible', e);
    return false;
  }
}
