// Retours tactiles sobres. Le réglage (profile.haptics) est passé en paramètre :
// le module ne dépend d'aucun contexte React. Sans effet sur le web.
import * as Haptics from 'expo-haptics';
import { Platform } from 'react-native';

function run(enabled: boolean, vibrate: () => Promise<void>) {
  if (!enabled || Platform.OS === 'web') return;
  try {
    vibrate().catch(() => {});
  } catch {
    // Pas de moteur de vibration : on ignore.
  }
}

/** Bonne réponse. */
export function good(enabled: boolean) {
  run(enabled, () => Haptics.notificationAsync(Haptics.NotificationFeedbackType.Success));
}

/** Mauvaise réponse : volontairement doux. */
export function bad(enabled: boolean) {
  run(enabled, () => Haptics.impactAsync(Haptics.ImpactFeedbackStyle.Light));
}

/** Montée de niveau, badge… */
export function celebrate(enabled: boolean) {
  run(enabled, () => Haptics.notificationAsync(Haptics.NotificationFeedbackType.Success));
}
