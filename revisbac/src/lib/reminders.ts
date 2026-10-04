// Rappel quotidien : notifications LOCALES planifiées sur le téléphone (aucun envoi depuis un serveur).
// La partie texte et calendrier est pure et testée. expo-notifications n'est chargé qu'à l'appel
// (import dynamique) pour que ce module reste importable dans les tests, hors React Native.
import { getTrack } from '../data/catalog';
import { addDays, dayKey, daysBetween, parseDay } from './dates';
import { effectiveStreak, todayStats, XP, type ProgressState } from './gamification';
import { DAILY_SIZE } from './quizBuilder';
import { activeMistakes } from './selectors';

/** Heures proposées (onboarding et Profil). */
export const REMINDER_HOURS = [7, 13, 17, 19, 21];
export const DEFAULT_REMINDER_HOUR = 19;
/** Les rappels sont planifiés pour 7 jours : sans retour dans l'appli, ils s'arrêtent d'eux-mêmes. */
export const REMINDER_DAYS = 7;
const ID_PREFIX = 'rappel-';
const CHANNEL_ID = 'rappels';

/** « 19 h 00 » */
export function formatHour(hour: number, minute = 0): string {
  return `${hour} h ${String(minute).padStart(2, '0')}`;
}

/**
 * Texte du rappel prévu pour le jour `day` (premier message applicable), null sans profil.
 * Jamais de menace ni de culpabilisation.
 */
export function reminderText(state: ProgressState, day: string, today = dayKey()): string | null {
  const profile = state.profile;
  if (!profile) return null;
  const track = getTrack(profile.track).label;
  const daysLeft = daysBetween(day, profile.examDate);

  if (daysLeft === 0) return `Bonne chance pour ton ${track} aujourd’hui, ${profile.name} ! 🍀`;
  if (daysLeft === 1) return `Demain c’est le ${track} : relis tes essentiels et dors tôt.`;

  const stats = todayStats(state, today);
  if (day === today && stats.xp > 0 && stats.xp < profile.dailyGoal) {
    return `Plus que ${profile.dailyGoal - stats.xp} XP pour ton objectif du jour 💪`;
  }

  const streak = effectiveStreak(state, day);
  if (streak >= 3) return `🔥 Ta série de ${streak} jours t’attend, ${profile.name} ! Une seule activité suffit.`;

  const due = activeMistakes(state, day).due.length;
  if (due >= 1) {
    return `🔁 ${due} question${due > 1 ? 's' : ''} à revoir aujourd’hui (≈ ${Math.max(1, Math.round(due / 2))} min)`;
  }

  // Rotation selon le jour.
  const options: string[] = [];
  if (!(day === today && stats.challengeDone)) options.push(`🎯 Le défi du jour est prêt : ${DAILY_SIZE} questions, +${XP.dailyChallenge} XP`);
  options.push('⚡ 5 minutes de révision express ?');
  if (daysLeft > 1 && daysLeft <= 60) options.push(`J-${daysLeft} avant le ${track} : une fiche de 3 minutes ?`);
  const n = daysBetween('2026-01-01', day);
  return options[((n % options.length) + options.length) % options.length];
}

export interface PlannedReminder {
  /** `rappel-AAAA-MM-JJ` */
  id: string;
  day: string;
  date: Date;
  body: string;
}

/**
 * Rappels à planifier (7 au plus) : aujourd'hui si l'heure n'est pas passée et l'objectif pas atteint,
 * puis les 6 jours suivants. Vide si le rappel est désactivé.
 */
export function planReminders(state: ProgressState, now: Date = new Date()): PlannedReminder[] {
  const reminder = state.profile?.reminder;
  if (!state.profile || !reminder?.enabled) return [];
  const today = dayKey(now);
  const out: PlannedReminder[] = [];
  for (let i = 0; i < REMINDER_DAYS; i++) {
    const day = addDays(today, i);
    const date = parseDay(day);
    date.setHours(reminder.hour, reminder.minute, 0, 0);
    if (i === 0 && (date.getTime() <= now.getTime() || todayStats(state, today).xp >= state.profile.dailyGoal)) continue;
    const body = reminderText(state, day, today);
    if (body) out.push({ id: `${ID_PREFIX}${day}`, day, date, body });
  }
  return out;
}

const isWeb = () => process.env.EXPO_OS === 'web';

async function notifications() {
  return import('expo-notifications');
}

async function ensureChannel(N: Awaited<ReturnType<typeof notifications>>) {
  if (process.env.EXPO_OS !== 'android') return;
  await N.setNotificationChannelAsync(CHANNEL_ID, { name: 'Rappels de révision', importance: N.AndroidImportance.DEFAULT });
}

export type ReminderPermission = 'granted' | 'denied' | 'undetermined';

/**
 * Traduit la réponse de getPermissionsAsync. Android 13 et plus : « denied » tant que la permission
 * n'a pas été demandée ; si on peut encore la demander, rien n'est bloqué.
 */
export function permissionState(status: string, canAskAgain: boolean): ReminderPermission {
  if (status === 'granted') return 'granted';
  return status === 'denied' && !canAskAgain ? 'denied' : 'undetermined';
}

/** État de la permission des notifications ('denied' aussi sur le web, où les rappels n'existent pas). */
export async function reminderPermission(): Promise<ReminderPermission> {
  if (isWeb()) return 'denied';
  try {
    const N = await notifications();
    const { status, canAskAgain } = await N.getPermissionsAsync();
    return permissionState(status, canAskAgain);
  } catch (e) {
    console.warn('Permission des notifications illisible', e);
    return 'denied';
  }
}

/** Demande la permission (à n'appeler qu'après le choix d'une heure). Vrai si elle est accordée. */
export async function requestReminderPermission(): Promise<boolean> {
  if (isWeb()) return false;
  try {
    const N = await notifications();
    // Android 13 et plus : la demande ne s'affiche qu'une fois un canal créé.
    await ensureChannel(N);
    const current = await N.getPermissionsAsync();
    if (current.granted) return true;
    return (await N.requestPermissionsAsync()).granted;
  } catch (e) {
    console.warn('Demande de permission impossible', e);
    return false;
  }
}

async function reschedule(state: ProgressState): Promise<void> {
  if (isWeb()) return;
  try {
    const N = await notifications();
    // On repart de zéro : les anciens rappels sont annulés, même si le rappel vient d'être désactivé.
    const scheduled = await N.getAllScheduledNotificationsAsync();
    await Promise.all(scheduled.filter((n) => n.identifier.startsWith(ID_PREFIX)).map((n) => N.cancelScheduledNotificationAsync(n.identifier)));
    if (!state.profile?.reminder?.enabled) return;
    if (!(await N.getPermissionsAsync()).granted) return;
    await ensureChannel(N);
    for (const r of planReminders(state)) {
      await N.scheduleNotificationAsync({
        identifier: r.id,
        content: { title: 'RéviBac', body: r.body },
        trigger: { type: N.SchedulableTriggerInputTypes.DATE, date: r.date, channelId: CHANNEL_ID },
      });
    }
  } catch (e) {
    console.warn('Rappels non planifiés', e);
  }
}

let queue: Promise<void> = Promise.resolve();

/**
 * Annule les rappels « rappel-… » et replanifie ceux des 7 prochains jours.
 * Sans effet sur le web ; rien n'est planifié si le rappel est désactivé ou la permission refusée.
 * Les appels passent l'un après l'autre.
 */
export function rescheduleReminders(state: ProgressState): Promise<void> {
  queue = queue.then(() => reschedule(state));
  return queue;
}
