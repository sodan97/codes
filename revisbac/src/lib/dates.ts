/** Date locale au format AAAA-MM-JJ. */
export function dayKey(date: Date = new Date()): string {
  const y = date.getFullYear();
  const m = String(date.getMonth() + 1).padStart(2, '0');
  const d = String(date.getDate()).padStart(2, '0');
  return `${y}-${m}-${d}`;
}

export function parseDay(key: string): Date {
  const [y, m, d] = key.split('-').map(Number);
  return new Date(y, m - 1, d);
}

/** Nombre de jours entre deux clés de date (b - a). */
export function daysBetween(a: string, b: string): number {
  return Math.round((parseDay(b).getTime() - parseDay(a).getTime()) / 86_400_000);
}

export function addDays(key: string, n: number): string {
  const d = parseDay(key);
  d.setDate(d.getDate() + n);
  return dayKey(d);
}

/** Ajoute n mois (négatif possible). Le jour est ramené au dernier jour du mois si besoin (31 janv. + 1 mois → 28 ou 29 févr.). */
export function addMonths(key: string, n: number): string {
  const [y, m, d] = key.split('-').map(Number);
  const target = new Date(y, m - 1 + n, 1);
  const lastDay = new Date(target.getFullYear(), target.getMonth() + 1, 0).getDate();
  target.setDate(Math.min(d, lastDay));
  return dayKey(target);
}

/** Vrai si la chaîne est une date au format AAAA-MM-JJ. */
export function isDayKey(value: unknown): value is string {
  return typeof value === 'string' && /^\d{4}-\d{2}-\d{2}$/.test(value);
}

const MONTHS = ['janv.', 'févr.', 'mars', 'avr.', 'mai', 'juin', 'juil.', 'août', 'sept.', 'oct.', 'nov.', 'déc.'];

export function formatDay(key: string): string {
  const d = parseDay(key);
  return `${d.getDate()} ${MONTHS[d.getMonth()]} ${d.getFullYear()}`;
}

/** Abréviation du jour de la semaine (Lu, Ma, Me, Je, Ve, Sa, Di) : mardi et mercredi restent distincts. */
export function weekdayLetter(key: string): string {
  return ['Di', 'Lu', 'Ma', 'Me', 'Je', 'Ve', 'Sa'][parseDay(key).getDay()];
}
