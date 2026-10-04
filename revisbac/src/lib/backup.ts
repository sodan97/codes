// Sauvegarde de la progression dans un code à garder (message, note…) : sans compte ni serveur.
// Format : « RB2: » + base64(JSON de l'état, sans l'historique des XP) + « : » + somme de contrôle FNV-1a (8 caractères hexa).
// RB2 range les flashcards et les erreurs suivies de façon compacte (voir packCards) pour que le code tienne
// dans un message WhatsApp ; les codes RB1 (état JSON tel quel) restent lisibles.
// Module pur et testé : base64 et UTF-8 sont codés ici pour donner le même résultat sous Hermes, Node et dans le navigateur.
import { getTrack } from '../data/tracks';
import { addDays, dayKey, daysBetween, isDayKey } from './dates';
import { effectiveStreak, migrate, type MistakeEntry, type ProgressState } from './gamification';
import { fnv1a } from './hash';
import { CARD_BOX_DAYS, type CardEntry } from './srs';

export const BACKUP_PREFIX = 'RB2';
const BACKUP_VERSION = 2;

export const BACKUP_ERRORS = {
  incomplete: 'Code incomplet : copie tout le message, du code qui commence par « RB » jusqu’au dernier caractère.',
  damaged: 'Code abîmé (somme de contrôle) : un caractère a changé. Recopie le message d’origine.',
  version: 'Version inconnue : ce code vient d’une version plus récente de RéviBac. Mets l’application à jour puis réessaie.',
} as const;

/** Au-delà, le code risque d'être coupé dans un message (WhatsApp : 65 536 caractères) : mieux vaut une note. */
export const BACKUP_SAFE_LENGTH = 60_000;

/** Message partagé avec le code. */
export function backupShareMessage(code: string): string {
  return `Garde ce message : il permet de récupérer ta progression RéviBac.\n${code}`;
}

// ---------- UTF-8 et base64 ----------

const B64 = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';

function utf8Bytes(text: string): number[] {
  const bytes: number[] = [];
  for (const ch of text) {
    const c = ch.codePointAt(0)!;
    if (c < 0x80) bytes.push(c);
    else if (c < 0x800) bytes.push(0xc0 | (c >> 6), 0x80 | (c & 63));
    else if (c < 0x10000) bytes.push(0xe0 | (c >> 12), 0x80 | ((c >> 6) & 63), 0x80 | (c & 63));
    else bytes.push(0xf0 | (c >> 18), 0x80 | ((c >> 12) & 63), 0x80 | ((c >> 6) & 63), 0x80 | (c & 63));
  }
  return bytes;
}

function utf8Text(bytes: number[]): string {
  let out = '';
  for (let i = 0; i < bytes.length; ) {
    const b = bytes[i];
    const size = b < 0x80 ? 1 : b >> 5 === 6 ? 2 : b >> 4 === 14 ? 3 : b >> 3 === 30 ? 4 : 0;
    if (size === 0 || i + size > bytes.length) throw new Error('UTF-8 invalide');
    let c = size === 1 ? b : b & (0xff >> (size + 1));
    for (let k = 1; k < size; k++) {
      const next = bytes[i + k];
      if (next >> 6 !== 2) throw new Error('UTF-8 invalide');
      c = (c << 6) | (next & 63);
    }
    out += String.fromCodePoint(c);
    i += size;
  }
  return out;
}

/** Texte → base64 (via UTF-8). */
export function encodeBase64(text: string): string {
  const bytes = utf8Bytes(text);
  let out = '';
  for (let i = 0; i < bytes.length; i += 3) {
    const [a, b, c] = [bytes[i], bytes[i + 1], bytes[i + 2]];
    const n = (a << 16) | ((b ?? 0) << 8) | (c ?? 0);
    out += B64[(n >> 18) & 63] + B64[(n >> 12) & 63] + (b === undefined ? '=' : B64[(n >> 6) & 63]) + (c === undefined ? '=' : B64[n & 63]);
  }
  return out;
}

/** base64 → texte (UTF-8). Lève une erreur si le code n'est pas du base64 valide. */
export function decodeBase64(b64: string): string {
  if (b64.length % 4 !== 0 || !/^[A-Za-z0-9+/]*={0,2}$/.test(b64)) throw new Error('base64 invalide');
  const bytes: number[] = [];
  for (let i = 0; i < b64.length; i += 4) {
    const chunk = b64.slice(i, i + 4);
    const n = [...chunk].reduce((acc, ch) => (acc << 6) | (ch === '=' ? 0 : B64.indexOf(ch)), 0);
    bytes.push((n >> 16) & 255);
    if (chunk[2] !== '=') bytes.push((n >> 8) & 255);
    if (chunk[3] !== '=') bytes.push(n & 255);
  }
  return utf8Text(bytes);
}

// ---------- Format compact (RB2) ----------

/**
 * Flashcards groupées par chapitre, dates en jours relatifs au jour de référence `d` :
 * { [chapitre]: { [reste de la clé après '#']: [boîte, échéance] ou [boîte, échéance, dernière révision] } }.
 * La dernière révision est omise quand elle se déduit de l'échéance (échéance − durée de la boîte).
 * Les clés sans '#' (identifiant de carte) vont dans le groupe "".
 */
type Packed = Record<string, Record<string, number[]>>;

/** Coupe une clé en (groupe, reste) au premier '#' (cartes) ou au dernier '-' (questions : « chapitre-q3 »). */
function split(key: string, sep: '#' | '-'): [string, string] {
  const cut = sep === '#' ? key.indexOf(sep) : key.lastIndexOf(sep);
  return cut > 0 ? [key.slice(0, cut), key.slice(cut + 1)] : ['', key];
}
const join = (group: string, rest: string, sep: '#' | '-') => (group ? `${group}${sep}${rest}` : rest);

function packCards(cards: Record<string, CardEntry>, d: string): Packed {
  const out: Packed = {};
  for (const [key, e] of Object.entries(cards)) {
    const [group, rest] = split(key, '#');
    const due = daysBetween(d, e.due);
    const implied = addDays(e.due, -CARD_BOX_DAYS[e.box]) === e.last;
    (out[group] ??= {})[rest] = implied ? [e.box, due] : [e.box, due, daysBetween(d, e.last)];
  }
  return out;
}

/** Erreurs suivies groupées de la même façon : [boîte, échéance, erreurs] ou [boîte, échéance, erreurs, dernier changement]. */
function packMistakes(mistakes: Record<string, MistakeEntry>, d: string): Packed {
  const out: Packed = {};
  for (const [id, e] of Object.entries(mistakes)) {
    const [group, rest] = split(id, '-');
    const row = [e.box, daysBetween(d, e.due), e.fails];
    (out[group] ??= {})[rest] = e.last ? [...row, daysBetween(d, e.last)] : row;
  }
  return out;
}

const isObj = (v: unknown): v is Record<string, unknown> => typeof v === 'object' && v !== null && !Array.isArray(v);
/** Jour relatif → date (undefined si le décalage est illisible : migrate() complète alors la valeur). */
const dayAt = (d: string, n: unknown) => (typeof n === 'number' && Number.isInteger(n) ? addDays(d, n) : undefined);

/** Remet cartes et erreurs à leur forme habituelle ; migrate() vérifie ensuite chaque entrée. */
function unpack(data: Record<string, unknown>): Record<string, unknown> {
  const { d, c, m, ...rest } = data;
  if (!isDayKey(d)) return rest;
  const cards: Record<string, unknown> = {};
  if (isObj(c)) {
    for (const [group, entries] of Object.entries(c)) {
      if (!isObj(entries)) continue;
      for (const [k, row] of Object.entries(entries)) {
        if (!Array.isArray(row)) continue;
        const [box, due, last] = row;
        const dueDay = dayAt(d, due);
        const implied = typeof box === 'number' && dueDay && box in CARD_BOX_DAYS ? addDays(dueDay, -CARD_BOX_DAYS[box as CardEntry['box']]) : undefined;
        cards[join(group, k, '#')] = { box, due: dueDay, last: row.length > 2 ? dayAt(d, last) : implied };
      }
    }
  }
  const mistakes: Record<string, unknown> = {};
  if (isObj(m)) {
    for (const [group, entries] of Object.entries(m)) {
      if (!isObj(entries)) continue;
      for (const [k, row] of Object.entries(entries)) {
        if (!Array.isArray(row)) continue;
        const [box, due, fails, last] = row;
        const lastDay = row.length > 3 ? dayAt(d, last) : undefined;
        mistakes[join(group, k, '-')] = { box, due: dayAt(d, due), fails, ...(lastDay ? { last: lastDay } : {}) };
      }
    }
  }
  return { ...rest, cards, mistakes };
}

// ---------- Export / import ----------

/** Code de sauvegarde de la progression (l'historique des XP des 30 derniers jours n'y est pas). */
export function exportCode(state: ProgressState): string {
  const { history: _history, cards, mistakes, ...rest } = state;
  // Jour de référence des dates relatives : le jour de l'état (toujours une date valide).
  const d = state.today.day;
  const payload = encodeBase64(JSON.stringify({ ...rest, d, c: packCards(cards, d), m: packMistakes(mistakes, d) }));
  return `${BACKUP_PREFIX}:${payload}:${fnv1a(payload)}`;
}

export interface BackupPreview {
  name: string;
  trackLabel: string;
  xp: number;
  streak: number;
  badges: number;
}

export type ImportResult = { ok: true; state: ProgressState; preview: BackupPreview } | { ok: false; error: string };

/**
 * Lit un code de sauvegarde, même entouré du message partagé ou coupé par des retours à la ligne.
 * La sauvegarde passe par migrate() : un code d'une version plus ancienne est mis à jour.
 */
export function importCode(text: string, today = dayKey()): ImportResult {
  const compact = text.replace(/\s+/g, '');
  const start = compact.search(/RB\d+:/);
  if (start < 0) return { ok: false, error: BACKUP_ERRORS.incomplete };
  const match = /^RB(\d+):([A-Za-z0-9+/=]*)(?::([0-9a-fA-F]{0,8}))?/.exec(compact.slice(start))!;
  // RB1 : état JSON tel quel ; RB2 : flashcards et erreurs compactes.
  const format = Number(match[1]);
  if (format < 1 || format > BACKUP_VERSION) return { ok: false, error: BACKUP_ERRORS.version };
  const [, , payload, checksum] = match;
  if (!payload || payload.length % 4 !== 0 || !checksum || checksum.length < 8) return { ok: false, error: BACKUP_ERRORS.incomplete };
  if (fnv1a(payload) !== checksum.toLowerCase()) return { ok: false, error: BACKUP_ERRORS.damaged };
  let state: ProgressState;
  try {
    const data: unknown = JSON.parse(decodeBase64(payload));
    state = migrate(format >= 2 && isObj(data) ? unpack(data) : data, today);
  } catch {
    return { ok: false, error: BACKUP_ERRORS.damaged };
  }
  // Examen inconnu de cette version de l'application.
  if (!state.profile) return { ok: false, error: BACKUP_ERRORS.version };
  return { ok: true, state, preview: backupPreview(state, today) };
}

export function backupPreview(state: ProgressState, today = dayKey()): BackupPreview {
  return {
    name: state.profile?.name ?? 'Champion',
    trackLabel: state.profile ? getTrack(state.profile.track).label : '',
    xp: state.xp,
    streak: effectiveStreak(state, today),
    badges: Object.keys(state.badges).length,
  };
}

/** Nombre à la française : 1 240 (espace insécable). */
function formatInt(n: number): string {
  return String(Math.round(n)).replace(/\B(?=(\d{3})+(?!\d))/g, ' ');
}

/** « Awa · Bac S · 1 240 XP · série 12 · 9 badges » */
export function previewText(p: BackupPreview): string {
  return [p.name, p.trackLabel, `${formatInt(p.xp)} XP`, `série ${p.streak}`, `${p.badges} badge${p.badges > 1 ? 's' : ''}`].filter(Boolean).join(' · ');
}
