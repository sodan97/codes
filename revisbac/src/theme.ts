import { Platform } from 'react-native';

export type Scheme = 'light' | 'dark';

// Palette inspirée du drapeau sénégalais (vert, jaune, rouge) sur fond clair.
export const lightColors = {
  primary: '#00853F',
  /** Vert plus contrasté : texte, liens. */
  primaryDark: '#006B32',
  primarySoft: '#E3F4EA',
  gold: '#F5B700',
  /** Or foncé pour le TEXTE sur fond clair (l'or vif est illisible en petit : contraste ≈ 1,8:1). */
  goldText: '#8A6100',
  goldSoft: '#FFF6D6',
  red: '#E31B23',
  redSoft: '#FDE7E8',
  bg: '#F6F7F4',
  card: '#FFFFFF',
  text: '#1B1F1D',
  muted: '#68706B',
  border: '#E3E6E1',
  /** Grand bandeau vert foncé (compte à rebours), texte blanc dans les deux thèmes. */
  hero: '#006B32',
  /** Voile derrière les fenêtres (récompenses). */
  backdrop: 'rgba(0,0,0,0.45)',
};

export type Colors = typeof lightColors;

// Même identité sur fond très sombre légèrement verdâtre : les couleurs du TEXTE sont éclaircies
// (contraste ≥ 4,5:1 sur les cartes), les fonds teintés sont assombris.
export const darkColors: Colors = {
  primary: '#3BB273',
  primaryDark: '#7FD6A4',
  primarySoft: '#16342A',
  gold: '#F5B700',
  goldText: '#F2C94C',
  goldSoft: '#3A3113',
  red: '#FF7A7A',
  redSoft: '#40201F',
  bg: '#0D1310',
  card: '#18211C',
  text: '#E7ECE8',
  muted: '#A3ADA6',
  border: '#2E3A33',
  hero: '#0E4F2D',
  backdrop: 'rgba(0,0,0,0.65)',
};

export const palettes: Record<Scheme, Colors> = { light: lightColors, dark: darkColors };

export const radius = { sm: 8, md: 14, lg: 20 };

export const shadow = {
  shadowColor: '#000',
  shadowOpacity: 0.06,
  shadowRadius: 8,
  shadowOffset: { width: 0, height: 2 },
  elevation: 2,
} as const;

/** Animations sur le moteur natif, absent dans le navigateur (version web). */
export const nativeDriver = Platform.OS !== 'web';

// --- Contrastes (WCAG) ---

const WHITE = '#FFFFFF';
/** Encre sombre posée sur les fonds vifs ou clairs. */
const INK = lightColors.text;

function rgb(hex: string): [number, number, number] {
  const h = hex.replace('#', '');
  const full = h.length === 3 ? h.replace(/./g, (c) => c + c) : h.slice(0, 6);
  const n = parseInt(full, 16);
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
}

function hex([r, g, b]: [number, number, number]): string {
  return '#' + [r, g, b].map((v) => Math.round(v).toString(16).padStart(2, '0')).join('').toUpperCase();
}

function luminance(color: string): number {
  const [r, g, b] = rgb(color).map((v) => {
    const c = v / 255;
    return c <= 0.03928 ? c / 12.92 : ((c + 0.055) / 1.055) ** 2.4;
  });
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}

/** Rapport de contraste entre deux couleurs opaques (#RRGGBB), de 1 à 21. */
export function contrast(a: string, b: string): number {
  const [hi, lo] = [luminance(a), luminance(b)].sort((x, y) => y - x);
  return (hi + 0.05) / (lo + 0.05);
}

/** Mélange `a` et `b` (t = 0 → a, t = 1 → b). */
function mix(a: string, b: string, t: number): string {
  const [ca, cb] = [rgb(a), rgb(b)];
  return hex([0, 1, 2].map((i) => ca[i] + (cb[i] - ca[i]) * t) as [number, number, number]);
}

/**
 * Couleur du texte posé sur un fond plein (bouton, bandeau de matière…) : blanc si lisible, sinon encre sombre.
 * À utiliser opaque : un suffixe alpha (textOn(x) + 'D9') ferait passer le contraste sous 4,5:1.
 */
export function textOn(fill: string): string {
  return contrast(WHITE, fill) >= 4.5 || contrast(WHITE, fill) >= contrast(INK, fill) ? WHITE : INK;
}

const accentCache = new Map<string, string>();

/**
 * Couleur d'accent (matière, vert, rouge…) utilisée pour du TEXTE : éclaircie en mode sombre,
 * assombrie en mode clair, juste assez pour un contraste ≥ 4,5:1 sur le fond et sur les cartes,
 * y compris teintés de cette couleur (jusqu'à ~12 % : color + '0F', '14', '1A', pastilles fg + '1A').
 */
export function accentText(color: string, palette: Colors): string {
  const key = `${color}|${palette.bg}`;
  const cached = accentCache.get(key);
  if (cached) return cached;
  const target = luminance(palette.bg) < 0.5 ? WHITE : '#000000';
  // Fond uni, ou teinté de la couleur d'origine ou de l'accent lui-même.
  const readable = (c: string) =>
    [palette.bg, palette.card].every((s) => [s, mix(color, s, 0.88), mix(c, s, 0.88)].every((f) => contrast(c, f) >= 4.5));
  let result = color;
  for (let t = 0.05; !readable(result) && t <= 1; t += 0.05) result = mix(color, target, t);
  accentCache.set(key, result);
  return result;
}
