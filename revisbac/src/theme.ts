import { Platform } from 'react-native';

// Palette inspirée du drapeau sénégalais (vert, jaune, rouge) sur fond clair.
export const colors = {
  primary: '#00853F',
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
};

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
