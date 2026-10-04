// Examens préparés par les élèves. Ce module ne dépend pas du contenu pédagogique :
// les règles du jeu (gamification.ts) peuvent l'utiliser sans charger toutes les matières.
import type { TrackId } from './types';

export interface TrackInfo {
  id: TrackId;
  label: string;
  description: string;
  emoji: string;
  /** Date d'examen indicative (AAAA-MM-JJ), modifiable par l'élève dans son profil. */
  defaultExamDate: string;
}

export const tracks: TrackInfo[] = [
  {
    id: 'bfm',
    label: 'BFEM',
    description: 'Brevet de Fin d’Études Moyennes — classe de 3e',
    emoji: '🎒',
    defaultExamDate: '2027-07-12',
  },
  {
    id: 'bac-s',
    label: 'Bac S',
    description: 'Terminale scientifique (S1, S2…)',
    emoji: '🔬',
    defaultExamDate: '2027-07-01',
  },
  {
    id: 'bac-l',
    label: 'Bac L',
    description: 'Terminale littéraire (L1, L2, L’…)',
    emoji: '📚',
    defaultExamDate: '2027-07-01',
  },
];

export function isTrackId(id: unknown): id is TrackId {
  return tracks.some((t) => t.id === id);
}

export function getTrack(id: TrackId): TrackInfo {
  return tracks.find((t) => t.id === id) ?? tracks[0];
}
