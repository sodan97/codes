import type { Subject } from '../types';

import mathsBfm from './maths-bfm';
import pcBfm from './pc-bfm';
import svtBfm from './svt-bfm';
import histoireBfm from './histoire-bfm';
import geoBfm from './geo-bfm';
import francaisBfm from './francais-bfm';
import anglaisBfm from './anglais-bfm';
import espagnolBfm from './espagnol-bfm';

import mathsS from './maths-s';
import mathsL from './maths-l';
import pcS from './pc-s';
import svtS from './svt-s';
import philoBac from './philo-bac';
import histoireBac from './histoire-bac';
import geoBac from './geo-bac';
import francaisBac from './francais-bac';
import anglaisBac from './anglais-bac';
import espagnolBac from './espagnol-bac';

/** Toutes les matières, dans l'ordre d'affichage. Pour ajouter une matière : créer le fichier puis l'ajouter ici. */
export const subjects: Subject[] = [
  mathsBfm, pcBfm, svtBfm, francaisBfm, histoireBfm, geoBfm, anglaisBfm, espagnolBfm,
  mathsS, pcS, svtS, mathsL, philoBac, francaisBac, histoireBac, geoBac, anglaisBac, espagnolBac,
];
