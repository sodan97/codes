/// <reference types="node" />
import assert from 'node:assert/strict';
import { test } from 'node:test';

import type { Subject } from '../src/data/types';
import { initialState } from '../src/lib/gamification';
import { nextStep, pace, phaseOf } from '../src/lib/plan';

const TODAY = '2026-10-04';

// Matières fictives : le plan ne dépend que de l'ordre et des identifiants des chapitres.
function subject(id: string, chapters: number): Subject {
  return {
    id,
    name: id.toUpperCase(),
    icon: '📘',
    color: '#000000',
    tracks: ['bac-s'],
    chapters: Array.from({ length: chapters }, (_, i) => ({
      id: `${id}-${i + 1}`,
      title: `${id} chapitre ${i + 1}`,
      summary: '',
      essentials: [],
      sections: [],
      flashcards: [],
      quiz: [],
    })),
  };
}

const maths = subject('maths', 3);
const svt = subject('svt', 2);
const subjects = [maths, svt];

function state(fichesRead: Record<string, string> = {}, quizBest: Record<string, number> = {}) {
  return { ...initialState(TODAY), fichesRead, quizBest };
}

test('nextStep a : la fiche lue le plus récemment sans quiz réussi', () => {
  const s = state({ 'maths-1': '2026-10-01', 'svt-1': '2026-10-03', 'maths-2': '2026-10-02' }, { 'svt-1': 80, 'maths-2': 40 });
  const step = nextStep(s, subjects, TODAY);
  assert.deepEqual(step, {
    kind: 'quiz',
    chapterId: 'maths-2',
    subjectId: 'maths',
    title: 'maths chapitre 2',
    reason: 'Fiche lue : teste-toi maintenant',
  });
});

test('nextStep b : premier chapitre non lu de la matière la moins avancée, dans l’ordre du programme', () => {
  // maths : 1/3 lu ; svt : 1/2 lu → maths est la moins avancée. Le chapitre 3 lu ne fait pas sauter le 2.
  const s = state({ 'maths-3': '2026-10-01', 'svt-1': '2026-10-01' }, { 'maths-3': 90, 'svt-1': 60 });
  const step = nextStep(s, subjects, TODAY);
  assert.equal(step?.kind, 'fiche');
  assert.equal(step?.chapterId, 'maths-1');
  assert.equal(step?.reason, 'Chapitre 1 de MATHS');

  // À égalité (rien de lu), l'ordre du catalogue départage.
  assert.equal(nextStep(state(), subjects, TODAY)?.chapterId, 'maths-1');
  assert.equal(nextStep(state(), [svt, maths], TODAY)?.chapterId, 'svt-1');
});

test('nextStep c puis d : tout est lu', () => {
  const read = Object.fromEntries([...maths.chapters, ...svt.chapters].map((c) => [c.id, '2026-09-01']));
  const best = { 'maths-1': 100, 'maths-2': 70, 'maths-3': 90, 'svt-1': 60, 'svt-2': 100 };
  const step = nextStep(state(read, best), subjects, TODAY);
  assert.equal(step?.kind, 'quiz');
  assert.equal(step?.chapterId, 'svt-1');
  assert.equal(step?.reason, 'Ton chapitre le plus fragile');

  const perfect = Object.fromEntries(Object.keys(read).map((id) => [id, 100]));
  assert.equal(nextStep(state(read, perfect), subjects, TODAY), null);
  assert.equal(nextStep(state(), [], TODAY), null);
});

test('nextStep ignore les matières absentes de la liste (masquées)', () => {
  const s = state({ 'svt-1': '2026-10-03' });
  assert.equal(nextStep(s, [maths], TODAY)?.chapterId, 'maths-1');
});

test('phases', () => {
  assert.equal(phaseOf(200), 'Découverte');
  assert.equal(phaseOf(91), 'Découverte');
  assert.equal(phaseOf(90), 'Consolidation');
  assert.equal(phaseOf(30), 'Consolidation');
  assert.equal(phaseOf(29), 'Sprint final');
  assert.equal(phaseOf(2), 'Sprint final');
  assert.equal(phaseOf(1), 'Veille');
  assert.equal(phaseOf(0), 'Jour J');
  assert.equal(phaseOf(-3), 'Passé');
});

test('pace : fiches restantes et rythme par semaine', () => {
  const s = state({ 'maths-1': '2026-10-01' });
  // 4 fiches restantes, 15 jours → 3 semaines → 2 par semaine.
  assert.deepEqual(pace(s, subjects, '2026-10-19', TODAY), {
    daysLeft: 15,
    remaining: 4,
    weeksLeft: 3,
    perWeek: 2,
    phase: 'Sprint final',
  });
  // Examen passé ou aujourd'hui : au moins une semaine.
  assert.equal(pace(s, subjects, '2026-09-01', TODAY).weeksLeft, 1);
  assert.equal(pace(s, subjects, '2026-09-01', TODAY).perWeek, 4);
  assert.equal(pace(s, subjects, TODAY, TODAY).phase, 'Jour J');
  // Tout est lu.
  const all = Object.fromEntries([...maths.chapters, ...svt.chapters].map((c) => [c.id, TODAY]));
  assert.equal(pace(state(all), subjects, '2027-07-01', TODAY).perWeek, 0);
});
