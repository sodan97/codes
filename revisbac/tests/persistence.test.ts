/// <reference types="node" />
import assert from 'node:assert/strict';
import { test } from 'node:test';

import { MAX_FREEZES, STATE_VERSION, initialState, migrate, rollDay } from '../src/lib/gamification';

const TODAY = '2026-10-04';

/** Sauvegarde telle qu'écrite par la version 1 de l'application. */
const V1 = {
  version: 1,
  profile: { name: 'Moussa', track: 'bfm', examDate: '2027-07-10', dailyGoal: 100 },
  xp: 420,
  streak: { current: 4, best: 9, lastDay: '2026-10-03', freezes: 1 },
  today: { day: '2026-10-03', xp: 60, goalBonusGiven: true, challengeDone: true },
  fichesRead: { 'maths-bfm-thales': '2026-09-20' },
  quizBest: { 'maths-bfm-thales': 80 },
  quizCount: 12,
  perfectCount: 2,
  challengesDone: 3,
  examBest: { 'maths-bfm': 13.5 },
  mistakes: { 'maths-bfm-thales-q1': 2, 'id-disparu': 1 },
  subjectsTouched: ['maths-bfm'],
  badges: { 'first-fiche': '2026-09-20' },
  history: { '2026-10-03': 60 },
};

test('migration v1 → v2 : tout est conservé, les erreurs deviennent des entrées Leitner', () => {
  const s = migrate(JSON.parse(JSON.stringify(V1)), TODAY);
  assert.equal(s.version, STATE_VERSION);
  assert.deepEqual(s.profile, {
    name: 'Moussa',
    track: 'bfm',
    examDate: '2027-07-10',
    dailyGoal: 100,
    hiddenSubjects: [],
    reminder: null,
    haptics: true,
  });
  assert.equal(s.xp, 420);
  assert.deepEqual(s.streak, V1.streak);
  assert.deepEqual(s.fichesRead, V1.fichesRead);
  assert.deepEqual(s.quizBest, V1.quizBest);
  assert.equal(s.quizCount, 12);
  assert.equal(s.examBest['maths-bfm'], 13.5);
  // Les ids inconnus du catalogue sont gardés (filtrés à l'affichage).
  assert.deepEqual(s.mistakes, {
    'maths-bfm-thales-q1': { box: 1, due: TODAY, fails: 2 },
    'id-disparu': { box: 1, due: TODAY, fails: 1 },
  });
  assert.deepEqual(s.badges, V1.badges);
  assert.deepEqual(s.dailyResults, {});
  assert.deepEqual(s.examHistory, {});
  assert.deepEqual(s.examLast, {});
  // Nouveau jour : compteurs du jour remis à zéro par rollDay.
  assert.deepEqual(s.today, { day: TODAY, xp: 0, goalBonusGiven: false, challengeDone: false, done: [] });
});

test('migration : une sauvegarde v2 relue donne le même état', () => {
  const once = migrate(V1, TODAY);
  once.mistakes.x = { box: 3, due: '2026-10-10', fails: 4 };
  once.dailyResults[TODAY] = { correct: 8, total: 10, grid: '✅✅❌✅✅✅✅❌✅✅' };
  once.today.done = ['daily', 'quiz:maths-bfm-thales'];
  once.profile!.reminder = { enabled: true, hour: 19, minute: 0 };
  const twice = migrate(JSON.parse(JSON.stringify(once)), TODAY);
  assert.deepEqual(twice, once);
});

test('migration : état partiel complété en profondeur', () => {
  const s = migrate(
    {
      version: 1,
      xp: 50,
      streak: { current: 2, best: 2, lastDay: '2026-10-03' }, // sans freezes
      profile: { name: 'Awa', track: 'bac-l', examDate: '2027-07-01' }, // sans dailyGoal
    },
    TODAY,
  );
  assert.equal(s.streak.freezes, 0);
  assert.equal(s.streak.current, 2);
  assert.equal(s.profile?.dailyGoal, 50);
  assert.deepEqual(s.today, initialState(TODAY).today);
  assert.deepEqual(s.fichesRead, {});
  assert.deepEqual(s.mistakes, {});
  assert.deepEqual(s.subjectsTouched, []);
});

test('migration : valeurs invalides bornées ou remplacées', () => {
  const s = migrate(
    {
      version: 2,
      xp: -12.7,
      quizCount: 'beaucoup',
      streak: { current: -3, best: 'x', lastDay: 'hier', freezes: 9 },
      profile: { name: '   ', track: 'bac-s', examDate: '1er juillet', dailyGoal: 42, reminder: { enabled: true, hour: 30, minute: -5 }, haptics: 'oui' },
      today: { day: TODAY, xp: 15.6, goalBonusGiven: 'non', challengeDone: true, done: ['daily', 3, 'daily'] },
      quizBest: { a: 140, b: 'x', c: -5 },
      examBest: { 'pc-s': 25 },
      examHistory: { 'pc-s': [{ day: '2026-10-01', note: 12 }, { day: 'x', note: 3 }, ...Array.from({ length: 6 }, (_, i) => ({ day: `2026-10-0${i + 2}`, note: i }))] },
      mistakes: { a: { box: 7, due: 'demain', fails: -1 }, b: 'cassé' },
      history: { [TODAY]: 30, bidon: 10 },
      dailyResults: { [TODAY]: { correct: 12, total: 10, grid: 7 } },
    },
    TODAY,
  );
  assert.equal(s.xp, 0);
  assert.equal(s.quizCount, 0);
  assert.deepEqual(s.streak, { current: 0, best: 0, lastDay: null, freezes: MAX_FREEZES });
  assert.equal(s.profile?.name, 'Champion');
  assert.equal(s.profile?.dailyGoal, 50);
  assert.equal(s.profile?.examDate, '2027-07-01'); // date indicative du Bac S
  assert.deepEqual(s.profile?.reminder, { enabled: true, hour: 23, minute: 0 });
  assert.equal(s.profile?.haptics, true);
  assert.deepEqual(s.today, { day: TODAY, xp: 15, goalBonusGiven: false, challengeDone: true, done: ['daily'] });
  assert.deepEqual(s.quizBest, { a: 100, c: 0 });
  assert.equal(s.examBest['pc-s'], 20);
  assert.equal(s.examHistory['pc-s'].length, 5);
  assert.deepEqual(s.mistakes, { a: { box: 1, due: TODAY, fails: 0 } });
  assert.deepEqual(s.history, { [TODAY]: 30 });
  assert.deepEqual(s.dailyResults, { [TODAY]: { correct: 10, total: 10, grid: '' } });
});

test('migration : examen inconnu → pas de profil (retour à l’onboarding)', () => {
  assert.equal(migrate({ ...V1, profile: { ...V1.profile, track: 'cfee' } }, TODAY).profile, null);
});

test('migration : une sauvegarde qui n’est pas un objet lève une erreur', () => {
  for (const raw of [null, 42, 'texte', [1, 2], undefined]) {
    assert.throws(() => migrate(raw, TODAY));
  }
});

test('rollDay vide today.done au changement de jour, sans rien toucher le même jour', () => {
  const s = initialState('2026-10-03');
  s.today.done = ['daily', 'flash:x'];
  assert.equal(rollDay(s, '2026-10-03'), s);
  assert.deepEqual(rollDay(s, TODAY).today.done, []);
});
