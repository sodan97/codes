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
  assert.deepEqual(s.cards, {});
  assert.deepEqual(s.flashcardsDone, {});
  assert.deepEqual(s.quiz80Days, {});
  // Nouveau jour : compteurs du jour remis à zéro par rollDay, et quêtes du jour tirées.
  const { quests, ...today } = s.today;
  assert.deepEqual(today, { day: TODAY, xp: 0, goalBonusGiven: false, challengeDone: false, done: [], chestOpened: false });
  assert.deepEqual(
    quests.map((q) => q.tier),
    ['facile', 'moyenne', 'variete'],
  );
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
  assert.deepEqual({ ...s.today, quests: [] }, initialState(TODAY).today);
  assert.equal(s.today.quests.length, 3);
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
  assert.deepEqual({ ...s.today, quests: [] }, { day: TODAY, xp: 15, goalBonusGiven: false, challengeDone: true, done: ['daily'], quests: [], chestOpened: false });
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

/** Sauvegarde telle qu'écrite par la version 2 de l'application (avant les étoiles et les quêtes), le même jour. */
const V2 = {
  version: 2,
  profile: { name: 'Awa', track: 'bac-s', examDate: '2027-07-01', dailyGoal: 50, hiddenSubjects: [], reminder: null, haptics: true },
  xp: 320,
  streak: { current: 3, best: 5, lastDay: TODAY, freezes: 1 },
  today: { day: TODAY, xp: 40, goalBonusGiven: false, challengeDone: true, done: ['daily', 'flash:maths-s-suites'] },
  fichesRead: { 'maths-s-suites': '2026-09-30' },
  quizBest: { 'maths-s-suites': 90 },
  quizCount: 4,
  perfectCount: 0,
  challengesDone: 2,
  examBest: {},
  examHistory: {},
  examLast: {},
  mistakes: { 'maths-s-suites-q1': { box: 2, due: '2026-10-06', fails: 1, last: '2026-10-03' } },
  dailyResults: {},
  subjectsTouched: ['maths-s'],
  badges: { 'first-fiche': '2026-09-30' },
  history: { [TODAY]: 40 },
};

test('migration v2 → v3 : progression conservée, suivi des cartes vide, quêtes du jour tirées le jour même', () => {
  const s = migrate(JSON.parse(JSON.stringify(V2)), TODAY);
  assert.equal(s.version, 3);
  assert.equal(STATE_VERSION, 3);
  assert.equal(s.xp, 320);
  assert.deepEqual(s.fichesRead, V2.fichesRead);
  assert.deepEqual(s.quizBest, V2.quizBest);
  assert.deepEqual(s.mistakes, V2.mistakes);
  assert.deepEqual(s.cards, {});
  // Paquet terminé le jour de la mise à jour (noté seulement dans today.done) : il compte pour la 3e étoile.
  assert.deepEqual(s.flashcardsDone, { 'maths-s-suites': TODAY });
  // Une sauvegarde v3 n'est pas relue de cette façon (flashcardsDone y fait foi).
  assert.deepEqual(migrate({ ...JSON.parse(JSON.stringify(V2)), version: 3 }, TODAY).flashcardsDone, {});
  // Un paquet d'un autre jour reste noté à sa date.
  assert.deepEqual(
    migrate({ ...JSON.parse(JSON.stringify(V2)), today: { ...V2.today, day: '2026-10-02' } }, TODAY).flashcardsDone,
    { 'maths-s-suites': '2026-10-02' },
  );
  // Le jour d'un quiz réussi n'était pas connu en v2 : rien n'est inventé.
  assert.deepEqual(s.quiz80Days, {});
  // Même jour : compteurs gardés, quêtes ajoutées (coffre fermé).
  assert.equal(s.today.xp, 40);
  assert.deepEqual(s.today.done, ['daily', 'flash:maths-s-suites']);
  assert.equal(s.today.chestOpened, false);
  assert.equal(s.today.quests.length, 3);
  assert.ok(s.today.quests.every((q) => q.progress === 0 && !q.done && q.id.startsWith(`${TODAY}:`)));
});

test('migration v3 : cartes, jours de quiz, quêtes et coffre relus et vérifiés', () => {
  const once = migrate(JSON.parse(JSON.stringify(V2)), TODAY);
  once.cards = { 'maths-s-suites#definition': { box: 3, due: '2026-10-11', last: TODAY } };
  once.flashcardsDone = { 'maths-s-suites': TODAY };
  once.quiz80Days = { 'maths-s-suites': ['2026-10-01', TODAY] };
  once.today.quests[0].progress = once.today.quests[0].target;
  once.today.quests[0].done = true;
  once.today.chestOpened = true;
  // Relue telle quelle : même état.
  assert.deepEqual(migrate(JSON.parse(JSON.stringify(once)), TODAY), once);

  const broken = migrate(
    {
      ...JSON.parse(JSON.stringify(once)),
      cards: { a: { box: 9, due: 'demain' }, b: 'x' },
      flashcardsDone: { a: 'hier', b: TODAY },
      quiz80Days: { a: [TODAY, TODAY, 'x', '2026-10-01'], b: 'x', c: [] },
      today: { ...once.today, quests: [{ id: 'q', kind: 'inconnue', tier: 'facile', target: 1 }, { id: 'r', kind: 'fix', tier: 'moyenne', target: 3, progress: 7 }] },
    },
    TODAY,
  );
  assert.deepEqual(broken.cards, { a: { box: 5, due: TODAY, last: TODAY } });
  assert.deepEqual(broken.flashcardsDone, { b: TODAY });
  assert.deepEqual(broken.quiz80Days, { a: ['2026-10-01', TODAY] });
  assert.deepEqual(broken.today.quests, [{ id: 'r', kind: 'fix', tier: 'moyenne', target: 3, progress: 3, done: true }]);
});

test('rollDay : nouvelles quêtes au changement de jour, coffre refermé ; rien ne change le même jour', () => {
  const s = migrate(JSON.parse(JSON.stringify(V2)), TODAY);
  s.today.quests.forEach((q) => {
    q.progress = q.target;
    q.done = true;
  });
  s.today.chestOpened = true;
  assert.equal(rollDay(s, TODAY), s);
  const next = rollDay(s, '2026-10-05');
  assert.equal(next.today.chestOpened, false);
  assert.equal(next.today.quests.length, 3);
  assert.ok(next.today.quests.every((q) => !q.done && q.id.startsWith('2026-10-05:')));
});
