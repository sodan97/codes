/// <reference types="node" />
import assert from 'node:assert/strict';
import { test } from 'node:test';

import { getSubjects } from '../src/data/catalog';
import { addDays, parseDay } from '../src/lib/dates';
import { createProfile, initialState, type ProgressState } from '../src/lib/gamification';
import { formatHour, permissionState, planReminders, reminderText } from '../src/lib/reminders';

const TODAY = '2026-10-04';
const EXAM = '2027-07-01';

function awa(patch: Partial<ProgressState> = {}): ProgressState {
  const s = initialState(TODAY);
  s.profile = createProfile({ name: 'Awa', track: 'bac-s', examDate: EXAM, dailyGoal: 50, reminder: { enabled: true, hour: 19, minute: 0 } });
  return { ...s, ...patch };
}

const question = getSubjects('bac-s')[0].chapters[0].quiz[0].id;

test('reminderText : jour de l’examen et veille passent avant tout', () => {
  const s = awa({ streak: { current: 10, best: 10, lastDay: TODAY, freezes: 0 } });
  s.today.xp = 20;
  assert.equal(reminderText(s, EXAM, TODAY), 'Bonne chance pour ton Bac S aujourd’hui, Awa ! 🍀');
  assert.equal(reminderText(s, addDays(EXAM, -1), TODAY), 'Demain c’est le Bac S : relis tes essentiels et dors tôt.');
});

test('reminderText : objectif partiel, seulement pour aujourd’hui', () => {
  const s = awa();
  s.today.xp = 20;
  assert.equal(reminderText(s, TODAY, TODAY), 'Plus que 30 XP pour ton objectif du jour 💪');
  assert.notEqual(reminderText(s, addDays(TODAY, 1), TODAY), 'Plus que 30 XP pour ton objectif du jour 💪');
  // 0 XP : pas de message d'objectif.
  assert.doesNotMatch(reminderText(awa(), TODAY, TODAY) ?? '', /Plus que/);
});

test('reminderText : série, puis erreurs à revoir, puis rotation', () => {
  const streak = awa({ streak: { current: 5, best: 5, lastDay: TODAY, freezes: 0 } });
  assert.equal(reminderText(streak, addDays(TODAY, 1), TODAY), '🔥 Ta série de 5 jours t’attend, Awa ! Une seule activité suffit.');
  // Le surlendemain sans activité, la série est perdue : on ne la promet pas.
  assert.doesNotMatch(reminderText(streak, addDays(TODAY, 2), TODAY) ?? '', /série/);

  const mistakes = awa({ mistakes: { [question]: { box: 1, due: addDays(TODAY, 1), fails: 1 } } });
  assert.equal(reminderText(mistakes, TODAY, TODAY)?.startsWith('🔁'), false); // pas encore due
  assert.equal(reminderText(mistakes, addDays(TODAY, 1), TODAY), '🔁 1 question à revoir aujourd’hui (≈ 1 min)');

  // Matière masquée : ses erreurs ne sont pas comptées.
  const hidden = awa({ mistakes: { [question]: { box: 1, due: TODAY, fails: 1 } } });
  hidden.profile = { ...hidden.profile!, hiddenSubjects: [getSubjects('bac-s')[0].id] };
  assert.doesNotMatch(reminderText(hidden, TODAY, TODAY) ?? '', /à revoir/);

  const texts = new Set(Array.from({ length: 6 }, (_, i) => reminderText(awa(), addDays(TODAY, i), TODAY)));
  assert.ok(texts.has('🎯 Le défi du jour est prêt : 10 questions, +50 XP'));
  assert.ok(texts.has('⚡ 5 minutes de révision express ?'));
  // Plus de 60 jours avant l'examen : pas de message J-n.
  assert.ok(![...texts].some((t) => t?.startsWith('J-')));
  const close = awa();
  close.profile = { ...close.profile!, examDate: addDays(TODAY, 20) };
  const closeTexts = Array.from({ length: 6 }, (_, i) => reminderText(close, addDays(TODAY, i), TODAY));
  const i = closeTexts.findIndex((t) => t?.startsWith('J-'));
  assert.ok(i >= 0);
  assert.equal(closeTexts[i], `J-${20 - i} avant le Bac S : une fiche de 3 minutes ?`);
});

test('reminderText : aucun texte sans profil', () => {
  assert.equal(reminderText(initialState(TODAY), TODAY, TODAY), null);
});

test('planReminders : 7 jours au plus, aujourd’hui seulement si l’heure n’est pas passée et l’objectif pas atteint', () => {
  const morning = parseDay(TODAY);
  morning.setHours(8, 0, 0, 0);
  const plan = planReminders(awa(), morning);
  assert.equal(plan.length, 7);
  assert.equal(plan[0].id, `rappel-${TODAY}`);
  assert.equal(plan[6].day, addDays(TODAY, 6));
  assert.equal(plan[1].date.getHours(), 19);

  const evening = parseDay(TODAY);
  evening.setHours(20, 0, 0, 0);
  assert.equal(planReminders(awa(), evening)[0].day, addDays(TODAY, 1));
  assert.equal(planReminders(awa(), evening).length, 6);

  const goalReached = awa();
  goalReached.today.xp = 50;
  assert.equal(planReminders(goalReached, morning)[0].day, addDays(TODAY, 1));

  const off = awa();
  off.profile = { ...off.profile!, reminder: { enabled: false, hour: 19, minute: 0 } };
  assert.deepEqual(planReminders(off, morning), []);
  assert.equal(formatHour(7), '7 h 00');
});

test('permission : « denied » seulement si on ne peut plus la demander (Android 13 et plus)', () => {
  assert.equal(permissionState('granted', true), 'granted');
  // Jamais demandée sur Android 13 et plus : le système répond « denied », mais rien n'est bloqué.
  assert.equal(permissionState('denied', true), 'undetermined');
  assert.equal(permissionState('denied', false), 'denied');
  assert.equal(permissionState('undetermined', true), 'undetermined');
});
