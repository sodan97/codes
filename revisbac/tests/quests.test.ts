/// <reference types="node" />
import assert from 'node:assert/strict';
import { test } from 'node:test';

import { getSubjects } from '../src/data/catalog';
import {
  CHEST_XP_MAX,
  CHEST_XP_MIN,
  MAX_FREEZES,
  XP,
  chestState,
  createProfile,
  ficheGain,
  flashcardsGain,
  initialState,
  migrate,
  openChest,
  quizGain,
  rollDay,
  withProfile,
  type AnswerResult,
  type ProgressState,
} from '../src/lib/gamification';
import { advanceQuests, generateQuests, questLink, questProgressLabel, questTitle, type Quest } from '../src/lib/quests';

const TODAY = '2026-10-04';

function bacS(examDate = '2027-07-01'): ProgressState {
  return withProfile(initialState(TODAY), createProfile({ name: 'Awa', track: 'bac-s', examDate, dailyGoal: 150 }), TODAY);
}

const quest = (kind: Quest['kind'], target: number, subjectId?: string): Quest => ({
  id: `${TODAY}:${kind}`,
  kind,
  tier: 'facile',
  target,
  progress: 0,
  done: false,
  ...(subjectId ? { subjectId } : {}),
});

const ok = (id: string, subjectId = 'maths-s'): AnswerResult => ({ questionId: id, subjectId, correct: true });
const ko = (id: string, subjectId = 'maths-s'): AnswerResult => ({ questionId: id, subjectId, correct: false });

test('quêtes du jour : une facile, une moyenne, une de variété, identiques pour un même jour et un même état', () => {
  const s = bacS();
  const quests = generateQuests(s, TODAY);
  assert.deepEqual(
    quests.map((q) => q.tier),
    ['facile', 'moyenne', 'variete'],
  );
  assert.deepEqual(generateQuests(s, TODAY), quests);
  assert.equal(new Set(quests.map((q) => q.kind)).size, 3);
  assert.ok(['correct', 'flashcards', 'fiche'].includes(quests[0].kind));
  assert.ok(['quiz80', 'combo', 'fix', 'exam'].includes(quests[1].kind));
  // Variété : un quiz dans une matière visible de l'examen.
  assert.equal(quests[2].kind, 'subject');
  assert.ok(getSubjects('bac-s').some((sub) => sub.id === quests[2].subjectId));
  assert.ok(quests.every((q) => q.progress === 0 && !q.done && q.id.startsWith(`${TODAY}:`)));
  // Sans profil : aucune quête.
  assert.deepEqual(generateQuests(initialState(TODAY), TODAY), []);
});

test('quêtes : seulement celles qui ont du sens (erreurs dues, fiches restantes, J-60)', () => {
  const far = bacS('2027-07-01');
  const near = bacS('2026-11-20');
  const subjects = getSubjects('bac-s');
  // Toutes les fiches lues : jamais « Lis une nouvelle fiche ».
  for (const sub of subjects) for (const c of sub.chapters) far.fichesRead[c.id] = TODAY;
  for (let d = 1; d <= 60; d++) {
    const day = `2026-${String(10 + Math.floor((d + 3) / 31)).padStart(2, '0')}-${String(((d + 3) % 31) + 1).padStart(2, '0')}`;
    const kinds = generateQuests(far, day).map((q) => q.kind);
    assert.ok(!kinds.includes('fiche'), day);
    assert.ok(!kinds.includes('fix'), day);
    assert.ok(!kinds.includes('exam'), day);
  }
  // À J-47 avec des erreurs dues, l'examen blanc et les erreurs à corriger finissent par être proposés.
  const qs = subjects[0].chapters[0].quiz.slice(0, 3).map((q) => q.id);
  for (const id of qs) near.mistakes[id] = { box: 1, due: TODAY, fails: 1 };
  const seen = new Set<string>();
  for (let d = 4; d <= 30; d++) for (const q of generateQuests(near, `2026-10-${String(d).padStart(2, '0')}`)) seen.add(q.kind);
  assert.ok(seen.has('exam'));
  assert.ok(seen.has('fix'));
});

test('advanceQuests : progression selon les événements, plafonnée à l’objectif', () => {
  const quests = [quest('correct', 10), quest('combo', 5), quest('quiz80', 1), quest('fix', 3), quest('subject', 1, 'pc-s'), quest('exam', 1)];
  let done = advanceQuests(quests, { kind: 'quiz', mode: 'express', correct: 4, maxCombo: 3, answered: true });
  assert.deepEqual(done, []);
  assert.deepEqual(
    quests.map(questProgressLabel),
    ['4/10', '3/5', '0/1', '0/3', '0/1', '0/1'],
  );
  done = advanceQuests(quests, { kind: 'quiz', mode: 'chapter', subjectId: 'pc-s', correct: 8, percent: 80, maxCombo: 6, mistakesFixed: 1, answered: true });
  assert.deepEqual(
    done.map((q) => q.kind),
    ['correct', 'combo', 'quiz80', 'subject'],
  );
  assert.equal(quests[0].progress, 10);
  assert.equal(questProgressLabel(quests[0]), '✓');
  // Une quête faite ne bouge plus.
  assert.deepEqual(advanceQuests(quests, { kind: 'quiz', mode: 'chapter', correct: 10, percent: 100, answered: true }), []);
  done = advanceQuests(quests, { kind: 'quiz', mode: 'exam', subjectId: 'pc-s', correct: 2, mistakesFixed: 2, answered: true });
  assert.deepEqual(
    done.map((q) => q.kind),
    ['fix', 'exam'],
  );
  // Fiches et flashcards.
  const other = [quest('fiche', 1), quest('flashcards', 1)];
  assert.deepEqual(advanceQuests(other, { kind: 'fiche' }).map((q) => q.kind), ['fiche']);
  assert.deepEqual(advanceQuests(other, { kind: 'flashcards' }).map((q) => q.kind), ['flashcards']);
});

test('intitulés et écrans des quêtes', () => {
  const s = bacS();
  const pc = getSubjects('bac-s').find((sub) => sub.id === 'pc-s')!;
  assert.equal(questTitle(quest('correct', 10)), 'Donne 10 bonnes réponses');
  assert.equal(questTitle(quest('subject', 1, 'pc-s')), `Fais un quiz en ${pc.name}`);
  assert.deepEqual(questLink(quest('fix', 3), s, TODAY), { pathname: '/quiz', params: { mode: 'review' } });
  assert.deepEqual(questLink(quest('subject', 1, 'pc-s'), s, TODAY), { pathname: '/matiere/[id]', params: { id: 'pc-s' } });
  assert.equal(questLink(quest('fiche', 1), s, TODAY).pathname, '/fiche/[id]');
  assert.deepEqual(questLink(quest('flashcards', 1), s, TODAY), { pathname: '/matieres' });
  // Examen blanc : la première matière jamais passée, sinon celle à la meilleure note la plus basse.
  const subjects = getSubjects('bac-s');
  assert.deepEqual(questLink(quest('exam', 1), s, TODAY), { pathname: '/quiz', params: { mode: 'exam', id: subjects[0].id } });
  const passed = { ...s, examBest: Object.fromEntries(subjects.map((sub, i) => [sub.id, i === 2 ? 8 : 14])) };
  assert.deepEqual(questLink(quest('exam', 1), passed, TODAY), { pathname: '/quiz', params: { mode: 'exam', id: subjects[2].id } });
});

test('+10 XP par quête accomplie, avec un message dans la récompense', () => {
  const s = bacS();
  s.today.quests = [quest('correct', 3), quest('flashcards', 1), quest('fiche', 1)];
  const r = quizGain(s, 'express', [ok('a'), ok('b'), ok('c'), ko('d')], {}, TODAY);
  assert.equal(r.reward.xp, 3 * XP.correctAnswer + XP.quest);
  assert.ok(r.reward.messages.includes(`🗺️ Quête accomplie : Donne 3 bonnes réponses : +${XP.quest} XP`));
  const f = flashcardsGain(r.state, 'maths-s', 'maths-s-suites', TODAY);
  assert.equal(f.reward.xp, XP.flashcards + XP.quest);
  assert.equal(chestState(f.state, TODAY), 'locked');
  const chapter = getSubjects('bac-s')[0].chapters[0];
  const g = ficheGain(f.state, chapter.id, 'maths-s', TODAY)!;
  assert.equal(g.reward.xp, XP.ficheRead + XP.quest);
  assert.ok(g.reward.messages.includes('🎁 Les 3 quêtes du jour sont faites : ouvre ton coffre !'));
  assert.equal(chestState(g.state, TODAY), 'ready');
});

test('combo : maxCombo du quiz, sinon plus longue suite de bonnes réponses', () => {
  const s = bacS();
  s.today.quests = [quest('combo', 3)];
  const results = [ok('a'), ok('b'), ko('c'), ok('d'), ok('e')];
  assert.equal(quizGain(s, 'express', results, {}, TODAY).state.today.quests[0].progress, 2);
  assert.equal(quizGain(s, 'express', results, { maxCombo: 3 }, TODAY).state.today.quests[0].done, true);
});

test('erreurs corrigées : une bonne réponse qui fait avancer une erreur suivie', () => {
  const s = bacS();
  s.today.quests = [quest('fix', 3)];
  s.mistakes = {
    a: { box: 1, due: TODAY, fails: 1, last: '2026-10-03' },
    b: { box: 3, due: TODAY, fails: 1, last: '2026-09-27' },
    c: { box: 1, due: TODAY, fails: 1, last: TODAY }, // créée aujourd'hui : n'avance pas
  };
  const r = quizGain(s, 'review', [ok('a'), ok('b'), ok('c'), ko('d')], {}, TODAY);
  assert.equal(r.state.today.quests[0].progress, 2);
});

test('coffre : fermé, prêt, ouvert une seule fois ; 20 à 40 XP ou un gel, déterministe par jour', () => {
  const s = bacS();
  assert.equal(openChest(s, TODAY), null);
  s.today.quests.forEach((q) => {
    q.progress = q.target;
    q.done = true;
  });
  const a = openChest(s, TODAY)!;
  const b = openChest(s, TODAY)!;
  assert.deepEqual(a.reward, b.reward);
  assert.equal(a.state.today.chestOpened, true);
  assert.equal(chestState(a.state, TODAY), 'opened');
  assert.equal(openChest(a.state, TODAY), null);
  const freeze = a.state.streak.freezes === s.streak.freezes + 1;
  if (freeze) assert.equal(a.reward.xp, 0);
  else assert.ok(a.reward.xp >= CHEST_XP_MIN && a.reward.xp <= CHEST_XP_MAX);
  // Sur plusieurs jours : des XP, et des gels tant que l'élève n'en a pas MAX_FREEZES.
  let xpDays = 0;
  let freezeDays = 0;
  for (let d = 1; d <= 28; d++) {
    const day = `2026-11-${String(d).padStart(2, '0')}`;
    const t = rollDay(s, day);
    t.today.quests.forEach((q) => {
      q.progress = q.target;
      q.done = true;
    });
    const r = openChest(t, day)!;
    if (r.state.streak.freezes > t.streak.freezes) freezeDays++;
    else {
      xpDays++;
      assert.ok(r.reward.xp >= CHEST_XP_MIN && r.reward.xp <= CHEST_XP_MAX);
    }
    // Avec déjà MAX_FREEZES gels, jamais de gel.
    assert.ok(openChest({ ...t, streak: { ...t.streak, freezes: MAX_FREEZES } }, day)!.reward.xp >= CHEST_XP_MIN);
  }
  assert.ok(xpDays > 0 && freezeDays > 0);
});

test('à minuit, les quêtes non faites disparaissent (aucun message) et le coffre se referme', () => {
  const s = bacS();
  s.today.quests[0].progress = 1;
  const next = rollDay(s, '2026-10-05');
  assert.ok(next.today.quests.every((q) => q.progress === 0 && q.id.startsWith('2026-10-05:')));
  const r = quizGain(s, 'express', [ok('a')], {}, '2026-10-05');
  assert.ok(!r.reward.messages.some((m) => /quête/i.test(m) && /échou|perdu|raté/i.test(m)));
});

test('changer d’examen tire de nouvelles quêtes ; une sauvegarde relue garde celles du jour', () => {
  const s = bacS();
  const bfm = withProfile(s, { ...s.profile!, track: 'bfm' }, TODAY);
  assert.ok(bfm.today.quests.every((q) => !q.subjectId || getSubjects('bfm').some((sub) => sub.id === q.subjectId)));
  const same = withProfile(s, { ...s.profile!, dailyGoal: 30 }, TODAY);
  assert.equal(same.today.quests, s.today.quests);
  assert.deepEqual(migrate(JSON.parse(JSON.stringify(s)), TODAY).today.quests, s.today.quests);
});

test('changer d’examen après avoir fait ses quêtes : rien n’est refait ni payé deux fois, le coffre reste prêt', () => {
  const s = bacS();
  s.today.quests = s.today.quests.map((q) => ({ ...q, progress: q.target, done: true }));
  assert.equal(chestState(s, TODAY), 'ready');
  const l = withProfile(s, { ...s.profile!, track: 'bac-l' }, TODAY);
  assert.deepEqual(l.today.quests, s.today.quests);
  assert.equal(chestState(l, TODAY), 'ready');
  const r = quizGain(l, 'express', Array.from({ length: 10 }, (_, i) => ok(`q${i}`)), {}, TODAY);
  assert.ok(!r.reward.messages.some((m) => m.startsWith('🗺️ Quête accomplie')));
});

test('quête de variété sur une matière masquée ou absente du nouvel examen : tirée à nouveau, le reste est gardé', () => {
  const s = bacS();
  const variety = s.today.quests.find((q) => q.kind === 'subject')!;
  assert.ok(variety.subjectId);
  s.today.quests[0] = { ...s.today.quests[0], progress: 1 };
  const hidden = withProfile(s, { ...s.profile!, hiddenSubjects: [variety.subjectId!] }, TODAY);
  const next = hidden.today.quests.find((q) => q.tier === 'variete');
  assert.ok(next && next.subjectId !== variety.subjectId);
  assert.deepEqual(hidden.today.quests.slice(0, 2), s.today.quests.slice(0, 2));
  // Changement d'examen avec une quête entamée : seule la quête de variété change de matière.
  const bfm = withProfile(s, { ...s.profile!, track: 'bfm' }, TODAY);
  assert.deepEqual(bfm.today.quests.slice(0, 2), s.today.quests.slice(0, 2));
  assert.ok(bfm.today.quests.every((q) => !q.subjectId || getSubjects('bfm').some((sub) => sub.id === q.subjectId)));
});
