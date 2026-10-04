/// <reference types="node" />
import assert from 'node:assert/strict';
import { test } from 'node:test';

import {
  XP,
  applyAnswerToMistakes,
  applyGain,
  createProfile,
  effectiveStreak,
  examNote,
  flashcardsGain,
  formatNote,
  initialState,
  levelInfo,
  mention,
  migrate,
  noteTrend,
  quizGain,
  rollDay,
  type AnswerResult,
  type ProgressState,
  type Reward,
} from '../src/lib/gamification';
import { dailyShareText, examShareText } from '../src/lib/shareText';

function withProfile(examDate = '2027-07-01'): ProgressState {
  const s = initialState('2026-10-01');
  s.profile = createProfile({ name: 'Awa', track: 'bac-s', examDate, dailyGoal: 50 });
  return s;
}

const ok = (id: string): AnswerResult => ({ questionId: id, subjectId: 'maths-s', correct: true });
const ko = (id: string): AnswerResult => ({ questionId: id, subjectId: 'maths-s', correct: false });
const skip = (id: string): AnswerResult => ({ questionId: id, subjectId: 'maths-s', correct: false, skipped: true });

const gain = (s: ProgressState, xp: number, day: string) => applyGain(s, xp, [], () => {}, day).state;
/** XP hors quêtes du jour (celles-ci sont testées dans quests.test.ts). */
const baseXp = (r: { reward: Reward }) => r.reward.xp - r.reward.messages.filter((m) => m.startsWith('🗺️')).length * XP.quest;

test('niveaux : seuils 0 / 100 / 300 / 600', () => {
  assert.equal(levelInfo(0).level, 1);
  assert.equal(levelInfo(99).level, 1);
  assert.equal(levelInfo(100).level, 2);
  assert.equal(levelInfo(300).level, 3);
  assert.equal(levelInfo(599).toNext, 1);
});

test('la série augmente sur des jours consécutifs et repart à 1 après un trou', () => {
  let s = withProfile();
  s = gain(s, 10, '2026-10-01');
  s = gain(s, 10, '2026-10-01');
  assert.equal(s.streak.current, 1);
  s = gain(s, 10, '2026-10-02');
  s = gain(s, 10, '2026-10-03');
  assert.equal(s.streak.current, 3);
  assert.equal(effectiveStreak(s, '2026-10-04'), 3);
  assert.equal(effectiveStreak(s, '2026-10-05'), 0);
  s = gain(s, 10, '2026-10-06');
  assert.equal(s.streak.current, 1);
  assert.equal(s.streak.best, 3);
});

test('7 jours d’affilée donnent un gel qui sauve la série après un jour manqué', () => {
  let s = withProfile();
  for (let d = 1; d <= 7; d++) s = gain(s, 10, `2026-10-0${d}`);
  assert.equal(s.streak.current, 7);
  assert.equal(s.streak.freezes, 1);
  assert.equal(effectiveStreak(s, '2026-10-09'), 7); // le 8 est manqué, couvert par le gel
  s = gain(s, 10, '2026-10-09');
  assert.equal(s.streak.current, 8);
  assert.equal(s.streak.freezes, 0);
});

test('bonus d’objectif du jour une seule fois, et compteurs remis à zéro le lendemain', () => {
  let s = withProfile();
  let r = applyGain(s, 40, [], () => {}, '2026-10-01');
  assert.equal(r.reward.xp, 40);
  r = applyGain(r.state, 20, [], () => {}, '2026-10-01');
  assert.equal(r.reward.xp, 20 + XP.dailyGoalReached);
  r = applyGain(r.state, 60, [], () => {}, '2026-10-01');
  assert.equal(r.reward.xp, 60);
  r = applyGain(r.state, 10, [], () => {}, '2026-10-02');
  assert.equal(r.state.today.xp, 10);
  assert.equal(r.state.today.goalBonusGiven, false);
});

test('quiz : XP, sans faute, défi du jour unique, erreurs mémorisées', () => {
  let r = quizGain(withProfile(), 'chapter', [ok('a'), ok('b'), ok('c'), ok('d'), ok('e')], { chapterId: 'maths-s-suites' }, '2026-10-01');
  assert.equal(baseXp(r), 5 * XP.correctAnswer + XP.perfectQuiz + XP.dailyGoalReached);
  assert.equal(r.state.quizBest['maths-s-suites'], 100);
  assert.ok(r.reward.newBadges.some((b) => b.id === 'perfect'));

  r = quizGain(r.state, 'daily', [ok('a'), ko('x')], {}, '2026-10-01');
  assert.equal(r.state.challengesDone, 1);
  assert.deepEqual(r.state.mistakes.x, { box: 1, due: '2026-10-02', fails: 1, last: '2026-10-01' });
});

test('examen blanc : note sur 20 et mention', () => {
  const results: AnswerResult[] = Array.from({ length: 20 }, (_, i) => ({ questionId: `q${i}`, subjectId: 'pc-s', correct: i < 15 }));
  const r = quizGain(withProfile(), 'exam', results, { subjectId: 'pc-s' }, '2026-10-01');
  assert.equal(r.state.examBest['pc-s'], 15);
  assert.deepEqual(r.state.examHistory['pc-s'], [{ day: '2026-10-01', note: 15 }]);
  assert.equal(r.state.examLast['pc-s'].length, 20);
  assert.equal(formatNote(13.5), '13,5');
  assert.equal(formatNote(12), '12');
  assert.equal(mention(15), 'Bien');
  assert.equal(mention(9.5), 'Insuffisant');
  assert.equal(mention(10), 'Passable');
});

// ---------- Leitner (P1-4) ----------

test('Leitner : boîte 1 → 2 → 3 → maîtrisée, avec des intervalles de 1, 3 et 7 jours', () => {
  let s = quizGain(withProfile(), 'chapter', [ko('q')], { chapterId: 'c' }, '2026-10-01').state;
  assert.deepEqual(s.mistakes.q, { box: 1, due: '2026-10-02', fails: 1, last: '2026-10-01' });
  s = quizGain(s, 'review', [ok('q')], {}, '2026-10-02').state;
  assert.deepEqual(s.mistakes.q, { box: 2, due: '2026-10-05', fails: 1, last: '2026-10-02' });
  s = quizGain(s, 'review', [ok('q')], {}, '2026-10-05').state;
  assert.deepEqual(s.mistakes.q, { box: 3, due: '2026-10-12', fails: 1, last: '2026-10-05' });
  const r = quizGain(s, 'review', [ok('q')], {}, '2026-10-12');
  assert.equal(r.state.mistakes.q, undefined);
  assert.ok(r.reward.messages.includes('✅ 1 question maîtrisée : elle sort de ta liste'));
});

test('Leitner : une erreur remet en boîte 1, une bonne réponse dans un autre mode fait avancer avant l’échéance', () => {
  const s = withProfile();
  s.mistakes.q = { box: 3, due: '2026-10-08', fails: 2 };
  applyAnswerToMistakes(s, ko('q'), 'daily', '2026-10-03');
  assert.deepEqual(s.mistakes.q, { box: 1, due: '2026-10-04', fails: 3, last: '2026-10-03' });
  // Entrée sans `last` (sauvegarde plus ancienne), échue le 10-08 : avance dès le 10-03.
  s.mistakes.r = { box: 1, due: '2026-10-08', fails: 1 };
  applyAnswerToMistakes(s, ok('r'), 'chapter', '2026-10-03');
  assert.deepEqual(s.mistakes.r, { box: 2, due: '2026-10-06', fails: 1, last: '2026-10-03' });
  // Une bonne réponse sur une question jamais ratée ne crée rien.
  assert.equal(applyAnswerToMistakes(s, ok('autre'), 'chapter', '2026-10-03'), 'none');
  assert.equal(s.mistakes.autre, undefined);
});

test('Leitner : au plus une avance par jour (rejouer la même série ne fait pas sortir une erreur)', () => {
  let s = quizGain(withProfile(), 'daily', [ko('q')], {}, '2026-10-01').state;
  for (let i = 0; i < 3; i++) s = quizGain(s, 'daily', [ok('q')], {}, '2026-10-01').state;
  assert.deepEqual(s.mistakes.q, { box: 1, due: '2026-10-02', fails: 1, last: '2026-10-01' });
  for (let i = 0; i < 3; i++) s = quizGain(s, 'chapter', [ok('q')], { chapterId: 'c' }, '2026-10-02').state;
  assert.equal(s.mistakes.q.box, 2);
  s = quizGain(s, 'chapter', [ok('q')], { chapterId: 'c' }, '2026-10-03').state;
  assert.equal(s.mistakes.q.box, 3);
  s = quizGain(s, 'chapter', [ok('q')], { chapterId: 'c' }, '2026-10-04').state;
  assert.equal(s.mistakes.q, undefined);
});

test('Leitner : le champ `last` est conservé à la lecture d’une sauvegarde', () => {
  const s = withProfile();
  s.mistakes.q = { box: 2, due: '2026-10-05', fails: 1, last: '2026-10-02' };
  s.mistakes.r = { box: 1, due: '2026-10-05', fails: 1 };
  const m = migrate(JSON.parse(JSON.stringify(s)), '2026-10-03').mistakes;
  assert.deepEqual(m.q, { box: 2, due: '2026-10-05', fails: 1, last: '2026-10-02' });
  assert.deepEqual(m.r, { box: 1, due: '2026-10-05', fails: 1 });
});

test('Leitner : en examen blanc, la boîte 3 reste en boîte 3 (jamais de suppression)', () => {
  const s = withProfile();
  s.mistakes.q = { box: 3, due: '2026-10-01', fails: 1 };
  assert.equal(applyAnswerToMistakes(s, ok('q'), 'exam', '2026-10-01'), 'promoted');
  assert.deepEqual(s.mistakes.q, { box: 3, due: '2026-10-08', fails: 1, last: '2026-10-01' });
});

test('Leitner : les questions non traitées (skipped) sont ignorées', () => {
  const s = withProfile();
  s.mistakes.q = { box: 2, due: '2026-10-01', fails: 1 };
  const r = quizGain(s, 'exam', [skip('q'), skip('n'), ok('a')], { subjectId: 'maths-s' }, '2026-10-01');
  assert.deepEqual(r.state.mistakes.q, { box: 2, due: '2026-10-01', fails: 1 });
  assert.equal(r.state.mistakes.n, undefined);
  assert.equal(r.state.examHistory['maths-s'][0].note, 6.5); // 1/3 → 6,67 arrondi au demi-point
});

test('Leitner : échéance plafonnée à examDate - 2, jamais avant demain', () => {
  const s = withProfile('2026-10-08');
  s.mistakes.q = { box: 2, due: '2026-10-01', fails: 1 };
  applyAnswerToMistakes(s, ok('q'), 'review', '2026-10-01'); // J+7 = 10-08 → plafonné au 10-06
  assert.equal(s.mistakes.q.due, '2026-10-06');
  const t = withProfile('2026-10-02');
  applyAnswerToMistakes(t, ko('q'), 'review', '2026-10-01'); // examDate - 2 = hier → demain
  assert.equal(t.mistakes.q.due, '2026-10-02');
});

// ---------- XP honnête (P1-6) ----------

test('défi du jour rejoué : 0 XP, compteurs et score officiel figés', () => {
  const first = quizGain(withProfile(), 'daily', [ok('a'), ko('b'), ok('c'), ok('d'), ok('e')], {}, '2026-10-01');
  assert.deepEqual(first.state.dailyResults['2026-10-01'], { correct: 4, total: 5, grid: '✅❌✅✅✅' });
  assert.ok(first.state.today.done.includes('daily'));
  const again = quizGain(first.state, 'daily', [ok('a'), ok('b'), ok('c'), ok('d'), ok('e')], {}, '2026-10-01');
  assert.equal(again.reward.xp, 0);
  assert.ok(again.reward.messages.includes('🔁 Entraînement : pas d’XP pour un défi déjà relevé'));
  assert.equal(again.state.quizCount, first.state.quizCount);
  assert.equal(again.state.perfectCount, first.state.perfectCount);
  assert.equal(again.state.challengesDone, 1);
  assert.deepEqual(again.state.dailyResults['2026-10-01'], first.state.dailyResults['2026-10-01']);
  // Erreur créée le jour même : elle n'avance pas avant demain (une avance par jour au plus).
  assert.equal(again.state.mistakes.b.box, 1);
  const tomorrow = quizGain(again.state, 'daily', [ok('a'), ok('b'), ok('c'), ok('d'), ok('e')], {}, '2026-10-02');
  assert.equal(tomorrow.state.mistakes.b.box, 2);
});

test('quiz de chapitre refait le même jour : XP divisés par 2, sans bonus ni compteurs', () => {
  const five = [ok('a'), ok('b'), ok('c'), ok('d'), ok('e')];
  const start = withProfile();
  start.profile!.dailyGoal = 150; // pas de bonus d'objectif dans ce test
  const first = quizGain(start, 'chapter', [ok('a'), ko('b'), ok('c'), ok('d'), ok('e')], { chapterId: 'ch' }, '2026-10-01');
  const again = quizGain(first.state, 'chapter', five, { chapterId: 'ch' }, '2026-10-01');
  assert.equal(baseXp(again), Math.floor((5 * XP.correctAnswer) / 2));
  assert.ok(again.reward.messages.includes('XP réduits de moitié : tu as déjà fait ce quiz aujourd’hui'));
  assert.equal(again.state.quizCount, first.state.quizCount);
  assert.equal(again.state.perfectCount, 0);
  assert.equal(again.state.quizBest.ch, 100);
  // Le lendemain, barème normal.
  const nextDay = quizGain(again.state, 'chapter', five, { chapterId: 'ch' }, '2026-10-02');
  assert.equal(baseXp(nextDay), 5 * XP.correctAnswer + XP.perfectQuiz);
});

test('examen blanc : bonus une fois par jour et par matière, historique limité à 5 notes', () => {
  const results: AnswerResult[] = Array.from({ length: 4 }, (_, i) => ({ questionId: `q${i}`, subjectId: 'pc-s', correct: i < 2 }));
  let s = withProfile();
  s.profile!.dailyGoal = 150; // pas de bonus d'objectif dans ce test
  const first = quizGain(s, 'exam', results, { subjectId: 'pc-s' }, '2026-10-01');
  assert.equal(first.reward.xp, 2 * XP.correctAnswer + XP.mockExam);
  const again = quizGain(first.state, 'exam', results, { subjectId: 'pc-s' }, '2026-10-01');
  assert.equal(again.reward.xp, 2 * XP.correctAnswer);
  assert.equal(again.state.quizCount, first.state.quizCount);
  s = again.state;
  for (let d = 2; d <= 6; d++) s = quizGain(s, 'exam', results, { subjectId: 'pc-s' }, `2026-10-0${d}`).state;
  assert.equal(s.examHistory['pc-s'].length, 5);
  assert.equal(s.examHistory['pc-s'][4].day, '2026-10-06');
});

test('flashcards : +5 XP une fois par jour et par chapitre', () => {
  const first = flashcardsGain(withProfile(), 'maths-s', 'ch', '2026-10-01');
  assert.equal(baseXp(first), XP.flashcards);
  const again = flashcardsGain(first.state, 'maths-s', 'ch', '2026-10-01');
  assert.equal(again.reward.xp, 0);
  assert.deepEqual(again.reward.messages, ['Paquet déjà terminé aujourd’hui : pas d’XP, mais bravo pour la révision !']);
  assert.equal(baseXp(flashcardsGain(again.state, 'maths-s', 'ch', '2026-10-02')), XP.flashcards);
});

test('toute activité terminée prolonge la série, même avec 0 bonne réponse', () => {
  let s = gain(withProfile(), 10, '2026-10-01');
  s = quizGain(s, 'chapter', [ko('a'), ko('b')], { chapterId: 'ch' }, '2026-10-02').state;
  assert.equal(s.streak.current, 2);
  assert.equal(s.streak.lastDay, '2026-10-02');
  // Sans XP ni activité terminée, la série n'avance pas.
  s = applyGain(s, 0, [], () => {}, '2026-10-03').state;
  assert.equal(s.streak.lastDay, '2026-10-02');
  s = applyGain(s, 0, [], () => {}, '2026-10-03', { activity: true }).state;
  assert.equal(s.streak.current, 3);
});

test('today.done est vidé au changement de jour', () => {
  const s = flashcardsGain(withProfile(), 'maths-s', 'ch', '2026-10-01').state;
  assert.deepEqual(s.today.done, ['flash:ch']);
  assert.deepEqual(rollDay(s, '2026-10-02').today.done, []);
});

test('score du défi : écrit une seule fois et gardé 30 jours', () => {
  let s = quizGain(withProfile(), 'daily', [ok('a')], {}, '2026-09-01').state;
  assert.ok(s.dailyResults['2026-09-01']);
  s = quizGain(s, 'daily', [ok('a')], {}, '2026-10-05').state;
  assert.equal(s.dailyResults['2026-09-01'], undefined);
  assert.ok(s.dailyResults['2026-10-05']);
});

// ---------- Examen blanc : questions non traitées, historique (P1-2) ----------

test('examen blanc : les questions non traitées comptent 0 et ne vont jamais dans les erreurs', () => {
  assert.equal(examNote(15, 20), 15);
  assert.equal(examNote(7, 12), 11.5); // arrondi au demi-point
  // 12 bonnes, 3 fausses, 5 non traitées (fin du temps) : note sur le total de 20.
  const results: AnswerResult[] = [
    ...Array.from({ length: 12 }, (_, i) => ok(`ok${i}`)),
    ...Array.from({ length: 3 }, (_, i) => ko(`ko${i}`)),
    ...Array.from({ length: 5 }, (_, i) => skip(`skip${i}`)),
  ];
  const s = quizGain(withProfile(), 'exam', results, { subjectId: 'maths-s' }, '2026-10-01').state;
  assert.equal(s.examBest['maths-s'], 12);
  assert.deepEqual(Object.keys(s.mistakes).sort(), ['ko0', 'ko1', 'ko2']);
  assert.equal(s.examLast['maths-s'].length, 20);
});

test('examen blanc sans aucune réponse (chrono expiré) : ni XP, ni série, ni note', () => {
  const results = Array.from({ length: 20 }, (_, i) => skip(`s${i}`));
  const none = quizGain(withProfile(), 'exam', results, { subjectId: 'maths-s' }, '2026-10-01');
  assert.equal(none.reward.xp, 0);
  assert.ok(none.reward.messages.includes('📝 Examen non traité : pas de note ni de bonus'));
  assert.deepEqual(none.state.streak, withProfile().streak);
  assert.equal(none.state.quizCount, 0);
  assert.equal(none.state.badges['first-quiz'], undefined);
  assert.deepEqual(none.state.examBest, {});
  assert.deepEqual(none.state.examHistory, {});
  assert.deepEqual(none.state.examLast, {});
  assert.deepEqual(none.state.mistakes, {});
  assert.ok(!none.state.today.done.includes('exam:maths-s'));
  assert.equal(none.state.history['2026-10-01'], undefined);
  // Un vrai examen ensuite le même jour reçoit encore le bonus.
  const real = quizGain(none.state, 'exam', [ok('a'), ...results.slice(1)], { subjectId: 'maths-s' }, '2026-10-01');
  assert.equal(real.reward.xp, XP.correctAnswer + XP.mockExam);
  assert.equal(real.state.examBest['maths-s'], 1);
});

test('défi du jour commencé la veille (avant minuit) : entraînement, le défi du jour reste à faire', () => {
  const r = quizGain(withProfile(), 'daily', [ok('a'), ok('b')], { day: '2026-10-04' }, '2026-10-05');
  assert.equal(r.reward.xp, 0); // comme un rejeu
  assert.deepEqual(r.reward.messages, ['🔁 Entraînement : ce défi date d’hier, pas d’XP']);
  assert.equal(r.state.dailyResults['2026-10-05'], undefined);
  assert.equal(r.state.dailyResults['2026-10-04'], undefined);
  assert.equal(r.state.today.challengeDone, false);
  assert.equal(r.state.challengesDone, 0);
  // La série avance quand même (des réponses ont été données).
  assert.equal(r.state.streak.lastDay, '2026-10-05');
  // Commencé et fini le même jour : défi officiel.
  const same = quizGain(withProfile(), 'daily', [ok('a')], { day: '2026-10-05' }, '2026-10-05');
  assert.ok(same.state.dailyResults['2026-10-05']);
});

test('historique : un jour d’activité à 0 XP est enregistré, pas un gain nul sans activité', () => {
  let s = gain(withProfile(), 10, '2026-10-01');
  s = quizGain(s, 'chapter', [ko('a')], { chapterId: 'ch' }, '2026-10-02').state;
  assert.equal(s.history['2026-10-02'], 0);
  assert.ok('2026-10-02' in s.history);
  s = applyGain(s, 0, [], () => {}, '2026-10-03').state;
  assert.ok(!('2026-10-03' in s.history));
});

test('noteTrend : « 8 → 11,5 → 13 » à partir de 2 notes', () => {
  assert.equal(noteTrend(undefined), null);
  assert.equal(noteTrend([{ day: '2026-10-01', note: 8 }]), null);
  assert.equal(
    noteTrend([
      { day: '2026-10-01', note: 8 },
      { day: '2026-10-02', note: 11.5 },
      { day: '2026-10-03', note: 13 },
    ]),
    '8 → 11,5 → 13',
  );
});

// ---------- Partage (P1-5) ----------

test('texte de partage du défi : ligne série seulement à partir de 2 jours, aucun lien', () => {
  const base = { day: '2026-10-04', trackLabel: 'Bac S', correct: 8, total: 10, grid: '✅✅❌✅✅✅✅❌✅✅' };
  const one = dailyShareText({ ...base, streak: 1 });
  assert.ok(one.startsWith('RéviBac · Défi du jour du '));
  assert.ok(one.includes('(Bac S)\n8/10 ✅✅❌✅✅✅✅❌✅✅\n'));
  assert.ok(!one.includes('Série'));
  assert.ok(one.endsWith('Fais le même défi que moi dans RéviBac !'));
  assert.ok(!/https?:/.test(one));
  assert.ok(dailyShareText({ ...base, streak: 3 }).includes('\n🔥 Série : 3 jours\n'));
});

test('texte de partage de l’examen blanc : note à la française et mention', () => {
  assert.equal(
    examShareText({ subjectName: 'Philosophie', note: 13.5 }),
    'RéviBac · Examen blanc de Philosophie : 13,5/20 (mention Assez-Bien).\nQuestions de révision RéviBac, pas un sujet officiel.',
  );
});

test('horloge qui recule : la série et les compteurs du jour sont conservés', () => {
  let s = withProfile();
  s = gain(s, 10, '2026-10-01');
  s = gain(s, 10, '2026-10-02');
  s = applyGain(s, 60, [], () => {}, '2026-10-02').state;
  assert.equal(s.streak.current, 2);
  assert.equal(s.today.goalBonusGiven, true);
  const back = applyGain(s, 10, [], () => {}, '2026-10-01').state;
  assert.equal(back.streak.current, 2);
  assert.equal(back.streak.lastDay, '2026-10-02');
  assert.equal(back.today.day, '2026-10-02');
  assert.equal(back.today.goalBonusGiven, true);
  assert.equal(rollDay(back, '2026-09-30'), back);
});
