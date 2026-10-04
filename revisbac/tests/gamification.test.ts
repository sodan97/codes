/// <reference types="node" />
import assert from 'node:assert/strict';
import { test } from 'node:test';

import {
  XP,
  applyGain,
  effectiveStreak,
  initialState,
  levelInfo,
  mention,
  quizGain,
  type AnswerResult,
  type ProgressState,
} from '../src/lib/gamification';

function withProfile(): ProgressState {
  const s = initialState();
  s.profile = { name: 'Awa', track: 'bac-s', examDate: '2027-07-01', dailyGoal: 50 };
  return s;
}

const gain = (s: ProgressState, xp: number, day: string) => applyGain(s, xp, [], () => {}, day).state;

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

test('quiz : XP, sans faute, défi du jour unique, erreurs mémorisées puis retirées en révision', () => {
  const ok = (id: string): AnswerResult => ({ questionId: id, subjectId: 'maths-s', correct: true });
  const ko = (id: string): AnswerResult => ({ questionId: id, subjectId: 'maths-s', correct: false });

  let r = quizGain(withProfile(), 'chapter', [ok('a'), ok('b'), ok('c'), ok('d'), ok('e')], { chapterId: 'maths-s-suites' }, '2026-10-01');
  assert.equal(r.reward.xp, 5 * XP.correctAnswer + XP.perfectQuiz + XP.dailyGoalReached);
  assert.equal(r.state.quizBest['maths-s-suites'], 100);
  assert.ok(r.reward.newBadges.some((b) => b.id === 'perfect'));

  r = quizGain(r.state, 'daily', [ok('a'), ko('x')], {}, '2026-10-01');
  assert.equal(r.state.challengesDone, 1);
  assert.equal(r.state.mistakes.x, 1);
  const again = quizGain(r.state, 'daily', [ok('a')], {}, '2026-10-01');
  assert.equal(again.state.challengesDone, 1);
  assert.equal(again.reward.xp, XP.correctAnswer);

  r = quizGain(r.state, 'review', [ok('x')], {}, '2026-10-01');
  assert.equal(r.state.mistakes.x, undefined);
});

test('examen blanc : note sur 20 et mention', () => {
  const results: AnswerResult[] = Array.from({ length: 20 }, (_, i) => ({ questionId: `q${i}`, subjectId: 'pc-s', correct: i < 15 }));
  const r = quizGain(withProfile(), 'exam', results, { subjectId: 'pc-s' }, '2026-10-01');
  assert.equal(r.state.examBest['pc-s'], 15);
  assert.equal(mention(15), 'Bien');
  assert.equal(mention(9.5), 'Insuffisant');
  assert.equal(mention(10), 'Passable');
});
