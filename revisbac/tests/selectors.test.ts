/// <reference types="node" />
import assert from 'node:assert/strict';
import { test } from 'node:test';

import { getSubjects } from '../src/data/catalog';
import { createProfile, initialState, type ProgressState } from '../src/lib/gamification';
import { buildQuiz, EXPRESS_SIZE, REVIEW_SIZE } from '../src/lib/quizBuilder';
import { activeMistakes, quizContext, readCount, visibleSubjects } from '../src/lib/selectors';

const TODAY = '2026-10-04';

function bacS(hiddenSubjects: string[] = []): ProgressState {
  const s = initialState(TODAY);
  s.profile = createProfile({ name: 'Awa', track: 'bac-s', examDate: '2027-07-01', dailyGoal: 50, hiddenSubjects });
  return s;
}

// Questions prises dans le catalogue réel (le contenu peut évoluer pendant la relecture).
const [first, second] = getSubjects('bac-s');
const q = (subjectIndex: 0 | 1, i: number) => (subjectIndex === 0 ? first : second).chapters.flatMap((c) => c.quiz)[i].id;
const bfmQuestion = getSubjects('bfm')[0].chapters[0].quiz[0].id;

test('visibleSubjects retire les matières masquées', () => {
  const all = visibleSubjects(bacS().profile).map((s) => s.id);
  assert.deepEqual(all, getSubjects('bac-s').map((s) => s.id));
  const visible = visibleSubjects(bacS([first.id]).profile).map((s) => s.id);
  assert.ok(!visible.includes(first.id));
  assert.equal(visible.length, all.length - 1);
  assert.deepEqual(visibleSubjects(null), []);
});

test('activeMistakes : seulement les échéances passées, du bon examen et des matières visibles', () => {
  const s = bacS();
  s.mistakes = {
    [q(0, 0)]: { box: 1, due: TODAY, fails: 1 },
    [q(0, 1)]: { box: 2, due: '2026-10-01', fails: 1 },
    [q(0, 2)]: { box: 1, due: '2026-10-01', fails: 3 },
    [q(1, 0)]: { box: 2, due: '2026-10-06', fails: 1 }, // pas encore à revoir
    [bfmQuestion]: { box: 1, due: TODAY, fails: 1 }, // autre examen
    'id-disparu': { box: 1, due: TODAY, fails: 1 }, // retiré du contenu
  };
  const { due, waiting } = activeMistakes(s, TODAY);
  // Tri : échéance la plus ancienne d'abord, puis le plus d'erreurs.
  assert.deepEqual(
    due.map((r) => r.question.id),
    [q(0, 2), q(0, 1), q(0, 0)],
  );
  assert.equal(waiting, 1);
  // Les erreurs d'une matière masquée sont conservées mais ni comptées ni proposées.
  const hidden = { ...s, profile: { ...s.profile!, hiddenSubjects: [first.id] } };
  assert.deepEqual(activeMistakes(hidden, TODAY), { due: [], waiting: 1 });
  assert.ok(hidden.mistakes[q(0, 0)]);
});

test('mode review : uniquement les erreurs dues, au plus REVIEW_SIZE', () => {
  const s = bacS();
  const ids = first.chapters.flatMap((c) => c.quiz).map((x) => x.id);
  ids.forEach((id, i) => (s.mistakes[id] = { box: 1, due: i % 2 ? '2026-10-05' : '2026-10-02', fails: 1 }));
  const session = buildQuiz('review', undefined, quizContext(s, TODAY)!);
  const dueIds = ids.filter((_, i) => i % 2 === 0);
  assert.equal(session.questions.length, Math.min(REVIEW_SIZE, dueIds.length));
  for (const r of session.questions) assert.ok(dueIds.includes(r.question.id));
});

test('mode express : 5 questions, au plus 2 par chapitre, erreurs dues en priorité', () => {
  const s = bacS();
  s.mistakes[q(1, 0)] = { box: 1, due: TODAY, fails: 1 };
  const session = buildQuiz('express', undefined, quizContext(s, TODAY)!);
  assert.equal(session.questions.length, EXPRESS_SIZE);
  assert.ok(session.questions.some((r) => r.question.id === q(1, 0)));
  const perChapter = new Map<string, number>();
  for (const r of session.questions) perChapter.set(r.chapter.id, (perChapter.get(r.chapter.id) ?? 0) + 1);
  assert.ok([...perChapter.values()].every((n) => n <= 2));
  assert.equal(new Set(session.questions.map((r) => r.question.id)).size, EXPRESS_SIZE);
});

test('readCount ne compte que les chapitres existants de l’examen', () => {
  const s = bacS();
  s.fichesRead = { [first.chapters[0].id]: TODAY, 'chapitre-disparu': TODAY, [getSubjects('bfm')[0].chapters[0].id]: TODAY };
  assert.equal(readCount(s), 1);
});
