/// <reference types="node" />
import assert from 'node:assert/strict';
import { test } from 'node:test';

import { getSubject, getSubjects, isOptional } from '../src/data/catalog';
import type { TrackId } from '../src/data/types';
import { createProfile, initialState, type ProgressState } from '../src/lib/gamification';
import { buildQuiz, DAILY_SIZE, EXAM_SIZE, examQuotas, EXPRESS_SIZE, type QuizSession } from '../src/lib/quizBuilder';
import { quizContext, visibleSubjects } from '../src/lib/selectors';

const TODAY = '2026-10-04';

function student(track: TrackId, hiddenSubjects: string[] = []): ProgressState {
  const s = initialState(TODAY);
  s.profile = createProfile({ name: 'Awa', track, examDate: '2027-07-01', dailyGoal: 50, hiddenSubjects });
  return s;
}

const ctxOf = (s: ProgressState, today = TODAY) => quizContext(s, today)!;
const ids = (session: QuizSession) => session.questions.map((r) => r.question.id);
const perChapter = (session: QuizSession) => {
  const m = new Map<string, number>();
  for (const r of session.questions) m.set(r.chapter.id, (m.get(r.chapter.id) ?? 0) + 1);
  return m;
};

// ---------- Examen blanc (P1-2) ----------

test('examQuotas : proportionnel, au moins 1 par chapitre, total = min(20, questions)', () => {
  assert.deepEqual(examQuotas([9, 9, 9]).reduce((a, b) => a + b), 20);
  assert.ok(examQuotas([9, 9, 9]).every((n) => n >= 6 && n <= 7));
  assert.deepEqual(examQuotas([30, 2]), [19, 1]);
  assert.deepEqual(examQuotas([3, 4]), [3, 4]); // petite banque : tout est pris
  // Plus de chapitres que de questions à tirer : 20 chapitres gardés, une question chacun.
  const many = examQuotas(Array.from({ length: 25 }, () => 6));
  assert.equal(many.reduce((a, b) => a + b), 20);
  assert.ok(many.every((n) => n <= 1));
});

test('examen blanc : chaque chapitre a au moins une question, 20 au total, sans doublon', () => {
  for (const track of ['bfm', 'bac-s', 'bac-l'] as const) {
    const ctx = ctxOf(student(track));
    for (const subject of getSubjects(track)) {
      const session = buildQuiz('exam', subject.id, ctx);
      const nQuestions = subject.chapters.reduce((n, c) => n + c.quiz.length, 0);
      assert.equal(session.questions.length, Math.min(EXAM_SIZE, nQuestions), subject.id);
      assert.equal(new Set(ids(session)).size, session.questions.length, subject.id);
      assert.equal(session.timeLimit, session.questions.length * 45);
      const counts = perChapter(session);
      if (subject.chapters.length <= EXAM_SIZE) for (const c of subject.chapters) assert.ok((counts.get(c.id) ?? 0) >= 1, `${c.id} absent`);
    }
  }
});

test('examen blanc : le tirage suivant privilégie les questions non vues au dernier examen', () => {
  const s = student('bfm');
  const subject = getSubject('espagnol-bfm')!;
  const all = subject.chapters.flatMap((c) => c.quiz.map((q) => q.id));
  const first = buildQuiz('exam', subject.id, ctxOf(s));
  s.examLast[subject.id] = ids(first);
  const unseen = all.filter((id) => !s.examLast[subject.id].includes(id));
  assert.equal(unseen.length, all.length - first.questions.length); // 27 - 20 = 7 aujourd'hui
  for (let i = 0; i < 10; i++) {
    const second = ids(buildQuiz('exam', subject.id, ctxOf(s)));
    for (const id of unseen) assert.ok(second.includes(id), `${id} non vue mais absente du second tirage`);
  }
});

// ---------- Défi du jour (P1-10) ----------

test('défi du jour : jamais de LV2, identique quelles que soient les matières masquées', () => {
  for (const track of ['bfm', 'bac-s', 'bac-l'] as const) {
    const optional = getSubjects(track).filter((s) => isOptional(s.id)).map((s) => s.id);
    const someHidden = [...optional, getSubjects(track).find((s) => !isOptional(s.id))!.id];
    for (const day of ['2026-10-04', '2026-10-05', '2026-11-20', '2027-03-01']) {
      const a = buildQuiz('daily', undefined, ctxOf(student(track), day));
      const b = buildQuiz('daily', undefined, ctxOf(student(track, someHidden), day));
      assert.equal(a.questions.length, DAILY_SIZE);
      assert.deepEqual(ids(a), ids(b));
      assert.ok(a.questions.every((r) => !r.subject.id.startsWith('espagnol-')));
    }
  }
});

// ---------- Révision express (P1-9) ----------

function checkExpress(session: QuizSession) {
  assert.equal(session.questions.length, EXPRESS_SIZE);
  assert.equal(new Set(ids(session)).size, EXPRESS_SIZE);
  assert.ok([...perChapter(session).values()].every((n) => n <= 2));
}

test('express : 5 questions, au plus 2 par chapitre, 2 erreurs dues au plus, les plus anciennes d’abord', () => {
  const s = student('bac-s');
  const [a, b] = getSubjects('bac-s');
  const q = (subject: typeof a, chapter: number, i: number) => subject.chapters[chapter].quiz[i].id;
  s.mistakes = {
    [q(a, 2, 0)]: { box: 1, due: '2026-10-03', fails: 1 }, // 3e erreur due : au-delà de la limite de 2
    [q(a, 1, 0)]: { box: 1, due: '2026-10-01', fails: 1 }, // la plus ancienne
    [q(b, 0, 0)]: { box: 2, due: '2026-10-02', fails: 2 },
    [q(b, 1, 0)]: { box: 1, due: '2026-10-09', fails: 1 }, // pas encore due
  };
  for (let i = 0; i < 10; i++) {
    const session = buildQuiz('express', undefined, ctxOf(s));
    checkExpress(session);
    const fromMistakes = ids(session).filter((id) => s.mistakes[id]);
    assert.deepEqual(fromMistakes.sort(), [q(a, 1, 0), q(b, 0, 0)].sort());
  }
});

test('express : chapitres lus et encore fragiles en priorité', () => {
  const s = student('bac-s');
  const [a] = getSubjects('bac-s');
  // Trois chapitres lus : deux fragiles, un bien maîtrisé.
  for (const c of a.chapters.slice(0, 3)) s.fichesRead[c.id] = '2026-10-01';
  s.quizBest[a.chapters[0].id] = 50;
  s.quizBest[a.chapters[2].id] = 100;
  for (let i = 0; i < 10; i++) {
    const session = buildQuiz('express', undefined, ctxOf(s));
    checkExpress(session);
    const counts = perChapter(session);
    assert.equal(counts.get(a.chapters[0].id), 2);
    assert.equal(counts.get(a.chapters[1].id), 2);
    assert.equal(counts.get(a.chapters[2].id), 1);
  }
});

test('express : seulement des matières visibles, même avec des erreurs dans une matière masquée', () => {
  const [hidden] = getSubjects('bac-s');
  const s = student('bac-s', [hidden.id]);
  s.mistakes[hidden.chapters[0].quiz[0].id] = { box: 1, due: TODAY, fails: 3 };
  s.fichesRead[hidden.chapters[0].id] = '2026-10-01';
  const visible = new Set(visibleSubjects(s.profile).map((x) => x.id));
  for (let i = 0; i < 10; i++) {
    const session = buildQuiz('express', undefined, ctxOf(s));
    checkExpress(session);
    assert.ok(session.questions.every((r) => visible.has(r.subject.id)));
  }
});

test('express : aucune fiche lue → questions des premiers chapitres de chaque matière', () => {
  const s = student('bfm');
  const firstChapters = new Set(getSubjects('bfm').map((x) => x.chapters[0].id));
  for (let i = 0; i < 10; i++) {
    const session = buildQuiz('express', undefined, ctxOf(s));
    checkExpress(session);
    assert.ok(session.questions.every((r) => firstChapters.has(r.chapter.id)));
  }
});
