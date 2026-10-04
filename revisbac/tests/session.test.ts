/// <reference types="node" />
import assert from 'node:assert/strict';
import { test } from 'node:test';

import { getSubjects } from '../src/data/catalog';
import { rebuildSession } from '../src/lib/quizBuilder';
import {
  EXAM_GRADE_MAX_AGE_MS,
  EXAM_GRACE_MS,
  SESSION_MAX_AGE_MS,
  ficheToContinue,
  gradeSaved,
  isPendingOtherExam,
  matchesSession,
  parseSavedSession,
  remainingMs,
  restoreSession,
  resumeButtonLabel,
  resumeCardText,
  resumeDecision,
  type SavedSession,
} from '../src/lib/resume';

const subject = getSubjects('bac-s')[0];
const chapter = subject.chapters[0];
const ids = chapter.quiz.slice(0, 4).map((q) => q.id);
// Midi heure locale : dayKey(T0) vaut '2026-10-04' quel que soit le fuseau du poste de test.
const T0 = new Date(2026, 9, 4, 12, 0).getTime();
const MIN = 60_000;

/** Quiz de chapitre : 2 réponses données sur 4, la 2e question passée (renvoyée au bout). */
function chapterSession(): SavedSession {
  return {
    mode: 'chapter',
    id: chapter.id,
    day: '2026-10-04',
    questionIds: ids,
    order: [2, 1],
    answers: [true, null, null, false],
    startedAt: T0,
    combo: 0,
    maxCombo: 1,
  };
}

function examSession(answered = 2): SavedSession {
  return {
    mode: 'exam',
    id: subject.id,
    day: '2026-10-04',
    questionIds: ids,
    order: [0, 1, 2, 3].slice(answered),
    answers: ids.map((_, i) => (i < answered ? i === 0 : null)),
    startedAt: T0,
    deadline: T0 + 3 * MIN,
    combo: 0,
    maxCombo: 1,
  };
}

test('parseSavedSession : une session relue est vérifiée', () => {
  const s = chapterSession();
  assert.deepEqual(parseSavedSession(JSON.parse(JSON.stringify(s))), s);
  assert.equal(parseSavedSession(null), null);
  assert.equal(parseSavedSession({ ...s, mode: 'duel' }), null);
  assert.equal(parseSavedSession({ ...s, id: null }), null); // quiz de chapitre sans chapitre
  assert.equal(parseSavedSession({ ...s, answers: [true] }), null);
  assert.equal(parseSavedSession({ ...s, order: [2, 2] }), null);
  assert.equal(parseSavedSession({ ...s, order: [0, 1] }), null); // question 0 déjà répondue
  assert.equal(parseSavedSession({ ...s, order: [2, 9] }), null);
  assert.equal(parseSavedSession({ ...s, startedAt: 'hier' }), null);
  assert.equal(parseSavedSession({ ...examSession(), deadline: undefined }), null);
  // Champs facultatifs complétés.
  const { combo: _c, maxCombo: _m, ...old } = s;
  assert.deepEqual(parseSavedSession(old), { ...s, combo: 0, maxCombo: 0 });
});

test('resumeDecision : reprendre pendant 24 h, sinon effacer', () => {
  const s = chapterSession();
  assert.equal(resumeDecision(s, T0 + 5 * MIN), 'resume');
  assert.equal(resumeDecision(s, T0 + SESSION_MAX_AGE_MS), 'resume');
  assert.equal(resumeDecision(s, T0 + SESSION_MAX_AGE_MS + 1), 'discard');
  // Rien de répondu : rien à reprendre.
  assert.equal(resumeDecision({ ...s, order: [0, 1, 2, 3], answers: [null, null, null, null] }, T0 + MIN), 'discard');
  // Tout est répondu (appli fermée avant l'écran de résultat) : on note.
  assert.equal(resumeDecision({ ...s, order: [], answers: [true, false, true, false] }, T0 + MIN), 'grade');
  assert.ok(matchesSession(s, 'chapter', chapter.id));
  assert.ok(!matchesSession(s, 'chapter', 'autre'));
  assert.ok(!matchesSession(s, 'exam', chapter.id));
  assert.ok(matchesSession({ ...s, mode: 'review', id: null }, 'review', undefined));
});

test('resumeDecision : examen blanc noté automatiquement quand le temps est écoulé', () => {
  const e = examSession();
  assert.equal(resumeDecision(e, T0 + MIN), 'resume');
  assert.equal(remainingMs(e, T0 + MIN), 2 * MIN);
  assert.equal(remainingMs(e, T0 + 10 * MIN), 0);
  assert.equal(remainingMs(chapterSession(), T0), null);
  // Temps écoulé appli fermée : noté à la réouverture…
  assert.equal(resumeDecision(e, e.deadline!), 'grade');
  assert.equal(resumeDecision(e, e.deadline! + EXAM_GRACE_MS), 'grade');
  // … au-delà de 30 min, noté sans proposition de reprise (même le surlendemain) ; au-delà de 7 jours, effacé.
  assert.equal(resumeDecision(e, e.deadline! + EXAM_GRACE_MS + 1), 'gradeSilently');
  assert.equal(resumeDecision(e, T0 + SESSION_MAX_AGE_MS + 1), 'gradeSilently');
  assert.equal(resumeDecision(e, T0 + 2 * SESSION_MAX_AGE_MS), 'gradeSilently');
  assert.equal(resumeDecision(e, T0 + EXAM_GRADE_MAX_AGE_MS + 1), 'discard');
  assert.equal(resumeDecision(examSession(0), T0 + SESSION_MAX_AGE_MS + 1), 'discard');
  // Examen commencé sans aucune réponse : reprise tant qu'il reste du temps, sinon rien à noter.
  assert.equal(resumeDecision(examSession(0), T0 + MIN), 'resume');
  assert.equal(resumeDecision(examSession(0), e.deadline! + MIN), 'discard');
});

test('resumeDecision : un défi du jour commencé la veille n\'est pas repris', () => {
  const daily: SavedSession = { ...chapterSession(), mode: 'daily', id: null, day: '2026-10-03', startedAt: T0 - 12 * 60 * MIN };
  // Le 04/10 à midi, le défi du 03 (commencé il y a 12 h, 2 réponses) laisse la place à celui du jour.
  assert.equal(resumeDecision(daily, T0), 'discard');
  assert.equal(resumeDecision({ ...daily, order: [], answers: [true, false, true, false] }, T0), 'discard');
  // Même jour : reprise proposée.
  assert.equal(resumeDecision({ ...daily, day: '2026-10-04', startedAt: T0 - 10 * MIN }, T0), 'resume');
  // Le jour peut être donné explicitement.
  assert.equal(resumeDecision(daily, T0, '2026-10-03'), 'resume');
});

test('isPendingOtherExam : un examen en cours n\'est pas écrasé par un autre quiz', () => {
  const e = examSession();
  assert.ok(isPendingOtherExam(e, 'daily', undefined, T0 + MIN));
  assert.ok(isPendingOtherExam(e, 'exam', 'autre-matiere', T0 + MIN));
  // Temps écoulé : toujours à noter.
  assert.ok(isPendingOtherExam(e, 'express', undefined, e.deadline! + EXAM_GRACE_MS + 1));
  // Le même examen : reprise habituelle.
  assert.ok(!isPendingOtherExam(e, 'exam', subject.id, T0 + MIN));
  // Aucune réponse, examen trop ancien ou simple quiz : rien à protéger.
  assert.ok(!isPendingOtherExam(examSession(0), 'daily', undefined, T0 + MIN));
  assert.ok(!isPendingOtherExam(e, 'daily', undefined, T0 + EXAM_GRADE_MAX_AGE_MS + 1));
  assert.ok(!isPendingOtherExam(chapterSession(), 'daily', undefined, T0 + MIN));
});

test('rebuildSession ignore les questions disparues du contenu', () => {
  const session = rebuildSession('chapter', chapter.id, [ids[0], 'id-disparu', ids[1]]);
  assert.equal(session.title, chapter.title);
  assert.equal(session.color, subject.color);
  assert.deepEqual(
    session.questions.map((q) => q.question.id),
    [ids[0], ids[1]],
  );
  const exam = rebuildSession('exam', subject.id, ids);
  assert.equal(exam.title, `Examen blanc · ${subject.name}`);
  assert.equal(exam.timeLimit, ids.length * 45);
  assert.equal(rebuildSession('review', undefined, ids).title, 'Revoir mes erreurs');
});

test('restoreSession : file et réponses recalées sur les questions restantes', () => {
  const s = { ...chapterSession(), questionIds: [ids[0], 'id-disparu', ids[2], ids[3]] };
  const r = restoreSession(s)!;
  assert.deepEqual(
    r.session.questions.map((q) => q.question.id),
    [ids[0], ids[2], ids[3]],
  );
  assert.deepEqual(r.answers, [true, null, false]);
  assert.deepEqual(r.order, [1]);
  assert.equal(r.answered, 2);
  assert.equal(r.total, 3);
  assert.equal(resumeCardText(r), `▶ Reprendre ton quiz : ${chapter.title}, 2/3`);
  assert.equal(resumeButtonLabel(r), 'Reprendre (2/3)');
  assert.equal(restoreSession({ ...s, questionIds: ['x', 'y', 'z', 'w'] }), null);
});

test('gradeSaved : les questions non traitées sont « skipped », options de finishQuiz', () => {
  const g = gradeSaved(restoreSession(examSession())!);
  assert.equal(g.mode, 'exam');
  assert.deepEqual(g.opts, { subjectId: subject.id, day: '2026-10-04', maxCombo: 1 });
  assert.deepEqual(
    g.results.map((r) => [r.correct, !!r.skipped]),
    [
      [true, false],
      [false, false],
      [false, true],
      [false, true],
    ],
  );
  assert.deepEqual(gradeSaved(restoreSession(chapterSession())!).opts, { chapterId: chapter.id, day: '2026-10-04', maxCombo: 1 });
});

test('ficheToContinue : dernière fiche ouverte, pas encore lue, visible et récente', () => {
  const visible = new Set([chapter.id]);
  const last = { chapterId: chapter.id, openedAt: T0 };
  assert.equal(ficheToContinue(last, { fichesRead: {} }, visible, T0 + MIN), chapter.id);
  assert.equal(ficheToContinue(last, { fichesRead: { [chapter.id]: '2026-10-04' } }, visible, T0 + MIN), null);
  assert.equal(ficheToContinue(last, { fichesRead: {} }, new Set(), T0 + MIN), null);
  assert.equal(ficheToContinue(last, { fichesRead: {} }, visible, T0 + 8 * 24 * 60 * MIN), null);
  assert.equal(ficheToContinue(null, { fichesRead: {} }, visible, T0), null);
});
