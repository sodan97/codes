/// <reference types="node" />
import assert from 'node:assert/strict';
import { test } from 'node:test';

import { getSubjects } from '../src/data/catalog';
import { BADGES, cardReview, createProfile, flashcardsGain, initialState, quizGain, type AnswerResult, type ProgressState } from '../src/lib/gamification';
import { seededRandom } from '../src/lib/random';
import { dueCards, dueCardsText } from '../src/lib/selectors';
import { capDue, cardKey, chapterDueCount, deckOrder, scheduleCard, slug, type CardEntry } from '../src/lib/srs';
import { chapterStars, subjectProgress } from '../src/lib/stats';

const TODAY = '2026-10-04';

function bacS(): ProgressState {
  const s = initialState(TODAY);
  s.profile = createProfile({ name: 'Awa', track: 'bac-s', examDate: '2027-07-01', dailyGoal: 50 });
  return s;
}

test('slug : minuscules sans accents, empreinte quand des symboles sont perdus', () => {
  assert.equal(slug('Définition d’une suite géométrique ?'), 'definition-d-une-suite-geometrique');
  assert.equal(slug('  Œuvre  À   CONNAÎTRE  '), 'oeuvre-a-connaitre');
  // Formules : sans empreinte, (a − b)² et (a + b)² auraient la même clé.
  assert.notEqual(slug('(a − b)² = ?'), slug('(a + b)² = ?'));
  assert.match(slug('(a − b)² = ?'), /^a-b-[0-9a-f]{6}$/);
  assert.match(slug('∫'), /^[0-9a-f]{6}$/);
  assert.ok(slug('x'.repeat(100)).length <= 60);
});

test('cardKey : id de la carte, sinon chapitre + slug du recto', () => {
  assert.equal(cardKey('maths-s-suites', { front: 'Suite arithmétique', back: '…' }), 'maths-s-suites#suite-arithmetique');
  assert.equal(cardKey('maths-s-suites', { id: 'maths-s-suites-c1', front: 'Suite arithmétique', back: '…' }), 'maths-s-suites-c1');
});

test('le contenu n’a pas deux cartes avec la même clé dans un chapitre', () => {
  for (const track of ['bfm', 'bac-s', 'bac-l'] as const) {
    for (const subject of getSubjects(track)) {
      for (const c of subject.chapters) {
        const keys = c.flashcards.map((card) => cardKey(c.id, card));
        assert.equal(new Set(keys).size, keys.length, c.id);
      }
    }
  }
});

test('Leitner 5 boîtes : intervalles 1, 3, 7, 14, 30 jours ; « À revoir » renvoie en boîte 1', () => {
  let e: CardEntry | undefined;
  const days = ['2026-10-04', '2026-10-07', '2026-10-14', '2026-10-28', '2026-11-27'];
  const expected = [
    { box: 2, due: '2026-10-07' },
    { box: 3, due: '2026-10-14' },
    { box: 4, due: '2026-10-28' },
    { box: 5, due: '2026-11-27' },
    { box: 5, due: '2026-12-27' },
  ];
  days.forEach((day, i) => {
    e = scheduleCard(e, true, day);
    assert.deepEqual({ box: e.box, due: e.due }, expected[i], day);
  });
  e = scheduleCard(e, false, '2026-12-27');
  assert.deepEqual(e, { box: 1, due: '2026-12-28', last: '2026-12-27' });
});

test('Leitner : au plus une avance par jour (carte ratée puis retrouvée en fin de paquet)', () => {
  const missed = scheduleCard(undefined, false, TODAY);
  assert.equal(scheduleCard(missed, true, TODAY), missed);
  const known = scheduleCard(undefined, true, TODAY);
  assert.equal(scheduleCard(known, true, TODAY), known);
  assert.equal(scheduleCard(missed, true, '2026-10-05').box, 2);
});

test('Leitner : échéance plafonnée à examDate - 2, jamais avant demain', () => {
  assert.equal(capDue('2026-12-01', TODAY, '2026-10-20'), '2026-10-18');
  assert.equal(capDue('2026-10-10', TODAY, '2026-10-05'), '2026-10-05');
  assert.equal(capDue('2026-12-01', TODAY, TODAY), '2026-12-01');
  assert.equal(scheduleCard({ box: 4, due: TODAY, last: '2026-09-20' }, true, TODAY, '2026-10-20').due, '2026-10-18');
});

test('paquet : cartes à revoir d’abord, puis les nouvelles, puis les autres', () => {
  const cards = ['A', 'B', 'C', 'D', 'E'].map((front) => ({ front, back: front }));
  const entries: Record<string, CardEntry> = {
    'ch#c': { box: 2, due: '2026-10-01', last: '2026-09-28' },
    'ch#e': { box: 1, due: TODAY, last: '2026-10-03' },
    'ch#a': { box: 3, due: '2026-10-09', last: '2026-10-02' },
  };
  const order = deckOrder('ch', cards, entries, TODAY, seededRandom('x')).map((c) => c.front);
  assert.deepEqual(order.slice(0, 2), ['C', 'E']);
  assert.deepEqual(new Set(order.slice(2, 4)), new Set(['B', 'D']));
  assert.equal(order[4], 'A');
  assert.equal(chapterDueCount('ch', cards, entries, TODAY), 2);
});

test('cardReview : suivi sans XP ; dueCards compte les cartes à revoir du contenu actuel', () => {
  const s = bacS();
  const chapter = getSubjects('bac-s')[0].chapters[0];
  const key = cardKey(chapter.id, chapter.flashcards[0]);
  const after = cardReview(s, key, false, TODAY);
  assert.equal(after.xp, 0);
  assert.deepEqual(after.cards[key], { box: 1, due: '2026-10-05', last: TODAY });
  assert.equal(dueCards(after, TODAY).total, 0);
  after.cards['carte-disparue'] = { box: 1, due: TODAY, last: '2026-10-01' };
  const due = dueCards(after, '2026-10-05');
  assert.equal(due.total, 1);
  assert.deepEqual(due.byChapter, [{ chapterId: chapter.id, subjectId: getSubjects('bac-s')[0].id, count: 1 }]);
  assert.equal(dueCardsText(1), '🧠 1 carte à revoir');
  assert.equal(dueCardsText(4), '🧠 4 cartes à revoir');
  assert.equal(dueCardsText(0), null);
});

const ok = (id: string): AnswerResult => ({ questionId: id, subjectId: 'maths-s', correct: true });
const ko = (id: string): AnswerResult => ({ questionId: id, subjectId: 'maths-s', correct: false });

test('étoiles : fiche lue, quiz ≥ 80 %, puis 80 % sur 2 jours et paquet terminé', () => {
  const ch = 'maths-s-suites';
  let s = bacS();
  assert.equal(chapterStars(ch, s, TODAY).stars, 0);
  assert.equal(chapterStars(ch, s, TODAY).next, 'Prochaine étoile : lis la fiche');
  s.fichesRead[ch] = TODAY;
  assert.equal(chapterStars(ch, s, TODAY).stars, 1);
  // 60 % : pas encore la 2e étoile.
  s = quizGain(s, 'chapter', [ok('a'), ok('b'), ok('c'), ko('d'), ko('e')], { chapterId: ch }, TODAY).state;
  assert.equal(chapterStars(ch, s, TODAY).stars, 1);
  assert.deepEqual(s.quiz80Days, {});
  s = quizGain(s, 'chapter', [ok('a'), ok('b'), ok('c'), ok('d'), ko('e')], { chapterId: ch }, TODAY).state;
  assert.deepEqual(s.quiz80Days[ch], [TODAY]);
  let stars = chapterStars(ch, s, TODAY);
  assert.equal(stars.stars, 2);
  assert.equal(stars.next, 'Prochaine étoile : termine le paquet de flashcards et refais le quiz à 80 % un autre jour');
  s = flashcardsGain(s, 'maths-s', ch, TODAY).state;
  assert.equal(s.flashcardsDone[ch], TODAY);
  assert.equal(chapterStars(ch, s, TODAY).next, 'Prochaine étoile : refais le quiz à 80 % un autre jour');
  s = quizGain(s, 'chapter', [ok('a'), ok('b'), ok('c'), ok('d'), ok('e')], { chapterId: ch }, '2026-10-06').state;
  stars = chapterStars(ch, s, '2026-10-06');
  assert.equal(stars.stars, 3);
  assert.equal(stars.next, null);
  assert.equal(stars.faded, false);
  // Badge « Chapitre en or ».
  assert.ok(s.badges['chapter-gold']);
  // Plus de 21 jours sans pratique : étoiles estompées, jamais retirées.
  const later = chapterStars(ch, s, '2026-10-28');
  assert.equal(later.stars, 3);
  assert.equal(later.faded, true);
  assert.equal(later.reminder, '🔄 petit rappel conseillé');
  assert.equal(chapterStars(ch, s, '2026-10-27').faded, false);
});

test('maîtrise d’une matière = étoiles / (3 × chapitres)', () => {
  const s = bacS();
  const subject = getSubjects('bac-s')[0];
  const [c1, c2] = subject.chapters;
  s.fichesRead[c1.id] = TODAY;
  s.fichesRead[c2.id] = TODAY;
  s.quizBest[c2.id] = 90;
  const p = subjectProgress(subject, s);
  assert.equal(p.stars, 3);
  assert.equal(p.maxStars, 3 * subject.chapters.length);
  assert.equal(p.mastery, 3 / (3 * subject.chapters.length));
  assert.equal(p.read, 2);
});

test('badges : avancement « 7/10 fiches » et « Matière maîtrisée »', () => {
  const s = bacS();
  const subject = getSubjects('bac-s')[0];
  for (const c of subject.chapters.slice(0, 7)) s.fichesRead[c.id] = TODAY;
  const fiches10 = BADGES.find((b) => b.id === 'fiches-10')!;
  assert.equal(fiches10.progress!(s).label, `${Math.min(7, subject.chapters.length)}/10 fiches`);
  const master = BADGES.find((b) => b.id === 'subject-master')!;
  assert.equal(master.earned(s), false);
  for (const c of subject.chapters) {
    s.fichesRead[c.id] = TODAY;
    s.quizBest[c.id] = 100;
    s.quiz80Days[c.id] = ['2026-10-01', TODAY];
    s.flashcardsDone[c.id] = TODAY;
  }
  assert.equal(master.earned(s), true);
  assert.equal(master.progress!(s).label, `${subject.chapters.length}/${subject.chapters.length} chapitres en or`);
  assert.equal(BADGES.find((b) => b.id === 'xp-1000')!.progress!({ ...s, xp: 1240 }).label, '1 000/1 000 XP');
});
