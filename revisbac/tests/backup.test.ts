/// <reference types="node" />
import assert from 'node:assert/strict';
import { test } from 'node:test';

import { getSubjects } from '../src/data/catalog';
import { addDays } from '../src/lib/dates';
import {
  BACKUP_ERRORS,
  backupShareMessage,
  decodeBase64,
  encodeBase64,
  exportCode,
  importCode,
  previewText,
} from '../src/lib/backup';
import { createProfile, initialState, STATE_VERSION, withProfile, type ProgressState } from '../src/lib/gamification';
import { fnv1a } from '../src/lib/hash';
import { CARD_BOX_DAYS, cardKey } from '../src/lib/srs';

const TODAY = '2026-10-04';

function awa(): ProgressState {
  const s = withProfile(initialState(TODAY), createProfile({ name: 'Awa Ndiaye', track: 'bac-s', examDate: '2027-07-01', dailyGoal: 100 }), TODAY);
  s.xp = 1240;
  s.streak = { current: 12, best: 15, lastDay: TODAY, freezes: 1 };
  s.fichesRead = { 'maths-s-suites': '2026-09-30' };
  s.quizBest = { 'maths-s-suites': 90 };
  s.mistakes = { 'maths-s-suites-q1': { box: 2, due: '2026-10-06', fails: 1, last: '2026-10-03' } };
  s.cards = { 'maths-s-suites#suite-arithmetique': { box: 3, due: '2026-10-11', last: TODAY } };
  s.badges = Object.fromEntries(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i'].map((id) => [id, TODAY]));
  s.history = { [TODAY]: 60 };
  return s;
}

test('base64 UTF-8 maison : accents, emoji et caractères hors BMP, identique à Buffer', () => {
  for (const text of ['', 'a', 'ab', 'abc', 'Élève · Révision ✅ 🦁 — « (a − b)² »', 'x'.repeat(1000)]) {
    const b64 = encodeBase64(text);
    assert.equal(b64, Buffer.from(text, 'utf8').toString('base64'));
    assert.equal(decodeBase64(b64), text);
  }
  assert.throws(() => decodeBase64('abc'));
  assert.throws(() => decodeBase64('ab$='));
  // Octets qui ne forment pas de l'UTF-8 valide.
  assert.throws(() => decodeBase64(Buffer.from([0xc3]).toString('base64')));
});

test('export puis import : même progression, sans l’historique des XP', () => {
  const s = awa();
  const code = exportCode(s);
  assert.match(code, /^RB2:[A-Za-z0-9+/=]+:[0-9a-f]{8}$/);
  const r = importCode(code, TODAY);
  assert.ok(r.ok);
  const { history, ...rest } = s;
  assert.deepEqual(history, { [TODAY]: 60 });
  assert.deepEqual(r.state, { ...rest, history: {} });
  assert.deepEqual(r.preview, { name: 'Awa Ndiaye', trackLabel: 'Bac S', xp: 1240, streak: 12, badges: 9 });
  assert.equal(previewText(r.preview), 'Awa Ndiaye · Bac S · 1 240 XP · série 12 · 9 badges');
});

test('import : le message partagé entier, coupé par des retours à la ligne, est accepté', () => {
  const code = exportCode(awa());
  const message = backupShareMessage(code);
  assert.ok(message.startsWith('Garde ce message : il permet de récupérer ta progression RéviBac.\n'));
  const wrapped = `${message.slice(0, 80)}\n${message.slice(80, 150)} \n${message.slice(150)}\n\nEnvoyé depuis mon téléphone`;
  const r = importCode(wrapped, TODAY);
  assert.ok(r.ok);
});

test('import : code incomplet, abîmé ou d’une version inconnue', () => {
  const code = exportCode(awa());
  const [, payload, sum] = code.split(':');
  assert.deepEqual(importCode('', TODAY), { ok: false, error: BACKUP_ERRORS.incomplete });
  assert.deepEqual(importCode('bonjour', TODAY), { ok: false, error: BACKUP_ERRORS.incomplete });
  assert.deepEqual(importCode(`RB2:${payload}`, TODAY), { ok: false, error: BACKUP_ERRORS.incomplete });
  assert.deepEqual(importCode(code.slice(0, -3), TODAY), { ok: false, error: BACKUP_ERRORS.incomplete });
  assert.deepEqual(importCode(`RB2:${payload.slice(0, 40)}:${sum}`, TODAY), { ok: false, error: BACKUP_ERRORS.damaged });
  // Un seul caractère changé.
  const i = 30;
  const altered = payload.slice(0, i) + (payload[i] === 'A' ? 'B' : 'A') + payload.slice(i + 1);
  assert.deepEqual(importCode(`RB2:${altered}:${sum}`, TODAY), { ok: false, error: BACKUP_ERRORS.damaged });
  assert.deepEqual(importCode(`RB3:${payload}:${sum}`, TODAY), { ok: false, error: BACKUP_ERRORS.version });
  assert.deepEqual(importCode(`RB0:${payload}:${sum}`, TODAY), { ok: false, error: BACKUP_ERRORS.version });
  // Somme de contrôle juste mais contenu illisible.
  const junk = encodeBase64('pas du json');
  assert.deepEqual(importCode(`RB1:${junk}:${fnv1a(junk)}`, TODAY), { ok: false, error: BACKUP_ERRORS.damaged });
  // Examen inconnu de cette version.
  const future = encodeBase64(JSON.stringify({ version: 9, profile: { name: 'Awa', track: 'cfee' } }));
  assert.deepEqual(importCode(`RB1:${future}:${fnv1a(future)}`, TODAY), { ok: false, error: BACKUP_ERRORS.version });
});

test('import d’un ancien format : la sauvegarde passe par migrate()', () => {
  const v1 = {
    version: 1,
    profile: { name: 'Moussa', track: 'bfm', examDate: '2027-07-10', dailyGoal: 100 },
    xp: 420,
    streak: { current: 4, best: 9, lastDay: '2026-10-03', freezes: 1 },
    mistakes: { 'maths-bfm-thales-q1': 2 },
    badges: { 'first-fiche': '2026-09-20' },
  };
  const payload = encodeBase64(JSON.stringify(v1));
  const r = importCode(`RB1:${payload}:${fnv1a(payload)}`, TODAY);
  assert.ok(r.ok);
  assert.equal(r.state.version, STATE_VERSION);
  assert.deepEqual(r.state.mistakes, { 'maths-bfm-thales-q1': { box: 1, due: TODAY, fails: 2 } });
  assert.deepEqual(r.state.cards, {});
  assert.equal(r.state.today.quests.length, 3);
  assert.equal(previewText(r.preview), 'Moussa · BFEM · 420 XP · série 4 · 1 badge');
});

test('import d’un code RB1 (état JSON tel quel) : toujours accepté', () => {
  const s = awa();
  const { history: _h, ...rest } = s;
  const payload = encodeBase64(JSON.stringify(rest));
  const r = importCode(`RB1:${payload}:${fnv1a(payload)}`, TODAY);
  assert.ok(r.ok);
  assert.deepEqual(r.state, { ...rest, history: {} });
});

/**
 * Élève de Bac S en fin d'année : toutes les flashcards revues (1 carte sur 10 à échéance plafonnée avant l'examen).
 * `full` : en plus, 1 question sur 4 en erreur suivie et tous les chapitres lus, réussis et révisés.
 */
function endOfYear(full: boolean): ProgressState {
  const s = awa();
  s.cards = {};
  s.mistakes = {};
  let i = 0;
  for (const sub of getSubjects('bac-s')) {
    for (const ch of sub.chapters) {
      for (const card of ch.flashcards) {
        i++;
        const box = ((i % 5) + 1) as 1 | 2 | 3 | 4 | 5;
        const last = addDays(TODAY, -(i % 9));
        const due = i % 10 === 0 ? addDays(last, 2) : addDays(last, CARD_BOX_DAYS[box]);
        s.cards[cardKey(ch.id, card)] = { box, due, last };
      }
      if (!full) continue;
      ch.quiz.forEach((q, k) => {
        if (k % 4 === 0) s.mistakes[q.id] = { box: ((k % 3) + 1) as 1 | 2 | 3, due: addDays(TODAY, k % 7), fails: 1 + (k % 3), ...(k % 2 ? { last: addDays(TODAY, -k) } : {}) };
      });
      s.fichesRead[ch.id] = TODAY;
      s.quizBest[ch.id] = 90;
      s.flashcardsDone[ch.id] = TODAY;
      s.quiz80Days[ch.id] = ['2026-09-30', TODAY];
    }
  }
  assert.equal(Object.keys(s.cards).length, 440);
  return s;
}

test('code compact : flashcards et erreurs reviennent à l’identique, et le code tient dans un message', () => {
  for (const full of [false, true]) {
    const s = endOfYear(full);
    const r = importCode(exportCode(s), TODAY);
    assert.ok(r.ok);
    assert.deepEqual(r.state.cards, s.cards);
    assert.deepEqual(r.state.mistakes, s.mistakes);
  }
  // Toutes les flashcards du Bac S (environ 72 000 caractères au format RB1).
  const cardsOnly = exportCode(endOfYear(false)).length;
  assert.ok(cardsOnly < 30_000, `flashcards seules : ${cardsOnly} caractères`);
  // Tout rempli : sous la limite d'un message WhatsApp (65 536 caractères).
  const all = exportCode(endOfYear(true)).length;
  assert.ok(all < 60_000, `progression complète : ${all} caractères`);
});
