// Vérifie l'intégrité du contenu pédagogique (ids uniques, réponses valides, trous cohérents).
// Usage : npm run validate
import { readdirSync } from 'node:fs';
import { join, dirname } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

import { cardKey } from '../src/lib/srs.ts';

const dir = join(dirname(fileURLToPath(import.meta.url)), '..', 'src', 'data', 'content');
const TRACKS = new Set(['bfm', 'bac-s', 'bac-l']);
const errors = [];
const warnings = [];
const ids = new Set();
const err = (where, msg) => errors.push(`${where}: ${msg}`);
const warn = (where, msg) => warnings.push(`${where}: ${msg}`);
const countOf = (list, word) => list.filter((w) => w === word).length;
const unique = (id, where) => (ids.has(id) ? err(where, `id dupliqué "${id}"`) : ids.add(id));

let nChapters = 0, nQuestions = 0, nCards = 0;
const files = readdirSync(dir).filter((f) => f.endsWith('.ts') && f !== 'index.ts');
for (const file of files) {
  const { default: s } = await import(pathToFileURL(join(dir, file)).href);
  const w = file;
  if (!s || !s.id) { err(w, 'export default manquant'); continue; }
  if (`${s.id}.ts` !== file) err(w, `le fichier doit s'appeler ${s.id}.ts`);
  unique(s.id, w);
  if (!s.tracks?.length || s.tracks.some((t) => !TRACKS.has(t))) err(w, 'tracks invalides');
  if (!/^#[0-9A-Fa-f]{6}$/.test(s.color)) err(w, 'couleur invalide');
  if (!s.chapters?.length) err(w, 'aucun chapitre');
  for (const c of s.chapters ?? []) {
    nChapters++;
    const wc = `${w} > ${c.id}`;
    unique(c.id, wc);
    if (!c.id.startsWith(s.id + '-')) err(wc, `l'id du chapitre doit commencer par "${s.id}-"`);
    if (!c.title || !c.summary) err(wc, 'titre/résumé manquant');
    if (!c.essentials?.length) err(wc, 'essentials vide');
    if (!c.sections?.length) err(wc, 'sections vides');
    for (const sec of c.sections ?? []) if (!sec.blocks?.length) err(wc, `section vide "${sec.title}"`);
    if ((c.flashcards?.length ?? 0) < 4) err(wc, 'moins de 4 flashcards');
    nCards += c.flashcards?.length ?? 0;
    // Chaque carte est suivie en répétition espacée par sa clé (id, sinon recto) : deux cartes d'un chapitre
    // ne doivent pas partager la même.
    const keys = new Set();
    for (const card of c.flashcards ?? []) {
      if (!card.front || !card.back) err(wc, 'flashcard sans recto ou verso');
      if (card.id !== undefined) {
        unique(card.id, wc);
        if (!String(card.id).startsWith(c.id + '-')) err(wc, `l'id de carte "${card.id}" doit commencer par "${c.id}-"`);
      }
      const key = cardKey(c.id, card);
      if (keys.has(key)) err(wc, `deux flashcards ont la même clé "${key}" : ajouter un id à l'une d'elles`);
      keys.add(key);
    }
    if ((c.quiz?.length ?? 0) < 6) err(wc, 'moins de 6 questions');
    for (const q of c.quiz ?? []) {
      nQuestions++;
      const wq = `${wc} > ${q.id}`;
      unique(q.id, wq);
      if (!q.id.startsWith(c.id + '-')) err(wq, `l'id doit commencer par "${c.id}-"`);
      if (!q.explanation) err(wq, 'explication manquante');
      if (q.type === 'qcm') {
        if (!(q.choices?.length >= 2)) err(wq, 'choix insuffisants');
        if (!Number.isInteger(q.answer) || q.answer < 0 || q.answer >= q.choices.length) err(wq, 'index de réponse invalide');
        if (new Set(q.choices).size !== q.choices.length) err(wq, 'choix dupliqués');
      } else if (q.type === 'vrai-faux') {
        if (typeof q.answer !== 'boolean') err(wq, 'réponse non booléenne');
      } else if (q.type === 'trous') {
        const blanks = q.prompt.split('___').length - 1;
        if (blanks < 1) err(wq, 'aucun trou "___"');
        if (blanks !== q.answers?.length) err(wq, `${blanks} trous mais ${q.answers?.length} réponses`);
        for (const a of new Set(q.answers ?? [])) {
          // Une réponse répétée doit figurer autant de fois dans la banque, sinon la question est impossible.
          const need = countOf(q.answers, a);
          const have = countOf(q.bank ?? [], a);
          if (have === 0) err(wq, `"${a}" absent de la banque de mots`);
          else if (have < need) err(wq, `"${a}" doit figurer ${need} fois dans la banque (${have} seulement)`);
        }
        // Mot en double sans réponse répétée : pas bloquant, mais signalé.
        for (const w of new Set(q.bank ?? [])) {
          if (countOf(q.bank, w) > Math.max(1, countOf(q.answers ?? [], w))) warn(wq, `mot "${w}" en double dans la banque`);
        }
        if ((q.bank?.length ?? 0) <= (q.answers?.length ?? 0)) err(wq, 'la banque doit contenir des distracteurs');
      } else err(wq, `type inconnu ${q.type}`);
    }
  }
}
if (warnings.length) console.warn(`⚠️ ${warnings.length} avertissement(s) :\n${warnings.join('\n')}`);
console.log(`${files.length} matières, ${nChapters} chapitres, ${nCards} flashcards, ${nQuestions} questions`);
if (errors.length) { console.error(errors.join('\n')); process.exit(1); }
console.log('Contenu valide ✔');
