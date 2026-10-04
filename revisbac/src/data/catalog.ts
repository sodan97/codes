import { subjects } from './content';
import type { Chapter, Question, Subject, TrackId } from './types';

export interface TrackInfo {
  id: TrackId;
  label: string;
  description: string;
  emoji: string;
  /** Date d'examen indicative (AAAA-MM-JJ), modifiable par l'élève dans son profil. */
  defaultExamDate: string;
}

export const tracks: TrackInfo[] = [
  {
    id: 'bfm',
    label: 'BFM',
    description: 'Brevet de Fin d’études Moyennes — classe de 3e',
    emoji: '🎒',
    defaultExamDate: '2027-07-12',
  },
  {
    id: 'bac-s',
    label: 'Bac S',
    description: 'Terminale scientifique (S1, S2…)',
    emoji: '🔬',
    defaultExamDate: '2027-07-01',
  },
  {
    id: 'bac-l',
    label: 'Bac L',
    description: 'Terminale littéraire (L1, L2, L’…)',
    emoji: '📚',
    defaultExamDate: '2027-07-01',
  },
];

export function getTrack(id: TrackId): TrackInfo {
  return tracks.find((t) => t.id === id) ?? tracks[0];
}

export function getSubjects(track: TrackId): Subject[] {
  return subjects.filter((s) => s.tracks.includes(track));
}

export function getSubject(id: string): Subject | undefined {
  return subjects.find((s) => s.id === id);
}

export interface ChapterRef {
  chapter: Chapter;
  subject: Subject;
}

export interface QuestionRef extends ChapterRef {
  question: Question;
}

const chapterIndex = new Map<string, ChapterRef>();
const questionIndex = new Map<string, QuestionRef>();
for (const subject of subjects) {
  for (const chapter of subject.chapters) {
    chapterIndex.set(chapter.id, { chapter, subject });
    for (const question of chapter.quiz) questionIndex.set(question.id, { question, chapter, subject });
  }
}

export function getChapter(id: string): ChapterRef | undefined {
  return chapterIndex.get(id);
}

export function getQuestion(id: string): QuestionRef | undefined {
  return questionIndex.get(id);
}

export function questionsOfSubject(subject: Subject): QuestionRef[] {
  return subject.chapters.flatMap((chapter) => chapter.quiz.map((question) => ({ question, chapter, subject })));
}
