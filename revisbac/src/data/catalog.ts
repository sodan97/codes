import { subjects } from './content';
import type { Chapter, Question, Subject, TrackId } from './types';

export { getTrack, isTrackId, tracks, type TrackInfo } from './tracks';

/** Matières facultatives (LV2) : masquables par l'élève et exclues du défi du jour commun. */
export const OPTIONAL_SUBJECTS: ReadonlySet<string> = new Set(['espagnol-bfm', 'espagnol-bac']);

export function isOptional(subjectId: string): boolean {
  return OPTIONAL_SUBJECTS.has(subjectId);
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
