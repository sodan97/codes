// Textes de partage (module pur, testé). Aucun lien : l'application n'est pas encore publiée,
// et jamais de comparaison chiffrée avec d'autres élèves.
import { formatDay } from './dates';
import { formatNote, mention } from './gamification';

export function dailyShareText({
  day,
  trackLabel,
  correct,
  total,
  grid,
  streak,
}: {
  day: string;
  trackLabel: string;
  correct: number;
  total: number;
  grid: string;
  streak: number;
}): string {
  const lines = [`RéviBac · Défi du jour du ${formatDay(day)} (${trackLabel})`, `${correct}/${total} ${grid}`];
  if (streak >= 2) lines.push(`🔥 Série : ${streak} jours`);
  lines.push('Fais le même défi que moi dans RéviBac !');
  return lines.join('\n');
}

export function examShareText({ subjectName, note }: { subjectName: string; note: number }): string {
  return `RéviBac · Examen blanc de ${subjectName} : ${formatNote(note)}/20 (mention ${mention(note)}).\nQuestions de révision RéviBac, pas un sujet officiel.`;
}
