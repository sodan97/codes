/** Empreinte FNV-1a 32 bits d'un texte, sur 8 caractères hexadécimaux (somme de contrôle, clés courtes). */
export function fnv1a(text: string): string {
  let h = 0x811c9dc5;
  for (let i = 0; i < text.length; i++) {
    h ^= text.charCodeAt(i);
    h = Math.imul(h, 0x01000193);
  }
  return (h >>> 0).toString(16).padStart(8, '0');
}
