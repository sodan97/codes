import { mkdirSync, writeFileSync } from "node:fs";
// Usage : node scripts/gen-emblem.mjs <dossier-de-sortie>
const out = process.argv[2] ?? ".";
mkdirSync(out, { recursive: true });
const r2 = (n) => Math.round(n * 100) / 100;
const polar = (cx, cy, r, deg) => {
  const a = ((deg - 90) * Math.PI) / 180;
  return [r2(cx + r * Math.cos(a)), r2(cy + r * Math.sin(a))];
};
const arc = (cx, cy, r, start, end) => {
  const [x1, y1] = polar(cx, cy, r, start);
  const [x2, y2] = polar(cx, cy, r, end);
  const large = (((end - start) % 360) + 360) % 360 > 180 ? 1 : 0;
  return `M${x1} ${y1}A${r} ${r} 0 ${large} 1 ${x2} ${y2}`;
};
const arcs = [
  arc(60, 60, 54, 8, 172),
  arc(60, 60, 54, 188, 352),
  arc(60, 60, 45, -64, -8),
  arc(60, 60, 45, 4, 64),
  arc(60, 60, 36, -56, 56),
  arc(60, 60, 27, -48, 10),
  arc(60, 60, 27, 22, 48),
  arc(60, 60, 18, -40, 40),
  arc(60, 93, 9, -62, 62),
  arc(60, 93, 16, -50, 50),
].join("");
const w = "M27 35 43.5 90 60 53 76.5 90 93 35";
writeFileSync(out + "/emblem.ts", `/** Géométrie de l'emblème WOOLYS (viewBox 0 0 120 120) — générée, ne pas éditer à la main. */
export const EMBLEM_ARCS =
  "${arcs}";

export const EMBLEM_W = "${w}";
`);
const svg = (bg) => `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 120 120">
  <defs>
    <linearGradient id="a" x1="8" y1="0" x2="112" y2="0" gradientUnits="userSpaceOnUse"><stop offset="0" stop-color="#E4CBA3"/><stop offset=".42" stop-color="#D4AF7C"/><stop offset=".62" stop-color="#3FB79F"/><stop offset="1" stop-color="#158F7A"/></linearGradient>
    <linearGradient id="g" x1="0" y1="30" x2="0" y2="95" gradientUnits="userSpaceOnUse"><stop offset="0" stop-color="#EFE0C6"/><stop offset=".55" stop-color="#D4AF7C"/><stop offset="1" stop-color="#B08452"/></linearGradient>
    <mask id="m" maskUnits="userSpaceOnUse" x="0" y="0" width="120" height="120"><rect width="120" height="120" fill="#fff"/><path d="${w}" fill="none" stroke="#000" stroke-width="15" stroke-linejoin="round" stroke-linecap="round"/></mask>
    <mask id="n" maskUnits="userSpaceOnUse" x="0" y="0" width="120" height="120"><rect width="120" height="120" fill="#fff"/><path d="${w}" fill="none" stroke="#000" stroke-width="2.2" stroke-linejoin="round" stroke-linecap="round"/></mask>
  </defs>${bg}
  <path d="${arcs}" fill="none" stroke="url(#a)" stroke-width="3.2" stroke-linecap="round" mask="url(#m)"/>
  <path d="${w}" fill="none" stroke="url(#g)" stroke-width="9" stroke-linejoin="round" stroke-linecap="round" mask="url(#n)"/>
</svg>
`;
writeFileSync(out + "/woolys-emblem.svg", svg(""));
writeFileSync(out + "/icon.svg", svg(`\n  <rect width="120" height="120" rx="26" fill="#0B3D33"/>`).replace('viewBox="0 0 120 120"', 'viewBox="-6 -6 132 132"').replace('<rect width="120" height="120" rx="26"', '<rect x="-6" y="-6" width="132" height="132" rx="28"'));
