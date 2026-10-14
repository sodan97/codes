import { Emblem } from "@/components/Logo";

const C = 260;

const ring = (count: number, radius: number, offset = 0) =>
  Array.from({ length: count }, (_, i) => {
    const angle = ((360 / count) * i + offset) * (Math.PI / 180);
    return { x: Math.round((C + radius * Math.cos(angle)) * 10) / 10, y: Math.round((C + radius * Math.sin(angle)) * 10) / 10 };
  });

const outer = ring(8, 215, 10);
const middle = ring(6, 155, 30);
const inner = ring(5, 95, 0);

const labels = [
  { text: "Gouvernance", className: "left-[2%] top-[18%]", delay: "0s" },
  { text: "Données", className: "right-[0%] top-[12%]", delay: "1.2s" },
  { text: "Cybersécurité", className: "right-[-2%] top-[56%]", delay: "2.4s" },
  { text: "Résilience", className: "left-[6%] bottom-[12%]", delay: "0.6s" },
  { text: "Intelligence", className: "left-[-4%] top-[52%]", delay: "1.8s" },
  { text: "Réseau", className: "right-[14%] bottom-[4%]", delay: "3s" },
];

/**
 * Visuel abstrait du hero : réseau, données, gouvernance, cybersécurité, intelligence, résilience.
 * 100 % SVG + CSS (aucun JavaScript), animations très lentes, désactivées si prefers-reduced-motion.
 */
export function NetworkVisual() {
  const orbitStyle = { transformOrigin: `${C}px ${C}px`, transformBox: "view-box" as const };
  return (
    <div className="relative mx-auto aspect-square w-full max-w-[34rem]" aria-hidden="true">
      <div className="absolute inset-[18%] rounded-full bg-emerald-500/25 blur-3xl" />
      <svg viewBox="0 0 520 520" className="relative h-full w-full" fill="none">
        <defs>
          <radialGradient id="nv-core" cx="50%" cy="50%" r="50%">
            <stop offset="0" stopColor="#158F7A" stopOpacity="0.35" />
            <stop offset="1" stopColor="#158F7A" stopOpacity="0" />
          </radialGradient>
          <linearGradient id="nv-line" x1="0" y1="0" x2="520" y2="520" gradientUnits="userSpaceOnUse">
            <stop offset="0" stopColor="#D4AF7C" />
            <stop offset="1" stopColor="#3FB79F" />
          </linearGradient>
        </defs>

        <circle cx={C} cy={C} r="250" fill="url(#nv-core)" />
        {[250, 215, 155, 95].map((r, i) => (
          <circle
            key={r}
            cx={C}
            cy={C}
            r={r}
            stroke="url(#nv-line)"
            strokeOpacity={i === 0 ? 0.15 : 0.28}
            strokeDasharray={i === 0 ? "2 9" : i === 2 ? "3 7" : undefined}
          />
        ))}

        {/* Orbite extérieure */}
        <g className="animate-spin-slow" style={orbitStyle}>
          {outer.map((p, i) => {
            const next = outer[(i + 1) % outer.length];
            const link = middle[i % middle.length];
            return (
              <g key={i}>
                <line x1={p.x} y1={p.y} x2={next.x} y2={next.y} stroke="#D4AF7C" strokeOpacity="0.22" />
                <line x1={p.x} y1={p.y} x2={link.x} y2={link.y} stroke="#3FB79F" strokeOpacity="0.18" />
                <circle cx={p.x} cy={p.y} r={i % 2 === 0 ? 4.5 : 3} fill={i % 2 === 0 ? "#D4AF7C" : "#3FB79F"} />
                <circle cx={p.x} cy={p.y} r="10" stroke="#D4AF7C" strokeOpacity="0.35" className="animate-glow" style={{ animationDelay: `${i * 0.5}s` }} />
              </g>
            );
          })}
        </g>

        {/* Orbite intermédiaire */}
        <g className="animate-spin-slower" style={orbitStyle}>
          {middle.map((p, i) => {
            const next = middle[(i + 2) % middle.length];
            const link = inner[i % inner.length];
            return (
              <g key={i}>
                <line x1={p.x} y1={p.y} x2={next.x} y2={next.y} stroke="#3FB79F" strokeOpacity="0.2" />
                <line x1={p.x} y1={p.y} x2={link.x} y2={link.y} stroke="#D4AF7C" strokeOpacity="0.2" />
                <circle cx={p.x} cy={p.y} r="3.5" fill="#3FB79F" />
              </g>
            );
          })}
        </g>

        {inner.map((p, i) => (
          <g key={i}>
            <line x1={C} y1={C} x2={p.x} y2={p.y} stroke="#D4AF7C" strokeOpacity="0.25" strokeDasharray="2 4" />
            <circle cx={p.x} cy={p.y} r="3" fill="#E4CBA3" className="animate-glow" style={{ animationDelay: `${i * 0.8}s` }} />
          </g>
        ))}
        <circle cx={C} cy={C} r="66" fill="#041E19" stroke="#D4AF7C" strokeOpacity="0.4" />
      </svg>

      <Emblem tone="onDark" className="absolute left-1/2 top-1/2 h-[21%] w-[21%] -translate-x-1/2 -translate-y-1/2" />

      {labels.map((label) => (
        <span
          key={label.text}
          className={`absolute hidden animate-float rounded-full border border-white/15 bg-forest-950/70 px-3.5 py-1.5 text-xs font-medium text-forest-100 shadow-lg backdrop-blur-md sm:inline-flex ${label.className}`}
          style={{ animationDelay: label.delay }}
        >
          <span className="mr-2 mt-[5px] inline-block h-1.5 w-1.5 rounded-full bg-gold-400" />
          {label.text}
        </span>
      ))}
    </div>
  );
}
