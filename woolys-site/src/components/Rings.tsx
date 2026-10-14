import { cn } from "@/lib/cn";

/** Anneaux concentriques décoratifs, écho graphique à l'emblème WOOLYS. */
export function Rings({ className, count = 7 }: { className?: string; count?: number }) {
  return (
    <svg viewBox="0 0 600 600" className={cn("pointer-events-none", className)} aria-hidden="true" fill="none">
      <defs>
        <linearGradient id="rings-gradient" x1="0" y1="0" x2="600" y2="0" gradientUnits="userSpaceOnUse">
          <stop offset="0" stopColor="#D4AF7C" />
          <stop offset="1" stopColor="#158F7A" />
        </linearGradient>
      </defs>
      {Array.from({ length: count }, (_, i) => (
        <circle
          key={i}
          cx="300"
          cy="300"
          r={60 + i * 36}
          stroke="url(#rings-gradient)"
          strokeWidth="1"
          strokeOpacity={0.5 - i * 0.05}
          strokeDasharray={i % 3 === 1 ? "2 10" : i % 3 === 2 ? "120 24" : undefined}
        />
      ))}
    </svg>
  );
}
