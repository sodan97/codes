import { useId } from "react";
import { EMBLEM_ARCS, EMBLEM_W } from "@/lib/emblem";
import { cn } from "@/lib/cn";

type Tone = "onDark" | "onLight" | "mono";

const gradients: Record<Exclude<Tone, "mono">, { arcs: string[]; w: string[] }> = {
  onDark: {
    arcs: ["#E4CBA3", "#D4AF7C", "#3FB79F", "#158F7A"],
    w: ["#EFE0C6", "#D4AF7C", "#B08452"],
  },
  onLight: {
    arcs: ["#C09764", "#A88450", "#158F7A", "#0F7464"],
    w: ["#D4AF7C", "#B08452", "#8A6838"],
  },
};

type EmblemProps = { className?: string; tone?: Tone; title?: string };

/**
 * Emblème WOOLYS : le « W » entouré d'empreintes concentriques (identité, confiance).
 * Recréation vectorielle de la charte graphique — à remplacer par le fichier vectoriel
 * officiel si le designer le fournit (public/brand/).
 */
export function Emblem({ className, tone = "onDark", title }: EmblemProps) {
  const uid = useId().replace(/[^a-zA-Z0-9_-]/g, "");
  const ids = { arcs: `a${uid}`, w: `g${uid}`, cut: `m${uid}`, line: `n${uid}` };
  const mono = tone === "mono";
  const palette = mono ? null : gradients[tone];

  return (
    <svg
      viewBox="0 0 120 120"
      className={className}
      role={title ? "img" : undefined}
      aria-hidden={title ? undefined : true}
      focusable="false"
    >
      {title ? <title>{title}</title> : null}
      <defs>
        {palette ? (
          <>
            <linearGradient id={ids.arcs} x1="8" y1="0" x2="112" y2="0" gradientUnits="userSpaceOnUse">
              <stop offset="0" stopColor={palette.arcs[0]} />
              <stop offset="0.42" stopColor={palette.arcs[1]} />
              <stop offset="0.62" stopColor={palette.arcs[2]} />
              <stop offset="1" stopColor={palette.arcs[3]} />
            </linearGradient>
            <linearGradient id={ids.w} x1="0" y1="30" x2="0" y2="95" gradientUnits="userSpaceOnUse">
              <stop offset="0" stopColor={palette.w[0]} />
              <stop offset="0.55" stopColor={palette.w[1]} />
              <stop offset="1" stopColor={palette.w[2]} />
            </linearGradient>
          </>
        ) : null}
        <mask id={ids.cut} maskUnits="userSpaceOnUse" x="0" y="0" width="120" height="120">
          <rect width="120" height="120" fill="#fff" />
          <path d={EMBLEM_W} fill="none" stroke="#000" strokeWidth="15" strokeLinejoin="round" strokeLinecap="round" />
        </mask>
        <mask id={ids.line} maskUnits="userSpaceOnUse" x="0" y="0" width="120" height="120">
          <rect width="120" height="120" fill="#fff" />
          <path d={EMBLEM_W} fill="none" stroke="#000" strokeWidth="2.2" strokeLinejoin="round" strokeLinecap="round" />
        </mask>
      </defs>
      <path
        d={EMBLEM_ARCS}
        fill="none"
        stroke={mono ? "currentColor" : `url(#${ids.arcs})`}
        strokeWidth="3.2"
        strokeLinecap="round"
        mask={`url(#${ids.cut})`}
      />
      <path
        d={EMBLEM_W}
        fill="none"
        stroke={mono ? "currentColor" : `url(#${ids.w})`}
        strokeWidth="9"
        strokeLinejoin="round"
        strokeLinecap="round"
        mask={`url(#${ids.line})`}
      />
    </svg>
  );
}

type LogoProps = {
  tone?: "onDark" | "onLight";
  showTagline?: boolean;
  className?: string;
};

/** Logo horizontal : emblème + logotype WOOLYS + signature. */
export function Logo({ tone = "onDark", showTagline = true, className }: LogoProps) {
  const dark = tone === "onDark";
  return (
    <span className={cn("inline-flex items-center gap-3", className)}>
      <Emblem tone={tone} className="h-10 w-10 shrink-0 sm:h-11 sm:w-11" />
      <span className="flex flex-col leading-none">
        <span
          className={cn(
            "font-display text-[1.35rem] font-semibold tracking-[0.14em] sm:text-[1.5rem]",
            dark ? "text-white" : "text-forest-800",
          )}
        >
          WOOLYS
        </span>
        {showTagline ? (
          <span
            className={cn(
              "mt-1 font-display text-[0.52rem] font-semibold uppercase tracking-[0.3em] sm:text-[0.58rem]",
              dark ? "text-gold-400" : "text-gold-700",
            )}
          >
            Cybersecurity &amp; Trust
          </span>
        ) : null}
      </span>
    </span>
  );
}
