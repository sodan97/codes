import type { ReactNode } from "react";
import { Breadcrumbs, type Crumb } from "./Breadcrumbs";
import { Reveal } from "./Reveal";
import { Rings } from "./Rings";
import { cn } from "@/lib/cn";

type PageHeroProps = {
  eyebrow?: string;
  title: ReactNode;
  subtitle?: ReactNode;
  breadcrumbs?: Crumb[];
  actions?: ReactNode;
  aside?: ReactNode;
  variant?: "default" | "labs";
};

/** En-tête sombre commun aux pages intérieures. */
export function PageHero({ eyebrow, title, subtitle, breadcrumbs, actions, aside, variant = "default" }: PageHeroProps) {
  const labs = variant === "labs";
  return (
    <section
      className={cn(
        "relative isolate overflow-hidden pb-20 pt-32 text-white sm:pb-24 sm:pt-40",
        labs ? "bg-[#03130F]" : "bg-forest-900",
      )}
    >
      <div className={cn("absolute inset-0 -z-10", labs ? "bg-grid opacity-80" : "bg-grid opacity-60")} />
      <div className="absolute inset-0 -z-10 bg-[radial-gradient(ellipse_at_top_right,rgba(21,143,122,0.28),transparent_55%)]" />
      <div className="absolute inset-0 -z-10 bg-[radial-gradient(ellipse_at_bottom_left,rgba(212,175,124,0.12),transparent_50%)]" />
      <Rings className="absolute -right-40 -top-24 -z-10 h-[38rem] w-[38rem] opacity-40 sm:-right-24" />
      <div className="container-site">
        {breadcrumbs ? <Breadcrumbs items={breadcrumbs} /> : null}
        <div className={cn("grid items-end gap-12", aside ? "lg:grid-cols-[1.4fr_1fr]" : null)}>
          <Reveal className="max-w-3xl">
            {eyebrow ? <p className="eyebrow mb-6 text-gold-400">{eyebrow}</p> : null}
            <h1 className="text-4xl font-bold leading-[1.1] sm:text-5xl lg:text-[3.5rem]">{title}</h1>
            {subtitle ? <div className="mt-6 max-w-2xl text-lg leading-relaxed text-forest-100/80">{subtitle}</div> : null}
            {actions ? <div className="mt-10 flex flex-col gap-4 sm:flex-row sm:flex-wrap">{actions}</div> : null}
          </Reveal>
          {aside ? <Reveal delay={0.15}>{aside}</Reveal> : null}
        </div>
      </div>
      <div className="absolute inset-x-0 bottom-0 h-px bg-gradient-to-r from-transparent via-gold-400/50 to-transparent" />
    </section>
  );
}
