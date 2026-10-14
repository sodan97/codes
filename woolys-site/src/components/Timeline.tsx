import { Reveal } from "./Reveal";
import { cn } from "@/lib/cn";

type Step = { number?: string; title: string; text: string };

/** Frise d'étapes : horizontale sur grand écran, verticale sur mobile. */
export function Timeline({ steps, tone = "light" }: { steps: Step[]; tone?: "light" | "dark" }) {
  const dark = tone === "dark";
  return (
    <ol className="relative grid gap-10 md:grid-cols-2 lg:grid-cols-4 lg:gap-8">
      <span
        aria-hidden="true"
        className={cn(
          "absolute left-[1.375rem] top-2 bottom-2 w-px lg:left-0 lg:right-0 lg:top-[1.375rem] lg:bottom-auto lg:h-px lg:w-auto md:hidden lg:block",
          dark ? "bg-gradient-to-b from-gold-400/60 via-white/15 to-transparent lg:bg-gradient-to-r" : "bg-gradient-to-b from-gold-500/70 via-forest-800/15 to-transparent lg:bg-gradient-to-r",
        )}
      />
      {steps.map((step, index) => (
        <Reveal as="li" key={step.title} delay={index * 0.1} className="relative pl-16 md:pl-0">
          <span
            className={cn(
              "absolute left-0 top-0 flex h-11 w-11 items-center justify-center rounded-full border font-display text-sm font-bold md:relative md:mb-6",
              dark ? "border-gold-400/60 bg-forest-950 text-gold-400" : "border-gold-500/60 bg-sand-50 text-forest-800",
            )}
          >
            {step.number ?? String(index + 1).padStart(2, "0")}
          </span>
          <h3 className={cn("text-xl font-semibold", dark ? "text-white" : "text-forest-800")}>{step.title}</h3>
          <p className={cn("mt-2 leading-relaxed", dark ? "text-forest-100/75" : "text-muted")}>{step.text}</p>
        </Reveal>
      ))}
    </ol>
  );
}
