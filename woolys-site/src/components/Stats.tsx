import { Reveal } from "./Reveal";
import { cn } from "@/lib/cn";

export type StatItem = { value: string; label: string };

/**
 * Chiffres clés. N'afficher QUE des informations vérifiables (ex. périmètre de l'offre).
 * Ne jamais y inscrire de nombre de clients, de collaborateurs ou d'années d'existence non confirmés.
 */
export function Stats({ items, tone = "light" }: { items: StatItem[]; tone?: "light" | "dark" }) {
  const dark = tone === "dark";
  return (
    <dl className={cn("grid grid-cols-2 gap-px overflow-hidden rounded-2xl lg:grid-cols-4", dark ? "bg-white/10" : "bg-forest-800/10")}>
      {items.map((item, index) => (
        <Reveal key={item.label} delay={index * 0.08} className={cn("flex flex-col p-6 sm:p-8", dark ? "bg-forest-900" : "bg-white")}>
          <dt className={cn("order-2 mt-2 text-sm leading-snug", dark ? "text-forest-100/70" : "text-muted")}>{item.label}</dt>
          <dd className={cn("order-1 font-display text-3xl font-bold sm:text-4xl", dark ? "text-gold-400" : "text-forest-800")}>
            {item.value}
          </dd>
        </Reveal>
      ))}
    </dl>
  );
}
