import { Reveal } from "@/components/Reveal";
import { SectionHeader } from "@/components/SectionHeader";
import type { TitledText } from "@/data/types";
import { cn } from "@/lib/cn";

/** Grille des offres d'un pôle. */
export function PoleOffers({
  eyebrow,
  title,
  intro,
  offers,
  tone = "light",
}: {
  eyebrow: string;
  title: string;
  intro?: string;
  offers: TitledText[];
  tone?: "light" | "dark";
}) {
  const dark = tone === "dark";
  return (
    <section className={cn("py-24 sm:py-28", dark ? "bg-[#03130F] text-white" : "bg-sand-50")} aria-labelledby="offres">
      <div className="container-site">
        <SectionHeader id="offres" tone={tone} eyebrow={eyebrow} title={title} intro={intro} />
        <ul className={cn("mt-14 grid gap-px overflow-hidden rounded-3xl sm:grid-cols-2 lg:grid-cols-3", dark ? "bg-emerald-400/15" : "bg-forest-800/10")}>
          {offers.map((offer, index) => (
            <Reveal
              as="li"
              key={offer.title}
              delay={(index % 3) * 0.06}
              className={cn("group p-8 transition-colors duration-500", dark ? "bg-[#03130F] hover:bg-forest-950" : "bg-white hover:bg-sand-50")}
            >
              <span className={cn("font-mono text-xs", dark ? "text-emerald-400" : "text-gold-700")}>{String(index + 1).padStart(2, "0")}</span>
              <h3 className={cn("mt-4 text-lg font-semibold", dark ? "text-white" : "text-forest-800")}>{offer.title}</h3>
              <p className={cn("mt-2 text-[0.95rem] leading-relaxed", dark ? "text-forest-100/70" : "text-muted")}>{offer.text}</p>
              <span className={cn("mt-6 block h-px w-8 transition-all duration-500 group-hover:w-16", dark ? "bg-emerald-400" : "bg-gold-400")} />
            </Reveal>
          ))}
        </ul>
      </div>
    </section>
  );
}
