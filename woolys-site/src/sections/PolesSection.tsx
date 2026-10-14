import Link from "next/link";
import { Icon } from "@/components/Icon";
import { Reveal } from "@/components/Reveal";
import { SectionHeader } from "@/components/SectionHeader";
import { poles } from "@/data/poles";

export function PolesSection() {
  return (
    <section className="relative overflow-hidden bg-forest-900 py-24 text-white sm:py-32" aria-labelledby="poles-title">
      <div className="absolute inset-0 bg-grid opacity-40" />
      <div className="container-site relative">
        <SectionHeader
          id="poles-title"
          tone="dark"
          eyebrow="Une architecture en trois pôles"
          title="Conseiller. Former. Innover."
          intro="Trois pôles complémentaires pour accompagner votre organisation, de la stratégie à la montée en compétences, jusqu'aux technologies émergentes."
        />
        <div className="mt-16 grid gap-6 lg:grid-cols-3">
          {poles.map((pole, index) => (
            <Reveal key={pole.slug} delay={index * 0.1} className="h-full">
              <article className="group relative flex h-full flex-col overflow-hidden rounded-3xl border border-white/10 bg-gradient-to-b from-white/[0.06] to-white/[0.01] p-8 transition-all duration-500 hover:border-gold-400/50 sm:p-10">
                <div className="absolute -right-16 -top-16 h-48 w-48 rounded-full border border-gold-400/10 transition-transform duration-700 group-hover:scale-125" />
                <div className="absolute -right-8 -top-8 h-32 w-32 rounded-full border border-emerald-400/15 transition-transform duration-700 group-hover:scale-125" />
                <div className="flex items-center justify-between">
                  <span className="font-display text-5xl font-bold text-white/10 transition-colors duration-500 group-hover:text-gold-400/30">
                    {pole.number}
                  </span>
                  <span className="flex h-12 w-12 items-center justify-center rounded-2xl bg-white/5 text-gold-400">
                    <Icon name={pole.icon} className="h-6 w-6" />
                  </span>
                </div>
                <h3 className="mt-10 font-display text-2xl font-bold tracking-wide">
                  <Link href={pole.href} className="after:absolute after:inset-0 focus-visible:outline-none">
                    {pole.name}
                  </Link>
                </h3>
                <p className="mt-4 flex-1 text-lg leading-relaxed text-forest-100/80">{pole.tagline}</p>
                <ul className="mt-8 flex flex-wrap gap-2">
                  {pole.offers.slice(0, 4).map((offer) => (
                    <li key={offer.title} className="rounded-full border border-white/10 px-3 py-1 text-xs text-forest-100/70">
                      {offer.title}
                    </li>
                  ))}
                </ul>
                <span className="mt-10 inline-flex items-center gap-2 text-sm font-semibold text-gold-400" aria-hidden="true">
                  {pole.cta}
                  <Icon name="arrowRight" className="h-4 w-4 transition-transform duration-300 group-hover:translate-x-1" />
                </span>
              </article>
            </Reveal>
          ))}
        </div>
      </div>
    </section>
  );
}
