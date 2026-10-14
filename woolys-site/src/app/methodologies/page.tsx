import type { Metadata } from "next";
import Link from "next/link";
import { CTA } from "@/components/CTA";
import { Icon } from "@/components/Icon";
import { PageHero } from "@/components/PageHero";
import { Reveal } from "@/components/Reveal";
import { getExpertise } from "@/data/expertises";
import { methodologies } from "@/data/methodologies";
import { pageMetadata } from "@/lib/seo";

export const metadata: Metadata = pageMetadata({
  title: "Méthodologies — EBIOS RM, ISO 27001, ISO 27005, NIST CSF, DORA…",
  description:
    "Les référentiels utilisés dans les missions WOOLYS : EBIOS RM, ISO 27001, ISO 27005, NIST CSF, CIS Controls, DORA, PCI-DSS et AWS Security. À quoi servent-ils et quand les utiliser ?",
  path: "/methodologies",
  keywords: ["EBIOS RM Sénégal", "ISO 27001 Sénégal", "NIST CSF", "DORA", "PCI-DSS"],
});

export default function MethodologiesPage() {
  return (
    <>
      <PageHero
        breadcrumbs={[{ label: "Méthodologies", href: "/methodologies" }]}
        eyebrow="Méthodologies"
        title={
          <>
            Des méthodes reconnues.
            <br />
            Des résultats concrets.
          </>
        }
        subtitle="Nos missions s'appuient sur des référentiels et des méthodes reconnus. Ils sont utilisés ou pris en compte selon votre contexte : WOOLYS ne se présente pas comme certifiée sur ces référentiels."
      />
      <nav aria-label="Sommaire des méthodologies" className="sticky top-[4.5rem] z-30 border-b border-forest-800/10 bg-sand-50/95 backdrop-blur">
        <ul className="container-site flex gap-2 overflow-x-auto py-3 [scrollbar-width:none]">
          {methodologies.map((item) => (
            <li key={item.slug} className="shrink-0">
              <a href={`#${item.slug}`} className="inline-flex rounded-full px-4 py-2 text-sm font-semibold text-forest-800 transition-colors hover:bg-white hover:text-emerald-600">
                {item.name}
              </a>
            </li>
          ))}
        </ul>
      </nav>
      <div className="bg-sand-50">
        {methodologies.map((item, index) => (
          <section key={item.slug} id={item.slug} aria-labelledby={`${item.slug}-title`} className={index % 2 === 1 ? "bg-white py-20 sm:py-24" : "py-20 sm:py-24"}>
            <div className="container-site grid gap-10 lg:grid-cols-[0.9fr_1.6fr] lg:gap-16">
              <Reveal>
                <p className="text-xs font-semibold uppercase tracking-[0.2em] text-gold-700">{item.category}</p>
                <h2 id={`${item.slug}-title`} className="mt-3 text-4xl font-bold text-forest-800">
                  {item.name}
                </h2>
                <p className="mt-3 text-sm text-muted">{item.fullName}</p>
                <ul className="mt-6 flex flex-wrap gap-2">
                  {item.expertises.map((slug) => {
                    const expertise = getExpertise(slug);
                    return expertise ? (
                      <li key={slug}>
                        <Link href={`/expertises/${slug}`} className="inline-flex rounded-full border border-forest-800/15 px-3 py-1 text-xs font-medium text-emerald-600 hover:border-emerald-600">
                          {expertise.title}
                        </Link>
                      </li>
                    ) : null;
                  })}
                </ul>
              </Reveal>
              <Reveal delay={0.08} className="grid gap-8 md:grid-cols-3">
                <div>
                  <h3 className="font-display text-base font-semibold text-forest-800">À quoi ça sert ?</h3>
                  <p className="mt-3 text-[0.95rem] leading-relaxed text-muted">{item.purpose}</p>
                </div>
                <div>
                  <h3 className="font-display text-base font-semibold text-forest-800">Quand l&apos;utiliser ?</h3>
                  <ul className="mt-3 space-y-2">
                    {item.whenToUse.map((entry) => (
                      <li key={entry} className="flex gap-2 text-[0.95rem] leading-relaxed text-muted">
                        <span className="mt-2.5 h-px w-3 shrink-0 bg-gold-500" aria-hidden="true" />
                        {entry}
                      </li>
                    ))}
                  </ul>
                </div>
                <div className="rounded-2xl bg-forest-800 p-6 text-white">
                  <h3 className="font-display text-base font-semibold text-gold-300">Comment WOOLYS peut accompagner ?</h3>
                  <ul className="mt-3 space-y-2">
                    {item.howWeHelp.map((entry) => (
                      <li key={entry} className="flex gap-2 text-sm leading-relaxed text-forest-100/85">
                        <Icon name="check" className="mt-0.5 h-4 w-4 shrink-0 text-gold-400" />
                        {entry}
                      </li>
                    ))}
                  </ul>
                </div>
              </Reveal>
            </div>
          </section>
        ))}
      </div>
      <CTA title="Quel référentiel pour votre organisation ?" text="Nous vous aidons à choisir et à mettre en œuvre les méthodes adaptées à votre contexte et à vos obligations." />
    </>
  );
}
