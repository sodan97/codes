import type { Metadata } from "next";
import { ButtonLink } from "@/components/Button";
import { CTA } from "@/components/CTA";
import { Icon, type IconName } from "@/components/Icon";
import { Emblem } from "@/components/Logo";
import { PageHero } from "@/components/PageHero";
import { Reveal } from "@/components/Reveal";
import { SectionHeader } from "@/components/SectionHeader";
import { t } from "@/lib/i18n";
import { pageMetadata } from "@/lib/seo";
import { ApproachSection } from "@/sections/ApproachSection";

export const metadata: Metadata = pageMetadata({
  title: "À propos — Construire une confiance numérique durable",
  description:
    "WOOLYS est un cabinet de conseil et de services numériques spécialisé dans la cybersécurité, la gestion des risques et la confiance numérique, au Sénégal, en Afrique et à l'international.",
  path: "/a-propos",
});

const values: { icon: IconName; title: string; text: string }[] = [
  { icon: "target", title: "Excellence", text: "Un niveau d'exigence élevé, des livrables rigoureux et des méthodes reconnues." },
  { icon: "governance", title: "Intégrité", text: "Indépendance de jugement, confidentialité et transparence dans nos recommandations." },
  { icon: "sparkle", title: "Innovation", text: "Une veille active sur le Cloud, l'IA et les nouvelles menaces, portée par WOOLYS Labs." },
  { icon: "people", title: "Proximité", text: "Une compréhension fine des réalités africaines et un accompagnement au plus près des équipes." },
  { icon: "book", title: "Transmission", text: "Le transfert de compétences au cœur de chaque mission, prolongé par WOOLYS Academy." },
];

export default function AboutPage() {
  return (
    <>
      <PageHero
        breadcrumbs={[{ label: "À propos", href: "/a-propos" }]}
        eyebrow="À propos de WOOLYS"
        title="Construire une confiance numérique durable."
        subtitle="WOOLYS est un cabinet de conseil et de services numériques spécialisé dans la cybersécurité, la gestion des risques et la confiance numérique."
        aside={
          <div className="flex justify-center">
            <Emblem className="h-56 w-56 sm:h-64 sm:w-64" title="Emblème WOOLYS" />
          </div>
        }
      />

      <section className="bg-sand-50 py-24 sm:py-28">
        <div className="container-site grid gap-8 lg:grid-cols-2">
          <Reveal className="rounded-3xl border border-forest-800/10 bg-white p-10">
            <p className="eyebrow text-gold-700">Mission</p>
            <h2 className="mt-5 text-2xl font-bold text-forest-800 sm:text-3xl">Aider les organisations à maîtriser leurs risques numériques.</h2>
            <p className="mt-5 leading-relaxed text-muted">
              Nous accompagnons les directions générales, les DSI, les RSSI et les fonctions risques et conformité pour structurer
              leur gouvernance, évaluer leurs risques, répondre à leurs obligations et développer les compétences de leurs équipes.
            </p>
          </Reveal>
          <Reveal delay={0.08} className="rounded-3xl bg-forest-800 p-10 text-white">
            <p className="eyebrow text-gold-400">Vision</p>
            <h2 className="mt-5 text-2xl font-bold sm:text-3xl">Un écosystème numérique africain plus sûr, plus résilient et plus compétitif.</h2>
            <p className="mt-5 leading-relaxed text-forest-100/80">
              Nous voulons contribuer à faire de la confiance numérique un levier de développement, au Sénégal, en Afrique de
              l&apos;Ouest et au-delà, en conjuguant standards internationaux et compréhension des réalités locales.
            </p>
          </Reveal>
        </div>
      </section>

      <section className="bg-white py-24 sm:py-28" aria-labelledby="valeurs">
        <div className="container-site">
          <SectionHeader id="valeurs" eyebrow="Valeurs" title="Ce qui guide chacune de nos missions." />
          <ul className="mt-14 grid gap-6 sm:grid-cols-2 lg:grid-cols-5">
            {values.map((value, index) => (
              <Reveal as="li" key={value.title} delay={index * 0.06} className="rounded-2xl border border-forest-800/10 bg-sand-50 p-7">
                <Icon name={value.icon} className="h-7 w-7 text-emerald-600" />
                <h3 className="mt-5 text-lg font-semibold text-forest-800">{value.title}</h3>
                <p className="mt-2 text-sm leading-relaxed text-muted">{value.text}</p>
              </Reveal>
            ))}
          </ul>
        </div>
      </section>

      <ApproachSection tone="dark" />

      {/* PLACEHOLDER : section équipe / fondateurs à ajouter avec des informations validées (noms, parcours, photos). */}
      <section className="bg-sand-100 py-20">
        <div className="container-site flex flex-col items-start justify-between gap-8 rounded-3xl border border-forest-800/10 bg-white p-10 lg:flex-row lg:items-center">
          <div className="max-w-2xl">
            <h2 className="text-2xl font-bold text-forest-800">Trois pôles, une même exigence.</h2>
            <p className="mt-3 text-muted">
              WOOLYS Advisory conseille, WOOLYS Academy transmet, WOOLYS Labs explore. Ensemble, ils couvrent tout le cycle de la
              confiance numérique.
            </p>
          </div>
          <ButtonLink href="/expertises" variant="outlineDark" arrow>
            {t.cta.discover}
          </ButtonLink>
        </div>
      </section>

      <CTA />
    </>
  );
}
