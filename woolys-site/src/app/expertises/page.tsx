import type { Metadata } from "next";
import { ButtonLink } from "@/components/Button";
import { CTA } from "@/components/CTA";
import { PageHero } from "@/components/PageHero";
import { Reveal } from "@/components/Reveal";
import { SectionHeader } from "@/components/SectionHeader";
import { ServiceCard } from "@/components/ServiceCard";
import { expertises } from "@/data/expertises";
import { getPole } from "@/data/poles";
import { t } from "@/lib/i18n";
import { pageMetadata } from "@/lib/seo";
import { ApproachSection } from "@/sections/ApproachSection";

export const metadata: Metadata = pageMetadata({
  title: "Nos expertises en cybersécurité, risques et conformité",
  description:
    "Gouvernance cyber, Risk Management, GRC & Compliance, Cloud Security, AI Security et sensibilisation : les expertises WOOLYS au Sénégal, en Afrique et à l'international.",
  path: "/expertises",
  keywords: ["expertises cybersécurité", "consultant cybersécurité Sénégal"],
});

export default function ExpertisesPage() {
  return (
    <>
      <PageHero
        breadcrumbs={[{ label: "Expertises", href: "/expertises" }]}
        eyebrow="Expertises"
        title="Six expertises pour une cybersécurité pilotée."
        subtitle="De la gouvernance à la sécurité de l'intelligence artificielle, nous couvrons les enjeux qui comptent pour les directions générales, les DSI et les RSSI."
        actions={
          <>
            <ButtonLink href="/contact" arrow>
              {t.cta.expert}
            </ButtonLink>
            <ButtonLink href="/contact?besoin=diagnostic" variant="outline">
              {t.cta.diagnostic}
            </ButtonLink>
          </>
        }
      />
      <section className="bg-sand-50 py-24 sm:py-28" aria-labelledby="liste-expertises">
        <div className="container-site">
          <SectionHeader
            id="liste-expertises"
            eyebrow="Domaines d'intervention"
            title="Choisissez votre point d'entrée"
            intro="Chaque expertise peut être mobilisée seule ou intégrée dans un programme global de transformation."
          />
          <div className="mt-14 grid gap-6 sm:grid-cols-2 lg:grid-cols-3">
            {expertises.map((item, index) => (
              <Reveal key={item.slug} delay={(index % 3) * 0.08} className="h-full">
                <ServiceCard
                  number={`${item.number} · ${getPole(item.pole)?.name ?? ""}`}
                  title={item.title}
                  description={item.summary}
                  icon={item.icon}
                  href={`/expertises/${item.slug}`}
                />
              </Reveal>
            ))}
          </div>
        </div>
      </section>
      <ApproachSection />
      <CTA />
    </>
  );
}
