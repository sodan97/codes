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
import { PoleOffers } from "@/sections/PoleOffers";

export const metadata: Metadata = pageMetadata({
  title: "WOOLYS Advisory — Conseil en cybersécurité, GRC et gestion des risques",
  description:
    "Conseil en cybersécurité au Sénégal et en Afrique : gouvernance, GRC, Risk Management, EBIOS RM, ISO 27001, conformité, DORA, PCI-DSS et accompagnement RSSI.",
  path: "/advisory",
  keywords: ["cabinet cybersécurité Sénégal", "conseil cybersécurité Afrique", "cybersecurity consulting Africa"],
});

const audiences = ["CEO & Directions générales", "DSI / CIO / CTO", "RSSI / CISO", "Risk Managers", "Compliance Officers", "Institutions publiques"];

export default function AdvisoryPage() {
  const pole = getPole("advisory")!;
  const related = expertises.filter((item) => item.pole === "advisory");
  return (
    <>
      <PageHero
        breadcrumbs={[{ label: "WOOLYS Advisory", href: "/advisory" }]}
        eyebrow="Pôle 01 · WOOLYS Advisory"
        title="Conseil, gouvernance et gestion des risques cyber."
        subtitle={pole.summary}
        actions={
          <>
            <ButtonLink href="/contact?besoin=advisory#formulaire" arrow>
              {t.cta.project}
            </ButtonLink>
            <ButtonLink href="/contact?besoin=diagnostic#formulaire" variant="outline">
              {t.cta.diagnostic}
            </ButtonLink>
          </>
        }
      />
      <section className="bg-white py-20" aria-labelledby="interlocuteurs">
        <div className="container-site grid gap-12 lg:grid-cols-[1fr_1.2fr] lg:items-center">
          <SectionHeader
            id="interlocuteurs"
            eyebrow="Pour qui ?"
            title="Un partenaire pour les décideurs."
            intro="Nous parlons le langage des directions comme celui des équipes techniques, pour transformer les enjeux cyber en décisions."
          />
          <Reveal>
            <ul className="grid grid-cols-2 gap-3 sm:grid-cols-3">
              {audiences.map((audience) => (
                <li key={audience} className="rounded-2xl border border-forest-800/10 bg-sand-50 p-5 text-center font-display text-sm font-semibold text-forest-800">
                  {audience}
                </li>
              ))}
            </ul>
          </Reveal>
        </div>
      </section>
      <PoleOffers eyebrow="Nos offres" title="Une expertise complète en conseil cyber" offers={pole.offers} />
      <section className="bg-sand-100 py-24 sm:py-28">
        <div className="container-site">
          <SectionHeader eyebrow="Expertises associées" title="Approfondir" />
          <div className="mt-12 grid gap-6 md:grid-cols-3">
            {related.map((item) => (
              <ServiceCard key={item.slug} number={item.number} title={item.title} description={item.summary} icon={item.icon} href={`/expertises/${item.slug}`} />
            ))}
          </div>
        </div>
      </section>
      <ApproachSection tone="dark" />
      <CTA title="Discutons de votre projet." primaryLabel={t.cta.project} primaryHref="/contact?besoin=advisory#formulaire" />
    </>
  );
}
