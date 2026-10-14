import type { Metadata } from "next";
import { CTA } from "@/components/CTA";
import { IndustryCard } from "@/components/IndustryCard";
import { PageHero } from "@/components/PageHero";
import { Reveal } from "@/components/Reveal";
import { sectors } from "@/data/sectors";
import { pageMetadata } from "@/lib/seo";

export const metadata: Metadata = pageMetadata({
  title: "Secteurs accompagnés — Finance, télécoms, assurance, secteur public…",
  description:
    "WOOLYS accompagne les secteurs les plus exposés au risque cyber au Sénégal et en Afrique : finance, télécoms, assurance, secteur public, énergie, industrie, technologie et infrastructures critiques.",
  path: "/secteurs",
});

export default function SectorsPage() {
  return (
    <>
      <PageHero
        breadcrumbs={[{ label: "Secteurs", href: "/secteurs" }]}
        eyebrow="Secteurs"
        title="Nous accompagnons les secteurs les plus exposés."
        subtitle="Chaque secteur a ses menaces, ses régulateurs et ses contraintes opérationnelles. Nos interventions en tiennent compte dès le premier jour."
      />
      <section className="bg-forest-950 py-24 sm:py-28" aria-label="Liste des secteurs">
        <div className="container-site grid gap-4 sm:grid-cols-2 lg:grid-cols-4">
          {sectors.map((sector, index) => (
            <Reveal key={sector.slug} delay={(index % 4) * 0.06} className="h-full">
              <IndustryCard name={sector.name} summary={sector.summary} icon={sector.icon} href={`/secteurs/${sector.slug}`} />
            </Reveal>
          ))}
        </div>
      </section>
      <CTA title="Votre secteur a ses spécificités. Parlons-en." />
    </>
  );
}
