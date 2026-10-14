import { IndustryCard } from "@/components/IndustryCard";
import { Reveal } from "@/components/Reveal";
import { SectionHeader } from "@/components/SectionHeader";
import { sectors } from "@/data/sectors";

export function SectorsSection() {
  return (
    <section className="relative overflow-hidden bg-forest-900 py-24 text-white sm:py-32" aria-labelledby="sectors-title">
      <div className="absolute inset-0 bg-[radial-gradient(ellipse_at_top_left,rgba(21,143,122,0.25),transparent_50%)]" />
      <div className="container-site relative">
        <SectionHeader
          id="sectors-title"
          tone="dark"
          eyebrow="Secteurs"
          title="Nous accompagnons les secteurs les plus exposés."
          intro="Chaque secteur a ses menaces, ses régulateurs et ses contraintes. Nos interventions en tiennent compte dès le premier jour."
        />
        <div className="mt-14 grid gap-4 sm:grid-cols-2 lg:grid-cols-4">
          {sectors.map((sector, index) => (
            <Reveal key={sector.slug} delay={(index % 4) * 0.06} className="h-full">
              <IndustryCard name={sector.name} summary={sector.summary} icon={sector.icon} href={`/secteurs/${sector.slug}`} />
            </Reveal>
          ))}
        </div>
      </div>
    </section>
  );
}
