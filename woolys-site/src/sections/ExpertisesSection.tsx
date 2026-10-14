import { ButtonLink } from "@/components/Button";
import { Reveal } from "@/components/Reveal";
import { SectionHeader } from "@/components/SectionHeader";
import { ServiceCard } from "@/components/ServiceCard";
import { expertises } from "@/data/expertises";

export function ExpertisesSection() {
  return (
    <section className="bg-sand-100 py-24 sm:py-32" aria-labelledby="expertises-title">
      <div className="container-site">
        <div className="flex flex-col justify-between gap-8 lg:flex-row lg:items-end">
          <SectionHeader
            id="expertises-title"
            eyebrow="Expertises"
            title="Nos expertises"
            intro="Six domaines d'intervention pour couvrir l'ensemble de vos enjeux de cybersécurité, de la gouvernance aux technologies émergentes."
          />
          <Reveal>
            <ButtonLink href="/expertises" variant="outlineDark" arrow>
              Toutes nos expertises
            </ButtonLink>
          </Reveal>
        </div>
        <div className="mt-14 grid gap-6 sm:grid-cols-2 lg:grid-cols-3">
          {expertises.map((item, index) => (
            <Reveal key={item.slug} delay={(index % 3) * 0.08} className="h-full">
              <ServiceCard
                number={item.number}
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
  );
}
