import { ButtonLink } from "@/components/Button";
import { MethodologyCard } from "@/components/MethodologyCard";
import { Reveal } from "@/components/Reveal";
import { SectionHeader } from "@/components/SectionHeader";
import { methodologies } from "@/data/methodologies";

export function MethodologiesSection() {
  return (
    <section className="relative overflow-hidden bg-white py-24 sm:py-32" aria-labelledby="methodo-title">
      <div className="absolute inset-0 bg-grid-light [mask-image:linear-gradient(to_bottom,black,transparent)]" />
      <div className="container-site relative">
        <SectionHeader
          id="methodo-title"
          align="center"
          eyebrow="Méthodologies"
          title={
            <>
              Des méthodes reconnues.
              <br />
              Des résultats concrets.
            </>
          }
          intro="Nos missions s'appuient sur des référentiels et méthodes reconnus internationalement, sélectionnés selon votre contexte, votre secteur et vos obligations."
        />
        <ul className="mx-auto mt-14 grid max-w-5xl grid-cols-2 gap-4 sm:grid-cols-4">
          {methodologies.map((item, index) => (
            <Reveal as="li" key={item.slug} delay={index * 0.05} className="h-full">
              <MethodologyCard name={item.name} category={item.category} href={`/methodologies#${item.slug}`} />
            </Reveal>
          ))}
        </ul>
        <Reveal className="mx-auto mt-10 max-w-3xl text-center">
          <p className="text-sm text-muted">
            Ces référentiels sont utilisés ou pris en compte dans nos missions. Les certifications éventuelles de votre
            organisation sont délivrées par des organismes accrédités indépendants.
          </p>
          <ButtonLink href="/methodologies" variant="ghost" arrow className="mt-4">
            Comprendre nos méthodologies
          </ButtonLink>
        </Reveal>
      </div>
    </section>
  );
}
