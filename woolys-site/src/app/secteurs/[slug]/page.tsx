import type { Metadata } from "next";
import Link from "next/link";
import { notFound } from "next/navigation";
import { ButtonLink } from "@/components/Button";
import { CTA } from "@/components/CTA";
import { Icon } from "@/components/Icon";
import { PageHero } from "@/components/PageHero";
import { Reveal } from "@/components/Reveal";
import { SectionHeader } from "@/components/SectionHeader";
import { ServiceCard } from "@/components/ServiceCard";
import { getExpertise } from "@/data/expertises";
import { getSector, sectors } from "@/data/sectors";
import { t } from "@/lib/i18n";
import { pageMetadata } from "@/lib/seo";

export const dynamicParams = false;

export function generateStaticParams() {
  return sectors.map((sector) => ({ slug: sector.slug }));
}

export async function generateMetadata({ params }: PageProps<"/secteurs/[slug]">): Promise<Metadata> {
  const { slug } = await params;
  const sector = getSector(slug);
  if (!sector) return {};
  return pageMetadata({
    title: `Cybersécurité ${sector.name}`,
    description: `${sector.hero} WOOLYS accompagne le secteur ${sector.name.toLowerCase()} au Sénégal et en Afrique : ${sector.summary.toLowerCase()}`,
    path: `/secteurs/${slug}`,
    keywords: [`cybersécurité ${sector.name.toLowerCase()} Sénégal`, `cybersécurité ${sector.name.toLowerCase()} Afrique`],
  });
}

export default async function SectorPage({ params }: PageProps<"/secteurs/[slug]">) {
  const { slug } = await params;
  const sector = getSector(slug);
  if (!sector) notFound();
  const related = sector.expertises.map(getExpertise).filter((item) => item !== undefined);
  const others = sectors.filter((other) => other.slug !== sector.slug);

  return (
    <>
      <PageHero
        breadcrumbs={[
          { label: "Secteurs", href: "/secteurs" },
          { label: sector.name, href: `/secteurs/${sector.slug}` },
        ]}
        eyebrow={`Secteur · ${sector.name}`}
        title={sector.hero}
        subtitle={sector.summary}
        actions={
          <ButtonLink href="/contact" arrow>
            {t.cta.expert}
          </ButtonLink>
        }
      />
      <section className="bg-sand-50 py-24 sm:py-28" aria-labelledby="enjeux-secteur">
        <div className="container-site grid gap-14 lg:grid-cols-2 lg:gap-20">
          <SectionHeader id="enjeux-secteur" eyebrow="Enjeux" title={`Les défis cyber du secteur ${sector.name.toLowerCase()}`} />
          <Reveal>
            <ul className="space-y-4">
              {sector.challenges.map((challenge) => (
                <li key={challenge} className="flex gap-4 rounded-2xl border border-forest-800/10 bg-white p-5">
                  <Icon name={sector.icon} className="mt-0.5 h-5 w-5 shrink-0 text-gold-700" />
                  <span className="text-forest-800">{challenge}</span>
                </li>
              ))}
            </ul>
          </Reveal>
        </div>
      </section>
      <section className="bg-forest-950 py-24 text-white sm:py-28" aria-labelledby="accompagnement">
        <div className="container-site">
          <SectionHeader id="accompagnement" tone="dark" eyebrow="Notre accompagnement" title="Comment WOOLYS peut vous aider" />
          <div className="mt-14 grid gap-6 md:grid-cols-3">
            {sector.support.map((item, index) => (
              <Reveal key={item.title} delay={index * 0.08} className="rounded-2xl border border-white/10 bg-white/[0.03] p-8">
                <span className="font-display text-sm font-semibold text-gold-400">{String(index + 1).padStart(2, "0")}</span>
                <h3 className="mt-4 text-xl font-semibold">{item.title}</h3>
                <p className="mt-3 leading-relaxed text-forest-100/75">{item.text}</p>
              </Reveal>
            ))}
          </div>
          <Reveal className="mt-12 flex flex-wrap items-center gap-3">
            <span className="text-sm text-forest-100/60">Référentiels souvent mobilisés :</span>
            {sector.frameworks.map((framework) => (
              <span key={framework} className="rounded-full border border-gold-400/40 px-4 py-1.5 text-sm text-gold-300">
                {framework}
              </span>
            ))}
          </Reveal>
        </div>
      </section>
      <section className="bg-sand-100 py-24 sm:py-28">
        <div className="container-site">
          <SectionHeader eyebrow="Expertises" title="Expertises recommandées" />
          <div className="mt-12 grid gap-6 md:grid-cols-3">
            {related.map((item) => (
              <ServiceCard key={item.slug} number={item.number} title={item.title} description={item.summary} icon={item.icon} href={`/expertises/${item.slug}`} />
            ))}
          </div>
          <nav aria-label="Autres secteurs" className="mt-16 border-t border-forest-800/10 pt-8">
            <p className="text-sm font-semibold text-forest-800">Autres secteurs</p>
            <ul className="mt-4 flex flex-wrap gap-2">
              {others.map((other) => (
                <li key={other.slug}>
                  <Link href={`/secteurs/${other.slug}`} className="inline-flex rounded-full border border-forest-800/15 bg-white px-4 py-2 text-sm text-forest-800 transition-colors hover:border-gold-400">
                    {other.name}
                  </Link>
                </li>
              ))}
            </ul>
          </nav>
        </div>
      </section>
      <CTA />
    </>
  );
}
