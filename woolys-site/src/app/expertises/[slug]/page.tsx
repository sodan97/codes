import type { Metadata } from "next";
import Link from "next/link";
import { notFound } from "next/navigation";
import { ButtonLink } from "@/components/Button";
import { CTA } from "@/components/CTA";
import { FAQ } from "@/components/FAQ";
import { Icon } from "@/components/Icon";
import { JsonLd } from "@/components/JsonLd";
import { PageHero } from "@/components/PageHero";
import { Reveal } from "@/components/Reveal";
import { SectionHeader } from "@/components/SectionHeader";
import { ServiceCard } from "@/components/ServiceCard";
import { Timeline } from "@/components/Timeline";
import { expertises, getExpertise } from "@/data/expertises";
import { getMethodology, methodologyLabel } from "@/data/methodologies";
import { getPole } from "@/data/poles";
import { t } from "@/lib/i18n";
import { pageMetadata } from "@/lib/seo";
import { absoluteUrl, bookingHref, siteConfig } from "@/lib/site";

export const dynamicParams = false;

export function generateStaticParams() {
  return expertises.map((item) => ({ slug: item.slug }));
}

export async function generateMetadata({ params }: PageProps<"/expertises/[slug]">): Promise<Metadata> {
  const { slug } = await params;
  const item = getExpertise(slug);
  if (!item) return {};
  return pageMetadata({ title: item.seo.title, description: item.seo.description, path: `/expertises/${slug}`, keywords: item.seo.keywords });
}

export default async function ExpertisePage({ params }: PageProps<"/expertises/[slug]">) {
  const { slug } = await params;
  const item = getExpertise(slug);
  if (!item) notFound();
  const pole = getPole(item.pole);
  const others = expertises.filter((other) => other.slug !== item.slug).slice(0, 3);

  return (
    <>
      <PageHero
        breadcrumbs={[
          { label: "Expertises", href: "/expertises" },
          { label: item.title, href: `/expertises/${item.slug}` },
        ]}
        eyebrow={`${item.number} · ${item.title}`}
        title={item.hero.title}
        subtitle={item.hero.subtitle}
        actions={
          <>
            <ButtonLink href={`/contact?besoin=${item.slug}#formulaire`} arrow>
              {t.cta.expert}
            </ButtonLink>
            <ButtonLink href="/contact?besoin=diagnostic" variant="outline">
              {t.cta.diagnostic}
            </ButtonLink>
          </>
        }
        aside={
          <div className="rounded-3xl border border-white/10 bg-white/[0.04] p-8 backdrop-blur">
            <span className="flex h-14 w-14 items-center justify-center rounded-2xl bg-gold-400 text-forest-950">
              <Icon name={item.icon} className="h-7 w-7" />
            </span>
            <p className="mt-6 text-sm uppercase tracking-[0.2em] text-gold-400">Pôle</p>
            <p className="mt-1 font-display text-xl font-semibold">{pole?.name}</p>
            <p className="mt-3 text-sm leading-relaxed text-forest-100/70">{pole?.tagline}</p>
          </div>
        }
      />

      {/* Problématique */}
      <section className="bg-sand-50 py-24 sm:py-28" aria-labelledby="problematique">
        <div className="container-site grid gap-14 lg:grid-cols-2 lg:gap-20">
          <SectionHeader id="problematique" eyebrow="Problématique" title={item.problem.title} intro={item.problem.intro} />
          <Reveal>
            <ul className="space-y-4">
              {item.problem.points.map((point) => (
                <li key={point} className="flex gap-4 rounded-2xl border border-forest-800/10 bg-white p-5">
                  <span className="mt-0.5 flex h-6 w-6 shrink-0 items-center justify-center rounded-full bg-gold-400/20 text-gold-700">
                    <Icon name="target" className="h-4 w-4" />
                  </span>
                  <span className="text-forest-800">{point}</span>
                </li>
              ))}
            </ul>
          </Reveal>
        </div>
      </section>

      {/* Approche */}
      <section className="bg-forest-950 py-24 text-white sm:py-28" aria-labelledby="approche">
        <div className="container-site">
          <SectionHeader id="approche" tone="dark" eyebrow="Notre approche" title="Une démarche structurée, adaptée à votre contexte." />
          <div className="mt-16">
            <Timeline tone="dark" steps={item.approach} />
          </div>
        </div>
      </section>

      {/* Services */}
      <section className="bg-sand-100 py-24 sm:py-28" aria-labelledby="services">
        <div className="container-site">
          <SectionHeader id="services" eyebrow="Services" title="Ce que nous faisons pour vous" />
          <div className="mt-14 grid gap-6 sm:grid-cols-2 lg:grid-cols-3">
            {item.services.map((service, index) => (
              <Reveal key={service.title} delay={(index % 3) * 0.06} className="h-full">
                <ServiceCard title={service.title} description={service.text} number={String(index + 1).padStart(2, "0")} />
              </Reveal>
            ))}
          </div>
        </div>
      </section>

      {/* Méthodologies & livrables */}
      <section className="bg-white py-24 sm:py-28">
        <div className="container-site grid gap-16 lg:grid-cols-2">
          <div>
            <SectionHeader eyebrow="Méthodologies" title="Référentiels mobilisés" />
            <Reveal>
              <ul className="mt-10 flex flex-wrap gap-3">
                {item.methodologies.map((entry) => {
                  const method = getMethodology(entry);
                  return (
                    <li key={entry}>
                      {method ? (
                        <Link
                          href={`/methodologies#${method.slug}`}
                          className="inline-flex rounded-full border border-forest-800/15 bg-sand-50 px-5 py-2.5 font-display text-sm font-semibold text-forest-800 transition-colors hover:border-gold-400 hover:bg-white"
                        >
                          {method.name}
                        </Link>
                      ) : (
                        <span className="inline-flex rounded-full border border-forest-800/15 bg-sand-50 px-5 py-2.5 font-display text-sm font-semibold text-forest-800">
                          {methodologyLabel(entry)}
                        </span>
                      )}
                    </li>
                  );
                })}
              </ul>
              <p className="mt-6 text-sm text-muted">
                Référentiels utilisés ou pris en compte dans nos missions, sélectionnés selon votre contexte.
              </p>
            </Reveal>
          </div>
          <div>
            <SectionHeader eyebrow="Livrables" title="Ce que vous obtenez" />
            <Reveal>
              <ul className="mt-10 space-y-3">
                {item.deliverables.map((deliverable) => (
                  <li key={deliverable} className="flex items-start gap-3 border-b border-forest-800/10 pb-3">
                    <Icon name="check" className="mt-0.5 h-5 w-5 shrink-0 text-emerald-600" />
                    <span className="text-forest-800">{deliverable}</span>
                  </li>
                ))}
              </ul>
            </Reveal>
          </div>
        </div>
      </section>

      {/* FAQ */}
      <section className="bg-sand-50 py-24 sm:py-28" aria-labelledby="faq">
        <div className="container-site grid gap-12 lg:grid-cols-[0.8fr_1.2fr]">
          <SectionHeader
            id="faq"
            eyebrow="FAQ"
            title="Questions fréquentes"
            intro={
              <>
                Une autre question ?{" "}
                <Link href={bookingHref} className="font-semibold text-emerald-600 underline underline-offset-4 hover:text-forest-800">
                  {t.cta.meeting}
                </Link>
              </>
            }
          />
          <Reveal>
            <FAQ items={item.faq} />
          </Reveal>
        </div>
      </section>

      {/* Autres expertises */}
      <section className="bg-sand-100 py-20" aria-labelledby="autres">
        <div className="container-site">
          <h2 id="autres" className="text-2xl font-bold text-forest-800">
            Expertises complémentaires
          </h2>
          <div className="mt-10 grid gap-6 md:grid-cols-3">
            {others.map((other) => (
              <ServiceCard key={other.slug} title={other.title} description={other.summary} icon={other.icon} href={`/expertises/${other.slug}`} />
            ))}
          </div>
        </div>
      </section>

      <CTA
        title="Parlons de votre projet."
        text={`Un besoin en ${item.title} ? Nos experts vous aident à définir une réponse adaptée à vos enjeux et à vos priorités.`}
        primaryLabel={t.cta.proposal}
        primaryHref="/contact?besoin=proposition"
      />

      <JsonLd
        data={{
          "@context": "https://schema.org",
          "@type": "Service",
          name: item.title,
          description: item.summary,
          serviceType: item.title,
          url: absoluteUrl(`/expertises/${item.slug}`),
          provider: { "@id": `${siteConfig.url}/#organization` },
          areaServed: siteConfig.areaServed,
        }}
      />
    </>
  );
}
