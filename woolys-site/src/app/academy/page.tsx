import type { Metadata } from "next";
import { ButtonLink } from "@/components/Button";
import { CTA } from "@/components/CTA";
import { Icon, type IconName } from "@/components/Icon";
import { Emblem } from "@/components/Logo";
import { PageHero } from "@/components/PageHero";
import { Reveal } from "@/components/Reveal";
import { SectionHeader } from "@/components/SectionHeader";
import { getPole } from "@/data/poles";
import { t } from "@/lib/i18n";
import { pageMetadata } from "@/lib/seo";
import { siteConfig } from "@/lib/site";
import { PoleOffers } from "@/sections/PoleOffers";

export const metadata: Metadata = pageMetadata({
  title: "WOOLYS Academy — Formation et sensibilisation à la cybersécurité",
  description:
    "Formations cybersécurité au Sénégal et en Afrique : sensibilisation, Executive Awareness, EBIOS RM, ISO 27001, GRC, Cloud Security, AI Security et plateforme e-learning.",
  path: "/academy",
  keywords: ["formation cybersécurité Sénégal", "formation EBIOS RM", "sensibilisation cybersécurité Afrique", "e-learning cybersécurité"],
});

const audiences: { icon: IconName; title: string; text: string }[] = [
  { icon: "governance", title: "Dirigeants", text: "Comprendre le risque cyber, piloter et décider en situation de crise." },
  { icon: "people", title: "Collaborateurs", text: "Adopter les bons réflexes au quotidien face aux menaces courantes." },
  { icon: "compliance", title: "Risk & Compliance", text: "Maîtriser les méthodes d'analyse de risques et les référentiels." },
  { icon: "chip", title: "Équipes IT & sécurité", text: "Monter en compétences sur le Cloud, l'IA et la sécurité opérationnelle." },
];

const formats = [
  { title: "Présentiel", text: "Sessions en intra-entreprise, dans vos locaux." },
  { title: "Distanciel", text: "Classes virtuelles interactives." },
  { title: "E-learning", text: "Parcours en ligne sur la plateforme WOOLYS Academy." },
  { title: "Ateliers & exercices", text: "Cas pratiques, simulations d'hameçonnage et exercices de crise." },
];

export default function AcademyPage() {
  const pole = getPole("academy")!;
  return (
    <>
      <PageHero
        breadcrumbs={[{ label: "WOOLYS Academy", href: "/academy" }]}
        eyebrow="Pôle 02 · WOOLYS Academy"
        title="Développer les compétences pour renforcer la sécurité."
        subtitle={pole.summary}
        actions={
          <>
            <ButtonLink href="/contact?besoin=formation#formulaire" variant="gold" arrow>
              {t.cta.training}
            </ButtonLink>
            {siteConfig.lmsUrl ? (
              <ButtonLink href={siteConfig.lmsUrl} variant="outline">
                {t.cta.lms}
              </ButtonLink>
            ) : (
              <ButtonLink href="#plateforme" variant="outline">
                Découvrir la plateforme
              </ButtonLink>
            )}
          </>
        }
      />

      <section className="bg-white py-24 sm:py-28" aria-labelledby="publics">
        <div className="container-site">
          <SectionHeader id="publics" eyebrow="Pour qui ?" title="Des parcours pour chaque public." />
          <div className="mt-14 grid gap-6 sm:grid-cols-2 lg:grid-cols-4">
            {audiences.map((audience, index) => (
              <Reveal key={audience.title} delay={index * 0.06} className="rounded-2xl border border-forest-800/10 bg-sand-50 p-7">
                <Icon name={audience.icon} className="h-7 w-7 text-emerald-600" />
                <h3 className="mt-5 text-lg font-semibold text-forest-800">{audience.title}</h3>
                <p className="mt-2 text-sm leading-relaxed text-muted">{audience.text}</p>
              </Reveal>
            ))}
          </div>
        </div>
      </section>

      <PoleOffers eyebrow="Catalogue" title="Formations, sensibilisation et ateliers" offers={pole.offers} />

      {/* Plateforme LMS */}
      <section id="plateforme" className="relative overflow-hidden bg-forest-900 py-24 text-white sm:py-28" aria-labelledby="lms-title">
        <div className="absolute inset-0 bg-grid opacity-50" />
        <div className="container-site relative grid items-center gap-14 lg:grid-cols-2">
          <Reveal>
            <p className="eyebrow mb-5 text-gold-400">Plateforme e-learning</p>
            <h2 id="lms-title" className="text-3xl font-bold leading-tight sm:text-4xl">
              Apprendre à son rythme, sur la plateforme WOOLYS Academy.
            </h2>
            <p className="mt-6 text-lg leading-relaxed text-forest-100/80">
              Notre plateforme de formation en ligne prolonge les sessions en présentiel : modules interactifs, parcours par
              public, quiz et suivi de progression pour les organisations.
            </p>
            <ul className="mt-8 space-y-3">
              {["Parcours de sensibilisation pour tous les collaborateurs", "Modules experts : EBIOS RM, GRC, Cloud, IA", "Suivi de la progression des équipes"].map((feature) => (
                <li key={feature} className="flex items-start gap-3 text-forest-100/85">
                  <Icon name="check" className="mt-0.5 h-5 w-5 shrink-0 text-gold-400" />
                  {feature}
                </li>
              ))}
            </ul>
            <div className="mt-10 flex flex-col gap-4 sm:flex-row">
              {siteConfig.lmsUrl ? (
                <ButtonLink href={siteConfig.lmsUrl} variant="gold" arrow>
                  {t.cta.lms}
                </ButtonLink>
              ) : (
                <p className="rounded-full border border-gold-400/40 px-5 py-3 text-sm text-gold-300">Accès en ligne bientôt disponible</p>
              )}
              <ButtonLink href="/contact?besoin=formation#formulaire" variant="outline">
                Demander un accès pour mon organisation
              </ButtonLink>
            </div>
          </Reveal>
          <Reveal delay={0.1}>
            <div className="relative rounded-3xl border border-white/10 bg-forest-950/80 p-6 shadow-2xl">
              <div className="flex items-center gap-2 border-b border-white/10 pb-4">
                <span className="h-2.5 w-2.5 rounded-full bg-white/20" />
                <span className="h-2.5 w-2.5 rounded-full bg-white/20" />
                <span className="h-2.5 w-2.5 rounded-full bg-white/20" />
                <span className="ml-3 text-xs text-forest-100/50">WOOLYS Academy</span>
              </div>
              <div className="mt-6 flex items-center gap-4">
                <Emblem className="h-12 w-12" />
                <div>
                  <p className="font-display font-semibold">Mon parcours</p>
                  <p className="text-xs text-forest-100/60">Aperçu illustratif</p>
                </div>
              </div>
              <ul className="mt-6 space-y-3">
                {[
                  { name: "Les fondamentaux de la cybersécurité", progress: 100 },
                  { name: "Reconnaître l'hameçonnage", progress: 70 },
                  { name: "Introduction à EBIOS RM", progress: 35 },
                  { name: "Utiliser l'IA générative en sécurité", progress: 0 },
                ].map((module) => (
                  <li key={module.name} className="rounded-xl bg-white/[0.04] p-4">
                    <div className="flex items-center justify-between gap-4 text-sm">
                      <span>{module.name}</span>
                      <span className="text-xs text-gold-300">{module.progress}%</span>
                    </div>
                    <div className="mt-3 h-1.5 rounded-full bg-white/10">
                      <div className="h-full rounded-full bg-gradient-to-r from-gold-400 to-emerald-400" style={{ width: `${module.progress}%` }} />
                    </div>
                  </li>
                ))}
              </ul>
            </div>
          </Reveal>
        </div>
      </section>

      <section className="bg-sand-100 py-24 sm:py-28" aria-labelledby="formats">
        <div className="container-site">
          <SectionHeader id="formats" eyebrow="Formats" title="Des formats adaptés à vos contraintes." />
          <div className="mt-12 grid gap-6 sm:grid-cols-2 lg:grid-cols-4">
            {formats.map((format, index) => (
              <Reveal key={format.title} delay={index * 0.06} className="rounded-2xl border-t-2 border-gold-400 bg-white p-7">
                <h3 className="text-lg font-semibold text-forest-800">{format.title}</h3>
                <p className="mt-2 text-sm leading-relaxed text-muted">{format.text}</p>
              </Reveal>
            ))}
          </div>
        </div>
      </section>

      <CTA
        title="Faisons de vos équipes votre première ligne de défense."
        text="Sensibilisation, formation des dirigeants ou parcours experts : construisons ensemble le programme adapté à votre organisation."
        primaryLabel={t.cta.training}
        primaryHref="/contact?besoin=formation#formulaire"
      />
    </>
  );
}
