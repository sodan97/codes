import type { Metadata } from "next";
import Link from "next/link";
import { ButtonLink } from "@/components/Button";
import { CTA } from "@/components/CTA";
import { Icon } from "@/components/Icon";
import { PageHero } from "@/components/PageHero";
import { Reveal } from "@/components/Reveal";
import { SectionHeader } from "@/components/SectionHeader";
import { getPole } from "@/data/poles";
import { t } from "@/lib/i18n";
import { pageMetadata } from "@/lib/seo";
import { PoleOffers } from "@/sections/PoleOffers";

export const metadata: Metadata = pageMetadata({
  title: "WOOLYS Labs — Cloud Security, AWS et AI Security",
  description:
    "WOOLYS Labs explore les nouvelles frontières de la cybersécurité : Cloud Security, AWS, AI Security, LLM, RAG, agents IA, identité numérique et automatisation.",
  path: "/labs",
  keywords: ["AI security Africa", "cloud security Africa", "LLM security", "AWS security Sénégal"],
});

const tracks = [
  {
    code: "CLOUD",
    title: "Cloud & AWS Security",
    text: "Landing zones sécurisées, identités, journalisation et conformité continue sur AWS.",
    href: "/expertises/cloud-security",
  },
  {
    code: "AI",
    title: "AI, LLM & RAG Security",
    text: "Gouvernance de l'IA, protection contre l'injection de prompt et cloisonnement des données.",
    href: "/expertises/ai-security",
  },
  {
    code: "AGENT",
    title: "Agentic AI Security",
    text: "Permissions, supervision humaine et traçabilité des agents autonomes.",
    href: "/expertises/ai-security",
  },
];

export default function LabsPage() {
  const pole = getPole("labs")!;
  return (
    <>
      <PageHero
        variant="labs"
        breadcrumbs={[{ label: "WOOLYS Labs", href: "/labs" }]}
        eyebrow="Pôle 03 · WOOLYS Labs"
        title={
          <>
            Explorer les nouvelles frontières <span className="text-gradient-gold">de la cybersécurité.</span>
          </>
        }
        subtitle={pole.summary}
        actions={
          <>
            <ButtonLink href="/contact?besoin=labs#formulaire" arrow>
              {t.cta.expert}
            </ButtonLink>
            <ButtonLink href="/insights" variant="outline">
              Lire nos analyses
            </ButtonLink>
          </>
        }
        aside={
          <div className="rounded-2xl border border-emerald-400/20 bg-black/40 p-6 font-mono text-xs leading-relaxed text-emerald-300/90 shadow-2xl">
            <p className="text-forest-100/40">{"// woolys-labs · posture.yaml"}</p>
            <p className="mt-3">
              <span className="text-gold-400">cloud</span>: least_privilege · logging · guardrails
            </p>
            <p>
              <span className="text-gold-400">ai</span>: prompt_injection_defense · data_isolation
            </p>
            <p>
              <span className="text-gold-400">agents</span>: human_in_the_loop · scoped_tools · audit_trail
            </p>
            <p>
              <span className="text-gold-400">identity</span>: mfa · zero_trust
            </p>
            <p className="mt-3 text-forest-100/40">
              status: <span className="text-emerald-300">exploring</span>
              <span className="ml-1 inline-block h-3 w-1.5 animate-glow bg-emerald-300 align-middle" />
            </p>
          </div>
        }
      />

      <section className="bg-[#03130F] pb-8 pt-24 text-white sm:pt-28" aria-labelledby="axes">
        <div className="container-site">
          <SectionHeader id="axes" tone="dark" eyebrow="Axes de recherche" title="Trois terrains d'exploration prioritaires." />
          <div className="mt-14 grid gap-6 lg:grid-cols-3">
            {tracks.map((track, index) => (
              <Reveal key={track.title} delay={index * 0.08} className="h-full">
                <Link
                  href={track.href}
                  className="group relative flex h-full flex-col overflow-hidden rounded-2xl border border-emerald-400/20 bg-gradient-to-br from-emerald-500/10 to-transparent p-8 transition-colors hover:border-gold-400/60"
                >
                  <span className="font-mono text-xs text-emerald-300">[{track.code}]</span>
                  <h3 className="mt-6 text-xl font-semibold">{track.title}</h3>
                  <p className="mt-3 flex-1 text-sm leading-relaxed text-forest-100/70">{track.text}</p>
                  <span className="mt-8 inline-flex items-center gap-2 text-sm font-semibold text-gold-400">
                    {t.cta.learnMore}
                    <Icon name="arrowRight" className="h-4 w-4 transition-transform group-hover:translate-x-1" />
                  </span>
                </Link>
              </Reveal>
            ))}
          </div>
        </div>
      </section>

      <PoleOffers tone="dark" eyebrow="Thématiques" title="Nos domaines d'innovation" offers={pole.offers} />

      <CTA
        title="Innover sans exposer votre organisation."
        text="Cloud, IA générative, agents autonomes : sécurisons ensemble vos projets les plus ambitieux."
        primaryHref="/contact?besoin=labs#formulaire"
      />
    </>
  );
}
