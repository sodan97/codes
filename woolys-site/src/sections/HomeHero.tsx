import { ButtonLink } from "@/components/Button";
import { Icon } from "@/components/Icon";
import { t } from "@/lib/i18n";
import { NetworkVisual } from "./NetworkVisual";

const pillars = ["Gouvernance", "Risk Management", "GRC & Conformité", "Cloud & AI Security", "Formation"];

export function HomeHero() {
  return (
    <section className="relative isolate overflow-hidden bg-forest-900 text-white" aria-labelledby="hero-title">
      <div className="absolute inset-0 -z-10 bg-grid opacity-60 [mask-image:radial-gradient(ellipse_at_center,black_40%,transparent_80%)]" />
      <div className="absolute inset-0 -z-10 bg-[radial-gradient(ellipse_at_75%_30%,rgba(21,143,122,0.35),transparent_55%)]" />
      <div className="absolute inset-0 -z-10 bg-[radial-gradient(ellipse_at_10%_90%,rgba(212,175,124,0.14),transparent_45%)]" />

      <div className="container-site grid min-h-[100svh] items-center gap-12 pb-20 pt-32 lg:grid-cols-[1.15fr_1fr] lg:gap-8 lg:pb-24 lg:pt-36">
        <div>
          <p className="eyebrow animate-fade-up text-gold-400">Cybersecurity &amp; Digital Trust Advisory</p>
          <h1
            id="hero-title"
            className="mt-7 animate-fade-up text-[2.5rem] font-bold leading-[1.06] [animation-delay:120ms] sm:text-6xl lg:text-[4.1rem]"
          >
            Maîtrisez vos risques.
            <br />
            <span className="text-gradient-gold">Renforcez votre cybersécurité.</span>
          </h1>
          <p className="mt-7 max-w-xl animate-fade-up text-lg leading-relaxed text-forest-100/80 [animation-delay:240ms] sm:text-xl">
            WOOLYS accompagne les organisations dans la gouvernance, la gestion des risques, la conformité et la cybersécurité pour
            construire une confiance numérique durable.
          </p>
          <div className="mt-10 flex animate-fade-up flex-col gap-4 [animation-delay:360ms] sm:flex-row">
            <ButtonLink href="/contact" variant="primary" arrow>
              {t.cta.expert}
            </ButtonLink>
            <ButtonLink href="/expertises" variant="outline">
              {t.cta.discover}
            </ButtonLink>
          </div>
          <ul className="mt-12 flex animate-fade-up flex-wrap gap-x-5 gap-y-3 [animation-delay:480ms]" aria-label="Domaines d'intervention">
            {pillars.map((pillar) => (
              <li key={pillar} className="flex items-center gap-2 text-sm text-forest-100/70">
                <Icon name="check" className="h-4 w-4 text-gold-400" />
                {pillar}
              </li>
            ))}
          </ul>
        </div>
        <div className="animate-fade-up [animation-delay:300ms]">
          <NetworkVisual />
        </div>
      </div>

      <a
        href="#enjeux"
        className="absolute bottom-8 left-1/2 hidden -translate-x-1/2 flex-col items-center gap-2 text-[0.65rem] uppercase tracking-[0.3em] text-forest-100/50 transition-colors hover:text-gold-300 lg:flex"
      >
        Découvrir
        <span className="h-10 w-px bg-gradient-to-b from-gold-400 to-transparent" />
      </a>
      <div className="absolute inset-x-0 bottom-0 h-px bg-gradient-to-r from-transparent via-gold-400/50 to-transparent" />
    </section>
  );
}
