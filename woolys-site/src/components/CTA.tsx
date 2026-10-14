import { ButtonLink } from "./Button";
import { Emblem } from "./Logo";
import { Reveal } from "./Reveal";
import { bookingHref } from "@/lib/site";
import { t } from "@/lib/i18n";

type CTAProps = {
  title?: string;
  text?: string;
  primaryLabel?: string;
  primaryHref?: string;
  secondaryLabel?: string;
  secondaryHref?: string;
};

/** Bloc d'appel à l'action final, fond sombre premium. */
export function CTA({
  title = "Construisons ensemble votre confiance numérique.",
  text = "Un projet de cybersécurité, de conformité, de gestion des risques ou de transformation numérique ?",
  primaryLabel = t.cta.expert,
  primaryHref = "/contact",
  secondaryLabel = t.cta.meeting,
  secondaryHref = bookingHref,
}: CTAProps) {
  return (
    <section className="relative isolate overflow-hidden bg-forest-950 py-24 text-white sm:py-32" aria-labelledby="cta-title">
      <div className="absolute inset-0 -z-10 bg-grid opacity-50" />
      <div className="absolute inset-0 -z-10 bg-[radial-gradient(circle_at_50%_120%,rgba(21,143,122,0.35),transparent_60%)]" />
      <Emblem tone="mono" className="absolute -bottom-40 left-1/2 -z-10 h-[34rem] w-[34rem] -translate-x-1/2 text-white/[0.03]" />
      <Reveal className="container-site text-center">
        <p className="eyebrow mb-6 justify-center text-gold-400">WOOLYS · Cybersecurity &amp; Trust</p>
        <h2 id="cta-title" className="mx-auto max-w-3xl text-3xl font-bold leading-tight sm:text-5xl">
          {title}
        </h2>
        <p className="mx-auto mt-6 max-w-2xl text-lg text-forest-100/80">{text}</p>
        <div className="mt-10 flex flex-col items-center justify-center gap-4 sm:flex-row">
          <ButtonLink href={primaryHref} variant="gold" arrow>
            {primaryLabel}
          </ButtonLink>
          <ButtonLink href={secondaryHref} variant="outline">
            {secondaryLabel}
          </ButtonLink>
        </div>
      </Reveal>
    </section>
  );
}
