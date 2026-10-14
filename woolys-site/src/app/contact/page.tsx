import type { Metadata } from "next";
import { Suspense } from "react";
import { ButtonLink } from "@/components/Button";
import { ContactForm, ContactFormFields } from "@/components/ContactForm";
import { Icon } from "@/components/Icon";
import { PageHero } from "@/components/PageHero";
import { Reveal } from "@/components/Reveal";
import { t } from "@/lib/i18n";
import { pageMetadata } from "@/lib/seo";
import { bookingHref, siteConfig } from "@/lib/site";

export const metadata: Metadata = pageMetadata({
  title: "Contact — Parler à un expert en cybersécurité",
  description:
    "Contactez WOOLYS pour un projet de cybersécurité, d'analyse de risques, d'audit, de conformité, de formation, de Cloud Security ou d'AI Security au Sénégal, en Afrique ou à l'international.",
  path: "/contact",
});

const steps = [
  "Nous analysons votre demande et revenons vers vous.",
  "Un premier échange permet de comprendre votre contexte et vos objectifs.",
  "Nous vous adressons une proposition adaptée à vos priorités.",
];

export default function ContactPage() {
  return (
    <>
      <PageHero
        breadcrumbs={[{ label: "Contact", href: "/contact" }]}
        eyebrow="Contact"
        title="Parlons de vos enjeux."
        subtitle="Un projet de cybersécurité, de conformité, de gestion des risques ou de transformation numérique ? Décrivez-nous votre besoin : un expert vous répondra."
      />
      <section className="bg-sand-50 py-20 sm:py-24">
        <div className="container-site grid gap-12 lg:grid-cols-[1.5fr_1fr] lg:gap-16">
          <Reveal className="rounded-3xl border border-forest-800/10 bg-white p-6 sm:p-10">
            <h2 id="formulaire" className="text-2xl font-bold text-forest-800">
              Votre demande
            </h2>
            <div className="mt-8">
              <Suspense fallback={<ContactFormFields />}>
                <ContactForm />
              </Suspense>
            </div>
          </Reveal>
          <aside className="space-y-6">
            <Reveal delay={0.08} className="rounded-3xl bg-forest-800 p-8 text-white">
              <Icon name="calendar" className="h-8 w-8 text-gold-400" />
              <h2 className="mt-5 text-xl font-bold">Vous préférez échanger directement ?</h2>
              <p className="mt-3 text-sm leading-relaxed text-forest-100/80">
                Réservez un créneau pour un premier échange, sans engagement, avec un expert WOOLYS.
              </p>
              {siteConfig.bookingUrl ? (
                <ButtonLink href={bookingHref} variant="gold" arrow className="mt-6">
                  {t.cta.meeting}
                </ButtonLink>
              ) : (
                <p className="mt-6 text-sm text-gold-300">
                  Indiquez « Planifier un échange » dans votre message : nous vous proposerons un créneau.
                </p>
              )}
              {siteConfig.contactEmail ? (
                <a href={`mailto:${siteConfig.contactEmail}`} className="mt-6 flex items-center gap-2 text-sm text-forest-100/80 hover:text-white">
                  <Icon name="mail" className="h-4 w-4" />
                  {siteConfig.contactEmail}
                </a>
              ) : null}
            </Reveal>
            <Reveal delay={0.16} className="rounded-3xl border border-forest-800/10 bg-white p-8">
              <h2 className="text-lg font-bold text-forest-800">Et ensuite ?</h2>
              <ol className="mt-5 space-y-5">
                {steps.map((step, index) => (
                  <li key={step} className="flex gap-4">
                    <span className="flex h-8 w-8 shrink-0 items-center justify-center rounded-full border border-gold-500/60 font-display text-xs font-bold text-forest-800">
                      {index + 1}
                    </span>
                    <span className="text-sm leading-relaxed text-muted">{step}</span>
                  </li>
                ))}
              </ol>
            </Reveal>
            <Reveal delay={0.24} className="rounded-3xl border border-forest-800/10 bg-white p-8">
              <Icon name="globe" className="h-7 w-7 text-emerald-600" />
              <p className="mt-4 text-sm leading-relaxed text-muted">
                Nous intervenons au <strong className="text-forest-800">Sénégal</strong>, en{" "}
                <strong className="text-forest-800">Afrique de l&apos;Ouest</strong>, en Afrique francophone et à l&apos;international.
              </p>
            </Reveal>
          </aside>
        </div>
      </section>
    </>
  );
}
