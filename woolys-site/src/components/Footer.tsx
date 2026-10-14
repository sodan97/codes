import Link from "next/link";
import { expertises } from "@/data/expertises";
import { poles } from "@/data/poles";
import { t } from "@/lib/i18n";
import { siteConfig } from "@/lib/site";
import { Icon } from "./Icon";
import { Logo } from "./Logo";

const companyLinks = [
  { label: "À propos", href: "/a-propos" },
  { label: "Secteurs", href: "/secteurs" },
  { label: "Méthodologies", href: "/methodologies" },
  { label: "Insights", href: "/insights" },
  { label: "Contact", href: "/contact" },
];

function FooterColumn({ title, links }: { title: string; links: { label: string; href: string; external?: boolean }[] }) {
  return (
    <div>
      <h2 className="font-display text-xs font-semibold uppercase tracking-[0.2em] text-gold-400">{title}</h2>
      <ul className="mt-5 space-y-3">
        {links.map((link) => (
          <li key={link.href}>
            {link.external ? (
              <a href={link.href} target="_blank" rel="noopener noreferrer" className="text-sm text-forest-100/70 transition-colors hover:text-white">
                {link.label}
                <span className="sr-only"> (s&apos;ouvre dans un nouvel onglet)</span>
              </a>
            ) : (
              <Link href={link.href} className="text-sm text-forest-100/70 transition-colors hover:text-white">
                {link.label}
              </Link>
            )}
          </li>
        ))}
      </ul>
    </div>
  );
}

export function Footer() {
  const year = new Date().getFullYear();
  const poleLinks = [
    ...poles.map((pole) => ({ label: pole.name, href: pole.href })),
    ...(siteConfig.lmsUrl ? [{ label: "Plateforme e-learning", href: siteConfig.lmsUrl, external: true }] : []),
  ];

  return (
    <footer className="relative overflow-hidden bg-forest-950 text-white">
      <div className="h-px bg-gradient-to-r from-transparent via-gold-400/60 to-transparent" />
      <div className="container-site grid gap-12 py-16 sm:py-20 lg:grid-cols-[1.3fr_2fr]">
        <div className="max-w-sm">
          <Logo tone="onDark" />
          <p className="mt-6 text-sm leading-relaxed text-forest-100/70">
            {siteConfig.positioning}
          </p>
          <p className="mt-4 font-display text-sm italic text-gold-300">{siteConfig.signature}</p>
          <div className="mt-8 flex flex-wrap gap-3">
            {siteConfig.contactEmail ? (
              <a
                href={`mailto:${siteConfig.contactEmail}`}
                className="inline-flex items-center gap-2 rounded-full border border-white/15 px-4 py-2 text-sm text-forest-100/80 transition-colors hover:border-gold-400 hover:text-white"
              >
                <Icon name="mail" className="h-4 w-4" />
                {siteConfig.contactEmail}
              </a>
            ) : null}
            {siteConfig.linkedinUrl ? (
              <a
                href={siteConfig.linkedinUrl}
                target="_blank"
                rel="noopener noreferrer"
                aria-label="WOOLYS sur LinkedIn (s'ouvre dans un nouvel onglet)"
                className="inline-flex h-10 w-10 items-center justify-center rounded-full border border-white/15 text-forest-100/80 transition-colors hover:border-gold-400 hover:text-white"
              >
                <Icon name="linkedin" className="h-4 w-4" />
              </a>
            ) : null}
          </div>
        </div>
        <div className="grid grid-cols-2 gap-10 sm:grid-cols-3">
          <FooterColumn title={t.footer.poles} links={poleLinks} />
          <FooterColumn
            title={t.footer.expertises}
            links={expertises.map((item) => ({ label: item.title, href: `/expertises/${item.slug}` }))}
          />
          <FooterColumn title={t.footer.company} links={companyLinks} />
        </div>
      </div>
      <div className="border-t border-white/10">
        <div className="container-site flex flex-col gap-4 py-6 text-xs text-forest-100/60 sm:flex-row sm:items-center sm:justify-between">
          <p>
            © {year} {siteConfig.legalName}. {t.footer.rights} · Sénégal · Afrique · International
          </p>
          <ul className="flex gap-6">
            <li>
              <Link href="/mentions-legales" className="hover:text-white">
                Mentions légales
              </Link>
            </li>
            <li>
              <Link href="/confidentialite" className="hover:text-white">
                Confidentialité
              </Link>
            </li>
          </ul>
        </div>
      </div>
    </footer>
  );
}
