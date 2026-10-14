/**
 * Configuration centrale du site WOOLYS.
 *
 * Les valeurs marquées PLACEHOLDER doivent être remplacées par les informations
 * officielles de l'entreprise avant la mise en production (via les variables
 * d'environnement décrites dans `.env.example`).
 */

const clean = (value: string | undefined) => {
  const trimmed = value?.trim();
  return trimmed ? trimmed : null;
};

export const siteConfig = {
  name: "WOOLYS",
  legalName: "WOOLYS", // PLACEHOLDER — raison sociale officielle
  tagline: "Cybersecurity & Trust",
  signature: "La confiance au cœur du numérique.",
  positioning:
    "Cabinet de conseil et de services numériques spécialisé dans la cybersécurité, la gestion des risques et la confiance numérique.",
  description:
    "WOOLYS accompagne les organisations du Sénégal, d'Afrique et de l'international en cybersécurité, gouvernance, gestion des risques, GRC, conformité, Cloud Security, AI Security et formation.",
  url: (clean(process.env.NEXT_PUBLIC_SITE_URL) ?? "https://www.woolys.example").replace(/\/$/, ""), // PLACEHOLDER — domaine définitif
  locale: "fr_FR",
  /** Plateforme de formation (LMS WOOLYS Academy). `null` tant qu'elle n'est pas configurée. */
  lmsUrl: clean(process.env.NEXT_PUBLIC_LMS_URL),
  /** E-mail de contact. `null` tant qu'il n'est pas configuré. */
  contactEmail: clean(process.env.NEXT_PUBLIC_CONTACT_EMAIL),
  /** Lien de prise de rendez-vous externe. `null` => renvoie vers le formulaire. */
  bookingUrl: clean(process.env.NEXT_PUBLIC_BOOKING_URL),
  linkedinUrl: clean(process.env.NEXT_PUBLIC_LINKEDIN_URL),
  areaServed: ["Sénégal", "Afrique de l'Ouest", "Afrique francophone", "Afrique", "Europe", "International"],
} as const;

export const bookingHref = siteConfig.bookingUrl ?? "/contact?besoin=echange#formulaire";

export const absoluteUrl = (path = "/") => `${siteConfig.url}${path.startsWith("/") ? path : `/${path}`}`;
