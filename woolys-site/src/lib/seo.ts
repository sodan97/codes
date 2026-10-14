import type { Metadata } from "next";
import { absoluteUrl, siteConfig } from "./site";

type PageMetaInput = {
  title: string;
  description: string;
  path: string;
  keywords?: string[];
  type?: "website" | "article";
  publishedTime?: string;
};

export const baseKeywords = [
  "cybersécurité Sénégal",
  "cabinet cybersécurité Sénégal",
  "consultant cybersécurité Sénégal",
  "GRC Sénégal",
  "analyse risques cyber Sénégal",
  "EBIOS RM Sénégal",
  "ISO 27001 Sénégal",
  "audit cybersécurité Sénégal",
  "cybersecurity Africa",
  "cybersecurity consulting Africa",
  "AI security Africa",
  "cloud security Africa",
];

/** Métadonnées cohérentes (canonical, Open Graph, Twitter) pour chaque page. */
export function pageMetadata({ title, description, path, keywords = [], type = "website", publishedTime }: PageMetaInput): Metadata {
  const url = absoluteUrl(path);
  return {
    title,
    description,
    keywords: [...keywords, ...baseKeywords.slice(0, 6)],
    alternates: {
      canonical: url,
      // Pour l'anglais : languages: { "fr-FR": url, "en": absoluteUrl(`/en${path}`) }
    },
    openGraph: {
      type,
      url,
      title: `${title} | ${siteConfig.name}`,
      description,
      siteName: siteConfig.name,
      locale: siteConfig.locale,
      ...(publishedTime ? { publishedTime } : {}),
    },
    twitter: {
      card: "summary_large_image",
      title: `${title} | ${siteConfig.name}`,
      description,
    },
  };
}

export const organizationJsonLd = () => ({
  "@context": "https://schema.org",
  "@graph": [
    {
      "@type": "Organization",
      "@id": `${siteConfig.url}/#organization`,
      name: siteConfig.name,
      legalName: siteConfig.legalName,
      url: siteConfig.url,
      logo: absoluteUrl("/brand/woolys-emblem.svg"),
      slogan: siteConfig.signature,
      description: siteConfig.description,
      ...(siteConfig.contactEmail ? { email: siteConfig.contactEmail } : {}),
      ...(siteConfig.linkedinUrl ? { sameAs: [siteConfig.linkedinUrl] } : {}),
    },
    {
      "@type": "ProfessionalService",
      "@id": `${siteConfig.url}/#service`,
      name: `${siteConfig.name} — ${siteConfig.tagline}`,
      url: siteConfig.url,
      image: absoluteUrl("/opengraph-image"),
      description: siteConfig.positioning,
      parentOrganization: { "@id": `${siteConfig.url}/#organization` },
      areaServed: siteConfig.areaServed.map((name) => ({ "@type": "Place", name })),
      knowsAbout: [
        "Cybersécurité",
        "Gouvernance de la sécurité",
        "Gestion des risques",
        "GRC",
        "Conformité",
        "EBIOS RM",
        "ISO 27001",
        "Cloud Security",
        "AI Security",
        "Sensibilisation à la cybersécurité",
      ],
      availableLanguage: ["fr"],
    },
    {
      "@type": "WebSite",
      "@id": `${siteConfig.url}/#website`,
      url: siteConfig.url,
      name: siteConfig.name,
      inLanguage: "fr-FR",
      publisher: { "@id": `${siteConfig.url}/#organization` },
    },
  ],
});
