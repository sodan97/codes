import type { IconName } from "@/components/Icon";

export type FaqItem = { question: string; answer: string };

export type TitledText = { title: string; text: string };

export type Expertise = {
  slug: string;
  number: string;
  title: string;
  icon: IconName;
  /** Pôle WOOLYS principalement concerné. */
  pole: PoleSlug;
  summary: string;
  hero: { title: string; subtitle: string };
  problem: { title: string; intro: string; points: string[] };
  approach: TitledText[];
  services: TitledText[];
  /** Slugs de `methodologies.ts` (ou libellés libres pour les référentiels sans page dédiée). */
  methodologies: string[];
  deliverables: string[];
  faq: FaqItem[];
  seo: { title: string; description: string; keywords: string[] };
};

export type PoleSlug = "advisory" | "academy" | "labs";

export type Pole = {
  slug: PoleSlug;
  number: string;
  name: string;
  tagline: string;
  summary: string;
  href: string;
  cta: string;
  icon: IconName;
  offers: TitledText[];
};

export type Sector = {
  slug: string;
  name: string;
  icon: IconName;
  summary: string;
  hero: string;
  challenges: string[];
  frameworks: string[];
  support: TitledText[];
  expertises: string[];
};

export type Methodology = {
  slug: string;
  name: string;
  fullName: string;
  category: string;
  summary: string;
  purpose: string;
  whenToUse: string[];
  howWeHelp: string[];
  expertises: string[];
};

export type ArticleBlock =
  | { type: "p"; text: string }
  | { type: "h2"; text: string }
  | { type: "h3"; text: string }
  | { type: "ul"; items: string[] };

export type Article = {
  slug: string;
  title: string;
  excerpt: string;
  category: string;
  /** Date ISO de publication. */
  date: string;
  readingMinutes: number;
  author: string;
  keywords: string[];
  content: ArticleBlock[];
};
