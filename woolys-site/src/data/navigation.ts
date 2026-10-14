import { expertises } from "./expertises";

export type NavItem = {
  label: string;
  href: string;
  children?: { label: string; href: string; description?: string }[];
};

export const mainNav: NavItem[] = [
  {
    label: "Expertises",
    href: "/expertises",
    children: expertises.map((item) => ({
      label: item.title,
      href: `/expertises/${item.slug}`,
      description: item.summary,
    })),
  },
  { label: "Secteurs", href: "/secteurs" },
  { label: "Méthodologies", href: "/methodologies" },
  { label: "Academy", href: "/academy" },
  { label: "Labs", href: "/labs" },
  { label: "Insights", href: "/insights" },
  { label: "À propos", href: "/a-propos" },
];

export const approachSteps = [
  { number: "01", title: "Comprendre", text: "Comprendre votre organisation, vos enjeux et votre exposition." },
  { number: "02", title: "Évaluer", text: "Identifier vos risques, vulnérabilités et écarts." },
  { number: "03", title: "Transformer", text: "Construire une feuille de route adaptée à vos priorités." },
  { number: "04", title: "Renforcer", text: "Mesurer, améliorer et développer durablement votre résilience." },
];

export const needTypes = [
  "Conseil",
  "Analyse de risques",
  "Audit",
  "Conformité",
  "Formation",
  "Cloud Security",
  "AI Security",
  "Autre",
] as const;

export const budgetRanges = [
  "Non défini",
  "Moins de 5 M FCFA",
  "5 à 15 M FCFA",
  "15 à 50 M FCFA",
  "Plus de 50 M FCFA",
] as const;
