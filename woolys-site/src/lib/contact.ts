import { budgetRanges, needTypes } from "@/data/navigation";

export const countries = [
  "Sénégal",
  "Bénin",
  "Burkina Faso",
  "Cameroun",
  "Côte d'Ivoire",
  "Gabon",
  "Gambie",
  "Guinée",
  "Guinée-Bissau",
  "Mali",
  "Maroc",
  "Mauritanie",
  "Niger",
  "République démocratique du Congo",
  "République du Congo",
  "Rwanda",
  "Tchad",
  "Togo",
  "Tunisie",
  "Autre pays d'Afrique",
  "Belgique",
  "Canada",
  "France",
  "Suisse",
  "Autre pays",
] as const;

export type ContactPayload = {
  firstName: string;
  lastName: string;
  company: string;
  jobTitle: string;
  email: string;
  phone: string;
  country: string;
  need: string;
  budget: string;
  message: string;
  consent: boolean;
};

export type ContactErrors = Partial<Record<keyof ContactPayload, string>>;

const EMAIL_RE = /^[^\s@]+@[^\s@]+\.[^\s@]{2,}$/;
const PHONE_RE = /^[+()\d\s.-]{6,24}$/;

const str = (value: unknown, max: number) => (typeof value === "string" ? value.trim().slice(0, max) : "");

/** Validation partagée entre le navigateur et l'API. */
export function validateContact(input: Record<string, unknown>): { value: ContactPayload; errors: ContactErrors } {
  const value: ContactPayload = {
    firstName: str(input.firstName, 80),
    lastName: str(input.lastName, 80),
    company: str(input.company, 120),
    jobTitle: str(input.jobTitle, 120),
    email: str(input.email, 160),
    phone: str(input.phone, 30),
    country: str(input.country, 80),
    need: str(input.need, 60),
    budget: str(input.budget, 60),
    message: str(input.message, 5000),
    consent: input.consent === true || input.consent === "on" || input.consent === "true",
  };

  const errors: ContactErrors = {};
  if (!value.firstName) errors.firstName = "Veuillez indiquer votre prénom.";
  if (!value.lastName) errors.lastName = "Veuillez indiquer votre nom.";
  if (!value.company) errors.company = "Veuillez indiquer votre entreprise ou organisation.";
  if (!value.jobTitle) errors.jobTitle = "Veuillez indiquer votre fonction.";
  if (!EMAIL_RE.test(value.email)) errors.email = "Veuillez saisir une adresse e-mail valide.";
  if (value.phone && !PHONE_RE.test(value.phone)) errors.phone = "Le numéro de téléphone semble invalide.";
  if (!value.country || !(countries as readonly string[]).includes(value.country)) errors.country = "Veuillez sélectionner un pays.";
  if (!(needTypes as readonly string[]).includes(value.need)) errors.need = "Veuillez sélectionner le type de besoin.";
  if (value.budget && !(budgetRanges as readonly string[]).includes(value.budget)) errors.budget = "Valeur de budget invalide.";
  if (value.message.length < 20) errors.message = "Merci de décrire votre besoin en quelques phrases (20 caractères minimum).";
  if (!value.consent) errors.consent = "Votre accord est nécessaire pour que nous puissions vous recontacter.";

  return { value, errors };
}

const slug = (text: string) =>
  text
    .normalize("NFD")
    .replace(/[̀-ͯ]/g, "")
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "-")
    .replace(/(^-|-$)/g, "");

const aliases: Record<string, string> = {
  diagnostic: "Audit",
  proposition: "Conseil",
  echange: "Conseil",
  academy: "Formation",
  sensibilisation: "Formation",
  "risk-management": "Analyse de risques",
  grc: "Conformité",
  "grc-compliance": "Conformité",
  "gouvernance-cyber": "Conseil",
  advisory: "Conseil",
  labs: "Cloud Security",
};

/** Convertit le paramètre d'URL `?besoin=` en type de besoin du formulaire. */
export function needFromParam(param: string | null | undefined): string {
  if (!param) return "";
  const key = slug(param);
  return needTypes.find((need) => slug(need) === key) ?? aliases[key] ?? "";
}
