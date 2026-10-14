import type { Methodology } from "./types";

/**
 * Référentiels et méthodes utilisés ou pris en compte dans les missions WOOLYS.
 * IMPORTANT : ils ne sont pas présentés comme des certifications détenues par WOOLYS.
 */
export const methodologies: Methodology[] = [
  {
    slug: "ebios-rm",
    name: "EBIOS RM",
    fullName: "Expression des Besoins et Identification des Objectifs de Sécurité — Risk Manager",
    category: "Analyse de risques",
    summary: "La méthode d'analyse de risques cyber publiée par l'ANSSI, structurée en cinq ateliers.",
    purpose:
      "EBIOS Risk Manager permet d'identifier et de hiérarchiser les scénarios de risque les plus réalistes, en partant des valeurs métier de l'organisation et des sources de risque qui pourraient les viser. Elle aboutit à une stratégie de traitement des risques compréhensible par la direction.",
    whenToUse: [
      "Avant un projet numérique sensible ou une transformation majeure",
      "Pour construire ou réviser la cartographie des risques cyber",
      "Dans le cadre d'une démarche d'homologation ou de conformité",
      "Pour justifier des choix de sécurité auprès d'une direction ou d'un régulateur",
    ],
    howWeHelp: [
      "Animation des cinq ateliers : cadrage et socle, sources de risque, scénarios stratégiques, scénarios opérationnels, traitement",
      "Mobilisation ciblée des métiers pour limiter la charge des équipes",
      "Production de livrables clairs et réutilisables",
      "Transfert de compétences à vos équipes via WOOLYS Academy",
    ],
    expertises: ["risk-management", "gouvernance-cyber"],
  },
  {
    slug: "iso-27001",
    name: "ISO 27001",
    fullName: "ISO/IEC 27001 — Systèmes de management de la sécurité de l'information",
    category: "Norme internationale",
    summary: "La norme de référence pour mettre en place un système de management de la sécurité de l'information (SMSI).",
    purpose:
      "ISO/IEC 27001 définit les exigences d'un système de management de la sécurité de l'information : une démarche d'amélioration continue fondée sur l'analyse des risques, des mesures de sécurité adaptées et un pilotage par la direction. Elle peut donner lieu à une certification délivrée par un organisme accrédité.",
    whenToUse: [
      "Pour structurer durablement la sécurité de l'information",
      "Lorsqu'un client, un partenaire ou un appel d'offres exige une certification",
      "Pour démontrer la maturité de votre organisation à l'international",
      "Comme socle commun pour répondre à d'autres exigences",
    ],
    howWeHelp: [
      "Analyse d'écarts par rapport à la norme",
      "Définition du périmètre et de la déclaration d'applicabilité",
      "Rédaction de la documentation et mise en place des processus",
      "Audit interne et préparation à l'audit de certification",
    ],
    expertises: ["grc-compliance", "gouvernance-cyber"],
  },
  {
    slug: "iso-27005",
    name: "ISO 27005",
    fullName: "ISO/IEC 27005 — Gestion des risques liés à la sécurité de l'information",
    category: "Norme internationale",
    summary: "Les lignes directrices pour gérer les risques de sécurité de l'information, en cohérence avec ISO 27001.",
    purpose:
      "ISO/IEC 27005 fournit un cadre pour identifier, analyser, évaluer et traiter les risques liés à la sécurité de l'information. Elle soutient l'exigence d'appréciation des risques d'ISO 27001 et s'articule avec des méthodes opérationnelles comme EBIOS RM.",
    whenToUse: [
      "Pour formaliser un processus de gestion des risques dans un SMSI",
      "Pour harmoniser les pratiques d'analyse de risques entre entités",
      "Pour intégrer les risques cyber à la gestion des risques d'entreprise",
    ],
    howWeHelp: [
      "Conception d'une méthodologie de gestion des risques adaptée",
      "Définition des échelles de vraisemblance et de gravité",
      "Conduite des appréciations de risques",
      "Mise en place du suivi et du registre des risques",
    ],
    expertises: ["risk-management", "grc-compliance"],
  },
  {
    slug: "nist-csf",
    name: "NIST CSF",
    fullName: "NIST Cybersecurity Framework 2.0",
    category: "Cadre de référence",
    summary: "Un cadre international pour évaluer et améliorer la posture de cybersécurité, organisé en six fonctions.",
    purpose:
      "Le NIST Cybersecurity Framework structure la cybersécurité autour de six fonctions : Gouverner, Identifier, Protéger, Détecter, Répondre et Rétablir. Il permet de mesurer une maturité, de définir un profil cible et de dialoguer avec la direction dans un langage commun.",
    whenToUse: [
      "Pour réaliser un diagnostic de maturité global",
      "Pour définir une trajectoire d'amélioration entre profil actuel et profil cible",
      "Pour comparer plusieurs entités d'un même groupe",
      "Comme langage commun avec des partenaires internationaux",
    ],
    howWeHelp: [
      "Évaluation de la maturité selon les six fonctions",
      "Construction des profils actuel et cible",
      "Feuille de route priorisée",
      "Reporting de maturité pour la direction",
    ],
    expertises: ["gouvernance-cyber", "risk-management", "cloud-security"],
  },
  {
    slug: "cis-controls",
    name: "CIS Controls",
    fullName: "CIS Critical Security Controls v8",
    category: "Mesures de sécurité",
    summary: "Un ensemble priorisé de mesures de sécurité concrètes, organisé en groupes de mise en œuvre.",
    purpose:
      "Les CIS Controls proposent dix-huit familles de mesures de sécurité concrètes et priorisées, réparties en trois groupes de mise en œuvre (IG1 à IG3) selon la maturité et les moyens de l'organisation. Ils constituent une base pragmatique pour réduire rapidement l'exposition aux attaques courantes.",
    whenToUse: [
      "Pour définir une hygiène de sécurité de base",
      "Lorsque les ressources sont limitées et qu'il faut prioriser",
      "Pour évaluer concrètement le niveau de sécurité technique",
    ],
    howWeHelp: [
      "Évaluation selon le groupe de mise en œuvre adapté",
      "Identification des mesures prioritaires",
      "Plan de remédiation réaliste",
      "Suivi de la progression",
    ],
    expertises: ["gouvernance-cyber", "cloud-security"],
  },
  {
    slug: "dora",
    name: "DORA",
    fullName: "Digital Operational Resilience Act — Règlement (UE) 2022/2554",
    category: "Réglementation européenne",
    summary: "Le règlement européen sur la résilience opérationnelle numérique du secteur financier.",
    purpose:
      "DORA impose aux entités financières de l'Union européenne un cadre de gestion des risques liés aux technologies de l'information et de la communication : gouvernance, gestion des incidents, tests de résilience et maîtrise des prestataires TIC. Il est applicable depuis janvier 2025.",
    whenToUse: [
      "Si votre organisation est une entité financière européenne ou une filiale de groupe européen",
      "Si vous êtes prestataire TIC d'acteurs financiers européens",
      "Comme référence de bonnes pratiques de résilience opérationnelle",
    ],
    howWeHelp: [
      "Analyse d'applicabilité et analyse d'écarts",
      "Cadre de gestion des risques TIC",
      "Processus de gestion et de notification des incidents",
      "Encadrement des prestataires TIC",
    ],
    expertises: ["grc-compliance", "risk-management"],
  },
  {
    slug: "pci-dss",
    name: "PCI-DSS",
    fullName: "Payment Card Industry Data Security Standard v4.0",
    category: "Standard sectoriel",
    summary: "Le standard de sécurité des données de cartes de paiement.",
    purpose:
      "PCI-DSS définit les exigences de sécurité applicables à toute organisation qui stocke, traite ou transmet des données de cartes de paiement. Sa version 4.0 renforce notamment l'authentification, la surveillance et l'approche fondée sur les risques.",
    whenToUse: [
      "Si vous traitez des paiements par carte (banques, commerçants, monétique)",
      "Lorsqu'un réseau de paiement ou un acquéreur exige une conformité",
      "Pour réduire le périmètre de données sensibles",
    ],
    howWeHelp: [
      "Définition et réduction du périmètre de conformité",
      "Analyse d'écarts selon la version 4.0",
      "Accompagnement à la remédiation",
      "Préparation à l'évaluation de conformité",
    ],
    expertises: ["grc-compliance"],
  },
  {
    slug: "aws-security",
    name: "AWS Security",
    fullName: "AWS Well-Architected Framework — Pilier Sécurité",
    category: "Sécurité Cloud",
    summary: "Les bonnes pratiques de sécurité des environnements Amazon Web Services.",
    purpose:
      "Le pilier sécurité du AWS Well-Architected Framework, associé au modèle de responsabilité partagée, définit les bonnes pratiques pour protéger les identités, les données, les réseaux et détecter les incidents dans un environnement AWS.",
    whenToUse: [
      "Avant ou pendant une migration vers AWS",
      "Pour évaluer la sécurité d'un environnement AWS existant",
      "Pour industrialiser la sécurité via l'automatisation",
    ],
    howWeHelp: [
      "Revue de sécurité de vos comptes et architectures",
      "Conception d'une landing zone sécurisée",
      "Mise en place de la journalisation et de la détection",
      "Automatisation des contrôles de conformité",
    ],
    expertises: ["cloud-security"],
  },
];

export const getMethodology = (slug: string) => methodologies.find((item) => item.slug === slug);

/** Retourne le libellé d'un référentiel, qu'il dispose d'une page dédiée ou non. */
export const methodologyLabel = (slugOrLabel: string) => getMethodology(slugOrLabel)?.name ?? slugOrLabel;
