import type { Pole } from "./types";

export const poles: Pole[] = [
  {
    slug: "advisory",
    number: "01",
    name: "WOOLYS Advisory",
    tagline: "Conseil, gouvernance et gestion des risques cyber.",
    summary:
      "Nous accompagnons les directions dans la structuration de leur gouvernance cyber, l'analyse de leurs risques et leur mise en conformité.",
    href: "/advisory",
    cta: "Découvrir Advisory",
    icon: "governance",
    offers: [
      { title: "Cybersecurity Advisory", text: "Conseil stratégique et opérationnel en cybersécurité." },
      { title: "GRC", text: "Gouvernance, risques et conformité réunis dans un dispositif cohérent." },
      { title: "Risk Management", text: "Identification, évaluation et traitement des risques cyber." },
      { title: "EBIOS RM", text: "Analyse de risques selon la méthode de l'ANSSI." },
      { title: "ISO 27001 / ISO 27005", text: "Système de management et gestion des risques de la sécurité de l'information." },
      { title: "Compliance", text: "Protection des données, régulateurs sectoriels et exigences contractuelles." },
      { title: "DORA", text: "Résilience opérationnelle numérique du secteur financier." },
      { title: "PCI-DSS", text: "Sécurité des données de cartes de paiement." },
      { title: "NIST CSF & CIS Controls", text: "Évaluation de maturité et mesures de sécurité priorisées." },
      { title: "Cyber Governance", text: "Organisation, politiques, comités et tableaux de bord." },
      { title: "Accompagnement RSSI / DSSI", text: "Appui ou externalisation de la fonction sécurité." },
      { title: "Cyber Risk Assessment", text: "Évaluation rapide de l'exposition pour éclairer une décision." },
    ],
  },
  {
    slug: "academy",
    number: "02",
    name: "WOOLYS Academy",
    tagline: "Former les équipes. Développer les compétences. Renforcer la culture cyber.",
    summary:
      "Des programmes de formation et de sensibilisation, en présentiel et sur notre plateforme e-learning, pour transmettre durablement les compétences.",
    href: "/academy",
    cta: "Découvrir Academy",
    icon: "graduation",
    offers: [
      { title: "Sensibilisation", text: "Les bons réflexes cyber pour tous les collaborateurs." },
      { title: "Executive Awareness", text: "Sessions dédiées aux dirigeants et aux instances de gouvernance." },
      { title: "Cyber Culture", text: "Programmes dans la durée pour ancrer les comportements." },
      { title: "Formation EBIOS RM", text: "Maîtriser la méthode et conduire ses propres analyses." },
      { title: "Formation GRC & Risk Management", text: "Gouvernance, conformité et gestion des risques en pratique." },
      { title: "Formation ISO 27001", text: "Comprendre et mettre en œuvre un SMSI." },
      { title: "Cloud Security", text: "Les fondamentaux de la sécurité des environnements Cloud." },
      { title: "AI Security", text: "Usages sûrs de l'IA et sécurité des applications LLM." },
      { title: "Ateliers pratiques", text: "Cas concrets, exercices de crise et mises en situation." },
    ],
  },
  {
    slug: "labs",
    number: "03",
    name: "WOOLYS Labs",
    tagline: "Explorer les nouvelles frontières de la cybersécurité, du Cloud et de l'IA.",
    summary:
      "Notre pôle d'innovation : veille, expérimentation et expertise sur les technologies qui transforment la sécurité.",
    href: "/labs",
    cta: "Découvrir Labs",
    icon: "flask",
    offers: [
      { title: "Cloud Security", text: "Architectures Cloud sécurisées dès la conception." },
      { title: "AWS Security", text: "Bonnes pratiques et automatisation de la sécurité sur AWS." },
      { title: "AI Security", text: "Gouvernance et sécurité des systèmes d'intelligence artificielle." },
      { title: "LLM Security", text: "Protection des applications fondées sur les modèles de langage." },
      { title: "RAG Security", text: "Contrôle des accès et cloisonnement des données dans les systèmes RAG." },
      { title: "Agentic AI Security", text: "Garde-fous, permissions et supervision des agents autonomes." },
      { title: "Digital Identity", text: "Identité numérique, authentification et gestion des accès." },
      { title: "Security Automation", text: "Automatisation des contrôles et de la réponse." },
      { title: "Cybersecurity Innovation", text: "Veille et expérimentations sur les menaces et défenses émergentes." },
    ],
  },
];

export const getPole = (slug: string) => poles.find((item) => item.slug === slug);
