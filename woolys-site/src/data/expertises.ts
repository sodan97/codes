import type { Expertise } from "./types";

export const expertises: Expertise[] = [
  {
    slug: "gouvernance-cyber",
    number: "01",
    title: "Cybersecurity Governance",
    icon: "governance",
    pole: "advisory",
    summary:
      "Structurer l'organisation, les rôles, les politiques et le pilotage de la cybersécurité au niveau de la direction.",
    hero: {
      title: "Gouvernance cyber : faire de la sécurité une décision de direction.",
      subtitle:
        "Nous aidons les directions générales, DSI et RSSI à définir une organisation de sécurité claire, des priorités assumées et un pilotage mesurable.",
    },
    problem: {
      title: "Une cybersécurité sans gouvernance reste une addition de projets.",
      intro:
        "Dans de nombreuses organisations, la sécurité progresse au rythme des incidents et des audits. Les responsabilités sont diffuses, les budgets difficiles à justifier et la direction manque d'une vision consolidée de l'exposition réelle.",
      points: [
        "Rôles et responsabilités cyber mal définis entre métiers, DSI et direction",
        "Absence de politique de sécurité à jour ou appliquée",
        "Difficulté à prioriser les investissements et à les justifier",
        "Reporting au comité de direction trop technique ou inexistant",
        "Fonction RSSI isolée ou à créer",
      ],
    },
    approach: [
      { title: "Diagnostic de maturité", text: "Évaluation de l'organisation, des processus et des mesures existantes au regard d'un référentiel reconnu (ISO 27001, NIST CSF)." },
      { title: "Cible et feuille de route", text: "Définition d'un modèle de gouvernance adapté à votre taille, votre secteur et vos obligations, décliné en chantiers priorisés." },
      { title: "Mise en œuvre", text: "Rédaction des politiques, mise en place des comités, des indicateurs et des rituels de pilotage." },
      { title: "Pilotage dans la durée", text: "Accompagnement du RSSI, tableaux de bord de direction et revues périodiques." },
    ],
    services: [
      { title: "Schéma directeur sécurité", text: "Une vision pluriannuelle alignée sur la stratégie de l'organisation." },
      { title: "Politique de sécurité (PSSI)", text: "Politiques, chartes et procédures rédigées pour être appliquées." },
      { title: "Organisation & comitologie", text: "Rôles, RACI, comités sécurité et circuits de décision." },
      { title: "RSSI externalisé / DSSI", text: "Un accompagnement opérationnel ou en appui de votre RSSI." },
      { title: "Tableaux de bord de direction", text: "Indicateurs lisibles pour la direction générale et le conseil." },
      { title: "Gestion des tiers", text: "Encadrement de la sécurité de vos fournisseurs et prestataires." },
    ],
    methodologies: ["iso-27001", "nist-csf", "cis-controls", "iso-27005"],
    deliverables: [
      "Rapport de maturité et cartographie des écarts",
      "Modèle de gouvernance cible et matrice RACI",
      "Politique de sécurité et documents associés",
      "Feuille de route priorisée et chiffrée",
      "Tableau de bord sécurité pour la direction",
    ],
    faq: [
      { question: "Faut-il obligatoirement un RSSI à temps plein ?", answer: "Pas nécessairement. Selon votre taille et votre exposition, une fonction RSSI partagée ou externalisée peut suffire, à condition que les responsabilités et le temps alloué soient clairement définis." },
      { question: "Par où commencer lorsqu'il n'existe rien de formalisé ?", answer: "Par un diagnostic de maturité court, qui permet d'identifier les priorités et de produire rapidement une première feuille de route et une politique de sécurité socle." },
      { question: "Comment présenter la cybersécurité à un comité de direction ?", answer: "En parlant de risques métier, d'impacts et de décisions à prendre plutôt que de technologies. Nous construisons avec vous un reporting adapté à ce niveau de lecture." },
    ],
    seo: {
      title: "Gouvernance cybersécurité — PSSI, RSSI, schéma directeur",
      description:
        "Gouvernance de la cybersécurité au Sénégal et en Afrique : diagnostic de maturité, PSSI, organisation, RSSI externalisé et tableaux de bord de direction.",
      keywords: ["gouvernance cybersécurité Sénégal", "RSSI externalisé Sénégal", "PSSI", "consultant cybersécurité Sénégal"],
    },
  },
  {
    slug: "risk-management",
    number: "02",
    title: "Risk Management",
    icon: "radar",
    pole: "advisory",
    summary:
      "Identifier, analyser et traiter les risques cyber avec des méthodes reconnues comme EBIOS RM et ISO 27005.",
    hero: {
      title: "Analyse de risques cyber : décider en connaissance de cause.",
      subtitle:
        "Nous cartographions vos risques numériques avec des méthodes éprouvées pour orienter vos investissements vers ce qui protège réellement votre activité.",
    },
    problem: {
      title: "On ne protège bien que ce que l'on a compris.",
      intro:
        "Sans analyse de risques, les mesures de sécurité sont choisies par habitude ou par imitation. Les scénarios les plus critiques pour l'activité restent parfois invisibles, alors que des moyens sont consacrés à des menaces secondaires.",
      points: [
        "Pas de vision consolidée des risques cyber et de leurs impacts métier",
        "Cartographie des risques d'entreprise qui ignore le numérique",
        "Projets lancés sans analyse de risques préalable",
        "Besoin de justifier des choix de sécurité auprès d'un régulateur ou d'un auditeur",
      ],
    },
    approach: [
      { title: "Cadrage", text: "Définition du périmètre, des valeurs métier, des biens supports et du socle de sécurité existant." },
      { title: "Analyse", text: "Identification des sources de risque, construction des scénarios stratégiques et opérationnels." },
      { title: "Évaluation", text: "Estimation de la vraisemblance et de la gravité, hiérarchisation des risques." },
      { title: "Traitement", text: "Plan de traitement, risques résiduels et suivi dans le temps." },
    ],
    services: [
      { title: "Analyse de risques EBIOS RM", text: "Une démarche structurée en cinq ateliers, centrée sur les scénarios réalistes." },
      { title: "Appréciation des risques ISO 27005", text: "Une gestion des risques intégrée à votre système de management." },
      { title: "Cyber Risk Assessment", text: "Une évaluation rapide de l'exposition pour éclairer une décision." },
      { title: "Intégration à l'ERM", text: "Raccordement des risques cyber à la cartographie des risques de l'entreprise." },
      { title: "Analyse de risques projet", text: "Sécurité intégrée dès la conception des projets numériques." },
      { title: "Risques tiers", text: "Évaluation des fournisseurs et partenaires critiques." },
    ],
    methodologies: ["ebios-rm", "iso-27005", "iso-27001", "nist-csf"],
    deliverables: [
      "Cartographie des risques cyber",
      "Scénarios stratégiques et opérationnels documentés",
      "Plan de traitement des risques",
      "Registre des risques résiduels",
      "Synthèse décisionnelle pour la direction",
    ],
    faq: [
      { question: "Quelle différence entre EBIOS RM et ISO 27005 ?", answer: "ISO 27005 fixe les principes et le processus de gestion des risques de sécurité de l'information. EBIOS RM est une méthode opérationnelle, publiée par l'ANSSI, qui permet de conduire concrètement une analyse compatible avec ces principes." },
      { question: "Combien de temps dure une analyse de risques ?", answer: "Cela dépend du périmètre. Une analyse ciblée sur un système ou un projet peut être menée en quelques semaines ; une analyse à l'échelle d'une organisation demande davantage d'ateliers." },
      { question: "Nos équipes doivent-elles être mobilisées ?", answer: "Oui, la qualité de l'analyse repose sur la connaissance des métiers. Nous organisons des ateliers courts et ciblés pour limiter la charge des participants." },
    ],
    seo: {
      title: "Analyse de risques cyber — EBIOS RM et ISO 27005",
      description:
        "Analyse de risques cyber au Sénégal et en Afrique avec EBIOS RM et ISO 27005 : cartographie des risques, scénarios, plan de traitement.",
      keywords: ["analyse risques cyber Sénégal", "EBIOS RM Sénégal", "ISO 27005", "risk management cybersécurité"],
    },
  },
  {
    slug: "grc-compliance",
    number: "03",
    title: "GRC & Compliance",
    icon: "compliance",
    pole: "advisory",
    summary:
      "Répondre aux exigences réglementaires et normatives tout en construisant un dispositif de conformité durable.",
    hero: {
      title: "GRC & conformité : transformer l'obligation en avantage.",
      subtitle:
        "Nous vous aidons à répondre aux exigences ISO 27001, PCI-DSS, DORA ou aux réglementations locales, avec un dispositif proportionné et maintenable.",
    },
    problem: {
      title: "La conformité ne doit pas être un exercice ponctuel.",
      intro:
        "Les exigences se multiplient : protection des données personnelles, régulateurs sectoriels, normes internationales, exigences contractuelles de partenaires européens. Traitées une par une, elles génèrent des efforts redondants et une conformité fragile.",
      points: [
        "Audit ou certification à préparer dans un délai contraint",
        "Exigences réglementaires multiples et difficiles à consolider",
        "Preuves de conformité dispersées ou inexistantes",
        "Partenaires ou maisons mères qui imposent des standards internationaux",
      ],
    },
    approach: [
      { title: "Analyse d'écarts", text: "Évaluation de votre situation face aux exigences applicables." },
      { title: "Plan de mise en conformité", text: "Actions priorisées, responsables et échéances réalistes." },
      { title: "Accompagnement", text: "Rédaction documentaire, mise en place des processus et collecte des preuves." },
      { title: "Préparation à l'audit", text: "Audit à blanc, revue des preuves et accompagnement pendant l'audit." },
    ],
    services: [
      { title: "Accompagnement ISO 27001", text: "De l'analyse d'écarts à la préparation de la certification." },
      { title: "Conformité PCI-DSS", text: "Réduction du périmètre et mise en conformité des environnements de paiement." },
      { title: "Préparation DORA", text: "Résilience opérationnelle numérique pour les acteurs financiers concernés." },
      { title: "Protection des données", text: "Mise en conformité avec les lois de protection des données personnelles." },
      { title: "Cadre de contrôle unifié", text: "Mutualisation des exigences de plusieurs référentiels." },
      { title: "Audit de conformité", text: "Évaluation indépendante de votre dispositif." },
    ],
    methodologies: ["iso-27001", "pci-dss", "dora", "nist-csf", "cis-controls"],
    deliverables: [
      "Rapport d'analyse d'écarts",
      "Plan de mise en conformité priorisé",
      "Déclaration d'applicabilité (SoA) et documentation du SMSI",
      "Matrice de conformité multi-référentiels",
      "Rapport d'audit à blanc",
    ],
    faq: [
      { question: "WOOLYS délivre-t-il des certifications ?", answer: "Non. Nous accompagnons votre organisation dans la préparation. La certification est délivrée par un organisme de certification accrédité, indépendant du conseil." },
      { question: "Peut-on traiter plusieurs référentiels en même temps ?", answer: "Oui. Un cadre de contrôle unifié permet de mutualiser les exigences communes et d'éviter les doublons documentaires." },
      { question: "DORA concerne-t-il les organisations africaines ?", answer: "DORA s'applique aux entités financières de l'Union européenne et à leurs prestataires TIC critiques. Des organisations africaines peuvent donc être concernées indirectement, via leurs relations avec des groupes ou partenaires européens." },
    ],
    seo: {
      title: "GRC et conformité — ISO 27001, PCI-DSS, DORA",
      description:
        "GRC et conformité cybersécurité au Sénégal et en Afrique : ISO 27001, PCI-DSS, DORA, protection des données, audit et préparation à la certification.",
      keywords: ["GRC Sénégal", "ISO 27001 Sénégal", "audit cybersécurité Sénégal", "conformité PCI-DSS Afrique"],
    },
  },
  {
    slug: "cloud-security",
    number: "04",
    title: "Cloud Security",
    icon: "cloud",
    pole: "labs",
    summary:
      "Sécuriser les migrations et les environnements Cloud, avec une expertise particulière sur AWS.",
    hero: {
      title: "Cloud Security : migrer sans transférer vos risques.",
      subtitle:
        "Nous sécurisons vos environnements Cloud dès la conception : architecture, identités, configuration, journalisation et gouvernance.",
    },
    problem: {
      title: "Le Cloud change les règles de la sécurité.",
      intro:
        "Le modèle de responsabilité partagée, la rapidité de déploiement et la multiplication des services exposent à de nouvelles erreurs, souvent liées à la configuration et aux identités plutôt qu'aux attaques sophistiquées.",
      points: [
        "Migration Cloud lancée sans cadre de sécurité défini",
        "Comptes et droits d'accès difficiles à maîtriser",
        "Ressources exposées publiquement par erreur",
        "Manque de visibilité et de journalisation centralisée",
        "Exigences de localisation ou de souveraineté des données",
      ],
    },
    approach: [
      { title: "Évaluer", text: "Revue de l'architecture et de la configuration existantes au regard des bonnes pratiques." },
      { title: "Concevoir", text: "Landing zone sécurisée, modèle d'identités, segmentation et chiffrement." },
      { title: "Automatiser", text: "Contrôles de conformité continus et sécurité intégrée aux pipelines." },
      { title: "Gouverner", text: "Politiques Cloud, responsabilités et indicateurs de suivi." },
    ],
    services: [
      { title: "Revue de sécurité AWS", text: "Analyse de la configuration selon le pilier sécurité du Well-Architected Framework." },
      { title: "Landing zone sécurisée", text: "Une fondation multi-comptes saine pour accueillir vos charges." },
      { title: "Identités & accès (IAM)", text: "Moindre privilège, fédération et gestion des accès à privilèges." },
      { title: "Sécurité des migrations", text: "Analyse de risques et exigences de sécurité dès le projet." },
      { title: "Détection & journalisation", text: "Centralisation des logs et des alertes de sécurité." },
      { title: "Gouvernance Cloud", text: "Politiques, garde-fous et conformité continue." },
    ],
    methodologies: ["aws-security", "cis-controls", "nist-csf", "iso-27001"],
    deliverables: [
      "Rapport de revue de sécurité Cloud",
      "Architecture cible et modèle d'identités",
      "Plan de remédiation priorisé",
      "Politiques et garde-fous Cloud",
      "Indicateurs de posture de sécurité",
    ],
    faq: [
      { question: "Intervenez-vous uniquement sur AWS ?", answer: "AWS est un domaine d'expertise privilégié. Les principes de sécurité Cloud (identités, configuration, journalisation, gouvernance) s'appliquent également aux autres fournisseurs, sur lesquels nous pouvons intervenir selon le contexte." },
      { question: "Le Cloud est-il moins sûr qu'une infrastructure interne ?", answer: "Pas intrinsèquement. Il offre des capacités de sécurité avancées, mais la responsabilité de leur configuration revient au client. La plupart des incidents Cloud proviennent d'erreurs de configuration évitables." },
      { question: "Peut-on sécuriser un environnement déjà en production ?", answer: "Oui. Nous commençons par une revue de l'existant puis proposons des remédiations progressives, priorisées selon le risque." },
    ],
    seo: {
      title: "Cloud Security et sécurité AWS",
      description:
        "Sécurité Cloud et AWS en Afrique : revue de configuration, landing zone sécurisée, IAM, sécurité des migrations et gouvernance Cloud.",
      keywords: ["cloud security Africa", "sécurité AWS Sénégal", "sécurité cloud Afrique", "AWS security"],
    },
  },
  {
    slug: "ai-security",
    number: "05",
    title: "AI Security",
    icon: "ai",
    pole: "labs",
    summary:
      "Encadrer et sécuriser les usages de l'intelligence artificielle : LLM, RAG, agents et données.",
    hero: {
      title: "AI Security : adopter l'IA sans perdre le contrôle.",
      subtitle:
        "Nous aidons les organisations à déployer l'intelligence artificielle de manière sûre : gouvernance des usages, sécurité des LLM, des architectures RAG et des agents.",
    },
    problem: {
      title: "L'IA crée de nouvelles surfaces d'attaque.",
      intro:
        "Les assistants conversationnels, les systèmes RAG et les agents autonomes accèdent à des données sensibles et peuvent agir sur des systèmes. Ils introduisent des risques spécifiques : injection de prompt, fuite de données, actions non maîtrisées.",
      points: [
        "Usages de l'IA générative non encadrés par les collaborateurs",
        "Données confidentielles exposées à des services externes",
        "Applications LLM vulnérables à l'injection de prompt",
        "Agents disposant de permissions trop larges",
        "Absence de cadre de gouvernance de l'IA",
      ],
    },
    approach: [
      { title: "Cartographier", text: "Inventaire des cas d'usage, des modèles, des données et des intégrations." },
      { title: "Analyser", text: "Analyse de risques dédiée aux systèmes d'IA et aux menaces spécifiques." },
      { title: "Sécuriser", text: "Mesures de protection des prompts, des données, des outils et des sorties." },
      { title: "Gouverner", text: "Politique d'usage, rôles, validation des cas d'usage et suivi." },
    ],
    services: [
      { title: "Gouvernance de l'IA", text: "Politique d'usage, comité IA et processus de validation." },
      { title: "Sécurité des LLM", text: "Évaluation des applications face aux menaces spécifiques aux modèles de langage." },
      { title: "Sécurité RAG", text: "Contrôle des accès aux données, cloisonnement et traçabilité." },
      { title: "Sécurité des agents", text: "Permissions, supervision humaine et garde-fous des agents autonomes." },
      { title: "Analyse de risques IA", text: "Identification des scénarios de risque propres à vos cas d'usage." },
      { title: "Sensibilisation IA", text: "Former les équipes aux bons usages de l'IA générative." },
    ],
    methodologies: ["OWASP Top 10 for LLM Applications", "NIST AI RMF", "ISO/IEC 42001", "ebios-rm"],
    deliverables: [
      "Cartographie des usages et des risques IA",
      "Politique d'usage de l'IA",
      "Rapport d'évaluation de sécurité d'une application LLM",
      "Recommandations d'architecture sécurisée (RAG, agents)",
      "Plan d'actions priorisé",
    ],
    faq: [
      { question: "Faut-il interdire l'IA générative en entreprise ?", answer: "Une interdiction totale est rarement tenable. Il est préférable d'encadrer les usages : outils autorisés, données interdites, règles de vérification et sensibilisation des équipes." },
      { question: "Qu'est-ce qu'une injection de prompt ?", answer: "C'est une technique consistant à insérer des instructions malveillantes dans les données traitées par un modèle de langage pour détourner son comportement, par exemple pour lui faire divulguer des informations ou exécuter une action non prévue." },
      { question: "Sécurisez-vous des solutions développées en interne ?", answer: "Oui. Nous intervenons sur les applications développées en interne comme sur les solutions du marché intégrées à vos systèmes." },
    ],
    seo: {
      title: "AI Security — sécurité des LLM, RAG et agents IA",
      description:
        "Sécurité de l'intelligence artificielle en Afrique : gouvernance de l'IA, sécurité des LLM, des architectures RAG et des agents autonomes.",
      keywords: ["AI security Africa", "sécurité IA Sénégal", "LLM security", "sécurité RAG"],
    },
  },
  {
    slug: "sensibilisation",
    number: "06",
    title: "Cybersecurity Awareness",
    icon: "people",
    pole: "academy",
    summary:
      "Faire de chaque collaborateur un acteur de la sécurité grâce à des programmes de sensibilisation et de formation.",
    hero: {
      title: "Sensibilisation : vos collaborateurs, première ligne de défense.",
      subtitle:
        "Nous concevons des programmes de sensibilisation adaptés à chaque population, des dirigeants aux équipes opérationnelles, pour ancrer durablement une culture cyber.",
    },
    problem: {
      title: "La technologie seule ne suffit pas.",
      intro:
        "Hameçonnage, fraude au président, mots de passe réutilisés, partage imprudent d'informations : une grande partie des incidents commence par une action humaine. Une sensibilisation ponctuelle et générique change peu les comportements.",
      points: [
        "Collaborateurs peu préparés à reconnaître une tentative d'hameçonnage",
        "Dirigeants ciblés par des fraudes et des usurpations d'identité",
        "Formations génériques, peu adaptées aux métiers",
        "Aucune mesure de l'évolution des comportements",
      ],
    },
    approach: [
      { title: "Comprendre", text: "Identification des populations, des risques métiers et du niveau de maturité." },
      { title: "Concevoir", text: "Parcours adaptés à chaque public, en présentiel et en ligne." },
      { title: "Déployer", text: "Sessions, ateliers, exercices et campagnes de communication." },
      { title: "Mesurer", text: "Indicateurs de participation et d'évolution des comportements." },
    ],
    services: [
      { title: "Executive Awareness", text: "Sessions dédiées aux dirigeants et aux comités de direction." },
      { title: "Sensibilisation collaborateurs", text: "Les bons réflexes, illustrés par des situations concrètes." },
      { title: "Exercices d'hameçonnage", text: "Campagnes simulées et accompagnement pédagogique." },
      { title: "Exercices de crise", text: "Simulation d'une cyberattaque pour entraîner les équipes de décision." },
      { title: "Parcours e-learning", text: "Modules en ligne via la plateforme WOOLYS Academy." },
      { title: "Formation des équipes IT", text: "Montée en compétences des équipes techniques et sécurité." },
    ],
    methodologies: ["iso-27001", "nist-csf"],
    deliverables: [
      "Plan de sensibilisation annuel",
      "Supports et modules adaptés à chaque public",
      "Rapport de campagne d'hameçonnage simulé",
      "Retour d'expérience d'exercice de crise",
      "Indicateurs de suivi de la culture cyber",
    ],
    faq: [
      { question: "Les formations sont-elles adaptées à notre secteur ?", answer: "Oui. Les contenus sont contextualisés selon votre secteur, vos métiers et les menaces que rencontrent réellement vos équipes." },
      { question: "Proposez-vous des formations en ligne ?", answer: "Oui, via la plateforme de formation WOOLYS Academy, en complément de sessions en présentiel ou à distance." },
      { question: "Comment mesurer l'efficacité d'un programme ?", answer: "Par des indicateurs combinés : taux de participation, résultats des exercices simulés, signalements d'incidents par les collaborateurs et évolution dans le temps." },
    ],
    seo: {
      title: "Sensibilisation et formation cybersécurité",
      description:
        "Programmes de sensibilisation à la cybersécurité au Sénégal et en Afrique : dirigeants, collaborateurs, exercices d'hameçonnage et de crise.",
      keywords: ["sensibilisation cybersécurité Sénégal", "formation cybersécurité Sénégal", "awareness cybersécurité Afrique"],
    },
  },
];

export const getExpertise = (slug: string) => expertises.find((item) => item.slug === slug);
