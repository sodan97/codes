import type { Article } from "./types";

/**
 * Articles de démarrage du blog « Insights ».
 * PLACEHOLDER éditorial : contenus à relire et enrichir par l'équipe WOOLYS avant publication.
 * Pour ajouter un article, ajouter un objet à ce tableau : la page, le sitemap et le flux
 * de la homepage sont générés automatiquement.
 */
export const articles: Article[] = [
  {
    slug: "5-priorites-programme-cybersecurite",
    title: "Les 5 priorités d'un programme de cybersécurité",
    excerpt:
      "Gouvernance, connaissance des actifs, gestion des identités, détection et préparation à la crise : les fondations à poser avant toute chose.",
    category: "Gouvernance",
    date: "2026-10-01",
    readingMinutes: 6,
    author: "Équipe WOOLYS",
    keywords: ["programme cybersécurité", "gouvernance cyber", "priorités RSSI"],
    content: [
      { type: "p", text: "Face à la multiplication des menaces et des exigences, beaucoup d'organisations se demandent par où commencer. La tentation est forte d'acheter une nouvelle solution technique. Pourtant, un programme de cybersécurité efficace repose d'abord sur quelques fondations simples, qui conditionnent la réussite de tout le reste." },
      { type: "h2", text: "1. Une gouvernance claire, portée par la direction" },
      { type: "p", text: "La cybersécurité est un risque d'entreprise. Elle doit donc être pilotée comme tel : un responsable identifié, un mandat clair, un comité qui arbitre et un reporting régulier à la direction générale. Sans cette gouvernance, les projets de sécurité restent dispersés et les arbitrages budgétaires difficiles." },
      { type: "ul", items: ["Nommer un responsable de la sécurité et définir son périmètre", "Formaliser une politique de sécurité validée par la direction", "Mettre en place un comité de pilotage et des indicateurs lisibles"] },
      { type: "h2", text: "2. Connaître ce que l'on protège" },
      { type: "p", text: "On ne peut pas protéger un actif que l'on ignore. L'inventaire des applications, des données sensibles, des serveurs, des postes, des comptes Cloud et des prestataires est la base de toute démarche. Il permet ensuite d'identifier les actifs critiques pour l'activité et de concentrer les efforts." },
      { type: "h2", text: "3. Maîtriser les identités et les accès" },
      { type: "p", text: "Une grande partie des incidents exploite des identifiants compromis ou des droits trop larges. L'authentification multifacteur sur les accès sensibles, la revue régulière des droits et la gestion stricte des comptes à privilèges figurent parmi les mesures les plus rentables." },
      { type: "h2", text: "4. Détecter et réagir" },
      { type: "p", text: "Aucune protection n'est infaillible. Il faut donc être capable de détecter un comportement anormal et de réagir rapidement : journalisation centralisée, surveillance des événements critiques et procédures de réponse documentées." },
      { type: "h2", text: "5. Se préparer à la crise" },
      { type: "p", text: "Une cyberattaque majeure est une crise d'entreprise. Les organisations qui s'en sortent le mieux sont celles qui s'y sont préparées : sauvegardes isolées et testées, plan de continuité, cellule de crise identifiée et exercices réguliers." },
      { type: "h2", text: "Et ensuite ?" },
      { type: "p", text: "Ces cinq priorités constituent un socle. Elles doivent ensuite être déclinées selon les risques propres à chaque organisation, idéalement à partir d'une analyse de risques. C'est cette démarche qui permet de passer d'une sécurité subie à une sécurité pilotée." },
    ],
  },
  {
    slug: "ebios-rm-comprendre-analyse-de-risques",
    title: "EBIOS RM : comprendre l'analyse de risques",
    excerpt:
      "Une méthode en cinq ateliers pour identifier les scénarios de risque les plus réalistes et décider des mesures de sécurité prioritaires.",
    category: "Risk Management",
    date: "2026-10-07",
    readingMinutes: 7,
    author: "Équipe WOOLYS",
    keywords: ["EBIOS RM", "analyse de risques", "ANSSI", "EBIOS RM Sénégal"],
    content: [
      { type: "p", text: "EBIOS Risk Manager est la méthode d'analyse des risques cyber publiée par l'ANSSI, l'agence nationale de la sécurité des systèmes d'information en France. Utilisée bien au-delà des frontières françaises, elle s'est imposée comme une référence dans l'espace francophone pour sa capacité à relier les enjeux métier aux mesures de sécurité." },
      { type: "h2", text: "Une approche centrée sur les scénarios réalistes" },
      { type: "p", text: "Plutôt que de dresser une liste exhaustive de vulnérabilités, EBIOS RM part de ce qui compte pour l'organisation (ses valeurs métier) et des acteurs qui pourraient chercher à lui nuire. L'objectif est de construire des scénarios plausibles et de concentrer les efforts sur les plus critiques." },
      { type: "h2", text: "Les cinq ateliers de la méthode" },
      { type: "h3", text: "Atelier 1 — Cadrage et socle de sécurité" },
      { type: "p", text: "Définition du périmètre, des missions, des valeurs métier et des biens supports, puis évaluation du socle de sécurité existant au regard des référentiels applicables." },
      { type: "h3", text: "Atelier 2 — Sources de risque" },
      { type: "p", text: "Identification des sources de risque (cybercriminels, concurrents, acteurs étatiques, initiés…) et de leurs objectifs visés, afin de retenir les couples les plus pertinents." },
      { type: "h3", text: "Atelier 3 — Scénarios stratégiques" },
      { type: "p", text: "Construction de scénarios de haut niveau décrivant comment une source de risque pourrait atteindre son objectif, y compris en passant par l'écosystème : fournisseurs, partenaires, clients." },
      { type: "h3", text: "Atelier 4 — Scénarios opérationnels" },
      { type: "p", text: "Déclinaison technique des scénarios stratégiques : modes opératoires, biens supports exploités, et évaluation de la vraisemblance." },
      { type: "h3", text: "Atelier 5 — Traitement du risque" },
      { type: "p", text: "Définition de la stratégie de traitement, des mesures de sécurité et du suivi des risques résiduels, présentés de manière compréhensible par la direction." },
      { type: "h2", text: "Pourquoi adopter EBIOS RM ?" },
      { type: "ul", items: ["Un langage commun entre métiers, techniciens et direction", "Une prise en compte explicite de l'écosystème et des tiers", "Des livrables directement exploitables pour décider et prioriser", "Une démarche compatible avec ISO 27001 et ISO 27005"] },
      { type: "p", text: "Bien menée, une analyse EBIOS RM n'est pas un exercice documentaire : c'est un outil de décision. Elle permet de justifier les investissements de sécurité et de les orienter vers ce qui protège réellement l'activité." },
    ],
  },
  {
    slug: "cybersecurite-ia-nouveaux-enjeux-entreprises",
    title: "Cybersécurité et IA : les nouveaux enjeux pour les entreprises",
    excerpt:
      "Assistants, systèmes RAG et agents autonomes ouvrent de nouvelles opportunités, mais aussi de nouvelles surfaces d'attaque à maîtriser.",
    category: "AI Security",
    date: "2026-10-12",
    readingMinutes: 6,
    author: "Équipe WOOLYS",
    keywords: ["sécurité IA", "AI security", "LLM security", "RAG", "agents IA"],
    content: [
      { type: "p", text: "L'intelligence artificielle générative s'est installée dans les entreprises à une vitesse inédite. Assistants conversationnels, recherche augmentée dans les documents internes, automatisation par agents : les usages se multiplient, souvent plus vite que les règles qui devraient les encadrer." },
      { type: "h2", text: "Premier enjeu : la maîtrise des données" },
      { type: "p", text: "Lorsqu'un collaborateur colle un contrat, un fichier client ou du code source dans un outil d'IA public, ces informations quittent le périmètre de l'organisation. La première mesure consiste donc à définir quels outils sont autorisés et quelles données ne doivent jamais y être partagées." },
      { type: "h2", text: "Deuxième enjeu : les vulnérabilités propres aux LLM" },
      { type: "p", text: "Les applications fondées sur des modèles de langage présentent des risques spécifiques, recensés notamment par l'OWASP Top 10 for LLM Applications :" },
      { type: "ul", items: ["L'injection de prompt, directe ou cachée dans un document traité", "La divulgation d'informations sensibles dans les réponses", "La confiance excessive accordée aux sorties du modèle", "Les permissions trop larges accordées aux outils connectés"] },
      { type: "h2", text: "Troisième enjeu : les systèmes RAG" },
      { type: "p", text: "Les architectures RAG (Retrieval-Augmented Generation) permettent à un modèle d'interroger les documents de l'entreprise. Si les droits d'accès ne sont pas respectés lors de la recherche, un utilisateur peut obtenir des informations auxquelles il n'aurait jamais dû avoir accès. Le cloisonnement et la traçabilité sont essentiels." },
      { type: "h2", text: "Quatrième enjeu : les agents autonomes" },
      { type: "p", text: "Les agents ne se contentent plus de répondre : ils agissent, envoient des messages, modifient des données ou déclenchent des processus. Il devient indispensable d'appliquer le principe du moindre privilège, de prévoir une validation humaine pour les actions sensibles et de journaliser chaque action." },
      { type: "h2", text: "Construire une gouvernance de l'IA" },
      { type: "p", text: "Des cadres comme le NIST AI Risk Management Framework ou la norme ISO/IEC 42001 offrent des repères pour structurer cette gouvernance : inventaire des cas d'usage, analyse de risques, règles d'usage, responsabilités et suivi. L'objectif n'est pas de freiner l'innovation, mais de permettre une adoption de l'IA durable et digne de confiance." },
    ],
  },
];

export const getArticle = (slug: string) => articles.find((item) => item.slug === slug);

export const sortedArticles = () => [...articles].sort((a, b) => b.date.localeCompare(a.date));
