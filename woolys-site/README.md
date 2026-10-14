# WOOLYS — Site vitrine

Site vitrine de **WOOLYS — Cybersecurity & Trust** (« La confiance au cœur du numérique »).
Cabinet de conseil et de services numériques spécialisé dans la cybersécurité, la gestion des risques et la confiance numérique.

Ce site est la **vitrine commerciale** (prospection, prise de contact). La formation en ligne est assurée par la plateforme
**LMS WOOLYS Academy** (projet séparé), reliée au site via la variable `NEXT_PUBLIC_LMS_URL`.

## Stack

- **Next.js 16** (App Router, Turbopack, pages statiques générées au build)
- **TypeScript**
- **Tailwind CSS v4** (design system dans `src/app/globals.css`)
- **Framer Motion** (chargé en `LazyMotion` pour limiter le poids du bundle, respecte `prefers-reduced-motion`)

Aucune autre dépendance d'exécution.

## Charte graphique

| Couleur        | Hex       | Usage                                   |
| -------------- | --------- | --------------------------------------- |
| Vert profond   | `#0B3D33` | Couleur institutionnelle, fonds sombres |
| Vert émeraude  | `#158F7A` | Innovation, boutons principaux          |
| Or champagne   | `#D4AF7C` | **Accent uniquement** (lignes, hover…)  |
| Beige clair    | `#F5F1E8` | Fonds clairs                            |

Typographies : **Montserrat** (titres, logotype) et **Inter** (texte courant), servies localement par `next/font`.
Des variantes accessibles (`emerald-600`, `gold-700`) sont utilisées pour les petits textes afin de respecter les contrastes WCAG AA.

L'emblème (le « W » entouré d'empreintes) est une **recréation vectorielle** de la charte : `src/components/Logo.tsx`,
`public/brand/woolys-emblem.svg` et `src/app/icon.svg` (favicon). Si le designer fournit le fichier vectoriel officiel,
remplacez ces fichiers. La géométrie est générée par `npm run gen:emblem` (sortie dans `scripts/out/`).

## Structure

```
woolys-site/
├── public/brand/            Logo / emblème SVG
├── scripts/gen-emblem.mjs   Générateur de la géométrie de l'emblème
└── src/
    ├── app/                 Routes (App Router — remplace l'ancien dossier pages/)
    │   ├── page.tsx                 Accueil
    │   ├── expertises/              /expertises + /expertises/[slug] (6 pages)
    │   ├── advisory/ academy/ labs/ Les 3 pôles
    │   ├── secteurs/                /secteurs + /secteurs/[slug] (8 pages)
    │   ├── methodologies/           8 référentiels
    │   ├── insights/                Blog + /insights/[slug]
    │   ├── a-propos/ contact/       À propos, formulaire de contact
    │   ├── mentions-legales/ confidentialite/
    │   ├── api/contact/route.ts     Réception du formulaire → webhook
    │   ├── sitemap.ts robots.ts manifest.ts opengraph-image.tsx icon.svg
    │   └── layout.tsx globals.css   Layout global + design system
    ├── components/          Composants réutilisables (Navbar, Footer, PageHero, SectionHeader,
    │                        ServiceCard, IndustryCard, MethodologyCard, CTA, ContactForm,
    │                        ArticleCard, Stats, Timeline, FAQ, Breadcrumbs, Logo, Icon…)
    ├── sections/            Sections de la page d'accueil et blocs de pages
    ├── data/                Tout le contenu éditorial (expertises, secteurs, méthodologies,
    │                        pôles, articles, navigation) — typé dans data/types.ts
    └── lib/                 Configuration (site.ts), SEO, i18n, validation du formulaire, carte d'Afrique
```

**Modifier le contenu** : tout se trouve dans `src/data/`. Ajouter un article = ajouter un objet dans
`src/data/insights.ts` (page, sitemap et homepage sont mis à jour automatiquement).

## Commandes

```bash
cd woolys-site
npm install          # installation (Node.js ≥ 20.9)
npm run dev          # développement → http://localhost:3000
npm run build        # build de production
npm run start        # serveur de production (après le build)
npm run lint         # ESLint
npm run typecheck    # génération des types de routes + TypeScript
```

## Variables d'environnement

Copier `.env.example` en `.env.local` (développement) ou les définir chez l'hébergeur.

| Variable                    | Obligatoire | Description                                                                 |
| --------------------------- | ----------- | --------------------------------------------------------------------------- |
| `NEXT_PUBLIC_SITE_URL`      | Oui (prod)  | URL publique (canonical, sitemap, Open Graph, Schema.org)                   |
| `CONTACT_WEBHOOK_URL`       | Oui (prod)  | Webhook recevant les demandes de contact en JSON (Make, Zapier, n8n, CRM…)  |
| `NEXT_PUBLIC_LMS_URL`       | Non         | URL de la plateforme LMS WOOLYS Academy (liens masqués si vide)             |
| `NEXT_PUBLIC_CONTACT_EMAIL` | Non         | E-mail affiché (footer, contact, repli en cas d'erreur)                     |
| `NEXT_PUBLIC_BOOKING_URL`   | Non         | Lien de prise de rendez-vous (« Planifier un échange »), sinon → /contact   |
| `NEXT_PUBLIC_LINKEDIN_URL`  | Non         | Page LinkedIn                                                               |

Les variables `NEXT_PUBLIC_*` sont intégrées **au moment du build** : relancer un build après modification.

Sans `CONTACT_WEBHOOK_URL`, l'API `/api/contact` affiche les demandes dans la console en développement et renvoie
une erreur 503 en production (le formulaire propose alors l'e-mail de contact) : aucune demande n'est perdue silencieusement.

## Déploiement

**Vercel (recommandé)**

1. Importer le dépôt, définir *Root Directory* = `woolys-site`.
2. Renseigner les variables d'environnement ci-dessus.
3. Déployer, puis brancher le domaine définitif et mettre à jour `NEXT_PUBLIC_SITE_URL`.

**Serveur Node (VPS, Docker…)**

```bash
npm ci && npm run build
NODE_ENV=production PORT=3000 npm run start   # derrière un reverse proxy HTTPS (Nginx, Caddy)
```

Toutes les pages sont générées statiquement ; seule `/api/contact` nécessite le runtime Node.js.

## Qualité

- **SEO** : metadata par page, Open Graph / Twitter, `sitemap.xml`, `robots.txt`, Schema.org
  (`Organization`, `ProfessionalService`, `WebSite`, `Service`, `Article`, `FAQPage`, `BreadcrumbList`), un seul H1 par page.
- **Accessibilité** : lien d'évitement, navigation clavier (menu mobile avec piège de focus et Échap), focus visibles,
  labels et messages d'erreur reliés aux champs (`aria-invalid`, `aria-describedby`), FAQ native `<details>`, contrastes AA.
- **Responsive** : vérifié à 360, 390, 768, 1024, 1440 et 1920 px (aucun débordement horizontal).
- **Performance** : pages statiques, polices auto-hébergées, aucune image bitmap (visuels en SVG/CSS), animations du hero en CSS pur.
- **Sécurité** : en-têtes HTTP (HSTS, X-Frame-Options, nosniff…), validation serveur, champ piège anti-spam, limitation de débit.

## À compléter avant la mise en ligne (PLACEHOLDERS)

- Raison sociale, mentions légales, hébergeur (`src/app/mentions-legales/page.tsx`, `siteConfig.legalName`)
- Durée de conservation des données et validation juridique (`src/app/confidentialite/page.tsx`)
- Domaine définitif, e-mail de contact, lien de rendez-vous, LinkedIn (variables d'environnement)
- Relecture éditoriale des articles de démarrage (`src/data/insights.ts`)
- Fichier vectoriel officiel du logo (optionnel)
- Tranches de budget du formulaire (`src/data/navigation.ts`)

Conformément au cahier des charges, **aucun client, partenaire, chiffre, certification, témoignage, bureau ni collaborateur
n'a été inventé**. Les référentiels (ISO 27001, EBIOS RM…) sont présentés comme méthodes utilisées, pas comme certifications.

## Internationalisation

Le site est en français. L'anglais est préparé : chaînes d'interface dans `src/lib/dictionaries/`, contenus typés dans
`src/data/`, et procédure décrite dans `src/lib/i18n.ts`.
