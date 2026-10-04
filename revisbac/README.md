# RéviBac 🇸🇳

Application mobile de révision gamifiée pour les élèves sénégalais qui préparent le **BFM** (3e) et le **Baccalauréat** (Terminale S et L).

Fiches synthétiques, flashcards, quiz (QCM, vrai/faux, textes à trous), défi du jour, examens blancs notés sur 20, séries de jours, niveaux et badges. Tout fonctionne **hors-ligne**.

➡️ Le concept détaillé, les améliorations proposées et la feuille de route : [docs/CONCEPT.md](docs/CONCEPT.md).

## Lancer l'application

Prérequis : Node.js ≥ 22 et l'application **Expo Go** sur le téléphone (Android ou iOS).

```bash
cd revisbac
npm install
npx expo start        # puis scanner le QR code avec Expo Go
npx expo start --web  # ou tester dans le navigateur
```

## Vérifications

```bash
npm run typecheck   # TypeScript
npm run validate    # cohérence du contenu (ids uniques, réponses valides, trous…)
npm test            # règles du jeu (XP, séries, gels, badges, notes)
```

## Organisation du code

```
src/
  app/                    écrans (Expo Router : chaque fichier = un écran)
    (tabs)/               Accueil, Réviser, Défis, Profil
    onboarding.tsx        choix du prénom, de l'examen et de l'objectif
    matiere/[id].tsx      chapitres d'une matière
    fiche/[id].tsx        fiche de révision
    flashcards/[id].tsx   cartes mémoire
    quiz.tsx              quiz de chapitre, défi du jour, examen blanc, révision des erreurs
  components/             composants d'interface (questions, blocs de fiche, récompenses…)
  data/
    types.ts              modèle du contenu pédagogique
    catalog.ts            examens, index des matières/chapitres/questions
    content/              une matière = un fichier (ex. maths-s.ts)
  lib/                    règles du jeu, construction des quiz, dates, aléatoire
  state/progress.tsx      progression de l'élève (sauvegardée sur le téléphone)
```

## Ajouter ou corriger du contenu

1. Ouvrir (ou créer) le fichier de la matière dans `src/data/content/`, en suivant les types de `src/data/types.ts`.
2. Pour une nouvelle matière, l'ajouter à la liste de `src/data/content/index.ts`.
3. Lancer `npm run validate`.

> ⚠️ Le contenu fourni est une première base rédigée d'après les grandes lignes du programme sénégalais. Il doit être **relu et validé par des enseignants** avant une diffusion auprès des élèves.

## Publier sur les stores

Avec [EAS](https://docs.expo.dev/eas/) : `npx eas-cli@latest build -p android` (APK/AAB pour le Play Store) et `-p ios`.
