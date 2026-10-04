# RéviBac — concept produit

Application mobile de **révision** (pas de cours complets) pour les élèves sénégalais qui préparent le **BFM** (3e) et le **Baccalauréat** (Terminale S et L). Elle se joue comme un jeu, pour donner envie de réviser un peu chaque jour.

## 1. L'idée de départ

- Toutes les matières du programme sénégalais : maths, sciences physiques, SVT, histoire, géographie, philosophie, français, anglais, espagnol…
- Pour chaque chapitre, seulement l'essentiel : résumés (histoire-géo), fiches de formules (maths, physique), règles et méthodes (langues, philo).
- Des exercices pour vérifier qu'on a compris : QCM, quiz, phrases à compléter.
- De la motivation : points bonus, défis, révision quotidienne, esprit « jeu ».

## 2. Ce que j'ai ajouté pour l'améliorer

| Idée | Pourquoi |
|---|---|
| **« L'essentiel en 30 secondes »** en tête de chaque fiche | Une révision express juste avant l'examen ou dans le car rapide. |
| **Blocs typés dans les fiches** : formule encadrée, définition, date clé, méthode 💡, piège à éviter ⚠️, exemple ✏️ | L'œil retrouve vite l'information ; les « pièges » ciblent les erreurs fréquentes des correcteurs. |
| **Flashcards** avec répétition (les cartes ratées reviennent) | La mémorisation active est plus efficace que la relecture. |
| **3 types d'exercices** : QCM, vrai/faux, texte à trous (banque de mots) | Variés et utilisables d'une main sur téléphone, sans clavier ni accents à taper. |
| **Correction expliquée** après chaque question | On apprend de ses erreurs au lieu de seulement voir « faux ». |
| **« Revoir mes erreurs »** à intervalles croissants (Leitner 3 boîtes) | Une question ratée revient le lendemain, puis 3 et 7 jours après chaque réussite, avant de sortir de la liste. Les réponses comptent dans tous les modes (quiz, défi, examen blanc), et tout est revu au plus tard 2 jours avant l'examen. |
| **Défi du jour** identique pour tous les candidats d'un même examen | On peut comparer son score avec ses camarades → émulation. |
| **Examen blanc** chronométré, noté sur 20 avec mention (Passable → Très Bien) | Se mettre en condition d'examen avec le barème sénégalais. |
| **Série de jours 🔥 + gels de série 🧊** | Habitude quotidienne ; le gel (gagné tous les 7 jours) évite de tout perdre pour un jour manqué (coupure, maladie…). |
| **Objectif quotidien** réglable (30 → 150 XP) | Chaque élève choisit son rythme. |
| **Niveaux** de « Débutant » à « Lauréat du Concours général » et **16 badges** (dont « Lion de la Teranga » pour 30 jours de série) | Progression visible et clin d'œil local. |
| **Compte à rebours J-xx** avant l'examen | Crée un sentiment d'urgence positif. |
| **Conseil méthodo du jour** | Apprendre à réviser, pas seulement quoi réviser. |
| **100 % hors-ligne**, aucune inscription | Les données mobiles coûtent cher : tout le contenu est embarqué, la progression est stockée sur le téléphone. |

## 3. Parcours de l'élève

1. **Accueil (onboarding)** : prénom → examen (BFM, Bac S, Bac L) → objectif quotidien.
2. **Accueil** : série, niveau, objectif du jour, défi du jour, J-xx, suggestion de fiche, erreurs à revoir, conseil.
3. **Réviser** : matières → chapitres → *Fiche* → *Flashcards* → *Quiz*.
4. **Défis** : défi du jour, revoir mes erreurs, statistiques de la semaine, examens blancs par matière.
5. **Profil** : statistiques, badges, réglages (objectif, date d'examen, changer d'examen).

## 4. Règles du jeu (XP)

| Action | XP |
|---|---|
| Bonne réponse | +10 |
| Quiz sans faute (≥ 5 questions) | +20 |
| Quiz de chapitre refait le même jour : XP ÷ 2 | +5 par bonne réponse, sans bonus |
| Première lecture d'une fiche | +10 |
| Flashcards : 1×/jour/chapitre | +5 |
| Défi du jour (1×/jour, rejeu sans XP) | +50 |
| Examen blanc : bonus 1×/jour/matière | +30 |
| Objectif quotidien atteint | +20 |

Toute activité terminée compte pour la série 🔥, même sans XP (un quiz raté, un paquet de flashcards refait). L'objectif du jour, lui, se calcule sur l'XP.

Niveau *n* atteint à 50 × (n−1) × n XP (100, 300, 600, 1 000…).

## 5. Feuille de route proposée

**V1 — prototype actuel (ce dépôt)** : application hors-ligne, 18 matières (BFM, Bac S, Bac L), fiches, flashcards, quiz, défis, badges.

**V2 — contenu et qualité**
- Faire **relire et compléter le contenu par des enseignants** sénégalais (un chapitre = une fiche validée).
- Couvrir tout le programme officiel, séries S1/S2/S3 et L1/L2/L’ distinguées quand les programmes diffèrent.
- **Annales corrigées** du BFM et du Bac (sujets des années précédentes, découpés par chapitre).
- Formules mathématiques mises en forme (rendu LaTeX) et schémas (SVT, physique, cartes en géographie).
- Matières supplémentaires : économie, éducation civique, arabe, portugais, allemand…

**V3 — social et suivi**
- Comptes élèves (facultatifs) et **classements** entre amis / par établissement / par région.
- Défis entre amis (« duel » de 5 questions).
- Notifications de rappel quotidien.
- Espace professeur : créer des quiz pour sa classe et suivre la progression.
- Fiches audio (y compris en wolof pour les explications) pour réviser en écoutant.

**Modèle économique possible** : gratuit pour l'essentiel, offre premium (annales corrigées, examens blancs illimités) payable par Orange Money / Wave, ou partenariats avec des écoles et ONG éducatives.
