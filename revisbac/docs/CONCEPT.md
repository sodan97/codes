# RéviBac — concept produit

Application mobile de **révision** (pas de cours complets) pour les élèves sénégalais qui préparent le **BFEM** (3e) et le **Baccalauréat** (Terminale S et L). Elle se joue comme un jeu, pour donner envie de réviser un peu chaque jour.

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
| **Flashcards** en répétition espacée (Leitner 5 boîtes : revues à 1, 3, 7, 14 puis 30 jours) | La mémorisation active est plus efficace que la relecture. « Je savais » fait monter la carte d'une boîte (une fois par jour au plus), « À revoir » la renvoie en boîte 1 ; le paquet propose d'abord les cartes à revoir, et l'accueil affiche « 🧠 n cartes à revoir ». Comme pour les erreurs, tout est revu au plus tard 2 jours avant l'examen. |
| **Étoiles de maîtrise** par chapitre : ⭐ fiche lue, ⭐⭐ quiz ≥ 80 %, ⭐⭐⭐ quiz ≥ 80 % sur 2 jours différents et paquet de flashcards terminé | Un objectif clair (« Prochaine étoile : … ») qui récompense la régularité. Après 21 jours sans pratique, les étoiles s'estompent avec un « 🔄 petit rappel conseillé », sans jamais être retirées. La maîtrise d'une matière vaut étoiles / (3 × chapitres). |
| **Quêtes du jour** (une facile, une moyenne, une de variété) et **coffre du jour** | Trois missions différentes chaque jour (bonnes réponses, nouvelle fiche, flashcards, quiz ≥ 80 %, erreurs à corriger, enchaînement, matière la moins travaillée, examen blanc à partir de J-60), tirées selon le jour, l'examen et ce qui a du sens pour l'élève. Les quêtes non faites disparaissent à minuit, sans message d'échec. |
| **3 types d'exercices** : QCM, vrai/faux, texte à trous (banque de mots) | Variés et utilisables d'une main sur téléphone, sans clavier ni accents à taper. |
| **Correction expliquée** après chaque question | On apprend de ses erreurs au lieu de seulement voir « faux ». |
| **« Revoir mes erreurs »** à intervalles croissants (Leitner 3 boîtes) | Une question ratée revient le lendemain, puis 3 et 7 jours après chaque réussite, avant de sortir de la liste. Les réponses comptent dans tous les modes (quiz, défi, examen blanc), et tout est revu au plus tard 2 jours avant l'examen. |
| **Défi du jour** identique pour tous les candidats d'un même examen | On peut comparer son score avec ses camarades → émulation. |
| **Examen blanc** chronométré, noté sur 20 avec mention (Passable → Très Bien) | Se mettre en condition d'examen avec le barème sénégalais. |
| **Série de jours 🔥 + gels de série 🧊** | Habitude quotidienne ; le gel (gagné tous les 7 jours) évite de tout perdre pour un jour manqué (coupure, maladie…). |
| **Objectif quotidien** réglable (30 → 150 XP) | Chaque élève choisit son rythme. |
| **Niveaux** de « Débutant » à « Lauréat du Concours général » et **18 badges** (dont « Lion de la Teranga » pour 30 jours de série, « Chapitre en or » et « Matière maîtrisée ») | Progression visible et clin d'œil local. Les badges verrouillés montrent où l'on en est (« 7/10 fiches »). |
| **Reprise d'un quiz interrompu** (appli fermée par le téléphone, appel, car rapide) | La session est gardée après chaque réponse pendant 24 h : « Reprendre (6/10) » ou « Recommencer ». Un examen blanc dont le temps s'est écoulé appli fermée est noté à la réouverture (questions non traitées comptées 0). |
| **Sauvegarde par code** (« RB2:… »), sans compte ni serveur | Changer de téléphone sans rien perdre : le code (avec une somme de contrôle) se garde dans un message et se colle dans « Restaurer une sauvegarde ». |
| **Compte à rebours J-xx** avant l'examen | Crée un sentiment d'urgence positif. |
| **Conseil méthodo du jour** | Apprendre à réviser, pas seulement quoi réviser. |
| **100 % hors-ligne**, aucune inscription | Les données mobiles coûtent cher : tout le contenu est embarqué, la progression est stockée sur le téléphone. |

## 3. Parcours de l'élève

1. **Accueil (onboarding)** : prénom → examen (BFEM, Bac S, Bac L) → objectif quotidien.
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
| Quête du jour accomplie (3 par jour) | +10 |
| Coffre du jour (les 3 quêtes faites) | +20 à +40, ou un gel de série 🧊 si l'élève en a moins de 2 |

Toute activité terminée compte pour la série 🔥, même sans XP (un quiz raté, un paquet de flashcards refait). L'objectif du jour, lui, se calcule sur l'XP. Les cartes revues une à une ne rapportent pas d'XP : c'est le paquet terminé qui compte.

Le contenu du coffre est tiré d'après le jour et l'examen (même coffre pour un même jour).

**Étoiles d'un chapitre** : ⭐ fiche lue → ⭐⭐ quiz de chapitre ≥ 80 % → ⭐⭐⭐ quiz ≥ 80 % sur 2 jours différents et paquet de flashcards terminé. Chaque étoile suppose la précédente.

**Sauvegarde** : `RB2:` + base64 (UTF-8) de la progression sans l'historique des XP + `:` + somme de contrôle FNV-1a (8 caractères hexadécimaux). Flashcards et erreurs suivies y sont rangées par chapitre, avec des dates relatives au jour de la sauvegarde, pour que le code tienne dans un message (environ 25 000 caractères pour toutes les flashcards du Bac S). Les codes `RB1:` (progression en JSON tel quel) restent lisibles. À l'import, la sauvegarde passe par la même migration que la progression du téléphone : un code d'une ancienne version reste lisible. Erreurs possibles : « Code incomplet », « Code abîmé (somme de contrôle) », « Version inconnue ».

**Reprise** : la session de quiz en cours est gardée sur le téléphone (24 h au plus). Un défi du jour commencé un autre jour n'est pas repris : le défi d'aujourd'hui l'emporte. Examen blanc : le temps continue de s'écouler appli fermée ; s'il est écoulé, l'examen est noté à la réouverture, et au-delà de 30 min après la fin du temps il est noté sans proposition de reprise (pendant 7 jours, ensuite il est effacé). Une seule session est gardée : si l'élève ouvre un autre quiz pendant un examen blanc commencé, l'écran lui propose d'abord de reprendre l'examen ou de le noter.

Niveau *n* atteint à 50 × (n−1) × n XP (100, 300, 600, 1 000…).

## 5. Feuille de route proposée

**V1 — prototype actuel (ce dépôt)** : application hors-ligne, 18 matières (BFEM, Bac S, Bac L), fiches, flashcards, quiz, défis, badges.

**V2 — contenu et qualité**
- Faire **relire et compléter le contenu par des enseignants** sénégalais (un chapitre = une fiche validée).
- Couvrir tout le programme officiel, séries S1/S2/S3 et L1/L2/L’ distinguées quand les programmes diffèrent.
- **Annales corrigées** du BFEM et du Bac (sujets des années précédentes, découpés par chapitre).
- Formules mathématiques mises en forme (rendu LaTeX) et schémas (SVT, physique, cartes en géographie).
- Matières supplémentaires : économie, éducation civique, arabe, portugais, allemand…

**V3 — social et suivi**
- Comptes élèves (facultatifs) et **classements** entre amis / par établissement / par région.
- Défis entre amis (« duel » de 5 questions).
- Notifications de rappel quotidien.
- Espace professeur : créer des quiz pour sa classe et suivre la progression.
- Fiches audio (y compris en wolof pour les explications) pour réviser en écoutant.

**Modèle économique possible** : gratuit pour l'essentiel, offre premium (annales corrigées, examens blancs illimités) payable par Orange Money / Wave, ou partenariats avec des écoles et ONG éducatives.
