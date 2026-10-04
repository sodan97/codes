import type { Subject } from '../types';

const subject: Subject = {
  id: "svt-s",
  name: "SVT",
  icon: "🧬",
  color: "#16A34A",
  tracks: ["bac-s"],
  chapters: [
    // ───────────────────────────── TISSU NERVEUX ─────────────────────────────
    {
      id: "svt-s-tissu-nerveux",
      title: "Le tissu nerveux : potentiels et synapse",
      summary: "Potentiel de repos, potentiel d'action, propagation du message nerveux et transmission synaptique.",
      essentials: [
        "Potentiel de repos ≈ −70 mV : l'intérieur de la fibre est négatif par rapport à l'extérieur.",
        "Potentiel d'action : dépolarisation (entrée de Na⁺) puis repolarisation (sortie de K⁺).",
        "Une fibre isolée obéit à la loi du tout ou rien ; un nerf présente un recrutement.",
        "Le message nerveux est codé en fréquence de potentiels d'action.",
        "À la synapse chimique, le message électrique devient chimique (neurotransmetteur), et la transmission est unidirectionnelle.",
      ],
      sections: [
        {
          title: "Potentiel de repos",
          blocks: [
            { kind: "definition", term: "Potentiel de repos", definition: "Différence de potentiel permanente (≈ −70 mV) entre l'intérieur et l'extérieur d'une cellule nerveuse au repos ; l'intérieur est négatif." },
            {
              kind: "list",
              title: "Origine",
              items: [
                "Répartition inégale des ions : K⁺ plus concentré à l'intérieur, Na⁺ plus concentré à l'extérieur.",
                "Au repos, la membrane est surtout perméable au K⁺ (fuite de K⁺ vers l'extérieur).",
                "La pompe Na⁺/K⁺, qui consomme de l'ATP, maintient ces différences : elle fait sortir 3 Na⁺ et entrer 2 K⁺.",
              ],
            },
          ],
        },
        {
          title: "Potentiel d'action",
          blocks: [
            {
              kind: "list",
              title: "Phases (fibre isolée)",
              items: [
                "Dépolarisation : ouverture des canaux Na⁺ voltage-dépendants, entrée massive de Na⁺ ; le potentiel passe de −70 mV à environ +30 mV.",
                "Repolarisation : fermeture des canaux Na⁺ et ouverture des canaux K⁺, sortie de K⁺.",
                "Hyperpolarisation : brève, puis retour au potentiel de repos.",
              ],
            },
            { kind: "definition", term: "Loi du tout ou rien", definition: "Pour une fibre isolée, une stimulation inférieure au seuil ne déclenche aucun potentiel d'action ; au-dessus du seuil, le potentiel d'action a toujours la même amplitude." },
            { kind: "definition", term: "Période réfractaire", definition: "Courte durée après un potentiel d'action pendant laquelle la fibre ne peut pas être de nouveau excitée (ou difficilement)." },
            { kind: "warning", text: "Un nerf n'obéit pas à la loi du tout ou rien : l'amplitude du potentiel global augmente avec l'intensité de stimulation (recrutement de fibres de seuils différents) jusqu'à un maximum." },
          ],
        },
        {
          title: "Propagation et codage",
          blocks: [
            {
              kind: "list",
              items: [
                "Fibre amyélinisée : propagation de proche en proche par courants locaux, lente.",
                "Fibre myélinisée : propagation saltatoire, de nœud de Ranvier en nœud de Ranvier, rapide.",
                "La vitesse augmente avec le diamètre de la fibre.",
                "Codage : plus la stimulation est intense, plus la fréquence des potentiels d'action est élevée (l'amplitude ne change pas).",
              ],
            },
          ],
        },
        {
          title: "La synapse chimique",
          blocks: [
            {
              kind: "list",
              title: "Étapes de la transmission",
              items: [
                "Arrivée du potentiel d'action dans la terminaison présynaptique.",
                "Entrée d'ions Ca²⁺ dans la terminaison.",
                "Exocytose des vésicules : libération du neurotransmetteur (ex. : acétylcholine) dans la fente synaptique.",
                "Fixation sur des récepteurs spécifiques de la membrane postsynaptique et ouverture de canaux ioniques.",
                "Naissance d'un potentiel postsynaptique ; le neurotransmetteur est ensuite dégradé (ex. : acétylcholinestérase) ou recapté.",
              ],
            },
            { kind: "definition", term: "PPSE / PPSI", definition: "Synapse excitatrice : potentiel postsynaptique excitateur (dépolarisation). Synapse inhibitrice : potentiel postsynaptique inhibiteur (hyperpolarisation)." },
            { kind: "text", text: "Le neurone postsynaptique intègre les messages : sommation spatiale (plusieurs synapses en même temps) et sommation temporelle (potentiels rapprochés sur une même synapse). Un potentiel d'action naît si le seuil est atteint." },
            { kind: "tip", text: "Au niveau de la synapse, le message est codé en concentration de neurotransmetteur. Il existe un délai synaptique d'environ 0,5 ms." },
          ],
        },
      ],
      flashcards: [
        { front: "Valeur du potentiel de repos ?", back: "Environ −70 mV (intérieur négatif)." },
        { front: "Ion responsable de la dépolarisation ?", back: "Na⁺ (entrée dans la cellule)." },
        { front: "Ion responsable de la repolarisation ?", back: "K⁺ (sortie de la cellule)." },
        { front: "Rôle de la pompe Na⁺/K⁺ ?", back: "Maintenir les différences de concentration : sort 3 Na⁺, entre 2 K⁺, avec consommation d'ATP." },
        { front: "Loi du tout ou rien ?", back: "Fibre isolée : pas de PA sous le seuil ; au-dessus, PA d'amplitude constante." },
        { front: "Codage du message nerveux sur une fibre ?", back: "En fréquence de potentiels d'action." },
        { front: "Propagation dans une fibre myélinisée ?", back: "Saltatoire, de nœud de Ranvier en nœud de Ranvier." },
        { front: "Ion qui déclenche l'exocytose du neurotransmetteur ?", back: "Ca²⁺." },
        { front: "Codage du message à la synapse ?", back: "En concentration de neurotransmetteur." },
      ],
      quiz: [
        {
          id: "svt-s-tissu-nerveux-q1",
          type: "qcm",
          prompt: "Quelle est la cause de la phase de dépolarisation du potentiel d'action ?",
          choices: ["Une sortie de K⁺", "Une entrée de Na⁺", "Une entrée de Cl⁻", "Le fonctionnement de la pompe Na⁺/K⁺"],
          answer: 1,
          explanation: "L'ouverture des canaux Na⁺ voltage-dépendants provoque une entrée massive de Na⁺ qui inverse la polarité de la membrane.",
        },
        {
          id: "svt-s-tissu-nerveux-q2",
          type: "vrai-faux",
          prompt: "Un nerf obéit à la loi du tout ou rien.",
          answer: false,
          explanation: "C'est la fibre isolée qui obéit à cette loi. Le nerf, formé de fibres de seuils différents, montre un recrutement progressif.",
        },
        {
          id: "svt-s-tissu-nerveux-q3",
          type: "trous",
          prompt: "Sur une fibre, le message nerveux est codé en ___ de potentiels d'action ; à la synapse, il est codé en ___ de neurotransmetteur.",
          answers: ["fréquence", "concentration"],
          bank: ["fréquence", "concentration", "amplitude", "vitesse", "durée"],
          explanation: "L'amplitude des PA est constante : seule leur fréquence varie. À la synapse, c'est la quantité de neurotransmetteur libéré qui varie.",
        },
        {
          id: "svt-s-tissu-nerveux-q4",
          type: "qcm",
          prompt: "Quel ion entre dans la terminaison présynaptique et déclenche l'exocytose ?",
          choices: ["Ca²⁺", "Na⁺", "K⁺", "Cl⁻"],
          answer: 0,
          explanation: "L'arrivée du PA ouvre des canaux Ca²⁺ ; l'entrée de Ca²⁺ déclenche la fusion des vésicules avec la membrane.",
        },
        {
          id: "svt-s-tissu-nerveux-q5",
          type: "vrai-faux",
          prompt: "La conduction est plus rapide dans une fibre myélinisée que dans une fibre amyélinisée de même diamètre.",
          answer: true,
          explanation: "Dans la fibre myélinisée, la conduction est saltatoire (d'un nœud de Ranvier à l'autre), donc beaucoup plus rapide.",
        },
        {
          id: "svt-s-tissu-nerveux-q6",
          type: "qcm",
          prompt: "Un potentiel postsynaptique inhibiteur (PPSI) correspond à :",
          choices: ["une dépolarisation", "un potentiel d'action", "une hyperpolarisation", "une absence de variation"],
          answer: 2,
          explanation: "Le PPSI éloigne le potentiel membranaire du seuil : c'est une hyperpolarisation.",
        },
        {
          id: "svt-s-tissu-nerveux-q7",
          type: "trous",
          prompt: "Au repos, l'ion ___ est plus concentré dans la cellule et l'ion ___ est plus concentré à l'extérieur.",
          answers: ["K⁺", "Na⁺"],
          bank: ["K⁺", "Na⁺", "myéline", "axone"],
          explanation: "Ces gradients, maintenus par la pompe Na⁺/K⁺, sont à l'origine du potentiel de repos et du potentiel d'action.",
        },
        {
          id: "svt-s-tissu-nerveux-q8",
          type: "vrai-faux",
          prompt: "La transmission synaptique chimique se fait dans les deux sens.",
          answer: false,
          explanation: "Elle est unidirectionnelle : les vésicules sont dans l'élément présynaptique et les récepteurs dans l'élément postsynaptique.",
        },
        {
          id: "svt-s-tissu-nerveux-q9",
          type: "qcm",
          prompt: "Deux PPSE arrivant presque en même temps par deux synapses différentes s'additionnent. On parle de :",
          choices: ["sommation temporelle", "période réfractaire", "recrutement", "sommation spatiale"],
          answer: 3,
          explanation: "Sommation spatiale : plusieurs synapses différentes en même temps. Sommation temporelle : potentiels rapprochés sur une même synapse.",
        },
      ],
    },

    // ───────────────────────────── MUSCLE ─────────────────────────────
    {
      id: "svt-s-muscle",
      title: "Le muscle squelettique et la contraction",
      summary: "Structure de la fibre musculaire, mécanisme moléculaire de la contraction et sources d'énergie.",
      essentials: [
        "Le sarcomère, compris entre deux stries Z, est l'unité de contraction.",
        "La contraction résulte du glissement des filaments d'actine entre les filaments de myosine.",
        "Le Ca²⁺ libéré par le réticulum sarcoplasmique déclenche la contraction.",
        "L'ATP est la seule source d'énergie directement utilisable par le muscle.",
        "L'ATP est régénéré par la phosphocréatine, la fermentation lactique et la respiration.",
      ],
      sections: [
        {
          title: "Organisation du muscle",
          blocks: [
            { kind: "text", text: "Muscle → faisceaux → fibres musculaires (cellules géantes plurinucléées) → myofibrilles → sarcomères." },
            {
              kind: "list",
              title: "Le sarcomère",
              items: [
                "Filaments fins : actine ; filaments épais : myosine.",
                "Bande claire I : actine seule.",
                "Bande sombre A : myosine (et actine dans les zones de chevauchement).",
                "Zone H : partie centrale de la bande A, myosine seule.",
              ],
            },
          ],
        },
        {
          title: "Mécanisme de la contraction",
          blocks: [
            {
              kind: "list",
              title: "Étapes",
              items: [
                "Le message nerveux arrive à la plaque motrice ; l'acétylcholine déclenche un potentiel d'action musculaire.",
                "Le réticulum sarcoplasmique libère des ions Ca²⁺.",
                "Le Ca²⁺ démasque les sites de fixation de la myosine sur l'actine.",
                "Les têtes de myosine se fixent sur l'actine (ponts actomyosine), pivotent et tirent l'actine vers le centre du sarcomère.",
                "La fixation d'une nouvelle molécule d'ATP détache la tête de myosine ; son hydrolyse « réarme » la tête.",
              ],
            },
            { kind: "warning", text: "Lors de la contraction, le sarcomère, la bande I et la zone H raccourcissent, mais la bande A garde la même longueur (les filaments ne raccourcissent pas, ils glissent)." },
            { kind: "tip", text: "Sans ATP, les têtes de myosine restent fixées à l'actine : c'est la rigidité cadavérique." },
          ],
        },
        {
          title: "Énergie de la contraction",
          blocks: [
            { kind: "formula", label: "Hydrolyse de l'ATP", formula: "ATP → ADP + Pi + énergie" },
            {
              kind: "list",
              title: "Régénération de l'ATP",
              items: [
                "Phosphocréatine : voie très rapide mais de courte durée (phosphocréatine + ADP → créatine + ATP).",
                "Fermentation lactique (sans O₂) : glucose → 2 acide lactique + 2 ATP ; rapide, peu rentable.",
                "Respiration (avec O₂) : glucose + 6 O₂ → 6 CO₂ + 6 H₂O + 36 à 38 ATP ; lente mais très rentable.",
              ],
            },
          ],
        },
        {
          title: "Le myogramme",
          blocks: [
            { kind: "definition", term: "Secousse musculaire", definition: "Réponse à une stimulation unique : temps de latence, phase de contraction, phase de relâchement." },
            { kind: "text", text: "Si les stimulations sont rapprochées, les secousses fusionnent : tétanos imparfait (plateau dentelé) puis tétanos parfait (plateau lisse) quand la fréquence augmente encore." },
          ],
        },
      ],
      flashcards: [
        { front: "Unité de contraction du muscle ?", back: "Le sarcomère (entre deux stries Z)." },
        { front: "Filaments fins et épais ?", back: "Fins : actine ; épais : myosine." },
        { front: "Quelles zones raccourcissent lors de la contraction ?", back: "Le sarcomère, la bande I et la zone H (la bande A ne change pas)." },
        { front: "Rôle du Ca²⁺ dans la contraction ?", back: "Il démasque les sites de fixation de la myosine sur l'actine." },
        { front: "Rôle de l'ATP dans la contraction ?", back: "Il permet le détachement des têtes de myosine et fournit l'énergie du pivotement." },
        { front: "Rendement de la respiration et de la fermentation ?", back: "Respiration : 36 à 38 ATP par glucose ; fermentation lactique : 2 ATP." },
        { front: "Voie de régénération de l'ATP la plus rapide ?", back: "La phosphocréatine." },
        { front: "Les 3 phases d'une secousse musculaire ?", back: "Latence, contraction, relâchement." },
      ],
      quiz: [
        {
          id: "svt-s-muscle-q1",
          type: "qcm",
          prompt: "Lors de la contraction, quelle partie du sarcomère garde la même longueur ?",
          choices: ["La bande I", "La zone H", "Le sarcomère entier", "La bande A"],
          answer: 3,
          explanation: "La bande A correspond à la longueur des filaments de myosine, qui ne raccourcissent pas.",
        },
        {
          id: "svt-s-muscle-q2",
          type: "vrai-faux",
          prompt: "Pendant la contraction, les filaments d'actine et de myosine raccourcissent.",
          answer: false,
          explanation: "Les filaments gardent leur longueur ; ils glissent les uns par rapport aux autres.",
        },
        {
          id: "svt-s-muscle-q3",
          type: "trous",
          prompt: "Les ions ___ sont libérés par le réticulum sarcoplasmique ; l'___ permet le détachement des têtes de myosine.",
          answers: ["Ca²⁺", "ATP"],
          bank: ["Ca²⁺", "ATP", "Na⁺", "ADP", "glucose"],
          explanation: "Le Ca²⁺ déclenche la formation des ponts actomyosine ; l'ATP est nécessaire à leur rupture.",
        },
        {
          id: "svt-s-muscle-q4",
          type: "qcm",
          prompt: "Combien d'ATP la fermentation lactique produit-elle par molécule de glucose ?",
          choices: ["2", "36", "6", "0"],
          answer: 0,
          explanation: "La fermentation lactique ne fournit que 2 ATP par glucose, contre 36 à 38 pour la respiration.",
        },
        {
          id: "svt-s-muscle-q5",
          type: "vrai-faux",
          prompt: "La rigidité cadavérique s'explique par l'absence d'ATP, qui empêche le détachement de la myosine.",
          answer: true,
          explanation: "Sans ATP, les ponts actomyosine ne peuvent plus se rompre : le muscle reste contracté et rigide.",
        },
        {
          id: "svt-s-muscle-q6",
          type: "qcm",
          prompt: "Quelle voie métabolique régénère l'ATP le plus rapidement au tout début d'un sprint ?",
          choices: ["La respiration", "La photosynthèse", "La phosphocréatine", "La digestion"],
          answer: 2,
          explanation: "La phosphocréatine cède directement son phosphate à l'ADP : c'est la voie la plus rapide, mais elle s'épuise en quelques secondes.",
        },
        {
          id: "svt-s-muscle-q7",
          type: "trous",
          prompt: "La bande claire I contient seulement de l'___ et la zone H contient seulement de la ___.",
          answers: ["actine", "myosine"],
          bank: ["actine", "myosine", "troponine", "myoglobine"],
          explanation: "Bande I : filaments fins d'actine ; zone H : partie centrale de la bande A, sans actine.",
        },
        {
          id: "svt-s-muscle-q8",
          type: "qcm",
          prompt: "Des stimulations de fréquence très élevée donnent un myogramme en plateau lisse. C'est :",
          choices: ["une secousse isolée", "le tétanos parfait", "le tétanos imparfait", "la fatigue musculaire"],
          answer: 1,
          explanation: "Quand les secousses fusionnent complètement, on obtient un plateau lisse : le tétanos parfait.",
        },
        {
          id: "svt-s-muscle-q9",
          type: "vrai-faux",
          prompt: "La fermentation lactique se déroule en absence de dioxygène.",
          answer: true,
          explanation: "C'est une voie anaérobie : elle produit de l'acide lactique et peu d'ATP.",
        },
      ],
    },

    // ───────────────────────────── GLYCÉMIE ─────────────────────────────
    {
      id: "svt-s-glycemie",
      title: "La régulation de la glycémie",
      summary: "Comment le pancréas, le foie et les hormones maintiennent la glycémie autour de 1 g/L, et ce qui se passe dans le diabète.",
      essentials: [
        "La glycémie normale à jeun est d'environ 1 g/L.",
        "L'insuline (cellules β des îlots de Langerhans) est la seule hormone hypoglycémiante.",
        "Le glucagon (cellules α) est hyperglycémiant.",
        "Le foie stocke le glucose sous forme de glycogène et peut le libérer dans le sang.",
        "Diabète de type 1 : absence d'insuline ; type 2 : insulinorésistance.",
      ],
      sections: [
        {
          title: "Le système de régulation",
          blocks: [
            { kind: "definition", term: "Glycémie", definition: "Concentration du glucose dans le sang : environ 1 g/L à jeun (entre 0,8 et 1,2 g/L environ)." },
            {
              kind: "list",
              items: [
                "Capteurs et centres : les cellules des îlots de Langerhans du pancréas détectent les variations de la glycémie.",
                "Messagers : l'insuline et le glucagon, hormones transportées par le sang.",
                "Effecteurs : foie, muscles, tissu adipeux.",
                "Rétrocontrôle négatif : le retour à la normale arrête la sécrétion de l'hormone.",
              ],
            },
          ],
        },
        {
          title: "Insuline et glucagon",
          blocks: [
            {
              kind: "list",
              title: "Insuline (hypoglycémiante), sécrétée en cas d'hyperglycémie",
              items: [
                "Favorise l'entrée du glucose dans les cellules (muscles, tissu adipeux).",
                "Stimule la glycogénogenèse (glucose → glycogène) dans le foie et les muscles.",
                "Stimule la lipogenèse dans le tissu adipeux.",
              ],
            },
            {
              kind: "list",
              title: "Glucagon (hyperglycémiant), sécrété en cas d'hypoglycémie",
              items: [
                "Stimule la glycogénolyse hépatique (glycogène → glucose).",
                "Stimule la néoglucogenèse (fabrication de glucose à partir d'autres molécules).",
              ],
            },
            { kind: "warning", text: "Les muscles stockent du glycogène mais ne peuvent pas libérer de glucose dans le sang : seul le foie est capable de le faire." },
          ],
        },
        {
          title: "Expériences historiques",
          blocks: [
            { kind: "example", title: "Ablation du pancréas (chien)", text: "Après pancréatectomie : hyperglycémie, glycosurie (glucose dans l'urine), amaigrissement puis mort. Une greffe de pancréas sous la peau (sans connexion nerveuse ni canal) corrige les troubles : le pancréas agit par voie sanguine (rôle endocrine)." },
            { kind: "example", title: "Le foie lavé de Claude Bernard", text: "Un foie lavé jusqu'à disparition du glucose en libère de nouveau quelques heures plus tard : il contient une réserve (le glycogène) transformée en glucose." },
          ],
        },
        {
          title: "Les diabètes",
          blocks: [
            {
              kind: "list",
              items: [
                "Diabète : glycémie à jeun ≥ 1,26 g/L (confirmée par deux mesures).",
                "Type 1 (insulinodépendant) : destruction auto-immune des cellules β ; souvent chez le sujet jeune ; traitement par injections d'insuline.",
                "Type 2 (non insulinodépendant) : les cellules cibles répondent mal à l'insuline (insulinorésistance) ; le plus fréquent ; favorisé par la sédentarité et le surpoids.",
              ],
            },
            { kind: "tip", text: "Glycosurie : le glucose apparaît dans les urines quand la glycémie dépasse le seuil rénal (environ 1,8 g/L)." },
          ],
        },
      ],
      flashcards: [
        { front: "Glycémie normale à jeun ?", back: "Environ 1 g/L." },
        { front: "Seule hormone hypoglycémiante ?", back: "L'insuline." },
        { front: "Cellules sécrétrices d'insuline ?", back: "Les cellules β des îlots de Langerhans." },
        { front: "Cellules sécrétrices de glucagon ?", back: "Les cellules α des îlots de Langerhans." },
        { front: "Glycogénolyse ?", back: "Hydrolyse du glycogène en glucose (stimulée par le glucagon)." },
        { front: "Glycogénogenèse ?", back: "Synthèse de glycogène à partir du glucose (stimulée par l'insuline)." },
        { front: "Seul organe capable de libérer du glucose dans le sang à partir du glycogène ?", back: "Le foie." },
        { front: "Différence diabète type 1 / type 2 ?", back: "Type 1 : pas d'insuline (cellules β détruites) ; type 2 : insulinorésistance." },
      ],
      quiz: [
        {
          id: "svt-s-glycemie-q1",
          type: "qcm",
          prompt: "Après un repas riche en glucides, quelle hormone est sécrétée en plus grande quantité ?",
          choices: ["Le glucagon", "L'adrénaline", "L'insuline", "La testostérone"],
          answer: 2,
          explanation: "L'hyperglycémie après le repas stimule les cellules β, qui sécrètent l'insuline pour faire baisser la glycémie.",
        },
        {
          id: "svt-s-glycemie-q2",
          type: "vrai-faux",
          prompt: "Les muscles peuvent libérer du glucose dans le sang à partir de leur glycogène.",
          answer: false,
          explanation: "Les muscles utilisent leur glycogène pour eux-mêmes. Seul le foie peut libérer du glucose dans le sang.",
        },
        {
          id: "svt-s-glycemie-q3",
          type: "trous",
          prompt: "L'insuline est sécrétée par les cellules ___ et le glucagon par les cellules ___ des îlots de Langerhans.",
          answers: ["β", "α"],
          bank: ["β", "α", "acineuses", "hépatiques"],
          explanation: "Cellules β → insuline (hypoglycémiante) ; cellules α → glucagon (hyperglycémiant).",
        },
        {
          id: "svt-s-glycemie-q4",
          type: "qcm",
          prompt: "Le glucagon augmente la glycémie principalement en stimulant :",
          choices: ["la glycogénolyse hépatique", "la lipogenèse", "la glycogénogenèse musculaire", "l'absorption intestinale"],
          answer: 0,
          explanation: "Le glucagon agit sur le foie : il stimule l'hydrolyse du glycogène en glucose, libéré dans le sang.",
        },
        {
          id: "svt-s-glycemie-q5",
          type: "vrai-faux",
          prompt: "Le diabète de type 1 est dû à une destruction des cellules β du pancréas.",
          answer: true,
          explanation: "C'est une maladie auto-immune : sans cellules β, il n'y a plus d'insuline, d'où le traitement par injections d'insuline.",
        },
        {
          id: "svt-s-glycemie-q6",
          type: "qcm",
          prompt: "Un chien dont on a retiré le pancréas reçoit une greffe de pancréas sous la peau, sans connexion nerveuse. Sa glycémie redevient normale. On en déduit que :",
          choices: ["le pancréas agit par voie nerveuse", "le pancréas n'intervient pas", "le foie remplace le pancréas", "le pancréas agit par voie sanguine (hormones)"],
          answer: 3,
          explanation: "Sans nerfs ni canal, seul le sang relie la greffe à l'organisme : le pancréas agit par des hormones.",
        },
        {
          id: "svt-s-glycemie-q7",
          type: "trous",
          prompt: "On parle de diabète quand la glycémie à jeun est supérieure ou égale à ___ ; le diabète de type ___ est le plus fréquent.",
          answers: ["1,26 g/L", "2"],
          bank: ["1,26 g/L", "2", "1", "0,8 g/L", "3"],
          explanation: "Le seuil de diagnostic est 1,26 g/L à jeun. Le diabète de type 2, lié à l'insulinorésistance, est le plus répandu.",
        },
        {
          id: "svt-s-glycemie-q8",
          type: "vrai-faux",
          prompt: "La régulation de la glycémie fait intervenir un rétrocontrôle négatif.",
          answer: true,
          explanation: "Quand la glycémie revient à la valeur normale, la sécrétion de l'hormone correctrice diminue.",
        },
        {
          id: "svt-s-glycemie-q9",
          type: "qcm",
          prompt: "La présence de glucose dans les urines s'appelle :",
          choices: ["glycogénolyse", "glycosurie", "néoglucogenèse", "hypoglycémie"],
          answer: 1,
          explanation: "La glycosurie apparaît quand la glycémie dépasse le seuil rénal (environ 1,8 g/L).",
        },
      ],
    },

    // ───────────────────────────── IMMUNOLOGIE ─────────────────────────────
    {
      id: "svt-s-immunologie",
      title: "Immunologie et VIH",
      summary: "Le soi et le non-soi, l'immunité non spécifique, les réponses spécifiques humorale et cellulaire, et l'infection par le VIH.",
      essentials: [
        "Le soi est défini par les marqueurs du CMH (système HLA chez l'homme).",
        "La réponse non spécifique (inflammation, phagocytose) est immédiate.",
        "Réponse humorale : LB → plasmocytes → anticorps ; réponse cellulaire : LT8 → LT cytotoxiques.",
        "Les LT4 (auxiliaires) coordonnent les deux réponses grâce aux interleukines.",
        "Le VIH détruit les LT4 : le système immunitaire s'effondre (sida).",
      ],
      sections: [
        {
          title: "Soi, non-soi et défenses non spécifiques",
          blocks: [
            { kind: "definition", term: "Soi", definition: "Ensemble des molécules propres à un individu, notamment les marqueurs du CMH (complexe majeur d'histocompatibilité), appelés HLA chez l'homme. Ils expliquent le rejet des greffes." },
            { kind: "definition", term: "Antigène", definition: "Molécule reconnue comme étrangère (non-soi) et capable de déclencher une réponse immunitaire spécifique." },
            {
              kind: "list",
              title: "Phagocytose (macrophages, granulocytes)",
              items: ["Adhésion.", "Ingestion (formation d'une vésicule de phagocytose).", "Digestion par les enzymes des lysosomes.", "Rejet des déchets."],
            },
          ],
        },
        {
          title: "Réponse spécifique humorale",
          blocks: [
            { kind: "text", text: "Les lymphocytes B reconnaissent l'antigène grâce à leurs anticorps membranaires. Activés (avec l'aide des LT4), ils se multiplient (sélection clonale, amplification) et se différencient en plasmocytes qui sécrètent des anticorps circulants." },
            { kind: "definition", term: "Anticorps (immunoglobuline)", definition: "Protéine en forme de Y formée de 2 chaînes lourdes et 2 chaînes légères. Ses deux sites de fixation (partie variable) se lient spécifiquement à l'antigène : complexe immun, ensuite éliminé par phagocytose." },
            { kind: "tip", text: "La réponse humorale est efficace contre les antigènes libres dans les liquides de l'organisme (toxines, bactéries, virus circulants)." },
          ],
        },
        {
          title: "Réponse spécifique cellulaire et coopération",
          blocks: [
            {
              kind: "list",
              items: [
                "Le macrophage, cellule présentatrice de l'antigène (CPA), présente des fragments d'antigène associés au CMH.",
                "Les LT4 reconnaissent l'antigène présenté, se multiplient et sécrètent des interleukines (dont l'IL-2).",
                "Les interleukines activent les LB et les LT8.",
                "Les LT8 deviennent des LT cytotoxiques (LTc) qui détruisent les cellules infectées (perforine).",
              ],
            },
            { kind: "text", text: "Origine : tous les lymphocytes naissent dans la moelle osseuse ; les LB y acquièrent leur maturité, les LT dans le thymus." },
            { kind: "definition", term: "Mémoire immunitaire", definition: "Après un premier contact, des lymphocytes mémoire persistent : la réponse secondaire est plus rapide, plus intense et plus durable. C'est la base de la vaccination." },
          ],
        },
        {
          title: "Le VIH et le sida",
          blocks: [
            {
              kind: "list",
              items: [
                "Le VIH est un rétrovirus : son génome est de l'ARN, copié en ADN par la transcriptase inverse puis intégré à l'ADN de la cellule hôte.",
                "Cibles : cellules portant le récepteur CD4, surtout les LT4, et les macrophages.",
                "Primo-infection : multiplication du virus, production d'anticorps anti-VIH (séropositivité).",
                "Phase asymptomatique : longue, le nombre de LT4 diminue lentement.",
                "Phase sida : LT4 très bas, apparition de maladies opportunistes (tuberculose, candidoses…).",
              ],
            },
            { kind: "warning", text: "Être séropositif signifie avoir des anticorps anti-VIH dans le sang : c'est la preuve d'une infection, pas d'une protection." },
            { kind: "tip", text: "Les antirétroviraux (trithérapie) bloquent des enzymes virales comme la transcriptase inverse ; ils ne guérissent pas mais contrôlent l'infection." },
          ],
        },
      ],
      flashcards: [
        { front: "Marqueurs du soi chez l'homme ?", back: "Les molécules du CMH (système HLA)." },
        { front: "Étapes de la phagocytose ?", back: "Adhésion, ingestion, digestion, rejet." },
        { front: "Cellules qui sécrètent les anticorps ?", back: "Les plasmocytes (issus des LB)." },
        { front: "Structure d'un anticorps ?", back: "En Y : 2 chaînes lourdes + 2 chaînes légères, 2 sites de fixation de l'antigène." },
        { front: "Rôle des LT4 ?", back: "Lymphocytes auxiliaires : sécrètent des interleukines qui activent LB et LT8." },
        { front: "Rôle des LT cytotoxiques ?", back: "Détruire les cellules infectées (perforine)." },
        { front: "Lieu de maturation des LT ?", back: "Le thymus." },
        { front: "Cellules cibles du VIH ?", back: "Les cellules CD4 : surtout les LT4, et les macrophages." },
        { front: "Que signifie séropositif ?", back: "Présence d'anticorps anti-VIH dans le sang." },
      ],
      quiz: [
        {
          id: "svt-s-immunologie-q1",
          type: "qcm",
          prompt: "Quelles cellules sécrètent les anticorps circulants ?",
          choices: ["Les LT8", "Les macrophages", "Les LT4", "Les plasmocytes"],
          answer: 3,
          explanation: "Les LB activés se différencient en plasmocytes, véritables usines à anticorps.",
        },
        {
          id: "svt-s-immunologie-q2",
          type: "vrai-faux",
          prompt: "Le VIH infecte principalement les lymphocytes T4.",
          answer: true,
          explanation: "Le VIH se fixe sur le récepteur CD4 présent surtout à la surface des LT4 (et des macrophages).",
        },
        {
          id: "svt-s-immunologie-q3",
          type: "trous",
          prompt: "Les LT4 sécrètent des ___ qui activent les LB et les LT8 ; les LT8 deviennent des lymphocytes T ___.",
          answers: ["interleukines", "cytotoxiques"],
          bank: ["interleukines", "cytotoxiques", "anticorps", "auxiliaires", "histamines"],
          explanation: "La coopération cellulaire passe par les interleukines des LT4 ; les LTc tuent les cellules infectées.",
        },
        {
          id: "svt-s-immunologie-q4",
          type: "qcm",
          prompt: "Où les lymphocytes T acquièrent-ils leur maturité ?",
          choices: ["Dans le thymus", "Dans la moelle osseuse", "Dans la rate", "Dans le foie"],
          answer: 0,
          explanation: "Tous les lymphocytes naissent dans la moelle osseuse ; les LT migrent et mûrissent dans le thymus.",
        },
        {
          id: "svt-s-immunologie-q5",
          type: "vrai-faux",
          prompt: "Une personne séropositive au VIH est protégée contre le virus grâce à ses anticorps.",
          answer: false,
          explanation: "La séropositivité indique une infection. Les anticorps anti-VIH ne suffisent pas à éliminer le virus.",
        },
        {
          id: "svt-s-immunologie-q6",
          type: "qcm",
          prompt: "Quelle enzyme permet au VIH de copier son ARN en ADN ?",
          choices: ["L'ADN polymérase humaine", "L'amylase", "La transcriptase inverse", "La perforine"],
          answer: 2,
          explanation: "La transcriptase inverse, enzyme virale, fabrique un ADN à partir de l'ARN viral. Certains antirétroviraux la bloquent.",
        },
        {
          id: "svt-s-immunologie-q7",
          type: "trous",
          prompt: "Un anticorps est formé de 2 chaînes ___ et de 2 chaînes ___.",
          answers: ["lourdes", "légères"],
          bank: ["lourdes", "légères", "doubles", "glucidiques"],
          explanation: "L'anticorps (immunoglobuline) a une forme de Y : 2 chaînes lourdes et 2 chaînes légères, reliées par des ponts disulfures.",
        },
        {
          id: "svt-s-immunologie-q8",
          type: "qcm",
          prompt: "La réponse secondaire, lors d'un deuxième contact avec le même antigène, est :",
          choices: ["plus lente et plus faible", "identique à la première", "absente", "plus rapide et plus intense"],
          answer: 3,
          explanation: "Grâce aux cellules mémoire, la réponse secondaire est plus rapide, plus forte et plus durable : c'est le principe du rappel de vaccin.",
        },
        {
          id: "svt-s-immunologie-q9",
          type: "vrai-faux",
          prompt: "La phagocytose est une réponse immunitaire non spécifique.",
          answer: true,
          explanation: "Les phagocytes ingèrent de nombreux types d'éléments étrangers sans reconnaissance spécifique d'un antigène précis.",
        },
      ],
    },

    // ───────────────────────────── REPRODUCTION ─────────────────────────────
    {
      id: "svt-s-reproduction",
      title: "Reproduction et régulation hormonale",
      summary: "Fonctionnement des testicules et des ovaires, cycles sexuels féminins et contrôle par l'axe hypothalamo-hypophysaire.",
      essentials: [
        "L'hypothalamus sécrète la GnRH qui stimule la sécrétion de FSH et LH par l'hypophyse antérieure.",
        "Chez l'homme : LH → testostérone (cellules de Leydig) ; FSH + testostérone → spermatogenèse.",
        "Chez la femme : phase folliculaire (œstrogènes), ovulation déclenchée par le pic de LH, phase lutéale (progestérone).",
        "Les hormones ovariennes exercent un rétrocontrôle négatif, sauf avant l'ovulation (rétrocontrôle positif).",
        "En cas de grossesse, l'HCG de l'embryon maintient le corps jaune.",
      ],
      sections: [
        {
          title: "Chez l'homme",
          blocks: [
            {
              kind: "list",
              items: [
                "Tubes séminifères : spermatogenèse, continue de la puberté à la fin de la vie.",
                "Cellules de Leydig (interstitielles) : sécrètent la testostérone.",
                "Cellules de Sertoli : nourrissent les cellules germinales et contrôlent la spermatogenèse.",
              ],
            },
            { kind: "formula", label: "Axe de commande", formula: "Hypothalamus (GnRH) → hypophyse (FSH, LH) → testicules", note: "LH agit sur les cellules de Leydig ; FSH agit sur les cellules de Sertoli." },
            { kind: "text", text: "La testostérone exerce un rétrocontrôle négatif sur l'hypothalamus et l'hypophyse : son taux reste à peu près constant." },
          ],
        },
        {
          title: "Chez la femme : les cycles",
          blocks: [
            {
              kind: "list",
              title: "Cycle ovarien (28 jours en moyenne)",
              items: [
                "Phase folliculaire (J1 à J14) : croissance des follicules ; les œstrogènes augmentent.",
                "Ovulation (vers J14) : libération de l'ovocyte II, déclenchée par le pic de LH.",
                "Phase lutéale (J14 à J28) : le follicule rompu devient le corps jaune, qui sécrète progestérone et œstrogènes.",
              ],
            },
            {
              kind: "list",
              title: "Cycle utérin",
              items: [
                "Menstruations : élimination de la couche superficielle de l'endomètre.",
                "Phase proliférative : l'endomètre s'épaissit sous l'effet des œstrogènes.",
                "Phase sécrétoire : formation de la dentelle utérine sous l'effet de la progestérone.",
              ],
            },
            { kind: "tip", text: "À la fin du cycle sans fécondation, le corps jaune régresse : la chute de la progestérone déclenche les règles." },
          ],
        },
        {
          title: "Régulation chez la femme",
          blocks: [
            {
              kind: "list",
              items: [
                "Rétrocontrôle négatif : des taux faibles ou moyens d'œstrogènes (et la progestérone en phase lutéale) freinent la sécrétion de GnRH, FSH et LH.",
                "Rétrocontrôle positif : un taux élevé d'œstrogènes, en fin de phase folliculaire, provoque le pic de LH, donc l'ovulation.",
              ],
            },
            { kind: "warning", text: "Ne pas confondre : c'est le pic de LH (et non la progestérone) qui déclenche l'ovulation." },
            { kind: "example", title: "La pilule contraceptive", text: "Elle contient des hormones de synthèse (œstrogènes et/ou progestatifs) qui exercent un rétrocontrôle négatif permanent : pas de pic de LH, donc pas d'ovulation." },
          ],
        },
        {
          title: "Gamétogenèse et grossesse",
          blocks: [
            { kind: "definition", term: "Méiose", definition: "Deux divisions successives qui donnent, à partir d'une cellule diploïde (2n), des cellules haploïdes (n) : les gamètes." },
            { kind: "text", text: "L'ovogenèse commence avant la naissance ; les ovocytes restent bloqués en prophase I. À l'ovulation, l'ovocyte II est libéré bloqué en métaphase II ; la méiose ne s'achève qu'en cas de fécondation." },
            { kind: "text", text: "Après la nidation, l'embryon sécrète l'HCG qui maintient le corps jaune actif : la progestérone reste élevée et il n'y a pas de règles. Les tests de grossesse détectent l'HCG." },
          ],
        },
      ],
      flashcards: [
        { front: "Hormone hypothalamique qui commande l'hypophyse ?", back: "La GnRH." },
        { front: "Hormones hypophysaires (gonadostimulines) ?", back: "FSH et LH." },
        { front: "Cellules sécrétrices de testostérone ?", back: "Les cellules de Leydig (interstitielles)." },
        { front: "Qu'est-ce qui déclenche l'ovulation ?", back: "Le pic de LH." },
        { front: "Hormone sécrétée surtout par le corps jaune ?", back: "La progestérone." },
        { front: "Quand y a-t-il rétrocontrôle positif ?", back: "En fin de phase folliculaire : un taux élevé d'œstrogènes provoque le pic de LH." },
        { front: "Hormone détectée par un test de grossesse ?", back: "L'HCG (sécrétée par l'embryon)." },
        { front: "Mode d'action de la pilule ?", back: "Rétrocontrôle négatif permanent : pas de pic de LH, donc pas d'ovulation." },
        { front: "Où est bloqué l'ovocyte libéré à l'ovulation ?", back: "En métaphase II de méiose." },
      ],
      quiz: [
        {
          id: "svt-s-reproduction-q1",
          type: "qcm",
          prompt: "Quelle hormone déclenche l'ovulation ?",
          choices: ["La progestérone", "Le pic de LH", "La testostérone", "L'HCG"],
          answer: 1,
          explanation: "Un pic de LH, provoqué par un taux élevé d'œstrogènes (rétrocontrôle positif), déclenche l'ovulation.",
        },
        {
          id: "svt-s-reproduction-q2",
          type: "vrai-faux",
          prompt: "La testostérone est sécrétée par les cellules de Sertoli.",
          answer: false,
          explanation: "La testostérone est sécrétée par les cellules de Leydig (interstitielles), sous l'action de la LH.",
        },
        {
          id: "svt-s-reproduction-q3",
          type: "trous",
          prompt: "L'hypothalamus sécrète la ___, qui stimule la sécrétion de FSH et de ___ par l'hypophyse.",
          answers: ["GnRH", "LH"],
          bank: ["GnRH", "LH", "HCG", "progestérone", "insuline"],
          explanation: "Axe hypothalamo-hypophysaire : GnRH → FSH et LH → gonades.",
        },
        {
          id: "svt-s-reproduction-q4",
          type: "qcm",
          prompt: "Pendant la phase lutéale, l'hormone ovarienne dominante est :",
          choices: ["la LH", "la FSH", "la GnRH", "la progestérone"],
          answer: 3,
          explanation: "Le corps jaune sécrète beaucoup de progestérone (et des œstrogènes), qui prépare l'endomètre à une éventuelle nidation.",
        },
        {
          id: "svt-s-reproduction-q5",
          type: "vrai-faux",
          prompt: "Un taux élevé d'œstrogènes en fin de phase folliculaire exerce un rétrocontrôle positif sur l'hypothalamus et l'hypophyse.",
          answer: true,
          explanation: "C'est ce rétrocontrôle positif qui provoque le pic de LH et l'ovulation.",
        },
        {
          id: "svt-s-reproduction-q6",
          type: "qcm",
          prompt: "Quelle hormone maintient le corps jaune en début de grossesse ?",
          choices: ["L'HCG", "La FSH", "La testostérone", "La GnRH"],
          answer: 0,
          explanation: "L'HCG sécrétée par l'embryon maintient le corps jaune, qui continue à produire de la progestérone : pas de règles.",
        },
        {
          id: "svt-s-reproduction-q7",
          type: "trous",
          prompt: "L'endomètre s'épaissit sous l'effet des ___ pendant la phase proliférative ; la dentelle utérine se forme sous l'effet de la ___.",
          answers: ["œstrogènes", "progestérone"],
          bank: ["œstrogènes", "progestérone", "testostérone", "LH", "insuline"],
          explanation: "Œstrogènes : prolifération de l'endomètre. Progestérone : développement de la dentelle utérine (phase sécrétoire).",
        },
        {
          id: "svt-s-reproduction-q8",
          type: "qcm",
          prompt: "La pilule contraceptive empêche l'ovulation car elle :",
          choices: ["détruit les ovocytes", "provoque un pic de LH", "exerce un rétrocontrôle négatif permanent", "bloque les trompes"],
          answer: 2,
          explanation: "Les hormones de synthèse freinent en permanence l'hypothalamus et l'hypophyse : il n'y a pas de pic de LH.",
        },
        {
          id: "svt-s-reproduction-q9",
          type: "vrai-faux",
          prompt: "La chute du taux de progestérone en fin de cycle déclenche les règles.",
          answer: true,
          explanation: "Sans fécondation, le corps jaune régresse ; la baisse de progestérone entraîne la desquamation de l'endomètre.",
        },
      ],
    },

    // ───────────────────────────── GÉNÉTIQUE ─────────────────────────────
    {
      id: "svt-s-genetique",
      title: "Génétique : lois de Mendel et hérédité",
      summary: "Monohybridisme, dihybridisme (gènes indépendants ou liés), hérédité liée au sexe et hérédité humaine.",
      essentials: [
        "1re loi de Mendel : uniformité des hybrides de 1re génération (F1).",
        "2e loi : pureté des gamètes, chaque gamète ne reçoit qu'un allèle de chaque gène.",
        "Monohybridisme avec dominance : F2 = 3/4 – 1/4 ; codominance : 1/4 – 1/2 – 1/4.",
        "Dihybridisme, gènes indépendants : F2 = 9/16 – 3/16 – 3/16 – 1/16 ; test-cross = 4 × 1/4.",
        "Hérédité liée à X : les résultats des croisements réciproques diffèrent.",
      ],
      sections: [
        {
          title: "Vocabulaire",
          blocks: [
            { kind: "definition", term: "Allèle", definition: "Version d'un gène. Un individu diploïde possède deux allèles de chaque gène (un sur chaque chromosome homologue)." },
            { kind: "definition", term: "Génotype / phénotype", definition: "Génotype : ensemble des allèles d'un individu. Phénotype : caractère observable." },
            { kind: "definition", term: "Homozygote / hétérozygote", definition: "Homozygote : deux allèles identiques. Hétérozygote : deux allèles différents." },
            { kind: "definition", term: "Dominance / codominance", definition: "Un allèle dominant s'exprime même à l'état hétérozygote. En cas de codominance, les deux allèles s'expriment (phénotype intermédiaire ou mixte)." },
          ],
        },
        {
          title: "Monohybridisme",
          blocks: [
            { kind: "text", text: "Croisement de deux lignées pures qui diffèrent par un seul caractère. F1 uniforme (1re loi). F1 × F1 donne la F2." },
            {
              kind: "list",
              items: [
                "Dominance complète : F2 = 3/4 phénotype dominant + 1/4 phénotype récessif.",
                "Codominance : F2 = 1/4 – 1/2 – 1/4 (trois phénotypes).",
                "Test-cross (F1 × homozygote récessif) : 1/2 – 1/2, ce qui prouve que F1 est hétérozygote.",
              ],
            },
            { kind: "tip", text: "Le croisement-test (test-cross) avec un individu homozygote récessif révèle directement les gamètes produits par l'individu testé." },
          ],
        },
        {
          title: "Dihybridisme",
          blocks: [
            {
              kind: "list",
              title: "Gènes indépendants (sur des chromosomes différents)",
              items: [
                "F2 : 9/16 – 3/16 – 3/16 – 1/16.",
                "Test-cross : 4 phénotypes équiprobables (1/4 chacun).",
              ],
            },
            {
              kind: "list",
              title: "Gènes liés (sur le même chromosome)",
              items: [
                "Liaison totale : test-cross → 2 phénotypes parentaux seulement (1/2 – 1/2).",
                "Liaison partielle (crossing-over) : 4 phénotypes, les parentaux majoritaires, les recombinés minoritaires (< 50 %).",
                "Distance entre deux gènes (en centimorgans, cM) = pourcentage de recombinés au test-cross.",
              ],
            },
            { kind: "warning", text: "Au test-cross, 4 phénotypes en proportions inégales (2 grandes classes, 2 petites) indiquent des gènes liés avec crossing-over, et non des gènes indépendants." },
          ],
        },
        {
          title: "Hérédité liée au sexe et hérédité humaine",
          blocks: [
            {
              kind: "list",
              items: [
                "Gène porté par le chromosome X (partie propre) : les mâles XY n'ont qu'un allèle (hémizygotes).",
                "Indice : les résultats des croisements réciproques sont différents.",
                "Récessif lié à X (hémophilie, daltonisme) : touche surtout les garçons ; un garçon reçoit son X de sa mère.",
                "Un père transmet son X à toutes ses filles et son Y à tous ses fils.",
              ],
            },
            { kind: "tip", text: "Arbre généalogique : deux parents sains ayant un enfant malade ⇒ l'allèle de la maladie est récessif." },
            { kind: "example", title: "La drépanocytose", text: "Maladie autosomique récessive fréquente en Afrique de l'Ouest. Deux parents hétérozygotes AS ont, à chaque naissance, un risque de 1/4 d'avoir un enfant SS (malade), 1/2 AS et 1/4 AA." },
          ],
        },
      ],
      flashcards: [
        { front: "1re loi de Mendel ?", back: "Uniformité des hybrides de première génération (F1)." },
        { front: "2e loi de Mendel ?", back: "Pureté des gamètes : chaque gamète reçoit un seul allèle de chaque gène." },
        { front: "Proportions en F2 (monohybridisme, dominance) ?", back: "3/4 – 1/4." },
        { front: "Proportions en F2 (monohybridisme, codominance) ?", back: "1/4 – 1/2 – 1/4." },
        { front: "Proportions en F2 (dihybridisme, gènes indépendants) ?", back: "9/16 – 3/16 – 3/16 – 1/16." },
        { front: "Test-cross avec deux gènes indépendants ?", back: "4 phénotypes à 1/4 chacun." },
        { front: "Signe d'une hérédité liée au sexe ?", back: "Les croisements réciproques donnent des résultats différents." },
        { front: "À qui un père transmet-il son chromosome X ?", back: "À toutes ses filles." },
        { front: "Risque d'enfant SS pour deux parents AS (drépanocytose) ?", back: "1/4 à chaque naissance." },
      ],
      quiz: [
        {
          id: "svt-s-genetique-q1",
          type: "qcm",
          prompt: "En monohybridisme avec dominance complète, la F2 donne les proportions phénotypiques :",
          choices: ["1/2 – 1/2", "1/4 – 1/2 – 1/4", "3/4 – 1/4", "9/16 – 7/16"],
          answer: 2,
          explanation: "F1 (Aa) × F1 (Aa) : 1/4 AA + 1/2 Aa + 1/4 aa, soit 3/4 de phénotype dominant et 1/4 de phénotype récessif.",
        },
        {
          id: "svt-s-genetique-q2",
          type: "vrai-faux",
          prompt: "Un test-cross donnant 4 phénotypes en proportions égales (1/4 chacun) indique deux gènes indépendants.",
          answer: true,
          explanation: "L'hybride produit 4 types de gamètes équiprobables, révélés directement par le test-cross.",
        },
        {
          id: "svt-s-genetique-q3",
          type: "trous",
          prompt: "En dihybridisme avec gènes indépendants, la F2 donne ___ ; avec gènes liés et crossing-over, le test-cross donne une majorité de phénotypes ___.",
          answers: ["9/16 – 3/16 – 3/16 – 1/16", "parentaux"],
          bank: ["9/16 – 3/16 – 3/16 – 1/16", "parentaux", "recombinés", "3/4 – 1/4", "1/4 – 1/2 – 1/4"],
          explanation: "Gènes indépendants : 9/3/3/1 en F2. Gènes liés : les recombinés (issus du crossing-over) sont minoritaires.",
        },
        {
          id: "svt-s-genetique-q4",
          type: "qcm",
          prompt: "Une femme conductrice de l'hémophilie (XᴴXʰ) et un homme sain (XᴴY) ont un fils. Quelle est la probabilité qu'il soit hémophile ?",
          choices: ["1/2", "1/4", "0", "1"],
          answer: 0,
          explanation: "Le fils reçoit le Y du père et un X de la mère : Xʰ avec une probabilité 1/2.",
        },
        {
          id: "svt-s-genetique-q5",
          type: "vrai-faux",
          prompt: "Un père atteint d'une maladie récessive liée à X la transmet directement à ses fils.",
          answer: false,
          explanation: "Le père donne son Y à ses fils, pas son X. Il transmet son X (porteur) à toutes ses filles, qui seront au moins conductrices.",
        },
        {
          id: "svt-s-genetique-q6",
          type: "qcm",
          prompt: "Au test-cross de deux gènes liés, on obtient 42 % + 42 % de parentaux et 8 % + 8 % de recombinés. La distance entre les deux gènes est :",
          choices: ["8 cM", "42 cM", "84 cM", "16 cM"],
          answer: 3,
          explanation: "Distance = pourcentage total de recombinés = 8 + 8 = 16 %, soit 16 cM.",
        },
        {
          id: "svt-s-genetique-q7",
          type: "trous",
          prompt: "La 1re loi de Mendel est la loi d'___ des hybrides de F1 ; la 2e est la loi de ___ des gamètes.",
          answers: ["uniformité", "pureté"],
          bank: ["uniformité", "pureté", "dominance", "diversité", "recombinaison"],
          explanation: "Uniformité de la F1 issue de deux lignées pures ; pureté des gamètes (un seul allèle par gène).",
        },
        {
          id: "svt-s-genetique-q8",
          type: "vrai-faux",
          prompt: "Deux parents sains qui ont un enfant malade indiquent que l'allèle de la maladie est récessif.",
          answer: true,
          explanation: "Les parents sont alors hétérozygotes : ils portent l'allèle sans l'exprimer, ce qui caractérise un allèle récessif.",
        },
        {
          id: "svt-s-genetique-q9",
          type: "qcm",
          prompt: "Deux parents hétérozygotes AS pour la drépanocytose attendent un enfant. Probabilité qu'il soit AS ?",
          choices: ["1/4", "1/2", "3/4", "0"],
          answer: 1,
          explanation: "AS × AS : 1/4 AA, 1/2 AS, 1/4 SS. La probabilité d'un enfant AS est 1/2.",
        },
        {
          id: "svt-s-genetique-q10",
          type: "qcm",
          prompt: "En cas de codominance, la F2 d'un monohybridisme présente :",
          choices: ["2 phénotypes : 3/4 – 1/4", "1 seul phénotype", "3 phénotypes : 1/4 – 1/2 – 1/4", "4 phénotypes égaux"],
          answer: 2,
          explanation: "Les hétérozygotes ont un phénotype propre : on observe 3 phénotypes dans les proportions 1/4 – 1/2 – 1/4.",
        },
      ],
    },
  ],
};

export default subject;
