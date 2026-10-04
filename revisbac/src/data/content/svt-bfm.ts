import type { Subject } from '../types';

const subject: Subject = {
  id: "svt-bfm",
  name: "SVT",
  icon: "🧬",
  color: "#16A34A",
  tracks: ["bfm"],
  chapters: [
    // ───────────────────────────── NUTRITION ─────────────────────────────
    {
      id: "svt-bfm-nutrition",
      title: "Digestion et absorption",
      summary: "Comment les aliments sont transformés en nutriments dans le tube digestif, puis passent dans le sang et la lymphe.",
      essentials: [
        "La digestion transforme les aliments en nutriments grâce à des actions mécaniques et chimiques.",
        "Les enzymes digestives sont spécifiques et agissent le mieux à 37 °C.",
        "Amidon → glucose ; protides → acides aminés ; lipides → acides gras + glycérol.",
        "L'absorption se fait surtout dans l'intestin grêle, au niveau des villosités.",
        "Une alimentation équilibrée apporte aliments énergétiques, bâtisseurs et fonctionnels.",
      ],
      sections: [
        {
          title: "Le tube digestif et les sucs",
          blocks: [
            { kind: "text", text: "Trajet des aliments : bouche → pharynx → œsophage → estomac → intestin grêle → gros intestin → anus." },
            {
              kind: "list",
              title: "Sucs digestifs",
              items: [
                "Salive (glandes salivaires) : amylase salivaire, amidon → maltose.",
                "Suc gastrique (estomac) : pepsine, protéines → peptides, en milieu acide.",
                "Suc pancréatique (pancréas) : amylase, lipase, trypsine.",
                "Bile (foie) : émulsionne les lipides ; elle ne contient pas d'enzyme.",
                "Suc intestinal : achève la digestion.",
              ],
            },
            { kind: "warning", text: "La bile n'est pas une enzyme : elle fragmente les graisses en fines gouttelettes pour faciliter l'action de la lipase." },
          ],
        },
        {
          title: "Les enzymes digestives",
          blocks: [
            { kind: "definition", term: "Enzyme", definition: "Protéine qui accélère une réaction chimique sans être consommée. Chaque enzyme agit sur une seule substance (spécificité)." },
            { kind: "example", title: "Expérience de l'amidon et de la salive", text: "Empois d'amidon + salive à 37 °C : après quelques minutes, l'eau iodée ne bleuit plus (l'amidon a disparu) et la liqueur de Fehling chauffée donne un précipité rouge brique (présence de sucre réducteur, le maltose)." },
            { kind: "tip", text: "Un témoin sans salive est indispensable : il montre que c'est bien la salive qui transforme l'amidon." },
          ],
        },
        {
          title: "Absorption des nutriments",
          blocks: [
            { kind: "definition", term: "Absorption intestinale", definition: "Passage des nutriments de l'intestin grêle vers le sang et la lymphe à travers la paroi des villosités." },
            {
              kind: "list",
              items: [
                "Glucose, acides aminés, eau, sels minéraux → sang (puis veine porte → foie).",
                "Acides gras et glycérol → surtout la lymphe.",
                "L'intestin grêle est long et plissé (villosités) : sa surface d'échange est très grande.",
                "Le gros intestin absorbe surtout l'eau.",
              ],
            },
          ],
        },
        {
          title: "Hygiène alimentaire",
          blocks: [
            {
              kind: "list",
              title: "Groupes d'aliments",
              items: [
                "Énergétiques : glucides et lipides (riz, mil, huile, sucre).",
                "Bâtisseurs : protides (poisson, viande, œufs, niébé, lait).",
                "Fonctionnels (protecteurs) : vitamines et sels minéraux (fruits, légumes).",
              ],
            },
            {
              kind: "list",
              title: "Maladies de carence",
              items: [
                "Kwashiorkor : manque de protéines.",
                "Marasme : manque global de nourriture.",
                "Goitre : manque d'iode.",
                "Scorbut : manque de vitamine C ; rachitisme : manque de vitamine D.",
                "Anémie : souvent un manque de fer.",
              ],
            },
          ],
        },
      ],
      flashcards: [
        { front: "Quelle enzyme de la salive agit sur l'amidon ?", back: "L'amylase salivaire (amidon → maltose)." },
        { front: "Où agit la pepsine ?", back: "Dans l'estomac, en milieu acide, sur les protéines." },
        { front: "Rôle de la bile ?", back: "Émulsionner les lipides (elle ne contient pas d'enzyme)." },
        { front: "Produit final de la digestion de l'amidon ?", back: "Le glucose." },
        { front: "Produits finaux de la digestion des lipides ?", back: "Acides gras et glycérol." },
        { front: "Où se fait l'essentiel de l'absorption ?", back: "Dans l'intestin grêle, au niveau des villosités." },
        { front: "Température optimale des enzymes digestives humaines ?", back: "Environ 37 °C." },
        { front: "Kwashiorkor : carence en quoi ?", back: "En protéines." },
      ],
      quiz: [
        {
          id: "svt-bfm-nutrition-q1",
          type: "qcm",
          prompt: "Quel est le produit final de la digestion des protéines ?",
          choices: ["Le glucose", "Les acides gras", "Le glycérol", "Les acides aminés"],
          answer: 3,
          explanation: "Les protéines sont découpées en peptides puis en acides aminés, qui sont absorbés dans le sang.",
        },
        {
          id: "svt-bfm-nutrition-q2",
          type: "vrai-faux",
          prompt: "La bile contient une enzyme qui digère les lipides.",
          answer: false,
          explanation: "La bile ne contient pas d'enzyme. Elle émulsionne les lipides ; c'est la lipase (suc pancréatique) qui les digère.",
        },
        {
          id: "svt-bfm-nutrition-q3",
          type: "trous",
          prompt: "L'amylase salivaire transforme l'___ en ___.",
          answers: ["amidon", "maltose"],
          bank: ["amidon", "maltose", "acides aminés", "glycérol", "protéine"],
          explanation: "L'amylase salivaire hydrolyse l'amidon en maltose, un sucre réducteur.",
        },
        {
          id: "svt-bfm-nutrition-q4",
          type: "qcm",
          prompt: "Dans quel organe se fait l'essentiel de l'absorption des nutriments ?",
          choices: ["L'estomac", "L'intestin grêle", "Le gros intestin", "L'œsophage"],
          answer: 1,
          explanation: "L'intestin grêle, grâce à ses nombreuses villosités, offre une très grande surface d'échange.",
        },
        {
          id: "svt-bfm-nutrition-q5",
          type: "vrai-faux",
          prompt: "Une enzyme est spécifique : elle agit sur un seul type de substance.",
          answer: true,
          explanation: "Par exemple, l'amylase agit sur l'amidon mais pas sur les protéines.",
        },
        {
          id: "svt-bfm-nutrition-q6",
          type: "qcm",
          prompt: "Après action de la salive sur l'empois d'amidon à 37 °C, la liqueur de Fehling chauffée donne :",
          choices: ["une coloration bleue", "aucun changement", "un précipité rouge brique", "une coloration violette"],
          answer: 2,
          explanation: "Le précipité rouge brique révèle un sucre réducteur (le maltose) formé à partir de l'amidon.",
        },
        {
          id: "svt-bfm-nutrition-q7",
          type: "trous",
          prompt: "Le poisson et le niébé sont des aliments ___ riches en protides ; le riz et l'huile sont des aliments ___.",
          answers: ["bâtisseurs", "énergétiques"],
          bank: ["bâtisseurs", "énergétiques", "fonctionnels", "toxiques"],
          explanation: "Les protides construisent et réparent le corps ; les glucides et les lipides fournissent l'énergie.",
        },
        {
          id: "svt-bfm-nutrition-q8",
          type: "vrai-faux",
          prompt: "Les acides gras et le glycérol passent principalement dans la lymphe.",
          answer: true,
          explanation: "Les produits de la digestion des lipides empruntent surtout la voie lymphatique, alors que le glucose et les acides aminés passent dans le sang.",
        },
        {
          id: "svt-bfm-nutrition-q9",
          type: "qcm",
          prompt: "Le goitre est dû à une carence en :",
          choices: ["iode", "fer", "vitamine C", "protéines"],
          answer: 0,
          explanation: "Le manque d'iode provoque un gonflement de la glande thyroïde : le goitre. On le prévient avec du sel iodé.",
        },
      ],
    },

    // ───────────────────────────── RESPIRATION ET CIRCULATION ─────────────────────────────
    {
      id: "svt-bfm-respiration-circulation",
      title: "Respiration et circulation sanguine",
      summary: "Les échanges gazeux dans les poumons, le rôle du cœur et le trajet du sang dans l'organisme.",
      essentials: [
        "Les échanges gazeux se font dans les alvéoles pulmonaires : O₂ vers le sang, CO₂ vers l'air.",
        "Les cellules utilisent le dioxygène pour dégrader le glucose et produire de l'énergie.",
        "Le cœur a 4 cavités : 2 oreillettes et 2 ventricules.",
        "Les artères partent du cœur, les veines y reviennent.",
        "Petite circulation : cœur ↔ poumons ; grande circulation : cœur ↔ organes.",
      ],
      sections: [
        {
          title: "La respiration",
          blocks: [
            {
              kind: "list",
              title: "Composition de l'air (valeurs approchées)",
              items: [
                "Air inspiré : environ 21 % de O₂ et 0,03 % de CO₂.",
                "Air expiré : environ 16 % de O₂ et 4 % de CO₂.",
              ],
            },
            { kind: "definition", term: "Alvéole pulmonaire", definition: "Petit sac à paroi très fine, entouré de capillaires sanguins, où se font les échanges gazeux entre l'air et le sang." },
            { kind: "text", text: "Inspiration : le diaphragme se contracte et s'abaisse, les côtes se soulèvent, le volume de la cage thoracique augmente et l'air entre. L'expiration est surtout passive." },
            { kind: "formula", label: "Respiration cellulaire", formula: "glucose + dioxygène → dioxyde de carbone + eau + énergie" },
          ],
        },
        {
          title: "Le sang",
          blocks: [
            {
              kind: "list",
              items: [
                "Plasma : liquide qui transporte nutriments, déchets, CO₂ et hormones.",
                "Globules rouges (hématies) : transportent le dioxygène grâce à l'hémoglobine.",
                "Globules blancs (leucocytes) : défendent l'organisme.",
                "Plaquettes : interviennent dans la coagulation.",
              ],
            },
          ],
        },
        {
          title: "Le cœur et la circulation",
          blocks: [
            { kind: "text", text: "Le cœur est un muscle creux qui fonctionne comme une double pompe. Le cœur droit contient du sang pauvre en O₂, le cœur gauche du sang riche en O₂. Les valvules empêchent le sang de revenir en arrière." },
            {
              kind: "list",
              title: "Les deux circulations",
              items: [
                "Grande circulation : ventricule gauche → aorte → organes → veines caves → oreillette droite.",
                "Petite circulation : ventricule droit → artère pulmonaire → poumons → veines pulmonaires → oreillette gauche.",
              ],
            },
            { kind: "warning", text: "Piège classique : l'artère pulmonaire transporte du sang pauvre en O₂ et les veines pulmonaires du sang riche en O₂. Une artère se définit par le fait qu'elle part du cœur, pas par la richesse en O₂." },
          ],
        },
        {
          title: "Hygiène",
          blocks: [
            { kind: "list", items: ["Le tabac apporte goudrons (cancers), nicotine (dépendance) et monoxyde de carbone (prend la place du O₂ sur l'hémoglobine).", "L'activité physique régulière renforce le cœur et les poumons.", "Une alimentation trop grasse et trop salée favorise les maladies cardiovasculaires."] },
            { kind: "tip", text: "Pendant un effort, le rythme cardiaque et le rythme respiratoire augmentent pour apporter plus de O₂ et de glucose aux muscles." },
          ],
        },
      ],
      flashcards: [
        { front: "Où se font les échanges gazeux respiratoires ?", back: "Dans les alvéoles pulmonaires." },
        { front: "Quel pigment transporte le dioxygène ?", back: "L'hémoglobine des globules rouges." },
        { front: "Combien de cavités possède le cœur ?", back: "4 : 2 oreillettes et 2 ventricules." },
        { front: "Définition d'une artère ?", back: "Vaisseau qui part du cœur et conduit le sang vers les organes." },
        { front: "Trajet de la petite circulation ?", back: "Ventricule droit → artère pulmonaire → poumons → veines pulmonaires → oreillette gauche." },
        { front: "Rôle des valvules cardiaques ?", back: "Empêcher le reflux du sang (circulation à sens unique)." },
        { front: "Équation de la respiration cellulaire ?", back: "Glucose + O₂ → CO₂ + H₂O + énergie." },
        { front: "Rôle des plaquettes ?", back: "La coagulation du sang." },
      ],
      quiz: [
        {
          id: "svt-bfm-respiration-circulation-q1",
          type: "qcm",
          prompt: "Quel vaisseau sort du ventricule gauche ?",
          choices: ["L'artère pulmonaire", "La veine cave", "L'aorte", "La veine pulmonaire"],
          answer: 2,
          explanation: "Le ventricule gauche envoie le sang riche en O₂ dans l'aorte, vers tous les organes (grande circulation).",
        },
        {
          id: "svt-bfm-respiration-circulation-q2",
          type: "vrai-faux",
          prompt: "L'artère pulmonaire transporte du sang riche en dioxygène.",
          answer: false,
          explanation: "Elle conduit du sang pauvre en O₂ du ventricule droit vers les poumons. C'est une artère car elle part du cœur.",
        },
        {
          id: "svt-bfm-respiration-circulation-q3",
          type: "trous",
          prompt: "Au niveau des alvéoles, le ___ passe de l'air vers le sang et le ___ passe du sang vers l'air.",
          answers: ["dioxygène", "dioxyde de carbone"],
          bank: ["dioxygène", "dioxyde de carbone", "glucose", "diazote", "monoxyde de carbone"],
          explanation: "Les gaz diffusent du milieu le plus concentré vers le moins concentré à travers la fine paroi alvéolaire.",
        },
        {
          id: "svt-bfm-respiration-circulation-q4",
          type: "qcm",
          prompt: "Quelle est la teneur approximative en dioxyde de carbone de l'air expiré ?",
          choices: ["0,03 %", "4 %", "21 %", "78 %"],
          answer: 1,
          explanation: "L'air expiré contient environ 4 % de CO₂, contre environ 0,03 % dans l'air inspiré.",
        },
        {
          id: "svt-bfm-respiration-circulation-q5",
          type: "vrai-faux",
          prompt: "Lors de l'inspiration, le diaphragme se contracte et s'abaisse.",
          answer: true,
          explanation: "La contraction du diaphragme augmente le volume de la cage thoracique : l'air entre dans les poumons.",
        },
        {
          id: "svt-bfm-respiration-circulation-q6",
          type: "qcm",
          prompt: "Quelles cellules sanguines transportent le dioxygène ?",
          choices: ["Les globules rouges", "Les globules blancs", "Les plaquettes", "Les cellules du plasma"],
          answer: 0,
          explanation: "Les globules rouges contiennent l'hémoglobine, qui fixe le dioxygène.",
        },
        {
          id: "svt-bfm-respiration-circulation-q7",
          type: "trous",
          prompt: "Le sang revient des organes au cœur par les ___ et arrive dans l'___ droite.",
          answers: ["veines caves", "oreillette"],
          bank: ["veines caves", "oreillette", "aorte", "ventricule", "artères pulmonaires"],
          explanation: "Grande circulation : les veines caves ramènent le sang pauvre en O₂ dans l'oreillette droite.",
        },
        {
          id: "svt-bfm-respiration-circulation-q8",
          type: "vrai-faux",
          prompt: "Le cœur droit contient du sang riche en dioxygène.",
          answer: false,
          explanation: "Le cœur droit contient du sang pauvre en O₂ qui revient des organes ; il l'envoie aux poumons.",
        },
        {
          id: "svt-bfm-respiration-circulation-q9",
          type: "qcm",
          prompt: "Quelle substance de la fumée de tabac prend la place du dioxygène sur l'hémoglobine ?",
          choices: ["La nicotine", "Les goudrons", "La vapeur d'eau", "Le monoxyde de carbone"],
          answer: 3,
          explanation: "Le monoxyde de carbone se fixe sur l'hémoglobine et réduit le transport du dioxygène.",
        },
      ],
    },

    // ───────────────────────────── REPRODUCTION ─────────────────────────────
    {
      id: "svt-bfm-reproduction",
      title: "Reproduction humaine et IST",
      summary: "Les appareils reproducteurs, la fécondation, la grossesse, la contraception et la prévention des IST dont le VIH/sida.",
      essentials: [
        "Testicules → spermatozoïdes ; ovaires → ovules.",
        "La fécondation a lieu dans la trompe et donne une cellule-œuf.",
        "La nidation se fait dans l'utérus environ une semaine après la fécondation.",
        "Le préservatif est le seul moyen de contraception qui protège aussi des IST.",
        "Le VIH se transmet par voie sexuelle, par le sang et de la mère à l'enfant.",
      ],
      sections: [
        {
          title: "Appareils reproducteurs",
          blocks: [
            {
              kind: "list",
              title: "Homme",
              items: [
                "Testicules : produisent les spermatozoïdes et la testostérone.",
                "Spermiductes, vésicules séminales, prostate, urètre, pénis.",
              ],
            },
            {
              kind: "list",
              title: "Femme",
              items: [
                "Ovaires : produisent les ovules et des hormones (œstrogènes, progestérone).",
                "Trompes : lieu de la fécondation.",
                "Utérus : accueille et nourrit l'embryon puis le fœtus.",
                "Vagin.",
              ],
            },
            { kind: "definition", term: "Puberté", definition: "Période de transformation (vers 10 à 16 ans) où les organes reproducteurs deviennent fonctionnels et où apparaissent les caractères sexuels secondaires." },
          ],
        },
        {
          title: "Cycle, fécondation et grossesse",
          blocks: [
            { kind: "text", text: "Le cycle menstruel dure en moyenne 28 jours. Dans un cycle de 28 jours, l'ovulation a lieu vers le 14e jour. S'il n'y a pas de fécondation, la muqueuse de l'utérus est éliminée : ce sont les règles." },
            { kind: "definition", term: "Fécondation", definition: "Union d'un spermatozoïde et d'un ovule dans la trompe. Elle donne une cellule-œuf." },
            { kind: "definition", term: "Nidation", definition: "Implantation de l'embryon dans la muqueuse utérine, environ 7 jours après la fécondation." },
            { kind: "text", text: "La grossesse dure environ 9 mois. Le placenta assure les échanges entre la mère et le fœtus (nutriments, O₂, déchets) sans mélange des sangs." },
            { kind: "warning", text: "L'alcool, le tabac et certains médicaments traversent le placenta et peuvent nuire au fœtus." },
          ],
        },
        {
          title: "Contraception",
          blocks: [
            {
              kind: "list",
              items: [
                "Préservatif (masculin ou féminin) : empêche la rencontre des gamètes et protège des IST.",
                "Pilule : bloque l'ovulation.",
                "Dispositif intra-utérin (stérilet) : empêche la nidation.",
              ],
            },
          ],
        },
        {
          title: "IST et VIH/sida",
          blocks: [
            {
              kind: "list",
              title: "Exemples d'IST",
              items: [
                "Dues à des bactéries : syphilis, gonococcie (blennorragie).",
                "Dues à des virus : VIH/sida, hépatite B, herpès génital.",
              ],
            },
            {
              kind: "list",
              title: "Transmission du VIH",
              items: [
                "Rapports sexuels non protégés.",
                "Sang : objets tranchants souillés, seringues partagées, transfusion non contrôlée.",
                "De la mère à l'enfant : pendant la grossesse, l'accouchement ou l'allaitement.",
              ],
            },
            { kind: "warning", text: "Le VIH ne se transmet PAS par une poignée de main, un repas partagé, les piqûres de moustiques ou l'utilisation des mêmes toilettes." },
            { kind: "tip", text: "Le dépistage se fait par un test sanguin. Les antirétroviraux (ARV) ne guérissent pas le sida mais contrôlent le virus et permettent de vivre longtemps." },
          ],
        },
      ],
      flashcards: [
        { front: "Où a lieu la fécondation ?", back: "Dans la trompe (de Fallope)." },
        { front: "Que donne la fécondation ?", back: "Une cellule-œuf (zygote)." },
        { front: "Jour de l'ovulation dans un cycle de 28 jours ?", back: "Vers le 14e jour." },
        { front: "Qu'est-ce que la nidation ?", back: "L'implantation de l'embryon dans la muqueuse utérine, environ 7 jours après la fécondation." },
        { front: "Rôle du placenta ?", back: "Échanges entre la mère et le fœtus (nutriments, O₂, déchets)." },
        { front: "Quel contraceptif protège aussi des IST ?", back: "Le préservatif." },
        { front: "Trois voies de transmission du VIH ?", back: "Sexuelle, sanguine, de la mère à l'enfant." },
        { front: "Deux IST d'origine bactérienne ?", back: "Syphilis et gonococcie." },
      ],
      quiz: [
        {
          id: "svt-bfm-reproduction-q1",
          type: "qcm",
          prompt: "Quel organe produit les spermatozoïdes ?",
          choices: ["La prostate", "Le testicule", "La vésicule séminale", "L'urètre"],
          answer: 1,
          explanation: "Les testicules produisent les spermatozoïdes et la testostérone.",
        },
        {
          id: "svt-bfm-reproduction-q2",
          type: "vrai-faux",
          prompt: "Le VIH peut se transmettre par une piqûre de moustique.",
          answer: false,
          explanation: "Le moustique ne transmet pas le VIH. Les voies de transmission sont sexuelle, sanguine et de la mère à l'enfant.",
        },
        {
          id: "svt-bfm-reproduction-q3",
          type: "trous",
          prompt: "La fécondation a lieu dans la ___ ; l'embryon s'implante ensuite dans l'___ : c'est la nidation.",
          answers: ["trompe", "utérus"],
          bank: ["trompe", "utérus", "ovaire", "vagin", "placenta"],
          explanation: "Spermatozoïde et ovule se rencontrent dans la trompe ; l'embryon migre ensuite vers l'utérus.",
        },
        {
          id: "svt-bfm-reproduction-q4",
          type: "qcm",
          prompt: "Quelle IST est causée par un virus ?",
          choices: ["La syphilis", "La gonococcie", "L'hépatite B", "Aucune"],
          answer: 2,
          explanation: "L'hépatite B est virale, comme le VIH et l'herpès génital. La syphilis et la gonococcie sont bactériennes.",
        },
        {
          id: "svt-bfm-reproduction-q5",
          type: "vrai-faux",
          prompt: "Le préservatif protège à la fois contre une grossesse et contre les IST.",
          answer: true,
          explanation: "C'est le seul moyen de contraception qui protège aussi des infections sexuellement transmissibles.",
        },
        {
          id: "svt-bfm-reproduction-q6",
          type: "qcm",
          prompt: "Que se passe-t-il à la fin du cycle s'il n'y a pas eu fécondation ?",
          choices: ["L'ovulation", "La nidation", "La puberté", "Les règles"],
          answer: 3,
          explanation: "La muqueuse utérine, épaissie pendant le cycle, est éliminée avec un peu de sang : ce sont les règles.",
        },
        {
          id: "svt-bfm-reproduction-q7",
          type: "trous",
          prompt: "Dans un cycle de ___ jours, l'ovulation a lieu vers le ___ jour.",
          answers: ["28", "14e"],
          bank: ["28", "14e", "7e", "21e", "40"],
          explanation: "Le cycle moyen dure 28 jours ; l'ovule est libéré vers le milieu du cycle, vers le 14e jour.",
        },
        {
          id: "svt-bfm-reproduction-q8",
          type: "vrai-faux",
          prompt: "Les médicaments antirétroviraux guérissent définitivement le sida.",
          answer: false,
          explanation: "Les ARV ne guérissent pas : ils empêchent le virus de se multiplier et permettent à la personne de vivre longtemps en bonne santé.",
        },
        {
          id: "svt-bfm-reproduction-q9",
          type: "qcm",
          prompt: "Quel organe assure les échanges entre la mère et le fœtus ?",
          choices: ["Le placenta", "L'ovaire", "La trompe", "Le vagin"],
          answer: 0,
          explanation: "Le placenta permet le passage du O₂ et des nutriments vers le fœtus, et des déchets vers la mère, sans mélange des sangs.",
        },
      ],
    },

    // ───────────────────────────── SYSTÈME NERVEUX ─────────────────────────────
    {
      id: "svt-bfm-systeme-nerveux",
      title: "Système nerveux et réflexes",
      summary: "L'organisation du système nerveux, le message nerveux, l'arc réflexe et le mouvement volontaire.",
      essentials: [
        "Système nerveux central : encéphale + moelle épinière ; périphérique : les nerfs.",
        "Le neurone est l'unité de base ; il transmet un message nerveux de nature électrique.",
        "Un réflexe est une réponse involontaire, rapide et stéréotypée.",
        "Arc réflexe : récepteur → nerf sensitif → centre nerveux → nerf moteur → effecteur.",
        "Le centre des mouvements volontaires est le cerveau ; chaque hémisphère commande le côté opposé du corps.",
      ],
      sections: [
        {
          title: "Organisation",
          blocks: [
            {
              kind: "list",
              items: [
                "Encéphale (protégé par le crâne) : cerveau, cervelet, tronc cérébral (bulbe rachidien).",
                "Moelle épinière : dans le canal de la colonne vertébrale.",
                "Nerfs : relient les centres nerveux aux organes.",
              ],
            },
            { kind: "definition", term: "Neurone", definition: "Cellule nerveuse formée d'un corps cellulaire, de dendrites et d'un axone. Elle conduit le message nerveux." },
            { kind: "definition", term: "Nerf rachidien", definition: "Nerf mixte relié à la moelle épinière par une racine dorsale (sensitive, avec un ganglion) et une racine ventrale (motrice)." },
          ],
        },
        {
          title: "Le réflexe",
          blocks: [
            { kind: "definition", term: "Réflexe", definition: "Réponse involontaire, rapide et toujours identique à une stimulation. Exemple : retirer la main d'un objet brûlant." },
            {
              kind: "list",
              title: "Les éléments de l'arc réflexe",
              items: [
                "Récepteur sensoriel (ex. : la peau).",
                "Conducteur sensitif (nerf sensitif).",
                "Centre nerveux (la moelle épinière pour les réflexes médullaires).",
                "Conducteur moteur (nerf moteur).",
                "Effecteur (un muscle).",
              ],
            },
            { kind: "example", title: "Expérience de la grenouille spinale", text: "On détruit l'encéphale d'une grenouille : elle réagit encore au pincement d'une patte (réflexe). Si on détruit ensuite la moelle épinière, plus aucun réflexe : la moelle est le centre du réflexe." },
          ],
        },
        {
          title: "Le mouvement volontaire",
          blocks: [
            { kind: "text", text: "Le mouvement volontaire est décidé par le cortex cérébral. Les aires motrices envoient un message vers la moelle épinière, puis vers les muscles par les nerfs moteurs." },
            { kind: "warning", text: "Les voies nerveuses se croisent : l'hémisphère gauche commande le côté droit du corps, et inversement. Une lésion du cerveau gauche peut paralyser le côté droit." },
          ],
        },
        {
          title: "Hygiène du système nerveux",
          blocks: [
            { kind: "list", items: ["Dormir suffisamment.", "Éviter l'alcool, le tabac et les drogues, qui perturbent le fonctionnement des neurones.", "Porter un casque à moto pour protéger l'encéphale."] },
            { kind: "tip", text: "Pour un schéma d'arc réflexe, flèche le trajet du message nerveux dans le bon sens : de la racine dorsale vers la racine ventrale." },
          ],
        },
      ],
      flashcards: [
        { front: "Que comprend le système nerveux central ?", back: "L'encéphale et la moelle épinière." },
        { front: "Unité de base du système nerveux ?", back: "Le neurone." },
        { front: "Caractères d'un réflexe ?", back: "Involontaire, rapide, stéréotypé (toujours identique)." },
        { front: "Les 5 éléments de l'arc réflexe ?", back: "Récepteur, nerf sensitif, centre nerveux, nerf moteur, effecteur." },
        { front: "Centre nerveux du réflexe de retrait de la main ?", back: "La moelle épinière." },
        { front: "Centre nerveux du mouvement volontaire ?", back: "Le cortex cérébral (cerveau)." },
        { front: "Quel hémisphère commande la main droite ?", back: "L'hémisphère gauche." },
        { front: "Racine sensitive du nerf rachidien ?", back: "La racine dorsale (porte le ganglion spinal)." },
      ],
      quiz: [
        {
          id: "svt-bfm-systeme-nerveux-q1",
          type: "qcm",
          prompt: "Quel est le centre nerveux du réflexe de retrait de la main ?",
          choices: ["Le cervelet", "Le cerveau", "La moelle épinière", "Le nerf moteur"],
          answer: 2,
          explanation: "Les réflexes médullaires ont pour centre la moelle épinière, ce qui les rend très rapides.",
        },
        {
          id: "svt-bfm-systeme-nerveux-q2",
          type: "vrai-faux",
          prompt: "Une grenouille dont l'encéphale est détruit peut encore réaliser des réflexes.",
          answer: true,
          explanation: "Tant que la moelle épinière est intacte, les réflexes médullaires restent possibles.",
        },
        {
          id: "svt-bfm-systeme-nerveux-q3",
          type: "trous",
          prompt: "Dans l'arc réflexe, le message va du récepteur au centre nerveux par le nerf ___, puis du centre à l'effecteur par le nerf ___.",
          answers: ["sensitif", "moteur"],
          bank: ["sensitif", "moteur", "optique", "crânien"],
          explanation: "Le nerf sensitif conduit le message vers le centre ; le nerf moteur le conduit vers le muscle.",
        },
        {
          id: "svt-bfm-systeme-nerveux-q4",
          type: "qcm",
          prompt: "Une personne a une lésion de l'hémisphère cérébral gauche. Quelle partie du corps risque d'être paralysée ?",
          choices: ["Le côté gauche", "Le côté droit", "Les deux côtés", "Aucune"],
          answer: 1,
          explanation: "Les voies motrices se croisent : l'hémisphère gauche commande la moitié droite du corps.",
        },
        {
          id: "svt-bfm-systeme-nerveux-q5",
          type: "vrai-faux",
          prompt: "Un réflexe est une réponse volontaire et réfléchie.",
          answer: false,
          explanation: "Un réflexe est involontaire, rapide et stéréotypé. Il se produit sans intervention de la volonté.",
        },
        {
          id: "svt-bfm-systeme-nerveux-q6",
          type: "qcm",
          prompt: "Dans l'arc réflexe, quel élément est l'effecteur ?",
          choices: ["La peau", "La moelle épinière", "Le nerf sensitif", "Le muscle"],
          answer: 3,
          explanation: "L'effecteur est l'organe qui réalise la réponse : ici le muscle qui se contracte.",
        },
        {
          id: "svt-bfm-systeme-nerveux-q7",
          type: "trous",
          prompt: "Le système nerveux central comprend l'___ et la ___.",
          answers: ["encéphale", "moelle épinière"],
          bank: ["encéphale", "moelle épinière", "nerf optique", "colonne vertébrale"],
          explanation: "L'encéphale (cerveau, cervelet, tronc cérébral) et la moelle épinière forment le système nerveux central.",
        },
        {
          id: "svt-bfm-systeme-nerveux-q8",
          type: "vrai-faux",
          prompt: "Le message nerveux est de nature électrique.",
          answer: true,
          explanation: "Le long des neurones, le message nerveux se propage sous forme de signaux électriques.",
        },
        {
          id: "svt-bfm-systeme-nerveux-q9",
          type: "qcm",
          prompt: "Quelle partie d'un neurone conduit le message nerveux loin du corps cellulaire ?",
          choices: ["L'axone", "Le noyau", "La dendrite", "Le ganglion"],
          answer: 0,
          explanation: "Les dendrites reçoivent les messages ; l'axone conduit le message vers d'autres cellules.",
        },
      ],
    },

    // ───────────────────────────── IMMUNITÉ ET SANTÉ ─────────────────────────────
    {
      id: "svt-bfm-immunite-sante",
      title: "Immunité, paludisme et vaccination",
      summary: "Comment l'organisme se défend contre les microbes, comment se transmet le paludisme et comment s'en protéger.",
      essentials: [
        "Les barrières naturelles (peau, muqueuses) empêchent l'entrée des microbes.",
        "La phagocytose est une défense rapide et non spécifique.",
        "Les anticorps, produits à partir des lymphocytes B, sont spécifiques d'un antigène.",
        "Paludisme : parasite Plasmodium, transmis par la piqûre de l'anophèle femelle.",
        "La vaccination est préventive et durable ; la sérothérapie est curative et temporaire.",
      ],
      sections: [
        {
          title: "Les défenses de l'organisme",
          blocks: [
            { kind: "definition", term: "Antigène", definition: "Substance ou élément étranger reconnu par l'organisme et qui déclenche une réaction immunitaire." },
            {
              kind: "list",
              items: [
                "Barrières naturelles : peau, muqueuses, larmes, mucus.",
                "Réaction inflammatoire : rougeur, chaleur, gonflement, douleur.",
                "Phagocytose : des globules blancs (phagocytes) englobent et digèrent les microbes.",
                "Réponse spécifique : lymphocytes B (anticorps) et lymphocytes T.",
              ],
            },
            { kind: "definition", term: "Anticorps", definition: "Protéine produite à partir des lymphocytes B, qui se fixe spécifiquement sur un antigène et aide à le neutraliser." },
          ],
        },
        {
          title: "Le paludisme",
          blocks: [
            {
              kind: "list",
              items: [
                "Agent : un parasite unicellulaire, le Plasmodium (P. falciparum est le plus dangereux).",
                "Vecteur : l'anophèle femelle, qui pique surtout la nuit.",
                "Chez l'homme, le parasite se multiplie d'abord dans le foie, puis dans les globules rouges.",
                "Signes : fièvre, frissons, maux de tête, vomissements.",
              ],
            },
            {
              kind: "list",
              title: "Prévention",
              items: [
                "Dormir sous une moustiquaire imprégnée d'insecticide.",
                "Supprimer les eaux stagnantes où pondent les moustiques.",
                "Traitement préventif chez la femme enceinte et les jeunes enfants.",
              ],
            },
            { kind: "tip", text: "Toute fièvre doit faire penser au paludisme : on fait un test de diagnostic rapide (TDR) au poste de santé avant de traiter." },
          ],
        },
        {
          title: "Vaccination et sérothérapie",
          blocks: [
            { kind: "definition", term: "Vaccination", definition: "Introduction d'un antigène rendu inoffensif (microbe tué ou atténué, toxine modifiée). L'organisme fabrique des anticorps et garde une mémoire : protection préventive, lente à s'installer mais durable." },
            { kind: "definition", term: "Sérothérapie", definition: "Injection d'un sérum contenant des anticorps tout prêts. Action curative, immédiate mais de courte durée." },
            { kind: "warning", text: "Les antibiotiques agissent sur les bactéries mais sont inefficaces contre les virus." },
          ],
        },
      ],
      flashcards: [
        { front: "Agent du paludisme ?", back: "Le Plasmodium, un parasite unicellulaire." },
        { front: "Vecteur du paludisme ?", back: "L'anophèle femelle (moustique)." },
        { front: "Meilleur moyen de prévention du paludisme la nuit ?", back: "Dormir sous une moustiquaire imprégnée d'insecticide." },
        { front: "Qu'est-ce que la phagocytose ?", back: "Englobement et digestion d'un microbe par un phagocyte (globule blanc)." },
        { front: "Les 4 signes de la réaction inflammatoire ?", back: "Rougeur, chaleur, gonflement, douleur." },
        { front: "Qu'est-ce qu'un antigène ?", back: "Un élément étranger qui déclenche une réaction immunitaire." },
        { front: "Vaccination : préventive ou curative ?", back: "Préventive, et durable grâce à la mémoire immunitaire." },
        { front: "Sérothérapie : préventive ou curative ?", back: "Curative, immédiate mais temporaire." },
      ],
      quiz: [
        {
          id: "svt-bfm-immunite-sante-q1",
          type: "qcm",
          prompt: "Quel est le vecteur du paludisme ?",
          choices: ["La mouche tsé-tsé", "L'anophèle femelle", "L'anophèle mâle", "Le Plasmodium"],
          answer: 1,
          explanation: "Le Plasmodium est l'agent (parasite) ; il est transmis par la piqûre de l'anophèle femelle, seule à piquer.",
        },
        {
          id: "svt-bfm-immunite-sante-q2",
          type: "vrai-faux",
          prompt: "Les antibiotiques sont efficaces contre les virus.",
          answer: false,
          explanation: "Les antibiotiques agissent sur les bactéries, pas sur les virus.",
        },
        {
          id: "svt-bfm-immunite-sante-q3",
          type: "trous",
          prompt: "La vaccination est ___ et durable, alors que la sérothérapie est ___ et temporaire.",
          answers: ["préventive", "curative"],
          bank: ["préventive", "curative", "inutile", "héréditaire"],
          explanation: "Le vaccin prépare l'organisme avant l'infection ; le sérum apporte des anticorps pour soigner une infection déjà présente.",
        },
        {
          id: "svt-bfm-immunite-sante-q4",
          type: "qcm",
          prompt: "Quelles cellules sont à l'origine de la production des anticorps ?",
          choices: ["Les globules rouges", "Les plaquettes", "Les neurones", "Les lymphocytes B"],
          answer: 3,
          explanation: "Les lymphocytes B, une fois activés, produisent des anticorps spécifiques de l'antigène.",
        },
        {
          id: "svt-bfm-immunite-sante-q5",
          type: "vrai-faux",
          prompt: "Dans l'organisme humain, le Plasmodium se multiplie dans le foie puis dans les globules rouges.",
          answer: true,
          explanation: "Après la piqûre, le parasite passe d'abord dans le foie, puis envahit et fait éclater les globules rouges, ce qui provoque les accès de fièvre.",
        },
        {
          id: "svt-bfm-immunite-sante-q6",
          type: "qcm",
          prompt: "Rougeur, chaleur, gonflement et douleur sont les signes de :",
          choices: ["la réaction inflammatoire", "la vaccination", "l'allergie au soleil", "la digestion"],
          answer: 0,
          explanation: "Ce sont les quatre signes de la réaction inflammatoire, première réponse de l'organisme après une blessure infectée.",
        },
        {
          id: "svt-bfm-immunite-sante-q7",
          type: "trous",
          prompt: "Les moustiques pondent dans les eaux ___ ; pour se protéger la nuit, on dort sous une ___ imprégnée.",
          answers: ["stagnantes", "moustiquaire"],
          bank: ["stagnantes", "moustiquaire", "courantes", "couverture", "salées"],
          explanation: "Supprimer les eaux stagnantes et utiliser des moustiquaires imprégnées réduit fortement les piqûres d'anophèles.",
        },
        {
          id: "svt-bfm-immunite-sante-q8",
          type: "vrai-faux",
          prompt: "Un anticorps peut neutraliser n'importe quel microbe.",
          answer: false,
          explanation: "Un anticorps est spécifique : il ne reconnaît qu'un antigène précis.",
        },
        {
          id: "svt-bfm-immunite-sante-q9",
          type: "qcm",
          prompt: "Un enfant est mordu par un chien suspect de rage et n'est pas vacciné. Quelle action rapide permet d'apporter immédiatement des anticorps ?",
          choices: ["Un antibiotique", "Une moustiquaire", "La sérothérapie", "Un régime alimentaire"],
          answer: 2,
          explanation: "La sérothérapie apporte des anticorps tout prêts et agit tout de suite ; on l'associe à la vaccination pour une protection durable.",
        },
      ],
    },
  ],
};

export default subject;
