import type { Subject } from '../types';

const subject: Subject = {
  id: 'geo-bac',
  name: 'Géographie',
  icon: '🌍',
  color: '#0891B2',
  tracks: ['bac-s', 'bac-l'],
  chapters: [
    // ───────────────────────── Chapitre 1 ─────────────────────────
    {
      id: 'geo-bac-mondialisation',
      title: 'La mondialisation et les grands ensembles économiques',
      summary:
        "La mondialisation intensifie les échanges entre les territoires, organisés autour de pôles dominants et de grands ensembles régionaux.",
      essentials: [
        'La mondialisation est la mise en relation croissante des économies et des sociétés par des flux de marchandises, de capitaux, de personnes et d’informations.',
        'Acteurs : firmes transnationales (FTN), États, organisations internationales (OMC, FMI, Banque mondiale), ONG.',
        "Pôles dominants : la Triade (Amérique du Nord, Europe, Asie orientale) ; montée des pays émergents (BRICS).",
        'Grands ensembles régionaux : UE, ACEUM (ex-ALENA), ASEAN, MERCOSUR, CEDEAO, ZLECAf.',
        "Les inégalités persistent : l'Afrique reste en marge (moins de 3 % du commerce mondial).",
      ],
      sections: [
        {
          title: 'Définitions clés',
          blocks: [
            {
              kind: 'definition',
              term: 'Mondialisation',
              definition: "Processus d'intégration des espaces du monde dans un système unique d'échanges (marchandises, capitaux, personnes, informations).",
            },
            {
              kind: 'definition',
              term: 'Firme transnationale (FTN)',
              definition: "Grande entreprise qui produit et vend dans plusieurs pays, avec une maison mère et des filiales à l'étranger.",
            },
            { kind: 'definition', term: 'IDE (investissement direct à l’étranger)', definition: "Investissement d'une entreprise pour créer ou contrôler une activité dans un autre pays." },
            { kind: 'definition', term: 'Délocalisation', definition: 'Transfert d’une activité de production vers un pays où les coûts (salaires surtout) sont plus faibles.' },
            { kind: 'definition', term: 'Pays émergent', definition: 'Pays en forte croissance qui s’intègre rapidement à l’économie mondiale (Chine, Inde, Brésil…).' },
          ],
        },
        {
          title: 'Les acteurs et les règles',
          blocks: [
            { kind: 'date', date: '1944', event: 'Accords de Bretton Woods : création du FMI et de la Banque mondiale.' },
            { kind: 'date', date: '1947', event: 'Signature du GATT (accord général sur les tarifs douaniers et le commerce).' },
            { kind: 'date', date: '1er janvier 1995', event: "L'OMC (Organisation mondiale du commerce) remplace le GATT ; siège à Genève." },
            {
              kind: 'list',
              title: 'Facteurs de la mondialisation',
              items: [
                'Révolution des transports : conteneurs, porte-conteneurs géants (maritimisation).',
                'Révolution des télécommunications : Internet, satellites.',
                'Libéralisation des échanges (baisse des droits de douane).',
              ],
            },
            {
              kind: 'example',
              title: "L'altermondialisme",
              text: "Mouvement qui veut une autre mondialisation, plus juste. Le Forum social mondial s'est tenu à Dakar en 2011.",
            },
          ],
        },
        {
          title: 'Un monde polycentrique',
          blocks: [
            {
              kind: 'definition',
              term: 'Triade',
              definition: "Les trois pôles dominants de l'économie mondiale : Amérique du Nord, Europe occidentale, Asie orientale (Japon puis Chine).",
            },
            {
              kind: 'list',
              title: 'Les puissances émergentes',
              items: [
                'BRICS : Brésil, Russie, Inde, Chine, Afrique du Sud (à partir de 2011) ; groupe élargi depuis 2024.',
                'G20 : forum des grandes économies ; l’Union africaine en est membre depuis 2023.',
              ],
            },
            {
              kind: 'list',
              title: 'Les grandes interfaces',
              items: [
                "Façades maritimes : Northern Range en Europe (Rotterdam), façade asiatique (Shanghai, 1er port mondial à conteneurs).",
                'Métropoles mondiales : New York, Londres, Tokyo, Paris, Shanghai.',
              ],
            },
          ],
        },
        {
          title: 'Les grands ensembles régionaux',
          blocks: [
            {
              kind: 'list',
              items: [
                'Union européenne (UE) : 27 États, marché unique, euro.',
                'ACEUM (accord Canada-États-Unis-Mexique), qui a remplacé l’ALENA (1994) en 2020.',
                'ASEAN (Asie du Sud-Est, 1967).',
                'MERCOSUR (Amérique du Sud, 1991).',
                'CEDEAO (Afrique de l’Ouest, 1975) et ZLECAf (Afrique, 2018).',
              ],
            },
            {
              kind: 'definition',
              term: 'Zone de libre-échange',
              definition: 'Espace où les droits de douane entre membres sont supprimés ; chaque pays garde ses tarifs vis-à-vis de l’extérieur.',
            },
            {
              kind: 'definition',
              term: 'Union douanière',
              definition: 'Zone de libre-échange + tarif extérieur commun.',
            },
            {
              kind: 'warning',
              text: "La mondialisation n'efface pas les inégalités : elle intègre fortement certains territoires (centres) et en marginalise d'autres (périphéries).",
            },
            {
              kind: 'tip',
              text: "Commentaire de document : identifie la nature, l'auteur et la date, puis relève les informations en les classant (acteurs, flux, espaces) avant de les expliquer avec tes connaissances.",
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Mondialisation', back: 'Intégration des espaces mondiaux par des flux de marchandises, capitaux, personnes, informations.' },
        { front: 'FTN', back: 'Firme transnationale : maison mère + filiales dans plusieurs pays.' },
        { front: 'OMC', back: 'Organisation mondiale du commerce, créée le 1er janvier 1995, siège à Genève.' },
        { front: 'Bretton Woods (1944)', back: 'Création du FMI et de la Banque mondiale.' },
        { front: 'Triade', back: 'Amérique du Nord, Europe occidentale, Asie orientale.' },
        { front: 'BRICS (membres d’origine)', back: 'Brésil, Russie, Inde, Chine, Afrique du Sud.' },
        { front: 'Délocalisation', back: 'Transfert de production vers des pays à coûts plus faibles.' },
        { front: 'ACEUM', back: "Accord Canada-États-Unis-Mexique, successeur de l'ALENA depuis 2020." },
        { front: 'Union douanière', back: 'Libre-échange entre membres + tarif extérieur commun.' },
      ],
      quiz: [
        {
          id: 'geo-bac-mondialisation-q1',
          type: 'qcm',
          prompt: "Quelle organisation a remplacé le GATT en 1995 ?",
          choices: ["L'OMC", 'Le FMI', 'La Banque mondiale', "L'OCDE"],
          answer: 0,
          explanation: "L'OMC est créée le 1er janvier 1995 ; elle fixe les règles du commerce mondial et règle les différends.",
        },
        {
          id: 'geo-bac-mondialisation-q2',
          type: 'trous',
          prompt: 'Les accords de ___ (1944) créent le FMI et la ___.',
          answers: ['Bretton Woods', 'Banque mondiale'],
          bank: ['Bretton Woods', 'Banque mondiale', 'Marrakech', 'OMC', 'Yalta'],
          explanation: "Bretton Woods organise l'économie d'après-guerre. Marrakech (1994) prépare la création de l'OMC.",
        },
        {
          id: 'geo-bac-mondialisation-q3',
          type: 'vrai-faux',
          prompt: 'Une firme transnationale possède des filiales dans plusieurs pays.',
          answer: true,
          explanation: "C'est sa définition : une maison mère et des filiales à l'étranger, qui produisent et vendent dans de nombreux pays.",
        },
        {
          id: 'geo-bac-mondialisation-q4',
          type: 'qcm',
          prompt: 'Quels sont les trois pôles de la Triade ?',
          choices: [
            'Amérique du Nord, Europe occidentale, Asie orientale',
            'Amérique du Sud, Afrique, Océanie',
            'Chine, Inde, Russie',
            'États-Unis, Brésil, Afrique du Sud',
          ],
          answer: 0,
          explanation: "La Triade désigne les trois centres d'impulsion de l'économie mondiale.",
        },
        {
          id: 'geo-bac-mondialisation-q5',
          type: 'vrai-faux',
          prompt: "L'ALENA existe toujours sous ce nom.",
          answer: false,
          explanation: "Depuis 2020, l'ALENA est remplacé par l'ACEUM (accord Canada-États-Unis-Mexique).",
        },
        {
          id: 'geo-bac-mondialisation-q6',
          type: 'trous',
          prompt: "Une union douanière est une zone de ___ dotée d'un tarif ___ commun.",
          answers: ['libre-échange', 'extérieur'],
          bank: ['libre-échange', 'extérieur', 'intérieur', 'monnaie unique', 'protectionnisme'],
          explanation: "L'UEMOA et l'UE sont des unions douanières : pas de droits de douane entre membres et un tarif commun face à l'extérieur.",
        },
        {
          id: 'geo-bac-mondialisation-q7',
          type: 'qcm',
          prompt: 'Quel pays a rejoint le groupe des BRIC, qui est alors devenu les BRICS ?',
          choices: ['Le Nigeria', 'Le Mexique', "L'Afrique du Sud", "L'Égypte"],
          answer: 2,
          explanation: "L'Afrique du Sud a rejoint le groupe (à partir de 2011, invitée fin 2010) ; le « S » vient de South Africa.",
        },
        {
          id: 'geo-bac-mondialisation-q8',
          type: 'vrai-faux',
          prompt: 'Le Forum social mondial, rendez-vous des altermondialistes, s’est tenu à Dakar en 2011.',
          answer: true,
          explanation: "Le FSM s'est tenu à Dakar en février 2011.",
        },
        {
          id: 'geo-bac-mondialisation-q9',
          type: 'qcm',
          prompt: 'Quel facteur a fortement fait baisser le coût du transport maritime ?',
          choices: ['La conteneurisation', 'Le protectionnisme', 'La fermeture des canaux', 'La hausse des droits de douane'],
          answer: 0,
          explanation: 'Le conteneur standardisé et les porte-conteneurs géants ont rendu le transport de marchandises très bon marché.',
        },
      ],
    },

    // ───────────────────────── Chapitre 2 ─────────────────────────
    {
      id: 'geo-bac-etats-unis',
      title: 'Les États-Unis, première puissance mondiale',
      summary:
        "Les États-Unis dominent le monde par leur économie, leur armée, leur monnaie et leur culture, malgré des faiblesses internes et la concurrence chinoise.",
      essentials: [
        "1re puissance économique (PIB), militaire et financière ; le dollar est la principale monnaie mondiale.",
        "Puissance complète : hard power (armée, économie) et soft power (culture, universités, numérique).",
        'Territoire immense (environ 9,8 millions de km², 50 États) et riche en ressources.',
        "Organisation de l'espace : Nord-Est (Mégalopolis), Sun Belt (Sud et Ouest dynamiques), littoraux et frontière mexicaine.",
        'Faiblesses : inégalités, dette et déficit commercial, concurrence de la Chine.',
      ],
      sections: [
        {
          title: 'Les fondements de la puissance',
          blocks: [
            {
              kind: 'list',
              items: [
                'Territoire d’environ 9,8 millions de km², ouvert sur deux océans.',
                'Ressources : pétrole, gaz (de schiste), charbon, terres agricoles.',
                'Population de plus de 330 millions d’habitants, enrichie par l’immigration.',
                "Agriculture très productive (agrobusiness) : l'un des premiers exportateurs agricoles mondiaux.",
              ],
            },
            {
              kind: 'definition',
              term: 'Agrobusiness',
              definition: "Ensemble des activités liées à l'agriculture, de l'amont (engrais, machines) à l'aval (transformation, distribution), contrôlé par de grandes firmes.",
            },
          ],
        },
        {
          title: 'Hard power et soft power',
          blocks: [
            {
              kind: 'definition',
              term: 'Hard power',
              definition: 'Puissance fondée sur la contrainte : force militaire et poids économique.',
            },
            {
              kind: 'definition',
              term: 'Soft power',
              definition: "Puissance d'influence et d'attraction (culture, valeurs, universités). Notion de Joseph Nye.",
            },
            {
              kind: 'list',
              title: 'Manifestations',
              items: [
                'Premier budget militaire du monde, bases sur tous les continents, OTAN.',
                'Le dollar, monnaie de réserve et d’échange internationale ; Wall Street.',
                'FTN : GAFAM (Google, Apple, Facebook/Meta, Amazon, Microsoft), Coca-Cola…',
                'Hollywood, musique, universités (Harvard, MIT), langue anglaise.',
                "Sièges de l'ONU (New York), du FMI et de la Banque mondiale (Washington).",
              ],
            },
          ],
        },
        {
          title: "L'organisation du territoire",
          blocks: [
            {
              kind: 'definition',
              term: 'Mégalopolis (BosWash)',
              definition: "Chapelet de métropoles du Nord-Est, de Boston à Washington en passant par New York : cœur politique et financier.",
            },
            {
              kind: 'definition',
              term: 'Manufacturing Belt',
              definition: "Ancienne région industrielle du Nord-Est et des Grands Lacs, touchée par la désindustrialisation (« Rust Belt »).",
            },
            {
              kind: 'definition',
              term: 'Sun Belt',
              definition: 'États du Sud et de l’Ouest (Californie, Texas, Floride…) : dynamisme démographique et économique, hautes technologies.',
            },
            {
              kind: 'example',
              title: 'Silicon Valley',
              text: 'Technopôle californien près de San Francisco : siège de grandes firmes du numérique.',
            },
          ],
        },
        {
          title: 'Limites et défis',
          blocks: [
            {
              kind: 'list',
              items: [
                'Inégalités sociales et tensions raciales.',
                'Déficit commercial et forte dette publique.',
                'Concurrence de la Chine.',
                'Contestation de leur rôle de « gendarme du monde » (Irak 2003, Afghanistan).',
              ],
            },
            { kind: 'date', date: '11 septembre 2001', event: 'Attentats terroristes contre le World Trade Center et le Pentagone.' },
            {
              kind: 'warning',
              text: "Ne dis pas que les États-Unis sont « en déclin » : ils restent la première puissance mondiale, mais leur domination est relative et contestée.",
            },
            {
              kind: 'tip',
              text: 'Plan type : I. Les fondements de la puissance ; II. Une puissance mondiale (hard et soft power) ; III. Des limites et des défis.',
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Hard power', back: 'Puissance par la contrainte : armée, économie.' },
        { front: 'Soft power', back: 'Puissance par l’influence et l’attraction : culture, valeurs, universités.' },
        { front: 'Mégalopolis (BosWash)', back: 'Chapelet de métropoles du Nord-Est, de Boston à Washington.' },
        { front: 'Sun Belt', back: 'États dynamiques du Sud et de l’Ouest.' },
        { front: 'Silicon Valley', back: 'Technopôle de Californie, cœur du numérique.' },
        { front: 'Manufacturing Belt / Rust Belt', back: 'Vieille région industrielle du Nord-Est et des Grands Lacs.' },
        { front: 'Superficie des États-Unis', back: 'Environ 9,8 millions de km² ; 50 États.' },
        { front: 'Siège de l’ONU / du FMI', back: 'New York / Washington.' },
      ],
      quiz: [
        {
          id: 'geo-bac-etats-unis-q1',
          type: 'qcm',
          prompt: 'Que désigne le soft power ?',
          choices: [
            'La puissance militaire',
            "La puissance d'influence et d'attraction",
            'La puissance agricole',
            "La puissance nucléaire",
          ],
          answer: 1,
          explanation: 'Le soft power (Joseph Nye) repose sur la culture, les valeurs, les universités, le cinéma.',
        },
        {
          id: 'geo-bac-etats-unis-q2',
          type: 'trous',
          prompt: 'La ___ est le chapelet de métropoles du Nord-Est ; la ___ désigne les États dynamiques du Sud et de l’Ouest.',
          answers: ['Mégalopolis', 'Sun Belt'],
          bank: ['Mégalopolis', 'Sun Belt', 'Rust Belt', 'Corn Belt', 'Silicon Valley'],
          explanation: "La Mégalopolis va de Boston à Washington ; la Sun Belt regroupe la Californie, le Texas, la Floride…",
        },
        {
          id: 'geo-bac-etats-unis-q3',
          type: 'vrai-faux',
          prompt: 'Le siège de l’ONU se trouve à Washington.',
          answer: false,
          explanation: "Le siège de l'ONU est à New York ; Washington accueille le FMI et la Banque mondiale.",
        },
        {
          id: 'geo-bac-etats-unis-q4',
          type: 'qcm',
          prompt: 'Où se situe la Silicon Valley ?',
          choices: ['Au Texas', 'En Floride', 'En Californie', 'Dans le Michigan'],
          answer: 2,
          explanation: 'La Silicon Valley se trouve en Californie, au sud de San Francisco.',
        },
        {
          id: 'geo-bac-etats-unis-q5',
          type: 'vrai-faux',
          prompt: 'Le dollar est la principale monnaie des échanges et des réserves internationales.',
          answer: true,
          explanation: "C'est un pilier de la puissance financière américaine.",
        },
        {
          id: 'geo-bac-etats-unis-q6',
          type: 'qcm',
          prompt: 'Combien d’États composent les États-Unis ?',
          choices: ['13', '48', '50', '52'],
          answer: 2,
          explanation: "Les États-Unis sont une fédération de 50 États (13 à l'indépendance en 1776).",
        },
        {
          id: 'geo-bac-etats-unis-q7',
          type: 'trous',
          prompt: 'La puissance par la contrainte est le ___ ; la notion de soft power a été proposée par Joseph ___.',
          answers: ['hard power', 'Nye'],
          bank: ['hard power', 'Nye', 'smart power', 'Monroe', 'Truman'],
          explanation: 'Joseph Nye a théorisé le soft power ; le hard power désigne la force militaire et économique.',
        },
        {
          id: 'geo-bac-etats-unis-q8',
          type: 'vrai-faux',
          prompt: 'La Manufacturing Belt a connu une forte désindustrialisation.',
          answer: true,
          explanation: "D'où son surnom de « Rust Belt » (ceinture de la rouille).",
        },
        {
          id: 'geo-bac-etats-unis-q9',
          type: 'qcm',
          prompt: "Lequel de ces éléments est une FAIBLESSE des États-Unis ?",
          choices: ['Les universités prestigieuses', 'La puissance militaire', 'Le dollar', 'Le déficit commercial'],
          answer: 3,
          explanation: 'Les États-Unis importent beaucoup plus qu’ils n’exportent : leur balance commerciale est déficitaire.',
        },
      ],
    },

    // ───────────────────────── Chapitre 3 ─────────────────────────
    {
      id: 'geo-bac-japon-chine',
      title: 'Les puissances asiatiques : le Japon et la Chine',
      summary:
        "Le Japon, puissance économique malgré un territoire contraignant, et la Chine, devenue la 2e économie mondiale, font de l'Asie orientale un pôle majeur.",
      essentials: [
        "Japon : archipel montagneux exposé aux risques (séismes, tsunamis, volcans, typhons), population d'environ 125 millions vieillissante.",
        'Japon : « miracle » économique (1950-1973), puissance industrielle et technologique, littoral Pacifique concentrant hommes et activités.',
        "Chine : environ 1,4 milliard d'habitants, 2e économie mondiale depuis 2010, « atelier du monde ».",
        "Chine : fortes inégalités entre le littoral et l'intérieur ; ZES, nouvelles routes de la soie.",
        "Les deux pays sont très présents en Afrique (TICAD pour le Japon, FOCAC pour la Chine).",
      ],
      sections: [
        {
          title: 'Le Japon : un territoire contraignant',
          blocks: [
            {
              kind: 'list',
              items: [
                'Archipel de quatre grandes îles : Hokkaido, Honshu, Shikoku, Kyushu.',
                'Relief montagneux sur environ les trois quarts du territoire : les hommes se concentrent sur les plaines littorales.',
                'Risques majeurs : séismes, tsunamis, volcans, typhons.',
                'Peu de matières premières : dépendance aux importations.',
              ],
            },
            { kind: 'date', date: '11 mars 2011', event: 'Séisme et tsunami dans le nord-est du Japon ; accident nucléaire de Fukushima.' },
            {
              kind: 'definition',
              term: 'Mégalopole japonaise',
              definition: 'Ensemble urbain continu sur la côte Pacifique, de Tokyo à Osaka et au-delà ; il concentre population et activités.',
            },
          ],
        },
        {
          title: 'Le Japon : une puissance économique',
          blocks: [
            { kind: 'date', date: '1950-1973', event: "« Miracle japonais » : croissance très forte, grâce à l'État (MITI), aux grands groupes et à l'exportation." },
            { kind: 'date', date: 'Années 1990', event: "« Décennie perdue » après l'éclatement d'une bulle financière." },
            {
              kind: 'list',
              items: [
                'Grandes firmes : Toyota, Sony, Honda, Panasonic…',
                'Toyotisme : production « juste-à-temps », qualité.',
                'Construction de terre-pleins sur la mer pour les zones industrialo-portuaires.',
                'Défis : vieillissement et baisse de la population, dette publique.',
              ],
            },
            {
              kind: 'example',
              title: 'Le Japon et l’Afrique',
              text: 'Depuis 1993, la TICAD (Conférence internationale de Tokyo sur le développement de l’Afrique) organise la coopération nippo-africaine.',
            },
          ],
        },
        {
          title: 'La Chine : une puissance émergée',
          blocks: [
            {
              kind: 'list',
              items: [
                "Environ 1,4 milliard d'habitants ; dépassée par l'Inde comme pays le plus peuplé vers 2023.",
                '2e économie mondiale (PIB) depuis 2010 ; premier exportateur mondial de marchandises.',
                "« Atelier du monde » : main-d'œuvre nombreuse, ZES, IDE.",
                'Politique de l’enfant unique (1979-2015) : vieillissement accéléré.',
              ],
            },
            {
              kind: 'definition',
              term: 'Mingong',
              definition: "Travailleurs migrants venus des campagnes de l'intérieur vers les villes du littoral, souvent sans les mêmes droits que les citadins.",
            },
            {
              kind: 'definition',
              term: 'Nouvelles routes de la soie',
              definition: "Projet lancé en 2013 par Xi Jinping : réseaux de ports, routes, voies ferrées financés par la Chine en Asie, en Europe et en Afrique.",
            },
          ],
        },
        {
          title: 'Les déséquilibres chinois',
          blocks: [
            {
              kind: 'list',
              items: [
                "Littoral riche (Shanghai, Shenzhen, Pékin) contre intérieur et Ouest plus pauvres.",
                'Inégalités villes / campagnes.',
                'Pollution : la Chine est le premier émetteur mondial de CO₂.',
                'Régime autoritaire ; tensions avec Taïwan et les États-Unis.',
              ],
            },
            {
              kind: 'warning',
              text: 'Ne confonds pas : le Japon est une puissance « ancienne » de la Triade ; la Chine est une puissance émergente devenue rivale des États-Unis.',
            },
            {
              kind: 'tip',
              text: "Croquis de la Chine : oppose clairement une Chine littorale (centres moteurs, ZES, grands ports) et une Chine intérieure (périphéries), et figure les flux de mingong.",
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Les 4 grandes îles du Japon', back: 'Hokkaido, Honshu, Shikoku, Kyushu.' },
        { front: 'Risques naturels au Japon', back: 'Séismes, tsunamis, volcans, typhons.' },
        { front: 'Fukushima', back: '11 mars 2011 : séisme, tsunami et accident nucléaire.' },
        { front: 'Miracle japonais', back: 'Croissance très forte de 1950 à 1973.' },
        { front: 'TICAD', back: 'Conférence de Tokyo sur le développement de l’Afrique (depuis 1993).' },
        { front: 'Atelier du monde', back: 'Surnom de la Chine, premier exportateur mondial de produits manufacturés.' },
        { front: 'Mingong', back: 'Migrants ruraux chinois travaillant dans les villes du littoral.' },
        { front: 'Nouvelles routes de la soie', back: 'Projet chinois d’infrastructures lancé en 2013.' },
      ],
      quiz: [
        {
          id: 'geo-bac-japon-chine-q1',
          type: 'qcm',
          prompt: 'Quelle est la plus grande île du Japon, où se trouvent Tokyo et Osaka ?',
          choices: ['Hokkaido', 'Kyushu', 'Honshu', 'Shikoku'],
          answer: 2,
          explanation: 'Honshu, l’île principale, accueille Tokyo, Osaka et la plus grande partie de la mégalopole.',
        },
        {
          id: 'geo-bac-japon-chine-q2',
          type: 'vrai-faux',
          prompt: 'Le Japon est riche en matières premières.',
          answer: false,
          explanation: 'Le Japon importe la plupart de son énergie et de ses matières premières : c’est une contrainte majeure.',
        },
        {
          id: 'geo-bac-japon-chine-q3',
          type: 'trous',
          prompt: 'Le « miracle japonais » va de 1950 à ___ ; les années 1990 sont appelées la « décennie ___ ».',
          answers: ['1973', 'perdue'],
          bank: ['1973', 'perdue', '1945', 'dorée', '2011'],
          explanation: 'Le premier choc pétrolier (1973) ralentit la croissance ; la crise des années 1990 vient de l’éclatement d’une bulle financière.',
        },
        {
          id: 'geo-bac-japon-chine-q4',
          type: 'qcm',
          prompt: 'Depuis quand la Chine est-elle la 2e économie mondiale ?',
          choices: ['1978', '2001', '2010', '2023'],
          answer: 2,
          explanation: "En 2010, le PIB chinois dépasse celui du Japon : la Chine devient la 2e économie mondiale.",
        },
        {
          id: 'geo-bac-japon-chine-q5',
          type: 'vrai-faux',
          prompt: "En Chine, les régions littorales sont plus riches que l'intérieur.",
          answer: true,
          explanation: "Le littoral concentre les ZES, les grands ports et les IDE ; l'intérieur et l'Ouest restent en retrait.",
        },
        {
          id: 'geo-bac-japon-chine-q6',
          type: 'trous',
          prompt: 'Les travailleurs migrants chinois sont appelés ___ ; le projet des nouvelles routes de la soie est lancé en ___.',
          answers: ['mingong', '2013'],
          bank: ['mingong', '2013', 'keiretsu', '1978', '2001'],
          explanation: 'Les mingong quittent les campagnes pour les villes côtières. Les routes de la soie sont lancées par Xi Jinping en 2013.',
        },
        {
          id: 'geo-bac-japon-chine-q7',
          type: 'qcm',
          prompt: 'Que désigne la TICAD ?',
          choices: [
            'La conférence de Tokyo sur le développement de l’Afrique',
            'Une zone économique spéciale chinoise',
            'Une firme japonaise',
            'Un traité entre la Chine et le Japon',
          ],
          answer: 0,
          explanation: 'Créée en 1993, la TICAD est le cadre de la coopération entre le Japon et l’Afrique.',
        },
        {
          id: 'geo-bac-japon-chine-q8',
          type: 'vrai-faux',
          prompt: "La population japonaise vieillit et diminue.",
          answer: true,
          explanation: 'Natalité faible et longue espérance de vie : la population baisse depuis la fin des années 2000.',
        },
        {
          id: 'geo-bac-japon-chine-q9',
          type: 'qcm',
          prompt: 'Quel pays a dépassé la Chine comme pays le plus peuplé du monde vers 2023 ?',
          choices: ["L'Inde", 'Les États-Unis', "L'Indonésie", 'Le Nigeria'],
          answer: 0,
          explanation: "Selon les estimations de l'ONU, l'Inde est devenue le pays le plus peuplé du monde en 2023.",
        },
      ],
    },

    // ───────────────────────── Chapitre 4 ─────────────────────────
    {
      id: 'geo-bac-union-europeenne',
      title: "L'Union européenne",
      summary:
        "Née de la construction européenne après 1945, l'UE est une grande puissance économique et commerciale, mais une puissance politique incomplète.",
      essentials: [
        'Étapes : CECA (1951), traité de Rome (1957, CEE à 6), traité de Maastricht (1992, UE), euro (1999-2002).',
        'Aujourd’hui 27 États membres, après la sortie du Royaume-Uni (Brexit, 2020).',
        'Puissance économique : marché unique d’environ 450 millions d’habitants, grand pôle commercial ; l’euro.',
        "Espace inégal : cœur riche (« dorsale européenne ») et périphéries.",
        "Liens avec l'Afrique : accords de Lomé (1975), Cotonou (2000), puis Samoa (2023).",
      ],
      sections: [
        {
          title: 'La construction européenne',
          blocks: [
            { kind: 'date', date: '1951', event: 'Traité de Paris : création de la CECA (charbon et acier) par 6 pays.' },
            { kind: 'date', date: '25 mars 1957', event: 'Traité de Rome : création de la CEE (marché commun) par la France, la RFA, l’Italie, la Belgique, les Pays-Bas et le Luxembourg.' },
            { kind: 'date', date: '1992', event: "Traité de Maastricht : naissance de l'Union européenne (en vigueur en 1993)." },
            { kind: 'date', date: '1999 / 2002', event: "Création de l'euro (1999) ; pièces et billets en circulation le 1er janvier 2002." },
            { kind: 'date', date: '2004', event: "Grand élargissement à 10 nouveaux pays, surtout d'Europe centrale et orientale." },
            { kind: 'date', date: '31 janvier 2020', event: "Sortie du Royaume-Uni (Brexit) : l'UE passe à 27 membres." },
          ],
        },
        {
          title: 'Les institutions',
          blocks: [
            { kind: 'definition', term: 'Commission européenne', definition: "Propose les lois et veille à l'application des traités (siège à Bruxelles)." },
            { kind: 'definition', term: 'Parlement européen', definition: 'Élu au suffrage universel direct depuis 1979 ; vote les lois avec le Conseil (siège à Strasbourg).' },
            { kind: 'definition', term: "Conseil de l'Union européenne", definition: 'Réunit les ministres des États ; vote les lois avec le Parlement.' },
            { kind: 'definition', term: 'BCE', definition: 'Banque centrale européenne : gère l’euro (siège à Francfort).' },
            {
              kind: 'warning',
              text: "Ne confonds pas le Conseil de l'Union européenne (ministres), le Conseil européen (chefs d'État et de gouvernement) et le Conseil de l'Europe (organisation distincte de l'UE).",
            },
          ],
        },
        {
          title: 'Une puissance économique',
          blocks: [
            {
              kind: 'list',
              items: [
                'Marché unique : libre circulation des marchandises, des capitaux, des services et des personnes.',
                'Espace Schengen : suppression des contrôles aux frontières intérieures.',
                "L'euro, deuxième monnaie de réserve du monde après le dollar.",
                'PAC (politique agricole commune, 1962) : puissance agricole.',
                'Grandes firmes et coopérations (Airbus, Ariane).',
              ],
            },
            {
              kind: 'definition',
              term: 'Dorsale européenne',
              definition: "Axe le plus riche et le plus peuplé de l'Europe, du sud de l'Angleterre au nord de l'Italie, en passant par le Benelux et la vallée du Rhin.",
            },
            { kind: 'example', title: 'Rotterdam', text: 'Premier port européen, au cœur de la façade maritime de la mer du Nord (Northern Range).' },
          ],
        },
        {
          title: "Limites et relations avec l'Afrique",
          blocks: [
            {
              kind: 'list',
              title: 'Limites',
              items: [
                "Pas d'armée commune ni de politique étrangère vraiment unifiée.",
                'Inégalités entre États et régions ; vieillissement.',
                'Euroscepticisme (Brexit), crise migratoire.',
              ],
            },
            {
              kind: 'list',
              title: "UE et Afrique",
              items: [
                'Convention de Lomé (1975) avec les pays ACP (Afrique, Caraïbes, Pacifique).',
                'Accord de Cotonou (2000), remplacé par l’accord de Samoa (2023).',
                'APE (accords de partenariat économique) : libéralisation des échanges, débattue en Afrique.',
                "L'UE est un partenaire commercial majeur et un grand bailleur d'aide au développement.",
              ],
            },
            {
              kind: 'tip',
              text: "Pour un sujet « L'UE, une puissance ? », montre ses atouts économiques (I), puis ses limites politiques et internes (II).",
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Traité de Rome', back: '25 mars 1957 : création de la CEE par 6 pays.' },
        { front: 'Les 6 pays fondateurs', back: 'France, RFA, Italie, Belgique, Pays-Bas, Luxembourg.' },
        { front: 'Traité de Maastricht', back: "1992 : naissance de l'Union européenne." },
        { front: "Mise en circulation de l'euro", back: 'Pièces et billets le 1er janvier 2002 (monnaie créée en 1999).' },
        { front: 'Nombre de membres de l’UE', back: '27 depuis le Brexit (2020).' },
        { front: 'Espace Schengen', back: 'Suppression des contrôles aux frontières intérieures.' },
        { front: 'Dorsale européenne', back: 'Axe le plus riche, du sud de l’Angleterre au nord de l’Italie.' },
        { front: 'Accords UE-ACP', back: 'Lomé (1975), Cotonou (2000), Samoa (2023).' },
      ],
      quiz: [
        {
          id: 'geo-bac-union-europeenne-q1',
          type: 'qcm',
          prompt: 'Quel traité crée la CEE en 1957 ?',
          choices: ['Le traité de Paris', 'Le traité de Maastricht', 'Le traité de Rome', 'Le traité de Lisbonne'],
          answer: 2,
          explanation: 'Le traité de Rome (25 mars 1957) crée la Communauté économique européenne. Le traité de Paris (1951) crée la CECA.',
        },
        {
          id: 'geo-bac-union-europeenne-q2',
          type: 'vrai-faux',
          prompt: "Le Royaume-Uni est toujours membre de l'Union européenne.",
          answer: false,
          explanation: "Le Royaume-Uni a quitté l'UE le 31 janvier 2020 (Brexit). L'UE compte 27 membres.",
        },
        {
          id: 'geo-bac-union-europeenne-q3',
          type: 'trous',
          prompt: "Le traité de ___ (1992) crée l'Union européenne ; l'euro est mis en circulation sous forme de billets en ___.",
          answers: ['Maastricht', '2002'],
          bank: ['Maastricht', '2002', 'Rome', '1999', 'Schengen'],
          explanation: "L'euro est créé en 1999 pour les banques et marchés ; pièces et billets arrivent le 1er janvier 2002.",
        },
        {
          id: 'geo-bac-union-europeenne-q4',
          type: 'qcm',
          prompt: 'Quelle institution propose les lois européennes ?',
          choices: ['Le Parlement européen', 'La BCE', 'La Commission européenne', 'La Cour de justice'],
          answer: 2,
          explanation: "La Commission a l'initiative des lois ; le Parlement et le Conseil de l'UE les votent.",
        },
        {
          id: 'geo-bac-union-europeenne-q5',
          type: 'vrai-faux',
          prompt: 'Le Parlement européen est élu au suffrage universel direct depuis 1979.',
          answer: true,
          explanation: 'Les citoyens européens élisent leurs députés tous les cinq ans depuis 1979.',
        },
        {
          id: 'geo-bac-union-europeenne-q6',
          type: 'qcm',
          prompt: "Quel accord a remplacé l'accord de Cotonou entre l'UE et les pays ACP ?",
          choices: ['La convention de Lomé', "L'accord de Samoa", 'Le traité de Lisbonne', "L'accord de Marrakech"],
          answer: 1,
          explanation: "Signé en 2023, l'accord de Samoa succède à l'accord de Cotonou (2000), lui-même successeur des conventions de Lomé.",
        },
        {
          id: 'geo-bac-union-europeenne-q7',
          type: 'trous',
          prompt: "La BCE, qui gère l'euro, siège à ___ ; la Commission européenne siège à ___.",
          answers: ['Francfort', 'Bruxelles'],
          bank: ['Francfort', 'Bruxelles', 'Strasbourg', 'Luxembourg', 'Paris'],
          explanation: "Strasbourg accueille le Parlement européen ; Luxembourg la Cour de justice de l'UE.",
        },
        {
          id: 'geo-bac-union-europeenne-q8',
          type: 'vrai-faux',
          prompt: "L'Union européenne dispose d'une armée commune.",
          answer: false,
          explanation: "L'UE n'a pas d'armée commune : c'est l'une de ses limites en tant que puissance politique.",
        },
        {
          id: 'geo-bac-union-europeenne-q9',
          type: 'qcm',
          prompt: 'Quel est le premier port européen ?',
          choices: ['Marseille', 'Le Havre', 'Hambourg', 'Rotterdam'],
          answer: 3,
          explanation: 'Rotterdam (Pays-Bas), sur la mer du Nord, est le premier port européen.',
        },
      ],
    },

    // ───────────────────────── Chapitre 5 ─────────────────────────
    {
      id: 'geo-bac-afrique-developpement',
      title: "Les défis du développement en Afrique",
      summary:
        "L'Afrique, continent jeune aux grandes ressources, fait face à de grands défis : démographie, urbanisation, pauvreté et intégration ; le Nigeria et l'Afrique du Sud en sont des exemples.",
      essentials: [
        "Population d'environ 1,4 milliard d'habitants (années 2020), la plus jeune et à la plus forte croissance du monde.",
        'Urbanisation rapide : mégapoles comme Lagos, Kinshasa, Le Caire ; quartiers précaires.',
        'Économies dépendantes des matières premières ; dette ; faible part du commerce mondial.',
        "L'intégration régionale (CEDEAO, SADC, ZLECAf) est une réponse aux défis.",
        "Nigeria (pays le plus peuplé, pétrole) et Afrique du Sud (puissance industrielle, BRICS) sont des puissances régionales.",
      ],
      sections: [
        {
          title: 'Le défi démographique',
          blocks: [
            {
              kind: 'list',
              items: [
                "Environ 1,4 milliard d'habitants au début des années 2020.",
                'Transition démographique inachevée : la mortalité a baissé, la natalité reste élevée.',
                "Population très jeune : grands besoins en éducation, santé, emplois.",
              ],
            },
            {
              kind: 'definition',
              term: 'Dividende démographique',
              definition:
                "Bonus de croissance possible quand la part des personnes en âge de travailler augmente, à condition qu'elles soient formées et employées.",
            },
            {
              kind: 'definition',
              term: 'IDH (indice de développement humain)',
              definition: 'Indicateur du PNUD combinant santé (espérance de vie), éducation et revenu par habitant.',
            },
          ],
        },
        {
          title: "Le défi urbain et économique",
          blocks: [
            {
              kind: 'list',
              title: 'Urbanisation',
              items: [
                'Croissance urbaine très rapide (exode rural + croissance naturelle).',
                'Mégapoles : Lagos, Kinshasa, Le Caire.',
                'Habitat précaire, manque d’eau, de transports, d’emplois formels (secteur informel).',
              ],
            },
            {
              kind: 'list',
              title: 'Économie',
              items: [
                'Économies de rente : exportation de pétrole, minerais, produits agricoles bruts.',
                'Faible industrialisation et faible transformation sur place.',
                "L'Afrique pèse moins de 3 % du commerce mondial.",
                'Dette, instabilité politique, conflits.',
              ],
            },
            {
              kind: 'definition',
              term: 'Économie de rente',
              definition: "Économie dépendant des revenus de l'exportation d'une ou de quelques matières premières, sensible aux variations des cours.",
            },
          ],
        },
        {
          title: 'Les réponses : intégration et projets',
          blocks: [
            {
              kind: 'list',
              items: [
                'Organisations régionales : CEDEAO, UEMOA, SADC (Afrique australe), EAC (Afrique de l’Est).',
                'ZLECAf : zone de libre-échange continentale (accord signé à Kigali en 2018, échanges lancés en 2021).',
                'Agenda 2063 de l’Union africaine.',
                'Grande Muraille verte : projet de reboisement à travers le Sahel, du Sénégal à Djibouti.',
              ],
            },
            {
              kind: 'example',
              title: 'Nouveaux partenaires',
              text: "Chine (FOCAC), Inde, Turquie, pays du Golfe s'ajoutent aux partenaires traditionnels (UE, États-Unis).",
            },
          ],
        },
        {
          title: 'Deux puissances régionales : Nigeria et Afrique du Sud',
          blocks: [
            {
              kind: 'list',
              title: 'Nigeria',
              items: [
                "Pays le plus peuplé d'Afrique (plus de 200 millions d'habitants).",
                'Fédération de 36 États ; capitale Abuja (depuis 1991), métropole économique Lagos.',
                'Grand producteur de pétrole (delta du Niger), membre de l’OPEP.',
                'Nollywood, grande industrie du cinéma.',
                'Défis : dépendance au pétrole, pauvreté, insécurité (Boko Haram), pollution du delta.',
              ],
            },
            {
              kind: 'list',
              title: 'Afrique du Sud',
              items: [
                'Puissance industrielle et minière (or, platine, diamants, charbon).',
                'Membre des BRICS et seul pays africain membre du G20 à titre national.',
                'Trois capitales : Pretoria (exécutif), Le Cap (législatif), Bloemfontein (judiciaire).',
                "Héritage de l'apartheid (1948-1991) : fortes inégalités, chômage.",
              ],
            },
            {
              kind: 'warning',
              text: "La capitale du Nigeria est Abuja, pas Lagos ; la capitale politique de l'Afrique du Sud est Pretoria, pas Johannesburg.",
            },
            {
              kind: 'tip',
              text: "Pour une étude de cas, suis le plan : atouts de puissance (I), limites et défis (II), rôle régional et mondial (III).",
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Population de l’Afrique (années 2020)', back: 'Environ 1,4 milliard d’habitants.' },
        { front: 'Dividende démographique', back: 'Bonus de croissance quand la population active augmente et trouve un emploi.' },
        { front: 'IDH', back: 'Indice du PNUD : santé, éducation, revenu.' },
        { front: 'Économie de rente', back: 'Dépendance à l’exportation de matières premières.' },
        { front: 'Capitale du Nigeria', back: 'Abuja (depuis 1991) ; Lagos est la métropole économique.' },
        { front: "Capitales de l'Afrique du Sud", back: 'Pretoria (exécutif), Le Cap (législatif), Bloemfontein (judiciaire).' },
        { front: 'Grande Muraille verte', back: 'Projet de reboisement à travers le Sahel, du Sénégal à Djibouti.' },
        { front: 'SADC', back: "Communauté de développement de l'Afrique australe." },
      ],
      quiz: [
        {
          id: 'geo-bac-afrique-developpement-q1',
          type: 'qcm',
          prompt: "Quel est le pays le plus peuplé d'Afrique ?",
          choices: ["L'Égypte", "L'Éthiopie", 'Le Nigeria', "L'Afrique du Sud"],
          answer: 2,
          explanation: "Le Nigeria compte plus de 200 millions d'habitants, devant l'Éthiopie et l'Égypte.",
        },
        {
          id: 'geo-bac-afrique-developpement-q2',
          type: 'vrai-faux',
          prompt: "La capitale du Nigeria est Lagos.",
          answer: false,
          explanation: 'Depuis 1991, la capitale est Abuja. Lagos reste la plus grande ville et le cœur économique.',
        },
        {
          id: 'geo-bac-afrique-developpement-q3',
          type: 'trous',
          prompt: "L'___ mesure le développement par la santé, l'éducation et le revenu ; une économie dépendant des matières premières est une économie de ___.",
          answers: ['IDH', 'rente'],
          bank: ['IDH', 'rente', 'PIB', 'marché', 'services'],
          explanation: "L'IDH est calculé par le PNUD. L'économie de rente est fragile face aux variations des cours.",
        },
        {
          id: 'geo-bac-afrique-developpement-q4',
          type: 'qcm',
          prompt: "Lequel de ces pays africains est membre des BRICS depuis 2011 ?",
          choices: ["L'Afrique du Sud", 'Le Nigeria', 'Le Kenya', 'Le Sénégal'],
          answer: 0,
          explanation: "L'Afrique du Sud a rejoint le groupe en 2011 (le « S » de BRICS).",
        },
        {
          id: 'geo-bac-afrique-developpement-q5',
          type: 'vrai-faux',
          prompt: "L'Afrique a la population la plus jeune du monde.",
          answer: true,
          explanation: 'La natalité y reste élevée : une grande partie de la population a moins de 20 ans.',
        },
        {
          id: 'geo-bac-afrique-developpement-q6',
          type: 'trous',
          prompt: 'Le Nigeria est un grand producteur de ___, extrait surtout dans le delta du ___.',
          answers: ['pétrole', 'Niger'],
          bank: ['pétrole', 'Niger', 'cuivre', 'Congo', 'Nil'],
          explanation: "Le pétrole du delta du Niger fournit l'essentiel des exportations nigérianes, mais pollue la région.",
        },
        {
          id: 'geo-bac-afrique-developpement-q7',
          type: 'qcm',
          prompt: "Quelle ville est la capitale exécutive de l'Afrique du Sud ?",
          choices: ['Johannesburg', 'Le Cap', 'Pretoria', 'Durban'],
          answer: 2,
          explanation: 'Pretoria est la capitale exécutive ; Le Cap la capitale législative ; Bloemfontein la capitale judiciaire.',
        },
        {
          id: 'geo-bac-afrique-developpement-q8',
          type: 'vrai-faux',
          prompt: 'La Grande Muraille verte est un projet de reboisement à travers le Sahel.',
          answer: true,
          explanation: "Ce projet panafricain, du Sénégal à Djibouti, vise à lutter contre la désertification.",
        },
        {
          id: 'geo-bac-afrique-developpement-q9',
          type: 'qcm',
          prompt: "Que signifie « dividende démographique » ?",
          choices: [
            'La baisse de la population',
            "L'argent envoyé par les émigrés",
            "Le bonus de croissance lié à l'augmentation de la population active, si elle est employée",
            'Une taxe sur les naissances',
          ],
          answer: 2,
          explanation: "Ce bonus n'est pas automatique : il faut former et employer les jeunes.",
        },
      ],
    },
  ],
};

export default subject;
