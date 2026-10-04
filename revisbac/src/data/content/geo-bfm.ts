import type { Subject } from '../types';

const subject: Subject = {
  id: 'geo-bfm',
  name: 'Géographie',
  icon: '🌍',
  color: '#0891B2',
  tracks: ['bfm'],
  chapters: [
    // ───────────────────────── Chapitre 1 ─────────────────────────
    {
      id: 'geo-bfm-milieu-physique',
      title: 'Le milieu physique du Sénégal',
      summary:
        "Le Sénégal est un pays plat au climat tropical, arrosé de plus en plus du nord au sud, drainé par quatre grands cours d'eau.",
      essentials: [
        "Superficie : 196 722 km². Le relief est bas et plat ; les seules hauteurs, modestes (pas plus de 650 m environ), sont au sud-est.",
        "Climat tropical à deux saisons : saison sèche et saison des pluies (hivernage), de juin-juillet à octobre.",
        'Les pluies augmentent du nord (moins de 400 mm) vers le sud (plus de 1 000 mm en Casamance).',
        'Trois vents : alizé maritime, harmattan (sec et chaud) et mousson (qui apporte la pluie).',
        'Principaux cours d’eau : Sénégal, Gambie, Casamance, Saloum.',
      ],
      sections: [
        {
          title: 'Situation et relief',
          blocks: [
            {
              kind: 'list',
              title: 'Situation',
              items: [
                "Pays le plus à l'ouest de l'Afrique continentale (pointe des Almadies, presqu'île du Cap-Vert).",
                'Superficie : 196 722 km².',
                'Voisins : Mauritanie, Mali, Guinée, Guinée-Bissau ; la Gambie forme une enclave.',
                "Façade maritime d'environ 700 km sur l'océan Atlantique.",
              ],
            },
            {
              kind: 'list',
              title: 'Relief',
              items: [
                'Un pays de plaines et de plateaux bas : la majeure partie du territoire est à moins de 100 m d’altitude.',
                'Au sud-est, les contreforts du Fouta-Djalon forment les seules hauteurs, modestes : elles ne dépassent pas 650 m environ (point culminant près de Népen Diakha, région de Kédougou).',
                "À Dakar, les collines des Mamelles sont d'anciens volcans.",
              ],
            },
          ],
        },
        {
          title: 'Le climat',
          blocks: [
            {
              kind: 'definition',
              term: 'Hivernage',
              definition: 'Nom donné à la saison des pluies au Sénégal (en gros de juin-juillet à octobre).',
            },
            {
              kind: 'list',
              title: 'Les trois vents',
              items: [
                "Alizé maritime : vent frais et humide venant de l'océan, sans pluie ; il souffle sur la côte nord.",
                "Harmattan (alizé continental) : vent chaud et sec venant du Sahara (nord-est), chargé de poussière.",
                'Mousson : vent humide venant du sud-ouest (océan), qui apporte les pluies d’hivernage.',
              ],
            },
            {
              kind: 'definition',
              term: 'Isohyète',
              definition: 'Ligne qui relie, sur une carte, les points recevant la même quantité de pluie par an.',
            },
            {
              kind: 'list',
              title: 'Les domaines climatiques (du nord au sud)',
              items: [
                'Sahélien : au nord, pluies faibles (moins de 500 mm).',
                'Soudanien : au centre et au sud-est.',
                'Subguinéen : en Casamance, le plus arrosé.',
                'Domaine côtier (subcanarien) : sur la Grande-Côte, températures adoucies par l’alizé maritime.',
              ],
            },
            {
              kind: 'warning',
              text: "Les pluies augmentent du NORD vers le SUD, pas d'ouest en est. Le nord est le plus sec.",
            },
          ],
        },
        {
          title: "Les cours d'eau",
          blocks: [
            {
              kind: 'list',
              items: [
                'Le fleuve Sénégal (plus de 1 700 km) : né en Guinée (Fouta-Djalon), il sert de frontière avec la Mauritanie et se jette dans l’océan près de Saint-Louis.',
                "La Gambie : née aussi au Fouta-Djalon, elle traverse le Sénégal oriental puis la Gambie.",
                'La Casamance : fleuve du sud.',
                "Le Saloum : un estuaire (bras de mer) où l'eau salée remonte loin à l'intérieur.",
              ],
            },
            {
              kind: 'definition',
              term: 'Lac de Guiers',
              definition: "Lac alimenté par le fleuve Sénégal ; il fournit une grande partie de l'eau potable de Dakar.",
            },
            {
              kind: 'date',
              date: '1986 et 1988',
              event: 'Mise en service des barrages de Diama (anti-sel, près de l’embouchure) et de Manantali (au Mali, électricité et régulation).',
            },
          ],
        },
        {
          title: 'Végétation, sols et grandes zones',
          blocks: [
            {
              kind: 'list',
              title: 'Végétation (du nord au sud)',
              items: [
                'Steppe à épineux au nord (Ferlo).',
                'Savane au centre.',
                'Forêt claire puis forêt plus dense en Casamance.',
                "Mangrove dans les estuaires (Saloum, Casamance).",
              ],
            },
            {
              kind: 'list',
              title: 'Sols',
              items: [
                "Sols « dior » : sableux, favorables à l'arachide (bassin arachidier).",
                'Sols « deck » : plus argileux, plus riches.',
                'Sols « hollaldé » : argileux, dans la vallée du fleuve Sénégal.',
              ],
            },
            {
              kind: 'example',
              title: 'Les Niayes',
              text: "Bande côtière entre Dakar et Saint-Louis, avec des cuvettes humides : c'est la grande zone de maraîchage (légumes).",
            },
            {
              kind: 'tip',
              text: 'Sur un croquis, utilise des couleurs du plus clair (nord sec) au plus foncé (sud humide) pour représenter la pluviométrie.',
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Superficie du Sénégal', back: '196 722 km².' },
        { front: 'Hivernage', back: 'La saison des pluies (de juin-juillet à octobre environ).' },
        { front: 'Harmattan', back: 'Vent chaud et sec venant du Sahara (nord-est).' },
        { front: 'Mousson', back: 'Vent humide du sud-ouest qui apporte les pluies.' },
        { front: 'Région la plus arrosée', back: 'La Casamance (sud).' },
        { front: 'Isohyète', back: 'Ligne reliant les points qui reçoivent la même quantité de pluie.' },
        { front: 'Où sont les seules hauteurs du Sénégal ?', back: 'Au sud-est, sur les contreforts du Fouta-Djalon.' },
        { front: 'Les Niayes', back: 'Zone côtière entre Dakar et Saint-Louis, spécialisée dans le maraîchage.' },
        { front: 'Lac de Guiers', back: "Lac alimenté par le fleuve Sénégal, qui approvisionne Dakar en eau." },
      ],
      quiz: [
        {
          id: 'geo-bfm-milieu-physique-q1',
          type: 'qcm',
          prompt: 'Quel vent apporte les pluies au Sénégal ?',
          choices: ["L'harmattan", "L'alizé maritime", 'La mousson'],
          answer: 2,
          explanation: "La mousson, vent humide venu du sud-ouest, apporte les pluies d'hivernage. L'harmattan est sec et l'alizé maritime ne donne pas de pluie.",
        },
        {
          id: 'geo-bfm-milieu-physique-q2',
          type: 'vrai-faux',
          prompt: 'Au Sénégal, les pluies augmentent du nord vers le sud.',
          answer: true,
          explanation: 'Le nord (Podor, Matam) reçoit moins de 400 mm par an, la Basse-Casamance plus de 1 000 mm.',
        },
        {
          id: 'geo-bfm-milieu-physique-q3',
          type: 'trous',
          prompt: "L'___ est un vent chaud et sec venu du Sahara ; la saison des pluies s'appelle l'___.",
          answers: ['harmattan', 'hivernage'],
          bank: ['harmattan', 'hivernage', 'mousson', 'alizé maritime', 'contre-saison'],
          explanation: "L'harmattan (alizé continental) souffle en saison sèche ; l'hivernage est la saison des pluies.",
        },
        {
          id: 'geo-bfm-milieu-physique-q4',
          type: 'qcm',
          prompt: 'Où se trouvent les seules hauteurs du Sénégal ?',
          choices: ['Au sud-est, près de la Guinée', 'Au nord, près du fleuve', 'Sur la côte', 'Au centre, dans le bassin arachidier'],
          answer: 0,
          explanation: 'Les contreforts du Fouta-Djalon, au sud-est (région de Kédougou), forment les seules hauteurs, modestes : elles ne dépassent pas 650 m environ.',
        },
        {
          id: 'geo-bfm-milieu-physique-q5',
          type: 'vrai-faux',
          prompt: 'Le fleuve Sénégal prend sa source au Sénégal.',
          answer: false,
          explanation: 'Le fleuve Sénégal naît en Guinée, dans le massif du Fouta-Djalon, puis traverse le Mali avant de longer le Sénégal.',
        },
        {
          id: 'geo-bfm-milieu-physique-q6',
          type: 'qcm',
          prompt: 'Quel barrage empêche la remontée de l’eau salée dans le fleuve Sénégal ?',
          choices: ['Manantali', 'Kandadji', 'Diama', 'Akosombo'],
          answer: 2,
          explanation: "Diama, près de l'embouchure, est un barrage anti-sel. Manantali (au Mali) produit de l'électricité et régule le débit.",
        },
        {
          id: 'geo-bfm-milieu-physique-q7',
          type: 'trous',
          prompt: 'Les sols ___, sableux, sont favorables à l’arachide ; la ___ pousse dans les estuaires du Saloum et de la Casamance.',
          answers: ['dior', 'mangrove'],
          bank: ['dior', 'mangrove', 'hollaldé', 'steppe', 'deck'],
          explanation: 'Les sols dior dominent dans le bassin arachidier. La mangrove (palétuviers) se développe dans les eaux saumâtres.',
        },
        {
          id: 'geo-bfm-milieu-physique-q8',
          type: 'vrai-faux',
          prompt: 'La Gambie forme une enclave à l’intérieur du territoire sénégalais.',
          answer: true,
          explanation: "La Gambie est entourée par le Sénégal, sauf sur sa petite façade atlantique : elle forme une enclave.",
        },
        {
          id: 'geo-bfm-milieu-physique-q9',
          type: 'qcm',
          prompt: 'Quel est le domaine climatique de la Casamance ?',
          choices: ['Sahélien', 'Subguinéen', 'Désertique', 'Subcanarien'],
          answer: 1,
          explanation: 'La Casamance a un climat subguinéen : c’est la région la plus arrosée, avec une végétation de forêt.',
        },
      ],
    },

    // ───────────────────────── Chapitre 2 ─────────────────────────
    {
      id: 'geo-bfm-population',
      title: 'La population du Sénégal',
      summary:
        "La population sénégalaise est jeune, en forte croissance, inégalement répartie et de plus en plus urbaine.",
      essentials: [
        'Environ 18 millions d’habitants au recensement de 2023, avec une croissance de près de 3 % par an.',
        'Une population très jeune : environ la moitié a moins de 20 ans.',
        "Répartition inégale : l'Ouest est très peuplé, l'Est presque vide ; la région de Dakar concentre près d'un quart de la population.",
        "Urbanisation rapide, alimentée par l'exode rural.",
        'Une population diverse : Wolof, Pulaar, Sérère, Diola, Mandingue, Soninké…',
      ],
      sections: [
        {
          title: 'Les notions démographiques',
          blocks: [
            { kind: 'definition', term: 'Taux de natalité', definition: "Nombre de naissances pour 1 000 habitants en un an." },
            { kind: 'definition', term: 'Taux de mortalité', definition: 'Nombre de décès pour 1 000 habitants en un an.' },
            {
              kind: 'formula',
              label: 'Accroissement naturel',
              formula: 'Taux d’accroissement naturel = taux de natalité − taux de mortalité',
              note: 'Exprimé en ‰ (pour mille) ou en %.',
            },
            { kind: 'formula', label: 'Densité', formula: 'Densité = population ÷ superficie (hab/km²)' },
            {
              kind: 'definition',
              term: 'Transition démographique',
              definition:
                'Passage d’une natalité et d’une mortalité fortes à une natalité et une mortalité faibles. Le Sénégal est en cours de transition : la mortalité a beaucoup baissé, la natalité baisse plus lentement.',
            },
          ],
        },
        {
          title: 'Une population nombreuse et jeune',
          blocks: [
            {
              kind: 'list',
              items: [
                'Environ 18 millions d’habitants (recensement général de 2023, ANSD).',
                'Croissance rapide : près de 3 % par an entre 2013 et 2023.',
                'Une pyramide des âges à base large : beaucoup d’enfants et de jeunes.',
              ],
            },
            {
              kind: 'list',
              title: 'Conséquences',
              items: [
                'Besoins importants en écoles, santé, emplois et logements.',
                'Une main-d’œuvre abondante : un atout si elle est formée.',
              ],
            },
            {
              kind: 'warning',
              text: 'Ne confonds pas croissance naturelle (naissances − décès) et croissance totale (qui inclut aussi les migrations).',
            },
          ],
        },
        {
          title: 'Une répartition inégale',
          blocks: [
            {
              kind: 'list',
              title: 'Zones très peuplées',
              items: [
                'La région de Dakar : près d’un quart des habitants sur 0,3 % du territoire.',
                'Thiès, Diourbel (Touba) et le bassin arachidier.',
                'La Basse-Casamance.',
              ],
            },
            {
              kind: 'list',
              title: 'Zones peu peuplées',
              items: [
                'Le Ferlo (zone d’élevage, sèche).',
                "Le Sénégal oriental (Tambacounda, Kédougou).",
              ],
            },
            {
              kind: 'text',
              text: "Densité moyenne : environ 90 hab/km² (2023), mais plusieurs milliers d'hab/km² dans la région de Dakar.",
            },
            {
              kind: 'list',
              title: 'Divisions administratives',
              items: ['14 régions, dont Dakar, Thiès, Diourbel, Saint-Louis, Kaolack, Ziguinchor, Tambacounda, Kédougou…'],
            },
          ],
        },
        {
          title: "L'urbanisation et l'exode rural",
          blocks: [
            {
              kind: 'definition',
              term: 'Exode rural',
              definition: 'Départ des habitants des campagnes vers les villes.',
            },
            {
              kind: 'definition',
              term: 'Urbanisation',
              definition: 'Augmentation de la part de la population qui vit en ville.',
            },
            {
              kind: 'list',
              title: "Causes de l'exode rural",
              items: [
                'Sécheresses (notamment dans les années 1970) et baisse des rendements agricoles.',
                'Pauvreté et manque d’emplois en saison sèche.',
                'Attraction de la ville : emplois, écoles, hôpitaux.',
              ],
            },
            {
              kind: 'list',
              title: 'Conséquences en ville',
              items: [
                'Quartiers spontanés (irréguliers), inondations en hivernage.',
                'Chômage et développement du secteur informel.',
                'Embouteillages, pollution, problèmes d’eau et d’assainissement.',
              ],
            },
            {
              kind: 'tip',
              text: 'Pour un sujet sur l’exode rural, organise ta réponse : causes (répulsion des campagnes / attraction des villes), puis conséquences (pour les villes / pour les campagnes), puis solutions.',
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Population du Sénégal (recensement 2023)', back: 'Environ 18 millions d’habitants.' },
        { front: 'Accroissement naturel', back: 'Taux de natalité − taux de mortalité.' },
        { front: 'Densité', back: 'Population ÷ superficie, en habitants par km².' },
        { front: 'Exode rural', back: 'Départ des ruraux vers les villes.' },
        { front: 'Région la plus peuplée', back: 'Dakar : près d’un quart de la population sur 0,3 % du territoire.' },
        { front: 'Zones peu peuplées', back: 'Le Ferlo et le Sénégal oriental.' },
        { front: 'Nombre de régions administratives', back: '14.' },
        { front: 'Transition démographique', back: 'Passage d’une forte natalité et mortalité à une faible natalité et mortalité.' },
      ],
      quiz: [
        {
          id: 'geo-bfm-population-q1',
          type: 'qcm',
          prompt: 'Comment calcule-t-on le taux d’accroissement naturel ?',
          choices: [
            'Natalité + mortalité',
            'Immigration − émigration',
            'Population ÷ superficie',
            'Natalité − mortalité',
          ],
          answer: 3,
          explanation: "L'accroissement naturel est la différence entre le taux de natalité et le taux de mortalité.",
        },
        {
          id: 'geo-bfm-population-q2',
          type: 'vrai-faux',
          prompt: 'La population du Sénégal est surtout composée de personnes âgées.',
          answer: false,
          explanation: 'Au contraire, elle est très jeune : environ la moitié des Sénégalais a moins de 20 ans.',
        },
        {
          id: 'geo-bfm-population-q3',
          type: 'trous',
          prompt: "Le départ des habitants des campagnes vers les villes s'appelle l'___ ; il accélère l'___.",
          answers: ['exode rural', 'urbanisation'],
          bank: ['exode rural', 'urbanisation', 'émigration', 'transhumance', 'natalité'],
          explanation: "L'exode rural nourrit la croissance des villes, donc l'urbanisation.",
        },
        {
          id: 'geo-bfm-population-q4',
          type: 'qcm',
          prompt: 'Quelle région est la plus densément peuplée ?',
          choices: ['Tambacounda', 'Matam', 'Dakar', 'Kédougou'],
          answer: 2,
          explanation: "La région de Dakar, la plus petite, concentre près d'un quart de la population : sa densité dépasse plusieurs milliers d'hab/km².",
        },
        {
          id: 'geo-bfm-population-q5',
          type: 'vrai-faux',
          prompt: 'Le Ferlo est une zone faiblement peuplée.',
          answer: true,
          explanation: "Le Ferlo, zone sèche du nord-est consacrée à l'élevage, a une faible densité de population.",
        },
        {
          id: 'geo-bfm-population-q6',
          type: 'qcm',
          prompt: 'Combien le Sénégal comptait-il d’habitants environ au recensement de 2023 ?',
          choices: ['8 millions', '12 millions', '18 millions', '30 millions'],
          answer: 2,
          explanation: 'Le 5e recensement général (2023) a dénombré environ 18 millions d’habitants.',
        },
        {
          id: 'geo-bfm-population-q7',
          type: 'trous',
          prompt: 'La densité se calcule en divisant la ___ par la ___.',
          answers: ['population', 'superficie'],
          bank: ['population', 'superficie', 'natalité', 'mortalité'],
          explanation: 'Densité = population ÷ superficie ; elle s’exprime en habitants par km².',
        },
        {
          id: 'geo-bfm-population-q8',
          type: 'vrai-faux',
          prompt: "Les sécheresses des années 1970 ont accéléré l'exode rural.",
          answer: true,
          explanation: 'Les mauvaises récoltes ont poussé de nombreux ruraux à partir vers les villes, surtout Dakar.',
        },
        {
          id: 'geo-bfm-population-q9',
          type: 'qcm',
          prompt: "Laquelle de ces conséquences n'est PAS liée à la croissance rapide des villes ?",
          choices: ['Les quartiers spontanés', 'Les embouteillages', 'Le développement du secteur informel', 'La hausse de la pluviométrie'],
          answer: 3,
          explanation: 'La pluviométrie dépend du climat, pas de l’urbanisation. Les trois autres sont des conséquences typiques de la croissance urbaine.',
        },
      ],
    },

    // ───────────────────────── Chapitre 3 ─────────────────────────
    {
      id: 'geo-bfm-agriculture-peche-elevage',
      title: "L'agriculture, la pêche et l'élevage",
      summary:
        "Le secteur primaire occupe une grande partie des Sénégalais, mais reste fragile car très dépendant des pluies.",
      essentials: [
        "L'agriculture est surtout pluviale : elle dépend de l'hivernage.",
        'Cultures vivrières (mil, sorgho, riz, maïs, niébé) et cultures commerciales (arachide, coton, horticulture).',
        "Le bassin arachidier est au centre-ouest ; la vallée du fleuve Sénégal pratique l'irrigation (riz, tomate, canne à sucre).",
        "Côtes très poissonneuses grâce à l'upwelling ; la pêche artisanale en pirogue domine.",
        "Élevage surtout extensif et transhumant, notamment dans le Ferlo.",
      ],
      sections: [
        {
          title: "L'agriculture",
          blocks: [
            { kind: 'definition', term: 'Culture vivrière', definition: 'Culture destinée à nourrir la famille ou le pays (mil, sorgho, riz, maïs, niébé, manioc).' },
            { kind: 'definition', term: 'Culture commerciale (de rente)', definition: "Culture destinée à la vente, souvent à l'exportation (arachide, coton, fruits et légumes)." },
            { kind: 'definition', term: 'Culture irriguée', definition: "Culture arrosée artificiellement grâce à des aménagements (canaux, pompes)." },
            {
              kind: 'list',
              title: 'Les grandes zones agricoles',
              items: [
                'Bassin arachidier (Kaolack, Kaffrine, Fatick, Diourbel, Thiès, Louga) : arachide et mil.',
                'Vallée du fleuve Sénégal : riz irrigué, tomate industrielle, canne à sucre (Richard-Toll).',
                'Niayes : maraîchage (oignon, chou, carotte…).',
                'Casamance : riz, fruits (mangues, anacarde/cajou).',
                'Sénégal oriental et Haute-Casamance : coton.',
              ],
            },
          ],
        },
        {
          title: "Les problèmes de l'agriculture",
          blocks: [
            {
              kind: 'list',
              items: [
                'Dépendance aux pluies irrégulières, sécheresses.',
                'Appauvrissement et salinisation des sols.',
                'Outils souvent traditionnels, faibles rendements.',
                'Ennemis des cultures (criquets, oiseaux granivores).',
                'Le pays importe beaucoup de riz pour se nourrir.',
              ],
            },
            {
              kind: 'list',
              title: 'Solutions',
              items: [
                'Maîtrise de l’eau : barrages, irrigation (SAED dans la vallée du fleuve).',
                'Semences améliorées, engrais, mécanisation.',
                'Formation des paysans, crédit agricole.',
              ],
            },
            {
              kind: 'warning',
              text: "L'arachide est une culture commerciale, mais elle est aussi consommée localement (huile, pâte) : ne la classe pas uniquement comme culture d'exportation.",
            },
          ],
        },
        {
          title: 'La pêche',
          blocks: [
            {
              kind: 'definition',
              term: 'Upwelling',
              definition:
                'Remontée d’eaux froides profondes, riches en plancton, le long des côtes : elle attire de nombreux poissons.',
            },
            {
              kind: 'list',
              items: [
                'Pêche artisanale en pirogue : la plus importante (Kayar, Joal, Mbour, Saint-Louis).',
                'Pêche industrielle : chalutiers, port de Dakar.',
                'Transformation : conserveries, fumage et séchage du poisson par les femmes.',
              ],
            },
            {
              kind: 'list',
              title: 'Problèmes',
              items: ['Surpêche et raréfaction du poisson.', 'Pêche illégale et chalutiers étrangers.', 'Accidents en mer et émigration irrégulière en pirogue.'],
            },
          ],
        },
        {
          title: "L'élevage",
          blocks: [
            { kind: 'definition', term: 'Élevage extensif', definition: 'Élevage sur de grands espaces, avec peu d’investissements ; les animaux se nourrissent dans les pâturages naturels.' },
            { kind: 'definition', term: 'Transhumance', definition: 'Déplacement saisonnier des troupeaux à la recherche d’eau et de pâturages.' },
            {
              kind: 'list',
              items: [
                'Bovins, ovins, caprins, volaille.',
                'Zone principale : le Ferlo, avec des forages pour abreuver le bétail.',
                'Forte demande en moutons pour la Tabaski.',
                'Problèmes : manque d’eau et de pâturages, maladies, vol de bétail.',
              ],
            },
            {
              kind: 'tip',
              text: "Pour chaque activité, retiens le schéma : localisation → atouts → problèmes → solutions. C'est le plan type d'une réponse au BFM.",
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Culture vivrière', back: 'Culture pour nourrir la population : mil, sorgho, riz, maïs, niébé.' },
        { front: 'Culture de rente', back: 'Culture destinée à la vente : arachide, coton, horticulture.' },
        { front: 'Bassin arachidier', back: 'Centre-ouest : Kaolack, Kaffrine, Fatick, Diourbel, Thiès, Louga.' },
        { front: 'Upwelling', back: 'Remontée d’eaux froides riches en plancton : côtes très poissonneuses.' },
        { front: 'Transhumance', back: 'Déplacement saisonnier des troupeaux vers l’eau et les pâturages.' },
        { front: 'Zone d’élevage principale', back: 'Le Ferlo.' },
        { front: 'Richard-Toll', back: 'Ville de la vallée du fleuve, centre de la culture de la canne à sucre.' },
        { front: 'Grands ports de pêche artisanale', back: 'Kayar, Joal, Mbour, Saint-Louis.' },
      ],
      quiz: [
        {
          id: 'geo-bfm-agriculture-peche-elevage-q1',
          type: 'qcm',
          prompt: 'Laquelle de ces cultures est une culture vivrière ?',
          choices: ['Le mil', 'Le coton', 'La canne à sucre pour l’exportation', 'Le tabac'],
          answer: 0,
          explanation: 'Le mil sert d’abord à nourrir les familles : c’est une culture vivrière.',
        },
        {
          id: 'geo-bfm-agriculture-peche-elevage-q2',
          type: 'vrai-faux',
          prompt: "L'agriculture sénégalaise dépend surtout des pluies.",
          answer: true,
          explanation: "La majorité des cultures sont pluviales : une mauvaise saison des pluies entraîne de mauvaises récoltes.",
        },
        {
          id: 'geo-bfm-agriculture-peche-elevage-q3',
          type: 'trous',
          prompt: "Le phénomène d'___ rend les côtes sénégalaises très poissonneuses ; la pêche ___ en pirogue est la plus importante.",
          answers: ['upwelling', 'artisanale'],
          bank: ['upwelling', 'artisanale', 'industrielle', 'harmattan', 'mousson'],
          explanation: "L'upwelling fait remonter des eaux froides riches en plancton. La pêche artisanale (Kayar, Joal, Mbour) domine.",
        },
        {
          id: 'geo-bfm-agriculture-peche-elevage-q4',
          type: 'qcm',
          prompt: 'Dans quelle zone se trouve surtout le maraîchage ?',
          choices: ['Le Ferlo', 'Le Sénégal oriental', 'Les Niayes', 'Le delta du Saloum'],
          answer: 2,
          explanation: 'Les Niayes, entre Dakar et Saint-Louis, ont des cuvettes humides idéales pour les légumes.',
        },
        {
          id: 'geo-bfm-agriculture-peche-elevage-q5',
          type: 'vrai-faux',
          prompt: 'Le Sénégal produit tout le riz qu’il consomme.',
          answer: false,
          explanation: 'Malgré la riziculture de la vallée du fleuve et de la Casamance, le pays importe encore beaucoup de riz.',
        },
        {
          id: 'geo-bfm-agriculture-peche-elevage-q6',
          type: 'qcm',
          prompt: 'Où cultive-t-on la canne à sucre au Sénégal ?',
          choices: ['À Kaolack', 'À Richard-Toll', 'À Ziguinchor', 'À Kayar'],
          answer: 1,
          explanation: 'La canne à sucre est cultivée en irrigué à Richard-Toll, dans la vallée du fleuve Sénégal.',
        },
        {
          id: 'geo-bfm-agriculture-peche-elevage-q7',
          type: 'trous',
          prompt: "L'élevage est surtout ___ ; les troupeaux pratiquent la ___ pour trouver eau et pâturages.",
          answers: ['extensif', 'transhumance'],
          bank: ['extensif', 'transhumance', 'intensif', 'irrigation', 'sédentarisation'],
          explanation: "L'élevage extensif utilise de grands espaces ; la transhumance est le déplacement saisonnier des troupeaux.",
        },
        {
          id: 'geo-bfm-agriculture-peche-elevage-q8',
          type: 'vrai-faux',
          prompt: 'La surpêche est un problème pour la pêche sénégalaise.',
          answer: true,
          explanation: 'Trop de bateaux (dont des chalutiers étrangers) et la pêche illégale font diminuer les stocks de poissons.',
        },
        {
          id: 'geo-bfm-agriculture-peche-elevage-q9',
          type: 'qcm',
          prompt: 'Quelle région fait partie du bassin arachidier ?',
          choices: ['Kaolack', 'Kédougou', 'Ziguinchor', 'Matam'],
          answer: 0,
          explanation: 'Kaolack est au cœur du bassin arachidier, avec Kaffrine, Fatick, Diourbel, Thiès et Louga.',
        },
      ],
    },

    // ───────────────────────── Chapitre 4 ─────────────────────────
    {
      id: 'geo-bfm-industrie-energie-transports',
      title: "L'industrie, les mines, l'énergie et les transports",
      summary:
        "Concentrée autour de Dakar, l'industrie sénégalaise transforme les produits agricoles et miniers ; l'énergie et les transports se modernisent.",
      essentials: [
        "L'industrie est concentrée dans la région de Dakar et sur l'axe Dakar-Thiès.",
        'Industries agroalimentaires (huileries, sucre, conserveries), chimiques (phosphates) et cimenteries.',
        'Ressources minières : phosphates, or (Sabodala, Kédougou), zircon ; pétrole et gaz exploités depuis 2024-2025.',
        "Énergie surtout thermique, mais développement du solaire et de l'éolien.",
        'Grands projets de transport : autoroute à péage, aéroport Blaise-Diagne (2017), TER (2021).',
      ],
      sections: [
        {
          title: "L'industrie",
          blocks: [
            {
              kind: 'list',
              title: 'Principales branches',
              items: [
                "Agroalimentaire : huileries d'arachide, sucre (Richard-Toll), conserveries de poisson, minoteries.",
                'Chimie : transformation des phosphates en acide phosphorique et engrais (ICS).',
                'Matériaux de construction : cimenteries (Rufisque, Pout…).',
              ],
            },
            {
              kind: 'list',
              title: 'Problèmes',
              items: [
                'Concentration à Dakar : déséquilibre régional.',
                "Coût élevé de l'énergie.",
                'Concurrence des produits importés, faible transformation sur place.',
              ],
            },
            { kind: 'definition', term: 'Pôle urbain de Diamniadio', definition: "Nouvelle ville en construction près de Dakar pour désengorger la capitale (administrations, université, industries)." },
          ],
        },
        {
          title: 'Les mines et les hydrocarbures',
          blocks: [
            {
              kind: 'list',
              items: [
                'Phosphates : Taïba et région de Thiès.',
                'Or : Sabodala, dans la région de Kédougou.',
                'Zircon : sables de la Grande-Côte.',
                'Calcaire et attapulgite pour la construction et l’industrie.',
              ],
            },
            {
              kind: 'list',
              title: 'Pétrole et gaz',
              items: [
                'Gisement pétrolier de Sangomar (au large) : production commencée en 2024.',
                'Gisement gazier GTA (Grand Tortue Ahmeyim), partagé avec la Mauritanie : production démarrée entre fin 2024 et 2025.',
              ],
            },
            {
              kind: 'warning',
              text: "Le gaz de GTA est partagé avec la MAURITANIE, pas avec la Gambie.",
            },
          ],
        },
        {
          title: "L'énergie",
          blocks: [
            {
              kind: 'list',
              items: [
                "Électricité surtout produite par des centrales thermiques (fioul, gaz) ; la SENELEC la distribue.",
                "Hydroélectricité partagée grâce au barrage de Manantali (OMVS).",
                'Énergies renouvelables en essor : centrales solaires, parc éolien de Taïba Ndiaye.',
                "Le bois et le charbon de bois restent utilisés pour la cuisine, ce qui favorise la déforestation.",
              ],
            },
            {
              kind: 'definition',
              term: 'Énergie renouvelable',
              definition: 'Énergie tirée d’une source qui ne s’épuise pas : soleil, vent, eau.',
            },
          ],
        },
        {
          title: 'Les transports',
          blocks: [
            {
              kind: 'list',
              items: [
                'Route : principal mode de transport ; autoroute à péage Dakar-Diamniadio-AIBD.',
                'Rail : TER (Train express régional) entre Dakar et Diamniadio ; ancien chemin de fer Dakar-Niger.',
                'Port autonome de Dakar : principal port du pays, débouché pour le Mali.',
                'Aéroport international Blaise-Diagne (AIBD), à Diass.',
                'BRT (bus rapides sur voie réservée) à Dakar depuis 2024.',
              ],
            },
            { kind: 'date', date: '7 décembre 2017', event: "Ouverture de l'aéroport international Blaise-Diagne (AIBD)." },
            { kind: 'date', date: 'Décembre 2021', event: 'Mise en service du TER entre Dakar et Diamniadio.' },
            {
              kind: 'tip',
              text: 'Sur une carte des transports, représente les axes par des lignes (épaisseur = importance) et les ports et aéroports par des symboles.',
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Où se concentre l’industrie ?', back: 'Dans la région de Dakar et sur l’axe Dakar-Thiès.' },
        { front: 'ICS', back: 'Industries chimiques du Sénégal : transforment les phosphates (acide phosphorique, engrais).' },
        { front: 'Mine d’or de Sabodala', back: 'Dans la région de Kédougou (sud-est).' },
        { front: 'Sangomar', back: 'Gisement de pétrole en mer, en production depuis 2024.' },
        { front: 'GTA (Grand Tortue Ahmeyim)', back: 'Gisement de gaz partagé avec la Mauritanie.' },
        { front: 'AIBD', back: "Aéroport international Blaise-Diagne, ouvert en décembre 2017." },
        { front: 'TER', back: 'Train express régional Dakar-Diamniadio, mis en service fin 2021.' },
        { front: 'Énergies renouvelables au Sénégal', back: 'Solaire, éolien (Taïba Ndiaye), hydroélectricité (Manantali).' },
      ],
      quiz: [
        {
          id: 'geo-bfm-industrie-energie-transports-q1',
          type: 'qcm',
          prompt: 'Où se trouve la mine d’or de Sabodala ?',
          choices: ['Dans la région de Thiès', 'Dans la région de Dakar', 'Dans la région de Saint-Louis', 'Dans la région de Kédougou'],
          answer: 3,
          explanation: 'Sabodala est dans le Sénégal oriental, région de Kédougou, la grande zone de l’or.',
        },
        {
          id: 'geo-bfm-industrie-energie-transports-q2',
          type: 'vrai-faux',
          prompt: "L'industrie sénégalaise est répartie de façon équilibrée sur tout le territoire.",
          answer: false,
          explanation: "Elle est très concentrée dans la région de Dakar et sur l'axe Dakar-Thiès : c'est un déséquilibre régional.",
        },
        {
          id: 'geo-bfm-industrie-energie-transports-q3',
          type: 'trous',
          prompt: "Le gisement gazier GTA est partagé avec la ___ ; le gisement pétrolier de ___ produit depuis 2024.",
          answers: ['Mauritanie', 'Sangomar'],
          bank: ['Mauritanie', 'Sangomar', 'Gambie', 'Taïba', 'Guinée-Bissau'],
          explanation: "GTA est à cheval sur la frontière maritime sénégalo-mauritanienne. Sangomar est un gisement pétrolier au large du Saloum.",
        },
        {
          id: 'geo-bfm-industrie-energie-transports-q4',
          type: 'qcm',
          prompt: 'Que transforment les ICS (Industries chimiques du Sénégal) ?',
          choices: ['Les phosphates', "L'arachide", 'Le poisson', "L'or"],
          answer: 0,
          explanation: 'Les ICS transforment les phosphates en acide phosphorique et en engrais.',
        },
        {
          id: 'geo-bfm-industrie-energie-transports-q5',
          type: 'vrai-faux',
          prompt: "L'aéroport international Blaise-Diagne a ouvert en 2017.",
          answer: true,
          explanation: "L'AIBD, situé à Diass, a ouvert le 7 décembre 2017, remplaçant l'aéroport Léopold-Sédar-Senghor de Yoff.",
        },
        {
          id: 'geo-bfm-industrie-energie-transports-q6',
          type: 'qcm',
          prompt: 'Laquelle de ces sources d’énergie est renouvelable ?',
          choices: ['Le fioul', 'Le gaz naturel', 'Le vent', 'Le charbon'],
          answer: 2,
          explanation: "Le vent (éolien) est inépuisable. Fioul, gaz et charbon sont des énergies fossiles.",
        },
        {
          id: 'geo-bfm-industrie-energie-transports-q7',
          type: 'trous',
          prompt: 'Le ___ relie Dakar à Diamniadio par le rail ; le ___ de Dakar est le principal port du pays.',
          answers: ['TER', 'port autonome'],
          bank: ['TER', 'port autonome', 'BRT', 'AIBD', 'Dakar-Niger'],
          explanation: "Le TER (train express régional) est en service depuis fin 2021. Le port autonome de Dakar sert aussi de débouché au Mali.",
        },
        {
          id: 'geo-bfm-industrie-energie-transports-q8',
          type: 'vrai-faux',
          prompt: "Le barrage de Manantali fournit de l'électricité au Sénégal.",
          answer: true,
          explanation: "Construit au Mali par l'OMVS, Manantali produit de l'hydroélectricité partagée entre le Mali, la Mauritanie et le Sénégal.",
        },
        {
          id: 'geo-bfm-industrie-energie-transports-q9',
          type: 'qcm',
          prompt: 'Quel est le principal mode de transport au Sénégal ?',
          choices: ['Le rail', 'Le transport fluvial', 'La route', "L'avion"],
          answer: 2,
          explanation: "La route transporte l'essentiel des voyageurs et des marchandises.",
        },
      ],
    },

    // ───────────────────────── Chapitre 5 ─────────────────────────
    {
      id: 'geo-bfm-integration-regionale',
      title: "Le Sénégal et l'intégration régionale",
      summary:
        "Le Sénégal coopère avec ses voisins au sein d'organisations régionales (CEDEAO, UEMOA, OMVS, OMVG) pour se développer.",
      essentials: [
        "L'intégration régionale permet un marché plus grand, la libre circulation et des projets communs.",
        "CEDEAO : créée en 1975 (traité de Lagos), siège à Abuja ; libre circulation des personnes.",
        'UEMOA : créée le 10 janvier 1994 à Dakar, siège à Ouagadougou ; 8 pays utilisant le franc CFA, émis par la BCEAO (dont le siège est à Dakar).',
        'OMVS : créée en 1972 (Mali, Mauritanie, Sénégal, rejoints par la Guinée) ; barrages de Diama et Manantali.',
        'Obstacles : conflits, mauvaises routes, faibles échanges, tensions politiques.',
      ],
      sections: [
        {
          title: "Qu'est-ce que l'intégration régionale ?",
          blocks: [
            {
              kind: 'definition',
              term: 'Intégration régionale',
              definition: 'Rapprochement de pays voisins qui mettent en commun leur économie, leurs projets ou leurs politiques.',
            },
            {
              kind: 'list',
              title: 'Avantages',
              items: [
                'Un marché plus vaste pour les entreprises.',
                'Libre circulation des personnes et des marchandises.',
                'Projets communs (barrages, routes, énergie).',
                'Plus de poids face aux grandes puissances.',
              ],
            },
          ],
        },
        {
          title: 'La CEDEAO et l’UEMOA',
          blocks: [
            { kind: 'date', date: '28 mai 1975', event: 'Création de la CEDEAO (Communauté économique des États de l’Afrique de l’Ouest) par le traité de Lagos. Siège : Abuja (Nigeria).' },
            {
              kind: 'list',
              title: 'La CEDEAO',
              items: [
                'Libre circulation des personnes (protocole de 1979) et passeport commun.',
                'Missions de paix (ECOMOG).',
                'Projet de monnaie commune : l’éco.',
                'Le Mali, le Burkina Faso et le Niger, réunis depuis 2023 dans l’Alliance des États du Sahel (AES), l’ont quittée (retrait effectif en janvier 2025).',
              ],
            },
            { kind: 'date', date: '10 janvier 1994', event: "Création de l'UEMOA (Union économique et monétaire ouest-africaine) à Dakar. Siège : Ouagadougou." },
            {
              kind: 'list',
              title: "L'UEMOA",
              items: [
                '8 membres : Bénin, Burkina Faso, Côte d’Ivoire, Guinée-Bissau, Mali, Niger, Sénégal, Togo.',
                'Monnaie commune : le franc CFA, émis par la BCEAO (siège à Dakar).',
                'Union douanière : tarif extérieur commun.',
              ],
            },
            {
              kind: 'warning',
              text: "Ne confonds pas : la CEDEAO (siège Abuja) regroupe des pays ayant des monnaies différentes ; l'UEMOA (siège Ouagadougou) regroupe les pays du franc CFA d'Afrique de l'Ouest.",
            },
          ],
        },
        {
          title: "L'OMVS et l'OMVG",
          blocks: [
            { kind: 'date', date: '11 mars 1972', event: "Création de l'OMVS (Organisation pour la mise en valeur du fleuve Sénégal) par le Mali, la Mauritanie et le Sénégal. La Guinée l'a rejointe en 2006. Siège : Dakar." },
            {
              kind: 'list',
              title: "Réalisations de l'OMVS",
              items: [
                'Barrage de Diama (anti-sel) : favorise l’irrigation.',
                'Barrage de Manantali (au Mali) : électricité et régulation du fleuve.',
                'Objectifs : agriculture irriguée, énergie, navigation.',
              ],
            },
            {
              kind: 'definition',
              term: 'OMVG',
              definition: "Organisation pour la mise en valeur du fleuve Gambie (Gambie, Guinée, Guinée-Bissau, Sénégal)." },
            { kind: 'example', title: 'Pont de la Sénégambie', text: "Inauguré en janvier 2019 sur le fleuve Gambie, il facilite la liaison entre le nord du Sénégal et la Casamance." },
          ],
        },
        {
          title: "Les obstacles à l'intégration",
          blocks: [
            {
              kind: 'list',
              items: [
                'Conflits et instabilité politique (coups d’État, terrorisme au Sahel).',
                'Infrastructures insuffisantes (routes, chemins de fer).',
                'Économies semblables qui échangent peu entre elles.',
                'Tracasseries aux frontières malgré la libre circulation.',
              ],
            },
            {
              kind: 'tip',
              text: "Pour chaque organisation, retiens 4 éléments : date de création, siège, membres, réalisations.",
            },
          ],
        },
      ],
      flashcards: [
        { front: 'CEDEAO : création et siège', back: '1975 (traité de Lagos) ; siège à Abuja.' },
        { front: 'UEMOA : création et siège', back: '10 janvier 1994 à Dakar ; siège à Ouagadougou.' },
        { front: "Nombre de pays de l'UEMOA", back: '8.' },
        { front: 'BCEAO', back: 'Banque centrale des États de l’Afrique de l’Ouest : émet le franc CFA ; siège à Dakar.' },
        { front: 'OMVS', back: 'Créée en 1972 par le Mali, la Mauritanie et le Sénégal (+ Guinée en 2006) ; siège à Dakar.' },
        { front: "Barrages de l'OMVS", back: 'Diama (anti-sel) et Manantali (électricité).' },
        { front: 'AES', back: 'Alliance des États du Sahel : Mali, Burkina Faso, Niger, sortis de la CEDEAO en janvier 2025.' },
        { front: "Monnaie commune prévue par la CEDEAO", back: "L'éco." },
      ],
      quiz: [
        {
          id: 'geo-bfm-integration-regionale-q1',
          type: 'qcm',
          prompt: 'Où se trouve le siège de la CEDEAO ?',
          choices: ['Dakar', 'Ouagadougou', 'Abuja', 'Lagos'],
          answer: 2,
          explanation: "Le siège est à Abuja (Nigeria). La CEDEAO a été créée à Lagos en 1975, mais son siège est à Abuja.",
        },
        {
          id: 'geo-bfm-integration-regionale-q2',
          type: 'trous',
          prompt: "L'UEMOA a été créée en ___ ; sa monnaie, le ___, est émise par la BCEAO.",
          answers: ['1994', 'franc CFA'],
          bank: ['1994', 'franc CFA', '1975', 'éco', 'naira'],
          explanation: "L'UEMOA naît le 10 janvier 1994 à Dakar. Ses 8 membres utilisent le franc CFA.",
        },
        {
          id: 'geo-bfm-integration-regionale-q3',
          type: 'vrai-faux',
          prompt: "L'OMVS a été fondée par le Mali, la Mauritanie et le Sénégal.",
          answer: true,
          explanation: 'L’OMVS est créée en 1972 par ces trois pays ; la Guinée la rejoint en 2006.',
        },
        {
          id: 'geo-bfm-integration-regionale-q4',
          type: 'qcm',
          prompt: "Lequel de ces pays n'est PAS membre de l'UEMOA ?",
          choices: ['Le Togo', 'La Guinée-Bissau', 'Le Ghana', 'Le Bénin'],
          answer: 2,
          explanation: "Le Ghana a sa propre monnaie (le cedi) : il est membre de la CEDEAO mais pas de l'UEMOA.",
        },
        {
          id: 'geo-bfm-integration-regionale-q5',
          type: 'vrai-faux',
          prompt: 'Le siège de la BCEAO se trouve à Abidjan.',
          answer: false,
          explanation: 'Le siège de la BCEAO est à Dakar.',
        },
        {
          id: 'geo-bfm-integration-regionale-q6',
          type: 'trous',
          prompt: "Les deux barrages de l'OMVS sont ___ (anti-sel) et ___ (électricité).",
          answers: ['Diama', 'Manantali'],
          bank: ['Diama', 'Manantali', 'Akosombo', 'Kossou', 'Sélingué'],
          explanation: "Diama, près de l'embouchure, bloque l'eau salée ; Manantali, au Mali, produit de l'électricité.",
        },
        {
          id: 'geo-bfm-integration-regionale-q7',
          type: 'qcm',
          prompt: 'Quel est un obstacle à l’intégration régionale en Afrique de l’Ouest ?',
          choices: ['La libre circulation', 'Le manque d’infrastructures de transport', 'La monnaie commune', 'Le tarif extérieur commun'],
          answer: 1,
          explanation: "Routes et voies ferrées insuffisantes freinent les échanges. Les autres réponses sont des outils de l'intégration.",
        },
        {
          id: 'geo-bfm-integration-regionale-q8',
          type: 'vrai-faux',
          prompt: 'En janvier 2025, le Mali, le Burkina Faso et le Niger ont quitté la CEDEAO.',
          answer: true,
          explanation: "Ces trois pays ont formé l'Alliance des États du Sahel (AES) et leur retrait de la CEDEAO est devenu effectif en janvier 2025.",
        },
        {
          id: 'geo-bfm-integration-regionale-q9',
          type: 'qcm',
          prompt: "Qu'est-ce que l'éco ?",
          choices: ["La monnaie de l'UEMOA aujourd'hui", 'Un barrage sur la Gambie', 'Une banque régionale', 'Le projet de monnaie commune de la CEDEAO'],
          answer: 3,
          explanation: "L'éco est un projet de monnaie unique pour la CEDEAO, pas encore en circulation.",
        },
      ],
    },
  ],
};

export default subject;
