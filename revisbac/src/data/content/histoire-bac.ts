import type { Subject } from '../types';

const subject: Subject = {
  id: 'histoire-bac',
  name: 'Histoire',
  icon: '🏛️',
  color: '#B45309',
  tracks: ['bac-s', 'bac-l'],
  chapters: [
    // ───────────────────────── Chapitre 1 ─────────────────────────
    {
      id: 'histoire-bac-guerre-froide',
      title: 'Les relations internationales depuis 1945 : la Guerre froide',
      summary:
        "De 1947 à 1991, les États-Unis et l'URSS, porteurs de deux modèles opposés, s'affrontent sans guerre directe, entre crises, détente et regain de tensions.",
      essentials: [
        '1947 : doctrine Truman (endiguement) contre doctrine Jdanov (deux camps) : le monde se divise en deux blocs.',
        'Grandes crises : blocus de Berlin (1948-1949), guerre de Corée (1950-1953), mur de Berlin (1961), crise de Cuba (1962).',
        'Coexistence pacifique puis détente (SALT I en 1972, Helsinki en 1975).',
        'Regain de tensions (Afghanistan, 1979), puis réformes de Gorbatchev (1985) et chute du mur (9 novembre 1989).',
        "Fin de la Guerre froide avec la disparition de l'URSS en décembre 1991.",
      ],
      sections: [
        {
          title: 'Deux modèles opposés',
          blocks: [
            {
              kind: 'list',
              title: 'Le modèle américain',
              items: [
                'Démocratie libérale : pluralisme, élections libres, libertés individuelles.',
                "Économie capitaliste (libérale) : propriété privée, libre entreprise, marché.",
                "Mode de vie (« American way of life ») diffusé par le cinéma et la consommation.",
              ],
            },
            {
              kind: 'list',
              title: 'Le modèle soviétique',
              items: [
                'Régime communiste à parti unique (PCUS).',
                'Économie planifiée (plans quinquennaux, Gosplan), collectivisation (kolkhozes, sovkhozes).',
                "Priorité à l'industrie lourde et à l'armement ; contrôle de la société.",
              ],
            },
          ],
        },
        {
          title: 'La naissance des blocs (1945-1949)',
          blocks: [
            { kind: 'date', date: 'Février 1945', event: 'Conférence de Yalta (Roosevelt, Churchill, Staline).' },
            { kind: 'date', date: '5 mars 1946', event: "Discours de Churchill à Fulton : un « rideau de fer » s'est abattu sur l'Europe." },
            { kind: 'date', date: '12 mars 1947', event: 'Doctrine Truman : les États-Unis veulent « endiguer » (containment) le communisme.' },
            { kind: 'date', date: 'Juin 1947', event: "Annonce du plan Marshall : aide économique américaine à l'Europe." },
            { kind: 'date', date: 'Septembre 1947', event: "Doctrine Jdanov : le monde est divisé en un camp « impérialiste » (américain) et un camp « anti-impérialiste » (soviétique). Création du Kominform." },
            { kind: 'date', date: 'Juin 1948 - mai 1949', event: "Blocus de Berlin-Ouest par l'URSS ; pont aérien américain. Première crise de la Guerre froide." },
            { kind: 'date', date: '1949', event: "Création de l'OTAN (avril), de la RFA et de la RDA ; l'URSS a la bombe atomique ; le CAEM (Comecon) est créé." },
            {
              kind: 'definition',
              term: 'Endiguement (containment)',
              definition: "Stratégie américaine visant à empêcher l'extension du communisme dans le monde.",
            },
          ],
        },
        {
          title: 'Les grandes crises',
          blocks: [
            { kind: 'date', date: '1950-1953', event: "Guerre de Corée : la Corée du Nord communiste envahit le Sud ; armistice de Panmunjom (1953), retour au 38e parallèle." },
            { kind: 'date', date: '1955', event: 'Création du pacte de Varsovie.' },
            { kind: 'date', date: '1956', event: "Khrouchtchev prône la « coexistence pacifique » ; l'URSS écrase l'insurrection de Budapest." },
            { kind: 'date', date: '13 août 1961', event: 'Construction du mur de Berlin.' },
            { kind: 'date', date: 'Octobre 1962', event: "Crise de Cuba : l'URSS installe des missiles à Cuba ; blocus américain ; Khrouchtchev retire les missiles. Le monde frôle la guerre nucléaire." },
            { kind: 'date', date: '1963', event: 'Installation du « téléphone rouge » entre Washington et Moscou.' },
            { kind: 'date', date: '1964-1973', event: "Engagement militaire massif des États-Unis au Vietnam ; accords de Paris (1973), puis prise de Saïgon par les communistes (1975)." },
            {
              kind: 'definition',
              term: 'Équilibre de la terreur',
              definition: "Situation où chaque superpuissance possède assez d'armes nucléaires pour détruire l'autre : cela dissuade d'une guerre directe.",
            },
          ],
        },
        {
          title: 'De la détente à la fin de la Guerre froide',
          blocks: [
            { kind: 'date', date: '1972', event: 'Accords SALT I : limitation des armements stratégiques.' },
            { kind: 'date', date: '1975', event: "Acte final d'Helsinki (CSCE) : respect des frontières et des droits de l'homme en Europe." },
            { kind: 'date', date: 'Décembre 1979', event: "Intervention soviétique en Afghanistan : fin de la détente (« guerre fraîche »)." },
            { kind: 'date', date: '1985', event: 'Gorbatchev arrive au pouvoir : perestroïka (restructuration) et glasnost (transparence).' },
            { kind: 'date', date: '9 novembre 1989', event: 'Chute du mur de Berlin.' },
            { kind: 'date', date: '3 octobre 1990', event: "Réunification de l'Allemagne." },
            { kind: 'date', date: 'Décembre 1991', event: "Disparition de l'URSS (démission de Gorbatchev le 25 décembre) : fin de la Guerre froide." },
            {
              kind: 'tip',
              text: "Composition : un plan chronologique fonctionne bien (1947-1953 : naissance et affrontement ; 1953-1975 : coexistence et détente ; 1975-1991 : tensions puis fin). Illustre chaque partie par une crise datée.",
            },
            {
              kind: 'warning',
              text: "Ne confonds pas le blocus de Berlin (1948-1949) et la construction du mur (1961) : ce sont deux crises différentes.",
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Doctrine Truman', back: '12 mars 1947 : endiguement (containment) du communisme.' },
        { front: 'Doctrine Jdanov', back: 'Septembre 1947 : le monde divisé en deux camps, « impérialiste » et « anti-impérialiste ».' },
        { front: 'Blocus de Berlin', back: 'Juin 1948 - mai 1949 ; pont aérien américain.' },
        { front: 'Guerre de Corée', back: '1950-1953 ; armistice de Panmunjom, frontière au 38e parallèle.' },
        { front: 'Construction du mur de Berlin', back: '13 août 1961.' },
        { front: 'Crise de Cuba', back: 'Octobre 1962 : missiles soviétiques à Cuba ; retrait après négociation.' },
        { front: 'Perestroïka et glasnost', back: 'Réformes de Gorbatchev (à partir de 1985) : restructuration et transparence.' },
        { front: "Fin de l'URSS", back: 'Décembre 1991.' },
        { front: 'OTAN / pacte de Varsovie', back: '1949 / 1955.' },
      ],
      quiz: [
        {
          id: 'histoire-bac-guerre-froide-q1',
          type: 'qcm',
          prompt: 'Que propose la doctrine Truman (1947) ?',
          choices: [
            "L'endiguement du communisme",
            "Le partage de l'Allemagne",
            'La coexistence pacifique',
            'La division du monde en deux camps selon Moscou',
          ],
          answer: 0,
          explanation: "La doctrine Truman (12 mars 1947) vise à « endiguer » l'expansion communiste. La division en deux camps est la thèse de Jdanov.",
        },
        {
          id: 'histoire-bac-guerre-froide-q2',
          type: 'trous',
          prompt: 'En 1947, la doctrine ___ répond à la doctrine Truman ; la même année, le plan ___ aide l’Europe.',
          answers: ['Jdanov', 'Marshall'],
          bank: ['Jdanov', 'Marshall', 'Brejnev', 'Monroe', 'Kennan'],
          explanation: "Jdanov théorise les deux camps (septembre 1947). Le plan Marshall, annoncé en juin 1947, est refusé par l'URSS pour son bloc.",
        },
        {
          id: 'histoire-bac-guerre-froide-q3',
          type: 'vrai-faux',
          prompt: 'Le mur de Berlin a été construit pendant le blocus de 1948-1949.',
          answer: false,
          explanation: "Le blocus date de 1948-1949 ; le mur est construit le 13 août 1961, lors d'une autre crise berlinoise.",
        },
        {
          id: 'histoire-bac-guerre-froide-q4',
          type: 'qcm',
          prompt: 'Quelle crise met le monde au bord de la guerre nucléaire en octobre 1962 ?',
          choices: ['La crise de Berlin', 'La guerre de Corée', 'La crise de Cuba', "L'invasion de l'Afghanistan"],
          answer: 2,
          explanation: "La découverte de missiles soviétiques à Cuba provoque un blocus américain ; Khrouchtchev finit par les retirer.",
        },
        {
          id: 'histoire-bac-guerre-froide-q5',
          type: 'vrai-faux',
          prompt: "L'armistice de Panmunjom (1953) rétablit la frontière entre les deux Corées autour du 38e parallèle.",
          answer: true,
          explanation: 'Après trois ans de guerre, la Corée reste divisée de part et d’autre du 38e parallèle.',
        },
        {
          id: 'histoire-bac-guerre-froide-q6',
          type: 'trous',
          prompt: "Gorbatchev lance la ___ (restructuration) et la ___ (transparence).",
          answers: ['perestroïka', 'glasnost'],
          bank: ['perestroïka', 'glasnost', 'détente', 'NEP', 'coexistence pacifique'],
          explanation: "Ces réformes, à partir de 1985, visent à sauver le système soviétique mais accélèrent sa chute.",
        },
        {
          id: 'histoire-bac-guerre-froide-q7',
          type: 'qcm',
          prompt: 'Quel événement met fin à la détente en décembre 1979 ?',
          choices: [
            "L'intervention soviétique en Afghanistan",
            'La crise de Cuba',
            "L'acte final d'Helsinki",
            'La chute du mur de Berlin',
          ],
          answer: 0,
          explanation: "L'invasion de l'Afghanistan par l'URSS ouvre la période de la « guerre fraîche ».",
        },
        {
          id: 'histoire-bac-guerre-froide-q8',
          type: 'vrai-faux',
          prompt: "Churchill a employé l'expression « rideau de fer » dans son discours de Fulton en 1946.",
          answer: true,
          explanation: "Le 5 mars 1946, à Fulton (États-Unis), Churchill dénonce le « rideau de fer » qui coupe l'Europe en deux.",
        },
        {
          id: 'histoire-bac-guerre-froide-q9',
          type: 'qcm',
          prompt: 'Quelle date marque la réunification de l’Allemagne ?',
          choices: ['9 novembre 1989', '25 décembre 1991', '3 octobre 1990', '13 août 1961'],
          answer: 2,
          explanation: "Le mur tombe le 9 novembre 1989 ; l'Allemagne est réunifiée le 3 octobre 1990 ; l'URSS disparaît en décembre 1991.",
        },
        {
          id: 'histoire-bac-guerre-froide-q10',
          type: 'qcm',
          prompt: "Qu'est-ce qui caractérise l'économie soviétique ?",
          choices: ['Le libre marché', 'La planification centralisée', 'La propriété privée des usines', 'La concurrence entre entreprises'],
          answer: 1,
          explanation: "L'économie soviétique est planifiée par l'État (Gosplan, plans quinquennaux) ; les moyens de production sont collectifs.",
        },
      ],
    },

    // ───────────────────────── Chapitre 2 ─────────────────────────
    {
      id: 'histoire-bac-decolonisation-asie',
      title: 'La décolonisation en Asie : Inde et Indochine',
      summary:
        "Après 1945, l'Asie se décolonise la première : de façon plutôt négociée en Inde (1947), par la guerre en Indochine (1946-1954).",
      essentials: [
        'Causes générales : affaiblissement des métropoles après 1945, Charte de l’ONU, nationalismes, soutien des deux Grands.',
        "Inde : lutte non violente de Gandhi et du Congrès ; indépendance le 15 août 1947 avec partition entre l'Inde et le Pakistan.",
        'Indochine : Hô Chi Minh proclame l’indépendance du Vietnam (2 septembre 1945) ; guerre contre la France (1946-1954).',
        'Défaite française à Diên Biên Phu (7 mai 1954) ; accords de Genève (juillet 1954) : Vietnam coupé au 17e parallèle.',
        "Indonésie : indépendance proclamée par Sukarno en 1945, reconnue par les Pays-Bas en 1949.",
      ],
      sections: [
        {
          title: 'Les causes de la décolonisation',
          blocks: [
            {
              kind: 'definition',
              term: 'Décolonisation',
              definition: "Processus par lequel les colonies accèdent à l'indépendance, par la négociation ou par la lutte armée.",
            },
            {
              kind: 'list',
              items: [
                "Affaiblissement des métropoles européennes après la Seconde Guerre mondiale ; défaites face au Japon en Asie.",
                'Charte de l’Atlantique (1941) et Charte de l’ONU (1945) : droit des peuples à disposer d’eux-mêmes.',
                'Anticolonialisme des États-Unis et de l’URSS.',
                'Montée des nationalismes portés par des élites formées (souvent en Occident).',
              ],
            },
          ],
        },
        {
          title: "L'indépendance de l'Inde",
          blocks: [
            { kind: 'date', date: '1885', event: 'Création du Congrès national indien.' },
            { kind: 'date', date: '1930', event: 'Marche du sel de Gandhi contre le monopole britannique sur le sel.' },
            { kind: 'date', date: '1942', event: 'Campagne « Quit India » (Quittez l’Inde) lancée par le Congrès.' },
            { kind: 'date', date: '15 août 1947', event: "Indépendance de l'Inde (Nehru Premier ministre) ; le Pakistan musulman (Jinnah) est créé le 14 août." },
            { kind: 'date', date: '30 janvier 1948', event: 'Assassinat de Gandhi par un extrémiste hindou.' },
            {
              kind: 'definition',
              term: 'Non-violence (ahimsa)',
              definition: 'Méthode de lutte de Gandhi : désobéissance civile, boycott, grèves, sans recours aux armes.',
            },
            {
              kind: 'warning',
              text: "L'indépendance négociée de l'Inde ne fut pas pacifique : la partition provoque des violences entre hindous et musulmans, des centaines de milliers de morts et des millions de réfugiés.",
            },
          ],
        },
        {
          title: "La guerre d'Indochine",
          blocks: [
            { kind: 'date', date: '1941', event: 'Hô Chi Minh fonde le Viet Minh (Ligue pour l’indépendance du Vietnam).' },
            { kind: 'date', date: '2 septembre 1945', event: 'À Hanoï, Hô Chi Minh proclame l’indépendance de la République démocratique du Vietnam.' },
            { kind: 'date', date: 'Fin 1946', event: 'Début de la guerre entre la France et le Viet Minh (bombardement de Haïphong, novembre 1946).' },
            { kind: 'date', date: '7 mai 1954', event: "Chute du camp retranché de Diên Biên Phu : défaite française face au général Giap." },
            { kind: 'date', date: 'Juillet 1954', event: 'Accords de Genève : fin de la guerre ; le Vietnam est divisé au 17e parallèle ; l’indépendance du Laos et du Cambodge est confirmée.' },
            {
              kind: 'text',
              text: "Avec la victoire communiste en Chine (1949), la guerre d'Indochine s'inscrit dans la Guerre froide : les États-Unis financent l'effort français.",
            },
          ],
        },
        {
          title: 'Méthode : comparer deux décolonisations',
          blocks: [
            {
              kind: 'example',
              title: 'Inde / Indochine',
              text: "Inde : décolonisation négociée, leader non violent, puissance coloniale britannique qui cède. Indochine : décolonisation par la guerre, mouvement communiste, puissance française qui refuse.",
            },
            {
              kind: 'tip',
              text: 'Pour comparer, utilise toujours les mêmes critères : acteurs, méthodes, attitude de la métropole, étapes, bilan.',
            },
            { kind: 'date', date: '17 août 1945', event: "Sukarno proclame l'indépendance de l'Indonésie (reconnue par les Pays-Bas en décembre 1949)." },
          ],
        },
      ],
      flashcards: [
        { front: "Indépendance de l'Inde", back: '15 août 1947, avec la partition Inde / Pakistan.' },
        { front: 'Méthode de Gandhi', back: 'La non-violence (ahimsa), la désobéissance civile.' },
        { front: 'Marche du sel', back: '1930 : action non violente de Gandhi contre le monopole britannique du sel.' },
        { front: 'Fondateur du Pakistan', back: 'Muhammad Ali Jinnah (Ligue musulmane).' },
        { front: 'Viet Minh', back: 'Mouvement nationaliste et communiste fondé par Hô Chi Minh en 1941.' },
        { front: 'Diên Biên Phu', back: '7 mai 1954 : défaite décisive de la France en Indochine.' },
        { front: 'Accords de Genève', back: 'Juillet 1954 : fin de la guerre d’Indochine ; Vietnam divisé au 17e parallèle.' },
        { front: "Indépendance de l'Indonésie", back: 'Proclamée par Sukarno le 17 août 1945, reconnue en 1949.' },
      ],
      quiz: [
        {
          id: 'histoire-bac-decolonisation-asie-q1',
          type: 'qcm',
          prompt: "Quand l'Inde devient-elle indépendante ?",
          choices: ['15 août 1947', '2 septembre 1945', '26 janvier 1950', '7 mai 1954'],
          answer: 0,
          explanation: "L'Inde est indépendante le 15 août 1947. Le 26 janvier 1950 est l'entrée en vigueur de sa Constitution (république).",
        },
        {
          id: 'histoire-bac-decolonisation-asie-q2',
          type: 'vrai-faux',
          prompt: 'La partition de l’Inde en 1947 s’est faite sans violence.',
          answer: false,
          explanation: 'Elle provoque des massacres entre communautés et d’immenses déplacements de population.',
        },
        {
          id: 'histoire-bac-decolonisation-asie-q3',
          type: 'trous',
          prompt: "La France est battue à ___ le 7 mai 1954 ; les accords de ___ divisent le Vietnam au 17e parallèle.",
          answers: ['Diên Biên Phu', 'Genève'],
          bank: ['Diên Biên Phu', 'Genève', 'Évian', 'Haïphong', 'Panmunjom'],
          explanation: "Diên Biên Phu met fin à la présence militaire française ; les accords de Genève (juillet 1954) règlent la paix. Évian concerne l'Algérie.",
        },
        {
          id: 'histoire-bac-decolonisation-asie-q4',
          type: 'qcm',
          prompt: 'Qui proclame l’indépendance du Vietnam le 2 septembre 1945 ?',
          choices: ['Le général Giap', 'Bao Dai', 'Hô Chi Minh', 'Sukarno'],
          answer: 2,
          explanation: "Hô Chi Minh, chef du Viet Minh, proclame l'indépendance à Hanoï. Giap est son chef militaire.",
        },
        {
          id: 'histoire-bac-decolonisation-asie-q5',
          type: 'vrai-faux',
          prompt: 'Gandhi a été assassiné en 1948.',
          answer: true,
          explanation: 'Gandhi est assassiné le 30 janvier 1948 par un extrémiste hindou qui lui reprochait sa modération envers les musulmans.',
        },
        {
          id: 'histoire-bac-decolonisation-asie-q6',
          type: 'qcm',
          prompt: 'Quel dirigeant est à l’origine de la création du Pakistan ?',
          choices: ['Nehru', 'Nasser', 'Gandhi', 'Jinnah'],
          answer: 3,
          explanation: 'Muhammad Ali Jinnah, chef de la Ligue musulmane, obtient la création d’un État pour les musulmans.',
        },
        {
          id: 'histoire-bac-decolonisation-asie-q7',
          type: 'trous',
          prompt: 'Gandhi prône la ___ ; en 1930, il organise la marche du ___.',
          answers: ['non-violence', 'sel'],
          bank: ['non-violence', 'sel', 'lutte armée', 'riz', 'coton'],
          explanation: 'La marche du sel (1930) est un acte de désobéissance civile contre le monopole britannique sur le sel.',
        },
        {
          id: 'histoire-bac-decolonisation-asie-q8',
          type: 'vrai-faux',
          prompt: 'La guerre d’Indochine est liée à la Guerre froide.',
          answer: true,
          explanation: "Après 1949, la Chine communiste soutient le Viet Minh et les États-Unis financent largement l'effort de guerre français.",
        },
        {
          id: 'histoire-bac-decolonisation-asie-q9',
          type: 'qcm',
          prompt: "Lequel de ces facteurs n'est PAS une cause de la décolonisation après 1945 ?",
          choices: [
            "L'affaiblissement des métropoles européennes",
            'La Charte des Nations unies',
            'Le soutien des puissances coloniales à l’indépendance immédiate',
            'La montée des nationalismes',
          ],
          answer: 2,
          explanation: "Les puissances coloniales (France surtout) ont souvent refusé l'indépendance ; elle a été arrachée.",
        },
      ],
    },

    // ───────────────────────── Chapitre 3 ─────────────────────────
    {
      id: 'histoire-bac-decolonisation-afrique',
      title: "La décolonisation en Afrique",
      summary:
        "Entre 1951 et 1975, l'Afrique accède à l'indépendance, de façon négociée (Ghana, Afrique noire française) ou par la guerre (Algérie, Kenya, colonies portugaises).",
      essentials: [
        'Ghana (ex-Gold Coast) : Nkrumah obtient l’indépendance le 6 mars 1957, première colonie britannique d’Afrique noire indépendante.',
        'Kenya : révolte des Mau Mau (années 1950), indépendance en 1963 avec Jomo Kenyatta.',
        'Algérie : guerre de 1954 (1er novembre) à 1962 ; accords d’Évian (18 mars 1962), indépendance le 5 juillet 1962.',
        'Afrique noire française : Brazzaville (1944), loi-cadre (1956), référendum (1958), indépendances de 1960.',
        "Colonies portugaises indépendantes après la révolution des Œillets (1974) ; Amilcar Cabral en Guinée-Bissau.",
      ],
      sections: [
        {
          title: 'Afrique britannique : Ghana et Kenya',
          blocks: [
            { kind: 'date', date: '1949', event: "Kwame Nkrumah fonde le CPP (Convention People's Party) en Gold Coast." },
            { kind: 'date', date: '6 mars 1957', event: "La Gold Coast devient indépendante sous le nom de Ghana, avec Nkrumah à sa tête." },
            { kind: 'date', date: '1952', event: "Au Kenya, début de la révolte des Mau Mau ; les Britanniques proclament l'état d'urgence et répriment durement." },
            { kind: 'date', date: '12 décembre 1963', event: 'Indépendance du Kenya ; Jomo Kenyatta, longtemps emprisonné, devient le dirigeant du pays.' },
            {
              kind: 'definition',
              term: 'Indirect rule',
              definition: "Administration britannique s'appuyant sur les chefs traditionnels ; elle a facilité des transitions plus négociées.",
            },
          ],
        },
        {
          title: "La guerre d'Algérie",
          blocks: [
            { kind: 'date', date: '8 mai 1945', event: 'Massacres de Sétif et Guelma : répression sanglante de manifestations nationalistes.' },
            { kind: 'date', date: '1er novembre 1954', event: 'Le FLN (Front de libération nationale) lance l’insurrection (« Toussaint rouge »).' },
            { kind: 'date', date: '1957', event: "Bataille d'Alger ; usage de la torture par l'armée française." },
            { kind: 'date', date: '1958', event: "Crise du 13 mai à Alger ; retour au pouvoir du général de Gaulle. Création du GPRA." },
            { kind: 'date', date: '18 mars 1962', event: "Accords d'Évian entre la France et le FLN (cessez-le-feu le 19 mars)." },
            { kind: 'date', date: '5 juillet 1962', event: "Indépendance de l'Algérie, après le référendum du 1er juillet." },
            {
              kind: 'text',
              text: "La présence d'environ un million d'Européens (« pieds-noirs ») explique le refus français et la violence du conflit.",
            },
            { kind: 'example', title: 'Maroc et Tunisie', text: 'Ces deux protectorats obtiennent leur indépendance de façon plus négociée en 1956.' },
          ],
        },
        {
          title: "L'Afrique noire française et belge",
          blocks: [
            { kind: 'date', date: '1944', event: 'Conférence de Brazzaville : réformes promises, indépendance refusée.' },
            { kind: 'date', date: '1946', event: 'Union française ; création du RDA à Bamako ; abolition du travail forcé.' },
            { kind: 'date', date: '1947', event: 'Insurrection à Madagascar, très durement réprimée.' },
            { kind: 'date', date: '1956', event: 'Loi-cadre Defferre : suffrage universel et autonomie interne.' },
            { kind: 'date', date: '28 septembre 1958', event: 'Référendum sur la Communauté : seule la Guinée de Sékou Touré vote « non » (indépendante le 2 octobre 1958).' },
            { kind: 'date', date: '1960', event: "« Année de l'Afrique » : 17 États africains deviennent indépendants, dont 14 anciennes colonies françaises." },
            { kind: 'date', date: '30 juin 1960', event: 'Indépendance du Congo belge ; Patrice Lumumba Premier ministre (assassiné en janvier 1961).' },
            {
              kind: 'example',
              title: 'Cameroun',
              text: "L'UPC de Ruben Um Nyobè mène une lutte armée contre la France ; le Cameroun devient indépendant le 1er janvier 1960.",
            },
          ],
        },
        {
          title: 'Les colonies portugaises et la fin des dominations blanches',
          blocks: [
            { kind: 'date', date: '20 janvier 1973', event: 'Assassinat à Conakry d’Amilcar Cabral, chef du PAIGC (Guinée-Bissau et Cap-Vert).' },
            { kind: 'date', date: '25 avril 1974', event: 'Révolution des Œillets au Portugal : elle ouvre la voie aux indépendances.' },
            { kind: 'date', date: '1974-1975', event: 'Indépendance de la Guinée-Bissau (reconnue en 1974), du Mozambique, du Cap-Vert et de l’Angola (1975).' },
            { kind: 'date', date: '1980 et 1990', event: 'Indépendance du Zimbabwe (1980), puis de la Namibie (1990).' },
            { kind: 'date', date: '1994', event: "Fin de l'apartheid : Nelson Mandela élu président de l'Afrique du Sud." },
            {
              kind: 'tip',
              text: "Pour une composition sur la décolonisation africaine, distingue les décolonisations négociées (Ghana, AOF) et les décolonisations par la guerre (Algérie, Kenya, colonies portugaises).",
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Indépendance du Ghana', back: '6 mars 1957, sous la direction de Kwame Nkrumah.' },
        { front: 'Mau Mau', back: 'Révolte contre la colonisation britannique au Kenya (années 1950).' },
        { front: 'Début de la guerre d’Algérie', back: '1er novembre 1954 (« Toussaint rouge »), lancée par le FLN.' },
        { front: "Accords d'Évian", back: '18 mars 1962 : fin de la guerre d’Algérie.' },
        { front: "Indépendance de l'Algérie", back: '5 juillet 1962.' },
        { front: 'Pays ayant voté « non » en 1958', back: 'La Guinée de Sékou Touré.' },
        { front: 'Amilcar Cabral', back: 'Chef du PAIGC (Guinée-Bissau, Cap-Vert), assassiné le 20 janvier 1973.' },
        { front: 'Révolution des Œillets', back: '25 avril 1974 au Portugal ; elle permet l’indépendance des colonies portugaises.' },
        { front: 'Patrice Lumumba', back: 'Premier ministre du Congo indépendant (1960), assassiné en janvier 1961.' },
      ],
      quiz: [
        {
          id: 'histoire-bac-decolonisation-afrique-q1',
          type: 'qcm',
          prompt: 'Quel était le nom colonial du Ghana ?',
          choices: ['La Gold Coast', 'La Côte d’Ivoire', 'Le Tanganyika', 'La Rhodésie'],
          answer: 0,
          explanation: "La Gold Coast (Côte-de-l'Or) devient le Ghana à l'indépendance, le 6 mars 1957.",
        },
        {
          id: 'histoire-bac-decolonisation-afrique-q2',
          type: 'trous',
          prompt: "La guerre d'Algérie commence le 1er novembre ___ ; elle se termine par les accords d'___ en 1962.",
          answers: ['1954', 'Évian'],
          bank: ['1954', 'Évian', '1945', 'Genève', 'Brazzaville'],
          explanation: "Le FLN lance l'insurrection en 1954 ; les accords d'Évian (18 mars 1962) mènent à l'indépendance (5 juillet 1962).",
        },
        {
          id: 'histoire-bac-decolonisation-afrique-q3',
          type: 'vrai-faux',
          prompt: 'Au Kenya, la révolte des Mau Mau était dirigée contre la colonisation britannique.',
          answer: true,
          explanation: "La révolte, liée à la question des terres, est durement réprimée ; le Kenya devient indépendant en 1963.",
        },
        {
          id: 'histoire-bac-decolonisation-afrique-q4',
          type: 'qcm',
          prompt: 'Quel événement ouvre la voie à l’indépendance des colonies portugaises ?',
          choices: ['La conférence de Bandung', 'Les accords de Genève', 'La révolution des Œillets', 'La création de l’OUA'],
          answer: 2,
          explanation: "La révolution des Œillets (25 avril 1974) renverse la dictature portugaise, qui menait des guerres coloniales.",
        },
        {
          id: 'histoire-bac-decolonisation-afrique-q5',
          type: 'vrai-faux',
          prompt: 'Le Sénégal a voté « non » au référendum de 1958.',
          answer: false,
          explanation: 'Seule la Guinée de Sékou Touré a voté « non ». Le Sénégal a voté « oui » et est entré dans la Communauté.',
        },
        {
          id: 'histoire-bac-decolonisation-afrique-q6',
          type: 'qcm',
          prompt: 'Qui dirigeait le PAIGC ?',
          choices: ['Agostinho Neto', 'Samora Machel', 'Amilcar Cabral', 'Ruben Um Nyobè'],
          answer: 2,
          explanation: "Amilcar Cabral dirigeait le PAIGC (Guinée-Bissau et Cap-Vert). Neto a dirigé le MPLA (Angola), Machel le FRELIMO (Mozambique), Um Nyobè l'UPC (Cameroun).",
        },
        {
          id: 'histoire-bac-decolonisation-afrique-q7',
          type: 'trous',
          prompt: "Le Congo belge devient indépendant le 30 juin 1960 avec ___ comme Premier ministre ; le Kenya devient indépendant en ___.",
          answers: ['Patrice Lumumba', '1963'],
          bank: ['Patrice Lumumba', '1963', 'Mobutu', '1957', 'Kenyatta'],
          explanation: "Lumumba est assassiné en janvier 1961. Le Kenya devient indépendant le 12 décembre 1963.",
        },
        {
          id: 'histoire-bac-decolonisation-afrique-q8',
          type: 'vrai-faux',
          prompt: "En 1960, 17 États africains accèdent à l'indépendance.",
          answer: true,
          explanation: "D'où le nom d'« année de l'Afrique » : 14 anciennes colonies françaises, plus le Nigeria, la Somalie et le Congo belge.",
        },
        {
          id: 'histoire-bac-decolonisation-afrique-q9',
          type: 'qcm',
          prompt: "Pourquoi la France a-t-elle refusé longtemps l'indépendance de l'Algérie ?",
          choices: [
            "Parce que l'Algérie n'avait pas de mouvement nationaliste",
            "Parce qu'elle était considérée comme une partie de la France, avec environ un million d'Européens",
            "Parce que l'ONU s'y opposait",
          ],
          answer: 1,
          explanation: "L'Algérie était divisée en départements français et comptait une forte population européenne (« pieds-noirs »).",
        },
      ],
    },

    // ───────────────────────── Chapitre 4 ─────────────────────────
    {
      id: 'histoire-bac-tiers-monde',
      title: 'Le Tiers-monde et le non-alignement',
      summary:
        "Les pays nouvellement indépendants cherchent à peser dans un monde bipolaire : Bandung (1955) puis le mouvement des non-alignés (1961), et la lutte pour un nouvel ordre économique.",
      essentials: [
        "L'expression « Tiers-monde » est créée par Alfred Sauvy en 1952, par analogie avec le Tiers-État.",
        "Conférence de Bandung (avril 1955) : 29 pays d'Asie et d'Afrique condamnent le colonialisme.",
        'Conférence de Belgrade (1961) : naissance du mouvement des non-alignés (Tito, Nehru, Nasser).',
        'Revendications économiques : CNUCED (1964), Groupe des 77, nouvel ordre économique international (1974).',
        "Diversification du Tiers-monde : pays pétroliers, pays émergents, pays les moins avancés ; crise de la dette dans les années 1980.",
      ],
      sections: [
        {
          title: "La naissance du Tiers-monde",
          blocks: [
            {
              kind: 'definition',
              term: 'Tiers-monde',
              definition:
                "Expression créée en 1952 par Alfred Sauvy pour désigner les pays pauvres, souvent anciennement colonisés, qui n'appartiennent ni au bloc occidental ni au bloc soviétique.",
            },
            {
              kind: 'list',
              title: 'Caractéristiques communes',
              items: [
                'Passé colonial.',
                'Pauvreté, forte croissance démographique.',
                "Économies dépendantes de l'exportation de matières premières.",
              ],
            },
            {
              kind: 'definition',
              term: "Détérioration des termes de l'échange",
              definition:
                "Baisse du prix des matières premières exportées par rapport au prix des produits industriels importés : le pays doit exporter plus pour importer autant.",
            },
          ],
        },
        {
          title: 'De Bandung à Belgrade',
          blocks: [
            {
              kind: 'date',
              date: '18-24 avril 1955',
              event: "Conférence de Bandung (Indonésie) : 29 pays afro-asiatiques (Nehru, Sukarno, Nasser, Zhou Enlai) condamnent le colonialisme et prônent la coexistence pacifique.",
            },
            { kind: 'date', date: '1956', event: 'Nasser nationalise le canal de Suez ; échec de l’expédition franco-britannique : victoire politique du Tiers-monde.' },
            {
              kind: 'date',
              date: 'Septembre 1961',
              event: 'Conférence de Belgrade : naissance du mouvement des non-alignés autour de Tito (Yougoslavie), Nehru (Inde) et Nasser (Égypte).',
            },
            {
              kind: 'definition',
              term: 'Non-alignement',
              definition: "Refus de s'aligner sur l'un des deux blocs de la Guerre froide.",
            },
            {
              kind: 'warning',
              text: "Bandung n'est pas encore le mouvement des non-alignés : celui-ci naît à Belgrade en 1961. Et beaucoup de non-alignés penchaient en fait vers un bloc.",
            },
          ],
        },
        {
          title: 'Les revendications économiques',
          blocks: [
            { kind: 'date', date: '1960', event: "Création de l'OPEP (Organisation des pays exportateurs de pétrole) à Bagdad." },
            { kind: 'date', date: '1964', event: 'Première CNUCED à Genève ; formation du Groupe des 77.' },
            { kind: 'date', date: '1973', event: 'Premier choc pétrolier : les pays de l’OPEP augmentent fortement le prix du pétrole.' },
            { kind: 'date', date: '1974', event: "L'Assemblée générale de l'ONU adopte la déclaration pour un nouvel ordre économique international (NOEI)." },
            {
              kind: 'example',
              title: 'Le slogan de la CNUCED',
              text: '« Trade, not aid » (du commerce, pas de l’aide) : les pays du Sud réclament des prix justes pour leurs produits.',
            },
          ],
        },
        {
          title: "L'éclatement du Tiers-monde",
          blocks: [
            {
              kind: 'list',
              items: [
                'Années 1980 : crise de la dette (à partir du Mexique en 1982).',
                'Programmes d’ajustement structurel (PAS) imposés par le FMI et la Banque mondiale.',
                'Dévaluation de 50 % du franc CFA le 12 janvier 1994.',
                "Divergences : pays pétroliers, nouveaux pays industrialisés d'Asie (« dragons »), pays les moins avancés (PMA).",
              ],
            },
            {
              kind: 'definition',
              term: 'PAS (programme d’ajustement structurel)',
              definition: "Ensemble de mesures (réduction des dépenses publiques, privatisations, libéralisation) exigées en échange de prêts.",
            },
            {
              kind: 'tip',
              text: 'Plan possible : I. L’émergence du Tiers-monde sur la scène internationale (Bandung, Belgrade) ; II. Le combat pour le développement (CNUCED, NOEI) ; III. Les difficultés et la diversification.',
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Qui a inventé l’expression « Tiers-monde » ?', back: 'Alfred Sauvy, en 1952.' },
        { front: 'Conférence de Bandung', back: 'Avril 1955, Indonésie : 29 pays afro-asiatiques condamnent le colonialisme.' },
        { front: 'Conférence de Belgrade', back: '1961 : naissance du mouvement des non-alignés.' },
        { front: 'Principaux fondateurs du non-alignement', back: 'Tito, Nehru, Nasser (avec aussi Sukarno et Nkrumah).' },
        { front: 'CNUCED', back: 'Conférence des Nations unies sur le commerce et le développement (1964).' },
        { front: 'NOEI', back: 'Nouvel ordre économique international, réclamé à l’ONU en 1974.' },
        { front: 'Dévaluation du franc CFA', back: '12 janvier 1994 : 50 %.' },
        { front: "Détérioration des termes de l'échange", back: 'Baisse du prix des matières premières par rapport aux produits industriels.' },
      ],
      quiz: [
        {
          id: 'histoire-bac-tiers-monde-q1',
          type: 'qcm',
          prompt: 'Où s’est tenue la conférence afro-asiatique de 1955 ?',
          choices: ['Belgrade', 'Le Caire', 'Bandung', 'Accra'],
          answer: 2,
          explanation: 'La conférence de Bandung (Indonésie), du 18 au 24 avril 1955, réunit 29 pays d’Asie et d’Afrique.',
        },
        {
          id: 'histoire-bac-tiers-monde-q2',
          type: 'vrai-faux',
          prompt: "L'expression « Tiers-monde » a été créée par Alfred Sauvy.",
          answer: true,
          explanation: 'En 1952, Sauvy compare ces pays au Tiers-État de 1789 : ignorés, exploités, ils veulent être quelque chose.',
        },
        {
          id: 'histoire-bac-tiers-monde-q3',
          type: 'trous',
          prompt: 'Le mouvement des non-alignés naît à ___ en 1961, autour de Tito, Nehru et ___.',
          answers: ['Belgrade', 'Nasser'],
          bank: ['Belgrade', 'Nasser', 'Bandung', 'Khrouchtchev', 'Mao Zedong'],
          explanation: "Bandung (1955) prépare le mouvement, mais celui-ci naît officiellement à Belgrade en septembre 1961.",
        },
        {
          id: 'histoire-bac-tiers-monde-q4',
          type: 'qcm',
          prompt: 'Que signifie « non-alignement » ?',
          choices: [
            "Refuser de s'aligner sur l'un des deux blocs",
            "Refuser d'adhérer à l'ONU",
            'Rejoindre le pacte de Varsovie',
            'Refuser toute aide extérieure',
          ],
          answer: 0,
          explanation: "Les non-alignés veulent rester indépendants des États-Unis et de l'URSS.",
        },
        {
          id: 'histoire-bac-tiers-monde-q5',
          type: 'vrai-faux',
          prompt: "Le franc CFA a été dévalué de 50 % le 12 janvier 1994.",
          answer: true,
          explanation: "Cette dévaluation, dans le contexte de l'ajustement structurel, a fortement renchéri les produits importés.",
        },
        {
          id: 'histoire-bac-tiers-monde-q6',
          type: 'qcm',
          prompt: 'Quel dirigeant nationalise le canal de Suez en 1956 ?',
          choices: ['Nehru', 'Nkrumah', 'Tito', 'Nasser'],
          answer: 3,
          explanation: "Le 26 juillet 1956, Nasser nationalise le canal ; l'intervention franco-britannique échoue sous la pression des deux Grands.",
        },
        {
          id: 'histoire-bac-tiers-monde-q7',
          type: 'trous',
          prompt: 'En 1964, la première ___ se tient à Genève ; en 1974, les pays du Sud réclament un ___.',
          answers: ['CNUCED', 'NOEI'],
          bank: ['CNUCED', 'NOEI', 'OPEP', 'PAS', 'OTAN'],
          explanation: "La CNUCED défend le commerce équitable pour le Sud ; le NOEI vise à réformer les règles de l'économie mondiale.",
        },
        {
          id: 'histoire-bac-tiers-monde-q8',
          type: 'vrai-faux',
          prompt: 'La conférence de Bandung a réuni surtout des pays européens.',
          answer: false,
          explanation: "Bandung réunit 29 pays d'Asie et d'Afrique, d'où le nom de conférence afro-asiatique.",
        },
        {
          id: 'histoire-bac-tiers-monde-q9',
          type: 'qcm',
          prompt: 'Que désigne la « détérioration des termes de l’échange » ?',
          choices: [
            'La baisse du prix des matières premières par rapport aux produits industriels',
            'La hausse des prix des matières premières',
            "L'arrêt du commerce entre le Nord et le Sud",
          ],
          answer: 0,
          explanation: 'Les pays du Sud doivent exporter toujours plus de matières premières pour acheter la même quantité de produits manufacturés.',
        },
      ],
    },

    // ───────────────────────── Chapitre 5 ─────────────────────────
    {
      id: 'histoire-bac-chine',
      title: 'La Chine depuis 1949',
      summary:
        "Communiste depuis 1949, la Chine passe du maoïsme (Grand Bond en avant, Révolution culturelle) aux réformes de Deng Xiaoping, qui en font une grande puissance.",
      essentials: [
        '1er octobre 1949 : Mao Zedong proclame la République populaire de Chine ; les nationalistes se réfugient à Taïwan.',
        'Grand Bond en avant (1958-1960) : échec et famine meurtrière.',
        'Révolution culturelle (1966-1976) : Gardes rouges, culte de Mao, violences.',
        "1971 : la Chine populaire entre à l'ONU ; 1972 : visite de Nixon.",
        'À partir de 1978, Deng Xiaoping lance les réformes : ouverture, zones économiques spéciales, « socialisme de marché ».',
      ],
      sections: [
        {
          title: 'La Chine de Mao (1949-1976)',
          blocks: [
            { kind: 'date', date: '1er octobre 1949', event: 'Mao Zedong proclame à Pékin la République populaire de Chine. Tchang Kaï-chek se replie à Taïwan.' },
            { kind: 'date', date: '1950', event: "Traité d'amitié avec l'URSS ; réforme agraire." },
            { kind: 'date', date: '1958-1960', event: "Grand Bond en avant : communes populaires, industrialisation forcée. Échec et famine faisant des dizaines de millions de morts selon la plupart des estimations." },
            { kind: 'date', date: 'Vers 1960', event: "Rupture sino-soviétique : retrait des experts soviétiques." },
            { kind: 'date', date: '1964', event: 'Première bombe atomique chinoise.' },
            { kind: 'date', date: '1966-1976', event: 'Révolution culturelle : les Gardes rouges, armés du « Petit Livre rouge », s’attaquent aux cadres et aux intellectuels.' },
            { kind: 'date', date: '9 septembre 1976', event: 'Mort de Mao Zedong.' },
            {
              kind: 'definition',
              term: 'Commune populaire',
              definition: 'Grande unité collective créée pendant le Grand Bond en avant, regroupant agriculture, industrie locale, école et milice.',
            },
          ],
        },
        {
          title: 'La Chine dans le monde',
          blocks: [
            { kind: 'date', date: '1950-1953', event: 'Des « volontaires » chinois combattent en Corée contre les forces de l’ONU.' },
            { kind: 'date', date: '1955', event: 'Zhou Enlai participe à la conférence de Bandung.' },
            { kind: 'date', date: '25 octobre 1971', event: "La Chine populaire obtient le siège de la Chine à l'ONU (et au Conseil de sécurité) à la place de Taïwan." },
            { kind: 'date', date: 'Février 1972', event: 'Visite du président américain Nixon en Chine.' },
            {
              kind: 'example',
              title: "La Chine et l'Afrique",
              text: "Dans les années 1970, la Chine construit le chemin de fer Tanzanie-Zambie (Tazara). Depuis 2000, le Forum sur la coopération sino-africaine (FOCAC) renforce ses liens avec l'Afrique.",
            },
          ],
        },
        {
          title: 'Les réformes de Deng Xiaoping',
          blocks: [
            { kind: 'date', date: 'Décembre 1978', event: "Deng Xiaoping lance la politique de réforme et d'ouverture ; les « Quatre modernisations » (agriculture, industrie, défense, sciences et techniques)." },
            { kind: 'date', date: '1980', event: 'Création des premières zones économiques spéciales (ZES), comme Shenzhen.' },
            { kind: 'date', date: '4 juin 1989', event: 'Répression sanglante des manifestants de la place Tian’anmen : pas de démocratisation politique.' },
            { kind: 'date', date: '1997', event: 'Rétrocession de Hong Kong à la Chine (Macao en 1999).' },
            { kind: 'date', date: '2001', event: "Entrée de la Chine à l'OMC." },
            {
              kind: 'definition',
              term: 'Socialisme de marché',
              definition: "Système chinois qui associe économie de marché et pouvoir politique monopolisé par le Parti communiste.",
            },
            {
              kind: 'definition',
              term: 'ZES (zone économique spéciale)',
              definition: 'Zone littorale ouverte aux capitaux étrangers avec des avantages fiscaux.',
            },
          ],
        },
        {
          title: 'Bilan et méthode',
          blocks: [
            {
              kind: 'list',
              items: [
                "Depuis 2010, la Chine est la 2e économie mondiale (PIB).",
                'Fortes inégalités entre littoral et intérieur, villes et campagnes.',
                'Régime autoritaire : Parti unique, contrôle de l’information.',
              ],
            },
            {
              kind: 'warning',
              text: "Les réformes de Deng sont économiques, pas politiques : le Parti communiste garde le monopole du pouvoir.",
            },
            {
              kind: 'tip',
              text: 'Plan chronologique conseillé : I. La Chine de Mao (1949-1976) ; II. La Chine de Deng et de ses successeurs (depuis 1978).',
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Proclamation de la République populaire de Chine', back: '1er octobre 1949, par Mao Zedong.' },
        { front: 'Grand Bond en avant', back: '1958-1960 : communes populaires, échec et famine.' },
        { front: 'Révolution culturelle', back: '1966-1976 : Gardes rouges, Petit Livre rouge, violences.' },
        { front: 'Entrée de la Chine populaire à l’ONU', back: '25 octobre 1971.' },
        { front: 'Les Quatre modernisations', back: 'Agriculture, industrie, défense, sciences et techniques.' },
        { front: 'ZES', back: 'Zones économiques spéciales, ouvertes aux capitaux étrangers (Shenzhen, 1980).' },
        { front: 'Tian’anmen', back: '4 juin 1989 : répression du mouvement démocratique.' },
        { front: 'Entrée de la Chine à l’OMC', back: '2001.' },
      ],
      quiz: [
        {
          id: 'histoire-bac-chine-q1',
          type: 'qcm',
          prompt: 'Qui proclame la République populaire de Chine en 1949 ?',
          choices: ['Tchang Kaï-chek', 'Deng Xiaoping', 'Mao Zedong', 'Zhou Enlai'],
          answer: 2,
          explanation: 'Mao Zedong la proclame le 1er octobre 1949 ; Tchang Kaï-chek, chef nationaliste vaincu, se réfugie à Taïwan.',
        },
        {
          id: 'histoire-bac-chine-q2',
          type: 'vrai-faux',
          prompt: 'Le Grand Bond en avant a été un succès économique.',
          answer: false,
          explanation: "C'est un échec : désorganisation de l'agriculture et famine qui fait des dizaines de millions de victimes selon la plupart des estimations.",
        },
        {
          id: 'histoire-bac-chine-q3',
          type: 'trous',
          prompt: 'La Révolution culturelle dure de 1966 à ___ ; les jeunes ___ y jouent un rôle central.',
          answers: ['1976', 'Gardes rouges'],
          bank: ['1976', 'Gardes rouges', '1989', 'Khmers rouges', '1958'],
          explanation: "Elle s'achève à la mort de Mao (1976). Les Khmers rouges sont cambodgiens.",
        },
        {
          id: 'histoire-bac-chine-q4',
          type: 'qcm',
          prompt: 'Quel dirigeant lance la politique de réforme et d’ouverture à partir de 1978 ?',
          choices: ['Deng Xiaoping', 'Mao Zedong', 'Xi Jinping', 'Lin Biao'],
          answer: 0,
          explanation: "Deng Xiaoping introduit l'économie de marché, les ZES et l'ouverture aux capitaux étrangers.",
        },
        {
          id: 'histoire-bac-chine-q5',
          type: 'vrai-faux',
          prompt: 'En 1971, la Chine populaire remplace Taïwan à l’ONU et au Conseil de sécurité.',
          answer: true,
          explanation: 'Le 25 octobre 1971, l’Assemblée générale attribue le siège de la Chine à la République populaire.',
        },
        {
          id: 'histoire-bac-chine-q6',
          type: 'trous',
          prompt: 'Les premières ___, comme Shenzhen, sont créées en ___.',
          answers: ['zones économiques spéciales', '1980'],
          bank: ['zones économiques spéciales', '1980', 'communes populaires', '1958', '1949'],
          explanation: 'Les ZES, sur le littoral, attirent les capitaux étrangers grâce à des avantages fiscaux.',
        },
        {
          id: 'histoire-bac-chine-q7',
          type: 'qcm',
          prompt: 'Que se passe-t-il place Tian’anmen le 4 juin 1989 ?',
          choices: [
            'La proclamation de la République populaire',
            'La rétrocession de Hong Kong',
            'La répression sanglante d’un mouvement démocratique',
            'Les funérailles de Mao',
          ],
          answer: 2,
          explanation: "L'armée écrase le mouvement des étudiants qui réclamaient plus de libertés : les réformes restent économiques.",
        },
        {
          id: 'histoire-bac-chine-q8',
          type: 'vrai-faux',
          prompt: 'Hong Kong a été rétrocédée à la Chine en 1997.',
          answer: true,
          explanation: 'Colonie britannique, Hong Kong revient à la Chine le 1er juillet 1997 ; Macao (portugaise) en 1999.',
        },
        {
          id: 'histoire-bac-chine-q9',
          type: 'qcm',
          prompt: "En quelle année la Chine entre-t-elle à l'OMC ?",
          choices: ['1971', '1978', '1997', '2001'],
          answer: 3,
          explanation: "La Chine rejoint l'Organisation mondiale du commerce en 2001, ce qui accélère son essor commercial.",
        },
      ],
    },

    // ───────────────────────── Chapitre 6 ─────────────────────────
    {
      id: 'histoire-bac-panafricanisme',
      title: "Le panafricanisme, l'OUA et l'UA",
      summary:
        "Né dans la diaspora noire, le panafricanisme inspire les luttes d'indépendance puis l'unité africaine, avec l'OUA (1963) puis l'UA (2002).",
      essentials: [
        'Le panafricanisme naît dans la diaspora : Sylvester Williams (1900), W. E. B. Du Bois, Marcus Garvey.',
        'Le congrès de Manchester (1945) réunit Nkrumah, Kenyatta et Du Bois : il réclame l’indépendance.',
        'Deux visions en 1961 : groupe de Casablanca (unité rapide, Nkrumah) et groupe de Monrovia (coopération progressive, Senghor).',
        "L'OUA est créée le 25 mai 1963 à Addis-Abeba par 32 États.",
        "L'Union africaine (UA) la remplace en 2002 ; elle compte 55 États membres.",
      ],
      sections: [
        {
          title: 'Les origines du panafricanisme',
          blocks: [
            {
              kind: 'definition',
              term: 'Panafricanisme',
              definition: "Mouvement politique et culturel qui affirme la solidarité des peuples africains et de la diaspora noire et vise l'unité de l'Afrique.",
            },
            { kind: 'date', date: '1900', event: 'Conférence panafricaine de Londres, organisée par Henry Sylvester Williams (Trinidad).' },
            { kind: 'date', date: '1919', event: 'Premier congrès panafricain, à Paris, organisé par W. E. B. Du Bois (avec l’appui de Blaise Diagne).' },
            { kind: 'definition', term: 'Marcus Garvey', definition: "Jamaïcain, fondateur de l'UNIA (1914) ; il prône le retour des Noirs en Afrique (« Back to Africa »)." },
            { kind: 'date', date: 'Octobre 1945', event: 'Ve congrès panafricain de Manchester : Nkrumah, Kenyatta, Du Bois, Padmore réclament l’indépendance.' },
            {
              kind: 'example',
              title: 'La négritude',
              text: "Mouvement culturel fondé dans les années 1930 par Senghor, Césaire et Damas : il revendique la valeur des cultures noires. Cheikh Anta Diop (« Nations nègres et culture », 1954) l'enrichit par l'histoire.",
            },
          ],
        },
        {
          title: 'Unité africaine : deux visions',
          blocks: [
            { kind: 'date', date: 'Avril 1958', event: 'Conférence des États africains indépendants à Accra (Nkrumah).' },
            {
              kind: 'list',
              title: 'Groupe de Casablanca (1961) : « progressistes »',
              items: ['Ghana, Guinée, Mali, Maroc, Égypte…', 'Unité politique rapide, gouvernement continental (Nkrumah : « Africa must unite »).'],
            },
            {
              kind: 'list',
              title: 'Groupe de Monrovia (1961) : « modérés »',
              items: ['Sénégal, Côte d’Ivoire, Nigeria, Liberia…', 'Coopération progressive entre États souverains (Senghor, Houphouët-Boigny).'],
            },
            {
              kind: 'example',
              title: 'Cheikh Anta Diop',
              text: "Dans « Les fondements économiques et culturels d'un État fédéral d'Afrique noire » (première édition en 1960 sous un titre un peu différent, titre actuel depuis 1974), il plaide pour une fédération africaine.",
            },
          ],
        },
        {
          title: "L'OUA (1963-2002)",
          blocks: [
            { kind: 'date', date: '25 mai 1963', event: "Création de l'OUA (Organisation de l'unité africaine) à Addis-Abeba (Éthiopie), par 32 États. Le 25 mai est la Journée de l'Afrique." },
            {
              kind: 'list',
              title: 'Principes',
              items: [
                'Souveraineté et égalité des États.',
                'Non-ingérence dans les affaires intérieures.',
                'Intangibilité des frontières héritées de la colonisation (1964).',
                'Libération des territoires encore colonisés et lutte contre l’apartheid.',
              ],
            },
            {
              kind: 'list',
              title: 'Bilan',
              items: [
                'Réussite : soutien aux mouvements de libération (colonies portugaises, Afrique australe).',
                'Limites : impuissance face aux coups d’État, guerres civiles et dictatures ; « syndicat de chefs d’État ».',
              ],
            },
          ],
        },
        {
          title: "L'Union africaine (depuis 2002)",
          blocks: [
            { kind: 'date', date: '9 septembre 1999', event: "Déclaration de Syrte (Libye) : décision de créer l'Union africaine." },
            { kind: 'date', date: 'Juillet 2002', event: "Lancement de l'Union africaine à Durban (Afrique du Sud). Siège : Addis-Abeba." },
            {
              kind: 'list',
              title: 'Organes et projets',
              items: [
                'Commission de l’UA, Parlement panafricain, Conseil de paix et de sécurité.',
                'NEPAD (2001) et Agenda 2063.',
                'ZLECAf (Zone de libre-échange continentale africaine), signée à Kigali en 2018.',
              ],
            },
            {
              kind: 'definition',
              term: 'Non-indifférence',
              definition: "Principe de l'UA : elle peut intervenir dans un État membre en cas de crimes graves (génocide, crimes de guerre), contrairement à la non-ingérence de l'OUA.",
            },
            {
              kind: 'tip',
              text: "Sujet type : « Le panafricanisme de 1945 à nos jours ». Plan : I. Un mouvement de libération ; II. L'OUA, unité inachevée ; III. L'UA, nouvelles ambitions et limites.",
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Panafricanisme', back: "Mouvement pour la solidarité des peuples noirs et l'unité de l'Afrique." },
        { front: 'Congrès de Manchester', back: 'Octobre 1945 : Nkrumah, Kenyatta, Du Bois réclament l’indépendance.' },
        { front: 'Marcus Garvey', back: 'Jamaïcain, fondateur de l’UNIA, partisan du « Back to Africa ».' },
        { front: 'Création de l’OUA', back: '25 mai 1963 à Addis-Abeba, par 32 États.' },
        { front: 'Groupe de Casablanca vs groupe de Monrovia', back: 'Unité politique rapide (Nkrumah) vs coopération progressive (Senghor).' },
        { front: "Lancement de l'UA", back: 'Juillet 2002 à Durban ; siège à Addis-Abeba.' },
        { front: 'Intangibilité des frontières', back: "Principe de l'OUA : on ne remet pas en cause les frontières coloniales." },
        { front: 'ZLECAf', back: 'Zone de libre-échange continentale africaine, accord signé à Kigali en 2018.' },
      ],
      quiz: [
        {
          id: 'histoire-bac-panafricanisme-q1',
          type: 'qcm',
          prompt: "Où et quand l'OUA a-t-elle été créée ?",
          choices: ['Accra, 1958', 'Addis-Abeba, 1963', 'Durban, 2002', 'Lagos, 1975'],
          answer: 1,
          explanation: "L'OUA naît le 25 mai 1963 à Addis-Abeba. Durban (2002) correspond au lancement de l'UA ; Lagos (1975) à la CEDEAO.",
        },
        {
          id: 'histoire-bac-panafricanisme-q2',
          type: 'vrai-faux',
          prompt: "Le Sénégal de Senghor faisait partie du groupe de Casablanca.",
          answer: false,
          explanation: "Le Sénégal appartenait au groupe de Monrovia, partisan d'une unité progressive. Le groupe de Casablanca rassemblait le Ghana, la Guinée, le Mali…",
        },
        {
          id: 'histoire-bac-panafricanisme-q3',
          type: 'trous',
          prompt: 'Le Ve congrès panafricain se tient à ___ en 1945 ; Kwame ___ y participe.',
          answers: ['Manchester', 'Nkrumah'],
          bank: ['Manchester', 'Nkrumah', 'Londres', 'Paris', 'Garvey'],
          explanation: "Le congrès de Manchester (1945) marque le passage du panafricanisme à la revendication d'indépendance.",
        },
        {
          id: 'histoire-bac-panafricanisme-q4',
          type: 'qcm',
          prompt: "Quel principe de l'OUA interdit de remettre en cause les frontières coloniales ?",
          choices: ['La non-indifférence', 'Le droit d’ingérence', 'La souveraineté populaire', "L'intangibilité des frontières"],
          answer: 3,
          explanation: "Adopté en 1964, ce principe visait à éviter des guerres en chaîne pour redessiner les frontières.",
        },
        {
          id: 'histoire-bac-panafricanisme-q5',
          type: 'vrai-faux',
          prompt: "L'Union africaine a remplacé l'OUA en 2002.",
          answer: true,
          explanation: "Décidée à Syrte (1999), l'UA est lancée à Durban en juillet 2002.",
        },
        {
          id: 'histoire-bac-panafricanisme-q6',
          type: 'qcm',
          prompt: "Qui a fondé l'UNIA et lancé le mot d'ordre « Back to Africa » ?",
          choices: ['Marcus Garvey', 'W. E. B. Du Bois', 'George Padmore', 'Aimé Césaire'],
          answer: 0,
          explanation: 'Le Jamaïcain Marcus Garvey fonde l’UNIA en 1914 et prône le retour des Noirs en Afrique.',
        },
        {
          id: 'histoire-bac-panafricanisme-q7',
          type: 'trous',
          prompt: "L'OUA est fondée par ___ États ; le ___ est la Journée de l'Afrique.",
          answers: ['32', '25 mai'],
          bank: ['32', '25 mai', '55', '4 avril', '1er juillet'],
          explanation: "32 États fondent l'OUA le 25 mai 1963. L'UA compte aujourd'hui 55 membres.",
        },
        {
          id: 'histoire-bac-panafricanisme-q8',
          type: 'vrai-faux',
          prompt: "Le principe de non-indifférence de l'UA permet d'intervenir en cas de génocide ou de crimes de guerre.",
          answer: true,
          explanation: "C'est une rupture avec la stricte non-ingérence de l'OUA, critiquée après le génocide des Tutsi au Rwanda (1994).",
        },
        {
          id: 'histoire-bac-panafricanisme-q9',
          type: 'qcm',
          prompt: "Quel ouvrage de Cheikh Anta Diop plaide pour une fédération africaine ?",
          choices: [
            "« L'Afrique doit s'unir »",
            '« Les Damnés de la terre »',
            "« Les fondements économiques et culturels d'un État fédéral d'Afrique noire »",
            '« Discours sur le colonialisme »',
          ],
          answer: 2,
          explanation: "Première édition en 1960 (sous un titre un peu différent), titre actuel depuis la réédition de 1974. « L'Afrique doit s'unir » est de Nkrumah, « Les Damnés de la terre » de Fanon, le « Discours sur le colonialisme » de Césaire.",
        },
      ],
    },
  ],
};

export default subject;
