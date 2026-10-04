import type { Subject } from '../types';

const subject: Subject = {
  id: 'histoire-bfm',
  name: 'Histoire',
  icon: '🏛️',
  color: '#B45309',
  tracks: ['bfm'],
  chapters: [
    // ───────────────────────── Chapitre 1 ─────────────────────────
    {
      id: 'histoire-bfm-premiere-guerre-mondiale',
      title: 'La Première Guerre mondiale et la révolution russe',
      summary:
        "De 1914 à 1918, une guerre totale ravage l'Europe, mobilise les colonies et provoque la révolution russe de 1917.",
      essentials: [
        "L'attentat de Sarajevo (28 juin 1914) déclenche la guerre entre deux systèmes d'alliances.",
        'La guerre devient une guerre totale : tranchées, Verdun (1916), mobilisation des civils et des colonies.',
        "Près de 200 000 soldats d'AOF (les tirailleurs sénégalais) sont recrutés pour la France ; Blaise Diagne organise le recrutement de 1918.",
        "En 1917, les bolcheviks de Lénine prennent le pouvoir en Russie ; l'URSS naît en 1922.",
        "L'armistice du 11 novembre 1918 et le traité de Versailles (28 juin 1919) mettent fin à la guerre.",
      ],
      sections: [
        {
          title: 'Les causes et le déclenchement',
          blocks: [
            {
              kind: 'list',
              title: 'Deux blocs d’alliances face à face',
              items: [
                'Triple Entente : France, Royaume-Uni, Russie.',
                'Triple Alliance : Allemagne, Autriche-Hongrie, Italie (qui rejoint finalement l’Entente en 1915).',
              ],
            },
            {
              kind: 'list',
              title: 'Causes profondes',
              items: [
                'Rivalités économiques et coloniales entre puissances européennes.',
                'Nationalismes (Balkans, Alsace-Lorraine perdue par la France en 1871).',
                'Course aux armements.',
              ],
            },
            {
              kind: 'date',
              date: '28 juin 1914',
              event: "Assassinat de l'archiduc François-Ferdinand, héritier d'Autriche-Hongrie, à Sarajevo : c'est la cause immédiate.",
            },
            {
              kind: 'date',
              date: 'Août 1914',
              event: "Le jeu des alliances entraîne l'Europe dans la guerre (l'Allemagne déclare la guerre à la Russie le 1er août, à la France le 3 août).",
            },
          ],
        },
        {
          title: 'Les grandes phases de la guerre',
          blocks: [
            {
              kind: 'date',
              date: 'Septembre 1914',
              event: "Bataille de la Marne : l'avancée allemande vers Paris est stoppée. Fin de la guerre de mouvement.",
            },
            {
              kind: 'definition',
              term: 'Guerre de position (des tranchées)',
              definition: 'De 1915 à 1917, les armées s’enterrent dans des tranchées et s’usent dans des offensives meurtrières.',
            },
            {
              kind: 'date',
              date: '1916',
              event: 'Bataille de Verdun (février-décembre) : symbole de la violence de masse et de la résistance française.',
            },
            {
              kind: 'date',
              date: 'Avril 1917',
              event: 'Les États-Unis entrent en guerre aux côtés de l’Entente.',
            },
            {
              kind: 'date',
              date: '11 novembre 1918',
              event: "Signature de l'armistice : fin des combats. Bilan : environ 10 millions de soldats tués.",
            },
            {
              kind: 'definition',
              term: 'Guerre totale',
              definition:
                "Guerre qui mobilise toutes les ressources d'un pays : hommes, économie, colonies, propagande, et qui touche aussi les civils.",
            },
          ],
        },
        {
          title: "L'Afrique et les tirailleurs sénégalais",
          blocks: [
            {
              kind: 'definition',
              term: 'Tirailleurs sénégalais',
              definition:
                "Soldats africains de l'armée coloniale française. Le premier bataillon est créé en 1857 par Faidherbe. Ils venaient de toute l'Afrique occidentale, pas seulement du Sénégal.",
            },
            {
              kind: 'date',
              date: '1914',
              event: "Blaise Diagne est élu député du Sénégal : premier Africain noir élu à la Chambre des députés française.",
            },
            {
              kind: 'date',
              date: '1918',
              event: "Blaise Diagne, nommé commissaire de la République, mène une grande campagne de recrutement en AOF.",
            },
            {
              kind: 'list',
              title: "Rôle de l'AOF",
              items: [
                'Fourniture de soldats (Verdun, la Somme, le Chemin des Dames…).',
                'Fourniture de produits (arachide, céréales) et de main-d’œuvre.',
                'Recrutement souvent forcé, provoquant des fuites et des révoltes.',
              ],
            },
            {
              kind: 'warning',
              text: "Ne dis pas que les tirailleurs étaient tous volontaires : beaucoup ont été recrutés de force.",
            },
          ],
        },
        {
          title: 'La révolution russe',
          blocks: [
            {
              kind: 'date',
              date: 'Février 1917',
              event: 'Première révolution : le tsar Nicolas II abdique. Un gouvernement provisoire continue la guerre.',
            },
            {
              kind: 'date',
              date: 'Octobre 1917',
              event: 'Les bolcheviks dirigés par Lénine prennent le pouvoir à Petrograd.',
            },
            {
              kind: 'date',
              date: 'Mars 1918',
              event: 'Traité de Brest-Litovsk : la Russie bolchevique sort de la guerre.',
            },
            {
              kind: 'date',
              date: '1918-1921',
              event: "Guerre civile entre l'Armée rouge (créée par Trotski) et les Blancs ; victoire des bolcheviks.",
            },
            {
              kind: 'date',
              date: 'Décembre 1922',
              event: "Création de l'URSS (Union des républiques socialistes soviétiques).",
            },
            {
              kind: 'definition',
              term: 'Bolcheviks',
              definition: 'Révolutionnaires communistes russes dirigés par Lénine, qui veulent instaurer la dictature du prolétariat.',
            },
            {
              kind: 'warning',
              text: "Les révolutions de « février » et « octobre » 1917 suivent l'ancien calendrier russe : elles ont eu lieu en mars et novembre dans notre calendrier.",
            },
          ],
        },
        {
          title: 'Les traités de paix',
          blocks: [
            {
              kind: 'date',
              date: '28 juin 1919',
              event: "Traité de Versailles : l'Allemagne, jugée responsable, perd l'Alsace-Lorraine et ses colonies, doit payer des réparations et limiter son armée.",
            },
            {
              kind: 'definition',
              term: 'SDN (Société des Nations)',
              definition: 'Organisation créée en 1919 pour maintenir la paix, siège à Genève. Elle échoue faute de moyens.',
            },
            {
              kind: 'tip',
              text: "Composition : pour expliquer les conséquences de la guerre, classe-les en bilan humain, bilan matériel/économique et bilan politique (nouvelles frontières, révolution russe, SDN).",
            },
          ],
        },
      ],
      flashcards: [
        { front: "Quel événement déclenche la Première Guerre mondiale ?", back: "L'assassinat de l'archiduc François-Ferdinand à Sarajevo, le 28 juin 1914." },
        { front: 'Pays de la Triple Entente', back: 'France, Royaume-Uni, Russie.' },
        { front: 'Date de l’armistice de la Première Guerre mondiale', back: '11 novembre 1918.' },
        { front: 'Bataille symbole de 1916', back: 'Verdun (février-décembre 1916).' },
        { front: 'Qui a créé les tirailleurs sénégalais, et quand ?', back: 'Le gouverneur Faidherbe, en 1857.' },
        { front: 'Blaise Diagne', back: "Premier Africain noir élu député (1914) ; chargé du recrutement de soldats en AOF en 1918." },
        { front: 'Qui dirige la révolution d’octobre 1917 ?', back: 'Lénine, à la tête des bolcheviks.' },
        { front: 'Traité de Versailles', back: "28 juin 1919 : l'Allemagne perd l'Alsace-Lorraine et ses colonies, et doit payer des réparations." },
        { front: "Création de l'URSS", back: 'Décembre 1922.' },
      ],
      quiz: [
        {
          id: 'histoire-bfm-premiere-guerre-mondiale-q1',
          type: 'qcm',
          prompt: 'Quel événement est la cause immédiate de la Première Guerre mondiale ?',
          choices: [
            "L'assassinat de François-Ferdinand à Sarajevo",
            "L'invasion de la Pologne",
            'La révolution russe',
            'Le traité de Versailles',
          ],
          answer: 0,
          explanation: "Le 28 juin 1914, l'héritier d'Autriche-Hongrie est assassiné à Sarajevo ; le jeu des alliances fait ensuite basculer l'Europe dans la guerre.",
        },
        {
          id: 'histoire-bfm-premiere-guerre-mondiale-q2',
          type: 'vrai-faux',
          prompt: "L'armistice de la Première Guerre mondiale a été signé le 11 novembre 1918.",
          answer: true,
          explanation: "Le 11 novembre 1918, l'armistice met fin aux combats. La paix est signée ensuite, à Versailles, le 28 juin 1919.",
        },
        {
          id: 'histoire-bfm-premiere-guerre-mondiale-q3',
          type: 'trous',
          prompt: 'La Triple Entente regroupe la France, le ___ et la ___.',
          answers: ['Royaume-Uni', 'Russie'],
          bank: ['Royaume-Uni', 'Russie', 'Allemagne', 'Autriche-Hongrie', 'Italie'],
          explanation: "L'Entente (France, Royaume-Uni, Russie) s'oppose à la Triple Alliance (Allemagne, Autriche-Hongrie, Italie).",
        },
        {
          id: 'histoire-bfm-premiere-guerre-mondiale-q4',
          type: 'qcm',
          prompt: 'Qui organise en 1918 une grande campagne de recrutement de soldats en AOF ?',
          choices: ['Lamine Guèye', 'Léopold Sédar Senghor', 'Faidherbe', 'Blaise Diagne'],
          answer: 3,
          explanation: "Blaise Diagne, député du Sénégal depuis 1914, est nommé commissaire de la République en 1918 pour recruter des soldats en AOF.",
        },
        {
          id: 'histoire-bfm-premiere-guerre-mondiale-q5',
          type: 'vrai-faux',
          prompt: 'Les tirailleurs sénégalais venaient uniquement du Sénégal.',
          answer: false,
          explanation: "Le nom est trompeur : ils étaient recrutés dans toute l'Afrique occidentale française (Soudan, Guinée, Haute-Volta, etc.).",
        },
        {
          id: 'histoire-bfm-premiere-guerre-mondiale-q6',
          type: 'qcm',
          prompt: 'Quelle bataille de 1916 est devenue le symbole de la guerre des tranchées ?',
          choices: ['Verdun', 'La Marne', 'Stalingrad', 'Dien Bien Phu'],
          answer: 0,
          explanation: 'Verdun (février-décembre 1916) symbolise la violence de la guerre de position. Stalingrad et Dien Bien Phu appartiennent à d’autres conflits.',
        },
        {
          id: 'histoire-bfm-premiere-guerre-mondiale-q7',
          type: 'trous',
          prompt: 'En octobre 1917, les ___ dirigés par ___ prennent le pouvoir en Russie.',
          answers: ['bolcheviks', 'Lénine'],
          bank: ['bolcheviks', 'Lénine', 'Staline', 'Nicolas II', 'mencheviks'],
          explanation: "Lénine et les bolcheviks renversent le gouvernement provisoire. Staline ne prend le pouvoir qu'après la mort de Lénine (1924).",
        },
        {
          id: 'histoire-bfm-premiere-guerre-mondiale-q8',
          type: 'vrai-faux',
          prompt: 'Par le traité de Brest-Litovsk (1918), la Russie bolchevique sort de la guerre.',
          answer: true,
          explanation: "En mars 1918, Lénine signe la paix avec l'Allemagne pour consolider la révolution.",
        },
        {
          id: 'histoire-bfm-premiere-guerre-mondiale-q9',
          type: 'trous',
          prompt: "Le traité de ___, signé le 28 juin 1919, oblige l'Allemagne à rendre l'___ à la France.",
          answers: ['Versailles', 'Alsace-Lorraine'],
          bank: ['Versailles', 'Alsace-Lorraine', 'Brest-Litovsk', 'Rome', 'Pologne'],
          explanation: "Le traité de Versailles impose aussi à l'Allemagne la perte de ses colonies, des réparations et la limitation de son armée.",
        },
        {
          id: 'histoire-bfm-premiere-guerre-mondiale-q10',
          type: 'qcm',
          prompt: 'Que signifie « guerre totale » ?',
          choices: [
            'Une guerre qui dure plus de quatre ans',
            'Une guerre qui se déroule sur tous les continents',
            'Une guerre qui mobilise toutes les ressources d’un pays, y compris les civils',
          ],
          answer: 2,
          explanation: "La guerre totale mobilise les hommes, l'économie, les colonies et l'opinion (propagande) : toute la société est engagée.",
        },
      ],
    },

    // ───────────────────────── Chapitre 2 ─────────────────────────
    {
      id: 'histoire-bfm-crise-totalitarismes',
      title: 'La crise de 1929 et les totalitarismes',
      summary:
        'Le krach de Wall Street (1929) provoque une crise mondiale qui favorise la montée des régimes totalitaires en Europe.',
      essentials: [
        'Le « jeudi noir » (24 octobre 1929) à Wall Street déclenche une crise économique mondiale.',
        'Aux États-Unis, Roosevelt lance le New Deal à partir de 1933.',
        'La crise touche aussi les colonies : en AOF, le prix de l’arachide s’effondre.',
        "Trois régimes totalitaires : l'Italie fasciste (Mussolini), l'Allemagne nazie (Hitler), l'URSS stalinienne.",
        'Un totalitarisme repose sur un parti unique, un chef tout-puissant, la propagande et la terreur.',
      ],
      sections: [
        {
          title: 'La crise de 1929',
          blocks: [
            {
              kind: 'date',
              date: '24 octobre 1929',
              event: 'Jeudi noir : effondrement des cours de la Bourse de New York (Wall Street).',
            },
            {
              kind: 'definition',
              term: 'Krach boursier',
              definition: "Effondrement brutal des cours (prix) des actions en Bourse.",
            },
            {
              kind: 'list',
              title: 'Enchaînement de la crise',
              items: [
                'Krach boursier → faillites de banques → faillites d’entreprises.',
                'Chômage massif : environ un actif américain sur quatre en 1933.',
                'Les capitaux américains quittent l’Europe : la crise devient mondiale (années 1930).',
              ],
            },
            {
              kind: 'date',
              date: '1933',
              event: "Franklin D. Roosevelt, président des États-Unis, lance le New Deal : grands travaux, aides sociales, intervention de l'État dans l'économie.",
            },
            {
              kind: 'example',
              title: 'La crise en AOF',
              text: "Le prix de l'arachide s'effondre : les paysans du Sénégal s'appauvrissent, alors que les impôts coloniaux restent à payer.",
            },
          ],
        },
        {
          title: "L'Italie fasciste et l'Allemagne nazie",
          blocks: [
            {
              kind: 'date',
              date: 'Octobre 1922',
              event: 'Marche sur Rome : Mussolini (le « Duce ») devient chef du gouvernement italien.',
            },
            {
              kind: 'date',
              date: '30 janvier 1933',
              event: 'Adolf Hitler, chef du parti nazi (NSDAP), est nommé chancelier en Allemagne.',
            },
            {
              kind: 'date',
              date: '1935',
              event: "Lois de Nuremberg : les Juifs d'Allemagne perdent leurs droits de citoyens. La même année, l'Italie envahit l'Éthiopie.",
            },
            {
              kind: 'list',
              title: 'Idéologie nazie',
              items: [
                'Racisme et antisémitisme (prétendue supériorité d’une « race aryenne »).',
                'Recherche d’un « espace vital » à l’Est.',
                'Rejet de la démocratie et du traité de Versailles.',
              ],
            },
          ],
        },
        {
          title: "L'URSS de Staline",
          blocks: [
            {
              kind: 'text',
              text: 'Après la mort de Lénine (1924), Staline élimine ses rivaux et gouverne seul.',
            },
            {
              kind: 'list',
              items: [
                'Plans quinquennaux (à partir de 1928) : priorité à l’industrie lourde.',
                'Collectivisation forcée des terres : kolkhozes et sovkhozes.',
                'Terreur : goulag (camps de travail forcé), Grandes Purges (1936-1938).',
                'Culte de la personnalité autour de Staline.',
              ],
            },
          ],
        },
        {
          title: 'Les caractères du totalitarisme',
          blocks: [
            {
              kind: 'definition',
              term: 'Totalitarisme',
              definition:
                "Régime politique où l'État, dirigé par un parti unique et un chef, contrôle toute la société et la vie des individus.",
            },
            {
              kind: 'list',
              title: 'Points communs',
              items: [
                'Parti unique et chef charismatique (culte de la personnalité).',
                'Propagande et embrigadement de la jeunesse.',
                'Police politique et terreur contre les opposants.',
                "Contrôle de l'économie et de la culture.",
              ],
            },
            {
              kind: 'warning',
              text: "Ne confonds pas les idéologies : le nazisme et le fascisme sont d'extrême droite et nationalistes ; le stalinisme se réclame du communisme.",
            },
            {
              kind: 'tip',
              text: 'Pour comparer les régimes, construis un tableau : chef, parti, idéologie, méthodes, victimes.',
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Jeudi noir', back: '24 octobre 1929 : krach de la Bourse de New York (Wall Street).' },
        { front: 'New Deal', back: "Politique de Roosevelt (à partir de 1933) : grands travaux, aides sociales, intervention de l'État." },
        { front: 'Mussolini arrive au pouvoir…', back: 'En octobre 1922, après la marche sur Rome.' },
        { front: 'Hitler devient chancelier…', back: 'Le 30 janvier 1933.' },
        { front: 'Lois de Nuremberg', back: '1935 : les Juifs allemands perdent leurs droits de citoyens.' },
        { front: 'Goulag', back: "Système de camps de travail forcé en URSS." },
        { front: 'Kolkhoze', back: 'Exploitation agricole collective en URSS.' },
        { front: '4 caractères du totalitarisme', back: 'Parti unique et chef, propagande, terreur, contrôle de toute la société.' },
      ],
      quiz: [
        {
          id: 'histoire-bfm-crise-totalitarismes-q1',
          type: 'qcm',
          prompt: 'Où commence la crise de 1929 ?',
          choices: ['À Londres', 'À Berlin', 'À la Bourse de New York', 'À Paris'],
          answer: 2,
          explanation: 'Le krach de Wall Street (Bourse de New York), le 24 octobre 1929, déclenche la crise.',
        },
        {
          id: 'histoire-bfm-crise-totalitarismes-q2',
          type: 'vrai-faux',
          prompt: 'Le New Deal est la politique de Roosevelt pour sortir de la crise.',
          answer: true,
          explanation: "À partir de 1933, Roosevelt fait intervenir l'État : grands travaux, aides aux chômeurs, contrôle des banques.",
        },
        {
          id: 'histoire-bfm-crise-totalitarismes-q3',
          type: 'trous',
          prompt: 'Hitler est nommé ___ le 30 janvier ___.',
          answers: ['chancelier', '1933'],
          bank: ['chancelier', '1933', 'président', '1929', '1922'],
          explanation: "Hitler, chef du parti nazi, devient chancelier le 30 janvier 1933, puis supprime rapidement la démocratie.",
        },
        {
          id: 'histoire-bfm-crise-totalitarismes-q4',
          type: 'qcm',
          prompt: 'Quel dirigeant est surnommé le « Duce » ?',
          choices: ['Hitler', 'Mussolini', 'Staline', 'Franco'],
          answer: 1,
          explanation: 'Mussolini, fondateur du fascisme, prend le pouvoir en Italie en 1922 et se fait appeler le Duce (le guide).',
        },
        {
          id: 'histoire-bfm-crise-totalitarismes-q5',
          type: 'vrai-faux',
          prompt: 'La crise de 1929 n’a eu aucune conséquence en Afrique occidentale.',
          answer: false,
          explanation: "La crise touche les colonies : en AOF, le prix de l'arachide s'effondre et les paysans s'appauvrissent.",
        },
        {
          id: 'histoire-bfm-crise-totalitarismes-q6',
          type: 'trous',
          prompt: 'En URSS, Staline impose la collectivisation des terres dans des ___ et envoie les opposants au ___.',
          answers: ['kolkhozes', 'goulag'],
          bank: ['kolkhozes', 'goulag', 'ghettos', 'communes populaires', 'soviets'],
          explanation: "Les kolkhozes sont des fermes collectives ; le goulag est le système de camps de travail forcé.",
        },
        {
          id: 'histoire-bfm-crise-totalitarismes-q7',
          type: 'qcm',
          prompt: "Lequel de ces éléments n'est PAS un caractère d'un régime totalitaire ?",
          choices: ['Le parti unique', 'La propagande', 'Le pluralisme des partis', 'La police politique'],
          answer: 2,
          explanation: 'Le pluralisme (plusieurs partis) caractérise la démocratie. Les régimes totalitaires imposent un parti unique.',
        },
        {
          id: 'histoire-bfm-crise-totalitarismes-q8',
          type: 'vrai-faux',
          prompt: 'Les lois de Nuremberg (1935) privent les Juifs allemands de leurs droits de citoyens.',
          answer: true,
          explanation: "Ces lois antisémites marquent une étape de la persécution qui mènera au génocide pendant la Seconde Guerre mondiale.",
        },
        {
          id: 'histoire-bfm-crise-totalitarismes-q9',
          type: 'qcm',
          prompt: 'Quelle priorité fixent les plans quinquennaux de Staline ?',
          choices: ["L'industrie lourde", 'Les biens de consommation', 'Le tourisme', 'Les exportations agricoles vers les États-Unis'],
          answer: 0,
          explanation: "Les plans quinquennaux (à partir de 1928) privilégient l'industrie lourde (acier, charbon, machines).",
        },
      ],
    },

    // ───────────────────────── Chapitre 3 ─────────────────────────
    {
      id: 'histoire-bfm-seconde-guerre-mondiale',
      title: 'La Seconde Guerre mondiale (1939-1945)',
      summary:
        "Guerre d'anéantissement mondiale, elle oppose l'Axe aux Alliés, mobilise l'Afrique et s'achève en 1945 ; à Thiaroye, des tirailleurs sont massacrés en 1944.",
      essentials: [
        "La guerre commence avec l'invasion de la Pologne par l'Allemagne (1er septembre 1939).",
        "1942-1943 est le tournant (Stalingrad) ; l'Allemagne capitule le 8 mai 1945, le Japon le 2 septembre 1945.",
        'La Shoah : génocide d’environ 6 millions de Juifs par les nazis.',
        "L'AOF fournit soldats, produits et main-d'œuvre ; les tirailleurs participent à la libération de la France.",
        "Le 1er décembre 1944, à Thiaroye, l'armée française tire sur des tirailleurs qui réclamaient leur dû.",
      ],
      sections: [
        {
          title: 'Les grandes étapes',
          blocks: [
            { kind: 'date', date: '1er septembre 1939', event: "L'Allemagne envahit la Pologne. Le 3 septembre, la France et le Royaume-Uni lui déclarent la guerre." },
            { kind: 'date', date: 'Juin 1940', event: "Défaite de la France. Appel du général de Gaulle à Londres (18 juin) ; armistice demandé par le gouvernement du maréchal Pétain et signé le 22 juin." },
            { kind: 'date', date: '22 juin 1941', event: "L'Allemagne attaque l'URSS." },
            { kind: 'date', date: '7 décembre 1941', event: 'Le Japon attaque Pearl Harbor : les États-Unis entrent en guerre.' },
            { kind: 'date', date: '1942-1943', event: "Bataille de Stalingrad : première grande défaite allemande (février 1943). C'est le tournant de la guerre." },
            { kind: 'date', date: '6 juin 1944', event: 'Débarquement allié en Normandie.' },
            { kind: 'date', date: '8 mai 1945', event: "Capitulation de l'Allemagne : fin de la guerre en Europe." },
            { kind: 'date', date: '6 et 9 août 1945', event: 'Bombes atomiques américaines sur Hiroshima et Nagasaki.' },
            { kind: 'date', date: '2 septembre 1945', event: 'Capitulation du Japon : fin de la Seconde Guerre mondiale.' },
          ],
        },
        {
          title: "Une guerre d'anéantissement",
          blocks: [
            {
              kind: 'definition',
              term: 'Shoah (génocide des Juifs)',
              definition:
                "Extermination organisée d'environ 6 millions de Juifs d'Europe par les nazis, notamment dans des centres de mise à mort comme Auschwitz-Birkenau.",
            },
            {
              kind: 'definition',
              term: 'Génocide',
              definition: "Destruction volontaire et organisée d'un peuple, d'un groupe ethnique ou religieux.",
            },
            {
              kind: 'list',
              title: 'Les deux camps',
              items: [
                "L'Axe : Allemagne, Italie, Japon.",
                'Les Alliés : Royaume-Uni, France libre, URSS (à partir de 1941), États-Unis (à partir de 1941), Chine…',
              ],
            },
            { kind: 'text', text: 'Bilan : plus de 50 millions de morts, dont une majorité de civils.' },
          ],
        },
        {
          title: "L'AOF et le Sénégal dans la guerre",
          blocks: [
            { kind: 'date', date: 'Septembre 1940', event: "Bataille de Dakar : les Britanniques et les Français libres de De Gaulle tentent en vain de prendre Dakar, fidèle au régime de Vichy." },
            { kind: 'date', date: 'Novembre 1942', event: "Après le débarquement allié en Afrique du Nord, l'AOF se rallie aux Alliés." },
            {
              kind: 'list',
              title: "L'effort de guerre de l'AOF",
              items: [
                'Recrutement de tirailleurs, qui combattent en 1940 puis dans la libération de la France (débarquement de Provence, août 1944).',
                "Réquisitions de produits (arachide, riz, caoutchouc) et travail forcé.",
                'Pénuries et souffrances des populations.',
              ],
            },
            { kind: 'date', date: 'Janvier-février 1944', event: "Conférence de Brazzaville : la France libre promet des réformes dans ses colonies, mais refuse l'idée d'indépendance." },
          ],
        },
        {
          title: 'Le massacre de Thiaroye',
          blocks: [
            {
              kind: 'date',
              date: '1er décembre 1944',
              event:
                "Au camp militaire de Thiaroye, près de Dakar, l'armée française tire sur des tirailleurs rapatriés (anciens prisonniers de guerre) qui réclamaient le paiement de leurs soldes et primes.",
            },
            {
              kind: 'text',
              text: "Le bilan officiel français était de 35 morts ; les historiens estiment qu'il y a eu bien davantage de victimes. Le nombre exact reste inconnu.",
            },
            {
              kind: 'example',
              title: 'Mémoire',
              text: 'Le film « Camp de Thiaroye » (1988) d’Ousmane Sembène a fait connaître ce drame. Il est aujourd’hui un symbole de l’injustice coloniale.',
            },
            {
              kind: 'warning',
              text: "Thiaroye a lieu APRÈS la libération de Paris (août 1944) mais AVANT la fin de la guerre (mai 1945). Ne le place pas en 1945.",
            },
            {
              kind: 'tip',
              text: "Dans une composition sur l'Afrique dans la guerre, montre que la guerre a renforcé la conscience politique africaine et préparé la décolonisation.",
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Début de la Seconde Guerre mondiale', back: "1er septembre 1939 : l'Allemagne envahit la Pologne." },
        { front: 'Appel du 18 juin 1940', back: 'Appel du général de Gaulle, depuis Londres, à continuer le combat.' },
        { front: 'Tournant de la guerre (1942-1943)', back: 'La bataille de Stalingrad.' },
        { front: 'Débarquement de Normandie', back: '6 juin 1944.' },
        { front: "Capitulation de l'Allemagne", back: '8 mai 1945.' },
        { front: 'Fin de la Seconde Guerre mondiale', back: '2 septembre 1945 : capitulation du Japon.' },
        { front: 'Massacre de Thiaroye', back: "1er décembre 1944 : l'armée française tire sur des tirailleurs qui réclamaient leur solde." },
        { front: 'Conférence de Brazzaville', back: '1944 : la France libre promet des réformes coloniales, sans indépendance.' },
        { front: 'Shoah', back: 'Génocide d’environ 6 millions de Juifs par les nazis.' },
      ],
      quiz: [
        {
          id: 'histoire-bfm-seconde-guerre-mondiale-q1',
          type: 'qcm',
          prompt: 'Quel événement déclenche la Seconde Guerre mondiale ?',
          choices: ['Pearl Harbor', 'Le jeudi noir', 'La bataille de Dakar', "L'invasion de la Pologne par l'Allemagne"],
          answer: 3,
          explanation: "Le 1er septembre 1939, l'Allemagne envahit la Pologne ; la France et le Royaume-Uni déclarent la guerre le 3 septembre.",
        },
        {
          id: 'histoire-bfm-seconde-guerre-mondiale-q2',
          type: 'trous',
          prompt: 'Le massacre de Thiaroye a eu lieu le 1er décembre ___, près de ___.',
          answers: ['1944', 'Dakar'],
          bank: ['1944', 'Dakar', '1945', 'Saint-Louis', '1940'],
          explanation: "Le camp militaire de Thiaroye se trouve près de Dakar. Le drame a lieu le 1er décembre 1944.",
        },
        {
          id: 'histoire-bfm-seconde-guerre-mondiale-q3',
          type: 'vrai-faux',
          prompt: 'À Thiaroye, les tirailleurs réclamaient le paiement de leurs soldes et primes.',
          answer: true,
          explanation: "Anciens prisonniers de guerre rapatriés, ils réclamaient l'argent qui leur était dû ; l'armée française a ouvert le feu.",
        },
        {
          id: 'histoire-bfm-seconde-guerre-mondiale-q4',
          type: 'qcm',
          prompt: 'Quelle bataille marque le tournant de la guerre en 1942-1943 ?',
          choices: ['Stalingrad', 'Verdun', 'La Marne', 'Dakar'],
          answer: 0,
          explanation: "La capitulation allemande à Stalingrad (février 1943) est la première grande défaite de l'Allemagne.",
        },
        {
          id: 'histoire-bfm-seconde-guerre-mondiale-q5',
          type: 'vrai-faux',
          prompt: "En septembre 1940, les Français libres de De Gaulle réussissent à prendre Dakar.",
          answer: false,
          explanation: "La bataille de Dakar (septembre 1940) est un échec. L'AOF reste fidèle à Vichy jusqu'en novembre 1942.",
        },
        {
          id: 'histoire-bfm-seconde-guerre-mondiale-q6',
          type: 'qcm',
          prompt: 'Quels pays forment l’Axe ?',
          choices: ['Allemagne, Italie, Japon', 'Allemagne, URSS, Japon', 'Allemagne, Italie, Royaume-Uni', 'Italie, Japon, Chine'],
          answer: 0,
          explanation: "L'Axe regroupe l'Allemagne, l'Italie et le Japon. L'URSS, le Royaume-Uni et la Chine sont du côté des Alliés.",
        },
        {
          id: 'histoire-bfm-seconde-guerre-mondiale-q7',
          type: 'trous',
          prompt: "L'Allemagne capitule le ___ 1945 ; le Japon capitule le ___ 1945.",
          answers: ['8 mai', '2 septembre'],
          bank: ['8 mai', '2 septembre', '11 novembre', '6 juin', '6 août'],
          explanation: "Le 8 mai 1945 marque la fin de la guerre en Europe ; la capitulation japonaise du 2 septembre 1945 met fin à la guerre.",
        },
        {
          id: 'histoire-bfm-seconde-guerre-mondiale-q8',
          type: 'vrai-faux',
          prompt: 'La conférence de Brazzaville (1944) promet l’indépendance aux colonies africaines.',
          answer: false,
          explanation: "Elle promet des réformes (représentation, fin progressive du travail forcé) mais écarte toute idée d'indépendance.",
        },
        {
          id: 'histoire-bfm-seconde-guerre-mondiale-q9',
          type: 'qcm',
          prompt: "Qui a réalisé le film « Camp de Thiaroye » (1988) ?",
          choices: ['Djibril Diop Mambéty', 'Safi Faye', 'Ousmane Sembène', 'Souleymane Cissé'],
          answer: 2,
          explanation: "Ousmane Sembène, grand cinéaste et écrivain sénégalais, a coréalisé ce film qui a fait connaître le massacre.",
        },
        {
          id: 'histoire-bfm-seconde-guerre-mondiale-q10',
          type: 'qcm',
          prompt: 'Quand a lieu le débarquement allié en Normandie ?',
          choices: ['6 juin 1944', '15 août 1944', '8 mai 1945', '18 juin 1940'],
          answer: 0,
          explanation: 'Le débarquement de Normandie a lieu le 6 juin 1944. Le 15 août 1944, c’est le débarquement de Provence, où combattent de nombreux tirailleurs.',
        },
      ],
    },

    // ───────────────────────── Chapitre 4 ─────────────────────────
    {
      id: 'histoire-bfm-independance-senegal',
      title: "La décolonisation et l'indépendance du Sénégal",
      summary:
        "Après 1945, le Sénégal passe par des réformes successives (1946, 1956, 1958) puis par la Fédération du Mali avant de devenir indépendant en 1960.",
      essentials: [
        "1946 : l'Union française ; la loi Lamine Guèye donne la citoyenneté aux habitants des colonies.",
        '1956 : la loi-cadre Defferre instaure le suffrage universel et des gouvernements locaux.',
        '28 septembre 1958 : le Sénégal vote « oui » au référendum (Communauté française) ; la Guinée vote « non ».',
        '4 avril 1960 : accords de transfert des compétences (fête nationale) ; 20 juin 1960 : indépendance de la Fédération du Mali.',
        '20 août 1960 : éclatement de la Fédération, le Sénégal proclame son indépendance ; Senghor devient président.',
      ],
      sections: [
        {
          title: 'Les réformes après 1945',
          blocks: [
            {
              kind: 'definition',
              term: 'Décolonisation',
              definition: "Processus par lequel une colonie obtient son indépendance, de façon pacifique ou par la guerre.",
            },
            { kind: 'date', date: '1946', event: "Création de l'Union française. Loi Houphouët-Boigny : abolition du travail forcé. Loi Lamine Guèye : citoyenneté pour les habitants des colonies." },
            { kind: 'date', date: '1946', event: 'Création à Bamako du RDA (Rassemblement démocratique africain), présidé par Félix Houphouët-Boigny.' },
            { kind: 'date', date: '1948', event: 'Senghor fonde le BDS (Bloc démocratique sénégalais) après avoir quitté le parti de Lamine Guèye.' },
            {
              kind: 'date',
              date: '23 juin 1956',
              event: "Loi-cadre Defferre : suffrage universel, collège unique, création de conseils de gouvernement dans chaque territoire.",
            },
            {
              kind: 'warning',
              text: "Les « quatre communes » (Saint-Louis, Gorée, Rufisque, Dakar) avaient déjà des citoyens français (les « originaires ») avant 1946.",
            },
          ],
        },
        {
          title: 'Le référendum de 1958',
          blocks: [
            {
              kind: 'date',
              date: '28 septembre 1958',
              event: "Référendum sur la Communauté française proposée par De Gaulle. Le Sénégal vote « oui » : il devient un État autonome membre de la Communauté.",
            },
            {
              kind: 'date',
              date: '2 octobre 1958',
              event: 'La Guinée de Sékou Touré, qui a voté « non », devient indépendante.',
            },
            {
              kind: 'definition',
              term: 'Communauté française',
              definition: "Association créée en 1958 entre la France et ses anciennes colonies devenues États autonomes ; la France garde la défense, la monnaie et la politique étrangère.",
            },
          ],
        },
        {
          title: 'La Fédération du Mali et l’indépendance',
          blocks: [
            { kind: 'date', date: 'Janvier 1959', event: 'Création de la Fédération du Mali entre le Sénégal et le Soudan français (actuel Mali).' },
            { kind: 'date', date: '4 avril 1960', event: "Signature à Paris des accords de transfert des compétences. C'est la date de la fête nationale du Sénégal." },
            { kind: 'date', date: '20 juin 1960', event: 'Proclamation de l’indépendance de la Fédération du Mali.' },
            { kind: 'date', date: '20 août 1960', event: "Éclatement de la Fédération : le Sénégal proclame son indépendance." },
            { kind: 'date', date: '5 septembre 1960', event: 'Léopold Sédar Senghor est élu premier président de la République du Sénégal ; Mamadou Dia est président du Conseil.' },
            { kind: 'date', date: '22 septembre 1960', event: 'Le Soudan proclame son indépendance sous le nom de République du Mali (Modibo Keïta).' },
            {
              kind: 'list',
              title: "Causes de l'éclatement",
              items: [
                'Désaccords politiques entre dirigeants sénégalais et soudanais.',
                'Rivalité pour les postes (présidence de la Fédération).',
                'Conceptions différentes : fédération souple pour le Sénégal, État plus centralisé pour le Soudan.',
              ],
            },
            {
              kind: 'tip',
              text: "Retiens le trio de dates de 1960 : 4 avril (transfert des compétences), 20 juin (indépendance de la Fédération), 20 août (indépendance du Sénégal seul).",
            },
          ],
        },
        {
          title: 'Les acteurs sénégalais',
          blocks: [
            { kind: 'definition', term: 'Lamine Guèye', definition: 'Avocat et député, auteur de la loi de 1946 sur la citoyenneté ; figure de la vie politique sénégalaise.' },
            { kind: 'definition', term: 'Léopold Sédar Senghor', definition: 'Poète de la négritude, député, fondateur du BDS (1948), premier président du Sénégal (1960-1980).' },
            { kind: 'definition', term: 'Mamadou Dia', definition: 'Président du Conseil du Sénégal ; arrêté lors de la crise de décembre 1962 qui l’oppose à Senghor.' },
            { kind: 'definition', term: 'Valdiodio Ndiaye', definition: "Ministre de l'Intérieur ; il accueille De Gaulle à Dakar en août 1958 par un discours réclamant l'indépendance." },
          ],
        },
      ],
      flashcards: [
        { front: 'Loi Lamine Guèye', back: '1946 : citoyenneté accordée aux habitants des colonies.' },
        { front: 'Loi-cadre Defferre', back: '23 juin 1956 : suffrage universel, collège unique, conseils de gouvernement.' },
        { front: 'Référendum sur la Communauté', back: '28 septembre 1958 : le Sénégal vote « oui », la Guinée « non ».' },
        { front: 'Fête nationale du Sénégal', back: '4 avril (accords de transfert des compétences, 1960).' },
        { front: 'Membres de la Fédération du Mali', back: 'Le Sénégal et le Soudan français (actuel Mali).' },
        { front: 'Indépendance de la Fédération du Mali', back: '20 juin 1960.' },
        { front: 'Éclatement de la Fédération du Mali', back: '20 août 1960 : le Sénégal proclame son indépendance.' },
        { front: 'Premier président du Sénégal', back: 'Léopold Sédar Senghor, élu le 5 septembre 1960.' },
        { front: 'Les quatre communes', back: 'Saint-Louis, Gorée, Rufisque, Dakar.' },
      ],
      quiz: [
        {
          id: 'histoire-bfm-independance-senegal-q1',
          type: 'qcm',
          prompt: 'Que commémore la fête nationale sénégalaise du 4 avril ?',
          choices: [
            "L'éclatement de la Fédération du Mali",
            "L'élection de Senghor",
            'La signature des accords de transfert des compétences en 1960',
            'Le référendum de 1958',
          ],
          answer: 2,
          explanation: "Le 4 avril 1960, les accords de transfert des compétences sont signés à Paris : c'est la fête nationale.",
        },
        {
          id: 'histoire-bfm-independance-senegal-q2',
          type: 'trous',
          prompt: 'La Fédération du Mali (1959-1960) regroupait le Sénégal et le ___ (actuel Mali).',
          answers: ['Soudan français'],
          bank: ['Soudan français', 'Dahomey', 'Haute-Volta', 'Guinée'],
          explanation: "Le Dahomey et la Haute-Volta devaient y entrer mais se sont retirés : seuls le Sénégal et le Soudan forment la Fédération.",
        },
        {
          id: 'histoire-bfm-independance-senegal-q3',
          type: 'vrai-faux',
          prompt: 'Au référendum du 28 septembre 1958, le Sénégal a voté « non ».',
          answer: false,
          explanation: "Le Sénégal vote « oui » et entre dans la Communauté. C'est la Guinée qui vote « non » et devient indépendante le 2 octobre 1958.",
        },
        {
          id: 'histoire-bfm-independance-senegal-q4',
          type: 'qcm',
          prompt: 'Quelle loi de 1956 instaure le suffrage universel dans les territoires d’outre-mer ?',
          choices: ['La loi Lamine Guèye', 'La loi Houphouët-Boigny', 'La loi-cadre Defferre', 'La loi Diagne'],
          answer: 2,
          explanation: "La loi-cadre Defferre (23 juin 1956) crée le suffrage universel, le collège unique et des conseils de gouvernement.",
        },
        {
          id: 'histoire-bfm-independance-senegal-q5',
          type: 'trous',
          prompt: 'Le Sénégal proclame son indépendance le ___ 1960 et ___ est élu président le 5 septembre.',
          answers: ['20 août', 'Senghor'],
          bank: ['20 août', 'Senghor', '4 avril', 'Mamadou Dia', 'Lamine Guèye', '20 juin'],
          explanation: "Après l'éclatement de la Fédération (20 août 1960), Léopold Sédar Senghor devient le premier président du Sénégal.",
        },
        {
          id: 'histoire-bfm-independance-senegal-q6',
          type: 'vrai-faux',
          prompt: "La loi Houphouët-Boigny de 1946 abolit le travail forcé dans les colonies.",
          answer: true,
          explanation: "Votée en avril 1946, cette loi supprime le travail forcé, une revendication majeure des Africains.",
        },
        {
          id: 'histoire-bfm-independance-senegal-q7',
          type: 'qcm',
          prompt: 'Quel parti Senghor fonde-t-il en 1948 ?',
          choices: ['Le RDA', 'Le BDS', 'Le PAI', 'L’UPC'],
          answer: 1,
          explanation: "Senghor fonde le Bloc démocratique sénégalais (BDS) en 1948, en s'appuyant sur le monde rural.",
        },
        {
          id: 'histoire-bfm-independance-senegal-q8',
          type: 'vrai-faux',
          prompt: 'La Fédération du Mali est devenue indépendante le 20 juin 1960.',
          answer: true,
          explanation: "La Fédération est indépendante le 20 juin 1960, mais elle éclate dès le 20 août 1960.",
        },
        {
          id: 'histoire-bfm-independance-senegal-q9',
          type: 'qcm',
          prompt: 'Laquelle de ces villes ne faisait PAS partie des « quatre communes » ?',
          choices: ['Gorée', 'Rufisque', 'Saint-Louis', 'Kaolack'],
          answer: 3,
          explanation: 'Les quatre communes sont Saint-Louis, Gorée, Rufisque et Dakar. Kaolack n’en faisait pas partie.',
        },
      ],
    },

    // ───────────────────────── Chapitre 5 ─────────────────────────
    {
      id: 'histoire-bfm-onu-monde-apres-1945',
      title: "L'ONU et le monde après 1945",
      summary:
        "Créée en 1945 pour préserver la paix, l'ONU agit dans un monde divisé par la Guerre froide et transformé par la décolonisation.",
      essentials: [
        "L'ONU est créée par la Charte de San Francisco (26 juin 1945) ; son siège est à New York.",
        'Le Conseil de sécurité compte 5 membres permanents avec droit de veto : États-Unis, Russie (ex-URSS), Chine, Royaume-Uni, France.',
        'La Déclaration universelle des droits de l’homme est adoptée le 10 décembre 1948.',
        'La Guerre froide (1947-1991) oppose le bloc américain (OTAN) et le bloc soviétique (pacte de Varsovie).',
        "L'ONU soutient la décolonisation ; le Sénégal y est admis en septembre 1960.",
      ],
      sections: [
        {
          title: "La création de l'ONU",
          blocks: [
            { kind: 'date', date: '26 juin 1945', event: 'Signature de la Charte des Nations unies à San Francisco (50 États ; la Pologne signe peu après : 51 membres fondateurs).' },
            { kind: 'date', date: '24 octobre 1945', event: "Entrée en vigueur de la Charte : naissance officielle de l'ONU (Journée des Nations unies)." },
            {
              kind: 'list',
              title: "Buts de l'ONU",
              items: [
                'Maintenir la paix et la sécurité internationales.',
                'Développer la coopération entre les nations.',
                'Défendre les droits de l’homme et le droit des peuples à disposer d’eux-mêmes.',
              ],
            },
            {
              kind: 'warning',
              text: "L'ONU remplace la SDN (créée en 1919), qui avait échoué à empêcher la Seconde Guerre mondiale.",
            },
          ],
        },
        {
          title: "Les organes de l'ONU",
          blocks: [
            { kind: 'definition', term: 'Assemblée générale', definition: 'Réunit tous les États membres ; chaque État a une voix.' },
            {
              kind: 'definition',
              term: 'Conseil de sécurité',
              definition: '15 membres dont 5 permanents disposant du droit de veto : États-Unis, Russie (ex-URSS), Chine, Royaume-Uni, France.',
            },
            { kind: 'definition', term: 'Secrétaire général', definition: "Dirige l'administration de l'ONU (le Ghanéen Kofi Annan l'a été de 1997 à 2006)." },
            { kind: 'definition', term: 'Droit de veto', definition: "Droit pour un membre permanent de bloquer une décision du Conseil de sécurité." },
            {
              kind: 'list',
              title: 'Quelques institutions spécialisées',
              items: [
                'UNESCO (éducation, science, culture) : dirigée par le Sénégalais Amadou-Mahtar M’Bow de 1974 à 1987.',
                'OMS (santé), FAO (alimentation et agriculture), UNICEF (enfance).',
              ],
            },
            { kind: 'date', date: '10 décembre 1948', event: "Adoption à Paris de la Déclaration universelle des droits de l'homme." },
          ],
        },
        {
          title: 'La Guerre froide',
          blocks: [
            {
              kind: 'definition',
              term: 'Guerre froide',
              definition:
                "Affrontement (1947-1991) entre les États-Unis et l'URSS, sans guerre directe entre eux, mais avec des crises et des conflits localisés.",
            },
            { kind: 'date', date: '1949', event: "Création de l'OTAN, alliance militaire du bloc occidental." },
            { kind: 'date', date: '1955', event: 'Création du pacte de Varsovie, alliance militaire du bloc soviétique.' },
            { kind: 'date', date: '1961', event: 'Construction du mur de Berlin, symbole de la division de l’Europe.' },
            { kind: 'date', date: '9 novembre 1989', event: 'Chute du mur de Berlin.' },
            { kind: 'date', date: 'Décembre 1991', event: "Disparition de l'URSS : fin de la Guerre froide." },
          ],
        },
        {
          title: "L'ONU, la décolonisation et l'Afrique",
          blocks: [
            {
              kind: 'text',
              text: "La Charte affirme le droit des peuples à disposer d'eux-mêmes. L'ONU devient une tribune pour les peuples colonisés.",
            },
            { kind: 'date', date: '1960', event: "« Année de l'Afrique » : de nombreux pays africains deviennent indépendants et entrent à l'ONU. Le Sénégal est admis le 28 septembre 1960." },
            { kind: 'example', title: 'Casques bleus', text: "Soldats de l'ONU chargés de maintenir la paix. Le Sénégal participe régulièrement à ces missions." },
            {
              kind: 'tip',
              text: "Pour juger l'action de l'ONU, présente ses réussites (aide humanitaire, décolonisation, santé) puis ses limites (veto, manque de moyens, guerres non empêchées).",
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Charte des Nations unies', back: 'Signée à San Francisco le 26 juin 1945.' },
        { front: "Siège de l'ONU", back: 'New York.' },
        { front: 'Les 5 membres permanents du Conseil de sécurité', back: 'États-Unis, Russie (ex-URSS), Chine, Royaume-Uni, France.' },
        { front: 'Droit de veto', back: 'Pouvoir d’un membre permanent de bloquer une décision du Conseil de sécurité.' },
        { front: "Déclaration universelle des droits de l'homme", back: '10 décembre 1948, à Paris.' },
        { front: 'Amadou-Mahtar M’Bow', back: "Sénégalais, directeur général de l'UNESCO de 1974 à 1987." },
        { front: 'OTAN / pacte de Varsovie', back: 'Alliance occidentale (1949) / alliance du bloc soviétique (1955).' },
        { front: 'Admission du Sénégal à l’ONU', back: '28 septembre 1960.' },
      ],
      quiz: [
        {
          id: 'histoire-bfm-onu-monde-apres-1945-q1',
          type: 'qcm',
          prompt: "Dans quelle ville la Charte des Nations unies a-t-elle été signée ?",
          choices: ['New York', 'Genève', 'San Francisco', 'Paris'],
          answer: 2,
          explanation: "La Charte est signée à San Francisco le 26 juin 1945. Le siège de l'ONU est ensuite installé à New York.",
        },
        {
          id: 'histoire-bfm-onu-monde-apres-1945-q2',
          type: 'vrai-faux',
          prompt: "L'Allemagne est membre permanent du Conseil de sécurité.",
          answer: false,
          explanation: 'Les 5 membres permanents sont les États-Unis, la Russie, la Chine, le Royaume-Uni et la France (les vainqueurs de 1945).',
        },
        {
          id: 'histoire-bfm-onu-monde-apres-1945-q3',
          type: 'trous',
          prompt: "La Déclaration universelle des droits de l'homme est adoptée le 10 décembre ___ à ___.",
          answers: ['1948', 'Paris'],
          bank: ['1948', 'Paris', '1945', 'New York', 'Genève'],
          explanation: "Elle est adoptée par l'Assemblée générale de l'ONU réunie à Paris, le 10 décembre 1948.",
        },
        {
          id: 'histoire-bfm-onu-monde-apres-1945-q4',
          type: 'qcm',
          prompt: 'Quelle organisation l’ONU remplace-t-elle ?',
          choices: ['La SDN', "L'OTAN", "L'OUA", 'La CEDEAO'],
          answer: 0,
          explanation: 'La Société des Nations (SDN), créée en 1919, avait échoué à empêcher la Seconde Guerre mondiale.',
        },
        {
          id: 'histoire-bfm-onu-monde-apres-1945-q5',
          type: 'vrai-faux',
          prompt: "Le Sénégalais Amadou-Mahtar M'Bow a dirigé l'UNESCO.",
          answer: true,
          explanation: "Il a été directeur général de l'UNESCO de 1974 à 1987, premier Africain à diriger une grande agence de l'ONU.",
        },
        {
          id: 'histoire-bfm-onu-monde-apres-1945-q6',
          type: 'trous',
          prompt: "Pendant la Guerre froide, l'alliance du bloc occidental est l'___ ; celle du bloc soviétique est le pacte de ___.",
          answers: ['OTAN', 'Varsovie'],
          bank: ['OTAN', 'Varsovie', 'ONU', 'Moscou', 'Berlin'],
          explanation: "L'OTAN est créée en 1949 et le pacte de Varsovie en 1955.",
        },
        {
          id: 'histoire-bfm-onu-monde-apres-1945-q7',
          type: 'qcm',
          prompt: 'Combien de membres compte le Conseil de sécurité ?',
          choices: ['5', '10', '15', '51'],
          answer: 2,
          explanation: '15 membres : 5 permanents avec droit de veto et 10 non permanents élus pour deux ans. 51 est le nombre d’États fondateurs de l’ONU.',
        },
        {
          id: 'histoire-bfm-onu-monde-apres-1945-q8',
          type: 'vrai-faux',
          prompt: 'Le mur de Berlin est tombé le 9 novembre 1989.',
          answer: true,
          explanation: "Construit en 1961, le mur tombe le 9 novembre 1989 ; l'URSS disparaît en décembre 1991.",
        },
        {
          id: 'histoire-bfm-onu-monde-apres-1945-q9',
          type: 'qcm',
          prompt: "Pourquoi appelle-t-on 1960 « l'année de l'Afrique » ?",
          choices: [
            "Parce que l'OUA est créée cette année-là",
            "Parce que l'ONU s'installe en Afrique",
            'Parce que de nombreux pays africains deviennent indépendants',
          ],
          answer: 2,
          explanation: "En 1960, de nombreux États africains (dont le Sénégal) accèdent à l'indépendance et entrent à l'ONU. L'OUA naît en 1963.",
        },
      ],
    },
  ],
};

export default subject;
