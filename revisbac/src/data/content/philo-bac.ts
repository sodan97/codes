import type { Subject } from '../types';

const subject: Subject = {
  id: 'philo-bac',
  name: 'Philosophie',
  icon: '🤔',
  color: '#DB2777',
  tracks: ['bac-s', 'bac-l'],
  chapters: [
    // ───────────────────────── Chapitre 1 ─────────────────────────
    {
      id: 'philo-bac-methodologie',
      title: 'Méthodologie : dissertation et commentaire de texte',
      summary:
        "Réussir l'épreuve de philosophie, c'est maîtriser deux exercices : la dissertation, qui discute une question, et le commentaire, qui explique et évalue un texte.",
      essentials: [
        'Tout commence par l’analyse du sujet : définir les termes et dégager un problème (la problématique).',
        "Introduction de dissertation : amorce, position du problème (problématique), annonce du plan.",
        'Le développement progresse en 2 ou 3 parties argumentées, illustrées d’exemples et de références.',
        "Commentaire : dégager thème, thèse, problème et structure du texte, puis faire l'étude ordonnée et l'intérêt philosophique.",
        'Pièges majeurs : hors-sujet, récitation de cours, paraphrase, catalogue d’auteurs.',
      ],
      sections: [
        {
          title: 'Analyser le sujet de dissertation',
          blocks: [
            {
              kind: 'list',
              items: [
                'Lire le sujet plusieurs fois ; souligner les mots clés.',
                'Définir chaque terme (sens courant, sens philosophique) et repérer les liens entre eux.',
                'Repérer les présupposés : ce que le sujet suppose sans le dire.',
                'Formuler le problème : une tension entre deux réponses possibles.',
              ],
            },
            {
              kind: 'definition',
              term: 'Problématique',
              definition:
                "Question centrale (ou ensemble de questions) qui fait apparaître la difficulté du sujet et guide toute la réflexion.",
            },
            {
              kind: 'example',
              title: '« Peut-on être libre sans loi ? »',
              text: "Paradoxe : la loi semble limiter la liberté, mais sans loi, la liberté de chacun est menacée par la force des autres. D'où le problème : la loi est-elle un obstacle ou une condition de la liberté ?",
            },
          ],
        },
        {
          title: 'Construire la dissertation',
          blocks: [
            {
              kind: 'list',
              title: 'Introduction',
              items: [
                'Amorce : une situation, un exemple ou une idée qui mène au sujet (éviter « Depuis la nuit des temps… »).',
                'Position du problème : définitions, paradoxe, problématique.',
                'Annonce du plan (claire et brève).',
              ],
            },
            {
              kind: 'list',
              title: 'Développement',
              items: [
                'Plan dialectique fréquent : thèse, antithèse, dépassement (synthèse).',
                'Chaque partie : une idée directrice, des arguments, des exemples, une référence d’auteur expliquée.',
                'Des transitions qui montrent pourquoi on passe à la partie suivante.',
              ],
            },
            {
              kind: 'list',
              title: 'Conclusion',
              items: ['Bilan du parcours.', 'Réponse nette au problème posé.', 'Ouverture éventuelle (pertinente, pas obligatoire).'],
            },
            {
              kind: 'warning',
              text: 'Une citation ne remplace pas un argument : il faut l’expliquer et montrer en quoi elle répond au problème.',
            },
          ],
        },
        {
          title: 'Le commentaire de texte',
          blocks: [
            { kind: 'definition', term: 'Thème', definition: 'Ce dont parle le texte (ex. : le travail).' },
            { kind: 'definition', term: 'Thèse', definition: 'Ce que l’auteur affirme sur ce thème (ex. : le travail libère l’homme).' },
            { kind: 'definition', term: 'Problème', definition: 'La question à laquelle le texte répond.' },
            { kind: 'definition', term: 'Enjeu', definition: 'Ce qui est en jeu, l’importance de la question (morale, politique, existentielle…).' },
            {
              kind: 'list',
              title: 'Les étapes',
              items: [
                'Introduction : thème, problème, thèse, structure (mouvements du texte).',
                'Étude ordonnée : expliquer le texte pas à pas, en suivant ses mouvements.',
                "Intérêt philosophique : évaluer la thèse (portée et limites) en la confrontant à d'autres auteurs.",
                'Conclusion : bilan et réponse.',
              ],
            },
            {
              kind: 'warning',
              text: "La paraphrase (répéter le texte avec d'autres mots) n'est pas une explication : il faut dire POURQUOI l'auteur affirme ce qu'il affirme.",
            },
          ],
        },
        {
          title: 'Conseils de gestion du temps',
          blocks: [
            {
              kind: 'tip',
              text: "Consacre un temps important au brouillon (analyse, problématique, plan détaillé) ; rédige directement au propre l'introduction et la conclusion soignées.",
            },
            {
              kind: 'tip',
              text: "Prépare des exemples variés : vie quotidienne, histoire, actualité, littérature africaine (Cheikh Hamidou Kane, Senghor…).",
            },
            {
              kind: 'warning',
              text: "Ne fais pas de « catalogue » d'auteurs : mieux vaut une référence bien expliquée que cinq noms cités.",
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Les 3 étapes de l’introduction de dissertation', back: 'Amorce, position du problème (problématique), annonce du plan.' },
        { front: 'Plan dialectique', back: 'Thèse, antithèse, dépassement (synthèse).' },
        { front: 'Problématique', back: 'Question centrale qui fait apparaître la difficulté du sujet.' },
        { front: 'Thème d’un texte', back: 'Ce dont parle le texte.' },
        { front: 'Thèse d’un texte', back: "Ce que l'auteur affirme à propos du thème." },
        { front: 'Les 2 parties du commentaire (développement)', back: 'Étude ordonnée, puis intérêt philosophique.' },
        { front: 'Paraphrase', back: 'Répéter le texte avec d’autres mots sans l’expliquer : à éviter.' },
        { front: 'Contenu de la conclusion', back: 'Bilan, réponse au problème, ouverture éventuelle.' },
      ],
      quiz: [
        {
          id: 'philo-bac-methodologie-q1',
          type: 'qcm',
          prompt: "Quel élément ne fait PAS partie de l'introduction d'une dissertation ?",
          choices: ['L’amorce', 'La problématique', 'L’annonce du plan', 'La réponse finale détaillée'],
          answer: 3,
          explanation: 'La réponse finale se trouve dans la conclusion. L’introduction pose le problème sans le résoudre.',
        },
        {
          id: 'philo-bac-methodologie-q2',
          type: 'vrai-faux',
          prompt: 'Dans un commentaire, la paraphrase est une bonne méthode d’explication.',
          answer: false,
          explanation: "La paraphrase répète sans expliquer. Il faut dégager les arguments, les concepts et leur enchaînement.",
        },
        {
          id: 'philo-bac-methodologie-q3',
          type: 'trous',
          prompt: 'Le ___ est ce dont parle le texte ; la ___ est ce que l’auteur affirme.',
          answers: ['thème', 'thèse'],
          bank: ['thème', 'thèse', 'problème', 'enjeu', 'amorce'],
          explanation: 'Thème = le sujet du texte ; thèse = la position de l’auteur ; problème = la question à laquelle il répond.',
        },
        {
          id: 'philo-bac-methodologie-q4',
          type: 'qcm',
          prompt: 'Quelle est la structure classique d’un plan dialectique ?',
          choices: [
            'Thèse, antithèse, dépassement',
            'Introduction, citation, conclusion',
            'Causes, conséquences, solutions',
            'Définition, exemple, résumé',
          ],
          answer: 0,
          explanation: 'Le plan dialectique examine une thèse, la discute, puis dépasse l’opposition.',
        },
        {
          id: 'philo-bac-methodologie-q5',
          type: 'vrai-faux',
          prompt: "Une citation d'auteur doit toujours être expliquée et reliée au problème.",
          answer: true,
          explanation: "Une citation non expliquée n'apporte rien : elle doit servir l'argumentation.",
        },
        {
          id: 'philo-bac-methodologie-q6',
          type: 'trous',
          prompt: "Le développement d'un commentaire comprend l'étude ___ puis l'intérêt ___.",
          answers: ['ordonnée', 'philosophique'],
          bank: ['ordonnée', 'philosophique', 'historique', 'littéraire', 'dialectique'],
          explanation: "L'étude ordonnée explique le texte pas à pas ; l'intérêt philosophique en évalue la portée et les limites.",
        },
        {
          id: 'philo-bac-methodologie-q7',
          type: 'qcm',
          prompt: 'Que doit contenir une bonne conclusion ?',
          choices: [
            'Une nouvelle partie avec de nouveaux arguments',
            'Uniquement une citation',
            'Un bilan et une réponse claire au problème',
            'La recopie de l’introduction',
          ],
          answer: 2,
          explanation: "La conclusion fait le bilan, répond au problème et peut proposer une ouverture pertinente.",
        },
        {
          id: 'philo-bac-methodologie-q8',
          type: 'vrai-faux',
          prompt: "Repérer les présupposés d'un sujet aide à construire la problématique.",
          answer: true,
          explanation: 'Les présupposés révèlent ce que le sujet tient pour acquis : les interroger fait naître le problème.',
        },
        {
          id: 'philo-bac-methodologie-q9',
          type: 'qcm',
          prompt: 'Qu’appelle-t-on « catalogue d’auteurs » ?',
          choices: [
            'Une bibliographie en fin de copie',
            'Une liste de noms et de citations sans véritable explication',
            'Un plan en trois parties',
          ],
          answer: 1,
          explanation: "C'est une erreur : les correcteurs attendent une réflexion personnelle appuyée sur quelques références bien comprises.",
        },
      ],
    },

    // ───────────────────────── Chapitre 2 ─────────────────────────
    {
      id: 'philo-bac-conscience-inconscient',
      title: "La conscience et l'inconscient",
      summary:
        "La conscience fait de l'homme un sujet qui se connaît et répond de ses actes ; l'hypothèse freudienne de l'inconscient montre qu'il n'est pas entièrement transparent à lui-même.",
      essentials: [
        "La conscience est la connaissance qu'a le sujet de ses états, de ses actes et du monde.",
        'Descartes : « Je pense, donc je suis » ; la conscience est la première certitude.',
        "Freud : le psychisme comporte un inconscient (refoulement, rêves, lapsus) ; « le moi n'est pas maître dans sa propre maison ».",
        "Sartre critique l'inconscient : il parle plutôt de mauvaise foi, pour préserver la responsabilité.",
        "Fanon analyse comment la domination coloniale aliène la conscience du colonisé.",
      ],
      sections: [
        {
          title: 'Définitions',
          blocks: [
            { kind: 'definition', term: 'Conscience psychologique', definition: "Connaissance immédiate qu'a le sujet de ses états et de ses actes, et du monde qui l'entoure." },
            { kind: 'definition', term: 'Conscience morale', definition: 'Capacité de juger le bien et le mal de ses propres actions.' },
            { kind: 'definition', term: 'Conscience réfléchie', definition: 'Retour de la pensée sur elle-même : je sais que je pense.' },
            { kind: 'definition', term: 'Inconscient (sens freudien)', definition: "Partie du psychisme formée de désirs refoulés, inaccessibles directement à la conscience, mais qui agissent sur nos comportements." },
          ],
        },
        {
          title: 'La conscience, fondement du sujet',
          blocks: [
            {
              kind: 'definition',
              term: 'Descartes (1596-1650)',
              definition:
                "Dans le Discours de la méthode (1637), après avoir tout mis en doute, il trouve une certitude : « Je pense, donc je suis ». Le sujet est une « chose qui pense ».",
            },
            {
              kind: 'definition',
              term: 'Pascal (1623-1662)',
              definition: "« L'homme n'est qu'un roseau, le plus faible de la nature ; mais c'est un roseau pensant » (Pensées) : la pensée fait sa dignité.",
            },
            {
              kind: 'definition',
              term: 'Husserl (1859-1938)',
              definition: "« Toute conscience est conscience de quelque chose » : la conscience est intentionnalité, toujours tournée vers un objet.",
            },
            {
              kind: 'definition',
              term: 'Leibniz (1646-1716)',
              definition: "Les « petites perceptions » : nous percevons sans nous en apercevoir (le bruit de chaque vague dans le bruit de la mer).",
            },
          ],
        },
        {
          title: "L'hypothèse de l'inconscient : Freud",
          blocks: [
            {
              kind: 'list',
              items: [
                "Freud (1856-1939), fondateur de la psychanalyse.",
                'Le refoulement : des désirs inacceptables sont repoussés hors de la conscience.',
                "Ils reviennent sous forme déguisée : rêves, lapsus, actes manqués, symptômes névrotiques.",
                "L'interprétation des rêves est « la voie royale » qui mène à la connaissance de l'inconscient (L'Interprétation des rêves, 1900).",
              ],
            },
            {
              kind: 'list',
              title: 'La seconde topique',
              items: [
                'Le ça : réservoir des pulsions.',
                'Le surmoi : intériorisation des interdits (parents, société).',
                'Le moi : instance qui tente de concilier le ça, le surmoi et la réalité.',
              ],
            },
            {
              kind: 'text',
              text: "Freud parle de trois « blessures narcissiques » infligées à l'homme : Copernic (la Terre n'est pas le centre), Darwin (l'homme descend de l'animal), la psychanalyse (« le moi n'est pas maître dans sa propre maison »).",
            },
          ],
        },
        {
          title: "Critiques et prolongements",
          blocks: [
            {
              kind: 'definition',
              term: 'Sartre et la mauvaise foi',
              definition:
                "Pour Sartre (L'Être et le Néant, 1943), l'inconscient risque d'excuser l'homme. Il parle de « mauvaise foi » : se mentir à soi-même pour fuir sa liberté et sa responsabilité.",
            },
            {
              kind: 'definition',
              term: 'Frantz Fanon',
              definition:
                "Psychiatre martiniquais engagé en Algérie. Dans Peau noire, masques blancs (1952), il analyse l'aliénation psychique produite par le racisme colonial.",
            },
            {
              kind: 'warning',
              text: "« Je pense, donc je suis » est de Descartes, pas de Pascal. Et « le moi n'est pas maître dans sa propre maison » est de Freud, pas de Sartre.",
            },
            {
              kind: 'tip',
              text: "Sujet : « Suis-je ce que j'ai conscience d'être ? » Plan possible : I. La conscience me fait connaître moi-même (Descartes) ; II. Mais l'inconscient limite cette connaissance (Freud) ; III. Je reste responsable de ce que je fais de moi (Sartre).",
            },
          ],
        },
      ],
      flashcards: [
        { front: '« Je pense, donc je suis »', back: 'Descartes, Discours de la méthode (1637).' },
        { front: '« Le moi n’est pas maître dans sa propre maison »', back: 'Freud.' },
        { front: '« Toute conscience est conscience de quelque chose »', back: 'Husserl : la conscience est intentionnalité.' },
        { front: '« L’homme n’est qu’un roseau… mais c’est un roseau pensant »', back: 'Pascal, Pensées.' },
        { front: 'Ça, moi, surmoi', back: 'Seconde topique de Freud.' },
        { front: 'Refoulement', back: 'Mécanisme qui repousse hors de la conscience des désirs inacceptables.' },
        { front: 'Mauvaise foi', back: 'Sartre : se mentir à soi-même pour fuir sa liberté.' },
        { front: 'Petites perceptions', back: 'Leibniz : perceptions dont nous n’avons pas conscience.' },
        { front: 'Peau noire, masques blancs', back: "Frantz Fanon (1952) : l'aliénation du colonisé." },
      ],
      quiz: [
        {
          id: 'philo-bac-conscience-inconscient-q1',
          type: 'qcm',
          prompt: 'Qui a écrit « Je pense, donc je suis » ?',
          choices: ['Pascal', 'Sartre', 'Kant', 'Descartes'],
          answer: 3,
          explanation: 'Descartes, dans le Discours de la méthode (1637), IVe partie.',
        },
        {
          id: 'philo-bac-conscience-inconscient-q2',
          type: 'vrai-faux',
          prompt: "Pour Freud, l'interprétation des rêves est « la voie royale » vers la connaissance de l'inconscient.",
          answer: true,
          explanation: "Dans L'Interprétation des rêves, Freud voit dans le rêve la réalisation déguisée d'un désir refoulé.",
        },
        {
          id: 'philo-bac-conscience-inconscient-q3',
          type: 'trous',
          prompt: 'Dans la seconde topique de Freud, le ___ est le réservoir des pulsions et le ___ intériorise les interdits.',
          answers: ['ça', 'surmoi'],
          bank: ['ça', 'surmoi', 'moi', 'cogito', 'préconscient'],
          explanation: 'Le moi essaie de concilier les exigences du ça, du surmoi et de la réalité.',
        },
        {
          id: 'philo-bac-conscience-inconscient-q4',
          type: 'qcm',
          prompt: "Quel philosophe oppose à l'inconscient la notion de « mauvaise foi » ?",
          choices: ['Leibniz', 'Freud', 'Sartre', 'Spinoza'],
          answer: 2,
          explanation: "Sartre (L'Être et le Néant) refuse que l'inconscient serve d'excuse : l'homme se ment à lui-même pour fuir sa responsabilité.",
        },
        {
          id: 'philo-bac-conscience-inconscient-q5',
          type: 'vrai-faux',
          prompt: "La phrase « l'homme est un roseau pensant » est de Descartes.",
          answer: false,
          explanation: "Elle est de Pascal (Pensées). Descartes est l'auteur du cogito.",
        },
        {
          id: 'philo-bac-conscience-inconscient-q6',
          type: 'qcm',
          prompt: "Pour Husserl, la conscience est…",
          choices: ['Toujours conscience de quelque chose', 'Une illusion', 'Une chose matérielle', 'Le surmoi'],
          answer: 0,
          explanation: "C'est l'intentionnalité : la conscience est toujours visée d'un objet.",
        },
        {
          id: 'philo-bac-conscience-inconscient-q7',
          type: 'trous',
          prompt: "Les lapsus et les actes ___ révèlent, selon Freud, des désirs ___.",
          answers: ['manqués', 'refoulés'],
          bank: ['manqués', 'refoulés', 'réussis', 'conscients', 'raisonnables'],
          explanation: "Ce sont des « ratés » qui laissent passer un désir que la conscience avait repoussé.",
        },
        {
          id: 'philo-bac-conscience-inconscient-q8',
          type: 'vrai-faux',
          prompt: 'Frantz Fanon a analysé les effets psychologiques de la domination coloniale.',
          answer: true,
          explanation: "Psychiatre, il montre dans Peau noire, masques blancs (1952) comment le racisme colonial aliène le colonisé.",
        },
        {
          id: 'philo-bac-conscience-inconscient-q9',
          type: 'qcm',
          prompt: 'À quel philosophe associe-t-on les « petites perceptions » ?',
          choices: ['Leibniz', 'Hegel', 'Husserl', 'Freud'],
          answer: 0,
          explanation: "Leibniz (Nouveaux essais sur l'entendement humain) : nous percevons sans nous en apercevoir, comme chaque vague dans le bruit de la mer.",
        },
      ],
    },

    // ───────────────────────── Chapitre 3 ─────────────────────────
    {
      id: 'philo-bac-liberte',
      title: 'La liberté',
      summary:
        "Sommes-nous libres ou déterminés ? La liberté peut être pensée comme libre arbitre, comme illusion, comme connaissance de la nécessité ou comme obéissance à la loi qu'on se donne.",
      essentials: [
        'Liberté au sens courant : faire ce que l’on veut ; au sens philosophique : se déterminer soi-même.',
        "Descartes défend le libre arbitre ; la liberté d'indifférence en est « le plus bas degré ».",
        "Spinoza : le libre arbitre est une illusion née de l'ignorance des causes ; être libre, c'est comprendre la nécessité.",
        "Sartre : « l'existence précède l'essence », l'homme est « condamné à être libre » et responsable.",
        "Rousseau et Kant : la vraie liberté est l'autonomie, obéir à la loi qu'on s'est prescrite.",
      ],
      sections: [
        {
          title: 'Définitions',
          blocks: [
            { kind: 'definition', term: 'Libre arbitre', definition: 'Pouvoir de choisir entre plusieurs possibles sans y être contraint.' },
            { kind: 'definition', term: 'Déterminisme', definition: 'Idée que tout événement a des causes qui le produisent nécessairement.' },
            { kind: 'definition', term: 'Fatalisme', definition: "Croyance que tout est écrit d'avance, quoi que l'on fasse." },
            { kind: 'definition', term: 'Autonomie', definition: 'Fait de se donner à soi-même sa propre loi (contraire : hétéronomie).' },
            {
              kind: 'warning',
              text: "Ne confonds pas déterminisme et fatalisme : le déterminisme permet d'agir en connaissant les causes ; le fatalisme décourage l'action.",
            },
          ],
        },
        {
          title: 'Le libre arbitre et sa critique',
          blocks: [
            {
              kind: 'definition',
              term: 'Descartes',
              definition:
                "La volonté est infinie : c'est par elle que l'homme ressemble à Dieu. Choisir sans raison (liberté d'indifférence) est « le plus bas degré de la liberté » (Méditations métaphysiques, IV).",
            },
            {
              kind: 'definition',
              term: 'Spinoza (1632-1677)',
              definition:
                "Les hommes se croient libres parce qu'ils sont conscients de leurs actions mais ignorent les causes qui les déterminent (Éthique). Exemple de la pierre qui, consciente, se croirait libre de rouler.",
            },
            {
              kind: 'text',
              text: 'Pour Spinoza, la vraie liberté consiste à agir selon la nécessité de sa propre nature, éclairée par la raison.',
            },
            {
              kind: 'list',
              title: 'Autres déterminismes',
              items: ['Psychique (Freud : l’inconscient).', 'Social (Durkheim, Bourdieu : la société façonne nos choix).', 'Biologique (gènes).'],
            },
          ],
        },
        {
          title: "L'existentialisme : Sartre",
          blocks: [
            {
              kind: 'definition',
              term: "« L'existence précède l'essence »",
              definition:
                "Sartre (L'existentialisme est un humanisme, 1946) : l'homme existe d'abord, puis se définit par ses actes. Il n'y a pas de nature humaine fixée à l'avance.",
            },
            {
              kind: 'definition',
              term: '« L’homme est condamné à être libre »',
              definition: "Sartre : il ne peut pas ne pas choisir ; il est donc responsable de ce qu'il fait.",
            },
            {
              kind: 'example',
              title: 'Mauvaise foi',
              text: "Dire « je n'avais pas le choix » pour fuir sa responsabilité est, pour Sartre, un acte de mauvaise foi.",
            },
          ],
        },
        {
          title: 'Liberté, loi et société',
          blocks: [
            {
              kind: 'definition',
              term: 'Rousseau (1712-1778)',
              definition:
                "« L'homme est né libre, et partout il est dans les fers » (Du contrat social, 1762). La liberté civile consiste dans « l'obéissance à la loi qu'on s'est prescrite ».",
            },
            {
              kind: 'definition',
              term: 'Kant (1724-1804)',
              definition: "La liberté est autonomie : la raison se donne à elle-même la loi morale.",
            },
            {
              kind: 'definition',
              term: 'Montesquieu',
              definition: "« La liberté est le droit de faire tout ce que les lois permettent » (De l'esprit des lois, 1748).",
            },
            {
              kind: 'example',
              title: 'Liberté et libération',
              text: "La négritude de Senghor et Césaire, comme les combats de Fanon et de Cabral, pensent la liberté comme libération collective face à la domination coloniale.",
            },
            {
              kind: 'tip',
              text: "Sujet : « La loi s'oppose-t-elle à la liberté ? » Utilise Rousseau et Montesquieu pour montrer que la loi peut être la condition de la liberté.",
            },
          ],
        },
      ],
      flashcards: [
        { front: '« L’existence précède l’essence »', back: "Sartre, L'existentialisme est un humanisme (1946)." },
        { front: '« L’homme est condamné à être libre »', back: 'Sartre.' },
        { front: '« L’homme est né libre, et partout il est dans les fers »', back: 'Rousseau, Du contrat social (1762).' },
        { front: "La liberté d'indifférence est « le plus bas degré de la liberté »", back: 'Descartes, Méditations métaphysiques (IV).' },
        { front: 'Le libre arbitre est une illusion due à l’ignorance des causes', back: 'Spinoza, Éthique.' },
        { front: '« La liberté est le droit de faire tout ce que les lois permettent »', back: "Montesquieu, De l'esprit des lois." },
        { front: 'Autonomie', back: 'Se donner à soi-même sa propre loi (Kant).' },
        { front: 'Déterminisme / fatalisme', back: 'Tout a une cause / tout est écrit d’avance.' },
      ],
      quiz: [
        {
          id: 'philo-bac-liberte-q1',
          type: 'qcm',
          prompt: '« L’existence précède l’essence » est une thèse de…',
          choices: ['Spinoza', 'Descartes', 'Sartre', 'Rousseau'],
          answer: 2,
          explanation: "Sartre, dans L'existentialisme est un humanisme : l'homme se fait par ses actes.",
        },
        {
          id: 'philo-bac-liberte-q2',
          type: 'vrai-faux',
          prompt: 'Pour Spinoza, le libre arbitre est une illusion.',
          answer: true,
          explanation: "Les hommes se croient libres parce qu'ils ignorent les causes de leurs actions.",
        },
        {
          id: 'philo-bac-liberte-q3',
          type: 'trous',
          prompt: '« L’homme est né ___, et partout il est dans les ___. »',
          answers: ['libre', 'fers'],
          bank: ['libre', 'fers', 'bon', 'chaînes', 'esclave'],
          explanation: 'Première phrase du Contrat social (1762) de Rousseau.',
        },
        {
          id: 'philo-bac-liberte-q4',
          type: 'qcm',
          prompt: "Pour Descartes, la liberté d'indifférence est…",
          choices: [
            'La plus haute forme de liberté',
            'Le plus bas degré de la liberté',
            'Une illusion totale',
            'Le fondement du droit',
          ],
          answer: 1,
          explanation: "Choisir sans raison témoigne d'un défaut de connaissance : c'est « le plus bas degré de la liberté » (IVe Méditation).",
        },
        {
          id: 'philo-bac-liberte-q5',
          type: 'vrai-faux',
          prompt: 'Déterminisme et fatalisme veulent dire la même chose.',
          answer: false,
          explanation: "Le déterminisme dit que tout a des causes, qu'on peut connaître pour agir ; le fatalisme dit que tout arrivera quoi qu'on fasse.",
        },
        {
          id: 'philo-bac-liberte-q6',
          type: 'qcm',
          prompt: "Qui définit la liberté comme « le droit de faire tout ce que les lois permettent » ?",
          choices: ['Hobbes', 'Kant', 'Montesquieu', 'Sartre'],
          answer: 2,
          explanation: "Montesquieu, De l'esprit des lois (1748), livre XI.",
        },
        {
          id: 'philo-bac-liberte-q7',
          type: 'trous',
          prompt: "Se donner à soi-même sa propre loi, c'est l'___ ; recevoir sa loi d'autrui, c'est l'___.",
          answers: ['autonomie', 'hétéronomie'],
          bank: ['autonomie', 'hétéronomie', 'anarchie', 'indifférence', 'nécessité'],
          explanation: "Pour Kant, la volonté morale est autonome : elle obéit à la loi que la raison se donne.",
        },
        {
          id: 'philo-bac-liberte-q8',
          type: 'vrai-faux',
          prompt: "Pour Sartre, dire « je n'avais pas le choix » peut être un acte de mauvaise foi.",
          answer: true,
          explanation: "On choisit toujours, même en refusant de choisir : se cacher derrière les circonstances, c'est fuir sa responsabilité.",
        },
        {
          id: 'philo-bac-liberte-q9',
          type: 'qcm',
          prompt: 'Selon Rousseau, en quoi consiste la liberté civile ?',
          choices: [
            'Faire tout ce que l’on désire',
            'La soumission au plus fort',
            "L'absence totale de lois",
            "L'obéissance à la loi qu'on s'est prescrite",
          ],
          answer: 3,
          explanation: "Dans le Contrat social, le citoyen est libre parce qu'il obéit à la volonté générale, à laquelle il participe.",
        },
      ],
    },

    // ───────────────────────── Chapitre 4 ─────────────────────────
    {
      id: 'philo-bac-travail-technique',
      title: 'Le travail et la technique',
      summary:
        "Par le travail et la technique, l'homme transforme la nature et se transforme lui-même ; mais le travail peut aliéner et la technique menacer.",
      essentials: [
        "Le travail est une activité consciente de transformation de la nature pour satisfaire des besoins.",
        "Hegel : dans la dialectique du maître et de l'esclave, c'est l'esclave qui se libère par le travail.",
        "Marx : le travail est propre à l'homme (l'architecte et l'abeille) mais il est aliéné dans le capitalisme.",
        'Technique : mythe de Prométhée (Platon), homo faber (Bergson), « maîtres et possesseurs de la nature » (Descartes).',
        'Hans Jonas : la puissance technique impose un « principe responsabilité » envers les générations futures.',
      ],
      sections: [
        {
          title: 'Définitions',
          blocks: [
            { kind: 'definition', term: 'Travail', definition: "Activité consciente et pénible par laquelle l'homme transforme la nature pour satisfaire ses besoins." },
            { kind: 'definition', term: 'Technique', definition: "Ensemble des procédés et outils permettant d'obtenir un résultat ; savoir-faire efficace." },
            { kind: 'definition', term: 'Aliénation', definition: "Fait de devenir étranger à soi-même, dépossédé de son activité et de son produit." },
            { kind: 'definition', term: 'Division du travail', definition: 'Répartition des tâches entre travailleurs ; Adam Smith (1776) en montre les gains de productivité.' },
          ],
        },
        {
          title: 'Le travail : malédiction ou libération ?',
          blocks: [
            {
              kind: 'list',
              title: 'Un travail dévalorisé',
              items: [
                "Genèse : « C'est à la sueur de ton visage que tu mangeras du pain » : le travail comme punition.",
                "Aristote : l'esclave est un « instrument animé » ; le citoyen libre se consacre à la politique et à la pensée.",
              ],
            },
            {
              kind: 'definition',
              term: "Hegel, la dialectique du maître et de l'esclave",
              definition:
                "Phénoménologie de l'esprit (1807) : le maître dépend du travail de l'esclave ; l'esclave, en transformant la nature, prend conscience de lui-même et se libère.",
            },
            {
              kind: 'definition',
              term: 'Marx (1818-1883)',
              definition:
                "Le travail distingue l'homme : « ce qui distingue dès l'abord le plus mauvais architecte de l'abeille la plus experte, c'est qu'il a construit la cellule dans sa tête avant de la construire dans la ruche » (Le Capital). Mais dans le capitalisme, l'ouvrier est aliéné.",
            },
            {
              kind: 'example',
              title: 'Hannah Arendt',
              text: "Dans Condition de l'homme moderne (1958), elle distingue le travail (besoins vitaux), l'œuvre (fabriquer des objets durables) et l'action (vie politique).",
            },
          ],
        },
        {
          title: 'La technique',
          blocks: [
            {
              kind: 'definition',
              term: 'Le mythe de Prométhée',
              definition:
                "Raconté par Platon (Protagoras) : Épiméthée a distribué toutes les qualités aux animaux et oublié l'homme ; Prométhée vole le feu et le savoir technique pour le lui donner. La technique compense la faiblesse naturelle de l'homme.",
            },
            {
              kind: 'definition',
              term: 'Homo faber',
              definition: "Bergson (L'Évolution créatrice, 1907) : l'intelligence humaine est d'abord la faculté de fabriquer des outils.",
            },
            {
              kind: 'definition',
              term: 'Descartes',
              definition: "La science doit nous rendre « comme maîtres et possesseurs de la nature » (Discours de la méthode, VIe partie).",
            },
            {
              kind: 'warning',
              text: "Descartes dit « COMME maîtres et possesseurs » : ce « comme » marque une limite, l'homme n'est pas Dieu.",
            },
          ],
        },
        {
          title: 'Les dangers de la technique',
          blocks: [
            {
              kind: 'list',
              items: [
                'Destruction de l’environnement, armes de destruction massive.',
                "Déshumanisation du travail (taylorisme, travail à la chaîne).",
                'Dépendance aux machines et au numérique.',
              ],
            },
            {
              kind: 'definition',
              term: 'Hans Jonas',
              definition:
                "Le Principe responsabilité (1979) : « Agis de façon que les effets de ton action soient compatibles avec la permanence d'une vie authentiquement humaine sur terre. »",
            },
            {
              kind: 'tip',
              text: "Sujet : « La technique libère-t-elle l'homme ? » I. Elle le libère des contraintes naturelles (Prométhée, Descartes) ; II. Elle peut l'asservir (aliénation, dangers) ; III. Tout dépend de l'usage responsable qu'on en fait (Jonas).",
            },
          ],
        },
      ],
      flashcards: [
        { front: "Dialectique du maître et de l'esclave", back: "Hegel, Phénoménologie de l'esprit (1807) : l'esclave se libère par le travail." },
        { front: "L'architecte et l'abeille", back: 'Marx, Le Capital : l’homme conçoit son ouvrage avant de le réaliser.' },
        { front: 'Aliénation (Marx)', back: "L'ouvrier est dépossédé de son travail et de son produit." },
        { front: 'Homo faber', back: "Bergson : l'homme est d'abord fabricateur d'outils." },
        { front: '« Comme maîtres et possesseurs de la nature »', back: 'Descartes, Discours de la méthode (VIe partie).' },
        { front: 'Mythe de Prométhée', back: "Platon, Protagoras : le feu et la technique compensent la faiblesse de l'homme." },
        { front: 'Le Principe responsabilité', back: 'Hans Jonas (1979) : responsabilité envers les générations futures.' },
        { front: 'Travail, œuvre, action', back: "Hannah Arendt, Condition de l'homme moderne (1958)." },
      ],
      quiz: [
        {
          id: 'philo-bac-travail-technique-q1',
          type: 'qcm',
          prompt: "Qui a développé la dialectique du maître et de l'esclave ?",
          choices: ['Hegel', 'Marx', 'Aristote', 'Rousseau'],
          answer: 0,
          explanation: "Hegel, dans la Phénoménologie de l'esprit (1807). Marx s'en inspirera.",
        },
        {
          id: 'philo-bac-travail-technique-q2',
          type: 'vrai-faux',
          prompt: "La comparaison entre l'architecte et l'abeille est de Marx.",
          answer: true,
          explanation: "Dans Le Capital : l'architecte construit la cellule dans sa tête avant de la construire réellement, ce que l'abeille ne fait pas.",
        },
        {
          id: 'philo-bac-travail-technique-q3',
          type: 'trous',
          prompt: 'Pour Descartes, la science doit nous rendre « comme ___ et ___ de la nature ».',
          answers: ['maîtres', 'possesseurs'],
          bank: ['maîtres', 'possesseurs', 'esclaves', 'gardiens', 'enfants'],
          explanation: 'Discours de la méthode, VIe partie (1637).',
        },
        {
          id: 'philo-bac-travail-technique-q4',
          type: 'qcm',
          prompt: "Quel philosophe a défini l'homme comme homo faber ?",
          choices: ['Bergson', 'Heidegger', 'Platon', 'Kant'],
          answer: 0,
          explanation: "Bergson, L'Évolution créatrice (1907) : l'intelligence est d'abord faculté de fabriquer des outils.",
        },
        {
          id: 'philo-bac-travail-technique-q5',
          type: 'vrai-faux',
          prompt: 'Le Principe responsabilité est un ouvrage de Karl Marx.',
          answer: false,
          explanation: "C'est un ouvrage de Hans Jonas (1979), qui réfléchit à la responsabilité face à la puissance technique.",
        },
        {
          id: 'philo-bac-travail-technique-q6',
          type: 'qcm',
          prompt: 'Dans quel dialogue de Platon trouve-t-on le mythe de Prométhée et Épiméthée ?',
          choices: ['La République', 'Le Phédon', 'Le Protagoras', 'Le Banquet'],
          answer: 2,
          explanation: "Le mythe est raconté par Protagoras dans le dialogue du même nom.",
        },
        {
          id: 'philo-bac-travail-technique-q7',
          type: 'trous',
          prompt: "Pour Marx, dans le capitalisme, le travail est ___ ; pour Aristote, l'esclave est un « instrument ___ ».",
          answers: ['aliéné', 'animé'],
          bank: ['aliéné', 'animé', 'libéré', 'inerte', 'divin'],
          explanation: "Marx dénonce l'aliénation de l'ouvrier ; Aristote (Politique) voit l'esclave comme un instrument animé.",
        },
        {
          id: 'philo-bac-travail-technique-q8',
          type: 'vrai-faux',
          prompt: "Hannah Arendt distingue le travail, l'œuvre et l'action.",
          answer: true,
          explanation: "Condition de l'homme moderne (1958) : trois activités fondamentales de la vie humaine.",
        },
        {
          id: 'philo-bac-travail-technique-q9',
          type: 'qcm',
          prompt: "Selon Hegel, qui finit par se libérer dans la relation maître-esclave ?",
          choices: ['Le maître, grâce à sa domination', "L'esclave, grâce à son travail", 'Aucun des deux'],
          answer: 1,
          explanation: "En transformant la nature, l'esclave se forme et prend conscience de lui-même ; le maître devient dépendant.",
        },
      ],
    },

    // ───────────────────────── Chapitre 5 ─────────────────────────
    {
      id: 'philo-bac-etat-droit-justice',
      title: "L'État, le droit et la justice",
      summary:
        "Pourquoi obéir à l'État ? Les théories du contrat social fondent le pouvoir sur le consentement ; le droit et la justice doivent limiter la force.",
      essentials: [
        "Weber : l'État détient le « monopole de la violence physique légitime ».",
        "Hobbes : l'état de nature est une guerre de chacun contre chacun ; les hommes se soumettent à un souverain absolu (Léviathan, 1651).",
        'Rousseau : le contrat social fonde la souveraineté du peuple et la volonté générale (1762).',
        "Montesquieu : « le pouvoir arrête le pouvoir » (séparation des pouvoirs) ; Marx : l'État est un instrument de domination de classe.",
        "Droit naturel / droit positif ; légal / légitime. Pascal : la justice sans la force est impuissante.",
      ],
      sections: [
        {
          title: 'Définitions',
          blocks: [
            { kind: 'definition', term: 'État', definition: "Institution qui exerce le pouvoir politique souverain sur une population et un territoire." },
            { kind: 'definition', term: 'Droit positif', definition: "Ensemble des lois effectivement en vigueur dans une société donnée." },
            { kind: 'definition', term: 'Droit naturel', definition: "Droits que l'homme posséderait par nature, avant toute loi écrite (vie, liberté…)." },
            { kind: 'definition', term: 'Légalité / légitimité', definition: 'Conformité à la loi / conformité à ce qui est juste. Une loi peut être légale sans être légitime.' },
            {
              kind: 'example',
              title: 'Antigone (Sophocle)',
              text: "Antigone enterre son frère malgré l'interdiction du roi Créon, au nom des lois non écrites : conflit entre loi de la cité et justice supérieure.",
            },
          ],
        },
        {
          title: "Les fondements de l'État",
          blocks: [
            { kind: 'definition', term: 'Aristote', definition: "« L'homme est par nature un animal politique » (Politique) : la cité est naturelle." },
            {
              kind: 'definition',
              term: 'Hobbes (1588-1679)',
              definition:
                "Léviathan (1651) : à l'état de nature règne la « guerre de chacun contre chacun ». Par peur de la mort, les hommes transfèrent leur pouvoir à un souverain absolu qui assure la paix.",
            },
            {
              kind: 'definition',
              term: 'Locke (1632-1704)',
              definition: "L'État doit protéger les droits naturels (vie, liberté, propriété) ; s'il les viole, le peuple peut résister.",
            },
            {
              kind: 'definition',
              term: 'Rousseau (1712-1778)',
              definition:
                "Du contrat social (1762) : chacun s'unit à tous et n'obéit qu'à la volonté générale. Le peuple est souverain.",
            },
            {
              kind: 'warning',
              text: "« L'homme est un loup pour l'homme » est une formule du poète latin Plaute, reprise par Hobbes (Du citoyen) : ne l'attribue pas à Rousseau, pour qui l'homme est naturellement bon.",
            },
          ],
        },
        {
          title: 'Limiter et critiquer le pouvoir',
          blocks: [
            {
              kind: 'definition',
              term: 'Montesquieu (1689-1755)',
              definition:
                "De l'esprit des lois (1748) : pour éviter l'abus de pouvoir, « il faut que, par la disposition des choses, le pouvoir arrête le pouvoir » : séparation des pouvoirs législatif, exécutif et judiciaire.",
            },
            {
              kind: 'definition',
              term: 'Max Weber (1864-1920)',
              definition: "L'État est la communauté humaine qui revendique avec succès le « monopole de la violence physique légitime » (Le Savant et le politique, 1919).",
            },
            {
              kind: 'definition',
              term: 'Marx et Engels',
              definition: "L'État est un instrument de domination de la classe dominante ; dans la société sans classes, il doit dépérir.",
            },
            {
              kind: 'example',
              title: 'Penseurs africains',
              text: "Nkrumah : « Cherchez d'abord le royaume politique » (l'indépendance politique avant tout). Le juriste sénégalais Kéba Mbaye a théorisé dès 1972 le « droit au développement ». La Charte africaine des droits de l'homme et des peuples est adoptée en 1981.",
            },
          ],
        },
        {
          title: 'La justice',
          blocks: [
            {
              kind: 'definition',
              term: 'Aristote',
              definition: 'Justice commutative (égalité stricte dans les échanges) et justice distributive (répartition proportionnelle au mérite).',
            },
            {
              kind: 'definition',
              term: 'Pascal',
              definition:
                "« La justice sans la force est impuissante ; la force sans la justice est tyrannique » (Pensées).",
            },
            {
              kind: 'definition',
              term: 'Rousseau',
              definition:
                "« Le plus fort n'est jamais assez fort pour être toujours le maître, s'il ne transforme sa force en droit et l'obéissance en devoir » : la force ne fait pas le droit.",
            },
            {
              kind: 'definition',
              term: 'John Rawls',
              definition: "Théorie de la justice (1971) : les principes justes sont ceux qu'on choisirait sous un « voile d'ignorance ».",
            },
            {
              kind: 'tip',
              text: "Sujet : « Doit-on toujours obéir aux lois ? » I. Oui : la loi garantit la paix et la liberté (Hobbes, Rousseau) ; II. Non : une loi injuste n'est pas légitime (Antigone, Locke) ; III. Désobéir de façon justifiée et responsable (désobéissance civile).",
            },
          ],
        },
      ],
      flashcards: [
        { front: "« Monopole de la violence physique légitime »", back: "Max Weber, définition de l'État." },
        { front: 'Léviathan (1651)', back: 'Hobbes : souverain absolu pour sortir de la guerre de chacun contre chacun.' },
        { front: 'Volonté générale', back: 'Rousseau, Du contrat social (1762).' },
        { front: '« Le pouvoir arrête le pouvoir »', back: "Montesquieu, De l'esprit des lois (1748) : séparation des pouvoirs." },
        { front: "« L'homme est par nature un animal politique »", back: 'Aristote, Politique.' },
        { front: '« La justice sans la force est impuissante ; la force sans la justice est tyrannique »', back: 'Pascal, Pensées.' },
        { front: 'Légal / légitime', back: 'Conforme à la loi / conforme à la justice.' },
        { front: 'Droit au développement', back: 'Notion théorisée par le juriste sénégalais Kéba Mbaye (1972).' },
        { front: 'Voile d’ignorance', back: 'John Rawls, Théorie de la justice (1971).' },
      ],
      quiz: [
        {
          id: 'philo-bac-etat-droit-justice-q1',
          type: 'qcm',
          prompt: "Qui définit l'État par le « monopole de la violence physique légitime » ?",
          choices: ['Hobbes', 'Max Weber', 'Marx', 'Rousseau'],
          answer: 1,
          explanation: 'Max Weber, Le Savant et le politique (1919).',
        },
        {
          id: 'philo-bac-etat-droit-justice-q2',
          type: 'vrai-faux',
          prompt: "Pour Hobbes, l'état de nature est une guerre de chacun contre chacun.",
          answer: true,
          explanation: "Dans le Léviathan, l'absence de pouvoir commun rend la vie dangereuse ; d'où le contrat qui crée le souverain.",
        },
        {
          id: 'philo-bac-etat-droit-justice-q3',
          type: 'trous',
          prompt: "Pour Rousseau, le citoyen obéit à la volonté ___ ; pour Montesquieu, « le pouvoir arrête le ___ ».",
          answers: ['générale', 'pouvoir'],
          bank: ['générale', 'pouvoir', 'particulière', 'peuple', 'roi'],
          explanation: "La volonté générale vise l'intérêt commun ; la séparation des pouvoirs empêche l'abus de pouvoir.",
        },
        {
          id: 'philo-bac-etat-droit-justice-q4',
          type: 'qcm',
          prompt: "Qui a écrit « L'homme est par nature un animal politique » ?",
          choices: ['Platon', 'Hobbes', 'Aristote', 'Locke'],
          answer: 2,
          explanation: "Aristote, Politique, livre I : l'homme ne s'accomplit que dans la cité.",
        },
        {
          id: 'philo-bac-etat-droit-justice-q5',
          type: 'vrai-faux',
          prompt: "La formule « l'homme est un loup pour l'homme » est de Rousseau.",
          answer: false,
          explanation: "Elle vient du poète latin Plaute et a été reprise par Hobbes. Rousseau pense au contraire que l'homme est naturellement bon.",
        },
        {
          id: 'philo-bac-etat-droit-justice-q6',
          type: 'qcm',
          prompt: "Pour Marx et Engels, l'État est…",
          choices: [
            "L'expression de la volonté générale",
            'Une institution naturelle',
            'La réalisation de la raison',
            "Un instrument de domination de la classe dominante",
          ],
          answer: 3,
          explanation: "Ils voient l'État comme un outil au service de la classe qui possède les moyens de production.",
        },
        {
          id: 'philo-bac-etat-droit-justice-q7',
          type: 'trous',
          prompt: "Le droit ___ est l'ensemble des lois en vigueur ; le droit ___ désigne les droits que l'homme posséderait par nature.",
          answers: ['positif', 'naturel'],
          bank: ['positif', 'naturel', 'divin', 'coutumier', 'pénal'],
          explanation: 'Cette distinction permet de juger une loi en vigueur au nom de principes supérieurs.',
        },
        {
          id: 'philo-bac-etat-droit-justice-q8',
          type: 'vrai-faux',
          prompt: "Pascal a écrit : « La justice sans la force est impuissante ; la force sans la justice est tyrannique. »",
          answer: true,
          explanation: 'Pensées : il faut mettre ensemble la justice et la force.',
        },
        {
          id: 'philo-bac-etat-droit-justice-q9',
          type: 'qcm',
          prompt: 'Quel juriste sénégalais a théorisé le « droit au développement » ?',
          choices: ['Kéba Mbaye', 'Cheikh Anta Diop', 'Mamadou Dia', 'Abdoulaye Wade'],
          answer: 0,
          explanation: "Kéba Mbaye, magistrat et juriste, a formulé cette notion en 1972 ; elle sera reconnue par l'ONU en 1986.",
        },
        {
          id: 'philo-bac-etat-droit-justice-q10',
          type: 'qcm',
          prompt: "Quelle œuvre met en scène le conflit entre la loi de la cité et les lois non écrites ?",
          choices: ['Antigone de Sophocle', 'Le Léviathan', 'Le Prince de Machiavel', 'La République de Platon'],
          answer: 0,
          explanation: "Antigone désobéit à Créon pour enterrer son frère : c'est l'exemple classique du conflit entre légalité et légitimité.",
        },
      ],
    },

    // ───────────────────────── Chapitre 6 ─────────────────────────
    {
      id: 'philo-bac-verite-raison',
      title: 'La vérité et la raison',
      summary:
        "Comment atteindre la vérité ? Par la raison, le doute et la méthode, en se méfiant de l'opinion ; la raison elle-même a plusieurs formes, comme l'ont montré aussi des penseurs africains.",
      essentials: [
        "La vérité est l'accord de la pensée avec le réel (adéquation) ou de la pensée avec elle-même (cohérence).",
        "Platon (allégorie de la caverne) : il faut sortir de l'opinion pour atteindre le vrai.",
        'Descartes : le doute méthodique et l’évidence ; la raison (le bon sens) est « la chose du monde la mieux partagée ».',
        "Bachelard (obstacle épistémologique) et Popper (réfutabilité) : la science progresse en corrigeant ses erreurs.",
        "Senghor (raison intuitive), Cheikh Anta Diop, Hountondji et Souleymane Bachir Diagne enrichissent le débat sur la raison et l'universel.",
      ],
      sections: [
        {
          title: 'Définitions',
          blocks: [
            { kind: 'definition', term: 'Vérité', definition: "Qualité d'un jugement conforme à la réalité (vérité matérielle) ou cohérent avec lui-même (vérité formelle)." },
            { kind: 'definition', term: 'Réalité', definition: "Ce qui existe. Le réel n'est ni vrai ni faux : seul un jugement sur le réel peut l'être." },
            { kind: 'definition', term: 'Raison', definition: 'Faculté de penser, de juger, de raisonner et de distinguer le vrai du faux.' },
            { kind: 'definition', term: 'Opinion (doxa)', definition: 'Croyance non fondée sur une démonstration ; s’oppose à la science (épistémè).' },
            {
              kind: 'warning',
              text: 'Ne confonds pas vérité et réalité : « la table existe » relève de la réalité ; « la table est en bois » est un jugement qui peut être vrai ou faux.',
            },
          ],
        },
        {
          title: "Sortir de l'opinion",
          blocks: [
            {
              kind: 'definition',
              term: 'Platon, allégorie de la caverne',
              definition:
                "République, livre VII : des prisonniers enchaînés prennent des ombres pour la réalité. Le philosophe est celui qui sort de la caverne vers la lumière du vrai.",
            },
            {
              kind: 'definition',
              term: 'Descartes, le doute méthodique',
              definition:
                "Douter de tout ce qui n'est pas absolument certain pour trouver une vérité indubitable (le cogito). Règle de l'évidence : n'accepter que les idées claires et distinctes.",
            },
            {
              kind: 'text',
              text: "« Le bon sens est la chose du monde la mieux partagée » (Descartes, Discours de la méthode) : tous les hommes ont la raison, mais il faut bien l'appliquer, avec méthode.",
            },
          ],
        },
        {
          title: 'Raison, expérience et science',
          blocks: [
            {
              kind: 'list',
              items: [
                'Rationalisme (Descartes, Leibniz) : la connaissance vient d’abord de la raison.',
                'Empirisme (Locke, Hume) : toute connaissance vient de l’expérience.',
                "Kant : « Des pensées sans contenu sont vides, des intuitions sans concepts sont aveugles » : il faut les deux.",
              ],
            },
            {
              kind: 'definition',
              term: 'Bachelard (1884-1962)',
              definition:
                "La Formation de l'esprit scientifique (1938) : la science se construit contre l'opinion et les « obstacles épistémologiques ». « L'opinion pense mal ; elle ne pense pas. »",
            },
            {
              kind: 'definition',
              term: 'Popper (1902-1994)',
              definition: "Une théorie est scientifique si elle est réfutable (falsifiable), c'est-à-dire si l'expérience peut la contredire.",
            },
            {
              kind: 'definition',
              term: 'Pascal',
              definition: "« Le cœur a ses raisons que la raison ne connaît point » (Pensées) : la raison a des limites.",
            },
          ],
        },
        {
          title: 'La raison et la pensée africaine',
          blocks: [
            {
              kind: 'definition',
              term: 'Senghor et la négritude',
              definition:
                "Pour Senghor, la négritude est l'ensemble des valeurs culturelles du monde noir. Il oppose une « raison-œil » analytique à une « raison-étreinte » intuitive, et écrit : « L'émotion est nègre, comme la raison hellène ». Formule très critiquée pour son essentialisme.",
            },
            {
              kind: 'definition',
              term: 'Cheikh Anta Diop (1923-1986)',
              definition:
                "Nations nègres et culture (1954) : l'Égypte ancienne était négro-africaine et a contribué à la science et à la philosophie grecques. Il applique la méthode scientifique (laboratoire de datation au carbone 14 de l'IFAN).",
            },
            {
              kind: 'definition',
              term: 'Paulin Hountondji',
              definition:
                "Sur la « philosophie africaine » (1976) : il critique l'« ethnophilosophie » (comme La Philosophie bantoue de Tempels) qui prête à tout un peuple une philosophie collective et implicite.",
            },
            {
              kind: 'definition',
              term: 'Souleymane Bachir Diagne',
              definition:
                "Philosophe sénégalais (spécialiste de Bergson, de la logique et de la pensée de Muhammad Iqbal) ; il défend un « universel latéral », construit par le dialogue et la traduction entre les cultures.",
            },
            {
              kind: 'tip',
              text: "Sujet : « La raison peut-elle tout connaître ? » I. La raison est l'instrument du vrai (Platon, Descartes) ; II. Ses limites (Pascal, Kant) ; III. Une raison plurielle et critique (Bachelard, Senghor, Diagne).",
            },
          ],
        },
      ],
      flashcards: [
        { front: 'Allégorie de la caverne', back: 'Platon, République, livre VII.' },
        { front: '« Le bon sens est la chose du monde la mieux partagée »', back: 'Descartes, Discours de la méthode.' },
        { front: '« Le cœur a ses raisons que la raison ne connaît point »', back: 'Pascal, Pensées.' },
        { front: 'Obstacle épistémologique', back: "Bachelard, La Formation de l'esprit scientifique (1938)." },
        { front: 'Réfutabilité (falsifiabilité)', back: 'Popper : critère de scientificité.' },
        { front: '« L’émotion est nègre, comme la raison hellène »', back: 'Senghor ; formule critiquée pour son essentialisme.' },
        { front: 'Nations nègres et culture', back: 'Cheikh Anta Diop (1954).' },
        { front: 'Critique de l’ethnophilosophie', back: 'Paulin Hountondji, Sur la « philosophie africaine » (1976).' },
        { front: 'Universel latéral', back: 'Notion défendue par Souleymane Bachir Diagne : un universel construit par le dialogue des cultures.' },
      ],
      quiz: [
        {
          id: 'philo-bac-verite-raison-q1',
          type: 'qcm',
          prompt: "Dans quel ouvrage de Platon se trouve l'allégorie de la caverne ?",
          choices: ['La République', 'Le Banquet', 'Le Phèdre', 'Le Gorgias'],
          answer: 0,
          explanation: "Elle ouvre le livre VII de la République.",
        },
        {
          id: 'philo-bac-verite-raison-q2',
          type: 'vrai-faux',
          prompt: '« Le cœur a ses raisons que la raison ne connaît point » est une phrase de Descartes.',
          answer: false,
          explanation: 'Elle est de Pascal (Pensées). Descartes a écrit : « Le bon sens est la chose du monde la mieux partagée ».',
        },
        {
          id: 'philo-bac-verite-raison-q3',
          type: 'trous',
          prompt: "Pour Popper, une théorie est scientifique si elle est ___ ; pour Bachelard, la science se construit contre les obstacles ___.",
          answers: ['réfutable', 'épistémologiques'],
          bank: ['réfutable', 'épistémologiques', 'irréfutable', 'métaphysiques', 'évidente'],
          explanation: "Une théorie que rien ne peut contredire n'est pas scientifique (Popper) ; l'opinion et l'expérience immédiate sont des obstacles (Bachelard).",
        },
        {
          id: 'philo-bac-verite-raison-q4',
          type: 'qcm',
          prompt: 'Qui a écrit Nations nègres et culture (1954) ?',
          choices: ['Léopold Sédar Senghor', 'Aimé Césaire', 'Cheikh Anta Diop', 'Kwame Nkrumah'],
          answer: 2,
          explanation: "Cheikh Anta Diop y soutient le caractère négro-africain de l'Égypte ancienne.",
        },
        {
          id: 'philo-bac-verite-raison-q5',
          type: 'vrai-faux',
          prompt: "Paulin Hountondji critique l'ethnophilosophie.",
          answer: true,
          explanation: "Il refuse l'idée d'une philosophie collective et implicite propre à tout un peuple, et défend une philosophie critique et individuelle.",
        },
        {
          id: 'philo-bac-verite-raison-q6',
          type: 'qcm',
          prompt: "Quelle formule de Senghor a été critiquée pour son essentialisme ?",
          choices: [
            "« L'émotion est nègre, comme la raison hellène »",
            '« Je pense, donc je suis »',
            "« L'opinion pense mal »",
            '« Cherchez d’abord le royaume politique »',
          ],
          answer: 0,
          explanation: "On lui reproche d'enfermer les peuples dans des natures figées. La dernière formule est de Nkrumah, la troisième de Bachelard.",
        },
        {
          id: 'philo-bac-verite-raison-q7',
          type: 'trous',
          prompt: 'Le ___ affirme que la connaissance vient de la raison ; l’___ affirme qu’elle vient de l’expérience.',
          answers: ['rationalisme', 'empirisme'],
          bank: ['rationalisme', 'empirisme', 'scepticisme', 'dogmatisme', 'existentialisme'],
          explanation: 'Descartes et Leibniz sont rationalistes ; Locke et Hume empiristes. Kant cherche à les concilier.',
        },
        {
          id: 'philo-bac-verite-raison-q8',
          type: 'vrai-faux',
          prompt: "La réalité et la vérité sont exactement la même chose.",
          answer: false,
          explanation: "La réalité est ce qui existe ; la vérité est la qualité d'un jugement sur cette réalité.",
        },
        {
          id: 'philo-bac-verite-raison-q9',
          type: 'qcm',
          prompt: 'Quel philosophe sénégalais défend un « universel latéral » fondé sur le dialogue et la traduction ?',
          choices: ['Souleymane Bachir Diagne', 'Cheikh Hamidou Kane', 'Mamadou Dia', 'Ousmane Sembène'],
          answer: 0,
          explanation: "Souleymane Bachir Diagne, spécialiste de Bergson et d'Iqbal, enseigne notamment à l'université Columbia (New York).",
        },
        {
          id: 'philo-bac-verite-raison-q10',
          type: 'qcm',
          prompt: '« Des pensées sans contenu sont vides, des intuitions sans concepts sont aveugles. » Qui l’a écrit ?',
          choices: ['Hume', 'Kant', 'Hegel', 'Descartes'],
          answer: 1,
          explanation: "Kant, Critique de la raison pure : la connaissance exige à la fois l'expérience (intuitions) et l'entendement (concepts).",
        },
      ],
    },
  ],
};

export default subject;
