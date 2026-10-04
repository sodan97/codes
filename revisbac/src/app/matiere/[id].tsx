import { router, Stack, useLocalSearchParams } from 'expo-router';
import { Pressable, StyleSheet, Text, View } from 'react-native';

import { Stars } from '../../components/Stars';
import { Button, Card, MAX_FONT_SCALE, Pill, ProgressBar, Screen, SectionTitle, useUi } from '../../components/ui';
import { getSubject } from '../../data/catalog';
import { dayKey } from '../../lib/dates';
import { formatNote, mention, noteTrend } from '../../lib/gamification';
import { EXAM_SIZE } from '../../lib/quizBuilder';
import { chapterStars, STAR_QUIZ_PERCENT, starsAccessibilityLabel, starsLabel, subjectProgress } from '../../lib/stats';
import { useProgress } from '../../state/progress';
import { useStyles, useTheme } from '../../state/theme';
import { accentText, radius, textOn, type Colors } from '../../theme';

export default function SubjectScreen() {
  const { id } = useLocalSearchParams<{ id: string }>();
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  const ui = useUi();
  const { state } = useProgress();
  const subject = getSubject(id);
  if (!subject) return <Text style={[ui.body, { padding: 20 }]}>Matière introuvable.</Text>;
  // Texte posé sur la couleur de la matière : blanc ou sombre selon le contraste.
  const onSubject = textOn(subject.color);

  const p = subjectProgress(subject, state);
  const examNote = state.examBest[subject.id];
  const examTrend = noteTrend(state.examHistory[subject.id]);
  const nQuestions = subject.chapters.reduce((n, c) => n + c.quiz.length, 0);
  const today = dayKey();

  return (
    <Screen edges={['bottom']}>
      <Stack.Screen options={{ title: subject.name }} />
      <Card style={{ backgroundColor: subject.color, gap: 8 }}>
        <View style={ui.row}>
          <Text style={{ fontSize: 40 }}>{subject.icon}</Text>
          <View style={{ flex: 1 }}>
            <Text style={[ui.h2, { color: onSubject }]}>{subject.name}</Text>
            <Text style={{ color: onSubject, fontSize: 14 }}>
              {p.total} chapitres · {nQuestions} questions
            </Text>
          </View>
        </View>
        <ProgressBar value={p.mastery} color={colors.gold} track={onSubject + '40'} height={10} accessibilityLabel="Maîtrise de la matière" />
        <Text
          style={{ color: onSubject, fontWeight: '700' }}
          accessibilityLabel={`Maîtrise : ${starsAccessibilityLabel(p.stars, p.maxStars)}, ${Math.round(p.mastery * 100)} %`}
        >
          Maîtrise : {starsLabel(p.stars, p.maxStars)} · {Math.round(p.mastery * 100)} %
        </Text>
      </Card>

      <Card style={{ gap: 8 }}>
        <Text style={styles.title}>📝 Examen blanc</Text>
        <Text style={ui.body}>
          {Math.min(EXAM_SIZE, nQuestions)} questions tirées de tous les chapitres, en temps limité, corrigées à la fin. Ta note est sur 20, avec la mention.
        </Text>
        {examNote !== undefined && (
          <Text style={ui.muted}>
            Meilleure note : {formatNote(examNote)}/20 ({mention(examNote)})
          </Text>
        )}
        {examTrend && <Text style={ui.muted}>Tes dernières notes : {examTrend}</Text>}
        <Button label="Lancer l’examen blanc" color={subject.color} onPress={() => router.push({ pathname: '/quiz', params: { mode: 'exam', id: subject.id } })} />
      </Card>

      <SectionTitle>Chapitres</SectionTitle>
      <Text style={ui.muted}>
        ★ fiche lue · ★★ quiz réussi à {STAR_QUIZ_PERCENT} % · ★★★ quiz réussi à {STAR_QUIZ_PERCENT} % deux jours différents et flashcards terminées
      </Text>
      {subject.chapters.map((c, i) => {
        const read = !!state.fichesRead[c.id];
        const best = state.quizBest[c.id];
        const s = chapterStars(c.id, state, today);
        return (
          <Card key={c.id} style={{ gap: 10 }}>
            <View style={[ui.row, { alignItems: 'flex-start' }]}>
              <View
                style={[styles.num, { backgroundColor: read ? subject.color : colors.border }]}
                accessible
                accessibilityLabel={read ? `Chapitre ${i + 1}, fiche lue` : `Chapitre ${i + 1}`}
              >
                <Text style={{ color: read ? onSubject : colors.text, fontWeight: '800' }} maxFontSizeMultiplier={MAX_FONT_SCALE}>
                  {read ? '✓' : i + 1}
                </Text>
              </View>
              <View style={{ flex: 1, gap: 4 }}>
                <Text style={styles.title}>{c.title}</Text>
                <Stars stars={s.stars} faded={s.faded} />
                <Text style={ui.muted}>{c.summary}</Text>
                {best !== undefined && (
                  <Pill label={`Meilleur quiz : ${best} %`} color={best >= STAR_QUIZ_PERCENT ? colors.primary : best >= 50 ? colors.goldText : colors.red} />
                )}
                {s.reminder && <Pill label={s.reminder} color={colors.goldText} />}
                {s.next && <Text style={styles.next}>{s.next}</Text>}
              </View>
            </View>
            <View style={styles.actions}>
              <Action
                label="📄 Fiche"
                accessibilityLabel={`Fiche : ${c.title}`}
                color={subject.color}
                onPress={() => router.push({ pathname: '/fiche/[id]', params: { id: c.id } })}
              />
              <Action
                label="🃏 Cartes"
                accessibilityLabel={`Flashcards : ${c.title}`}
                color={subject.color}
                onPress={() => router.push({ pathname: '/flashcards/[id]', params: { id: c.id } })}
              />
              <Action
                label="✅ Quiz"
                accessibilityLabel={`Quiz : ${c.title}`}
                color={subject.color}
                filled
                onPress={() => router.push({ pathname: '/quiz', params: { mode: 'chapter', id: c.id } })}
              />
            </View>
          </Card>
        );
      })}
    </Screen>
  );
}

function Action({
  label,
  accessibilityLabel,
  color,
  onPress,
  filled,
}: {
  label: string;
  accessibilityLabel: string;
  color: string;
  onPress: () => void;
  filled?: boolean;
}) {
  const { colors } = useTheme();
  const styles = useStyles(makeStyles);
  return (
    <Pressable
      onPress={onPress}
      accessibilityRole="button"
      accessibilityLabel={accessibilityLabel}
      style={({ pressed }) => [
        styles.action,
        { borderColor: color, backgroundColor: filled ? color : colors.card },
        pressed && { opacity: 0.7 },
      ]}
    >
      <Text style={{ fontWeight: '700', color: filled ? textOn(color) : accentText(color, colors) }}>{label}</Text>
    </Pressable>
  );
}

const makeStyles = (colors: Colors) =>
  StyleSheet.create({
    title: { fontSize: 16, fontWeight: '800', color: colors.text },
    next: { fontSize: 13, fontWeight: '700', color: colors.primaryDark },
    num: { minWidth: 30, minHeight: 30, borderRadius: 15, paddingHorizontal: 4, alignItems: 'center', justifyContent: 'center' },
    actions: { flexDirection: 'row', gap: 8 },
    action: { flex: 1, borderWidth: 2, borderRadius: radius.sm, paddingVertical: 10, alignItems: 'center' },
  });
