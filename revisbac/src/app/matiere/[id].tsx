import { router, Stack, useLocalSearchParams } from 'expo-router';
import { Pressable, StyleSheet, Text, View } from 'react-native';

import { Button, Card, Pill, ProgressBar, Screen, SectionTitle, styles as ui } from '../../components/ui';
import { getSubject } from '../../data/catalog';
import { formatNote, mention, noteTrend } from '../../lib/gamification';
import { EXAM_SIZE } from '../../lib/quizBuilder';
import { subjectProgress } from '../../lib/stats';
import { useProgress } from '../../state/progress';
import { colors, radius } from '../../theme';

export default function SubjectScreen() {
  const { id } = useLocalSearchParams<{ id: string }>();
  const { state } = useProgress();
  const subject = getSubject(id);
  if (!subject) return <Text style={{ padding: 20 }}>Matière introuvable.</Text>;

  const p = subjectProgress(subject, state);
  const examNote = state.examBest[subject.id];
  const examTrend = noteTrend(state.examHistory[subject.id]);
  const nQuestions = subject.chapters.reduce((n, c) => n + c.quiz.length, 0);

  return (
    <Screen edges={['bottom']}>
      <Stack.Screen options={{ title: subject.name }} />
      <Card style={{ backgroundColor: subject.color, gap: 8 }}>
        <View style={ui.row}>
          <Text style={{ fontSize: 40 }}>{subject.icon}</Text>
          <View style={{ flex: 1 }}>
            <Text style={[ui.h2, { color: '#fff' }]}>{subject.name}</Text>
            <Text style={{ color: '#ffffffd0' }}>
              {p.total} chapitres · {nQuestions} questions
            </Text>
          </View>
        </View>
        <ProgressBar value={p.mastery} color={colors.gold} track="#ffffff40" height={10} />
        <Text style={{ color: '#fff', fontWeight: '700' }}>Maîtrise : {Math.round(p.mastery * 100)} %</Text>
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
      {subject.chapters.map((c, i) => {
        const read = !!state.fichesRead[c.id];
        const best = state.quizBest[c.id];
        return (
          <Card key={c.id} style={{ gap: 10 }}>
            <View style={[ui.row, { alignItems: 'flex-start' }]}>
              <View style={[styles.num, { backgroundColor: read ? subject.color : colors.border }]}>
                <Text style={{ color: read ? '#fff' : colors.muted, fontWeight: '800' }}>{read ? '✓' : i + 1}</Text>
              </View>
              <View style={{ flex: 1, gap: 4 }}>
                <Text style={styles.title}>{c.title}</Text>
                <Text style={ui.muted}>{c.summary}</Text>
                {best !== undefined && (
                  <Pill label={`Meilleur quiz : ${best} %`} color={best >= 80 ? colors.primary : best >= 50 ? colors.goldText : colors.red} />
                )}
              </View>
            </View>
            <View style={styles.actions}>
              <Action label="📄 Fiche" color={subject.color} onPress={() => router.push({ pathname: '/fiche/[id]', params: { id: c.id } })} />
              <Action label="🃏 Cartes" color={subject.color} onPress={() => router.push({ pathname: '/flashcards/[id]', params: { id: c.id } })} />
              <Action
                label="✅ Quiz"
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

function Action({ label, color, onPress, filled }: { label: string; color: string; onPress: () => void; filled?: boolean }) {
  return (
    <Pressable
      onPress={onPress}
      style={({ pressed }) => [
        styles.action,
        { borderColor: color, backgroundColor: filled ? color : colors.card },
        pressed && { opacity: 0.7 },
      ]}
    >
      <Text style={{ fontWeight: '700', color: filled ? '#fff' : color }}>{label}</Text>
    </Pressable>
  );
}

const styles = StyleSheet.create({
  title: { fontSize: 16, fontWeight: '800', color: colors.text },
  num: { width: 30, height: 30, borderRadius: 15, alignItems: 'center', justifyContent: 'center' },
  actions: { flexDirection: 'row', gap: 8 },
  action: { flex: 1, borderWidth: 2, borderRadius: radius.sm, paddingVertical: 10, alignItems: 'center' },
});
