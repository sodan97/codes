import { router } from 'expo-router';
import { StyleSheet, Text, View } from 'react-native';

import { Card, Pill, ProgressBar, Screen, styles as ui } from '../../components/ui';
import { getSubjects, getTrack } from '../../data/catalog';
import { visibleSubjects } from '../../lib/selectors';
import { subjectProgress } from '../../lib/stats';
import { useProgress } from '../../state/progress';
import { colors } from '../../theme';

export default function Subjects() {
  const { state } = useProgress();
  const profile = state.profile!;
  const track = getTrack(profile.track);
  const subjects = visibleSubjects(profile);
  const hiddenCount = getSubjects(track.id).length - subjects.length;
  const progress = new Map(subjects.map((s) => [s.id, subjectProgress(s, state)]));

  // Priorité : la matière commencée la moins maîtrisée, dès que deux matières sont commencées.
  const started = subjects.filter((s) => progress.get(s.id)!.read > 0 || s.chapters.some((c) => state.quizBest[c.id] !== undefined));
  const priority =
    started.length >= 2 ? started.reduce((min, s) => (progress.get(s.id)!.mastery < progress.get(min.id)!.mastery ? s : min)).id : null;

  return (
    <Screen>
      <Text style={ui.h1}>Réviser</Text>
      <Text style={ui.muted}>
        {track.emoji} {track.label} · {subjects.length} matières. Pour chaque chapitre : une fiche, des flashcards et un quiz.
      </Text>
      {subjects.map((s) => {
        const p = progress.get(s.id)!;
        return (
          <Card key={s.id} onPress={() => router.push({ pathname: '/matiere/[id]', params: { id: s.id } })}>
            <View style={ui.row}>
              <View style={[styles.icon, { backgroundColor: s.color + '1A' }]}>
                <Text style={{ fontSize: 26 }}>{s.icon}</Text>
              </View>
              <View style={{ flex: 1, gap: 4 }}>
                <View style={[ui.row, { gap: 8, flexWrap: 'wrap' }]}>
                  <Text style={styles.name}>{s.name}</Text>
                  {s.id === priority && <Pill label="Priorité" color={colors.red} />}
                </View>
                <Text style={ui.muted}>
                  {p.read}/{p.total} fiches lues · maîtrise {Math.round(p.mastery * 100)} %
                </Text>
                <ProgressBar value={p.mastery} color={s.color} />
              </View>
              <Text style={styles.chevron}>›</Text>
            </View>
          </Card>
        );
      })}
      {hiddenCount > 0 && (
        <Text style={[ui.muted, { textAlign: 'center' }]} onPress={() => router.navigate('/profil')}>
          {hiddenCount} matière{hiddenCount > 1 ? 's' : ''} masquée{hiddenCount > 1 ? 's' : ''} · modifier dans Profil
        </Text>
      )}
    </Screen>
  );
}

const styles = StyleSheet.create({
  icon: { width: 52, height: 52, borderRadius: 14, alignItems: 'center', justifyContent: 'center' },
  name: { fontSize: 17, fontWeight: '800', color: colors.text },
  chevron: { fontSize: 28, color: colors.muted },
});
