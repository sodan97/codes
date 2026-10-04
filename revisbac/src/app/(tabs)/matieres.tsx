import { router } from 'expo-router';
import { StyleSheet, Text, View } from 'react-native';

import { Card, ProgressBar, Screen, styles as ui } from '../../components/ui';
import { getSubjects, getTrack } from '../../data/catalog';
import { subjectProgress } from '../../lib/stats';
import { useProgress } from '../../state/progress';
import { colors } from '../../theme';

export default function Subjects() {
  const { state } = useProgress();
  const track = getTrack(state.profile!.track);
  const subjects = getSubjects(track.id);

  return (
    <Screen>
      <Text style={ui.h1}>Réviser</Text>
      <Text style={ui.muted}>
        {track.emoji} {track.label} · {subjects.length} matières. Pour chaque chapitre : une fiche, des flashcards et un quiz.
      </Text>
      {subjects.map((s) => {
        const p = subjectProgress(s, state);
        return (
          <Card key={s.id} onPress={() => router.push({ pathname: '/matiere/[id]', params: { id: s.id } })}>
            <View style={ui.row}>
              <View style={[styles.icon, { backgroundColor: s.color + '1A' }]}>
                <Text style={{ fontSize: 26 }}>{s.icon}</Text>
              </View>
              <View style={{ flex: 1, gap: 4 }}>
                <Text style={styles.name}>{s.name}</Text>
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
    </Screen>
  );
}

const styles = StyleSheet.create({
  icon: { width: 52, height: 52, borderRadius: 14, alignItems: 'center', justifyContent: 'center' },
  name: { fontSize: 17, fontWeight: '800', color: colors.text },
  chevron: { fontSize: 28, color: colors.muted },
});
