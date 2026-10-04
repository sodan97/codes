import { Stack } from 'expo-router';
import { StatusBar } from 'expo-status-bar';
import { SafeAreaProvider } from 'react-native-safe-area-context';

import { ProgressProvider } from '../state/progress';
import { colors } from '../theme';

export default function RootLayout() {
  return (
    <SafeAreaProvider>
      <ProgressProvider>
        <StatusBar style="dark" />
        <Stack
          screenOptions={{
            headerTintColor: colors.text,
            headerTitleStyle: { fontWeight: '800' },
            headerStyle: { backgroundColor: colors.bg },
            headerShadowVisible: false,
            contentStyle: { backgroundColor: colors.bg },
            headerBackTitle: 'Retour',
          }}
        >
          <Stack.Screen name="(tabs)" options={{ headerShown: false }} />
          <Stack.Screen name="onboarding" options={{ headerShown: false }} />
          <Stack.Screen name="matiere/[id]" options={{ title: '' }} />
          <Stack.Screen name="fiche/[id]" options={{ title: 'Fiche de révision' }} />
          <Stack.Screen name="flashcards/[id]" options={{ title: 'Flashcards' }} />
          <Stack.Screen name="quiz" options={{ title: 'Quiz', gestureEnabled: false }} />
        </Stack>
      </ProgressProvider>
    </SafeAreaProvider>
  );
}
