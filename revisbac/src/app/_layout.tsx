import { DarkTheme, DefaultTheme, ThemeProvider as NavigationThemeProvider, router, Stack, type ErrorBoundaryProps, type Theme } from 'expo-router';
import * as SplashScreen from 'expo-splash-screen';
import { StatusBar } from 'expo-status-bar';
import { useEffect, useMemo } from 'react';
import { Text, View } from 'react-native';
import { SafeAreaProvider } from 'react-native-safe-area-context';

import { Button, useUi } from '../components/ui';
import { ProgressProvider, useProgress } from '../state/progress';
import { ThemeProvider, useTheme } from '../state/theme';

// L'écran de démarrage reste affiché jusqu'à la lecture du thème puis de la progression (voir ThemeProvider et SplashGate).
SplashScreen.preventAutoHideAsync().catch(() => {});

/** Filet de sécurité : une page qui plante (contenu mal formé…) n'emporte pas toute l'application. */
export function ErrorBoundary({ error, retry }: ErrorBoundaryProps) {
  const { colors } = useTheme();
  const ui = useUi();
  useEffect(() => {
    console.error(error);
  }, [error]);

  const home = () => {
    try {
      router.dismissTo('/');
    } catch {
      // Pas de navigation disponible (la mise en page elle-même a planté) : retry suffit.
    }
    void retry();
  };

  return (
    <View style={{ flex: 1, alignItems: 'center', justifyContent: 'center', gap: 16, padding: 24, backgroundColor: colors.bg }}>
      <Text style={{ fontSize: 48 }}>🛠️</Text>
      <Text style={[ui.h2, { textAlign: 'center' }]}>Oups, cette page a un problème. Ta progression est sauvegardée.</Text>
      <Button label="Réessayer" onPress={() => void retry()} style={{ alignSelf: 'stretch' }} />
      <Button label="Retour à l’accueil" variant="secondary" onPress={home} style={{ alignSelf: 'stretch' }} />
    </View>
  );
}

/** Cache l'écran de démarrage une fois la progression chargée. */
function SplashGate() {
  const { loaded } = useProgress();
  useEffect(() => {
    if (loaded) SplashScreen.hideAsync().catch(() => {});
  }, [loaded]);
  return null;
}

export default function RootLayout() {
  return (
    <SafeAreaProvider>
      <ThemeProvider>
        <ProgressProvider>
          <SplashGate />
          <ThemedStack />
        </ProgressProvider>
      </ThemeProvider>
    </SafeAreaProvider>
  );
}

/** Navigation aux couleurs du thème : en-têtes, fond pendant les transitions, barre d'état. */
function ThemedStack() {
  const { colors, scheme } = useTheme();
  const navigationTheme = useMemo<Theme>(() => {
    const base = scheme === 'dark' ? DarkTheme : DefaultTheme;
    return {
      ...base,
      colors: {
        ...base.colors,
        primary: colors.primary,
        background: colors.bg,
        card: colors.card,
        text: colors.text,
        border: colors.border,
        notification: colors.red,
      },
    };
  }, [scheme, colors]);
  return (
    <NavigationThemeProvider value={navigationTheme}>
      <StatusBar style={scheme === 'dark' ? 'light' : 'dark'} />
      <Stack
        // Chaque écran a son propre filet : la navigation reste en place pour revenir à l'accueil.
        unstable_screenErrorBoundary={ErrorBoundary}
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
    </NavigationThemeProvider>
  );
}
