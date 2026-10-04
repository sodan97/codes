import { Redirect, Tabs } from 'expo-router';
import { ActivityIndicator, Text, View, type ColorValue } from 'react-native';

import { Button, MAX_FONT_SCALE, useUi } from '../../components/ui';
import { useProgress } from '../../state/progress';
import { useTheme } from '../../state/theme';

/** Emoji décoratif : TalkBack lit seulement le nom de l'onglet. */
function TabIcon({ emoji, focused }: { emoji: string; focused: boolean }) {
  return (
    <Text
      importantForAccessibility="no-hide-descendants"
      accessibilityElementsHidden
      maxFontSizeMultiplier={MAX_FONT_SCALE}
      style={{ fontSize: 22, opacity: focused ? 1 : 0.5 }}
    >
      {emoji}
    </Text>
  );
}

/** Nom de l'onglet, agrandi avec la police du téléphone sans déborder de la barre. */
function TabLabel({ color, position, children }: { color: ColorValue; position: 'below-icon' | 'beside-icon'; children: string }) {
  return (
    <Text
      maxFontSizeMultiplier={MAX_FONT_SCALE}
      numberOfLines={1}
      style={[{ color, fontWeight: '700', fontSize: 11, textAlign: 'center' }, position === 'beside-icon' && { marginLeft: 16, fontSize: 13 }]}
    >
      {children}
    </Text>
  );
}

export default function TabsLayout() {
  const { loaded, loadError, retryLoad, state } = useProgress();
  const { colors } = useTheme();
  const ui = useUi();
  if (!loaded) {
    return (
      <View style={{ flex: 1, alignItems: 'center', justifyContent: 'center', backgroundColor: colors.bg }}>
        <ActivityIndicator color={colors.primary} />
      </View>
    );
  }
  if (loadError) {
    // Surtout pas d'onboarding ici : la sauvegarde existe peut-être et ne doit pas être écrasée.
    return (
      <View style={{ flex: 1, alignItems: 'center', justifyContent: 'center', gap: 16, padding: 24, backgroundColor: colors.bg }}>
        <Text style={{ fontSize: 48 }}>💾</Text>
        <Text style={[ui.h2, { textAlign: 'center' }]}>Impossible de lire ta progression. Elle n’a pas été effacée.</Text>
        <Button label="Réessayer" onPress={retryLoad} style={{ alignSelf: 'stretch' }} />
      </View>
    );
  }
  if (!state.profile) return <Redirect href="/onboarding" />;

  return (
    <Tabs
      screenOptions={{
        headerShown: false,
        tabBarActiveTintColor: colors.primary,
        tabBarInactiveTintColor: colors.muted,
        tabBarLabel: TabLabel,
        tabBarStyle: { backgroundColor: colors.card, borderTopColor: colors.border },
      }}
    >
      <Tabs.Screen name="index" options={{ title: 'Accueil', tabBarIcon: ({ focused }) => <TabIcon emoji="🏠" focused={focused} /> }} />
      <Tabs.Screen name="matieres" options={{ title: 'Réviser', tabBarIcon: ({ focused }) => <TabIcon emoji="📚" focused={focused} /> }} />
      <Tabs.Screen name="defis" options={{ title: 'Défis', tabBarIcon: ({ focused }) => <TabIcon emoji="🎯" focused={focused} /> }} />
      <Tabs.Screen name="profil" options={{ title: 'Profil', tabBarIcon: ({ focused }) => <TabIcon emoji="🏅" focused={focused} /> }} />
    </Tabs>
  );
}
