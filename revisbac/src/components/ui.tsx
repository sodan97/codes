import type { ReactNode } from 'react';
import { Pressable, ScrollView, StyleSheet, Text, View, type StyleProp, type ViewStyle } from 'react-native';
import { SafeAreaView, type Edge } from 'react-native-safe-area-context';

import { colors, radius, shadow } from '../theme';

export function Screen({ children, edges = ['top'], scroll = true }: { children: ReactNode; edges?: Edge[]; scroll?: boolean }) {
  return (
    <SafeAreaView style={styles.screen} edges={edges}>
      {scroll ? (
        <ScrollView contentContainerStyle={styles.screenContent} showsVerticalScrollIndicator={false}>
          {children}
        </ScrollView>
      ) : (
        <View style={[styles.screenContent, { flex: 1 }]}>{children}</View>
      )}
    </SafeAreaView>
  );
}

export function Card({ children, style, onPress }: { children: ReactNode; style?: StyleProp<ViewStyle>; onPress?: () => void }) {
  if (onPress) {
    return (
      <Pressable onPress={onPress} style={({ pressed }) => [styles.card, style, pressed && { opacity: 0.85, transform: [{ scale: 0.99 }] }]}>
        {children}
      </Pressable>
    );
  }
  return <View style={[styles.card, style]}>{children}</View>;
}

type ButtonVariant = 'primary' | 'secondary' | 'ghost' | 'gold';

export function Button({
  label,
  onPress,
  variant = 'primary',
  disabled,
  color,
  style,
}: {
  label: string;
  onPress: () => void;
  variant?: ButtonVariant;
  disabled?: boolean;
  color?: string;
  style?: StyleProp<ViewStyle>;
}) {
  const bg =
    variant === 'primary' ? (color ?? colors.primary) : variant === 'gold' ? colors.gold : variant === 'secondary' ? colors.card : 'transparent';
  const fg = variant === 'primary' ? '#fff' : variant === 'gold' ? colors.text : (color ?? colors.primary);
  return (
    <Pressable
      accessibilityRole="button"
      disabled={disabled}
      onPress={onPress}
      style={({ pressed }) => [
        styles.button,
        { backgroundColor: bg },
        variant === 'secondary' && { borderWidth: 2, borderColor: color ?? colors.primary },
        disabled && { opacity: 0.4 },
        pressed && { opacity: 0.8 },
        style,
      ]}
    >
      <Text style={[styles.buttonLabel, { color: fg }]}>{label}</Text>
    </Pressable>
  );
}

export function ProgressBar({ value, color = colors.primary, height = 8, track = colors.border }: { value: number; color?: string; height?: number; track?: string }) {
  const pct = Math.max(0, Math.min(1, value)) * 100;
  return (
    <View style={{ height, borderRadius: height, backgroundColor: track, overflow: 'hidden' }}>
      <View style={{ width: `${pct}%`, height: '100%', backgroundColor: color, borderRadius: height }} />
    </View>
  );
}

export function Pill({ label, color = colors.primary, bg }: { label: string; color?: string; bg?: string }) {
  return (
    <View style={[styles.pill, { backgroundColor: bg ?? color + '1A' }]}>
      <Text style={[styles.pillText, { color }]}>{label}</Text>
    </View>
  );
}

export function SectionTitle({ children, right }: { children: ReactNode; right?: ReactNode }) {
  return (
    <View style={styles.sectionTitleRow}>
      <Text style={styles.sectionTitle}>{children}</Text>
      {right}
    </View>
  );
}

export const styles = StyleSheet.create({
  screen: { flex: 1, backgroundColor: colors.bg },
  screenContent: { padding: 16, paddingBottom: 40, gap: 14, maxWidth: 720, width: '100%', alignSelf: 'center' },
  card: { backgroundColor: colors.card, borderRadius: radius.md, padding: 16, ...shadow },
  button: { paddingVertical: 14, paddingHorizontal: 18, borderRadius: radius.md, alignItems: 'center', justifyContent: 'center' },
  buttonLabel: { fontSize: 16, fontWeight: '700' },
  pill: { paddingHorizontal: 10, paddingVertical: 4, borderRadius: 999, alignSelf: 'flex-start' },
  pillText: { fontSize: 12, fontWeight: '700' },
  sectionTitleRow: { flexDirection: 'row', alignItems: 'center', justifyContent: 'space-between', marginTop: 6 },
  sectionTitle: { fontSize: 18, fontWeight: '800', color: colors.text },
  h1: { fontSize: 26, fontWeight: '800', color: colors.text },
  h2: { fontSize: 20, fontWeight: '800', color: colors.text },
  body: { fontSize: 15, lineHeight: 22, color: colors.text },
  muted: { fontSize: 13, color: colors.muted },
  row: { flexDirection: 'row', alignItems: 'center', gap: 10 },
});
