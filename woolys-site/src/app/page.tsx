import type { Metadata } from "next";
import { CTA } from "@/components/CTA";
import { AfricaSection } from "@/sections/AfricaSection";
import { ApproachSection } from "@/sections/ApproachSection";
import { ChallengeSection } from "@/sections/ChallengeSection";
import { ExpertisesSection } from "@/sections/ExpertisesSection";
import { HomeHero } from "@/sections/HomeHero";
import { InsightsSection } from "@/sections/InsightsSection";
import { IssuesSection } from "@/sections/IssuesSection";
import { MethodologiesSection } from "@/sections/MethodologiesSection";
import { PolesSection } from "@/sections/PolesSection";
import { SectorsSection } from "@/sections/SectorsSection";
import { siteConfig } from "@/lib/site";

export const metadata: Metadata = {
  title: { absolute: `${siteConfig.name} — Cabinet de cybersécurité, GRC et gestion des risques au Sénégal et en Afrique` },
  alternates: { canonical: siteConfig.url },
};

export default function HomePage() {
  return (
    <>
      <HomeHero />
      <ChallengeSection />
      <PolesSection />
      <ExpertisesSection />
      <IssuesSection />
      <MethodologiesSection />
      <SectorsSection />
      <ApproachSection />
      <AfricaSection />
      <InsightsSection />
      <CTA />
    </>
  );
}
