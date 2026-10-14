import Link from "next/link";
import { Icon } from "@/components/Icon";
import { Reveal } from "@/components/Reveal";
import { SectionHeader } from "@/components/SectionHeader";

const issues = [
  { question: "Vous devez structurer votre gouvernance cyber ?", answer: "WOOLYS Advisory", href: "/advisory" },
  { question: "Vous devez cartographier vos risques ?", answer: "Risk Management / EBIOS RM", href: "/expertises/risk-management" },
  { question: "Vous devez répondre à une exigence réglementaire ?", answer: "GRC & Compliance", href: "/expertises/grc-compliance" },
  { question: "Vous migrez vers le Cloud ?", answer: "Cloud Security", href: "/expertises/cloud-security" },
  { question: "Vous déployez des solutions d'IA ?", answer: "AI Security", href: "/expertises/ai-security" },
  {
    question: "Vos collaborateurs doivent devenir votre première ligne de défense ?",
    answer: "WOOLYS Academy",
    href: "/academy",
  },
];

export function IssuesSection() {
  return (
    <section className="bg-sand-50 py-24 sm:py-32" aria-labelledby="issues-title">
      <div className="container-site grid gap-14 lg:grid-cols-[0.8fr_1.2fr] lg:gap-20">
        <div className="lg:sticky lg:top-32 lg:self-start">
          <SectionHeader
            id="issues-title"
            eyebrow="Votre situation"
            title="Vos enjeux. Notre expertise."
            intro="Chaque organisation part d'une question concrète. Identifiez la vôtre : nous vous orientons vers la bonne réponse."
          />
        </div>
        <ul className="divide-y divide-forest-800/10 border-y border-forest-800/10">
          {issues.map((issue, index) => (
            <Reveal as="li" key={issue.question} delay={index * 0.05}>
              <Link
                href={issue.href}
                className="group grid gap-3 py-7 transition-colors sm:grid-cols-[1fr_auto] sm:items-center sm:gap-8"
              >
                <span className="font-display text-lg font-semibold leading-snug text-forest-800 sm:text-xl">{issue.question}</span>
                <span className="inline-flex items-center gap-3 text-sm font-semibold text-emerald-600 transition-colors group-hover:text-gold-700">
                  <span className="h-px w-6 bg-current transition-all duration-300 group-hover:w-10" aria-hidden="true" />
                  {issue.answer}
                  <Icon name="arrowRight" className="h-4 w-4 transition-transform duration-300 group-hover:translate-x-1" />
                </span>
              </Link>
            </Reveal>
          ))}
        </ul>
      </div>
    </section>
  );
}
