import { Icon, type IconName } from "@/components/Icon";
import { Reveal } from "@/components/Reveal";
import { Stats } from "@/components/Stats";

const drivers: { icon: IconName; title: string; text: string }[] = [
  { icon: "radar", title: "Menaces cyber croissantes", text: "Rançongiciels, fraudes, hameçonnage ciblé : les attaques se professionnalisent." },
  { icon: "compliance", title: "Exigences renforcées", text: "Régulateurs, protection des données et partenaires internationaux élèvent le niveau attendu." },
  { icon: "cloud", title: "Transformation Cloud", text: "Les migrations déplacent le périmètre de sécurité et les responsabilités." },
  { icon: "ai", title: "Accélération de l'IA", text: "L'intelligence artificielle ouvre de nouvelles opportunités et de nouvelles surfaces d'attaque." },
];

export function ChallengeSection() {
  return (
    <section id="enjeux" className="relative bg-sand-50 py-24 sm:py-32" aria-labelledby="enjeux-title">
      <div className="container-site">
        <div className="grid gap-14 lg:grid-cols-[1fr_1.1fr] lg:gap-20">
          <Reveal>
            <p className="eyebrow mb-5 text-gold-700">Le contexte</p>
            <h2 id="enjeux-title" className="text-3xl font-bold leading-[1.15] text-forest-800 sm:text-4xl lg:text-[2.75rem]">
              Les risques numériques évoluent.
              <br />
              Votre stratégie de sécurité aussi.
            </h2>
            <p className="mt-6 text-lg leading-relaxed text-muted">
              Les organisations doivent aujourd&apos;hui composer avec des menaces cyber croissantes, des exigences réglementaires
              renforcées, la transformation Cloud et l&apos;accélération de l&apos;intelligence artificielle.
            </p>
            <p className="mt-6 border-l-2 border-gold-400 pl-5 font-display text-xl font-semibold leading-snug text-forest-800">
              Nous vous aidons à transformer ces enjeux en décisions concrètes.
            </p>
          </Reveal>
          <ul className="grid gap-5 sm:grid-cols-2">
            {drivers.map((driver, index) => (
              <Reveal as="li" key={driver.title} delay={index * 0.08} className="rounded-2xl border border-forest-800/10 bg-white p-6">
                <Icon name={driver.icon} className="h-7 w-7 text-emerald-600" />
                <h3 className="mt-5 text-base font-semibold text-forest-800">{driver.title}</h3>
                <p className="mt-2 text-sm leading-relaxed text-muted">{driver.text}</p>
              </Reveal>
            ))}
          </ul>
        </div>
        <div className="mt-20">
          <Stats
            items={[
              { value: "3", label: "pôles complémentaires : Advisory, Academy, Labs" },
              { value: "6", label: "domaines d'expertise" },
              { value: "8", label: "référentiels internationaux pris en compte" },
              { value: "8", label: "secteurs prioritaires accompagnés" },
            ]}
          />
        </div>
      </div>
    </section>
  );
}
