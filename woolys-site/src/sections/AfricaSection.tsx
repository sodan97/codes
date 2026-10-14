import { Reveal } from "@/components/Reveal";
import { AfricaMap } from "./AfricaMap";

const zones = [
  { name: "Sénégal", text: "Notre marché prioritaire et notre ancrage.", color: "bg-gold-400" },
  { name: "Afrique de l'Ouest", text: "Une compréhension des réalités régionales et des cadres communautaires.", color: "bg-emerald-400" },
  { name: "Afrique francophone", text: "Une expertise portée dans la langue de nos clients.", color: "bg-emerald-500" },
  { name: "International", text: "Des standards internationaux, pour des organisations ouvertes sur le monde.", color: "bg-white" },
];

export function AfricaSection() {
  return (
    <section className="relative overflow-hidden bg-forest-950 py-24 text-white sm:py-32" aria-labelledby="africa-title">
      <div className="absolute inset-0 bg-grid opacity-40" />
      <div className="absolute inset-0 bg-[radial-gradient(ellipse_at_30%_50%,rgba(21,143,122,0.25),transparent_55%)]" />
      <div className="container-site relative grid items-center gap-16 lg:grid-cols-2">
        <Reveal className="order-2 lg:order-1">
          <AfricaMap />
        </Reveal>
        <div className="order-1 lg:order-2">
          <Reveal>
            <p className="eyebrow mb-5 text-gold-400">Ancrage</p>
            <h2 id="africa-title" className="text-3xl font-bold leading-[1.15] sm:text-4xl lg:text-[2.75rem]">
              Une expertise internationale.
              <br />
              <span className="text-gradient-gold">Une compréhension des réalités africaines.</span>
            </h2>
            <p className="mt-6 text-lg leading-relaxed text-forest-100/80">
              WOOLYS ambitionne de contribuer à la construction d&apos;un écosystème numérique africain plus sûr, plus résilient et plus
              compétitif.
            </p>
          </Reveal>
          <ul className="mt-10 grid gap-4 sm:grid-cols-2">
            {zones.map((zone, index) => (
              <Reveal as="li" key={zone.name} delay={index * 0.08} className="rounded-2xl border border-white/10 bg-white/[0.03] p-5">
                <p className="flex items-center gap-3 font-display font-semibold">
                  <span className={`h-2.5 w-2.5 rounded-full ${zone.color}`} aria-hidden="true" />
                  {zone.name}
                </p>
                <p className="mt-2 text-sm leading-relaxed text-forest-100/70">{zone.text}</p>
              </Reveal>
            ))}
          </ul>
        </div>
      </div>
    </section>
  );
}
