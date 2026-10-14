import { africaDots, DAKAR, MAP_HEIGHT, MAP_WIDTH } from "@/lib/africa";

const fills = {
  senegal: "#D4AF7C",
  west: "#3FB79F",
  francophone: "#158F7A",
  africa: "#FFFFFF",
} as const;

const opacities = { senegal: 1, west: 0.85, francophone: 0.45, africa: 0.16 } as const;

/** Carte stylisée de l'Afrique (points), mettant en avant le Sénégal et l'Afrique de l'Ouest. */
export function AfricaMap() {
  const dots = africaDots();
  return (
    <figure className="relative mx-auto w-full max-w-xl">
      <svg
        viewBox={`0 0 ${MAP_WIDTH} ${MAP_HEIGHT}`}
        className="h-auto w-full"
        role="img"
        aria-labelledby="africa-map-title"
      >
        <title id="africa-map-title">
          Carte stylisée de l&apos;Afrique mettant en évidence le Sénégal, l&apos;Afrique de l&apos;Ouest et l&apos;Afrique francophone
        </title>
        {dots.map((dot) => (
          <circle key={`${dot.x}-${dot.y}`} cx={dot.x} cy={dot.y} r={2.9} fill={fills[dot.region]} fillOpacity={opacities[dot.region]} />
        ))}
        {/* Rayonnement depuis Dakar */}
        <g>
          {[0, 1, 2].map((i) => (
            <circle
              key={i}
              cx={DAKAR.x}
              cy={DAKAR.y}
              r="10"
              fill="none"
              stroke="#D4AF7C"
              strokeWidth="1.5"
              className="animate-pulse-ring"
              style={{ transformOrigin: `${DAKAR.x}px ${DAKAR.y}px`, transformBox: "view-box", animationDelay: `${i * 0.85}s` }}
            />
          ))}
          <circle cx={DAKAR.x} cy={DAKAR.y} r="6" fill="#D4AF7C" />
          <circle cx={DAKAR.x} cy={DAKAR.y} r="2.5" fill="#041E19" />
        </g>
        {/* Liaisons vers l'international */}
        <g fill="none" stroke="#D4AF7C" strokeOpacity="0.5" strokeDasharray="3 6">
          <path d={`M${DAKAR.x} ${DAKAR.y} Q ${DAKAR.x + 40} ${DAKAR.y - 120} ${DAKAR.x + 150} ${-10}`} />
          <path d={`M${DAKAR.x} ${DAKAR.y} Q ${DAKAR.x - 30} ${DAKAR.y - 60} ${-10} ${DAKAR.y - 110}`} />
        </g>
      </svg>
    </figure>
  );
}
