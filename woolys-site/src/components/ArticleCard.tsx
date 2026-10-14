import Link from "next/link";
import type { Article } from "@/data/types";
import { formatDate } from "@/lib/format";
import { Icon } from "./Icon";

/** Visuel abstrait généré (pas de photo stock) : variante selon l'index. */
function ArticleVisual({ seed }: { seed: number }) {
  const rings = [0, 1, 2, 3, 4];
  return (
    <svg viewBox="0 0 400 220" className="h-full w-full" aria-hidden="true" preserveAspectRatio="xMidYMid slice">
      <rect width="400" height="220" fill="#0B3D33" />
      <g fill="none" strokeWidth="1.2">
        {rings.map((ring) => (
          <circle
            key={ring}
            cx={seed % 2 === 0 ? 320 : 80}
            cy={seed === 1 ? 40 : 190}
            r={40 + ring * 34}
            stroke={ring % 2 === 0 ? "#D4AF7C" : "#158F7A"}
            strokeOpacity={0.55 - ring * 0.08}
            strokeDasharray={ring === 2 ? "4 8" : undefined}
          />
        ))}
      </g>
      <g fill="#D4AF7C">
        {Array.from({ length: 9 }, (_, i) => (
          <circle key={i} cx={40 + ((i * 47 + seed * 31) % 330)} cy={30 + ((i * 61 + seed * 17) % 160)} r={i % 3 === 0 ? 2.2 : 1.3} opacity={0.7} />
        ))}
      </g>
    </svg>
  );
}

export function ArticleCard({ article, index = 0 }: { article: Article; index?: number }) {
  return (
    <article className="group relative flex h-full flex-col overflow-hidden rounded-2xl border border-forest-800/10 bg-white transition-all duration-500 hover:-translate-y-1 hover:shadow-[0_24px_48px_-28px_rgba(11,61,51,0.45)]">
      <div className="relative aspect-[16/8] overflow-hidden">
        <div className="h-full w-full transition-transform duration-700 group-hover:scale-105">
          <ArticleVisual seed={index} />
        </div>
        <span className="absolute left-5 top-5 rounded-full bg-forest-950/70 px-3 py-1 text-xs font-semibold text-gold-300 backdrop-blur">
          {article.category}
        </span>
      </div>
      <div className="flex flex-1 flex-col p-7">
        <p className="text-xs font-medium uppercase tracking-wider text-muted">
          <time dateTime={article.date}>{formatDate(article.date)}</time> · {article.readingMinutes} min de lecture
        </p>
        <h3 className="mt-3 text-lg font-semibold leading-snug text-forest-800">
          <Link href={`/insights/${article.slug}`} className="after:absolute after:inset-0 focus-visible:outline-none">
            {article.title}
          </Link>
        </h3>
        <p className="mt-3 flex-1 text-[0.95rem] leading-relaxed text-muted">{article.excerpt}</p>
        <span className="mt-6 inline-flex items-center gap-2 text-sm font-semibold text-emerald-600" aria-hidden="true">
          Lire l&apos;article
          <Icon name="arrowRight" className="h-4 w-4 transition-transform duration-300 group-hover:translate-x-1" />
        </span>
      </div>
    </article>
  );
}
