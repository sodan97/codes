import Link from "next/link";
import { cn } from "@/lib/cn";

type MethodologyCardProps = {
  name: string;
  category: string;
  summary?: string;
  href?: string;
  tone?: "light" | "dark";
};

export function MethodologyCard({ name, category, summary, href, tone = "light" }: MethodologyCardProps) {
  const dark = tone === "dark";
  return (
    <article
      className={cn(
        "group relative flex h-full flex-col rounded-2xl border p-6 transition-all duration-500",
        dark
          ? "border-white/10 bg-forest-900/60 hover:border-gold-400/60"
          : "border-forest-800/10 bg-white hover:border-gold-400/70 hover:shadow-[0_20px_40px_-28px_rgba(11,61,51,0.4)]",
      )}
    >
      <p className={cn("text-[0.7rem] font-semibold uppercase tracking-[0.18em]", dark ? "text-gold-400" : "text-gold-700")}>
        {category}
      </p>
      <h3 className={cn("mt-3 font-display text-xl font-bold", dark ? "text-white" : "text-forest-800")}>
        {href ? (
          <Link href={href} className="after:absolute after:inset-0 after:rounded-2xl focus-visible:outline-none">
            {name}
          </Link>
        ) : (
          name
        )}
      </h3>
      {summary ? (
        <p className={cn("mt-3 text-sm leading-relaxed", dark ? "text-forest-100/70" : "text-muted")}>{summary}</p>
      ) : null}
    </article>
  );
}
