import Link from "next/link";
import { Icon, type IconName } from "./Icon";
import { cn } from "@/lib/cn";

type ServiceCardProps = {
  number?: string;
  title: string;
  description: string;
  icon?: IconName;
  href?: string;
  linkLabel?: string;
  className?: string;
};

/** Carte d'expertise / de service. Toute la carte est cliquable lorsqu'un lien est fourni. */
export function ServiceCard({ number, title, description, icon, href, linkLabel = "En savoir plus", className }: ServiceCardProps) {
  return (
    <article
      className={cn(
        "group relative flex h-full flex-col rounded-2xl border border-forest-800/10 bg-white p-7 transition-all duration-500 hover:-translate-y-1 hover:border-gold-400/60 hover:shadow-[0_24px_48px_-28px_rgba(11,61,51,0.45)]",
        className,
      )}
    >
      <div className="mb-8 flex items-start justify-between">
        {icon ? (
          <span className="flex h-12 w-12 items-center justify-center rounded-xl bg-forest-50 text-emerald-600 transition-colors duration-500 group-hover:bg-forest-800 group-hover:text-gold-400">
            <Icon name={icon} className="h-6 w-6" />
          </span>
        ) : (
          <span />
        )}
        {number ? <span className="font-display text-sm font-semibold text-gold-700">{number}</span> : null}
      </div>
      <h3 className="text-lg font-semibold leading-snug text-forest-800">
        {href ? (
          <Link href={href} className="after:absolute after:inset-0 after:rounded-2xl focus-visible:outline-none">
            {title}
          </Link>
        ) : (
          title
        )}
      </h3>
      <p className="mt-3 flex-1 text-[0.95rem] leading-relaxed text-muted">{description}</p>
      {href ? (
        <span className="mt-6 inline-flex items-center gap-2 text-sm font-semibold text-emerald-600 transition-colors group-hover:text-forest-800" aria-hidden="true">
          {linkLabel}
          <Icon name="arrowRight" className="h-4 w-4 transition-transform duration-300 group-hover:translate-x-1" />
        </span>
      ) : null}
      <span className="pointer-events-none absolute inset-x-7 bottom-0 h-px origin-left scale-x-0 bg-gold-400 transition-transform duration-500 group-hover:scale-x-100" />
    </article>
  );
}
