import Link from "next/link";
import { Icon, type IconName } from "./Icon";

type IndustryCardProps = { name: string; summary: string; icon: IconName; href: string };

export function IndustryCard({ name, summary, icon, href }: IndustryCardProps) {
  return (
    <article className="group relative flex h-full flex-col justify-between overflow-hidden rounded-2xl border border-white/10 bg-white/[0.03] p-6 transition-all duration-500 hover:border-gold-400/50 hover:bg-white/[0.06]">
      <span className="flex h-11 w-11 items-center justify-center rounded-full border border-white/15 text-gold-400 transition-colors duration-500 group-hover:border-gold-400 group-hover:bg-gold-400 group-hover:text-forest-950">
        <Icon name={icon} className="h-5 w-5" />
      </span>
      <div className="mt-10">
        <h3 className="text-lg font-semibold text-white">
          <Link href={href} className="after:absolute after:inset-0 after:rounded-2xl focus-visible:outline-none">
            {name}
          </Link>
        </h3>
        <p className="mt-2 text-sm leading-relaxed text-forest-100/70">{summary}</p>
      </div>
      <Icon
        name="arrowUpRight"
        className="absolute right-6 top-6 h-5 w-5 text-white/30 transition-all duration-500 group-hover:-translate-y-0.5 group-hover:translate-x-0.5 group-hover:text-gold-400"
      />
    </article>
  );
}
