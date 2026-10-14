import Link from "next/link";
import type { ReactNode } from "react";
import { Icon } from "./Icon";
import { cn } from "@/lib/cn";

export type ButtonVariant = "primary" | "gold" | "outline" | "outlineDark" | "ghost";

const base =
  "group relative inline-flex items-center justify-center gap-2 rounded-full px-6 py-3.5 text-sm font-semibold tracking-wide transition-all duration-300 ease-out focus-visible:outline-2 focus-visible:outline-offset-4 active:scale-[0.98] disabled:cursor-not-allowed disabled:opacity-60";

const variants: Record<ButtonVariant, string> = {
  primary:
    "bg-emerald-600 text-white shadow-[0_8px_24px_-12px_rgba(15,116,100,0.8)] hover:bg-emerald-700 hover:shadow-[0_12px_28px_-12px_rgba(15,116,100,0.9)]",
  gold: "bg-gold-400 text-forest-950 hover:bg-gold-300 shadow-[0_8px_24px_-12px_rgba(212,175,124,0.7)]",
  outline: "border border-white/25 text-white hover:border-gold-400 hover:text-gold-300 hover:bg-white/5",
  outlineDark: "border border-forest-800/25 text-forest-800 hover:border-forest-800 hover:bg-forest-800 hover:text-white",
  ghost: "px-0 py-2 text-emerald-600 hover:text-forest-800",
};

type CommonProps = {
  children: ReactNode;
  variant?: ButtonVariant;
  className?: string;
  arrow?: boolean;
};

const isExternal = (href: string) => /^https?:\/\//.test(href) || href.startsWith("mailto:");

export function ButtonLink({ href, children, variant = "primary", className, arrow = false }: CommonProps & { href: string }) {
  const content = (
    <>
      <span>{children}</span>
      {arrow ? (
        <Icon name="arrowRight" className="h-4 w-4 transition-transform duration-300 group-hover:translate-x-1" />
      ) : null}
    </>
  );

  if (isExternal(href)) {
    const newTab = href.startsWith("http");
    return (
      <a
        href={href}
        className={cn(base, variants[variant], className)}
        {...(newTab ? { target: "_blank", rel: "noopener noreferrer" } : {})}
      >
        {content}
        {newTab ? <span className="sr-only"> (s&apos;ouvre dans un nouvel onglet)</span> : null}
      </a>
    );
  }

  return (
    <Link href={href} className={cn(base, variants[variant], className)}>
      {content}
    </Link>
  );
}

export function Button({
  children,
  variant = "primary",
  className,
  arrow = false,
  ...props
}: CommonProps & React.ButtonHTMLAttributes<HTMLButtonElement>) {
  return (
    <button className={cn(base, variants[variant], className)} {...props}>
      <span>{children}</span>
      {arrow ? (
        <Icon name="arrowRight" className="h-4 w-4 transition-transform duration-300 group-hover:translate-x-1" />
      ) : null}
    </button>
  );
}

/** Lien texte avec flèche animée (« En savoir plus → »). */
export function ArrowLink({
  href,
  children,
  className,
  tone = "light",
}: {
  href: string;
  children: ReactNode;
  className?: string;
  tone?: "light" | "dark";
}) {
  return (
    <Link
      href={href}
      className={cn(
        "group inline-flex items-center gap-2 text-sm font-semibold transition-colors",
        tone === "light" ? "text-emerald-600 hover:text-forest-800" : "text-gold-400 hover:text-gold-200",
        className,
      )}
    >
      {children}
      <Icon name="arrowRight" className="h-4 w-4 transition-transform duration-300 group-hover:translate-x-1" />
    </Link>
  );
}
