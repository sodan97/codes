import type { ReactNode } from "react";
import { Reveal } from "./Reveal";
import { cn } from "@/lib/cn";

type SectionHeaderProps = {
  eyebrow?: string;
  title: ReactNode;
  intro?: ReactNode;
  align?: "left" | "center";
  tone?: "light" | "dark";
  as?: "h1" | "h2";
  className?: string;
  id?: string;
};

export function SectionHeader({
  eyebrow,
  title,
  intro,
  align = "left",
  tone = "light",
  as: Heading = "h2",
  className,
  id,
}: SectionHeaderProps) {
  const dark = tone === "dark";
  return (
    <Reveal className={cn("max-w-3xl", align === "center" && "mx-auto text-center", className)}>
      {eyebrow ? (
        <p className={cn("eyebrow mb-5", dark ? "text-gold-400" : "text-gold-700", align === "center" && "justify-center")}>
          {eyebrow}
        </p>
      ) : null}
      <Heading
        id={id}
        className={cn(
          "text-3xl font-bold leading-[1.15] sm:text-4xl lg:text-[2.75rem]",
          dark ? "text-white" : "text-forest-800",
        )}
      >
        {title}
      </Heading>
      {intro ? (
        <div className={cn("mt-5 text-base leading-relaxed sm:text-lg", dark ? "text-forest-100/80" : "text-muted")}>
          {intro}
        </div>
      ) : null}
    </Reveal>
  );
}
