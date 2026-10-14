"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";
import { useEffect, useState } from "react";
import { t } from "@/lib/i18n";
import { cn } from "@/lib/cn";
import { Icon } from "./Icon";

/** CTA collant en bas d'écran sur mobile, affiché après le premier écran. */
export function MobileStickyCTA() {
  const pathname = usePathname();
  const [visible, setVisible] = useState(false);

  useEffect(() => {
    const onScroll = () => setVisible(window.scrollY > window.innerHeight * 0.7);
    onScroll();
    window.addEventListener("scroll", onScroll, { passive: true });
    return () => window.removeEventListener("scroll", onScroll);
  }, []);

  if (pathname.startsWith("/contact")) return null;

  return (
    <div
      className={cn(
        "fixed inset-x-0 bottom-0 z-40 border-t border-white/10 bg-forest-950/95 px-4 pb-[calc(0.75rem+env(safe-area-inset-bottom))] pt-3 backdrop-blur-xl transition-transform duration-500 sm:hidden",
        visible ? "translate-y-0" : "translate-y-full",
      )}
      aria-hidden={!visible}
    >
      <Link
        href="/contact"
        tabIndex={visible ? undefined : -1}
        className="group flex w-full items-center justify-center gap-2 rounded-full bg-gold-400 py-3.5 text-sm font-semibold text-forest-950 active:scale-[0.98]"
      >
        {t.cta.expert}
        <Icon name="arrowRight" className="h-4 w-4 transition-transform group-hover:translate-x-1" />
      </Link>
    </div>
  );
}
