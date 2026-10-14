"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";
import { useCallback, useEffect, useRef, useState } from "react";
import { mainNav } from "@/data/navigation";
import { t } from "@/lib/i18n";
import { siteConfig } from "@/lib/site";
import { cn } from "@/lib/cn";
import { ButtonLink } from "./Button";
import { Icon } from "./Icon";
import { Logo } from "./Logo";

const isActive = (pathname: string, href: string) => (href === "/" ? pathname === "/" : pathname === href || pathname.startsWith(`${href}/`));

export function Navbar() {
  const pathname = usePathname();
  const [scrolled, setScrolled] = useState(false);
  const [mobileOpen, setMobileOpen] = useState(false);
  const [dropdownOpen, setDropdownOpen] = useState(false);
  const toggleRef = useRef<HTMLButtonElement>(null);
  const panelRef = useRef<HTMLDivElement>(null);
  const dropdownRef = useRef<HTMLLIElement>(null);
  const [lastPathname, setLastPathname] = useState(pathname);

  // Fermeture des menus à chaque changement de page.
  if (lastPathname !== pathname) {
    setLastPathname(pathname);
    setMobileOpen(false);
    setDropdownOpen(false);
  }

  useEffect(() => {
    const onScroll = () => setScrolled(window.scrollY > 24);
    onScroll();
    window.addEventListener("scroll", onScroll, { passive: true });
    return () => window.removeEventListener("scroll", onScroll);
  }, []);

  const closeMobile = useCallback((restoreFocus = false) => {
    setMobileOpen(false);
    if (restoreFocus) toggleRef.current?.focus();
  }, []);

  // Menu mobile : verrouillage du scroll, focus initial, Échap et piège de focus.
  useEffect(() => {
    if (!mobileOpen) return;
    const previousOverflow = document.body.style.overflow;
    document.body.style.overflow = "hidden";
    const panel = panelRef.current;
    panel?.querySelector<HTMLElement>("a, button")?.focus();

    const onKeyDown = (event: KeyboardEvent) => {
      if (event.key === "Escape") {
        closeMobile(true);
        return;
      }
      if (event.key !== "Tab" || !panel) return;
      // Le bouton de fermeture (dans le header) fait partie du cycle de focus.
      const focusables = [
        ...(toggleRef.current ? [toggleRef.current] : []),
        ...Array.from(panel.querySelectorAll<HTMLElement>("a[href], button:not([disabled])")),
      ];
      if (focusables.length === 0) return;
      const first = focusables[0];
      const last = focusables[focusables.length - 1];
      if (event.shiftKey && document.activeElement === first) {
        event.preventDefault();
        last.focus();
      } else if (!event.shiftKey && document.activeElement === last) {
        event.preventDefault();
        first.focus();
      }
    };

    document.addEventListener("keydown", onKeyDown);
    return () => {
      document.body.style.overflow = previousOverflow;
      document.removeEventListener("keydown", onKeyDown);
    };
  }, [mobileOpen, closeMobile]);

  // Menu déroulant desktop : fermeture au clic extérieur et à Échap.
  useEffect(() => {
    if (!dropdownOpen) return;
    const onPointer = (event: PointerEvent) => {
      if (!dropdownRef.current?.contains(event.target as Node)) setDropdownOpen(false);
    };
    const onKey = (event: KeyboardEvent) => {
      if (event.key === "Escape") {
        setDropdownOpen(false);
        dropdownRef.current?.querySelector<HTMLButtonElement>("button")?.focus();
      }
    };
    document.addEventListener("pointerdown", onPointer);
    document.addEventListener("keydown", onKey);
    return () => {
      document.removeEventListener("pointerdown", onPointer);
      document.removeEventListener("keydown", onKey);
    };
  }, [dropdownOpen]);

  const solid = scrolled || mobileOpen;

  return (
    <>
      <header
        className={cn(
          "fixed inset-x-0 top-0 z-50 transition-all duration-500",
          solid ? "border-b border-white/10 bg-forest-950/90 shadow-[0_10px_30px_-20px_rgba(0,0,0,0.6)] backdrop-blur-xl" : "bg-transparent",
        )}
      >
        <div className={cn("container-site flex items-center justify-between gap-6 transition-all duration-500", solid ? "h-[4.5rem]" : "h-20 sm:h-24")}>
          <Link href="/" aria-label="WOOLYS — Accueil" className="shrink-0 rounded-lg" onClick={() => closeMobile()}>
            <Logo tone="onDark" />
          </Link>

          <nav aria-label={t.nav.main} className="hidden xl:block">
            <ul className="flex items-center gap-1">
              {mainNav.map((item) =>
                item.children ? (
                  <li
                    key={item.href}
                    ref={dropdownRef}
                    className="relative"
                    onMouseEnter={() => setDropdownOpen(true)}
                    onMouseLeave={() => setDropdownOpen(false)}
                  >
                    <button
                      type="button"
                      aria-expanded={dropdownOpen}
                      aria-controls="menu-expertises"
                      onClick={() => setDropdownOpen((open) => !open)}
                      className={cn(
                        "flex items-center gap-1 rounded-full px-3.5 py-2 text-sm font-medium transition-colors",
                        isActive(pathname, item.href) ? "text-gold-300" : "text-white/85 hover:text-white",
                      )}
                    >
                      {item.label}
                      <Icon name="chevronDown" className={cn("h-4 w-4 transition-transform duration-300", dropdownOpen && "rotate-180")} />
                    </button>
                    <div
                      id="menu-expertises"
                      className={cn(
                        "absolute left-1/2 top-full w-[40rem] -translate-x-1/2 pt-3 transition-all duration-300",
                        dropdownOpen ? "visible translate-y-0 opacity-100" : "invisible -translate-y-2 opacity-0",
                      )}
                    >
                      <div className="overflow-hidden rounded-2xl border border-white/10 bg-forest-950/95 p-3 shadow-2xl backdrop-blur-xl">
                        <ul className="grid grid-cols-2 gap-1">
                          {item.children.map((child) => (
                            <li key={child.href}>
                              <Link
                                href={child.href}
                                onClick={() => setDropdownOpen(false)}
                                className="block rounded-xl p-4 transition-colors hover:bg-white/5 focus-visible:bg-white/5"
                              >
                                <span className="block text-sm font-semibold text-white">{child.label}</span>
                                <span className="mt-1 block text-xs leading-relaxed text-forest-100/60">{child.description}</span>
                              </Link>
                            </li>
                          ))}
                        </ul>
                        <Link
                          href={item.href}
                          onClick={() => setDropdownOpen(false)}
                          className="mt-2 flex items-center justify-between rounded-xl border-t border-white/10 px-4 py-3 text-sm font-semibold text-gold-400 hover:text-gold-200"
                        >
                          Toutes nos expertises
                          <Icon name="arrowRight" className="h-4 w-4" />
                        </Link>
                      </div>
                    </div>
                  </li>
                ) : (
                  <li key={item.href}>
                    <Link
                      href={item.href}
                      aria-current={isActive(pathname, item.href) ? "page" : undefined}
                      className={cn(
                        "relative rounded-full px-3.5 py-2 text-sm font-medium transition-colors",
                        isActive(pathname, item.href) ? "text-gold-300" : "text-white/85 hover:text-white",
                      )}
                    >
                      {item.label}
                    </Link>
                  </li>
                ),
              )}
            </ul>
          </nav>

          <div className="flex items-center gap-3">
            <div className="hidden sm:block">
              <ButtonLink href="/contact" variant="primary" className="px-5 py-3">
                {t.cta.expert}
              </ButtonLink>
            </div>
            <button
              ref={toggleRef}
              type="button"
              className="flex h-11 w-11 items-center justify-center rounded-full border border-white/15 text-white transition-colors hover:border-gold-400 hover:text-gold-300 xl:hidden"
              aria-expanded={mobileOpen}
              aria-controls="menu-mobile"
              aria-label={mobileOpen ? t.nav.closeMenu : t.nav.openMenu}
              onClick={() => setMobileOpen((open) => !open)}
            >
              <Icon name={mobileOpen ? "close" : "menu"} className="h-5 w-5" />
            </button>
          </div>
        </div>

      </header>

      {/* Menu mobile / tablette — hors du <header> : son backdrop-filter casserait le positionnement fixed. */}
      <div
        id="menu-mobile"
        ref={panelRef}
        role="dialog"
        aria-modal="true"
        aria-label="Menu"
        hidden={!mobileOpen}
        className="fixed inset-x-0 bottom-0 top-[4.5rem] z-50 overflow-y-auto border-t border-white/10 bg-forest-950 xl:hidden"
      >
        <nav aria-label="Navigation mobile" className="container-site flex min-h-full flex-col py-8">
          <ul className="flex flex-col">
            {mainNav.map((item) => (
              <li key={item.href} className="border-b border-white/10">
                <Link
                  href={item.href}
                  onClick={() => closeMobile()}
                  aria-current={isActive(pathname, item.href) ? "page" : undefined}
                  className={cn(
                    "flex items-center justify-between py-4 font-display text-xl font-semibold",
                    isActive(pathname, item.href) ? "text-gold-300" : "text-white",
                  )}
                >
                  {item.label}
                  <Icon name="arrowRight" className="h-5 w-5 text-gold-400" />
                </Link>
                {item.children ? (
                  <ul className="grid gap-1 pb-4 sm:grid-cols-2">
                    {item.children.map((child) => (
                      <li key={child.href}>
                        <Link
                          href={child.href}
                          onClick={() => closeMobile()}
                          className="block rounded-lg py-2 pl-4 text-sm text-forest-100/75 hover:text-white"
                        >
                          {child.label}
                        </Link>
                      </li>
                    ))}
                  </ul>
                ) : null}
              </li>
            ))}
            <li className="border-b border-white/10">
              <Link
                href="/advisory"
                onClick={() => closeMobile()}
                className="flex items-center justify-between py-4 font-display text-xl font-semibold text-white"
              >
                Advisory
                <Icon name="arrowRight" className="h-5 w-5 text-gold-400" />
              </Link>
            </li>
          </ul>
          <div className="mt-auto flex flex-col gap-3 pt-10">
            <ButtonLink href="/contact" variant="gold" arrow className="w-full">
              {t.cta.expert}
            </ButtonLink>
            {siteConfig.lmsUrl ? (
              <ButtonLink href={siteConfig.lmsUrl} variant="outline" className="w-full">
                {t.cta.lms}
              </ButtonLink>
            ) : null}
          </div>
        </nav>
      </div>
    </>
  );
}
