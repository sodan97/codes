import type { SVGProps } from "react";

/** Icônes linéaires maison (24×24, trait 1.5) pour un rendu sobre et cohérent. */
const paths = {
  governance: (
    <>
      <path d="M12 3 4 6.5v5.2c0 4.5 3.3 8.3 8 9.3 4.7-1 8-4.8 8-9.3V6.5L12 3Z" />
      <path d="M8.5 12h7M12 8.5v7" />
    </>
  ),
  radar: (
    <>
      <circle cx="12" cy="12" r="8.5" />
      <circle cx="12" cy="12" r="4.5" />
      <path d="M12 12 18 6" />
      <circle cx="15.5" cy="9.5" r="0.9" fill="currentColor" />
    </>
  ),
  compliance: (
    <>
      <rect x="5" y="3.5" width="14" height="17" rx="2" />
      <path d="M9 3.5h6v2.5H9z" />
      <path d="m8.5 13 2.2 2.2L15.5 10.5" />
    </>
  ),
  cloud: (
    <>
      <path d="M7 18.5h10a4 4 0 0 0 .6-7.96A5.5 5.5 0 0 0 7.1 9.6 4.5 4.5 0 0 0 7 18.5Z" />
      <path d="M12 12.5v3M10.5 14h3" />
    </>
  ),
  ai: (
    <>
      <rect x="6" y="6" width="12" height="12" rx="2.5" />
      <path d="M9.5 2.5v3.5M14.5 2.5v3.5M9.5 18v3.5M14.5 18v3.5M2.5 9.5H6M2.5 14.5H6M18 9.5h3.5M18 14.5h3.5" />
      <path d="m10 14.5 2-5 2 5M10.7 13h2.6" />
    </>
  ),
  people: (
    <>
      <circle cx="9" cy="8" r="3.2" />
      <path d="M3.5 19.5a5.5 5.5 0 0 1 11 0" />
      <circle cx="17" cy="9" r="2.4" />
      <path d="M15.5 14.2a4.5 4.5 0 0 1 5 4.3" />
    </>
  ),
  graduation: (
    <>
      <path d="M12 4 2.5 9 12 14l9.5-5L12 4Z" />
      <path d="M6.5 11.2v4.3c0 1.5 2.5 3 5.5 3s5.5-1.5 5.5-3v-4.3M21.5 9v5.5" />
    </>
  ),
  flask: (
    <>
      <path d="M9.5 3h5M10.5 3v6.2L5 18.2A1.8 1.8 0 0 0 6.6 21h10.8a1.8 1.8 0 0 0 1.6-2.8L13.5 9.2V3" />
      <path d="M7.5 15h9" />
    </>
  ),
  bank: (
    <>
      <path d="M3 9.5 12 4l9 5.5H3Z" />
      <path d="M5 9.5v8M9.5 9.5v8M14.5 9.5v8M19 9.5v8M3 20h18" />
    </>
  ),
  antenna: (
    <>
      <path d="M12 10v11M8.5 21h7" />
      <circle cx="12" cy="8.5" r="1.5" />
      <path d="M8 4.5a6 6 0 0 0 0 8M16 4.5a6 6 0 0 1 0 8M5.2 2a10 10 0 0 0 0 13M18.8 2a10 10 0 0 1 0 13" />
    </>
  ),
  umbrella: (
    <>
      <path d="M2.5 12a9.5 9.5 0 0 1 19 0h-19Z" />
      <path d="M12 12v6.5a2 2 0 0 1-4 0M12 2.5V2" />
    </>
  ),
  landmark: (
    <>
      <path d="M12 3v2M7 8.5a5 5 0 0 1 10 0" />
      <path d="M4 8.5h16M6 8.5v9M10 8.5v9M14 8.5v9M18 8.5v9M3 20.5h18M4 17.5h16" />
    </>
  ),
  bolt: <path d="M13 2.5 4.5 13.5H11l-1 8 8.5-11H12l1-8Z" />,
  factory: (
    <>
      <path d="M3 20.5V10l5 3V10l5 3V10l5 3V4h3v16.5H3Z" />
      <path d="M7 17h2M11 17h2M15 17h2" />
    </>
  ),
  chip: (
    <>
      <rect x="6.5" y="6.5" width="11" height="11" rx="1.5" />
      <rect x="9.5" y="9.5" width="5" height="5" rx="0.5" />
      <path d="M9.5 2.5v4M14.5 2.5v4M9.5 17.5v4M14.5 17.5v4M2.5 9.5h4M2.5 14.5h4M17.5 9.5h4M17.5 14.5h4" />
    </>
  ),
  infrastructure: (
    <>
      <circle cx="12" cy="5" r="2" />
      <circle cx="5" cy="18" r="2" />
      <circle cx="19" cy="18" r="2" />
      <circle cx="12" cy="13" r="2" />
      <path d="M12 7v4M10.4 14.3 6.5 16.8M13.6 14.3l3.9 2.5M7 18h10" />
    </>
  ),
  arrowRight: <path d="M4.5 12h15M13.5 6l6 6-6 6" />,
  arrowUpRight: <path d="M7 17 17 7M8.5 7H17v8.5" />,
  check: <path d="m5 12.5 4.5 4.5L19 7.5" />,
  plus: <path d="M12 5v14M5 12h14" />,
  menu: <path d="M4 7h16M4 12h16M4 17h10" />,
  close: <path d="M6 6l12 12M18 6 6 18" />,
  chevronDown: <path d="m6 9 6 6 6-6" />,
  calendar: (
    <>
      <rect x="3.5" y="5" width="17" height="15.5" rx="2" />
      <path d="M3.5 10h17M8 3v4M16 3v4" />
    </>
  ),
  mail: (
    <>
      <rect x="3" y="5" width="18" height="14" rx="2" />
      <path d="m3.5 6.5 8.5 7 8.5-7" />
    </>
  ),
  globe: (
    <>
      <circle cx="12" cy="12" r="9" />
      <path d="M3 12h18M12 3c2.5 2.6 3.8 5.6 3.8 9s-1.3 6.4-3.8 9c-2.5-2.6-3.8-5.6-3.8-9S9.5 5.6 12 3Z" />
    </>
  ),
  compass: (
    <>
      <circle cx="12" cy="12" r="9" />
      <path d="m15.5 8.5-2 5-5 2 2-5 5-2Z" />
    </>
  ),
  target: (
    <>
      <circle cx="12" cy="12" r="9" />
      <circle cx="12" cy="12" r="5" />
      <circle cx="12" cy="12" r="1.2" fill="currentColor" />
    </>
  ),
  layers: (
    <>
      <path d="m12 3 9 5-9 5-9-5 9-5Z" />
      <path d="m3 12.5 9 5 9-5M3 16.5l9 5 9-5" />
    </>
  ),
  fingerprint: (
    <>
      <path d="M6.5 5.5A8 8 0 0 1 20 11.5c0 2-.2 4-.7 5.8" />
      <path d="M4 9.5a8 8 0 0 0-.2 2c0 1.6.3 3.1.8 4.5" />
      <path d="M8.5 19.5c.9-2.4 1.2-5 1.2-8a2.3 2.3 0 0 1 4.6 0c0 3.2-.4 6.2-1.4 9" />
      <path d="M16.6 9a5 5 0 0 1 .4 2.5c0 2.5-.3 4.9-1 7.2M7 11.5a5 5 0 0 1 2.5-4.4" />
    </>
  ),
  linkedin: (
    <>
      <rect x="3" y="3" width="18" height="18" rx="2.5" />
      <path d="M8 10.5V17M8 7.2v.1M11.5 17v-6.5M11.5 13.2c0-1.6 1.2-2.7 2.6-2.7s2.4 1 2.4 2.7V17" />
    </>
  ),
  sparkle: <path d="M12 3c.6 4.6 3.4 7.4 9 9-5.6 1.6-8.4 4.4-9 9-.6-4.6-3.4-7.4-9-9 5.6-1.6 8.4-4.4 9-9Z" />,
  book: (
    <>
      <path d="M4 5.5A2.5 2.5 0 0 1 6.5 3H20v15H6.5A2.5 2.5 0 0 0 4 20.5v-15Z" />
      <path d="M4 20.5A2.5 2.5 0 0 1 6.5 18H20v3H6.5" />
    </>
  ),
  key: (
    <>
      <circle cx="8" cy="15" r="4.5" />
      <path d="m11.2 11.8 8.3-8.3M16.5 6.5l2.5 2.5M14 9l2 2" />
    </>
  ),
} as const;

export type IconName = keyof typeof paths;

type IconProps = SVGProps<SVGSVGElement> & { name: IconName; title?: string };

export function Icon({ name, title, className, ...props }: IconProps) {
  return (
    <svg
      viewBox="0 0 24 24"
      fill="none"
      stroke="currentColor"
      strokeWidth={1.5}
      strokeLinecap="round"
      strokeLinejoin="round"
      className={className ?? "h-6 w-6"}
      aria-hidden={title ? undefined : true}
      role={title ? "img" : undefined}
      focusable="false"
      {...props}
    >
      {title ? <title>{title}</title> : null}
      {paths[name]}
    </svg>
  );
}
