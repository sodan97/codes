import type { MetadataRoute } from "next";
import { expertises } from "@/data/expertises";
import { articles } from "@/data/insights";
import { sectors } from "@/data/sectors";
import { absoluteUrl } from "@/lib/site";

const staticRoutes: Array<{ path: string; priority: number }> = [
  { path: "/", priority: 1 },
  { path: "/expertises", priority: 0.9 },
  { path: "/advisory", priority: 0.9 },
  { path: "/academy", priority: 0.9 },
  { path: "/labs", priority: 0.8 },
  { path: "/secteurs", priority: 0.8 },
  { path: "/methodologies", priority: 0.8 },
  { path: "/insights", priority: 0.7 },
  { path: "/a-propos", priority: 0.7 },
  { path: "/contact", priority: 0.9 },
  { path: "/mentions-legales", priority: 0.2 },
  { path: "/confidentialite", priority: 0.2 },
];

export default function sitemap(): MetadataRoute.Sitemap {
  const now = new Date();
  return [
    ...staticRoutes.map(({ path, priority }) => ({
      url: absoluteUrl(path),
      lastModified: now,
      changeFrequency: "monthly" as const,
      priority,
    })),
    ...expertises.map((item) => ({
      url: absoluteUrl(`/expertises/${item.slug}`),
      lastModified: now,
      changeFrequency: "monthly" as const,
      priority: 0.8,
    })),
    ...sectors.map((item) => ({
      url: absoluteUrl(`/secteurs/${item.slug}`),
      lastModified: now,
      changeFrequency: "monthly" as const,
      priority: 0.7,
    })),
    ...articles.map((item) => ({
      url: absoluteUrl(`/insights/${item.slug}`),
      lastModified: new Date(`${item.date}T00:00:00Z`),
      changeFrequency: "yearly" as const,
      priority: 0.6,
    })),
  ];
}
