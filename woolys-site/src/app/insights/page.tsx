import type { Metadata } from "next";
import { ArticleCard } from "@/components/ArticleCard";
import { CTA } from "@/components/CTA";
import { PageHero } from "@/components/PageHero";
import { Reveal } from "@/components/Reveal";
import { sortedArticles } from "@/data/insights";
import { pageMetadata } from "@/lib/seo";

export const metadata: Metadata = pageMetadata({
  title: "Insights — Analyses cybersécurité, risques et IA",
  description:
    "Analyses et décryptages WOOLYS : gouvernance cyber, EBIOS RM, gestion des risques, conformité, Cloud et sécurité de l'intelligence artificielle.",
  path: "/insights",
});

export default function InsightsPage() {
  const articles = sortedArticles();
  const categories = Array.from(new Set(articles.map((article) => article.category)));
  return (
    <>
      <PageHero
        breadcrumbs={[{ label: "Insights", href: "/insights" }]}
        eyebrow="Insights"
        title="Analyses & décryptages."
        subtitle="Des contenus clairs pour aider les dirigeants et les équipes à comprendre les enjeux de cybersécurité, de risques et de confiance numérique."
      />
      <section className="bg-sand-50 py-20 sm:py-24" aria-label="Articles">
        <div className="container-site">
          <ul className="flex flex-wrap gap-2" aria-label="Thématiques">
            {categories.map((category) => (
              <li key={category} className="rounded-full border border-forest-800/15 bg-white px-4 py-1.5 text-sm text-forest-800">
                {category}
              </li>
            ))}
          </ul>
          <div className="mt-10 grid gap-6 md:grid-cols-2 lg:grid-cols-3">
            {articles.map((article, index) => (
              <Reveal key={article.slug} delay={(index % 3) * 0.08} className="h-full">
                <ArticleCard article={article} index={index} />
              </Reveal>
            ))}
          </div>
        </div>
      </section>
      <CTA />
    </>
  );
}
