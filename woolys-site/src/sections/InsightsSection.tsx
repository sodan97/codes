import { ArticleCard } from "@/components/ArticleCard";
import { ButtonLink } from "@/components/Button";
import { Reveal } from "@/components/Reveal";
import { SectionHeader } from "@/components/SectionHeader";
import { sortedArticles } from "@/data/insights";

export function InsightsSection() {
  const latest = sortedArticles().slice(0, 3);
  return (
    <section className="bg-sand-50 py-24 sm:py-32" aria-labelledby="insights-title">
      <div className="container-site">
        <div className="flex flex-col justify-between gap-8 lg:flex-row lg:items-end">
          <SectionHeader
            id="insights-title"
            eyebrow="Analyses & décryptages"
            title="Insights"
            intro="Nos analyses pour aider les dirigeants et les équipes à comprendre les enjeux de cybersécurité, de risques et de confiance numérique."
          />
          <Reveal>
            <ButtonLink href="/insights" variant="outlineDark" arrow>
              Tous les articles
            </ButtonLink>
          </Reveal>
        </div>
        <div className="mt-14 grid gap-6 md:grid-cols-2 lg:grid-cols-3">
          {latest.map((article, index) => (
            <Reveal key={article.slug} delay={index * 0.08} className="h-full">
              <ArticleCard article={article} index={index} />
            </Reveal>
          ))}
        </div>
      </div>
    </section>
  );
}
