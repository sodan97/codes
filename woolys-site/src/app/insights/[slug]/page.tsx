import type { Metadata } from "next";
import { notFound } from "next/navigation";
import { ArticleCard } from "@/components/ArticleCard";
import { ArrowLink } from "@/components/Button";
import { CTA } from "@/components/CTA";
import { JsonLd } from "@/components/JsonLd";
import { PageHero } from "@/components/PageHero";
import type { ArticleBlock } from "@/data/types";
import { articles, getArticle, sortedArticles } from "@/data/insights";
import { formatDate } from "@/lib/format";
import { pageMetadata } from "@/lib/seo";
import { absoluteUrl, siteConfig } from "@/lib/site";

export const dynamicParams = false;

export function generateStaticParams() {
  return articles.map((article) => ({ slug: article.slug }));
}

export async function generateMetadata({ params }: PageProps<"/insights/[slug]">): Promise<Metadata> {
  const { slug } = await params;
  const article = getArticle(slug);
  if (!article) return {};
  return pageMetadata({
    title: article.title,
    description: article.excerpt,
    path: `/insights/${slug}`,
    keywords: article.keywords,
    type: "article",
    publishedTime: article.date,
  });
}

function Block({ block }: { block: ArticleBlock }) {
  switch (block.type) {
    case "h2":
      return <h2>{block.text}</h2>;
    case "h3":
      return <h3>{block.text}</h3>;
    case "ul":
      return (
        <ul>
          {block.items.map((item) => (
            <li key={item}>{item}</li>
          ))}
        </ul>
      );
    default:
      return <p>{block.text}</p>;
  }
}

export default async function ArticlePage({ params }: PageProps<"/insights/[slug]">) {
  const { slug } = await params;
  const article = getArticle(slug);
  if (!article) notFound();
  const related = sortedArticles().filter((other) => other.slug !== article.slug).slice(0, 2);

  return (
    <>
      <PageHero
        breadcrumbs={[
          { label: "Insights", href: "/insights" },
          { label: article.title, href: `/insights/${article.slug}` },
        ]}
        eyebrow={article.category}
        title={article.title}
        subtitle={
          <p className="text-base">
            <time dateTime={article.date}>{formatDate(article.date)}</time> · {article.readingMinutes} min de lecture · {article.author}
          </p>
        }
      />
      <article className="bg-sand-50 py-16 sm:py-24">
        <div className="container-site">
          <div className="mx-auto max-w-3xl">
            <p className="font-display text-xl font-medium leading-relaxed text-forest-800">{article.excerpt}</p>
            <div className="mt-10 h-px bg-gradient-to-r from-gold-400 to-transparent" />
            <div className="prose-woolys mt-10">
              {article.content.map((block, index) => (
                <Block key={index} block={block} />
              ))}
            </div>
            <div className="mt-14 rounded-2xl border border-forest-800/10 bg-white p-8">
              <p className="font-display text-lg font-semibold text-forest-800">Ce sujet concerne votre organisation ?</p>
              <p className="mt-2 text-muted">Échangez avec un expert WOOLYS pour en évaluer l&apos;impact sur vos activités.</p>
              <ArrowLink href="/contact" className="mt-4">
                Parler à un expert
              </ArrowLink>
            </div>
          </div>
        </div>
      </article>
      <section className="bg-sand-100 py-20" aria-labelledby="related-title">
        <div className="container-site">
          <h2 id="related-title" className="text-2xl font-bold text-forest-800">
            À lire aussi
          </h2>
          <div className="mt-10 grid gap-6 md:grid-cols-2">
            {related.map((other, index) => (
              <ArticleCard key={other.slug} article={other} index={index + 1} />
            ))}
          </div>
        </div>
      </section>
      <CTA />
      <JsonLd
        data={{
          "@context": "https://schema.org",
          "@type": "Article",
          headline: article.title,
          description: article.excerpt,
          datePublished: article.date,
          inLanguage: "fr-FR",
          author: { "@type": "Organization", name: siteConfig.name },
          publisher: { "@id": `${siteConfig.url}/#organization` },
          mainEntityOfPage: absoluteUrl(`/insights/${article.slug}`),
          keywords: article.keywords.join(", "),
        }}
      />
    </>
  );
}
