import type { FaqItem } from "@/data/types";
import { Icon } from "./Icon";
import { JsonLd } from "./JsonLd";

/** FAQ accessible basée sur <details>/<summary> (clavier et lecteurs d'écran natifs). */
export function FAQ({ items }: { items: FaqItem[] }) {
  return (
    <>
      <div className="divide-y divide-forest-800/10 border-y border-forest-800/10">
        {items.map((item) => (
          <details key={item.question} className="group py-2">
            <summary className="flex cursor-pointer list-none items-center justify-between gap-6 rounded-md py-4 text-left font-display text-base font-semibold text-forest-800 sm:text-lg [&::-webkit-details-marker]:hidden">
              {item.question}
              <span className="flex h-8 w-8 shrink-0 items-center justify-center rounded-full border border-forest-800/15 text-emerald-600 transition-transform duration-300 group-open:rotate-45">
                <Icon name="plus" className="h-4 w-4" />
              </span>
            </summary>
            <p className="pb-5 pr-12 leading-relaxed text-muted">{item.answer}</p>
          </details>
        ))}
      </div>
      <JsonLd
        data={{
          "@context": "https://schema.org",
          "@type": "FAQPage",
          mainEntity: items.map((item) => ({
            "@type": "Question",
            name: item.question,
            acceptedAnswer: { "@type": "Answer", text: item.answer },
          })),
        }}
      />
    </>
  );
}
