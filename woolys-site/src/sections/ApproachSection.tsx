import { SectionHeader } from "@/components/SectionHeader";
import { Timeline } from "@/components/Timeline";
import { approachSteps } from "@/data/navigation";

export function ApproachSection({ tone = "light" }: { tone?: "light" | "dark" }) {
  const dark = tone === "dark";
  return (
    <section className={dark ? "bg-forest-950 py-24 text-white sm:py-32" : "bg-sand-100 py-24 sm:py-32"} aria-labelledby="approach-title">
      <div className="container-site">
        <SectionHeader
          id="approach-title"
          tone={tone}
          eyebrow="Notre approche"
          title={
            <>
              Une approche simple.
              <br />
              Une exigence élevée.
            </>
          }
        />
        <div className="mt-16">
          <Timeline steps={approachSteps} tone={tone} />
        </div>
      </div>
    </section>
  );
}
