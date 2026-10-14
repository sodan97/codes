import type { Metadata } from "next";
import { PageHero } from "@/components/PageHero";
import { pageMetadata } from "@/lib/seo";
import { siteConfig } from "@/lib/site";

export const metadata: Metadata = {
  ...pageMetadata({ title: "Mentions légales", description: "Mentions légales du site WOOLYS.", path: "/mentions-legales" }),
  robots: { index: false, follow: true },
};

/* PLACEHOLDER : toutes les informations ci-dessous doivent être complétées avec les données officielles de l'entreprise. */
export default function LegalPage() {
  return (
    <>
      <PageHero breadcrumbs={[{ label: "Mentions légales", href: "/mentions-legales" }]} title="Mentions légales" />
      <section className="bg-sand-50 py-16 sm:py-20">
        <div className="container-site">
          <div className="prose-woolys mx-auto max-w-3xl">
            <h2>Éditeur du site</h2>
            <p>
              <strong>{siteConfig.legalName}</strong> — [Forme juridique, capital social, RCCM, NINEA, adresse du siège : à compléter].
            </p>
            <p>Directeur de la publication : [à compléter].</p>
            <p>Contact : {siteConfig.contactEmail ?? "[adresse e-mail à compléter]"}.</p>
            <h2>Hébergement</h2>
            <p>[Nom, adresse et coordonnées de l&apos;hébergeur : à compléter].</p>
            <h2>Propriété intellectuelle</h2>
            <p>
              L&apos;ensemble des contenus de ce site (textes, logo, éléments graphiques) est la propriété de {siteConfig.name}, sauf
              mention contraire. Toute reproduction sans autorisation préalable est interdite.
            </p>
            <h2>Marques et référentiels cités</h2>
            <p>
              Les noms de normes, méthodes et services cités (ISO, EBIOS RM, NIST, CIS, PCI-DSS, AWS…) appartiennent à leurs
              propriétaires respectifs. Leur mention n&apos;implique aucune certification, affiliation ou partenariat.
            </p>
          </div>
        </div>
      </section>
    </>
  );
}
