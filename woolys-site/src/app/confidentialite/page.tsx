import type { Metadata } from "next";
import { PageHero } from "@/components/PageHero";
import { pageMetadata } from "@/lib/seo";
import { siteConfig } from "@/lib/site";

export const metadata: Metadata = {
  ...pageMetadata({
    title: "Politique de confidentialité",
    description: "Comment WOOLYS traite les données personnelles collectées sur ce site.",
    path: "/confidentialite",
  }),
  robots: { index: false, follow: true },
};

/* PLACEHOLDER : texte de base à faire valider juridiquement (loi sénégalaise n° 2008-12, RGPD le cas échéant). */
export default function PrivacyPage() {
  return (
    <>
      <PageHero breadcrumbs={[{ label: "Confidentialité", href: "/confidentialite" }]} title="Politique de confidentialité" />
      <section className="bg-sand-50 py-16 sm:py-20">
        <div className="container-site">
          <div className="prose-woolys mx-auto max-w-3xl">
            <h2>Données collectées</h2>
            <p>
              Via le formulaire de contact, nous collectons : prénom, nom, entreprise, fonction, e-mail professionnel, téléphone
              (facultatif), pays, type de besoin, budget indicatif (facultatif) et message.
            </p>
            <h2>Finalité</h2>
            <p>Ces données sont utilisées uniquement pour répondre à votre demande et, le cas échéant, préparer une proposition.</p>
            <h2>Durée de conservation</h2>
            <p>[Durée de conservation à définir — par exemple 3 ans après le dernier contact].</p>
            <h2>Cadre légal</h2>
            <p>
              Les traitements sont réalisés conformément à la loi sénégalaise n° 2008-12 sur la protection des données à caractère
              personnel et, lorsqu&apos;il s&apos;applique, au Règlement général sur la protection des données (RGPD).
            </p>
            <h2>Vos droits</h2>
            <p>
              Vous disposez d&apos;un droit d&apos;accès, de rectification, d&apos;opposition et de suppression de vos données. Pour
              l&apos;exercer : {siteConfig.contactEmail ?? "[adresse e-mail à compléter]"}.
            </p>
            <h2>Cookies</h2>
            <p>Ce site n&apos;utilise pas de cookies publicitaires ni de traceurs tiers.</p>
          </div>
        </div>
      </section>
    </>
  );
}
