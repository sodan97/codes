import { ButtonLink } from "@/components/Button";
import { Rings } from "@/components/Rings";

export default function NotFound() {
  return (
    <section className="relative isolate flex min-h-[80vh] items-center overflow-hidden bg-forest-900 pt-28 text-white">
      <div className="absolute inset-0 -z-10 bg-grid opacity-60" />
      <Rings className="absolute left-1/2 top-1/2 -z-10 h-[40rem] w-[40rem] -translate-x-1/2 -translate-y-1/2 opacity-40" />
      <div className="container-site text-center">
        <p className="eyebrow justify-center text-gold-400">Erreur 404</p>
        <h1 className="mx-auto mt-6 max-w-2xl text-4xl font-bold sm:text-5xl">Cette page est introuvable.</h1>
        <p className="mx-auto mt-6 max-w-xl text-lg text-forest-100/80">
          Le lien est peut-être obsolète ou la page a été déplacée.
        </p>
        <div className="mt-10 flex flex-col justify-center gap-4 sm:flex-row">
          <ButtonLink href="/" variant="gold" arrow>
            Retour à l&apos;accueil
          </ButtonLink>
          <ButtonLink href="/contact" variant="outline">
            Parler à un expert
          </ButtonLink>
        </div>
      </div>
    </section>
  );
}
