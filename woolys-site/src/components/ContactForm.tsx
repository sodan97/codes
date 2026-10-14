"use client";

import { useSearchParams } from "next/navigation";
import { useId, useRef, useState, type FormEvent, type ReactNode } from "react";
import { budgetRanges, needTypes } from "@/data/navigation";
import { countries, needFromParam, validateContact, type ContactErrors, type ContactPayload } from "@/lib/contact";
import { siteConfig } from "@/lib/site";
import { cn } from "@/lib/cn";
import { Button } from "./Button";
import { Icon } from "./Icon";

type Status = "idle" | "submitting" | "success" | "error";

const fieldBase =
  "mt-2 block w-full rounded-xl border bg-white px-4 py-3 text-[0.95rem] text-ink placeholder:text-muted/60 transition-colors focus:border-emerald-600 focus:outline-none focus:ring-2 focus:ring-emerald-600/20";

function Field({
  id,
  label,
  required,
  optional,
  error,
  hint,
  children,
  className,
}: {
  id: string;
  label: string;
  required?: boolean;
  optional?: boolean;
  error?: string;
  hint?: string;
  children: ReactNode;
  className?: string;
}) {
  return (
    <div className={className}>
      <label htmlFor={id} className="text-sm font-semibold text-forest-800">
        {label}
        {required ? (
          <span className="text-emerald-600" aria-hidden="true">
            {" "}
            *
          </span>
        ) : null}
        {optional ? <span className="font-normal text-muted"> (optionnel)</span> : null}
      </label>
      {children}
      {hint && !error ? (
        <p id={`${id}-hint`} className="mt-1.5 text-xs text-muted">
          {hint}
        </p>
      ) : null}
      {error ? (
        <p id={`${id}-error`} className="mt-1.5 text-xs font-medium text-red-700">
          {error}
        </p>
      ) : null}
    </div>
  );
}

/** Formulaire de contact. `defaultNeed` permet de présélectionner le type de besoin. */
export function ContactFormFields({ defaultNeed = "" }: { defaultNeed?: string }) {
  const uid = useId();
  const id = (name: string) => `${uid}-${name}`;
  const formRef = useRef<HTMLFormElement>(null);
  const [status, setStatus] = useState<Status>("idle");
  const [errors, setErrors] = useState<ContactErrors>({});
  const [serverMessage, setServerMessage] = useState("");

  const describedBy = (name: keyof ContactPayload, hint = false) =>
    errors[name] ? `${id(name)}-error` : hint ? `${id(name)}-hint` : undefined;

  const fieldClass = (name: keyof ContactPayload) =>
    cn(fieldBase, errors[name] ? "border-red-600/70" : "border-forest-800/15");

  async function onSubmit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    const form = event.currentTarget;
    const data = Object.fromEntries(new FormData(form).entries());
    const { errors: found } = validateContact({ ...data, consent: data.consent === "on" });
    setErrors(found);

    const firstInvalid = Object.keys(found)[0];
    if (firstInvalid) {
      form.querySelector<HTMLElement>(`[name="${firstInvalid}"]`)?.focus();
      return;
    }

    setStatus("submitting");
    setServerMessage("");
    try {
      const response = await fetch("/api/contact", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ ...data, consent: true }),
      });
      const result = (await response.json().catch(() => ({}))) as { ok?: boolean; message?: string; errors?: ContactErrors };
      if (!response.ok || !result.ok) {
        if (result.errors) setErrors(result.errors);
        setServerMessage(result.message ?? "L'envoi a échoué. Merci de réessayer.");
        setStatus("error");
        return;
      }
      setStatus("success");
      formRef.current?.reset();
    } catch {
      setServerMessage("Connexion impossible. Vérifiez votre réseau et réessayez.");
      setStatus("error");
    }
  }

  if (status === "success") {
    return (
      <div role="status" className="rounded-2xl border border-emerald-600/30 bg-white p-10 text-center">
        <span className="mx-auto flex h-14 w-14 items-center justify-center rounded-full bg-emerald-600 text-white">
          <Icon name="check" className="h-7 w-7" />
        </span>
        <h3 className="mt-6 text-2xl font-bold text-forest-800">Merci, votre demande a bien été envoyée.</h3>
        <p className="mx-auto mt-3 max-w-md text-muted">
          Un expert WOOLYS reviendra vers vous rapidement pour échanger sur votre projet.
        </p>
        <Button variant="outlineDark" className="mt-8" onClick={() => setStatus("idle")}>
          Envoyer une autre demande
        </Button>
      </div>
    );
  }

  return (
    <form ref={formRef} onSubmit={onSubmit} noValidate className="grid gap-6 sm:grid-cols-2" aria-describedby={`${uid}-required`}>
      <p id={`${uid}-required`} className="text-sm text-muted sm:col-span-2">
        Les champs marqués d&apos;un <span className="text-emerald-600">*</span> sont obligatoires.
      </p>

      <Field id={id("firstName")} label="Prénom" required error={errors.firstName}>
        <input id={id("firstName")} name="firstName" type="text" autoComplete="given-name" required aria-invalid={!!errors.firstName} aria-describedby={describedBy("firstName")} className={fieldClass("firstName")} />
      </Field>
      <Field id={id("lastName")} label="Nom" required error={errors.lastName}>
        <input id={id("lastName")} name="lastName" type="text" autoComplete="family-name" required aria-invalid={!!errors.lastName} aria-describedby={describedBy("lastName")} className={fieldClass("lastName")} />
      </Field>
      <Field id={id("company")} label="Entreprise" required error={errors.company}>
        <input id={id("company")} name="company" type="text" autoComplete="organization" required aria-invalid={!!errors.company} aria-describedby={describedBy("company")} className={fieldClass("company")} />
      </Field>
      <Field id={id("jobTitle")} label="Fonction" required error={errors.jobTitle}>
        <input id={id("jobTitle")} name="jobTitle" type="text" autoComplete="organization-title" required placeholder="DSI, RSSI, Directeur des risques…" aria-invalid={!!errors.jobTitle} aria-describedby={describedBy("jobTitle")} className={fieldClass("jobTitle")} />
      </Field>
      <Field id={id("email")} label="Email professionnel" required error={errors.email}>
        <input id={id("email")} name="email" type="email" autoComplete="email" inputMode="email" required aria-invalid={!!errors.email} aria-describedby={describedBy("email")} className={fieldClass("email")} />
      </Field>
      <Field id={id("phone")} label="Téléphone" optional error={errors.phone}>
        <input id={id("phone")} name="phone" type="tel" autoComplete="tel" inputMode="tel" placeholder="+221 …" aria-invalid={!!errors.phone} aria-describedby={describedBy("phone")} className={fieldClass("phone")} />
      </Field>
      <Field id={id("country")} label="Pays" required error={errors.country}>
        <select id={id("country")} name="country" required defaultValue="" aria-invalid={!!errors.country} aria-describedby={describedBy("country")} className={fieldClass("country")}>
          <option value="" disabled>
            Sélectionnez un pays
          </option>
          {countries.map((country) => (
            <option key={country}>{country}</option>
          ))}
        </select>
      </Field>
      <Field id={id("need")} label="Type de besoin" required error={errors.need}>
        <select id={id("need")} name="need" required defaultValue={defaultNeed} aria-invalid={!!errors.need} aria-describedby={describedBy("need")} className={fieldClass("need")}>
          <option value="" disabled>
            Sélectionnez votre besoin
          </option>
          {needTypes.map((need) => (
            <option key={need}>{need}</option>
          ))}
        </select>
      </Field>
      <Field id={id("budget")} label="Budget indicatif" optional error={errors.budget} className="sm:col-span-2">
        <select id={id("budget")} name="budget" defaultValue="" aria-invalid={!!errors.budget} aria-describedby={describedBy("budget")} className={fieldClass("budget")}>
          <option value="">Je ne sais pas encore</option>
          {budgetRanges.map((range) => (
            <option key={range}>{range}</option>
          ))}
        </select>
      </Field>
      <Field id={id("message")} label="Message" required error={errors.message} hint="Contexte, objectifs, échéances : quelques lignes suffisent." className="sm:col-span-2">
        <textarea id={id("message")} name="message" rows={6} required aria-invalid={!!errors.message} aria-describedby={describedBy("message", true)} className={cn(fieldClass("message"), "resize-y")} />
      </Field>

      {/* Champ piège anti-spam, invisible pour les humains */}
      <div aria-hidden="true" className="absolute -left-[9999px] h-px w-px overflow-hidden">
        <label htmlFor={id("website")}>Ne pas remplir</label>
        <input id={id("website")} name="website" type="text" tabIndex={-1} autoComplete="off" />
      </div>

      <div className="sm:col-span-2">
        <div className="flex items-start gap-3">
          <input
            id={id("consent")}
            name="consent"
            type="checkbox"
            required
            aria-invalid={!!errors.consent}
            aria-describedby={describedBy("consent")}
            className="mt-1 h-5 w-5 shrink-0 rounded border-forest-800/30 accent-emerald-600"
          />
          <label htmlFor={id("consent")} className="text-sm leading-relaxed text-muted">
            J&apos;accepte que WOOLYS utilise ces informations pour répondre à ma demande, conformément à la{" "}
            <a href="/confidentialite" className="font-medium text-emerald-600 underline underline-offset-2 hover:text-forest-800">
              politique de confidentialité
            </a>
            .<span className="text-emerald-600"> *</span>
          </label>
        </div>
        {errors.consent ? (
          <p id={`${id("consent")}-error`} className="ml-8 mt-1.5 text-xs font-medium text-red-700">
            {errors.consent}
          </p>
        ) : null}
      </div>

      <div className="flex flex-col gap-4 sm:col-span-2 sm:flex-row sm:items-center">
        <Button type="submit" arrow disabled={status === "submitting"} className="w-full sm:w-auto">
          {status === "submitting" ? "Envoi en cours…" : "Envoyer ma demande"}
        </Button>
        <p className="text-xs text-muted">Vos informations sont traitées de manière confidentielle.</p>
      </div>

      <div aria-live="assertive" className="sm:col-span-2">
        {status === "error" ? (
          <p className="rounded-xl border border-red-600/30 bg-red-50 p-4 text-sm text-red-800">
            {serverMessage}
            {siteConfig.contactEmail ? (
              <>
                {" "}Vous pouvez aussi nous écrire directement à{" "}
                <a className="font-semibold underline" href={`mailto:${siteConfig.contactEmail}`}>
                  {siteConfig.contactEmail}
                </a>
                .
              </>
            ) : null}
          </p>
        ) : null}
      </div>
    </form>
  );
}

/** Variante qui lit `?besoin=` dans l'URL pour présélectionner le type de besoin. */
export function ContactForm() {
  const params = useSearchParams();
  return <ContactFormFields defaultNeed={needFromParam(params.get("besoin"))} />;
}
