import { NextResponse } from "next/server";
import { validateContact } from "@/lib/contact";

const WINDOW_MS = 10 * 60 * 1000;
const MAX_REQUESTS = 5;
// Limitation de débit simple, en mémoire (best effort, par instance).
const hits = new Map<string, number[]>();

function rateLimited(ip: string) {
  const now = Date.now();
  const recent = (hits.get(ip) ?? []).filter((time) => now - time < WINDOW_MS);
  recent.push(now);
  hits.set(ip, recent);
  if (hits.size > 5000) hits.clear();
  return recent.length > MAX_REQUESTS;
}

export async function POST(request: Request) {
  const ip = request.headers.get("x-forwarded-for")?.split(",")[0]?.trim() || "unknown";
  if (rateLimited(ip)) {
    return NextResponse.json({ ok: false, message: "Trop de demandes. Merci de réessayer dans quelques minutes." }, { status: 429 });
  }

  let body: Record<string, unknown>;
  try {
    body = await request.json();
  } catch {
    return NextResponse.json({ ok: false, message: "Requête invalide." }, { status: 400 });
  }

  // Champ piège anti-spam : rempli uniquement par les robots.
  if (typeof body.website === "string" && body.website.length > 0) {
    return NextResponse.json({ ok: true });
  }

  const { value, errors } = validateContact(body);
  if (Object.keys(errors).length > 0) {
    return NextResponse.json({ ok: false, errors, message: "Certains champs sont invalides." }, { status: 422 });
  }

  const webhook = process.env.CONTACT_WEBHOOK_URL?.trim();
  const submission = { ...value, source: "site-vitrine", submittedAt: new Date().toISOString() };

  if (!webhook) {
    if (process.env.NODE_ENV !== "production") {
      console.info("[contact] CONTACT_WEBHOOK_URL non configuré — demande reçue en développement :", submission);
      return NextResponse.json({ ok: true });
    }
    console.error("[contact] CONTACT_WEBHOOK_URL non configuré : la demande n'a pas pu être transmise.");
    return NextResponse.json(
      { ok: false, message: "Le service d'envoi est momentanément indisponible." },
      { status: 503 },
    );
  }

  try {
    const response = await fetch(webhook, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(submission),
      signal: AbortSignal.timeout(10_000),
    });
    if (!response.ok) throw new Error(`Webhook status ${response.status}`);
  } catch (error) {
    console.error("[contact] Échec de transmission au webhook :", error);
    return NextResponse.json({ ok: false, message: "L'envoi a échoué. Merci de réessayer." }, { status: 502 });
  }

  return NextResponse.json({ ok: true });
}
