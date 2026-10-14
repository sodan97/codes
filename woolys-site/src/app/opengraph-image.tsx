import { readFile } from "node:fs/promises";
import path from "node:path";
import { ImageResponse } from "next/og";

export const alt = "WOOLYS — Cybersecurity & Trust. La confiance au cœur du numérique.";
export const size = { width: 1200, height: 630 };
export const contentType = "image/png";

export default async function OpengraphImage() {
  const svg = await readFile(path.join(process.cwd(), "public/brand/woolys-emblem.svg"), "utf8");
  const emblem = `data:image/svg+xml;base64,${Buffer.from(svg).toString("base64")}`;

  return new ImageResponse(
    (
      <div
        style={{
          width: "100%",
          height: "100%",
          display: "flex",
          flexDirection: "column",
          alignItems: "center",
          justifyContent: "center",
          background: "radial-gradient(circle at 75% 20%, #0F4A3E 0%, #072C25 45%, #041E19 100%)",
          color: "#FFFFFF",
          fontFamily: "sans-serif",
        }}
      >
        <img src={emblem} width={210} height={210} alt="" />
        <div style={{ marginTop: 28, fontSize: 92, fontWeight: 700, letterSpacing: 14 }}>WOOLYS</div>
        <div style={{ marginTop: 8, fontSize: 26, letterSpacing: 10, color: "#D4AF7C" }}>CYBERSECURITY &amp; TRUST</div>
        <div style={{ marginTop: 22, width: 160, height: 2, background: "#D4AF7C" }} />
        <div style={{ marginTop: 22, fontSize: 32, color: "#D5E6E1" }}>La confiance au cœur du numérique.</div>
      </div>
    ),
    size,
  );
}
