#!/usr/bin/env node
/**
 * Capture real domain-specific ISM & Sub-ISM screenshots for `tanmay-starslab/multiphase-ism`:
 *   - docs/readme/hero-spectra.png  (1440x900 desktop 12-ion Voigt absorption spectra + fitted components + live 3D mini)
 *   - docs/readme/hero-3d.png       (1440x900 desktop 3D WebGL2 volume render of the 5-level multiphase ISM / mixing layer)
 *   - docs/readme/sim-picker.png    (1440x900 desktop with #sim-pop open showing all 6 ISM & sub-ISM simulation regimes)
 *   - docs/readme/compare-both.png  (1440x900 desktop Layout "Compare both" 2R x 2S/N across multiphase ISM ions)
 *   - docs/readme/mobile-3d.png     (390x844 mobile 3D view with peek bottom sheet)
 *   - docs/ism/img/og/explorer.jpg  (1200x630 OpenGraph preview card)
 */
import {createServer} from "node:http";
import {readFileSync, writeFileSync, existsSync, statSync, mkdirSync} from "node:fs";
import {join, extname, resolve, dirname} from "node:path";
import {launch, sleep, READY, ROOT} from "../tests/js/lib/browser.mjs";

const ISM_DIR = join(ROOT, "docs/ism");
const BUILD_DIR = join(ROOT, "build_5lvl");
const README_DIR = join(ROOT, "docs/readme");
const OG_FILE = join(ISM_DIR, "img/og/explorer.jpg");

const TYPES = {
  ".html": "text/html; charset=utf-8",
  ".js": "text/javascript",
  ".mjs": "text/javascript",
  ".json": "application/json",
  ".css": "text/css",
  ".webp": "image/webp",
  ".png": "image/png",
  ".jpg": "image/jpeg",
  ".svg": "image/svg+xml",
  ".gz": "application/octet-stream",
  ".bin": "application/octet-stream",
  ".wasm": "application/wasm",
  ".woff2": "font/woff2",
};

function fastManifestJson(filePath) {
  // For fast headless software-WebGL screenshot rendering on a 2-core node, expose L0+L1 and S0+S1
  // so 3D volume textures bind completely in ~3s with zero "loading L2, L3, L4..." overlay banner.
  const m = JSON.parse(readFileSync(filePath, "utf8"));
  const keepG = ["L0", "L1"].filter(L => m.gas?.levels?.includes(L));
  const keepS = ["S0", "S1"].filter(S => m.stars?.levels?.includes(S));
  const keepAll = new Set([...keepG, ...keepS]);
  m.levels = Object.fromEntries(Object.entries(m.levels || {}).filter(([k]) => keepAll.has(k)));
  if (m.gas) {
    m.gas.levels = keepG;
    for (const f of Object.values(m.gas.fields || {})) {
      f.cubes = Object.fromEntries(Object.entries(f.cubes || {}).filter(([k]) => keepG.includes(k)));
    }
  }
  if (m.stars) {
    m.stars.levels = keepS;
    m.stars.cubes = Object.fromEntries(Object.entries(m.stars.cubes || {}).filter(([k]) => keepS.includes(k)));
  }
  return JSON.stringify(m);
}

async function startLocalServer() {
  const srv = createServer((req, res) => {
    const u = new URL(req.url, "http://127.0.0.1");
    const p = decodeURIComponent(u.pathname);
    let f = null;
    if (p === "/" || p === "/ism/" || p === "/ism/index.html") {
      const rawHtml = readFileSync(join(ISM_DIR, "index.html"), "utf8");
      const localHtml = rawHtml.replaceAll(
        "https://huggingface.co/datasets/tanmay-starslab/ism-atlas-3d/resolve/main/",
        "/hf/"
      );
      res.writeHead(200, {
        "Content-Type": "text/html; charset=utf-8",
        "Cache-Control": "no-store",
        "Access-Control-Allow-Origin": "*",
      });
      res.end(localHtml);
      return;
    } else if (p.startsWith("/hf/")) {
      f = resolve(join(BUILD_DIR, p.slice(4)));
    } else if (p.startsWith("/ism/data/viz/")) {
      const rel = p.slice("/ism/data/viz/".length);
      const cand = resolve(join(BUILD_DIR, "RTML_50pc/viz", rel));
      f = existsSync(cand) ? cand : resolve(join(ISM_DIR, "data/viz", rel));
    } else if (p.startsWith("/ism/data/gal/")) {
      const rel = p.slice("/ism/data/gal/".length);
      const cand = resolve(join(BUILD_DIR, rel));
      f = existsSync(cand) ? cand : resolve(join(ISM_DIR, "data/gal", rel));
    } else if (p.startsWith("/ism/")) {
      f = resolve(join(ISM_DIR, p.slice(5)));
    }
    if (!f || !existsSync(f) || statSync(f).isDirectory()) {
      res.writeHead(404);
      res.end("not found: " + p);
      return;
    }
    if (f.endsWith("manifest_v2.json")) {
      res.writeHead(200, {
        "Content-Type": "application/json",
        "Cache-Control": "no-store",
        "Access-Control-Allow-Origin": "*",
      });
      res.end(fastManifestJson(f));
      return;
    }
    res.writeHead(200, {
      "Content-Type": TYPES[extname(f)] || "application/octet-stream",
      "Cache-Control": "no-store",
      "Access-Control-Allow-Origin": "*",
    });
    res.end(readFileSync(f));
  });
  await new Promise((r, j) => {
    srv.once("error", j);
    srv.listen(0, "127.0.0.1", r);
  });
  const port = srv.address().port;
  return {
    port,
    base: `http://127.0.0.1:${port}/ism/?bg=off&tier=high`,
    close: () => new Promise(r => { srv.closeAllConnections?.(); srv.close(r); }),
  };
}

async function waitFor3dReady(P, timeoutMs = 45000) {
  await P.waitFor(
    `window.tngViewer3d && window.tngViewer3d.species && !JSON.stringify(window.tngViewer3d.levels()).includes("loading") && !document.querySelector(".v3d-msg")?.textContent`,
    timeoutMs,
    200
  );
  await P.eval(`window.tngViewer3d.render && window.tngViewer3d.render()`);
  await sleep(500);
}

async function waitForThumbReady(P, timeoutMs = 30000) {
  await P.waitFor(
    `document.querySelector("#v3d-mini")?.classList.contains("live") && !document.querySelector("#v3d-mini .v3d-msg")?.textContent`,
    timeoutMs,
    200
  );
  await sleep(400);
}

async function waitForFitsReady(P, timeoutMs = 25000) {
  await P.waitFor(
    `document.querySelector("#fitstatus")?.dataset?.state === "done" || (window.tngApp && window.tngApp._debug.fitCache() >= 4)`,
    timeoutMs,
    150
  );
  await sleep(500);
}

async function captureJpg(P, file, quality = 92) {
  mkdirSync(dirname(file), {recursive: true});
  const r = await P.send("Page.captureScreenshot", {format: "jpeg", quality});
  writeFileSync(file, Buffer.from(r.result.data, "base64"));
  return file;
}

async function main() {
  const srv = await startLocalServer();
  const edgeBin = process.env.TNG_EDGE || "/scratch/tsingh65/cgm_work/bin/chrome/headless-rm";
  const P = await launch({
    edge: edgeBin,
    width: 1440,
    height: 900,
    hover: true,
  });

  try {
    await P.scheme(true);
    await P.viewport(1440, 900, false);

    // 1. Hero Spectra on RTML_50pc with Fitted components enabled and live 3D mini volume
    console.log("[shots] 1/6 Loading RTML_50pc Spectra view at", srv.base);
    await P.goto(srv.base);
    await P.waitFor(READY, 30000);
    await P.eval(`(() => {
      const h = document.querySelector("#story-row .sh-x");
      if (h && !h.hidden) h.click();
      document.activeElement?.blur?.();
      // Pick the richest multi-ion sightline and enable fitted Voigt components
      const best = window.tngApp.entries().find(e => (e.showcase_ions || []).length >= 3) || window.tngApp.entries()[0];
      if (best) window.tngApp.select(best.id);
      window.tngApp.setInstrument({fit: true, ions: ["HI", "MgII", "SiIII", "CIV", "OVI"], vhw: 250});
    })()`);
    await waitForFitsReady(P);
    await waitForThumbReady(P);
    await P.eval(`document.activeElement?.blur?.()`);
    await P.hover(5, 5);
    await sleep(300);
    await P.shot(join(README_DIR, "hero-spectra.png"));
    console.log("[shots] Saved docs/readme/hero-spectra.png");

    // 2. Open Simulation regime picker (#top-sim -> #sim-pop) showing all 6 flagship ISM & sub-ISM regimes
    await P.click("#top-sim");
    await P.waitFor(`!document.querySelector("#sim-pop")?.hidden`, 5000);
    await P.hover(5, 5);
    await sleep(350);
    await P.shot(join(README_DIR, "sim-picker.png"));
    console.log("[shots] Saved docs/readme/sim-picker.png");

    // 3. Compare both (2R x 2S/N) on CloudCrush_200pc with rich cometary wake sightline & fitted components
    console.log("[shots] 3/6 Loading CloudCrush_200pc Compare both view");
    await P.goto(srv.base + "#sim=cloudcrush&g=CloudCrush_200pc");
    await P.waitFor(READY, 30000);
    await P.eval(`(() => {
      const h = document.querySelector("#story-row .sh-x");
      if (h && !h.hidden) h.click();
      document.activeElement?.blur?.();
      const best = window.tngApp.entries().find(e => (e.showcase_ions || []).length >= 3) || window.tngApp.entries()[0];
      if (best) window.tngApp.select(best.id);
      window.tngApp.setInstrument({
        layout: "4",
        fit: true,
        ions: ["HI", "MgII", "SiIII", "CIV", "OVI"],
        R: {uv: [45000, 18000]},
        snr: [30, 10],
        vhw: 250
      });
    })()`);
    await waitForFitsReady(P);
    await waitForThumbReady(P);
    await P.eval(`document.activeElement?.blur?.()`);
    await P.hover(5, 5);
    await sleep(350);
    await P.shot(join(README_DIR, "compare-both.png"));
    console.log("[shots] Saved docs/readme/compare-both.png");

    // 4. Hero 3D view on RTML_50pc (Volume mode, Gas phases THOR transfer function + starlight)
    console.log("[shots] 4/6 Loading RTML_50pc 3D volume view");
    await P.goto(srv.base);
    await P.waitFor(READY, 30000);
    const res3d = await P.eval(`window.tngApp.show3d({species: "T", mode: "volume", tf: "thor"})`, 90000);
    console.log("[shots] show3d RTML_50pc result:", JSON.stringify(res3d));
    await waitFor3dReady(P, 45000);
    await P.eval(`document.activeElement?.blur?.()`);
    await P.hover(5, 5);
    await sleep(400);
    await P.shot(join(README_DIR, "hero-3d.png"));
    console.log("[shots] Saved docs/readme/hero-3d.png");

    // 5. OpenGraph 1200x630 preview card (docs/ism/img/og/explorer.jpg)
    console.log("[shots] 5/6 Capturing 1200x630 OpenGraph card");
    await P.viewport(1200, 630, false);
    await sleep(500);
    await P.eval(`window.tngViewer3d && window.tngViewer3d.render && window.tngViewer3d.render()`);
    await sleep(300);
    await captureJpg(P, OG_FILE, 92);
    console.log("[shots] Saved docs/ism/img/og/explorer.jpg");

    // 6. Mobile 3D view (390x844) with peek bottom sheet on GMC_STARFORGE_20pc
    console.log("[shots] 6/6 Capturing Mobile 3D view on GMC_STARFORGE_20pc");
    await P.viewport(390, 844, true);
    await P.goto(srv.base + "#sim=starforge&g=GMC_STARFORGE_20pc");
    await P.waitFor(READY, 30000);
    const resMob = await P.eval(`window.tngApp.show3d({species: "nH", mode: "volume", tf: "thor"})`, 90000);
    console.log("[shots] show3d mobile GMC_STARFORGE_20pc result:", JSON.stringify(resMob));
    await waitFor3dReady(P, 45000);
    await sleep(400);
    await P.shot(join(README_DIR, "mobile-3d.png"));
    console.log("[shots] Saved docs/readme/mobile-3d.png");
  } finally {
    await P.close();
    await srv.close();
  }
}

main().catch(err => {
  console.error(err);
  process.exit(1);
});
