import assert from "node:assert/strict";
import {readFileSync} from "node:fs";
import {join, dirname} from "node:path";
import {fileURLToPath} from "node:url";
import {
  CLOUDY_LINES,
  DEFAULT_EMISSION_LINE,
  EPS_REF,
  LOG10_4PI_ARCSEC2_PER_SR,
  CM_PER_KPC,
  getEmissionLine,
  trilinearLogEps,
  integrateLogSB
} from "../../docs/ism/js/cloudy_emissivity.js";

const ROOT = join(dirname(fileURLToPath(import.meta.url)), "..", "..");
const meta = JSON.parse(readFileSync(join(ROOT, "docs/ism/data/cloudy_emissivity_meta.json"), "utf8"));
const raw = new Uint8Array(readFileSync(join(ROOT, "docs/ism/data/cloudy_emissivity_z0.bin")));

assert.equal(CLOUDY_LINES.length, 12, "all 12 CLOUDY lines present");
assert.equal(DEFAULT_EMISSION_LINE, "Ha_6563", "default ISM emission line is Halpha 6563 A");
assert.equal(getEmissionLine("Ha_6563").label, "H\u03b1 6563 \u00c5");
assert.equal(getEmissionLine("OIII_5007").label, "[O III] 5007 \u00c5");
assert.equal(raw.length, 12 * 56 * 6 * 61, "binary LUT exact byte length");

// Physical sanity checks at nH=1 cm^-3, Z=Zsun, T=10^4 K
const haWarm = trilinearLogEps(raw, meta, 0, 0.0, 0.0, 4.0);
assert.ok(haWarm > -30 && haWarm < -22, `Halpha at 1e4 K should be bright (got ${haWarm})`);

// O III 5007 peaks in warm-hot gas (~3e4 - 1e5 K)
const oiiiWarm = trilinearLogEps(raw, meta, 3, 0.0, 0.0, 4.8);
const oiiiCold = trilinearLogEps(raw, meta, 3, 0.0, 0.0, 3.2);
assert.ok(oiiiWarm > oiiiCold + 5, `[O III] 5007 should be much brighter at 10^4.8 K (${oiiiWarm}) than 10^3.2 K (${oiiiCold})`);

// O VII 21.6 A X-ray line peaks near 10^6.3 K
const oviiHot = trilinearLogEps(raw, meta, 9, -2.0, 0.0, 6.3);
const oviiCold = trilinearLogEps(raw, meta, 9, -2.0, 0.0, 4.0);
assert.ok(oviiHot > oviiCold + 10, `O VII 21.6 should peak in hot X-ray gas (${oviiHot} vs ${oviiCold})`);

// Exact SB formula check: uniform slab with logEps = -24, L = 10 pc = 0.01 kpc
const dsCm = 0.01 * CM_PER_KPC;
const sbZ0 = integrateLogSB([{logEps: -24.0, dsCm}], 0.0);
const expectedZ0 = -24.0 + Math.log10(dsCm) - LOG10_4PI_ARCSEC2_PER_SR;
assert.ok(Math.abs(sbZ0 - expectedZ0) < 1e-9, `SB z=0 exact match: ${sbZ0} vs ${expectedZ0}`);

// Cosmological (1+z)^4 dimming check at z=1 (-4*log10(2) = -1.20412)
const sbZ1 = integrateLogSB([{logEps: -24.0, dsCm}], 1.0);
assert.ok(Math.abs((sbZ0 - sbZ1) - 4 * Math.log10(2)) < 1e-9, "exact (1+z)^4 dimming");

console.log(`PASS cloudy_emissivity_unit (Ha@1e4K=${haWarm.toFixed(2)}, OIII@1e4.8K=${oiiiWarm.toFixed(2)}, SB=${sbZ0.toFixed(2)} erg/s/cm^2/arcsec^2)`);
