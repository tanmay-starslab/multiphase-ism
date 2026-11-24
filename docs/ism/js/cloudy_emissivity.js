// CLOUDY c17.03 (FG11 UVB, GASS10 solar, self-shielded) Line Emissivity LUT loader and evaluator.
// Lazy-imported by viewer3d.js ONLY when Emission (SB) mode is activated, so default page loads incur zero bytes.
// Physical output: log10 SB [erg s^-1 cm^-2 arcsec^-2] =
//   EPS_REF + log10(sum(10^(logEps - EPS_REF) * ds_cm)) - 11.7280614585 - 4*log10(1 + z).

export const DEFAULT_EMISSION_LINE = "Ha_6563";

export const EPS_REF = -25.0;
export const LOG10_4PI_ARCSEC2_PER_SR = 11.7280614585;
export const CM_PER_PC = 3.08567758e18;
export const CM_PER_KPC = 3.08567758e21;
export const LOG10_CM_PER_KPC = 21.489350545189;

export const CLOUDY_META_DEFAULTS = Object.freeze({
  source: "CLOUDY c17.03 (FG11)",
  redshift: 0.0,
  units: "log10(epsilon [erg s^-1 cm^-3])",
  sb_units: "log10(SB [erg s^-1 cm^-2 arcsec^-2])",
  log10_4pi_arcsec2_per_sr: LOG10_4PI_ARCSEC2_PER_SR,
  eps_min: -48.0,
  eps_max: -12.0,
  shape: [12, 56, 6, 61],
  nH_min: -7.0,
  nH_max: 4.0,
  nH_len: 56,
  Z_min: -3.0,
  Z_max: 1.0,
  Z_len: 6,
  T_min: 3.0,
  T_max: 9.0,
  T_len: 61
});

export const CLOUDY_LINES = Object.freeze([
  {id: "Ha_6563", label: "H\u03b1 6563 \u00c5", short: "H\u03b1", wave_A: 6562.81, band: "Optical", sb_range: [-21.5, -14.5]},
  {id: "Lya_1216", label: "Ly\u03b1 1216 \u00c5", short: "Ly\u03b1", wave_A: 1215.67, band: "FUV", sb_range: [-21.0, -14.0]},
  {id: "OII_3727", label: "[O II] 3727/29 \u00c5", short: "[O II]", wave_A: 3727.42, band: "Optical", sb_range: [-22.0, -15.0]},
  {id: "OIII_5007", label: "[O III] 5007 \u00c5", short: "[O III]", wave_A: 5006.84, band: "Optical", sb_range: [-22.0, -15.0]},
  {id: "NII_6583", label: "[N II] 6583 \u00c5", short: "[N II]", wave_A: 6583.45, band: "Optical", sb_range: [-22.5, -15.5]},
  {id: "MgII_2798", label: "Mg II 2796/03 \u00c5", short: "Mg II", wave_A: 2798.0, band: "NUV", sb_range: [-22.5, -15.5]},
  {id: "CIII_977", label: "C III 977 \u00c5", short: "C III", wave_A: 977.02, band: "FUV", sb_range: [-22.5, -16.0]},
  {id: "CIV_1548", label: "C IV 1548/51 \u00c5", short: "C IV", wave_A: 1549.06, band: "FUV", sb_range: [-23.0, -16.5]},
  {id: "OVI_1032", label: "O VI 1032/38 \u00c5", short: "O VI", wave_A: 1033.82, band: "FUV", sb_range: [-23.0, -16.5]},
  {id: "OVII_21.6", label: "O VII 21.6 \u00c5", short: "O VII", wave_A: 21.602, band: "X-ray", sb_range: [-24.0, -17.5]},
  {id: "OVIII_18.97", label: "O VIII 18.97 \u00c5", short: "O VIII", wave_A: 18.971, band: "X-ray", sb_range: [-24.0, -17.0]},
  {id: "FeXXV_1.85", label: "Fe XXV 1.85 \u00c5", short: "Fe XXV", wave_A: 1.8504, band: "X-ray", sb_range: [-24.5, -16.5]}
]);

const LINE_BY_ID = new Map(CLOUDY_LINES.map((l, idx) => [l.id, {...l, index: idx}]));

export function getEmissionLine(id) {
  return LINE_BY_ID.get(id) || LINE_BY_ID.get(DEFAULT_EMISSION_LINE);
}

/**
 * Trilinearly interpolate log10(epsilon [erg s^-1 cm^-3]) from the uint8 CLOUDY LUT
 * for a given line index and physical (lognH [cm^-3], logZ [Zsun], logT [K]).
 */
export function trilinearLogEps(raw, meta, lineIndex, lognH, logZ, logT) {
  const nNH = meta.nH_len, nZ = meta.Z_len, nT = meta.T_len;
  const stride = nNH * nZ * nT;
  const base = lineIndex * stride;
  const clamp = (x, lo, hi) => (x < lo ? lo : x > hi ? hi : x);
  const fnH = clamp((lognH - meta.nH_min) / (meta.nH_max - meta.nH_min), 0, 1) * (nNH - 1);
  const fZ = clamp((logZ - meta.Z_min) / (meta.Z_max - meta.Z_min), 0, 1) * (nZ - 1);
  const fT = clamp((logT - meta.T_min) / (meta.T_max - meta.T_min), 0, 1) * (nT - 1);

  const i0 = Math.min(nNH - 1, Math.floor(fnH)), i1 = Math.min(nNH - 1, i0 + 1), wnH = fnH - i0;
  const j0 = Math.min(nZ - 1, Math.floor(fZ)), j1 = Math.min(nZ - 1, j0 + 1), wZ = fZ - j0;
  const k0 = Math.min(nT - 1, Math.floor(fT)), k1 = Math.min(nT - 1, k0 + 1), wT = fT - k0;

  const at = (i, j, k) => raw[base + (i * nZ + j) * nT + k] / 255;
  const c00 = at(i0, j0, k0) * (1 - wT) + at(i0, j0, k1) * wT;
  const c01 = at(i0, j1, k0) * (1 - wT) + at(i0, j1, k1) * wT;
  const c10 = at(i1, j0, k0) * (1 - wT) + at(i1, j0, k1) * wT;
  const c11 = at(i1, j1, k0) * (1 - wT) + at(i1, j1, k1) * wT;
  const c0 = c00 * (1 - wZ) + c01 * wZ;
  const c1 = c10 * (1 - wZ) + c11 * wZ;
  const q = c0 * (1 - wnH) + c1 * wnH;
  return meta.eps_min + q * (meta.eps_max - meta.eps_min);
}

/**
 * Compute exact log10 SB [erg s^-1 cm^-2 arcsec^-2] from an array of ray samples
 * using float32-safe shifted log-accumulation (EPS_REF = -25.0).
 * @param {Array<{logEps: number, dsCm: number}>} samples
 * @param {number} redshift
 * @returns {number}
 */
export function integrateLogSB(samples, redshift = 0.0) {
  let acc = 0.0;
  for (let i = 0; i < samples.length; i++) {
    const s = samples[i];
    if (!(s.dsCm > 0)) continue;
    acc += Math.pow(10, s.logEps - EPS_REF) * s.dsCm;
  }
  if (!(acc > 0)) return -1e30;
  const dimZ = 4.0 * Math.log10(1.0 + Math.max(0.0, redshift));
  return EPS_REF + Math.log10(acc) - LOG10_4PI_ARCSEC2_PER_SR - dimZ;
}

/**
 * Build a Three.js Data3DTexture (width=nT, height=nZ, depth=nNH, RedFormat, LinearFilter)
 * for a single CLOUDY emission line slice so WebGL2 performs exact hardware trilinear lookup.
 */
export function buildLineTexture(THREE, raw, meta, lineId) {
  const line = getEmissionLine(lineId);
  const nNH = meta.nH_len, nZ = meta.Z_len, nT = meta.T_len;
  const stride = nNH * nZ * nT;
  const offset = line.index * stride;
  const slice = raw.slice(offset, offset + stride);
  const tex = new THREE.Data3DTexture(slice, nT, nZ, nNH);
  tex.format = THREE.RedFormat;
  tex.type = THREE.UnsignedByteType;
  tex.minFilter = THREE.LinearFilter;
  tex.magFilter = THREE.LinearFilter;
  tex.wrapS = THREE.ClampToEdgeWrapping;
  tex.wrapT = THREE.ClampToEdgeWrapping;
  tex.wrapR = THREE.ClampToEdgeWrapping;
  tex.unpackAlignment = 1;
  tex.generateMipmaps = false;
  tex.needsUpdate = true;
  return tex;
}

let cachedPromise = null;
let cachedBase = null;

export function loadCloudyEmissivity({baseUrl = "data/", fetchFn = globalThis.fetch} = {}) {
  const normBase = String(baseUrl || "data/").replace(/\/?$/, "/");
  if (cachedPromise && cachedBase === normBase) return cachedPromise;
  cachedBase = normBase;
  cachedPromise = (async () => {
    const [metaRes, binRes] = await Promise.all([
      fetchFn(`${normBase}cloudy_emissivity_meta.json`),
      fetchFn(`${normBase}cloudy_emissivity_z0.bin`)
    ]);
    if (!metaRes.ok) throw new Error(`CLOUDY meta HTTP ${metaRes.status}`);
    if (!binRes.ok) throw new Error(`CLOUDY LUT HTTP ${binRes.status}`);
    const meta = await metaRes.json();
    const raw = new Uint8Array(await binRes.arrayBuffer());
    const expected = meta.shape.reduce((a, b) => a * b, 1);
    if (raw.length !== expected) {
      throw new Error(`CLOUDY LUT size mismatch: got ${raw.length}, expected ${expected}`);
    }
    const lineIndexMap = new Map(meta.lines.map((l, i) => [l.id, i]));
    const indexOf = id => (lineIndexMap.has(id) ? lineIndexMap.get(id) : 0);
    return {
      meta,
      raw,
      lines: meta.lines,
      getLine(id) {
        const idx = indexOf(id);
        return {...meta.lines[idx], index: idx};
      },
      getLineSlice(id) {
        const idx = indexOf(id);
        const stride = meta.nH_len * meta.Z_len * meta.T_len;
        return raw.subarray(idx * stride, (idx + 1) * stride);
      },
      evalLogEps(id, lognH, logZ, logT) {
        return trilinearLogEps(raw, meta, indexOf(id), lognH, logZ, logT);
      },
      makeLineTexture(THREE, id) {
        return buildLineTexture(THREE, raw, meta, id);
      }
    };
  })().catch(err => {
    cachedPromise = null;
    cachedBase = null;
    throw err;
  });
  return cachedPromise;
}
