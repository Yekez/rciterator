// using and quoting RC Concrete book from Ersoy, Ozcebe and Tankut
// page 100

/*
Initial values (defaults used in solver when not provided):

  Section (concrete):
    width  = 400 mm
    height = 400 mm
    layers = []   (steel layers set by caller; each layer: barCount, diameter mm, depth mm from top)

  Materials:
    fck   = 30 MPa          (UI default; solver uses fc = 0.85*fck if fck given, else fc = 25)
    eps0  = 0.002           (Hognestad strain at peak)
    epsCu = 0.0038          (ultimate concrete strain; Eurocode 0.0035)
    Es    = 200000 MPa
    fyk    = 420 MPa         

  M–κ sweep:
    epsMin  = 0.0002
    epsMax  = 0.02
    epsStep = 0.0002

  Equilibrium (bisection):
    c range   = 1 to height*1.5 mm
    errorLimit = max(|N_app|*0.001, 1) N
    maxIter  = 80
    nStrip   = 80 (concrete integration)

  Curve stop:
    when M < 0.9 * M_max (10% drop from max moment, as I remember something similar we discussed in the class)
*/


/*
Algorithm summary:
  For a given axial force N_app, we build the M–κ curve by sweeping the top-fibre concrete strain ε_ci
  from epsMin to epsMax in steps of epsStep. For each ε_ci:
    1. Find neutral axis depth c by bisection so that Fc + Fs = N_app (force equilibrium).
    2. Fc = integral of Hognestad σ(ε) over the compression zone (y = 0 to c); Fs = sum of steel
       layer forces from compatibility ε(d) = ε_ci*(c−d)/c and elastic–perfectly plastic σ(ε).
    3. Compute M = Mc + Ms about the section centroid; κ = ε_ci / c.
    4. Append (κ, M) to the curve. Use previous c to narrow the next bisection bracket (hint: max = c*(1+Δε/ε)).
  Stop when moment has dropped 10% from the maximum. Concrete: Hognestad (parabolic + descending).
  Steel: elastic–perfectly plastic. No confinement or strain hardening.
*/
/*
1. Choose the extreme fiber concrete strain
2. Assume a neutral axis depth, assume "c"
3. From the compatibility equation, compute the steel strains
4. Compute the steel stresses
5. Compute the steel forces F_si = A_si * sigma_si
6. Compute the concrete force. The resultant concrete force is equal to the volume under the stress distribution.
   (Hognestad curve integrated over compression zone; area × b_w.)
7. Check the force equilibrium using the equilibrium equation: F_c + sum(F_si) = N_app
8. If it is satisfied, continue. If it is not satisfied, update the neutral axis depth and repeat until it is satisfied.
9. Compute moment of internal forces about the centroid of the gross concrete section.
10. Compute the curvature from K = -epsilon_ci / c
11. Go to step 1 and take another value of epsilon_ci
*/

/*
alt codes:
σ (sigma): 229
ε (epsilon): 238
*/
(function (global) {
  'use strict';

  const RCSolver = {
    section: null,
    materials: null,

    setSection(opts) {
      this.section = opts;
      return this;
    },

    setMaterials(opts) {
      this.materials = opts;
      return this;
    },

    // Hognestad concrete stress (MPa). Compression positive.
    sigmaConcrete(eps) {
      const { fc, eps0, epsCu } = this._concreteParams();
      if (eps <= 0) return 0;
      if (eps >= epsCu) return 0; // beyond ultimate
      if (eps <= eps0) {
        // ascending (parabolic): σ = 0.85·fck · [2·(ε_c/ε_0) − (ε_c/ε_0)²]
        const r = eps / eps0;
        return fc * (2 * r - r * r);
      }
      // descending branch: σ = 0.85·fck − (0.15·0.85·fck)/(ε_cu−ε_0)·(ε_c−ε_0) = fc·(1 − 0.15·(ε−ε_0)/(ε_cu−ε_0))
      const denom = epsCu - eps0;
      if (denom <= 0) return fc;
      return fc * (1 - 0.15 * (eps - eps0) / denom);
    },

    _concreteParams() {
      const m = this.materials || {};
      const fc = m.fc != null ? m.fc : (m.fck != null ? 0.85 * m.fck : 25);
      const eps0 = m.eps0 != null ? m.eps0 : 0.002;
      const epsCu = m.epsCu != null ? m.epsCu : 0.0038; // Ersoy/Özcebe/Tankut; Eurocode uses 0.0035
      return { fc, eps0, epsCu };
    },

    // Steel stress (MPa). Tension positive, elastic-perfectly plastic.
    sigmaSteel(eps) {
      const m = this.materials || {};
      const Es = m.Es != null ? m.Es : 200000;
      const fy = m.fy != null ? m.fy : (m.fyk != null ? m.fyk : 420);
      const sig = Es * eps;
      if (sig >= fy) return fy;
      if (sig <= -fy) return -fy;
      return sig;
    },

    // Strain at depth y (from top, mm) given NA depth c and extreme top strain epsCi.
    strainAtDepth(y, c, epsCi) {
      if (c <= 0) return 0;
      return epsCi * (c - y) / c;
    },

    // Internal axial force (N) for given c and epsCi. Returns N_int = Fc + Fs (compression positive).
    // Used in equilibrium: we find c so that N_int = N_app.
    axialForce(epsCi, c) {
      const s = this.section;
      if (!s) return 0;
      const b = s.width || 400;
      const h = s.height || 400;

      // --- Concrete force Fc (compression zone only, y from 0 to c) ---
      // Integrate stress over depth: Fc = ∫ σ(ε(y)) * b dy. Strain is linear: ε(y) = epsCi * (c - y) / c.
      const nStrip = 80;
      const dy = h / nStrip;
      let Fc = 0;
      for (let i = 0; i < nStrip; i++) {
        const y = (i + 0.5) * dy;  // strip centre
        if (y >= c) break;         // only above neutral axis (compression)
        const eps = this.strainAtDepth(y, c, epsCi);
        const sig = this.sigmaConcrete(eps);  // Hognestad
        Fc += sig * b * dy;
      }

      // --- Steel force Fs (each layer) ---
      // Compatibility: strain at layer depth d = epsCi * (c - d) / c. Then F_si = σ(ε) * As (tension + or compression -).
      const layers = s.layers || [];
      let Fs = 0;
      for (const L of layers) {
        const n = L.barCount || 0;
        const d = L.diameter || 0;
        const depth = L.depth != null ? L.depth : 0;
        const As = n * (Math.PI * d * d / 4);
        if (As <= 0) continue;
        const eps = this.strainAtDepth(depth, c, epsCi);
        const sig = this.sigmaSteel(eps);  // elastic–perfectly plastic
        Fs += sig * As;
      }

      return Fc + Fs;
    },

    // Internal moment (N·mm) about section centroid. Concrete from Hognestad integration so post-peak softening gives descending M–κ.
    moment(epsCi, c) {
      const s = this.section;
      if (!s) return 0;
      const b = s.width || 400;
      const h = s.height || 400;
      const yRef = h / 2;
      const nStrip = 80;
      const dy = h / nStrip;
      let Mc = 0;
      for (let i = 0; i < nStrip; i++) {
        const y = (i + 0.5) * dy;
        if (y >= c) break;
        const eps = this.strainAtDepth(y, c, epsCi);
        const sig = this.sigmaConcrete(eps);
        const arm = yRef - y;
        Mc += sig * b * dy * arm;
      }
      const layers = s.layers || [];
      let Ms = 0;
      for (const L of layers) {
        const n = L.barCount || 0;
        const d = L.diameter || 0;
        const depth = L.depth != null ? L.depth : 0;
        const As = n * (Math.PI * d * d / 4);
        if (As <= 0) continue;
        const eps = this.strainAtDepth(depth, c, epsCi);
        const sig = this.sigmaSteel(eps);
        const arm = yRef - depth;
        Ms += sig * As * arm;
      }
      return Mc + Ms;
    },

    // Find c so that F_c + sum(F_si) = N_app. cPrev and nextMax optional 
    // (hint: use previous c and max = c*(1+Δε/ε) for next bracket). Returns { ok, c, M } or { ok: false }.
    equilibrium(N_app, epsCi, cPrev, nextMax) {
      const s = this.section;
      if (!s) return { ok: false };
      const h = (s.height || 400) * 1.5;
      let cLo = 1;
      let cHi = h;
      if (nextMax != null && nextMax > 0) {
        cHi = Math.min(h, nextMax);
      } else if (cPrev != null && cPrev > 0) {
        cLo = Math.max(1, cPrev * 0.25);
        cHi = Math.min(h, cPrev * 2.5);
      }
      const errorLimit = Math.max(Math.abs(N_app) * 0.001, 1); // hint: max(N × %error, 1 N)
      const maxIter = 80;
      for (let it = 0; it < maxIter; it++) {
        const c = 0.5 * (cLo + cHi);
        const N_int = this.axialForce(epsCi, c);
        const err = N_app - N_int; // hint: error = N - (Fc + Fs)
        if (Math.abs(err) <= errorLimit) {
          const M = this.moment(epsCi, c);
          return { ok: true, c, M };
        }
        if (err < 0) cHi = c;
        else cLo = c;
      }
      return { ok: false };
    },

    // Build M–κ curve per hint: ε_c^top 0.0002 to 0.02 step 0.0002; bisection for c; then max = c*(1+Δε/ε) for next bracket.
    // Stop when moment has dropped 10% from max (considerable drop → kill curve and reduce numerical noise).
    momentCurvatureCurve(N_app, options) {
      const opts = options || {};
      const epsMin = opts.epsMin ?? 0.0002;
      const epsMax = opts.epsMax ?? 0.02;
      const epsStep = opts.epsStep ?? 0.0002;
      const curve = [];
      let cPrev = null;
      let nextMax = null; // hint: max = c*(1 + Δε_c^top/ε_c^top) for next iteration
      let M_max = 0;
      for (let epsCi = epsMin; epsCi <= epsMax + epsStep * 0.5; epsCi += epsStep) {
        // Try to find c so that Fc + Fs = N_app, using previous c and nextMax to narrow the search.
        let res = this.equilibrium(N_app, epsCi, cPrev, nextMax);
        // If that failed, retry with full bisection range (no bracket) in case we were too narrow.
        if (!res.ok && (cPrev != null || nextMax != null)) res = this.equilibrium(N_app, epsCi, null, null);
        // Still no solution for this strain step; skip and go to next ε_ci.
        if (!res.ok) continue;
        // Track the maximum moment we've seen so far (for the 10% drop stop).
        if (res.M > M_max) M_max = res.M;
        // Stop adding points once moment has dropped 10% from the peak.
        if (M_max > 0 && res.M < 0.9 * M_max) break; // 10% drop from max moment -> stop
        // Remember this c so the next iteration can start its bisection near here.
        cPrev = res.c;
        // Upper bound for next bisection: c*(1 + Δε/ε) so we stay on the same branch.
        nextMax = res.c * (1 + epsStep / epsCi); // hint: set max for next ε_c^top
        // Curvature from compatibility: κ = ε_ci / c.
        const kappa = epsCi / res.c; // κ = ε_c^top / C
        // Append this (κ, M) point to the curve; we'll sort by κ later.
        curve.push({ kappa, M: res.M, epsCi, c: res.c });
      }
      curve.sort((a, b) => a.kappa - b.kappa);
      return curve;
    }
  };

  if (typeof module !== 'undefined' && module.exports) {
    module.exports = RCSolver;
  } else {
    global.RCSolver = RCSolver;
  }
})(typeof self !== 'undefined' ? self : this);
