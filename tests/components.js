// Component checks: each piece against something it must get right in isolation, on small hand-built
// worlds. Seconds, not minutes. Numbered as the verification tests in docs/crust_ale_method.tex.
// Run: npm run check
require('./harness').inPage(() => {
  const res = {};
  // A two-plate world split at the equator (plate 0 north), uniform zero-age ocean so that Φ̃ = 0
  // everywhere: the only load is the body force a check applies.
  const twoPlates = (seamType, N = 1200) => makeWorld({ seed: 9, N, contFrac: 0, plates: 2, kappaM: 0.1, meshMode: 'lagrangian',
    plateOf: (x, y, z) => z > 0 ? 0 : 1, seamOf: () => ({ type: seamType, sub: seamType === TRENCH ? 1 : -1 }),
    stateOf: () => ({ h: PHYS.h_oc, rho: PHYS.rho_oc, Th: 0 }) });
  const tight = (w, f) => { const t = RHEO.tol; RHEO.tol = 1e-9; w.U = null; f(w); RHEO.tol = t; };
  // seam pairs: closing rate (u1 − u2)·n and slip, in cm/yr
  const seamRates = w => { const X = w.X, out = [];
    for (let e = 0; e < w.sA.length; e++) { const a = w.sA[e], b = w.sB[e], mx = X[3*a] + X[3*b], my = X[3*a+1] + X[3*b+1], mz = X[3*a+2] + X[3*b+2];
      const ex = X[3*b] - X[3*a], ey = X[3*b+1] - X[3*a+1], ez = X[3*b+2] - X[3*a+2], l = Math.hypot(ex, ey, ez);
      let nx = ey*mz - ez*my, ny = ez*mx - ex*mz, nz = ex*my - ey*mx; const nl = Math.hypot(nx, ny, nz); nx /= nl; ny /= nl; nz /= nl;
      for (let k = 0; k < 2; k++) { const c1 = w.sC[4*e+k], c2 = w.sC[4*e+2+k], dx = w.U[3*c1] - w.U[3*c2], dy = w.U[3*c1+1] - w.U[3*c2+1], dz = w.U[3*c1+2] - w.U[3*c2+2];
        out.push({ closing: (dx*nx + dy*ny + dz*nz) * R_KM / 10, slip: (dx*ex + dy*ey + dz*ez) / l * R_KM / 10 }); } }
    return out; };
  const mean = a => a.reduce((s, v) => s + v, 0) / a.length, rms = a => Math.sqrt(mean(a.map(v => v * v)));
  // body forces, N/m²: toward or away from the equator, or shear along it
  const s0 = 3e5;
  const meridional = sign => (x, y, z, p) => { const k = sign * (p === 0 ? -1 : 1) * s0; return [-k * z * x, -k * z * y, k * (1 - z * z)]; };
  const zonal = (x, y, z, p) => { const k = (p === 0 ? 1 : -1) * s0; return [-k * y, k * x, 0]; };

  // One-sided normal coupling (§1.8): a trench absorbs convergence and resists opening; a ridge absorbs
  // opening and resists closing; a transform resists both. Closing rates are relative to the same push
  // with no coupling at all.
  {
    const free = {};
    for (const [name, sign] of [['converge', 1], ['diverge', -1]]) {
      const w = twoPlates(RIDGE), cn = RHEO.Cn_tf; RHEO.Cn_tf = 0; w.bodyForce = meridional(sign);
      tight(w, w => solveVelocity(w, 5000)); RHEO.Cn_tf = cn; free[name] = mean(seamRates(w).map(r => r.closing));
    }
    res.uncoupled_closing_cm_yr = { converge: +free.converge.toFixed(2), diverge: +free.diverge.toFixed(2) };
    res.oneSided = [];
    for (const [tname, type] of [['trench', TRENCH], ['ridge', RIDGE], ['transform', TRANSFORM]])
      for (const [name, sign] of [['converge', 1], ['diverge', -1]]) {
        const w = twoPlates(type); w.bodyForce = meridional(sign);
        tight(w, w => solveVelocity(w, 5000));
        const stiffFrac = mean([...w.contact]);
        res.oneSided.push({ seam: tname, push: name, closing_over_uncoupled: +(mean(seamRates(w).map(r => r.closing)) / free[name]).toFixed(4), stiff_fraction: +stiffFrac.toFixed(2), passes: w.solve.passes, unsettled: w.solve.changed });
      }
  }

  // Test 6 (transform): opening relative to slip under a pure shear push, as C_n grows.
  res.transform = [];
  for (const k of [1, 10, 100]) {
    const w = twoPlates(TRANSFORM), cn = RHEO.Cn_tf; w.bodyForce = zonal; RHEO.Cn_tf = cn * k;
    tight(w, w => solveVelocity(w, 5000)); RHEO.Cn_tf = cn;
    const r = seamRates(w);
    res.transform.push({ Cn_factor: k, rms_opening_over_rms_slip: +(rms(r.map(q => q.closing)) / rms(r.map(q => q.slip))).toExponential(2) });
  }

  // Test 3 on the sphere: a rigid rotation produces no strain in any triangle.
  {
    const w = makeWorld({ seed: 7, N: 2000, contFrac: 0.3, plates: 6, kappaM: 0.1, meshMode: 'lagrangian' });
    const om = [0.3, -0.2, 0.9], U = new Float64Array(3 * w.CM.length);
    for (let c = 0; c < w.CM.length; c++) { const m = w.CM[c], x = w.X[3*m], y = w.X[3*m+1], z = w.X[3*m+2]; U[3*c] = om[1]*z - om[2]*y; U[3*c+1] = om[2]*x - om[0]*z; U[3*c+2] = om[0]*y - om[1]*x; }
    res.rigidRotation = { max_strain_rate_over_omega: +(Math.max(...strainRates(w, U)) / Math.hypot(...om)).toExponential(1) };
  }

  // Drag-dominated limit (instantaneous form of test 4): negligible strength, one plate, uniform crust,
  // smooth heat deficit; the solve must return u = −∇Φ/β.
  res.dragLimit = [];
  for (const N of [1000, 4000]) {
    const ax = [0.48, 0.6, 0.64], Th0 = 4e4, f = (x, y, z) => Th0 * (1 + 0.5 * (ax[0]*x + ax[1]*y + ax[2]*z));
    const w = makeWorld({ seed: 3, N, contFrac: 0, plates: 1, kappaM: 0.1, meshMode: 'lagrangian', stateOf: (x, y, z) => ({ h: 7, rho: 2900, Th: f(x, y, z) }) });
    const st = RHEO.strength; RHEO.strength = 1e-9; tight(w, w => solveVelocity(w, 5000)); RHEO.strength = st;
    let num = 0, den = 0;
    for (let c = 0; c < w.CM.length; c++) {
      const m = w.CM[c], x = w.X[3*m], y = w.X[3*m+1], z = w.X[3*m+2], Th = w.Z[c] / w.A[c], ad = ax[0]*x + ax[1]*y + ax[2]*z;
      const dPhi = (phiCol(7, 2900, Th * 1.0001) - phiCol(7, 2900, Th * 0.9999)) / (Th * 0.0002);
      const k = -dPhi * (0.5 * Th0 / R_M) / dragBeta(7, Th) * SEC_MYR / R_M;
      const ex = [k * (ax[0] - ad*x), k * (ax[1] - ad*y), k * (ax[2] - ad*z)];
      for (let d = 0; d < 3; d++) { num += (w.U[3*c+d] - ex[d]) ** 2; den += ex[d] ** 2; }
    }
    res.dragLimit.push({ N, relative_rms_error: +Math.sqrt(num / den).toExponential(2) });
  }

  // GPE: ridge-to-ocean drop in Φ against classical ridge push (eq. ridgepush).
  res.ridgePush = [10, 40, 80].map(tau => { const Th = 2 * PHYS.dT * Math.sqrt(PHYS.kappa * tau / Math.PI);
    const Frp = G * PHYS.rho_m * PHYS.alpha * PHYS.dT * (1 + 2 * PHYS.rho_m * PHYS.alpha * PHYS.dT / (Math.PI * (PHYS.rho_m - PHYS.rho_w))) * PHYS.kappa * tau * 1e6;
    return { age_Myr: tau, dPhi_over_Frp: +((PHI_R - phiCol(PHYS.h_oc, PHYS.rho_oc, Th)) / Frp).toFixed(3) }; });

  // Test 1 (budgets): a small default world under the force balance for 10 Myr, through every mesh operation.
  {
    const sum = (w, k) => w[k].reduce((s, v) => s + v, 0);
    const w = makeWorld({ seed: 11, N: 1500, contFrac: 0.3, plates: 5, kappaM: 0.1, meshMode: 'lagrangian' });
    let worst = 0;
    while (w.t < 10 && !w.exhausted) { advance(w, 0.5); const B = w.budget;
      for (const e of [(sum(w, 'V') - (B.V0 + B.Vc - B.Vs)) / B.V0, (sum(w, 'M') - (B.M0 + B.Mc - B.Ms)) / B.M0, (sum(w, 'Z') - (B.Z0 + B.Zsrc - B.Zs)) / B.Z0]) worst = Math.max(worst, Math.abs(e)); }
    res.budgets = { Myr: +w.t.toFixed(1), stopped: w.exhausted ? w.fail : null, worst_relative_error: +worst.toExponential(1), ops: w.ops };
  }

  // Test 2 (free-stream preservation): uniform state, u = 0, smooth imposed motion of interior markers.
  res.freeStream = [];
  for (const maintain of [false, true]) {
    const w = makeWorld({ seed: 5, N: 2000, contFrac: 0, plates: 4, kappaM: 0, meshMode: 'lagrangian' });
    w.kinematic = true; w.omega = w.omega.map(() => [0, 0, 0]); w.noSources = true; w.noMaintain = !maintain;
    for (let c = 0; c < w.CM.length; c++) { const A = w.A[c]; w.V[c] = 7 * A; w.M[c] = 7 * 2900 * A; w.Z[c] = 5e4 * A; w.S[c] = 7 * A; }
    const a = [0.3, -0.5, 0.8], b = [0.6, 0.2, -0.1], s = 0.05 * w.L0;
    w.extraW = (x, y, z) => { const f = Math.sin(3 * (b[0]*x + b[1]*y + b[2]*z)); return [s * f * (a[1]*z - a[2]*y), s * f * (a[2]*x - a[0]*z), s * f * (a[0]*y - a[1]*x)]; };
    let maxErr = 0;
    while (w.t < 4 - 1e-9) { advance(w, Math.min(0.5, 4 - w.t)); for (let c = 0; c < w.CM.length; c++) maxErr = Math.max(maxErr, Math.abs(w.V[c] / w.A[c] - 7) / 7); }
    res.freeStream.push({ maintenance: maintain, max_relative_h_error: +maxErr.toExponential(2) });
  }
  return res;
});
