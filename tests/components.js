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
  const tight = (w, f) => { const t = RHEO.tol, p = RHEO.passes; RHEO.tol = 1e-9; RHEO.passes = 60; w.U = null; f(w); RHEO.tol = t; RHEO.passes = p; };
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

  // Test 7 (slab feedback): a stiff polar cap (plate 0) with a trench on its leading half and a ridge on its
  // trailing half, driven along x̂ by a body force plus its own ridge push, cooling frozen so the arriving
  // heat deficit stays Θ(60 Myr). Lumped rigid-plate prediction: v (βA_p + ξℓ − Gℓ) = F, with
  // ℓ = Σ ℓ_e (n_e·x̂)² over trench edges and G = τ_a g ρm α Θ_tr. Γ = Gℓ / (ξℓ + βA_p).
  {
    const capCos = Math.cos(25 * Math.PI / 180), Th60 = 2 * PHYS.dT * Math.sqrt(PHYS.kappa * 60 / Math.PI), f0 = 2e6;
    // Rigid plate, frictionless boundary, as the lumped prediction assumes: plates 100× stronger; seams stay fully
    // damaged (a near-zero fault-zone width keeps seams at s ≈ 0, tiny s_min) so that their shear coupling C_t = Ĥ(s)/w vanishes.
    const st = RHEO.strength, si = RHEO.stepIter, ta = RHEO.tau_a, sm = RHEO.smin, kd = RHEO.kd, wf = RHEO.wFault; RHEO.strength = 100; RHEO.stepIter = 400; RHEO.smin = 1e-5; RHEO.wFault = 1e-4;
    const make = () => { const w = makeWorld({ seed: 4, N: 1500, contFrac: 0, plates: 2, kappaM: 0.1, meshMode: 'lagrangian',
        plateOf: (x, y, z) => z > capCos ? 0 : 1,
        seamOf: (p1, p2, mx) => mx > 0 ? { type: TRENCH, sub: 0 } : { type: RIDGE, sub: -1 },
        stateOf: () => ({ h: PHYS.h_oc, rho: PHYS.rho_oc, Th: Th60 }) });
      w.noSources = true; w.noEvents = true; w.bodyForce = (x, y, z, p) => p === 0 ? [f0 * (1 - x * x), -f0 * x * y, -f0 * x * z] : [0, 0, 0]; return w; };
    const w0 = make(), X = w0.X, phiCap = phiCol(PHYS.h_oc, PHYS.rho_oc, Th60) - PHI_R; let bA = 0, F = 0, ell = 0;
    for (let c = 0; c < w0.CM.length; c++) if (w0.CP[c] === 0) { const m = w0.CM[c], x = X[3*m], a = w0.A[c] * R_M * R_M; bA += dragBeta(PHYS.h_oc, Th60) * a; F += f0 * (1 - x * x) * a; }
    for (let e = 0; e < w0.sA.length; e++) { const a = w0.sA[e], b = w0.sB[e];
      const mx = X[3*a] + X[3*b], my = X[3*a+1] + X[3*b+1], mz = X[3*a+2] + X[3*b+2], ex = X[3*b] - X[3*a], ey = X[3*b+1] - X[3*a+1], ez = X[3*b+2] - X[3*a+2];
      const nx = (ey*mz - ez*my) / Math.hypot(ey*mz - ez*my, ez*mx - ex*mz, ex*my - ey*mx), len = Math.hypot(ex, ey, ez) * R_M;
      if (w0.sT[e] === TRENCH) ell += len * nx * nx;
      else F += len * phiCap * (w0.CP[w0.sC[4*e]] === 0 ? nx : -nx); } // the cap's own free edge at its ridge: Φ̃ n_out
    const Gc = G * PHYS.rho_m * PHYS.alpha * Th60 * 1000;
    const capSpeed = w => { let v = 0, a = 0; for (let c = 0; c < w.CM.length; c++) if (w.CP[c] === 0) { const m = w.CM[c], x = w.X[3*m];
      v += w.A[c] * (w.U[3*c] - x * (w.U[3*c]*x + w.U[3*c+1]*w.X[3*m+1] + w.U[3*c+2]*w.X[3*m+2])); a += w.A[c] * (1 - x * x); } return v / a * R_KM / 10; };
    res.slabFeedback = [];
    for (const tau of [5, 15]) {
      RHEO.tau_a = tau; const Gl = tau * SEC_MYR * Gc * ell, Gam = Gl / (RHEO.xi * ell + bA);
      const w = make(); solveVelocity(w, 5000, 60); const v0 = capSpeed(w), sp = [];
      while (w.t < 40 && !w.exhausted) { advance(w, 5); sp.push(+capSpeed(w).toFixed(2)); }
      res.slabFeedback.push({ tau_a: tau, Gamma: +Gam.toFixed(2),
        v0_cm_yr: +v0.toFixed(2), v0_pred: +(F / (bA + RHEO.xi * ell) * SEC_MYR / 1e4).toFixed(2),
        v_terminal_pred: Gam < 1 ? +(F / (bA + RHEO.xi * ell - Gl) * SEC_MYR / 1e4).toFixed(2) : 'runaway',
        v_every_5_Myr: sp.join(' '), stopped: w.exhausted ? w.fail : null });
      if (tau === 5) { // the reservoir is a line integral owned by copies: mesh maintenance must not change its total
        const tot = () => w.Wl.reduce((q, v) => q + v, 0), before = tot(), o = { ...w.ops };
        w.L0 *= 3; for (let i = 0; i < 5; i++) maintain(w); w.L0 /= 3; // coarser target spacing forces chain merges along the trench
        res.slabReservoirAcrossMaintenance = { seam_merges: w.ops.seam - o.seam, relative_change: +((tot() - before) / before).toExponential(1) };
      }
    }
    RHEO.strength = st; RHEO.stepIter = si; RHEO.tau_a = ta; RHEO.smin = sm; RHEO.kd = kd; RHEO.wFault = wf;
  }

  // Power-law creep (eq. visc): one plate, uniform state, negligible drag, a torque-free load (the
  // tangential gradient of a·x). Doubling the load must multiply strain rates by 2ⁿ.
  {
    const ax = [0.48, 0.6, 0.64], b0 = RHEO.beta0; res.powerLaw = [];
    for (const n of [1, 3]) {
      const n0 = RHEO.n; RHEO.n = n; RHEO.beta0 = 1e10; const rmsE = [], info = [];
      for (const k of [1, 2]) {
        const w = makeWorld({ seed: 3, N: 1500, contFrac: 0, plates: 1, kappaM: 0.1, meshMode: 'lagrangian', stateOf: () => ({ h: 7, rho: 2900, Th: 6e4 }) });
        const f0 = 1e6 * k; w.bodyForce = (x, y, z) => { const d = ax[0]*x + ax[1]*y + ax[2]*z; return [f0 * (ax[0] - d*x), f0 * (ax[1] - d*y), f0 * (ax[2] - d*z)]; };
        tight(w, w => solveVelocity(w, 5000, 60));
        const e = strainRates(w, w.U); let q = 0; for (const v of e) q += v * v; rmsE.push(Math.sqrt(q / e.length));
        info.push({ passes: w.solve.passes, du: +w.solve.du.toExponential(0), median_edot_per_s: +([...e].sort((a, b) => a - b)[e.length >> 1] / SEC_MYR).toExponential(1) });
      }
      RHEO.n = n0; RHEO.beta0 = b0;
      res.powerLaw.push({ n, strain_rate_ratio: +(rmsE[1] / rmsE[0]).toFixed(3), expected: 2 ** n, picard: info });
    }
  }

  // Damage integrator (Alg. 6): uneven steps at fixed strain rate must match the exact solution
  // s(t) = s_eq + (s0 − s_eq) e^{−(kd ε̇ + kh) t}, and healing at zero strain rate.
  {
    const w = makeWorld({ seed: 3, N: 400, contFrac: 0, plates: 1, kappaM: 0.1, meshMode: 'lagrangian' });
    for (let c = 0; c < w.CM.length; c++) w.S[c] = 0.5 * w.V[c];
    const nc = w.CM.length, ed = Float64Array.from({ length: nc }, (_, c) => c % 2 ? 0 : 0.3); // 1/Myr
    let t = 0; for (const dt of [0.1, 2.3, 0.7, 5, 0.01, 3]) { damageSources(w, dt, ed); t += dt; }
    let err = 0; for (let c = 0; c < nc; c++) { const r = RHEO.kd * ed[c] + RHEO.kh, seq = RHEO.kh / r, ex = seq + (0.5 - seq) * Math.exp(-r * t); err = Math.max(err, Math.abs(w.S[c] / w.V[c] - ex)); }
    res.damageIntegrator = { Myr: t, max_abs_error: +err.toExponential(1) };
  }

  // Localization: one plate (no seam), zonal shear push (east in the north, west in the south). With damage
  // the shear should concentrate near the equator over time. Index = mean strain rate within 6° of the
  // equator over the global mean.
  res.localization = [];
  for (const damage of [false, true]) {
    const w = makeWorld({ seed: 6, N: 1500, contFrac: 0, plates: 1, kappaM: 0.1, meshMode: 'lagrangian', stateOf: () => ({ h: 7, rho: 2900, Th: 6e4 }) });
    w.noSources = true; w.noEvents = true; w.noDamage = !damage; const f0 = 4e5;
    w.bodyForce = (x, y, z) => { const k = z > 0 ? f0 : -f0; return [-k * y, k * x, 0]; };
    solveVelocity(w, RHEO.maxIter, RHEO.firstPasses);
    const index = () => { const e = copyStrainRates(w, w.U); let a = 0, na = 0, b = 0;
      for (let c = 0; c < w.CM.length; c++) { const z = w.X[3*w.CM[c]+2]; b += e[c]; if (Math.abs(z) < 0.1) { a += e[c]; na++; } } return (a / na) / (b / w.CM.length); };
    const idx = [+index().toFixed(2)]; let smin = 1;
    while (w.t < 20 && !w.exhausted) { advance(w, 5); idx.push(+index().toFixed(2)); }
    for (let c = 0; c < w.CM.length; c++) smin = Math.min(smin, w.S[c] / w.V[c]);
    res.localization.push({ damage, index_every_5_Myr: idx.join(' '), min_integrity: +smin.toFixed(2), stopped: w.exhausted ? w.fail : null });
  }

  // Test 8 (seam events), on two hemispheres split at the equator unless noted. Budgets are checked across
  // every event; Θ is frozen (no cooling) so each scenario keeps the buoyancy it was built with.
  {
    const sum = (w, k) => w[k].reduce((q, v) => q + v, 0);
    const budgetErr = w => { const B = w.budget; return Math.max(Math.abs((sum(w, 'V') - (B.V0 + B.Vc - B.Vs)) / B.V0), Math.abs((sum(w, 'M') - (B.M0 + B.Mc - B.Ms)) / B.M0), Math.abs((sum(w, 'Z') - (B.Z0 + B.Zsrc - B.Zs)) / B.Z0)); };
    const types = w => { const n = [0, 0, 0]; for (let e = 0; e < w.sA.length; e++) n[w.sT[e]]++; return { ridge: n[0], transform: n[1], trench: n[2] }; };
    const hemis = (seamType, stateOf) => { const w = makeWorld({ seed: 9, N: 1200, contFrac: 0, plates: 2, kappaM: 0.1, meshMode: 'lagrangian',
      plateOf: (x, y, z) => z > 0 ? 0 : 1, seamOf: () => ({ type: seamType, sub: seamType === TRENCH ? 1 : -1 }), stateOf }); w.noSources = true; return w; };
    const ThAge = a => 2 * PHYS.dT * Math.sqrt(PHYS.kappa * a / Math.PI), ocean = a => () => ({ h: PHYS.h_oc, rho: PHYS.rho_oc, Th: ThAge(a) });
    const run = (w, T) => { let worst = 0; while (w.t < T && !w.exhausted) { advance(w, 0.5); worst = Math.max(worst, budgetErr(w)); } return worst; };
    res.events = {};
    // Healing: an idle transform whose damage recovers (k_h raised to keep the test short) must weld shut.
    { const kh = RHEO.kh; RHEO.kh = 0.5; const w = hemis(TRANSFORM, ocean(60)); const s0 = types(w); const err = run(w, 10); RHEO.kh = kh;
      let hdev = 0; for (let c = 0; c < w.CM.length; c++) hdev = Math.max(hdev, Math.abs(w.V[c] / w.A[c] - PHYS.h_oc) / PHYS.h_oc);
      res.events.healing = { seams_before: s0, seams_after: types(w), copies_minus_markers: w.CM.length - w.X.length / 3, plates_after: new Set(w.CP).size, max_h_deviation: +hdev.toExponential(1), worst_budget: +err.toExponential(1), events: w.events, ran_to_Myr: +w.t.toFixed(1), stopped: w.exhausted ? w.fail : null }; }
    // Nucleation: one plate (no seam) sheared along the equator, with fast damage. Transforms should appear near the equator.
    { const kd = RHEO.kd; RHEO.kd = 40; const w = makeWorld({ seed: 6, N: 1500, contFrac: 0, plates: 1, kappaM: 0.1, meshMode: 'lagrangian', stateOf: ocean(60) });
      w.noSources = true; const f0 = 4e5; w.bodyForce = (x, y, z) => { const k = z > 0 ? f0 : -f0; return [-k * y, k * x, 0]; }; solveVelocity(w, RHEO.maxIter, RHEO.firstPasses);
      const err = run(w, 14); RHEO.kd = kd; // the window in which the band nucleates (from about 9 Myr)

      let near = 0; for (let e = 0; e < w.sA.length; e++) { const z = 0.5 * (w.X[3*w.sA[e]+2] + w.X[3*w.sB[e]+2]); if (Math.abs(z) < Math.sin(10 * Math.PI / 180)) near++; }
      let split = 0; for (let m = 0; m < w.X.length / 3; m++) if (w.cOff[m+1] - w.cOff[m] > 1) split++;
      // chains: connected components of the seam-edge graph; and the band's own width for scale
      const par = new Map(), find = m => { while (par.get(m) !== m) m = par.get(m); return m; };
      for (let e = 0; e < w.sA.length; e++) for (const m of [w.sA[e], w.sB[e]]) if (!par.has(m)) par.set(m, m);
      for (let e = 0; e < w.sA.length; e++) { const r1 = find(w.sA[e]), r2 = find(w.sB[e]); if (r1 !== r2) par.set(r1, r2); }
      const chains = new Set([...par.keys()].map(find)).size;
      let lat = 0; for (let e = 0; e < w.sA.length; e++) lat = Math.max(lat, Math.abs(Math.asin(0.5 * (w.X[3*w.sA[e]+2] + w.X[3*w.sB[e]+2]))) * 180 / Math.PI);
      let branch = 0; for (let m = 0; m < w.X.length / 3; m++) if (w.sCount[m] > 2) branch++;
      res.events.nucleation = { seam_edges: w.sA.length, ring_edges_estimate: Math.round(2 * Math.PI / w.L0), chains, branch_markers: branch, max_latitude_deg: +lat.toFixed(1), within_10deg_of_equator: near, split_markers: split, worst_budget: +err.toExponential(1), events: w.events, ran_to_Myr: +w.t.toFixed(1), stopped: w.exhausted ? w.fail : null }; }
    // Subduction initiation: a transform pushed shut between old north (80 Myr) and young south (5 Myr). The old side must subduct.
    { const w = hemis(TRANSFORM, (x, y, z, p) => ({ h: PHYS.h_oc, rho: PHYS.rho_oc, Th: ThAge(p === 0 ? 80 : 5) })); w.bodyForce = meridional(1);
      // initiation fires at τ = 2 Myr; the push never relents, so later evolution is not this rule's test
      const err = run(w, 3); let right = 0, tr = 0;
      for (let e = 0; e < w.sA.length; e++) if (w.sT[e] === TRENCH) { tr++; const o = w.sSub[e] === 2 ? 2 : 0, c = w.sC[4*e + o]; if (w.Z[c] / w.A[c] > ThAge(40)) right++; }
      res.events.subductionInitiation = { seams_after: types(w), trench_edges_with_old_side_subducting: `${right}/${tr}`, subducted_Mkm3: +(w.budget.Vs * R_KM * R_KM / 1e6).toFixed(2), worst_budget: +err.toExponential(1), events: w.events, ran_to_Myr: +w.t.toFixed(1), stopped: w.exhausted ? w.fail : null }; }
    // Rifting: a transform pulled apart. Thin (2 Myr) lithosphere must open a ridge and create crust; thick (80 Myr) must not.
    for (const age of [2, 80]) { const w = hemis(TRANSFORM, ocean(age)); w.bodyForce = meridional(-1); const err = run(w, 8);
      res.events['rifting_' + age + 'Myr'] = { seams_after: types(w), created_Mkm3: +(w.budget.Vc * R_KM * R_KM / 1e6).toFixed(2), worst_budget: +err.toExponential(1), events: w.events, ran_to_Myr: +w.t.toFixed(1), stopped: w.exhausted ? w.fail : null }; }
    // Ridge death: a ridge pushed shut (young, so it cannot start subducting) must become a transform.
    { const w = hemis(RIDGE, ocean(5)); w.bodyForce = meridional(1); const err = run(w, 6);
      res.events.ridgeDeath = { seams_after: types(w), worst_budget: +err.toExponential(1), events: w.events, ran_to_Myr: +w.t.toFixed(1), stopped: w.exhausted ? w.fail : null }; }
    // Breakoff: a stiff polar cap subducting at its leading edge with a continent riding toward the centre of the
    // trench, driven by a steady body force. Expected: breakoff where the continent arrives; with the push still
    // on and old ocean on the other side, subduction may restart with the opposite polarity.
    { const capCos = Math.cos(25 * Math.PI / 180), Th60 = ThAge(60), f0 = 2e6, st = RHEO.strength, si = RHEO.stepIter, sm = RHEO.smin, wf = RHEO.wFault;
      RHEO.strength = 100; RHEO.stepIter = 400; RHEO.smin = 1e-5; RHEO.wFault = 1e-4;
      const w = makeWorld({ seed: 4, N: 4000, contFrac: 0, plates: 2, kappaM: 0.1, meshMode: 'lagrangian', plateOf: (x, y, z) => z > capCos ? 0 : 1,
        seamOf: (p1, p2, mx) => mx > 0 ? { type: TRENCH, sub: 0 } : { type: RIDGE, sub: -1 },
        stateOf: (x, y, z, p) => p === 0 && x > 0.25 && x < 0.37 && Math.abs(y) < 0.12 ? { h: PHYS.h_cc, rho: PHYS.rho_cc, Th: PHYS.Theta_ss } : { h: PHYS.h_oc, rho: PHYS.rho_oc, Th: Th60 } });
      w.noSources = true; w.bodyForce = (x, y, z, p) => p === 0 ? [f0 * (1 - x * x), -f0 * x * y, -f0 * x * z] : [0, 0, 0]; solveVelocity(w, 5000, 30);
      const hist = []; let worst = 0;
      while (w.t < 12 && !w.exhausted) { advance(w, 3); worst = Math.max(worst, budgetErr(w));
        let front = 0; for (let c = 0; c < w.CM.length; c++) if (w.V[c] / w.A[c] > 20) front = Math.max(front, w.X[3*w.CM[c]]);
        const cen = []; for (let e = 0; e < w.sA.length; e++) { const a = w.sA[e], b = w.sB[e]; if (0.5 * (w.X[3*a] + w.X[3*b]) > 0.2 && Math.abs(0.5 * (w.X[3*a+1] + w.X[3*b+1])) < 0.12) {
          if (w.sT[e] !== TRENCH) { cen.push(['ridge', 'transform'][w.sT[e]]); continue; }
          const o = w.sSub[e] === 2 ? 2 : 0, c = w.sC[4*e + o]; cen.push(w.V[c] / w.A[c] > 20 ? 'trench(continent subducting)' : 'trench'); } }
        const counts = {}; for (const k of cen) counts[k] = (counts[k] || 0) + 1;
        hist.push(`${Math.round(w.t)} Myr: continent front x = ${front.toFixed(2)} (trench near 0.41); central seam edges ${JSON.stringify(counts)}`); }
      RHEO.strength = st; RHEO.stepIter = si; RHEO.smin = sm; RHEO.wFault = wf;
      res.events.breakoff = { history: hist, worst_budget: +worst.toExponential(1), events: w.events, stopped: w.exhausted ? w.fail : null }; }
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

  // Test 10 (mesh update), prescribed rigid rotations about the pole. Plate 1 is a 30° wedge of longitude
  // subducting at both its edges under plates 0 (west, moving east) and 2 (east, moving west), at 5 cm/yr
  // each at the equator; 0 and 2 spread apart at a ridge on the far side. Trenches move with their
  // overriding plates, so the wedge closes at 2ω and must be gone at 30°/2ω = 33.4 Myr: pinched first at
  // the poles, where it is narrowest, then zipped shut. Uniform h throughout (the ridge makes h_oc too), so
  // h must stay uniform through every zip and snap. Known defect, present before the repair existed: copies at
  // and next to the polar junctions drift by about 1.5% per Myr (a junction's mesh velocity cannot equal every
  // overriding plate's velocity, so a trench face declared material sweeps area). Uniformity is therefore
  // reported for interior copies, and for all copies separately; both until the wedge is gone, because
  // afterwards plates 0 and 2 converge rigidly across a transform, which piles crust up by construction.
  {
    const sum = (w, k) => w[k].reduce((q, v) => q + v, 0);
    const budgetErr = w => { const B = w.budget; return Math.max(Math.abs((sum(w, 'V') - (B.V0 + B.Vc - B.Vs)) / B.V0), Math.abs((sum(w, 'M') - (B.M0 + B.Mc - B.Ms)) / B.M0), Math.abs((sum(w, 'Z') - (B.Z0 + B.Zsrc - B.Zs)) / B.Z0)); };
    const half = 15 * Math.PI / 180, om = 50 / R_KM;
    const make = () => { const w = makeWorld({ seed: 3, N: 2000, contFrac: 0, plates: 3, kappaM: 0.1, meshMode: 'lagrangian',
        plateOf: (x, y) => { const f = Math.atan2(y, x); return f < -half ? 0 : f <= half ? 1 : 2; },
        seamOf: (p1, p2) => (p1 === 1 || p2 === 1) ? { type: TRENCH, sub: 1 } : { type: RIDGE, sub: -1 },
        stateOf: () => ({ h: PHYS.h_oc, rho: PHYS.rho_oc, Th: 1000 }) });
      w.kinematic = true; w.noSources = true; w.noEvents = true; w.noDamage = true;
      w.omega = [[0, 0, om], [0, 0, 0], [0, 0, -om]]; return w; };
    const w = make(); let opsRep = null;
    const area1 = () => { let a = 0; for (let c = 0; c < w.CM.length; c++) if (w.CP[c] !== 0 && w.CP[c] !== 2) a += w.A[c]; return a; };
    const A0 = area1(); let worst = 0, hdev = 0, hseam = 0, gone = null, wrongPol = 0;
    while (w.t < 45 && !w.exhausted) { advance(w, 1); worst = Math.max(worst, budgetErr(w));
      if (gone === null) { // before the wedge is gone; interior copies only (see below)
        for (let c = 0; c < w.CM.length; c++) if (w.cOff[w.CM[c] + 1] - w.cOff[w.CM[c]] === 1) hdev = Math.max(hdev, Math.abs(w.V[c] / w.A[c] - PHYS.h_oc) / PHYS.h_oc);
        for (let c = 0; c < w.CM.length; c++) hseam = Math.max(hseam, Math.abs(w.V[c] / w.A[c] - PHYS.h_oc) / PHYS.h_oc); }
      for (let e = 0; e < w.sA.length; e++) if (w.sT[e] === TRENCH) { const p = w.CP[w.sC[4*e + (w.sSub[e] === 1 ? 0 : 2)]]; if (p === 0 || p === 2) wrongPol++; }
      if (gone === null && area1() < 1e-6 * A0) gone = +w.t.toFixed(1); }
    // Repair alone: a second run in which the state is reset to uniform before every maintenance call and
    // measured right after it. This is the claim the remap makes, free of transport errors.
    let hrep = 0; { const w2 = make(); w2.noMaintain = true; let gone2 = false;
      while (w2.t < 45 && !w2.exhausted && !gone2) {
        for (let c = 0; c < w2.CM.length; c++) { const A = w2.A[c]; w2.V[c] = PHYS.h_oc * A; w2.M[c] = PHYS.rho_oc * PHYS.h_oc * A; w2.Z[c] = 1000 * A; w2.S[c] = PHYS.h_oc * A; }
        maintain(w2);
        for (let c = 0; c < w2.CM.length; c++) hrep = Math.max(hrep, Math.abs(w2.V[c] / w2.A[c] - PHYS.h_oc) / PHYS.h_oc, Math.abs(w2.Z[c] / w2.A[c] - 1000) / 1000);
        if (!step(w2, 1)) break; let a = 0; for (let c = 0; c < w2.CM.length; c++) if (w2.CP[c] !== 0 && w2.CP[c] !== 2) a += w2.A[c]; gone2 = a < 1e-6 * A0; }
      opsRep = w2.ops; }
    const n = [0, 0, 0]; for (let e = 0; e < w.sA.length; e++) n[w.sT[e]]++;
    res.meshUpdate = { wedge_gone_at_Myr: gone, predicted_Myr: +(2 * half / (2 * om)).toFixed(1), ran_to_Myr: +w.t.toFixed(1), stopped: w.exhausted ? w.fail : null,
      seams_after: { ridge: n[0], transform: n[1], trench: n[2] }, plates_after: new Set(w.CP).size, trench_edges_with_wrong_side_subducting: wrongPol,
      max_h_deviation_interior: +hdev.toExponential(1), max_h_deviation_any: +hseam.toExponential(1), repair_alone_max_deviation: +hrep.toExponential(1), repair_alone_ops: opsRep, worst_budget: +worst.toExponential(1), ops: w.ops };
  }
  // Seam refinement: a trench world, polarity known, a linear h, a slab reservoir on the subducting side;
  // halve the target spacing so that maintenance splits seam edges.
  {
    const w = makeWorld({ seed: 9, N: 1200, contFrac: 0, plates: 2, kappaM: 0.1, meshMode: 'lagrangian',
      plateOf: (x, y, z) => z > 0 ? 0 : 1, seamOf: () => ({ type: TRENCH, sub: 1 }), stateOf: x => ({ h: PHYS.h_oc * (1 + 0.3 * x), rho: PHYS.rho_oc, Th: 1000 }) });
    for (let c = 0; c < w.CM.length; c++) w.Wl[c] = w.CP[c] === 1 && w.cOff[w.CM[c] + 1] - w.cOff[w.CM[c]] > 1 ? 1e12 : 0;
    const sum = k => w[k].reduce((q, v) => q + v, 0), V0 = sum('V'), Z0 = sum('Z'), W0 = sum('Wl'), e0 = w.sA.length;
    w.L0 *= 0.55; maintain(w); w.L0 /= 0.55;
    let wrong = 0; for (let e = 0; e < w.sA.length; e++) { const c = w.sC[4*e + (w.sSub[e] === 1 ? 0 : 2)]; if (w.sT[e] !== TRENCH || w.CP[c] !== 1) wrong++; }
    res.seamSplit = { seam_edges: `${e0} → ${w.sA.length}`, splits: w.ops.seamSplit, rel_V: +(sum('V') / V0 - 1).toExponential(1), rel_Z: +(sum('Z') / Z0 - 1).toExponential(1), rel_slab: +(sum('Wl') / W0 - 1).toExponential(1), wrong_polarity: wrong };
  }

  // Test 2 (free-stream preservation): uniform state, u = 0, smooth imposed motion of interior markers.
  res.freeStream = [];
  // With seams too, the imposed motion is applied to seam and junction markers as well: every seam face then
  // moves through material that does not move, and must carry exactly the flux of the area it sweeps.
  for (const [maintain, seams] of [[false, false], [true, false], [false, true]]) {
    const w = makeWorld({ seed: 5, N: 2000, contFrac: 0, plates: 4, kappaM: 0, meshMode: 'lagrangian' });
    w.kinematic = true; w.omega = w.omega.map(() => [0, 0, 0]); w.noSources = true; w.noEvents = true; w.noMaintain = !maintain; w.extraWAll = seams; w.noDamage = true;
    for (let c = 0; c < w.CM.length; c++) { const A = w.A[c]; w.V[c] = 7 * A; w.M[c] = 7 * 2900 * A; w.Z[c] = 5e4 * A; w.S[c] = 7 * A; }
    const a = [0.3, -0.5, 0.8], b = [0.6, 0.2, -0.1], s = 0.05 * w.L0;
    w.extraW = (x, y, z) => { const f = Math.sin(3 * (b[0]*x + b[1]*y + b[2]*z)); return [s * f * (a[1]*z - a[2]*y), s * f * (a[2]*x - a[0]*z), s * f * (a[0]*y - a[1]*x)]; };
    let maxErr = 0;
    while (w.t < 4 - 1e-9) { advance(w, Math.min(0.5, 4 - w.t)); for (let c = 0; c < w.CM.length; c++) maxErr = Math.max(maxErr, Math.abs(w.V[c] / w.A[c] - 7) / 7); }
    res.freeStream.push({ maintenance: maintain, seam_markers_moved: seams, max_relative_h_error: +maxErr.toExponential(2) });
  }
  return res;
});
