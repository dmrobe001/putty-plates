// Headless checks of the simulation inside index.html.
// Run:  npm install && npx playwright install chromium && npm run check
// Prints measurements. The verification tests are numbered as in crust_ale_method.tex.
const { chromium } = require('playwright');
const fs = require('fs'), path = require('path');

const root = path.join(__dirname, '..');
const page = '<!doctype html><html><head><meta charset="utf-8"></head><body>' +
  fs.readFileSync(path.join(root, 'index.html'), 'utf8').replace('refreshLabels(); regenerate(); drawLegend(); requestAnimationFrame(frame);', 'refreshLabels(); drawLegend();') +
  '</body></html>';
const delaunator = fs.readFileSync(require.resolve('delaunator/delaunator.min.js'), 'utf8');

(async () => {
  let browser; try { browser = await chromium.launch(); } catch (e) { browser = await chromium.launch({ executablePath: '/opt/pw-browsers/chromium' }); }
  const p = await browser.newPage();
  const errors = []; p.on('pageerror', e => errors.push(e.message));
  await p.route('**/delaunator.min.js', r => r.fulfill({ body: delaunator, contentType: 'application/javascript' }));
  await p.route('https://fonts.**', r => r.abort());
  await p.route('http://local/', r => r.fulfill({ body: page, contentType: 'text/html' }));
  await p.goto('http://local/'); await p.waitForTimeout(300);

  const out = await p.evaluate(() => {
    const res = {};
    const sum = (w, k) => w[k].reduce((s, v) => s + v, 0);
    const budget = w => { const B = w.budget; return {
      V: (sum(w, 'V') - (B.V0 + B.Vc - B.Vs)) / B.V0,
      M: (sum(w, 'M') - (B.M0 + B.Mc - B.Ms)) / B.M0,
      Z: (sum(w, 'Z') - (B.Z0 + B.Zsrc - B.Zs)) / B.Z0 }; };

    // Test 1 (budgets) and mesh lifetime: force-balance velocity, several seeds, until the mesh needs a topology change.
    res.runs = [];
    for (const seed of [11, 22, 33, 44]) {
      const w = makeWorld({ seed, N: 8000, contFrac: 0.3, plates: 6, kappaM: 0.1, meshMode: 'lagrangian' });
      const t0 = performance.now(); let worst = { V: 0, M: 0, Z: 0 }, iters = 0, steps = 0; const speed = [];
      while (!w.exhausted && w.t < 300) {
        const n = advance(w, 1); steps += n; iters += w.solve.iters * n;
        const b = budget(w); for (const k in b) worst[k] = Math.max(worst[k], Math.abs(b[k]));
        if (Math.abs(w.t - 10) < 0.5 && !speed.length) { for (let c = 0; c < w.CM.length; c++) speed.push(Math.hypot(w.U[3*c], w.U[3*c+1], w.U[3*c+2]) * R_KM / 10); speed.sort((a, b) => a - b); }
      }
      res.runs.push({ seed, lasted_Myr: +w.t.toFixed(1), stop: w.exhausted ? w.fail : 'ran to 300 Myr', ms_per_Myr: Math.round((performance.now() - t0) / w.t),
        mean_pcg_iters: Math.round(iters / steps), speed_at_10Myr_cm_yr: speed.length ? { median: +speed[speed.length >> 1].toFixed(1), p95: +speed[Math.floor(0.95 * speed.length)].toFixed(1), max: +speed[speed.length - 1].toFixed(1) } : null,
        worst_budget_error: Object.fromEntries(Object.entries(worst).map(([k, v]) => [k, +v.toExponential(1)])), ops: w.ops });
    }

    // Test 2 (free-stream preservation): uniform state, u = 0, a smooth imposed motion of interior markers,
    // no sources. Exact in the plane by the geometric conservation law. Reported separately without and
    // with mesh maintenance, whose remap is first order.
    res.freeStream = [];
    for (const maintain of [false, true]) {
      const w = makeWorld({ seed: 5, N: 8000, contFrac: 0, plates: 4, kappaM: 0, meshMode: 'lagrangian' });
      w.kinematic = true; w.omega = w.omega.map(() => [0, 0, 0]); w.noSources = true; w.noMaintain = !maintain;
      for (let c = 0; c < w.CM.length; c++) { const A = w.A[c]; w.V[c] = 7 * A; w.M[c] = 7 * 2900 * A; w.Z[c] = 5e4 * A; w.S[c] = 7 * A; }
      const a = [0.3, -0.5, 0.8], b = [0.6, 0.2, -0.1], s = 0.05 * w.L0; // up to 5% of a cell per Myr
      w.extraW = (x, y, z) => { const f = Math.sin(3 * (b[0]*x + b[1]*y + b[2]*z)); return [s * f * (a[1]*z - a[2]*y), s * f * (a[2]*x - a[0]*z), s * f * (a[0]*y - a[1]*x)]; };
      let maxErr = 0;
      while (w.t < 4 - 1e-9) { advance(w, Math.min(0.5, 4 - w.t)); for (let c = 0; c < w.CM.length; c++) maxErr = Math.max(maxErr, Math.abs(w.V[c] / w.A[c] - 7) / 7); }
      res.freeStream.push({ maintenance: maintain, Myr: 4, max_relative_h_error: +maxErr.toExponential(2), ops: w.ops });
    }

    // Test 3 on the sphere: a rigid rotation must produce no strain in any triangle (each flat triangle
    // sees a rigid motion of 3-D space projected onto its plane).
    {
      const w = makeWorld({ seed: 7, N: 8000, contFrac: 0.3, plates: 6, kappaM: 0.1, meshMode: 'lagrangian' });
      const om = [0.3, -0.2, 0.9], U = new Float64Array(3 * w.CM.length);
      for (let c = 0; c < w.CM.length; c++) { const m = w.CM[c], x = w.X[3*m], y = w.X[3*m+1], z = w.X[3*m+2]; U[3*c] = om[1]*z - om[2]*y; U[3*c+1] = om[2]*x - om[0]*z; U[3*c+2] = om[0]*y - om[1]*x; }
      res.rigidRotation = { max_strain_rate_over_omega: +(Math.max(...strainRates(w, U)) / Math.hypot(...om)).toExponential(1) };
    }

    // Drag-dominated limit (the instantaneous form of test 4): negligible strength, one plate, uniform crust,
    // smooth heat deficit. The solve must return u = −∇Φ/β.
    res.dragLimit = [];
    for (const N of [2000, 8000, 32000]) {
      const w = makeWorld({ seed: 3, N, contFrac: 0, plates: 1, kappaM: 0.1, meshMode: 'lagrangian' });
      const ax = [0.48, 0.6, 0.64], Th0 = 4e4, f = (x, y, z) => Th0 * (1 + 0.5 * (ax[0]*x + ax[1]*y + ax[2]*z));
      for (let c = 0; c < w.CM.length; c++) { const m = w.CM[c], A = w.A[c]; w.V[c] = 7 * A; w.M[c] = 7 * 2900 * A; w.S[c] = 7 * A; w.Z[c] = A * f(w.X[3*m], w.X[3*m+1], w.X[3*m+2]); }
      const s0 = RHEO.strength; RHEO.strength = 1e-9; RHEO.tol = 1e-10; w.U = null; solveVelocity(w, 5000); RHEO.strength = s0; RHEO.tol = 1e-4;
      let num = 0, den = 0;
      for (let c = 0; c < w.CM.length; c++) {
        const m = w.CM[c], x = w.X[3*m], y = w.X[3*m+1], z = w.X[3*m+2], Th = w.Z[c] / w.A[c], ad = ax[0]*x + ax[1]*y + ax[2]*z;
        const dPhi = (phiCol(7, 2900, Th * 1.0001) - phiCol(7, 2900, Th * 0.9999)) / (Th * 0.0002); // ∂Φ/∂Θ
        const g = 0.5 * Th0 / R_M, k = -dPhi * g / dragBeta(7, Th) * SEC_MYR / R_M; // u in rad/Myr
        const ex = [k * (ax[0] - ad*x), k * (ax[1] - ad*y), k * (ax[2] - ad*z)];
        for (let d = 0; d < 3; d++) { num += (w.U[3*c+d] - ex[d]) ** 2; den += ex[d] ** 2; }
      }
      res.dragLimit.push({ N, relative_rms_error: +Math.sqrt(num / den).toExponential(2) });
    }

    // Test 6 (transform): normal opening relative to tangential slip at transform copies, as C_n grows.
    res.transform = [];
    for (const k of [1, 10, 100]) {
      const w = makeWorld({ seed: 11, N: 8000, contFrac: 0.3, plates: 6, kappaM: 0.1, meshMode: 'lagrangian' });
      const cn = RHEO.Cn_tf; RHEO.Cn_tf = cn * k; RHEO.tol = 1e-9; solveVelocity(w, 5000); RHEO.Cn_tf = cn; RHEO.tol = 1e-4;
      let sn = 0, st = 0;
      for (let e = 0; e < w.sA.length; e++) { if (w.sT[e] !== TRANSFORM) continue;
        const a = w.sA[e], b = w.sB[e], X = w.X, mx = X[3*a] + X[3*b], my = X[3*a+1] + X[3*b+1], mz = X[3*a+2] + X[3*b+2];
        const ex = X[3*b] - X[3*a], ey = X[3*b+1] - X[3*a+1], ez = X[3*b+2] - X[3*a+2], l = Math.hypot(ex, ey, ez);
        let nx = ey*mz - ez*my, ny = ez*mx - ex*mz, nz = ex*my - ey*mx; const nl = Math.hypot(nx, ny, nz); nx /= nl; ny /= nl; nz /= nl;
        for (const [p, q] of [[0, 2], [1, 3]]) { const c1 = w.sC[4*e+p], c2 = w.sC[4*e+q];
          const dx = w.U[3*c2] - w.U[3*c1], dy = w.U[3*c2+1] - w.U[3*c1+1], dz = w.U[3*c2+2] - w.U[3*c1+2], dn = dx*nx + dy*ny + dz*nz, dt = (dx*ex + dy*ey + dz*ez) / l;
          sn += dn * dn; st += dt * dt; } }
      res.transform.push({ Cn_factor: k, rms_opening_over_rms_slip: +Math.sqrt(sn / st).toExponential(2) });
    }

    // Gravitational potential energy: ridge-to-ocean drop in Φ against the classical ridge push (eq. ridgepush).
    res.ridgePush = [10, 40, 80].map(tau => { const Th = 2 * PHYS.dT * Math.sqrt(PHYS.kappa * tau / Math.PI);
      const Frp = G * PHYS.rho_m * PHYS.alpha * PHYS.dT * (1 + 2 * PHYS.rho_m * PHYS.alpha * PHYS.dT / (Math.PI * (PHYS.rho_m - PHYS.rho_w))) * PHYS.kappa * tau * 1e6;
      return { age_Myr: tau, dPhi_over_Frp: +((PHI_R - phiCol(PHYS.h_oc, PHYS.rho_oc, Th)) / Frp).toFixed(3) }; });
    return res;
  });
  console.log(JSON.stringify(out, null, 1));
  console.log(errors.length ? 'page errors: ' + errors.join('; ') : 'no page errors');
  await browser.close();
})();
