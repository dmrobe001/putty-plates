// Whole-planet smoke test: budgets (test 1), speeds, solver effort and how long the mesh lasts, for a few
// seeds. Slow headless (software-emulated graphics); run occasionally, don't iterate on it.
// Run: npm run check:planet
require('./harness').inPage(() => {
  const sum = (w, k) => w[k].reduce((s, v) => s + v, 0);
  const budget = w => { const B = w.budget; return {
    V: (sum(w, 'V') - (B.V0 + B.Vc - B.Vs)) / B.V0, M: (sum(w, 'M') - (B.M0 + B.Mc - B.Ms)) / B.M0, Z: (sum(w, 'Z') - (B.Z0 + B.Zsrc - B.Zs)) / B.Z0 }; };
  const runs = [];
  for (const seed of [11, 22, 33, 44]) {
    const w = makeWorld({ seed, N: 8000, contFrac: 0.3, plates: 6, kappaM: 0.1, meshMode: 'lagrangian' });
    const t0 = performance.now(); let worst = { V: 0, M: 0, Z: 0 }, iters = 0, steps = 0; const speed = [];
    while (!w.exhausted && w.t < 300) {
      const n = advance(w, 1); steps += n; iters += w.solve.iters * n;
      const b = budget(w); for (const k in b) worst[k] = Math.max(worst[k], Math.abs(b[k]));
      if (Math.abs(w.t - 10) < 0.5 && !speed.length) { for (let c = 0; c < w.CM.length; c++) speed.push(Math.hypot(w.U[3*c], w.U[3*c+1], w.U[3*c+2]) * R_KM / 10); speed.sort((a, b) => a - b); }
    }
    runs.push({ seed, lasted_Myr: +w.t.toFixed(1), stop: w.exhausted ? w.fail : 'ran to 300 Myr', ms_per_Myr: Math.round((performance.now() - t0) / w.t),
      mean_pcg_iters_per_step: Math.round(iters / steps), speed_at_10Myr_cm_yr: speed.length ? { median: +speed[speed.length >> 1].toFixed(1), p95: +speed[Math.floor(0.95 * speed.length)].toFixed(1), max: +speed[speed.length - 1].toFixed(1) } : null,
      worst_budget_error: Object.fromEntries(Object.entries(worst).map(([k, v]) => [k, +v.toExponential(1)])), ops: w.ops });
  }
  return runs;
});
