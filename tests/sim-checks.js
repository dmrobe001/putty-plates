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
  const browser = await chromium.launch();
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

    // Test 1 (budgets) and mesh lifetime: default rules, several seeds, until the mesh needs a topology change.
    res.runs = [];
    for (const seed of [11, 22, 33, 44]) {
      const w = makeWorld({ seed, N: 8000, contFrac: 0.3, plates: 6, kappaM: 0.1, meshMode: 'lagrangian' });
      const t0 = performance.now(); let worst = { V: 0, M: 0, Z: 0 };
      while (!w.exhausted && w.t < 300) { advance(w, 1); const b = budget(w); for (const k in b) worst[k] = Math.max(worst[k], Math.abs(b[k])); }
      res.runs.push({ seed, lasted_Myr: +w.t.toFixed(1), stop: w.exhausted ? w.fail : 'ran to 300 Myr', ms_per_Myr: Math.round((performance.now() - t0) / w.t),
        worst_budget_error: Object.fromEntries(Object.entries(worst).map(([k, v]) => [k, +v.toExponential(1)])), ops: w.ops });
    }

    // Test 2 (free-stream preservation): uniform state, u = 0, a smooth imposed motion of interior markers,
    // no sources. Exact in the plane by the geometric conservation law. Reported separately without and
    // with mesh maintenance, whose remap is first order.
    res.freeStream = [];
    for (const maintain of [false, true]) {
      const w = makeWorld({ seed: 5, N: 8000, contFrac: 0, plates: 4, kappaM: 0, meshMode: 'lagrangian' });
      w.omega = w.omega.map(() => [0, 0, 0]); w.noSources = true; w.noMaintain = !maintain;
      for (let c = 0; c < w.CM.length; c++) { const A = w.A[c]; w.V[c] = 7 * A; w.M[c] = 7 * 2900 * A; w.Z[c] = 5e4 * A; w.S[c] = 7 * A; }
      const a = [0.3, -0.5, 0.8], b = [0.6, 0.2, -0.1], s = 0.05 * w.L0; // up to 5% of a cell per Myr
      w.extraW = (x, y, z) => { const f = Math.sin(3 * (b[0]*x + b[1]*y + b[2]*z)); return [s * f * (a[1]*z - a[2]*y), s * f * (a[2]*x - a[0]*z), s * f * (a[0]*y - a[1]*x)]; };
      let maxErr = 0;
      while (w.t < 4 - 1e-9) { advance(w, Math.min(0.5, 4 - w.t)); for (let c = 0; c < w.CM.length; c++) maxErr = Math.max(maxErr, Math.abs(w.V[c] / w.A[c] - 7) / 7); }
      res.freeStream.push({ maintenance: maintain, Myr: 4, max_relative_h_error: +maxErr.toExponential(2), ops: w.ops });
    }
    return res;
  });
  console.log(JSON.stringify(out, null, 1));
  console.log(errors.length ? 'page errors: ' + errors.join('; ') : 'no page errors');
  await browser.close();
})();
