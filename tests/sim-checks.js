// Headless checks of the simulation inside index.html.
// Run:  npm install && npx playwright install chromium && npm run check
// Prints measurements; it does not assert, because the rules are still provisional.
const { chromium } = require('playwright');
const fs = require('fs'), path = require('path');

const root = path.join(__dirname, '..');
const page = '<!doctype html><html><head><meta charset="utf-8"></head><body>' + fs.readFileSync(path.join(root, 'index.html'), 'utf8') + '</body></html>';
const delaunator = fs.readFileSync(require.resolve('delaunator/delaunator.min.js'), 'utf8');

(async () => {
  const browser = await chromium.launch({ args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'] });
  const p = await browser.newPage({ viewport: { width: 800, height: 600 } });
  const errors = []; p.on('pageerror', e => errors.push(e.message));
  await p.route('**/delaunator.min.js', r => r.fulfill({ body: delaunator, contentType: 'application/javascript' }));
  await p.route('https://fonts.**', r => r.abort());
  await p.route('http://local/', r => r.fulfill({ body: page, contentType: 'text/html' }));
  await p.goto('http://local/'); await p.waitForTimeout(300);

  const out = await p.evaluate(() => {
    setRunning(false);
    const res = {};
    // 1. Mesh: a closed triangulation of N points has exactly 2N - 4 triangles.
    res.mesh = { N: W.N, triangles: W.M.T.length / 3, expected: 2 * W.N - 4 };

    // 2. Rigid rotation: give every plate the same Euler vector; no rifts or sutures.
    //    The continental centroid should turn by |omega| t and continental area should not change.
    const base = { ...RULE }; RULE.riftTau = Infinity; RULE.sutureLen = 1e9;
    res.rigid = [];
    for (const vcm of [1, 5]) {
      seed = 12345; regenerate();
      const om = vcm * 10 / 6371, ax = [0.5, 0, Math.sqrt(0.75)];
      for (const q of W.plates) if (q) q.w = ax.map(v => v * om);
      const lon = () => { let x = 0, y = 0, z = 0, c = 0; for (let i = 0; i < W.N; i++) if (W.cont[i]) { x += W.M.P[3*i]; y += W.M.P[3*i+1]; z += W.M.P[3*i+2]; c++; }
        const d = x*ax[0] + y*ax[1] + z*ax[2]; x -= d*ax[0]; y -= d*ax[1]; z -= d*ax[2];
        const e1 = [ax[2], 0, -ax[0]], e2 = [ax[1]*e1[2] - ax[2]*e1[1], ax[2]*e1[0] - ax[0]*e1[2], ax[0]*e1[1] - ax[1]*e1[0]];
        return [Math.atan2(x*e2[0] + y*e2[1] + z*e2[2], x*e1[0] + y*e1[1] + z*e1[2]), c]; };
      const [a0, c0] = lon(); for (let s = 0; s < 50; s++) stepWorld(W, 2); resample(W, 0); const [a1, c1] = lon();
      let da = a1 - a0; if (da < 0) da += 2 * Math.PI;
      res.rigid.push({ speed_cm_per_yr: vcm, rotation_ratio_measured_over_expected: +(da / (om * 100)).toFixed(3), continental_area_change_pct: +(100 * (c1 - c0) / c0).toFixed(2) });
    }
    Object.assign(RULE, base);

    // 3. Long runs with default rules: plate count, continental share, continental elevation.
    res.history = [];
    for (const sd of [7919, 15838]) {
      seed = sd; regenerate(); const row = { seed: sd, samples: [] };
      for (let s = 0; s < 750; s++) { stepWorld(W, 2); if (s % 125 === 124) {
        let n = 0, c = 0; for (const q of W.plates) if (q && q.alive) { n++; c += q.cont; }
        const e = []; for (let i = 0; i < W.N; i++) if (W.cont[i]) e.push(elevOf(W, i)); e.sort((a, b) => a - b);
        row.samples.push({ t_Myr: Math.round(W.t), plates: n, continental: +(c / W.N).toFixed(3), cont_elev_p95_km: +e[Math.floor(e.length * .95)].toFixed(2) }); } }
      res.history.push(row);
    }
    return res;
  });
  console.log(JSON.stringify(out, null, 1));
  console.log(errors.length ? 'page errors: ' + errors.join('; ') : 'no page errors');
  await browser.close();
})();
