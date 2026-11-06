// Loads index.html in headless Chromium (no animation loop) and evaluates a function inside it.
const { chromium } = require('playwright');
const fs = require('fs'), path = require('path');

const root = path.join(__dirname, '..');
const page = '<!doctype html><html><head><meta charset="utf-8"></head><body>' +
  fs.readFileSync(path.join(root, 'index.html'), 'utf8').replace('refreshLabels(); regenerate(); drawLegend(); requestAnimationFrame(frame);', 'refreshLabels(); drawLegend();') +
  '</body></html>';
const delaunator = fs.readFileSync(require.resolve('delaunator/delaunator.min.js'), 'utf8');

async function inPage(fn) {
  let browser;
  try { browser = await chromium.launch(); } catch (e) { browser = await chromium.launch({ executablePath: '/opt/pw-browsers/chromium' }); }
  const p = await browser.newPage();
  const errors = []; p.on('pageerror', e => errors.push(e.message));
  await p.route('**/delaunator.min.js', r => r.fulfill({ body: delaunator, contentType: 'application/javascript' }));
  await p.route('https://fonts.**', r => r.abort());
  await p.route('http://local/', r => r.fulfill({ body: page, contentType: 'text/html' }));
  await p.goto('http://local/'); await p.waitForTimeout(300);
  const out = await p.evaluate(fn);
  console.log(JSON.stringify(out, null, 1));
  console.log(errors.length ? 'page errors: ' + errors.join('; ') : 'no page errors');
  await browser.close();
}
module.exports = { inPage };
