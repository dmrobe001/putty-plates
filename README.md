# Putty Plates

A plate-tectonics toy for generating Earth-adjacent worlds for a game. The goal is diversity that comes from a few rules acting over a long history, not from layers of procedural noise.

The whole app is one file, `index.html`. It loads one script from a content delivery network (CDN), [Delaunator](https://github.com/mapbox/delaunator), and needs a browser with WebGL2. It is written as a page fragment (no `<!doctype>`), because it is published as a claude.ai Artifact. Opened directly from disk it renders in quirks mode, which this layout tolerates.

## Controls

| Control | Effect |
| --- | --- |
| Time step | Million years (Myr) simulated per animation frame; takes effect immediately |
| Vertices | Mesh resolution; regenerates the world when you release the slider |
| Continental area fraction | Initial share of the surface that is continental crust; regenerates |
| Plates | Initial plate count; regenerates |
| Run / Pause (Space) | Pause or resume time |
| Reset (R) | Restart from t = 0 with the same seed |
| New seed (N) | A different random world |
| Drag, scroll, pinch | Orbit and zoom |

## Checks

`tests/sim-checks.js` runs the page in headless Chromium and prints measurements: mesh closure, a rigid-rotation test, and long-run statistics.

```
npm install
npx playwright install chromium
npm run check
```

Bit-for-bit reproducibility is not a goal: the simulation is seeded, but results may change with floating-point details and, later, with graphics processing unit (GPU) code paths.

See `PLAN.md` for what the current rules are, what has been measured, and what comes next.
