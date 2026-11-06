# Putty Plates

A plate-tectonics simulator for generating Earth-adjacent worlds for a game. The goal is diversity that comes from a few rules acting over a long history, not from layers of procedural noise.

The numerical method is described in `docs/crust_ale_method.tex`: a conservative vertex-centred arbitrary Lagrangian–Eulerian (ALE) scheme for a thin viscous lithospheric sheet. `PLAN.md` maps it onto slices and records what each slice has measured.

The whole app is one file, `index.html`. It loads one script from a content delivery network (CDN), [Delaunator](https://github.com/mapbox/delaunator), and needs a browser with WebGL2. It is written as a page fragment (no `<!doctype>`), because it is published as a claude.ai Artifact. Opened directly from disk it renders in quirks mode, which this layout tolerates.

## Controls

| Control | Effect |
| --- | --- |
| Time step | Million years (Myr) simulated per animation frame; takes effect immediately |
| Vertices | Mesh resolution; regenerates the world when you release the slider |
| Continental area fraction | Initial share of the surface that is continental crust; regenerates |
| Plates | Initial plate count; regenerates |
| Plate strength | Multiplies the lithosphere's depth-integrated viscosity; takes effect immediately. The hint shows the drag length under old ocean |
| Keel drag c_k | How much basal drag grows with lithospheric thickness (crust + thermal lithosphere) |
| Mesh velocity | Markers follow the material (the method's default) or stay put apart from smoothing; seams move either way |
| Smoothing κm | Mesh-velocity relaxation per step |
| View | Elevation, plates and seams, thermal age, crustal thickness, speed, or effective strain rate |
| Show triangulation | Draws the current triangles and seams |
| Show velocity | One segment per marker copy: where it would be after 5 Myr at the current velocity |
| Run / Pause (Space) | Pause or resume time |
| Reset (R) | Restart from t = 0 with the same seed |
| New seed (N) | A different random world |
| Drag, scroll, pinch | Orbit and zoom |

## Checks

`tests/sim-checks.js` runs the page in headless Chromium and prints measurements, numbered as the verification tests in the method document: budgets over whole runs (test 1), free-stream preservation (test 2), rigid rotations being strain-free on the sphere (test 3), the drag-dominated limit u = −∇Φ/β (the instantaneous form of test 4), transform opening against C_n (test 6), the GPE drop against classical ridge push, plate speeds, solver iterations, and how long the mesh lasts before seams would need to change topology.

```
npm install
npx playwright install chromium
npm run check
```

Bit-for-bit reproducibility is not a goal: the simulation is seeded, but results may change with floating-point details and, later, with graphics processing unit (GPU) code paths.

See `PLAN.md` for what the current rules are, what has been measured, and what comes next.
