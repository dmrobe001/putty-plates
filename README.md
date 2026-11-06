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
| Slab assimilation τ_a | How long a slab keeps pulling after it sinks; longer means stronger slab pull. The hint gives the steady pull per cm/yr of subduction |
| Mesh velocity | Markers follow the material (the method's default) or stay put apart from smoothing; seams move either way |
| Smoothing κm | Mesh-velocity relaxation per step |
| View | Elevation, plates and seams, thermal age, crustal thickness, speed, or effective strain rate |
| Show triangulation | Draws the current triangles and seams |
| Show velocity | One segment per marker copy: where it would be after 5 Myr at the current velocity |
| Run / Pause (Space) | Pause or resume time |
| Reset (R) | Restart from t = 0 with the same seed |
| New seed (N) | A different random world |
| Drag, scroll, pinch | Orbit and zoom |

## How we work

- **This is for a game, not for science.** Low resolution, loose solver tolerances, poor convergence and departures from Earth's statistics are all acceptable, as long as the simulation runs fairly stably and produces something interesting.
- **Build the physics, not stopgaps.** When the planet misbehaves because a phenomenon isn't implemented yet (slab pull, the seam life cycle, power-law creep), the answer is to implement it. Don't patch in a rule that imitates it, and don't tune parameters to make an incomplete model look normal. Odd whole-planet behaviour is expected until the main phenomena are in place.
- **Check components, not planets.** Test each piece against something it must get right in isolation: an exact solution, a conservation law, a limiting case, a small hand-built configuration. Run a full planet only when a question can only be answered by the whole system, and don't iterate on full runs. Headless browsers emulate the graphics processing unit (GPU) in software, so full runs are much slower there than in a real browser.

## Checks

Two scripts run the page in headless Chromium and print measurements, numbered as the verification tests in the method document.

- `npm run check` (`tests/components.js`, a few seconds): component checks on small hand-built worlds. They cover one-sided seam coupling, slab feedback (test 7) and the slab reservoir surviving mesh maintenance, transform opening against C_n (test 6), rigid rotations being strain-free on the sphere (test 3), the drag-dominated limit u = −∇Φ/β (the instantaneous form of test 4), the drop in gravitational potential energy (GPE) against classical ridge push, budgets (test 1) and free-stream preservation (test 2). Iterate against these.
- `npm run check:planet` (`tests/planet.js`, tens of minutes headless): a whole-planet smoke test of budgets, speeds, solver effort and mesh lifetime for four seeds. Run it occasionally.

`makeWorld` accepts `plateOf`, `seamOf` and `stateOf` to build test worlds by hand, and a world's `bodyForce` hook applies a test load.

```
npm install
npx playwright install chromium
npm run check
npm run check:planet
```

Bit-for-bit reproducibility is not a goal: the simulation is seeded, but results may change with floating-point details and, later, with graphics processing unit (GPU) code paths.

See `PLAN.md` for what the current rules are, what has been measured, and what comes next.
