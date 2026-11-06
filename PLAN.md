# Plan

The method is `docs/crust_ale_method.tex`: a conservative, vertex-centred arbitrary Lagrangian–Eulerian (ALE) scheme for a thin viscous lithospheric sheet. Markers own dual-cell integrals, triangles carry no state, plate boundaries are seams of split markers, and velocity comes from a depth-integrated force balance. Below, the document is cut into slices that each end in something visible. Section numbers refer to the document.

Two kinds of statement appear here. **Measured** means reproduced by `tests/sim-checks.js`. **Conjecture** or **choice** means not yet tested.

## Slice 0 (superseded)

A rule-based toy (Eulerian nearest-neighbour resampling, scripted rifting and suturing), used only to stand up the renderer and controls. Kept in git history (`53fd7ab`).

## Slice 1: the data model and conservative transport, on a sphere, with prescribed velocity (done)

Implemented from the document:

| Document | Code (`index.html`) |
| --- | --- |
| Column state h, ρc, Θ, s; integrals V, M, Z, S per marker copy (§1.2, §2.2) | `makeWorld`, `w.V/M/Z/S` |
| Thermal closure, exact step update with basal flux q_b (§1.3, eq. coolexact) | `step` → sources |
| Isostasy with water loading, e₀ from the ridge depth (§1.4) | `elevation` |
| Seams as split markers, one copy per angular sector (§2.1) | union-find in `makeWorld`, `deriveCopies` |
| Median dual cells, internal segments, seam face segments (§2.1, Alg. 2) | `computeAreas`, `fluxes` |
| Mesh velocity: material + smoothing; ridge = mean, trench = overriding copy (§2.6) | `meshVelocity` |
| Flux update at midpoint geometry, upwind, first order (φ = 0) (§2.7, Alg. 5) | `fluxes` |
| Ridge inflow, trench outflow, transform exchange (§2.7) | `fluxes` |
| Time-step halving on inversion or over-draining (§2.7) | `step` |
| Budgets V, M, Z tracked against created / subducted / cooling source (§2.2) | `w.budget`, panel |

Substituted for now:

- **Velocity is prescribed**: u = ω_plate × x, with random Euler poles. There is no force balance, Φ, viscosity or slab yet.
- **Sphere**: markers live on the unit sphere. Dual cells use flat triangles, and markers are re-projected after each move. (Document §5, last row.)
- **Mesh adaptation**: local operations stand in for Adapt plus a constrained Delaunay triangulation (CDT) rebuild. Interior markers collapse when an altitude falls below 0.3 cell, interior edges split above 1.6 cells, and non-seam edges flip toward Delaunay. Chain markers merge into a chain neighbour (the document's "merge along the chain"). Seams are never crossed by a flip or split.
- **Remap**: one first-order rule for every operation, applied per seam side. A removed copy hands over everything, losers hand over cell average × area lost, and gainers share the pool by area gained. It is conservative to round-off and keeps a uniform state uniform, but it is not the exact Sutherland–Hodgman intersection remap of §2.9.
- **Deviation**: merging a chain marker moves the seam. The swept sliver's crust is transferred to the side that gains it, at the losing side's cell average. The document's seam refinement never moves the seam, so it has no rule for this.
- **Junctions** (three or more copies) move with the least-squares velocity that best satisfies each incident seam's own rule (a McKenzie–Morgan construction). The document uses the plain average of copies.
- **Initial state**: plates are additively weighted Voronoi regions; continents grow from random seeds to the requested area; seam types come from the initial relative motion (ridge, transform or trench) and stay fixed. The less continental plate subducts. Initial Θ is half-space cooling at an age set by distance from the nearest ridge.

### Measured (8,000 markers, 6 plates, seeds 11/22/33/44)

- **Test 1, budgets**: the worst relative error in ΣV, ΣM and ΣZ over whole runs is below 10⁻¹⁴, across ridges, trenches and every mesh operation.
- **Test 2, free stream** (uniform state, u = 0, smooth imposed mesh motion of up to 5% of a cell per Myr, 4 Myr):
  - Without maintenance: max |Δh|/h = 1.7 × 10⁻⁴. This doesn't shrink with dt, so it is geometric: triangles tilt as markers slide over a curved surface, and the planar geometric conservation law isn't exact there.
  - With maintenance: 1 × 10⁻², from the first-order remap.
- **Mesh lifetime**: runs stop after 2–23 Myr with these seeds (earlier variants reached 30–96 Myr). Every stop is the same kind of event: a triangle whose vertices are all seam markers collapses near one or more triple junctions. Examples are two triple junctions one edge apart, or a four-plate trench junction present from t = 0. These are places where seams need to change topology; slice 1 has no seam events yet.
- **Cost**: about 100–130 ms per simulated Myr on one CPU thread. Mesh maintenance and the rebuilds after it are more than half of that.

## Candidate slices from here

2. **Velocity solve, n = 1** (§1.5–1.7, §2.5, Alg. 4): Φ from eq. (phi), P1 stiffness in a per-triangle tangent frame, lumped drag with the thickness-dependent β of eq. (drag), ridge free faces, seam coupling C_n and C_t. Solve with preconditioned conjugate gradient (PCG) in double precision on the CPU. Check: test 3 (patch test, planar sub-problem), test 4 (porous-medium spreading), drag-dominated limit u ≈ −∇Φ/β.
3. **Slab reservoir and power-law creep** (§1.9, eq. visc): W and L per trench copy, Picard iteration. Check: test 7 (Γ < 1 terminal velocity, Γ > 1 runaway).
4. **Second-order transport** (§2.3): least-squares gradients per seam side and the Barth–Jespersen limiter. Check: sharpness of continents after 100 Myr, and the free-stream test.
5. **Exact remap** (§2.9): Sutherland–Hodgman intersection of old and new dual cells, replacing the area-share rule. Check: the free-stream error with maintenance drops to the no-maintenance level.
6. **Seam life cycle** (§1.10, §2.8 "Seam events"): the buoyancy 𝒩 and thermal thickness δ per copy, persistence clocks per edge, then events in this order: slab breakoff (trench → transform), healing (transform → interior), ridge death, subduction initiation, rifting, and transform nucleation from damage. Check: test 8 (budgets across events; continent arriving at a trench; ridge opening in thinned lithosphere). Slab breakoff and healing alone should already remove most of the stops measured in slice 1.
7. **WebGPU kernels** (§3): flux, column-state and element kernels in WGSL with gather instead of scatter and two-sum accumulation. Keep this JavaScript code as the double-precision reference (§3.6).
8. **Sphere-exact geometry**: a geometric conservation law for moving flat triangles on a sphere (for example, the mesh-motion flux taken from swept areas projected on the tangent plane). Check: free-stream error at round-off without maintenance.
