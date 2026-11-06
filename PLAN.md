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

- **Velocity was prescribed** in this slice: u = ω_plate × x, with random Euler poles. Slice 2 replaces it; the Euler poles survive only to classify the initial seams (and as a test hook, `w.kinematic`).
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

## Slice 2: velocity from the force balance, n = 1 (done)

Implemented from the document:

| Document | Code (`index.html`) |
| --- | --- |
| Φ in closed form, Φ̃ = Φ − Φ_r (§1.5, eq. phi) | `phiCol`, `PHI_R` |
| Ĥ = ½ f(s)(B_c h + B_l δ) with n = 1 (§1.6, eq. visc) | `solveVelocity`, per triangle |
| β = β₀[1 + c_k(H/H₀)^m], H = h + δ (§1.7, eq. drag) | `dragBeta` |
| P1 stiffness 2A Ĥ BᵀDB, buoyancy load A Φ̃ ∇φ, lumped drag (§2.5) | `solveVelocity`, `elementFrame` |
| Free ridge faces as the natural condition; extra traction −Φ̃_f n on faces whose facing column is not the ridge column (§2.5) | seam loop in `solveVelocity` |
| Pair coupling C_n, C_t lumped onto copy pairs (§2.5, eq. interface): transform stiff normal, weak shear; trench weak normal, moderate shear | seam loop in `solveVelocity` |
| PCG, warm-started from the previous velocity (§3.4, stage 2, in double precision on the CPU) | `solveVelocity` |

Choices made here:

- **Sphere**: each triangle is treated in its own plane. Each marker carries a tangent basis, and each copy has two unknowns in it. A 2×2 map per vertex takes them into the triangle's plane. A rigid rotation of the sphere is then strain-free to round-off (measured below), so this is the sphere's version of the patch test.
- **Preconditioner**: block Jacobi with one block per marker (all of its copies at once, so the stiff transform penalty between copies is inverted exactly), plus a coarse correction from each plate's three rigid rotations. The coarse space cuts iterations by only about 20%. At these parameters the drag length (about 2,000 km) is shorter than a plate, so the slow modes are below plate scale. That is a job for multigrid (with slice 7).
- **Iteration budget**: the first solve runs to a relative residual of 10⁻⁴. After that each step does at most 40 iterations, warm-started, and convergence carries across steps.
- **Opening trenches**: with no slab pull, some trenches open. The subducting face of an opening trench now accretes new lithosphere, as a ridge face does. Without this, the cells at the trench drained to 0.3 km of crust and moved at more than 100 cm/yr. This stands in for the life cycle's rifting rule (§1.10). The accreted volume counts as created, so the budgets still close.
- **Parameters** (`RHEO`): B_c = 10²⁶ and B_l = 2 × 10²⁵ N·s/m per km; β₀ = 1.5 × 10¹⁴ Pa·s/m, c_k = 3, H₀ = 150 km, m = 3; transform C_n = 10²³ and C_t = 10²⁰ N·s/m²; trench C_n = 10¹⁹ and C_t = 10²¹. With these, old ocean has a drag length √(Ĥ/β) ≈ 2,100 km. The UI exposes a strength multiplier and c_k.

### Measured (8,000 markers, 6 plates)

- **Ridge push from Φ**: the drop in Φ from the ridge to ocean of age 10, 40 and 80 Myr is 0.92, 0.96 and 0.97 of the classical F_rp (eq. ridgepush). Isostatic columns have equal mass above −D to round-off once the thermal term is counted as added density.
- **Continents**: Φ̃ is −0.33 × 10¹² N/m for a 35 km continent with a saturated heat deficit. These continents sit slightly below the ridge column in GPE (gravitational potential energy), so they don't collapse under their own weight; ridges push them.
- **Rigid rotation (test 3 on the sphere)**: max ε̇_E/|ω| = 4 × 10⁻¹⁴.
- **Drag-dominated limit** (negligible strength, one plate, smooth Θ, compared with u = −∇Φ/β): relative RMS error 2.4%, 1.2% and 0.58% at 2,000, 8,000 and 32,000 markers. First order.
- **Transform (test 6)**: RMS opening over RMS slip is 9.3%, 1.9% and 0.9% at C_n × 1, 10 and 100. The floor is geometric: the initial seams zig-zag, and a kinked transform can't slip without opening at the kink. Excluding junctions doesn't change it.
- **Speeds** at 10 Myr (seeds 11, 33, 44): median 2.2–2.6 cm/yr, 95th percentile 5–6 cm/yr, all driven only by ridge push and GPE contrasts.
- **Solver**: a cold solve takes about 220 PCG iterations. With the 40-iteration cap, the velocity is within 0.1% of a tight solve. Converging each step to 10⁻⁴ instead gives 6 × 10⁻⁵, at about 1.6× the cost.
- **Budgets (test 1)**: below 10⁻¹⁴ relative, as in slice 1.
- **Necking**: where a plate is pulled apart, it thins, weakens and thins faster. A Newtonian sheet whose strength scales with its thickness, under a fixed pull, reaches zero thickness in finite time. Crust behind opening trenches falls to 3–5 km and speeds there grow from 7 to about 50 cm/yr over 20–50 Myr. The model is making rifts on its own. Turning them into ridges is the life cycle's job (slice 6). Slab pull (slice 3) would remove the pull that opens these trenches.
- **Mesh lifetime**: seeds 11, 22, 33 and 44 stop at 51, 2.4, 13 and 37 Myr. With prescribed rotations they stop at 22.7, 12.7, 19.8 and 2.1 Myr. Stops are the same all-seam slivers at junctions as in slice 1. Seed 22's "drain" is a near-zero-area copy at such a sliver.
- **Cost**: about 30 ms of assembly and 40 ms of maintenance per step, plus about 0.8 ms per PCG iteration. That is about 150 ms per simulated Myr while speeds are moderate. Necking zones shrink the time step, which pushes it to 400–900 ms per Myr.

## Candidate slices from here

3. **Slab reservoir and power-law creep** (§1.9, eq. visc): W and L per trench copy, Picard iteration. Check: test 7 (Γ < 1 terminal velocity, Γ > 1 runaway).
4. **Second-order transport** (§2.3): least-squares gradients per seam side and the Barth–Jespersen limiter. Check: sharpness of continents after 100 Myr, and the free-stream test.
5. **Exact remap** (§2.9): Sutherland–Hodgman intersection of old and new dual cells, replacing the area-share rule. Check: the free-stream error with maintenance drops to the no-maintenance level.
6. **Seam life cycle** (§1.10, §2.8 "Seam events"): the buoyancy 𝒩 and thermal thickness δ per copy, persistence clocks per edge, then events in this order: slab breakoff (trench → transform), healing (transform → interior), ridge death, subduction initiation, rifting, and transform nucleation from damage. Check: test 8 (budgets across events; continent arriving at a trench; ridge opening in thinned lithosphere). Slab breakoff and healing alone should already remove most of the stops measured in slice 1.
7. **WebGPU kernels** (§3): flux, column-state and element kernels in WGSL with gather instead of scatter and two-sum accumulation. Keep this JavaScript code as the double-precision reference (§3.6).
8. **Sphere-exact geometry**: a geometric conservation law for moving flat triangles on a sphere (for example, the mesh-motion flux taken from swept areas projected on the tangent plane). Check: free-stream error at round-off without maintenance.
