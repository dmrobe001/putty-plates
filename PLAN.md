# Plan

The method is `docs/crust_ale_method.tex`: a conservative, vertex-centred arbitrary Lagrangian–Eulerian (ALE) scheme for a thin viscous lithospheric sheet. Markers own dual-cell integrals, triangles carry no state, plate boundaries are seams of split markers, and velocity comes from a depth-integrated force balance. Below, the document is cut into slices that each end in something visible. Section numbers refer to the document.

Two kinds of statement appear here. **Measured** means reproduced by `npm run check` (component checks, `tests/components.js`) or `npm run check:planet` (whole-planet smoke test, `tests/planet.js`); measurements quoted for an earlier slice are as of that slice. **Conjecture** or **choice** means not yet tested.

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
  - With maintenance: 1 × 10⁻² at the time. Most of it was the sliver-transfer bug fixed in slice 2 (see there).
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
| Pair coupling C_n, C_t lumped onto copy pairs (§2.5, eq. interface), one-sided in the normal direction (§1.8): stiff C_n for normal motion the seam cannot absorb (any at a transform, opening at a trench, closing at a ridge), free otherwise | `seamContact`, seam loop in `solvePass` |
| Contact states updated from each iterate (Alg. 4) | `solveVelocity`, up to 4 passes |
| PCG, warm-started from the previous velocity (§3.4, stage 2, in double precision on the CPU) | `solveVelocity` |

Choices made here:

- **Sphere**: each triangle is treated in its own plane. Each marker carries a tangent basis, and each copy has two unknowns in it. A 2×2 map per vertex takes them into the triangle's plane. A rigid rotation of the sphere is then strain-free to round-off (measured below), so this is the sphere's version of the patch test.
- **Preconditioner**: block Jacobi with one block per marker (all of its copies at once, so the stiff transform penalty between copies is inverted exactly), plus a coarse correction from each plate's three rigid rotations. The coarse space cuts iterations by only about 20%. At these parameters the drag length (about 2,000 km) is shorter than a plate, so the slow modes are below plate scale. That is a job for multigrid (with slice 7).
- **Iteration budget**: the first solve runs to a relative residual of 10⁻⁴. After that each step does at most 40 iterations, warm-started, and convergence carries across steps.
- **One-sided normal coupling** (added to the document, §1.8): a trench consumes convergence and a ridge fills opening; motion in the other direction is resisted like a transform. This replaces an earlier stopgap that let opening trenches accrete crust.
- **Parameters** (`RHEO`): B_c = 10²⁶ and B_l = 2 × 10²⁵ N·s/m per km; β₀ = 1.5 × 10¹⁴ Pa·s/m, c_k = 3, H₀ = 150 km, m = 3; transform C_n = 10²³ and C_t = 10²⁰ N·s/m²; trench C_n = 10¹⁹ and C_t = 10²¹. With these, old ocean has a drag length √(Ĥ/β) ≈ 2,100 km. The UI exposes a strength multiplier and c_k.

### Fixes to slice 1 found by the component checks

- **Sliver transfer at chain merges**: when a seam marker merges into its chain neighbour, the seam shifts and a sliver changes sides. The transfer took the sliver's crust only from the merging marker's copy, capped at that copy's content. When the sliver was larger than that cell, the excess stayed on the losing side. A uniform state then came out 24% too thick there and 10% too thin across the seam, at 2,000 markers. Now the losing side hands the sliver over out of its whole remap pool (`remapPatch`, `xfer`).
- **Spherical dual areas**: dual cells use spherical-triangle areas, which tile the sphere exactly. With flat triangles, inserting a marker onto the sphere increases the total area, and the remap, being conservative, thins the new cell by 0.4–1% at 2,000 markers.
- Together, free stream with maintenance (test 2) went from 22% to 0.10% at 2,000 markers, and from 1% to 0.06% at 8,000. That is now the same as without maintenance: 0.08% and 0.02%.

### Measured: components (`npm run check`, about 5 s)

- **One-sided coupling** (two hemispheres split at the equator, Φ̃ = 0, pushed together or apart by a body force; closing rate relative to the same push with no coupling):

  | Seam | Converge | Diverge |
  | --- | --- | --- |
  | trench | 0.88 | 0.002 |
  | ridge | 0.002 | 1.00 |
  | transform | 0.002 | 0.002 |

  Contact states settle in one or two passes. The trench's 0.88 comes from its shear coupling acting on the zig-zag seam, where part of the closing motion is shear along individual edges.
- **Transform (test 6)**, same world under a shear push: RMS opening over RMS slip is 2.1%, 0.77% and 0.35% at C_n × 1, 10 and 100. It falls more slowly than 1/C_n because the seam zig-zags, and a kinked transform can't slip without opening.
- **Rigid rotation (test 3 on the sphere)**: max ε̇_E/|ω| = 3 × 10⁻¹⁴.
- **Drag-dominated limit** (negligible strength, one plate, smooth Θ, against u = −∇Φ/β): relative RMS error 3.5% at 1,000 markers and 1.7% at 4,000. First order.
- **Ridge push from Φ**: the drop in Φ from the ridge to ocean of age 10, 40 and 80 Myr is 0.92, 0.96 and 0.97 of the classical F_rp. Columns are isostatic to round-off.
- **Continents**: Φ̃ is −0.33 × 10¹² N/m for a 35 km continent with a saturated heat deficit, slightly below the ridge column. They don't spread under their own weight.
- **Budgets (test 1)**: 1,500 markers, 10 Myr, every mesh operation: 3 × 10⁻¹⁵.

### Measured: whole planet (`npm run check:planet`), before one-sided coupling

- **Speeds** at 10 Myr: median 2.2–2.6 cm/yr, 95th percentile 5–6 cm/yr, from ridge push and GPE contrasts alone.
- **Solver**: a cold solve takes about 220 PCG iterations. The 40-iteration cap keeps velocity within 0.1% of a tight solve.
- **Necking**: a stretched Newtonian sheet whose strength scales with thickness thins faster where it is already thin. Behind opening trenches crust fell to 3–5 km and speeds rose to about 50 cm/yr over 20–50 Myr.
- **Mesh lifetime**: 51, 2.4, 13 and 37 Myr for seeds 11, 22, 33 and 44, stopped by all-seam slivers at junctions, as in slice 1.
- **Cost**: about 150 ms per simulated Myr at 8,000 markers while speeds are moderate.

## Slice 3: slab reservoir (done)

Implemented from the document (§1.9, eq. slab and slabforce; §2.5 seam terms; Alg. 6):

| Document | Code (`index.html`) |
| --- | --- |
| W and L per trench copy, remapped by length along the trench | `w.Wl` = ∫W dℓ (N), `w.Ll` = ∫L dℓ (m²) per subducting copy; W = ∫W dℓ / ℓ with ℓ from `trenchLengths` |
| dW/dt = gρmα Θ_tr v_s − W/τ_a, dL/dt = v_s, exact decay per step | `slabSources`, fed by the outflow that `fluxes` records per subducting copy |
| Optional input cutoff at L > L_max | `RHEO.Lmax` (0 = off) |
| Slab force ℓ W n and drag ℓ ξ n nᵀ on the subducting copy | seam loop in `solvePass` |

Choices:

- **The reservoir as line integrals owned by copies**: this is the same idea as V = ∫h. A seam merge adds the removed copy's integrals to the copy that absorbs it, so the totals survive maintenance exactly.
- **Parameters**: τ_a = 5 Myr (UI slider, 1–20 Myr), ξ = 2 × 10²² N·s/m² (about 4πη for a mantle viscosity η ≈ 2 × 10²¹ Pa·s). For 80 Myr ocean, τ_a = 5 Myr gives a steady pull of about 12 TN/m per cm/yr of subduction against a slab drag of about 6 TN/m per cm/yr. Slab drag acts on every subducting copy, whether or not it has a slab yet, as in the document.
- **No trench suction yet** (P_ts = 0).

### Measured (`npm run check`)

- **Test 7, slab feedback**: a polar cap of radius 25° made 100× stiffer (so it is nearly rigid, as the test assumes), with a trench on its leading half and a ridge on its trailing half. It is driven by a body force plus its own ridge push, with cooling frozen at Θ(60 Myr). The lumped prediction uses ℓ = Σ ℓ_e(n_e·x̂)²:
  - Γ = 0.47: the speed starts at 2.00 cm/yr (2.22 predicted) and levels off at 4.46 (4.19 predicted), after about 4 τ_a.
  - Γ = 1.41: the speed grows from 2.0 to 14.3 cm/yr over 40 Myr without settling.
- **At normal strength** the same cap is not rigid. Its interior moves near the drag-limited f/β, and only a rim about one drag length wide feels the trench. The rigid-plate Γ doesn't describe it.
- **Reservoir across maintenance**: six seam merges along the trench leave Σ∫W dℓ unchanged.

## Slice 4: power-law creep and damage (done)

Implemented from the document:

| Document | Code (`index.html`) |
| --- | --- |
| Ĥ = ½ f(s) B̄ max(ε̇_E, ε̇_min)^(1/n−1) (eq. visc) | `solvePass`; B̄ calibrated so that B_c, B_l give the n = 1 strength at ε̇_ref = 10⁻¹⁵ s⁻¹; ε̇_min = 10⁻¹⁷ s⁻¹ |
| Picard iteration, with seam contact states updated in the same passes (Alg. 4) | `solveVelocity` |
| G = −k_d ε̇_E s + k_h(1 − s), f(s) = s_min + (1 − s_min)s (§1.8, Alg. 6) | `damageSources`, integrated exactly per step at that step's strain rate; S = V s |
| C_t(s) (eq. interface) | seam loop in `solvePass`: C_t × f(s̄) of the copy pair |

Choices:

- **Passes**: the first solve iterates up to 30 Picard passes (to a 1% velocity change). After that, each step does one pass, so viscosity and contact states lag by one step and keep converging from step to step. On the slab-feedback and localization checks, one pass per step matches two passes to within a few percent.
- **Contact stays stiff**: damage weakens shear coupling C_t, but not the stiff normal C_n that enforces one-sided contact. A damaged transform slides more easily but still doesn't open.
- **Parameters**: n = 3 (UI slider 1–4), k_d = 2 (strain to damage), k_h = 0.05 Myr⁻¹ (20 Myr healing), s_min = 0.1. At these values damage is mild at plate-interior strain rates and strong at boundary rates (s_eq ≈ 0.07 at 10⁻¹⁴ s⁻¹).

### Measured (`npm run check`)

- **Power law**: one plate, negligible drag, a torque-free load (the tangential gradient of a·x). Doubling the load multiplies strain rates by 2.000 at n = 1 and by 7.94 at n = 3 (8 expected). From a cold start, Picard takes 9–12 passes to reach a 1% velocity change.
- **Damage integrator**: six uneven steps at fixed strain rate match the exact solution to 3 × 10⁻¹⁶.
- **Localization**: one plate with no seam, sheared along the equator by a zonal push. The strain rate within 6° of the equator, relative to the global mean, goes 3.68 → 3.87 over 20 Myr with damage, and 3.68 → 3.73 without. Integrity falls to 0.88. So localization happens, but slowly at these k_d and k_h. How strongly it should localize is a design choice that matters for transform nucleation in step 5.
- **One-sided coupling** at n = 3: wrong-direction motion is still 0.005 of the uncoupled rate. A converging trench now passes 0.50 of the uncoupled closing rate (it was 0.88 at n = 1). The response is nonlinear, so these ratios depend on the push.
- **Slab feedback** at n = 3 with a stiff cap: terminal 4.4 cm/yr against 4.19 predicted. At Γ = 1.41 it runs away.

## Slice 5: seam life cycle (done)

The seam is the resolved limit of a damaged band (document §1.8, rewritten): a band one cell wide with damaged strength Ĥ(s) carries the same traction as a seam with C_t = Ĥ(s, ε̇)/w. Damage decides where seams exist; motion and buoyancy decide their type; the transform is the hub between types.

| Document | Code (`index.html`) |
| --- | --- |
| C_t = Ĥ(s, ε̇_w)/w, w = one cell, ε̇_w = slip/w (eq. interface) | seam loop in `solvePass` |
| Seam damage from a fault zone of physical width w_f (default 20 km) | `seamStrainRates`, fed to `damageSources` |
| Seams as an explicit edge set; tips allowed; polarity stored per edge relative to its orientation | `w.seamInfo`, `isSeam`, `deriveCopies` |
| Events with persistence clocks (§1.10, table) | `seamEvents`, parameters in `LIFE` |
| Copies rebuilt by sector, integrals split by corner area (§2.8) | `rebuildCopies` |
| Plates relabelled as the regions seams separate | `relabelPlates` |

Choices:

- **Loading of a transform** is read from its contact force (penalty × normal velocity): compression above P_init with 𝒩 > 𝒩_on on the denser side starts subduction (the denser side subducts); tension above T_rift with δ < 40 km on both sides opens a ridge.
- **Hysteresis**: the trench stall threshold v_stall = 0.05 cm/yr sits below P_init/ξ ≈ 0.08 cm/yr, the convergence at which the same load would restart it. With v_stall = 0.3 a small planet chattered: 60 subduction starts and 32 stalls in 15 Myr. With 0.05 it dropped to 46 and 10.
- **Initial seams** start fully damaged (s = 0).
- **Parameters** (`LIFE`): persistence 2 Myr; s_nuc = 0.15; s_heal = 0.85 with slip below 0.2 cm/yr; 𝒩_on = 1,500 t/m² (ocean of about 20 Myr), 𝒩_off = 0; P_init = T_rift = 5 × 10¹¹ N/m.

### Measured (`npm run check`, about 15 s)

- **Test 8, events** (two hemispheres split at the equator unless noted; cooling frozen; budgets close to about 5 × 10⁻¹⁵ in every scenario):
  - **Healing**: an idle transform (k_h raised to keep the test short) welds shut: 59 of 61 edges heal and the rest go by mesh merges. The result is one plate, no split markers, and h uniform to 0.13%.
  - **Nucleation**: a seamless sheet sheared along the equator, with fast damage (k_d = 40), grows 501 transform edges, all within 10° of the equator. It is a wide swath, not a single fault: every damaged edge in the band nucleates.
  - **Subduction initiation**: a transform pushed shut between 80 Myr and 5 Myr ocean becomes a trench on 57 of 61 edges. All 57 subduct the old side.
  - **Rifting**: pulled apart, 2 Myr lithosphere opens a ridge (59 of 60 edges) and creates crust; 80 Myr lithosphere stays a transform.
  - **Ridge death**: a ridge pushed shut becomes a transform.
  - **Breakoff** (stiff cap at 4,000 markers, continent riding to the centre of its trench, steady external push): the central trench subducts until the continent arrives (about 6 Myr), then breaks off to a transform. With the push still on and 60 Myr ocean on the other side, subduction restarts with the opposite polarity by 9 Myr. That is polarity reversal after collision, from the rules alone.
- **Slab feedback (test 7)** now needs a frictionless megathrust to match the rigid-plate prediction. The check uses a near-zero fault width, so seams stay at s ≈ 0. At Γ = 0.47 it levels off at 3.6–3.9 cm/yr (4.19 predicted); at Γ = 1.41 it runs away. With 100× plate strength and the default fault width, the megathrust's own C_t = Ĥ(s)/w welds the cap.

### Measured: small planets (2,000 markers, seeds 11 and 22; run once)

- Events fire throughout: in seed 11's first 6 Myr there are 18 nucleations, 38 subduction starts and 34 breakoffs (many at continental margins, followed by reversed subduction), plus rifts.
- Both runs still stop early, at 5.8 and 16.3 Myr, with the mesh tangling at all-seam triangles at junctions (markers with 3–4 copies).
- **The conjecture from slice 1 is false**: events don't remove junction slivers. Frequent type changes add junction complexity.

## Slice 6: the mesh update (done)

Document section "Mesh update" (§2.10). Seam segments are triangle edges, so "the mesh never tangles" and "seams never cross" are the same statement. The update is a move that never inverts a triangle (Δt halving, as before), then a repair between steps that restores an invariant: every altitude ≥ 0.3 L0, every edge < 1.6 L0. Repair never moves a marker and never takes one off a seam. Each violation is classified by its thin direction (vertex v too close to edge ab) and by what v and ab are, which picks one local, conservative operation.

| Document | Code (`index.html`) |
| --- | --- |
| Long seam edge → split, one copy per side, slab reservoir split by length | `splitSeam` |
| Seam marker near a seam edge, not near its ends → snap: flip ab, reroute the seam a–v–b | `snapSeam` (falls back to `splitSeam` at the foot) |
| Seam marker near another seam marker → zip: collapse, seams moved and merged | `collapseAny`, `mergeSeams` |
| Candidate search, priorities, up to three rounds per step | `seamRepair`, called from `maintain` before the flips |
| Remap when area changes sides: net area moves between the groups (sides) of a patch at the donor's average; vanished copies hand over everything | `remapCavity` |

Choices:

- **Zip and snap are seam events.** All the other operations are numerics. Two seams closer than a cell are one damaged band, and a band one cell wide is one seam. That is the resolution limit of the slice-5 equivalence. Merge rule for coinciding seam edges: a trench wins and keeps its polarity, unless the other is a trench of the opposite polarity; then the two outer plates collide and the edge becomes a transform for the life-cycle rules to reclassify. Equal types stay; any other pair becomes a transform.
- **Junctions keep their markers.** A collapse keeps the endpoint with more seam edges. Junction migration and plate consumption then follow from chain merges and zips, with no special junction code.

### Measured (`npm run check`, about 20 s)

- **Test 10, plate consumption** (prescribed rotations about the pole, 2,000 markers): a 30° wedge subducting under both neighbours at 5 cm/yr each pinches off at the poles first, then zips shut along its length.
  - It is gone at 35 Myr, against 33.4 predicted. The last slivers go once they are thinner than 0.3 L0.
  - Its two trenches become 38 transform edges between its former neighbours (opposite polarities → collision).
  - No trench ever had the wrong side subducting; budgets close to 4 × 10⁻¹⁵; the run continues to 45 Myr.
- **Repair alone keeps a uniform state uniform** to 6 × 10⁻¹⁵ over the whole consumption (state reset to uniform before each maintenance, measured after it), through 24 zips, 16 snaps and 14 seam splits.
- **Seam splits** conserve V, Z and the slab reservoir to about 1 × 10⁻¹⁵ and keep every trench's polarity.
- **The event scenarios used to stop on tangles, and the checks didn't say so.** Nucleation stopped at 14.1 Myr and subduction initiation at 2.6 Myr. Their slice-5 numbers were measured just before the runs died. Every scenario now reports `ran_to_Myr` and `stopped`. The runs no longer stop, so each now runs only as long as its rule needs: nucleation to 14 Myr, initiation to 3 Myr.
  - Results in those windows: nucleation 499 edges, all within 10° of the equator; initiation 61 of 61 trench edges subducting the old side.

### Found here (slice 7 below fixes the first)

- **Transport at junctions.** Copies at and next to a junction drifted away from a uniform state by about 1.5% per Myr. The old code shows the same drift to the third digit.
- **Beyond the windows:**
  - Under the nucleation push, the swath zips into fewer seams after 14 Myr (499 → 258 edges by 16 Myr). Then new seams nucleate outside the band: 494 edges, only 228 near the equator, by 18 Myr.
  - Under the unrelenting initiation push, interior damage spreads and nucleation runs away after 5 Myr.

## Slice 7 (step 6): seam-face fluxes consistent with how faces move (done)

**Cause, measured.** The old flux declared some seam faces flux-free: a trench's overriding face, and a ridge face the plate moves toward. Faces move with the marker velocity w, and at a junction w cannot equal every overriding plate's velocity. Those faces therefore swept area without material.
- In a uniform state, each drifting copy's content-versus-area mismatch over one step equalled the area its flux-free faces swept, to two digits (e.g. −1.06% against −1.05%).
- Only 19 copies drifted, all at seams next to junctions.

**Fix** (document, "Seam faces"; `fluxes`). Every face carries the flux of its own (u − w)·n. Facing segments first exchange: what one side sends fills the other's gap, upwind. The remainders sum to the physical convergence (u₁ − u₂)·n, independent of w. Only that remainder depends on the seam type:
- trench: excess is subducted (only the subducting side feeds the slab);
- ridge: an unfilled gap is new lithosphere;
- elsewhere: thickening or thinning.

Along ordinary seams this reduces to the old one-sided rules.

### Measured (`npm run check`, about 16 s)

- **Free stream, with seam and junction markers moved too** (new variant):

  | | Old flux | New flux |
  | --- | --- | --- |
  | Worst cell | 153% | 1.0% |

  The remaining errors are in sliver copies of 0.7–4.5% of a typical cell's area. They do not grow linearly with the imposed motion: 0.65% at a quarter of the amplitude, and exactly zero with no motion. One possibility is the known flat-versus-spherical geometry floor (0.1% in ordinary cells) divided by a tiny area. That doesn't explain the shape of the amplitude dependence, which I haven't explained.
- **Free stream, interior motion only, no maintenance:** 7.2% → 0.10%. The old figure was mostly this defect too.
- **Test 10 wedge, before it is gone:**
  - interior copies 3.8% → 1.7%;
  - copies anywhere 230% → 140%. What remains is real convergence that the prescribed rotations force across ridge edges at the polar junctions (both faces close). Ridges have no sink, and in the dynamic model the one-sided coupling would stiffen a closing ridge.
- **Everything else is unchanged:** one-sided coupling, slab feedback, every event scenario, and budgets at about 4 × 10⁻¹⁵.
- **The runaway nucleation under the unrelenting initiation push is not caused by the junction drift.** It is essentially unchanged: 695 nucleations by 8 Myr, against 739 before the fix.

## Slice 8 (step 7): nucleation along creases of integrity (done)

A fault forms along the weakest line of a damaged band, not across its width. The old rule nucleated every edge whose two ends averaged below s_nuc, which filled a band several cells wide with seams. The new rule (document §1.10 table; `seamEvents`) is a crease test:
- both ends of the edge are below s_nuc;
- both markers flanking it (the apexes of its two triangles) are less damaged than either end.

Edges across a band fail, because one end is up the flank. Edges beside a seam fail, because the seam's markers (integrity ≈ 0, a marker's integrity being that of its most damaged copy) are lower. The persistence time is unchanged.

### Measured (`npm run check`, about 15 s)

- **Shear-band nucleation** (one sheet, zonal shear, k_d = 40):

  | | Old rule | Crease rule |
  | --- | --- | --- |
  | Seam edges by 14 Myr | 501 (the whole band) | 76 |
  | Chains | | 2 |
  | Branch markers | | 2 |
  | Farthest from the equator | 10° | 2.6° |

  A full ring at this resolution is about 62 edges.
- **The post-window runaways are gone:**
  - Shear band, extended to 20 Myr: 74–84 edges throughout, all within 10° (old rule: 1,512 edges, two-thirds outside 10°).
  - Unrelenting initiation push, to 8 Myr: 3 nucleations (old rule: 695). Interior damage still goes deep (minimum integrity 0.008), but forms no creases.
  - So the runaway came from the old rule, not from the junction drift or the push.
- **Everything else is unchanged:** subduction initiation 59/59 old side, rifting, ridge death, breakoff and reversal, budgets about 4 × 10⁻¹⁵.

### Noticed, not investigated

- **Shear-only push:** 37 transform edges become trenches by 20 Myr. Possibly local compression on zig-zag edges.
- **Initiation push:** seam edges grow from 61 to 132 by 8 Myr without nucleation. Possibly seam splits along the trench.

## Slice 9: responsive page, and chasing the spike (in progress)

### The page

- **Responsiveness.** The simulation runs in a background worker built from the page's own script, and the page only draws. The page asks for the next update once the last has arrived, so a slow step never freezes drawing or input. A banner shows how long the current update has taken. If the browser refuses the worker, it runs in the page one step per task.
  - Measured headless (software graphics): drawing stays at about 23 fps while the simulation computes.
  - The simulation itself ran faster once off the drawing thread: 14 Myr in 10 s at 7,900 markers.
- **Garbled mesh on a new seed.** The renderer re-uploads triangles when a topology counter changes, and each world started its counter afresh, so a new world could match the old number. The counter is now global.
- **Also added:**
  - a seed field;
  - a floating banner for stop reasons and errors;
  - "Extremes" (highest point and thickest crust, with what is there, and a button to face it);
  - on phones, the globe is centred above the open controls.

### The spike: a velocity runaway, with crust piling up downstream

**Reproduced** headless at the default settings (seed 1, 7,900 markers): speeds of 3–5 cm/yr jump to hundreds or thousands of cm/yr within about a Myr, then subside. Crust piles up to hundreds of km where the runaway converges.

**Defects found and fixed along the way:**

| Fix | Blowup time (seed 1) |
| --- | --- |
| None | 15 Myr |
| The overriding plate is no longer subducted | 22.8 Myr |
| Polarity runs and tears | 32.6 Myr |
| Leaky transforms, decompression melting, slab-pull cap | 32.6 Myr |

1. **The overriding plate was subducted at its own trench.** The step-6 rule put convergent excess from either face of a trench down the trench. Only the subducting face's excess goes down now; overriding excess stays.
2. **Trench polarity flipped at single markers.** One copy was overriding on one trench edge and subducting on the next, pulled under from both sides.
   - Initiation now gives one polarity to each run of edges starting together. A run inherits an adjacent trench's polarity if that pushes down dense enough lithosphere.
   - Remaining flips tear: the less dense side's edge becomes a transform and releases its slab. A tear changes one edge, so it can't cascade.
   - Breakoff then polarity reversal still works, now with tears at the segment ends.
3. **Gaps at transforms stayed unfilled,** thinning cells at ridge–transform corners by about 2% per Myr. Any unfilled gap is now new ocean lithosphere (a leaky transform).
4. **Decompression melting.** Thin lithosphere that stretches adds oceanic crust; this is the distributed form of the ridge rule.
5. **Slab-pull saturation.** W ≤ g ρm α Θ_ss L_max (49 TN/m at 660 km), replacing the old cumulative cutoff, and now on by default. In the slab-feedback check, the gain-1.41 case levels off near 9 cm/yr instead of running away.

**Still open.** In seed 1 at 32.6 Myr the runaway starts with slab pull at only 21.5 TN/m, so the cap doesn't stop it.
- **Where it starts:** in overriding copies at a ridge–trench junction that are fully damaged (integrity 0) and young.
- **Delayed but not prevented** by a higher floor on damaged strength (s_min 0.3: 36.5 Myr).
- **Reading, not proven:** damaged cells decouple from their plate. Edge loads (the free-face potential across the trench) then act on a single cell resisted only by its own basal drag.
- **Unexplored option:** treat a sector bounded by seams and fully damaged lithosphere as below resolution, like a zip.

## Next, in order

Physics first (see "How we work" in the README). Each step comes with component checks on small hand-built worlds (`makeWorld` with `plateOf`, `seamOf`, `stateOf`; `w.bodyForce`).

Next: settle the open runaway above with you. It touches modelling choices: how weak fully damaged lithosphere is, and whether grid-scale microplates should exist.

Numerics, only when stability or looks need them: second-order transport (§2.3), the exact intersection remap (§2.9), a geometric conservation law exact on the sphere, WebGPU kernels and multigrid (§3).
