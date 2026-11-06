# Plan

The algorithm document (LaTeX) has not been received yet. Slice 0 was built without it, so the rules below are placeholders, chosen to get a moving, inspectable planet on screen. Once the document arrives, slices 1 and later get rewritten around it.

## Slice 0: kinematic skeleton (done)

**Mesh.** N jittered Fibonacci points on the unit sphere. The spherical Delaunay triangulation comes from projecting stereographically from one point and triangulating in the plane (circles map to circles, so the Delaunay property carries over). The adjacency is stored in compressed sparse row (CSR) form.

**State on the fixed mesh (Eulerian).** Per vertex: plate id, crust type (continental or oceanic), crust age in Myr, and an "uplift" field in km.

**Plate motion.** Each plate has an Euler (angular-velocity) vector ω and accumulates a pending rotation. Once a plate's pending rotation reaches 1.5 cell spacings, its material is resampled onto the mesh (semi-Lagrangian). Each vertex is traced back through the plate's rotation, and its nearest vertex is found by a greedy walk on the Delaunay graph. Plates below the threshold keep accumulating motion, and the renderer draws them at their pending position, so motion looks continuous.

**Interaction rules (placeholders).**

1. *Overlap* (several plates claim one vertex): the most buoyant crust wins. Continental beats oceanic, and younger ocean beats older.
   - Continental winner over continental loser: collision. The winner gains uplift, and the pair accumulates suture contact.
   - Any other loser: subduction. The winner gains arc uplift.
2. *Gap* (no plate claims a vertex): new oceanic crust with age 0, i.e. a ridge.
3. *Ocean depth* = −2.6 − 0.33·√age km, floored at −6.4 (half-space-cooling-like).
4. *Erosion*: uplift decays with a 200 Myr time constant and diffuses (κ = 900 km²/Myr).
5. *Arc accretion*: oceanic crust with more than 1.5 km of uplift, built up above −2 km, becomes continental. Elevation stays continuous through the conversion.
6. *Suturing*: once two plates' continental contact reaches 1000 km × 250 km of convergence (with a 25 Myr memory), they weld. The merged plate's ω is the area-weighted mean.
7. *Rifting*: each plate rifts as a Poisson process with rate (A_k/A)² · (0.5 + its continental share) / 2 Myr. The cut follows a low-noise flood fill between two distant seeds, and the halves get opposite rotations about the axis between the seeds (2–5 cm/yr opening).

Rule 7's area-squared hazard is what keeps the plate count from collapsing. There is no target count. Welding and rifting balance at a count that depends on the configuration.

### Measured (from `tests/sim-checks.js`, 20,000 vertices, 2 Myr steps)

- The mesh closes: exactly 2N − 4 triangles.
- Rigid rotation (every plate given the same ω): material turns 92.5% (1 cm/yr) and 96.3% (5 cm/yr) as far as it should over 100 Myr. Continental area changes by +0.03% and +2.4%.
- Default rules, 2 seeds, 1500 Myr: the plate count stays between 11 and 25. Continental share drifts from 30% to about 20%. The 95th percentile of continental elevation is about 2 km, with peaks of 4–6 km.

### Known artifacts

- **Near-pole stall.** Nearest-vertex resampling rounds sub-cell motion away. Material within roughly 20° of a plate's Euler pole moves too little, which is the source of the rotation deficit above.
- **Continental area drift.** Continental collisions destroy crust faster than arc accretion makes it. Separately, buoyancy-wins resampling slightly grows continents at boundaries.
- **Ragged boundaries.** Coastlines and plate edges are stair-stepped at the cell scale (about 160 km at 20,000 vertices).
- The CPU does all the stepping, at about 5–9 ms per step for 20,000 vertices.

## Candidate next slices (to be replaced by the document's structure)

1. **Read the document and map it onto the code**: which state variables, which rules, and which time integrator. Write it as a table of document section → function.
2. **Advection with no stall.** Options: keep sub-cell offsets per vertex, use Lagrangian particles with periodic remeshing, or stochastic barycentric sampling. Judge each by the rigid-rotation check (target: ratio ≥ 0.99, area drift < 0.5% per 100 Myr).
3. **Physics for plate velocities** (driving forces, if the document defines them) in place of fixed Euler vectors.
4. **Crust budget.** Thickness as a conserved field rather than a continental/oceanic flag, so that collisions conserve volume.
5. **GPU stepping** (WebGL2 render-to-texture passes, or WebGPU compute) once the update rule is stable enough to be worth porting.
6. **World parameters**: radius, sea level (or ocean volume), and planet-scale outputs (land fraction and hypsometry) reported as emergent statistics.
