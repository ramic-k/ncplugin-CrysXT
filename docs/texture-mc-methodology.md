# Monte-Carlo treatment of textured coherent-elastic scattering

This note documents the methodology behind the Monte-Carlo (MC) scattering support added to the
CrysXT texture model. It is intended for review: it states the model, the sampling procedure, the
exact consistency relations it must satisfy, and the validation that checks them.

## 1. Texture model

CrysXT describes preferred orientation with the **modified March–Dollase** model: a fibre
(axially-symmetric) texture about a preferred-orientation axis `a`, with one or two components of
strength `R_i` and weight `f_i` (`f1 + f2 = 1`). The pole density of a reflection whose plane
normal makes angle `beta` with the axis is

```
P_MD(beta; R) = ( R^2 cos^2(beta) + sin^2(beta)/R ) ^ (-3/2)
```

normalized so its average over the sphere is 1 (`R = 1` is the untextured / random limit).

The powder coherent-elastic cross section is the usual sum over Bragg-active reflections
(`lambda <= 2 d_hkl`):

```
sigma_powder(lambda) = (lambda^2 / 2 V0) * sum_hkl  d_hkl |F_hkl|^2 m_hkl
```

Texture re-weights each reflection by its pole density evaluated at the relevant geometry.

## 2. Monte-Carlo treatment

A textured polycrystal is **anisotropic**, so the correct MC treatment depends on whether a texture
orientation is supplied.

### 2a. Oriented (anisotropic) — the new MC path

When the material is oriented (single-crystal cfg with `dir1`/`dir2`), the texture axes are fixed in
the lab frame. For incident direction `k`:

- **Cross section.** Each open reflection `hkl` contributes `strength_hkl * E_hkl * W_hkl(k)`, where
  `E_hkl` is the optional extinction factor and `W_hkl(k)` is the **cone-averaged pole density** — the
  pole density averaged over the Debye–Scherrer cone of scattering vectors `Q` consistent with `hkl`
  at the Bragg angle. For a single component,

  ```
  W_hkl(k) = < P_MD( angle(Q, a) ) >_psi ,   Q(psi) = sinT*k + cosT*(cos psi e1 + sin psi e2)
  ```

  where `sinT = lambda/2 d_hkl`, `(e1,e2)` span the plane perpendicular to `k`, and the average is over
  the azimuth `psi`. Projecting onto the axis gives `Q.a = sinT*cos(gamma) + cosT*sin(gamma)*cos(psi)`
  with `gamma = angle(k, a)`, so **the cone average depends only on `(gamma, sinT, cos alpha)`**
  (`cos alpha = |normal.a|`). Two components add linearly: `W_hkl = f1 W_hkl^(1) + f2 W_hkl^(2)`.

- **Scattering angle.** A reflection is selected by inverse-CDF with probability proportional to its
  directional partial cross section. The Bragg polar angle is fixed (`cos 2theta = 1 - 2 E_hkl_edge/E`);
  the azimuth `psi` around the Debye cone is drawn from the reflection pole density `P_MD(angle(Q(psi), a))`
  by inverse-CDF (no rejection — guaranteed termination). The measure `d psi` carries no extra
  Lorentz/Jacobian factor.

### 2b. Non-oriented (isotropic)

When no orientation is supplied, the cross section is the powder-with-preferred-orientation
azimuthal average (the original cross-section-only CrysXT value), and scattering is sampled with the
ordinary uniform Debye cone. This equals the oriented value for a beam parallel to the axis (see §3),
and is the correct, fast choice for cross-section and transmission use where the scattering direction
is not tracked.

## 3. Exact consistency relations (anchors)

Two relations follow analytically from the model and are used as validation gates.

**1. Sphere average = random powder.** The cross section is linear in the per-reflection pole weights,
so the average over incident directions is `<sigma(k)>_k = sum_hkl strength_hkl E_hkl <W_hkl(k)>_k`.
Since `P_MD` is normalized to unit sphere-average, `<W_hkl>_k = 1` for every reflection, giving
`<sigma>_k = sigma_powder` (for `E_hkl = 1`). Texture redistributes scattering but conserves the
angle-integrated total.

**2. Beam ‖ axis = the non-oriented value.** For `k = a` the cone vectors satisfy `Q(psi).a = sinT`
for all `psi` (since `e1, e2` are perpendicular to `a`): every cone vector makes the same constant
angle `90°−theta` with the axis. The azimuthal average therefore collapses, `W_hkl(a) =
P_MD(90°−theta)` — exactly the quantity the non-oriented (Sato azimuthal-average) model evaluates.
Hence `sigma(k ‖ a) = sigma_non-oriented`. This is *why* a normal-incidence plate transmission
(beam ‖ axis) correctly uses the fast non-oriented value, while a reflector (all incidence directions;
the scattered direction sets the albedo) requires the oriented path.

## 3a. Implementation in NCrystal

The oriented process is exposed through NCrystal's `ScatterAnisotropicMat` interface
(`crossSection(neutron, dir)` and `sampleScatter(neutron, dir)`); the plugin factory selects it when
the configuration is oriented and the isotropic `ScatterIsotropicMat` otherwise. The Debye-cone azimuth
is sampled by building the discretized cumulative distribution of the selected reflection's pole
density over `psi` and inverting it (histogram inverse-CDF, no rejection — guaranteed termination). For
performance the cross sections are evaluated at runtime from edge-aware tables built once at
construction (validated results-preserving); the exact per-reflection evaluation is retained for
scattering-event sampling and as a reference. The texture model reuses the plugin's existing
March–Dollase machinery; only the direction-dependent evaluation and the angular sampler are new.

## 4. Validation

| check | result |
|---|---|
| Oriented `sigma(k ‖ axis)` vs. non-oriented (Sato) [anchor 2] | ~1e-4 |
| Sphere-average of `sigma(k)` vs. random powder [anchor 1] | ~1e-5 (exact path) |
| **R → 1 limit: textured xs vs. random powder** | **0.0000%** |
| **Isotropic incidence: oriented sampler scattering-angle distribution vs. powder** | **correlation 1.0000, max bin 0.4%** |
| Debye-cone azimuth sampler vs. March–Dollase pole density | correlation 0.997 |
| Extinction-only OpenMC transmission vs. analytic | 0.29% |
| Full-range textured OpenMC transmission vs. analytic | 0.30% |
| VENUS Cu/Ni/Fe textured-plate experimental transmission | 0.4–0.9% at the Bragg edges |

The two bold rows directly answer the natural first question — *does isotropic sampling recover the
powder average?* — at both the cross-section level (R → 1) and the **sampling** level (averaging the
oriented sampler over isotropic incidence reproduces the powder scattering-angle distribution).
Reproduced by `validation/powder_recovery.py`.

## 5. Scope and limitations

- This work provides the **MC-correct sampling of the existing March–Dollase texture model**. It does
  **not** introduce a new texture model.
- March–Dollase is a **fibre (axially-symmetric)** model: the pole density depends only on the angle to
  a single axis (or two, for the two-component form). It cannot represent a general, non-axially-symmetric
  texture (e.g. orthotropic rolled sheet with distinct RD/TD/ND), which requires a full orientation
  distribution function (ODF) / spherical-harmonic pole-density representation. Generalizing the texture
  model to relax axial symmetry is separate, future work.
- The sampling machinery is **model-agnostic**: it samples the Debye-cone azimuth from *whatever*
  per-reflection pole density is supplied. A harmonic/ODF pole density can replace March–Dollase with no
  change to the sampler — this note's framework is the foundation a general-texture model would plug into.
- The texture is an **input** that must be characterized (e.g. by transmission/diffraction); where it is
  unknown or spatially variable (bulk components), the sensitivity can be **bounded** rather than assumed.

## References

- W. A. Dollase, *Correction of intensities for preferred orientation in powder diffractometry:
  application of the March model*, J. Appl. Cryst. **19**, 267–272 (1986).
- X. Cai & T. Kittelmann, *NCrystal: a library for thermal neutron transport*, Comput. Phys. Commun.
  **246**, 107184 (2020).
- Modified March–Dollase texture model and extinction models as implemented in the upstream plugins
  (and references therein): [ncplugin-CrysText](https://github.com/highness-eu/ncplugin-CrysText)
  (texture) and [ncplugin-CrysExtn](https://github.com/XuShuqi7/ncplugin-CrysExtn) (extinction).
