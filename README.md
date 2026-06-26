# NCrystal plugin CrysXT

Combined crystalline extinction and texture correction plugin for NCrystal, based on [ncplugin-CrysExtn](https://github.com/XuShuqi7/ncplugin-CrysExtn) and [ncplugin-CrysText](https://github.com/highness-eu/ncplugin-CrysText). Please refer to those repositories for full model descriptions and parameter details. Only the Sabine (uncorr/corr) and Becker & Coppens (BC_pure, BC_mix, BC_mod) extinction models are supported.

BC models support an optional `rec=` flag selecting the y(x,θ) recipe: `cls` (BC1974 classic formula, Becker, P. J. & Coppens, P. (1974). Acta Cryst. A30, 129-147), `std` (BC2025 standard precision, default), or `lux` (BC2025 luxury precision). Refs: T. Kittelmann et al., Acta Cryst. (2026). A82, 163-178, https://doi.org/10.1107/S2053273326001245; Thomas Kittelmann. (2026). tkittel/bc-extinction-paper: Supporting material for BC extinction recipe updates v1.0.1 (v1.0.1). Zenodo. https://doi.org/10.5281/zenodo.18493059

## Installation

Create a conda environment, install NCrystal via `conda install ncrystal`, then:

```
pip install "git+https://github.com/dddijulio/ncplugin-CrysXT"
```

## Syntax

Append a `@CUSTOM_CRYSXT` section to any `.ncmat` file using the `Extinction` and/or `Texture` keywords. Both are optional but at least one must be present.

```
@CUSTOM_CRYSXT
  Extinction Sabine_uncorr 30000 1000 100000 rect
  Texture 1 1 1 0.5 0.5
  Texture 1 1 1 0.5 0.5
```

```
@CUSTOM_CRYSXT
  Extinction BC_mix 30000 0.001 100000 Gauss rec=lux
```

- `Extinction <model> <l[Å]> <g[1/rad]> <L[Å]> [dist] [rec=cls|std|lux]` — see [ncplugin-CrysExtn](https://github.com/XuShuqi7/ncplugin-CrysExtn) for parameter details. `Sabine_corr` takes no `dist`; `Sabine_uncorr` uses `rect`/`tri`; BC models use `Gauss`/`Lorentz`/`Fresnel`. The optional `rec=` flag (BC models only, default `std`) selects the extinction recipe: `cls` = BC1974 classic formula, `std` = BC2025 standard precision, `lux` = BC2025 luxury precision.
- `Texture <px> <py> <pz> <R> <f>` — exactly two lines required, with `f1 + f2 = 1`. See [ncplugin-CrysText](https://github.com/highness-eu/ncplugin-CrysText) for details on the modified March-Dollase model.

## Monte Carlo use (texture)

The texture model now supports **Monte-Carlo scattering** when the material is used as
an **oriented** single crystal (i.e. an orientation is supplied in the cfg via
`dir1`/`dir2`, as for any oriented NCrystal crystal). In that case the textured
coherent-elastic process is **anisotropic**: the cross section depends on the incident
neutron direction relative to the lab-fixed texture axis, and the scattering azimuth
around each Debye–Scherrer cone is sampled from the modified March-Dollase pole density.
This is the physically correct treatment for a textured (fibre/plate) polycrystal.

Example (texture axis [111] placed along the +x lab axis):

```
Fe_sg229_CrysXT.ncmat;dcutoff=0.5;mos=0.0005deg;dir1=@crys_hkl:1,1,1@lab:1,0,0;dir2=@crys_hkl:1,-1,0@lab:0,1,0
```

If the material is **not** oriented, the texture cross section is the orientation-
averaged value (equivalent to a beam parallel to the texture axis) and the scattering is
sampled isotropically — appropriate for cross-section/transmission use but not for a
direction-resolved MC.

**Extinction** needs no special MC treatment. It is a per-reflection *scalar* magnitude
reduction (`E_hkl ≤ 1`, isotropic, independent of sample orientation): it lowers the
strength of each Bragg edge but does not alter the scattering-angle distribution, so the
effective cross section the Monte-Carlo code already transports is correct. The plugin
weights hkl selection by `strength·E_hkl`, consistent with the extinction-reduced cross
section. This was verified by an extinction-only uncollided transmission in OpenMC
(MPI, 1e8 neutrons): the MC-extracted Σ matches the analytic extinction-reduced cross
section to ~0.3%, and at the strongly-extinguished (110) edge tracks the reduced curve
(not the kinematic peak) to <1%.

The texture sampler has been validated to reproduce: the orientation-averaged cross
section equal to the powder cross section (texture conserves the angle-integrated total);
the beam-parallel-to-axis cross section equal to the previous (cross-section-only) value;
and the sampled Debye-cone azimuth distribution equal to the March-Dollase pole density.

## Performance

Both the non-oriented and the oriented textured cross sections are tabulated at construction, so a
`crossSection` query is an interpolation rather than a per-reflection Bragg sum.

The **non-oriented** (isotropic) cross section uses an edge-aware energy table of `g(E)=xs·E` (the
energy-integrated cross section is preserved to <0.001%; scattering sampling still uses the exact
per-reflection weights).

The **oriented** (single-crystal) cross section depends only on the angle to each texture axis, so it
is tabulated per component as `σ_i(cos γ_i, E)`. The Debye-cone-averaged pole density — itself a
function of `(cos γ, sin θ, cos α)` — is precomputed in a small 3-D table that also accelerates the
sampler's plane selection. The tabulated oriented cross section reproduces the exact per-reflection
value to ~0.1% (the linear pole-density tables assume a moderate March-Dollase `R`, as in typical
fibre/plate textures).
