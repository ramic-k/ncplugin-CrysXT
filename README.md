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

## Note

This plugin is intended for cross-section calculations only and not for use in Monte Carlo simulations.
