# NCrystal plugin CrysXT

Combined crystalline extinction and texture correction plugin for NCrystal, based on [ncplugin-CrysExtn](https://github.com/XuShuqi7/ncplugin-CrysExtn) and [ncplugin-CrysText](https://github.com/highness-eu/ncplugin-CrysText). Please refer to those repositories for full model descriptions and parameter details. Only the Sabine (uncorr/corr) and Becker & Coppens (BC_pure, BC_mix, BC_mod) extinction models are supported.

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

- `Extinction <model> <l[Å]> <g[1/rad]> <L[Å]> [dist]` — see [ncplugin-CrysExtn](https://github.com/XuShuqi7/ncplugin-CrysExtn) for parameter details. `Sabine_corr` takes no `dist`; `Sabine_uncorr` uses `rect`/`tri`; BC models use `Gauss`/`Lorentz`/`Fresnel`.
- `Texture <px> <py> <pz> <R> <f>` — exactly two lines required, with `f1 + f2 = 1`. See [ncplugin-CrysText](https://github.com/highness-eu/ncplugin-CrysText) for details on the modified March-Dollase model.

## Note

This plugin is intended for cross-section calculations only and not for use in Monte Carlo simulations.
