# DarkAxionPortal

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![DOI](https://zenodo.org/badge/439945365.svg)](https://zenodo.org/badge/latestdoi/439945365)

**Code for studying the cosmology and detection of the Dark Axion Portal Model for Dark Matter.**

The core of the code is in [`code/DarkAxionPortal.py`](code/DarkAxionPortal.py) which includes a class defining the model, as well as functions to calculate couplings, relic abundances, decay rates and other quantities. All calculations are performed in the units defined in [`code/Units.py`](code/Units.py). The remaining scripts are used to generate various plots mapping out the parameter space for cosmological producton and detection of axions and dark photons.

**Authors:** Bradley J Kavanagh, Juan Cortabitarte Gutiérrez

### Getting started

Many of the plotting routines require digitised limits from [github.com/cajohare/AxionLimits](github.com/cajohare/AxionLimits). Download the `AxionLimits` repo and edit the file [`code/Dirs.py`](code/Dirs.py) to point towards your local `AxionLimits/` folder. Then, to generate the plots appearing in [arXiv:2112.XXXXX](https://arxiv.org/abs/2112.XXXXX) (*'Cosmology and direct detection of the Dark Axion Portal'*), simply run the script:
```bash
./GeneratePlots.sh
```

### Acknowledgments

Some code (especially in the `code/PlotFuncs*.py` scripts) has been adapted from [github.com/cajohare/AxionLimits](github.com/cajohare/AxionLimits) (MIT License), which we are grateful for.
