# geos-chem-coupling-tools
Tools for coupling GEOS-Chem to the CESM and WRF (WRF-GC)

To use YAML support for GEOS-Chem 14.0.0+, install `pyyaml` using `pip install pyyaml` first.

```
    Copyright (C) 2022 Haipeng Lin

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License version 2,
    as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
```

## generate_species_coupling
This script will generate species-related information lists for CESM and WRF.

**For CESM:** Generates `solsym`, `adv_mass` in `mo_sim_dat.F90`; dry and wet deposition lists for compset xml files; instructions to update `bld/configure`, `chem_mods.F90` for species counts.

**For WRF:** Generates `WRFGC_Get_WRF`, `WRFGC_Set_WRF` two-way species conversion lists and index definitions; `registry.chem` species registry lists; `mozbc` configuration files for initial/boundary conditions.

### Instructions
1. Install and run corresponding version of GEOS-Chem "Classic". Run the model for any amount of time to generate an output species database `.yml` file in `OutputDir`.
2. Move that `.yml` file to project directory, configure options in `generate_species_coupling.py` (file path, target model) then run `python generate_species_coupling.py`.
