My own copy of working with GOTM simulation, forked from [Fox-Kemper Lab](https://gitlab.com/fox-kemper_lab/gotmwork.git).

* To make the tools and paths accessible everywhere regardless of `pwd`:

  * `ln -s <path to gotmwork repo>/set_tools.sh /usr/local/bin/set_tools_gotm`
  * `source set_tools_gotm` will then activate the paths and tools.

A work directory of [GOTM](http://gotm.net), which contains a set of scripts and tools to preprocess the input data, set up runs, and analyze and visualize the output data.

__Latest version__: lt_v1.5, compatible with [GOTM version lt_v1.5.x](https://gitlab.com/fox-kemper_lab/gotm) and [CVMix version lt_v1.3](https://gitlab.com/fox-kemper_lab/CVMix-src).

---
## Quick Start

### Setup
Use `setup_gotmwork.sh` to check Python version (3.x) and required Python module, and set up the necessary environment variables. It only needs to be used once. It will generate a file `.gotmwork_env.sh` in the HOME directory which saves all the necessary environment variables.
You may optionally add `source $HOME/.gotmwork_env.sh` in the file `$HOME/.bashrc` (bash shell) to automatically set up these environment variables when opening a new terminal.

### How to build
CVMix has to be built before building GOTM.See [CVMix Homepage](http://cvmix.github.io) for more information on how to build CVMix.
Use `build_src.sh -build` or simply `build_src.sh` to compile GOTM.
Use `build_src.sh -clean` to clean the old build.
Use `build_src.sh -clean -build` to do a clean build.

### Preprocess data

### Run a case
Change directory to a test case (see [Test Cases](#Test-Cases) for more detail) and run `case_run`.

### Analysis & Visualization

---
## Test Cases

Test cases. In each case, `case_test` sets up the namelist, preprocess the input data and run the simulation, whereas `case_run` (if exists) runs the simulation from preprocessed input data. `case_postproc.sh` contains steps to postprocess output data and is used by `case_run` and `case_test`, either for visulaization or manipulation of the data.

- __COREII__
Run one set of simulations in each 4 by 4 degree box to cover the global ocean, foced by CORE-II.
  - `case_test_multi` sets up multiple runs under CORE-II forcing.
  - `case_run_multi` is similar to `case_test_multi`, but uses preprocessed CORE-II data and is therefore significantly faster.
  - `do_parallel` automatically submit parallel jobs to multiple cores.
  - `kill_all` kills all the jobs.
  - `preproc_data` preprocesses the CORE-II data.

- __JRA55do__
Run one set of simulations in each 4 by 4 degree box to cover the global ocean, foced by JRA55-do.
  - `case_run_multi` sets up multiple runs under JRA55-do forcing using preprocessed data.
  - `do_parallel` automatically submit parallel jobs to multiple cores.
  - `kill_all` kills all the jobs.
  - `preproc_data` preprocesses the JRA55-do data.

- __OCSKEO__

- __OCSPapa__
Single site simulation forced by [Ocean Station Papa](https://www.pmel.noaa.gov/ocs/Papa) data.
  - `case_postproc.sh` is a Bash script to postprocess a single run, used by `case_run`.
- __OSMOSIS__
Single site simulation forced by [OSMOSIS](https://www.bodc.ac.uk/projects/data_management/uk/osmosis/) data.
- __TEST_RES__
Sensitivity test of different boundary layer schemes to different vertical resolutions and time steps.
   - `case_loop.sh` loops over different turbulent methods, vertical resolutions and time steps.
   - `OCSPapa` runs test case using OCS Papa data
   - `OSMOSIS` runs test case using OSMOSIS data
   - `COREII` runs test case using selected COREII data

- __Idealized_Tests__
- __Idealized_Hurricane__

### Preprocessed Data

The preprocessed input data and namelists for [Test Cases](#Test-Cases) are in the directrory `./data/`:

- OCSPapa_20120101-20131204
- OSMOSIS_winter
- OSMOSIS_spring
- COREII_LAT-54_LON254_20080601-20091231
- COREII_LAT10_LON86_20080601-20091231
- COREII_LAT2_LON234_20080601-20091231
- Idealized
- Idealized_Hurricane

In each directory the tool `update_nml` can be used to update the namelist from `./data/namelist/` in the case where new entries are added.

Also included in this directory are the data description files in XML format, which are used by `case_preproc` to preprocess the input data.

### Namelist

All namelists are in the directory `./data/namelist/`. Use `init_namelist` to generate namelist from schemas in the source code.

---
## A List of Tools

A list of tools in the directory `./tools/`.
Most of the tools listed below are written in Python3, some in Bash script. The file `gotmtool.py` contains some shared Python3 functions used by many of the tools. Option `-h` can be used with all tools to get the usage.

| Tool name                  | Description |
| -------------------------- |:----------- |
| `case_preproc`             | Preprocess the input data for GOTM and modify the namelist according to the input xml file. |
| `gotm_archive_data`        | Compress and archive GOTM output data. |
| `gotm_extract_data`        | Extract data from archive generated by `gotm_archive_data`. |
| `gotm_map_quality_control` | Remove runs with NaNs in the output data. |
| `is_sea`                   | Check if the point given by latitude and longitude is a sea point. |
| `init_namelist`            | Initialize namelist from schemas according to the type of turbulence closure. Require Python tool `editscenario`, which can be installed using the script `install_python_tools.sh`. |
| `nc2dat`                   | Convert observational data (OCS etc.) in netCDF format to formatted text file for GOTM input. |
| `nc2dat_argo`              | Convert Argo profile data in netCDF format to formatted text file for GOTM input. |
| `nc2dat_cdip_spec`         | Convert CDIP wave spectrum data in netCDF format to formatted text file for GOTM input. |
| `nc2dat_core2swr`          | Read the daily maximum shortwave radiation from COREII in netCDF format, add an idealized diurnal cycle, and output the hourly shortwave radiation into formatted text file for GOTM input. |
|  `nc2dat_latlon`           | Select CORE-II/JRA55-do/CESM data (in netCDF format) at given latitude and longitude and output into formatted text file for GOTM input. |
| `nc2dat_ww3`               | Convert WW3 wave variables and partitioned surface Stokes drift data in netCDF format to formatted text file for GOTM input. |
| `nmlchange`                | Change the entry value of a namelist. |
| `nmlquery`                 | Query the value of en entry in a namelist. |
| `plotpfl`                  | Plot Hovmoller diagram (time-depth) from GOTM output.
| `plotts`                   | Plot time series from GOTM output. Accept multiple variables. |

### Other Tools

- `gotmanalysis.py` contains classes and functions for data analysis and visualization.
- `windwave.py` contains tools to estimate and test Stokes drift computed from empirical wind wave spectrum.

---
## A List of Scripts

A list of scripts in the directory `./scripts/`.
Bash scripts to setup the tools and runs, and Matlab and NCL scripts to preprocess the input data.

| Script name                | Description |
| -------------------------- |:----------- |
| `argo_mat2nc.m`            | Matlab script to convert Argo profile data from MAT to netCDF. |
| `case_turbmethod.sh`       | Bash script to set the namelist according to the turbulent methods. |
| `cesm_prep_fluxes.ncl`     | NCL script to prepare the surface fluxes data from CESM output. |
| `cesm_prep_profiles.ncl`   | NCL script to prepare the temperature and salinity profiles from CESM output. |
| `cesm_ww3a_to_gx1v6.ncl`   | NCL script to interpolate the WW3 output data onto POP grid gx1v6. |
| `core2_prep_meteo.ncl`     | NCL script to prepare meteorology data from CORE-II. |
| `install_python_tools.sh`  | Bash script to download and install Python tools for GOTM from Github, including `editscenario`, `xmlstore`, `xmlplot` and `gotmgui`. |
| `jar55do_prep_meteo.ncl`   | NCL script to prepare meteorology data from JRA55-do. |
| `ocs_heatflux.ncl`         | NCL script to prepare the net heat flux (excluding shortwave) data from longwave, sensible and latent heat fluxes for GOTM.
| `roms_dz.m`                | Matlab script to generate ROMS style stretching vertical grid. |

### Other Scripts

Scripts in the root directory:

- `build_src.sh` is a Bash script to build GOTM.
- `set_tools.sh` is a Bash script to set up paths and tools, used by `case_run`.
- `setup_gotmwork.sh` is a Bash script to check Python version (3.x) and required Python module, and set up the necessary environment variables.
