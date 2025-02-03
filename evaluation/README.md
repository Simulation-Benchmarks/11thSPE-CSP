<!--SPDX-FileCopyrightText: 2025 Bernd Flemisch <bernd@iws.uni-stuttgart.de-->
<!--SPDX-License-Identifier: CC-BY-4.0-->
# Evaluation

This folder contains scripts for evaluating results submitted to the 11th SPE CSP.

## Checking the format

Given a folder and case letter, the script [check_format.py](https://github.com/Simulation-Benchmarks/11thSPE-CSP/blob/main/evaluation/check_format.py) checks if the result files present in the folder are conforming to the requirements as given in the [description](https://doi.org/10.2118/218015-PA). For example,
```bash
python3 ./check_format.py -f path_to_folder -c A
```
will check the files in `path_to_folder` for the case CSP 11A. If an archive `spe11a.zip` is found, the files in that archive will be checked.

## Visualizing results

For each case, scripts are provided for visualizing the sparse data / time series and dense data / spatial maps as well as corresponding performance data. All scripts offer the same mechanism for specifying groups and folders via the following options:

* `-g`, `--groups`: Specify group names. For example, `-g GroupX GroupY-1` will specify two groups with the indicated names.

* `-gf`, `--groupfolders`: Specify for each group a folder where the result file(s) to be evaluated have to be located. Amending the example with `-gf path_to_folder1 path_to_folder2` specifies that `GroupX`'s result files are located in `path_to_folder1` and `GroupY-2`'s result files are located in `path_to_folder2`. Relative and absolute paths may be used. The paths should not contain any filename.

* `-f`, `--folder`: As an alternative to `-gf`, one single folder can be specified which contains the results of the groups as subfolders. For example, if `-f path_to_folder` is used instead of `-gf ...`, result files of `GroupX` are expected in `path_to_folder/groupx`, and result files of `GroupY-1` are expected in `path_to_folder/groupy/result1`.

### Sparse data / time series

The scripts `spe11x_assemble_time_series.py` visualize sparse data. They assume that each of the specified folders contains a file `spe11x_time_series.csv` in the required format. Image files `spe11b_time_series_{boxA, boxB, boxC, pressure, seal}.png` are created.

Analogously, the scripts `spe11x_assemble_performance_time_series.py` visualize performance sparse data.

### Dense data / spatial maps

The scripts `spe11x_visualize_spatial_maps.py` visualize dense data. In addition to the options discussed above, a specific reporting time step has to be specified with the option `-t`, `--time`. For example, `-t 100` will analyze the files `spe11a_spatial_map_100h.csv` in case of Case A and `spe11b_spatial_map_100y.csv` in case of Case B. If only one group name `GroupX` is provided via `-g`, the spatial maps of all relevant quantities are visualized together in one file, e.g., `spe11a_groupx_100h.csv`. If more than one group is specified, one file per quantity if produced, amounting to, e.g., files `spe11b_{pressure, saturation, mco2, mh2o, rhog, rhol, tmco2, temp}_100y.png`.

Analogously, the scripts `spe11x_visualize_performance_spatial_maps.py` visualize performance dense data.

## Calculating convection from spatial maps

The scripts `spe11x_calculate_convection.py` calculate the convection integral ((17) in the description) based on given spatial maps. They can be used to double-check the results required as part of the sparse data. They also take options `-g`, `-gf` and `-f` as discussed above.

## Calculating Wasserstein distances

The script `emd.py` calculates the Wasserstein distance between two CO2 mass distributions given via spatial maps. The following parameters are to be passed:

* `-in1`, `-in2`: Two spatial map files.

* `-nx`, `-ny`: Number of reporting cells in x and y direction. Defaults to `280` and `120`, corresponding to Case A. For Case B, `840` and `120` should be used.

* `-sf`: An integer shrinking factor `>= 2` needs to be provided since the underlying algorithm consumes too much memory. Defaults to `2`, corresponding to Case A. For Case B, at least `3` should be used.

The script `spe11a_calculate_wasserstein_distances.py` calculates distances for arbitrary number of groups and reporting time steps. It takes the same options `-g`, `-gf` and `-f` as discussed above, plus an option `--hours` by which the time steps can be passed. The result is a distance matrix stored in a file `distances.csv`.
