<!--SPDX-FileCopyrightText: 2025 Bernd Flemisch <bernd@iws.uni-stuttgart.de-->
<!--SPDX-License-Identifier: CC-BY-4.0-->
# Thermodynamics

The two Python scripts located in this folder create CSV tables for the pure phase properties and solubilities.

In addition, a MATLAB/Octave script can be used to generate Eclipse-type tables from the CSV files.

## make_component_table.py

Given a temperature and pressure range, the script queries the [NIST database](https://doi.org/10.18434/T4D303) for the properties of either CO2 or H2O and generates a csv table.
It can be invoked by, for example,
```bash
python3 ./make_component_table.py -c CO2 -t1 10 -t2 20 -nt 3 -p1 1e5 -p2 1.2e5 -np 4
```
for generating a table for CO2 properties in a temperature range [10, 20] °C with 3 sampling points and a pressure range [1e5, 1.2e5] Pa with 4 sampling points.

Currently, the table contains density, viscosity, enthalpy, thermal conductivity, isochoric heat capacity and isobaric heat capacity values.

##  make_solubility_table.py

The script implements the formulas for the solubility from [Spycher et al. 2003](https://doi.org/10.1016/S0016-7037(03)00273-4). It can be invoked similar to the other one by, for example,
```bash
python3 ./make_solubility_table.py -t1 10 -t2 20 -nt 3 -p1 1e5 -p2 1.2e5 -np 4
```

## make_ecl_tables.m

Convert the results from the above into Eclipse-type keywords. Limited to isothermal tables at the moment.