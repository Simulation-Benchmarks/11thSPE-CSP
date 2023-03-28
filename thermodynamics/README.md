# Thermodynamics

The two scripts located in this folder create tables for the pure phase properties and solubilities.

## make_component_table.py

Given a temperature and pressure range, the script queries the [NIST database](https://doi.org/10.18434/T4D303) for the properties of either CO2 or H2O and generates a csv table.
It can be invoked by, for example,
```bash
python3 ./make_component_table.py -c CO2 -t1 10 -t2 20 -nt 3 -p1 1e5 -p2 1.2e5 -np 4
```
for generating a table for CO2 properties in a temperature range [10, 20] °C with 3 sampling points and a pressure range [1e5, 1.2e5] Pa with 4 sampling points.

Currently, the table contains density, viscosity and enthalpy values.

##  make_solubility_table.py

The script implements the formulas for the solubility from [Spycher et al. 2003](https://doi.org/10.1016/S0016-7037(03)00273-4). It can be invoked similar to the other one by, for example,
```bash
python3 ./make_solubility_table.py -t1 10 -t2 20 -nt 3 -p1 1e5 -p2 1.2e5 -np 4
```
