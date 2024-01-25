## Eclipse input deck

This input deck for Eclipse consists of a simple grid together with porosity and intrinsic permeability.
Facies identifiers are associated with grid cells via a region array.

The grid data can be used to construct "black oil" (with PVT tables) or compositional (EOS) simulation cases.
For more information on the data, please read the comments at the top of the text files.

For obtaining PVT tables which are consistent with the CSP description, we refer to the script
[make_ecl_tables.m](https://github.com/Simulation-Benchmarks/11thSPE-CSP/blob/main/thermodynamics/make_ecl_tables.m).

A repository containing complete input files with some simplifications can be found in [the SPE11 decks repository](https://github.com/sintefmath/spe11-decks).
