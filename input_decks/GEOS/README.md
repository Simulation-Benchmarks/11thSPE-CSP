<!--SPDX-FileCopyrightText: 2025 Jacques Franc <jacquesfrancdev@gmail.com-->
<!--SPDX-License-Identifier: CC-BY-4.0-->
## GEOS input deck

The *SPE11b* input deck for GEOS open-source simulator can be found in the repository under [/inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b](https://github.com/GEOS-DEV/GEOS/tree/develop/inputFiles/compositionalMultiphaseFlow/benchmarks/SPE11/b).

It consists in:

1. Two deck files, one for the specific discretization, *spe11b_vti_source_00840x00120.xml*, one with the base XML tags shared by all cases *spe11b_vti_source_base.xml*.
2. an *include* folder that gathered all sub-XML related to boundary conditions, relperms per facies and other per facies properties
3. a *mesh* folder with the mesh file under VTI vtk's format
4. a *table* folder with all CSV-like data read from to fill in relative permeabilities and capillary pressures

A full tutorial can be found in [the GEOS documentation](https://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/advancedExamples/validationStudies/carbonStorage/spe11b/Example.html)
