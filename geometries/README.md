# SPE11 Geometry Files

This folder contains geometry files for the SPE11 variants in the [gmsh](https://gmsh.info/) geometry format.
__Note__: all files have been tested with `gmsh` version 4.11.1. In case you experience issues, consider upgrading
your `gmsh` version. This may be required, in particular, for the Python scripts using the `gmsh` API (see below).


You can open the geometry files with `gmsh`, or let it create a mesh directly via the command line with

```bash
gmsh -2 FILENAME.geo  # creates a 2d mesh in the file FILENAME.msh
gmsh -3 FILENAME.geo  # creates a 3d mesh in the file FILENAME.msh
```

This produces the gmsh file format (`.msh`), however, gmsh can export the mesh into a variety of different formats.
The easiest way to achieve this is by opening the `.msh` file with `gmsh`, and then clicking `file > export`. This
lets you choose from a list of supported file formats. Note that you can also open the geometry in `gmsh` and export
that into different geometry formats, if needed. For instance, in order to use a different meshing tool.

Note that the geometry files contain information about which facies the individual geometries belong to. This information
also ends up in the `.msh` files such that for each element of the discretization, a corresponding facies index can be
retrieved. In `gmsh` this is done via the mechanism of "physical indices". When you export into a different geometry or
mesh file format, you will have to check if or how this information is forwarded.

## SPE11-A (`spe11a.geo`)

This file contains the geometry of the SPE11 variant A. At the beginning of the file you can see a large list of mesh size
definitions around individual points of the geometry. You can modify these individually to achieve the desired mesh,
or, you can conveniently scale all mesh sizes with the provided command-line argument `refinement_factor`. That is, in
order to generate a mesh with a characteristic mesh size twice as coarse as the default, you may invoke gmsh as follows:

```bash
# produces spe11a.msh with mesh sizes twice as coarse as the default
gmsh -2 spe11a.geo -setnumber refinement_factor 2.0
```

## SPE11-B (`spe11b.geo`)

This file contains the geometry of the SPE11 variant B. Internally, it reuses `spe11a.geo` and scales it to the required
dimensions. The default mesh sizes of `spe11a.geo` are also scaled in order to fit to the extended domain sizes. However,
you can still scale the default mesh size with a runtime argument:

```bash
# produces spe11b.msh with mesh sizes twice as coarse as the default for variant B
gmsh -2 spe11b.geo -setnumber refinement_factor 2.0
```

## SPE11-C (`make_spe11c.geo`)

This is a generator script for producing the geometry of the SPE11 variant C. It makes use of the Python API of `gmsh`, and thus,
you need to have `python` and the `gmsh` Python package available. The latter can simply be installed via `pip` with
`pip install gmsh`.

The script takes the `spe11b.geo` file and rotates and extrudes the geometry according to the description of variant C. After
successful execution, the script produces the file `spe11c.geo` (and `spe11c.brep`, which is included internally).
The script also takes a runtime argument `mesh-size` that controls the mesh size to be used at all points of the resulting geometry:

```bash
# generates spe11c.geo with a characteristic mesh size of 250m to be used around all points
python3 make_spe11c_geo.py --mesh-size 250
```


## Structured grid generation (`make_structured_mesh.py`)

If you need a structured mesh for your simulator, you may use this script to generate a `.msh` file containing a structured mesh
including the facies indices that all elements of the mesh belong to. The script takes the SPE variant and the resolution as
runtime arguments. That is, to generate structured meshes for all three variants, you can use the script in the following way:

```bash
# generate a 200x200 mesh for variant A in the file spe11a_structured.msh
python3 make_structured_mesh.py --variant A -nx 200 -ny 200
# generate a 200x200 mesh for variant B in the file spe11b_structured.msh
python3 make_structured_mesh.py --variant B -nx 200 -ny 200
# generate a 200x200x200 mesh for variant C in the file spe11c_structured.msh
python3 make_structured_mesh.py --variant C -nx 200 -ny 200 -nz 200
```

Note that this script also requires the Python API of `gmsh`. Furthermore, passing the flag `--remove-cells-in-seal` creates
mesh files in which the cells in the seal layers are removed.
