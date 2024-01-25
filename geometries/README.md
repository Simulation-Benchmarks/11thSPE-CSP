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

__Important__: the geometry is defined in the x-y plane, that is, all point coordinates are given in the form of `(x, y, 0.0)`.
The description defines the geometry in the x-z plane (to be in agreement with variant C), however, we have chosen the x-y plane
here in order to facilitate the parsing of two-dimensional meshes created with this geometry file.


This file contains the geometry of the SPE11 variant A. At the beginning of the file you can see a large list of mesh size
definitions around individual points of the geometry. You can modify these individually to achieve the desired mesh,
or, you can conveniently scale all mesh sizes with the provided command-line argument `refinement_factor`. That is, in
order to generate a mesh with a characteristic mesh size twice as coarse as the default, you may invoke gmsh as follows:

```bash
# produces spe11a.msh with mesh sizes twice as coarse as the default
gmsh -2 spe11a.geo -setnumber refinement_factor 2.0
```

Besides this, you can use the `with_facies_7` option (default: `1`), to remove all elements in facies 7:

```bash
# produces spe11a.msh with all elements in facies 7 removed
gmsh -2 spe11a.geo -setnumber with_facies_7 0
```

## SPE11-B (`spe11b.geo`)

__Important__: as for variant A, the geometry is defined in the x-y plane, while the description defines the geometry in the
x-z plane in order to be in agreement with variant C. See the above section for more details.

This file contains the geometry of the SPE11 variant B. Internally, it reuses `spe11a.geo` and scales it to the required
dimensions. The default mesh sizes of `spe11a.geo` are also scaled in order to fit to the extended domain sizes. However,
you can still scale the default mesh size with a runtime argument:

```bash
# produces spe11b.msh with mesh sizes twice as coarse as the default for variant B
gmsh -2 spe11b.geo -setnumber refinement_factor 2.0
```

You can also again decide to remove all elements inside facies 7:

```bash
# produces spe11b.msh with all elements in facies 7 removed
gmsh -2 spe11b.geo -setnumber with_facies_7 0
```

## SPE11-C (`make_spe11c.py`)

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

__Important__: for the variants A & B, keep in mind that the geometries are defined in the x-y plane instead of the x-z plane
used in the description. See the above sections for more details.

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

Note that this script also requires the Python API of `gmsh`. Furthermore, passing the flag `--remove-cells-in-facies-7` creates
mesh files in which the cells in all regions of facies 7 are removed.


## Extrusion of 2D meshes to one cell thick 3D meshes

The script `extrude_and_rotate.py` can be used to extrude 2D generated meshes to one cell thick 3D vtk
versions for simulators that are not accepting such 2D meshes, .

Leveraging [vtk](https://pypi.org/project/vtk/), it reads from a [gmsh](https://gmsh.info/) generated `vtk` mesh, extrudes it
with one cell in the 3rd dimension and rotates it in the x-z plane.

__Important__: The script is not detecting either if it is a structured mesh (based on *quads*) or unstructured
(based on *triangles*), so we have to provide it with the explicit option `--quad` or `--tri`.

The full procedure is then, for instance for spe11-a:

```bash
gmsh -2 spe11a.geo -format vtk
python3 extrude_and_rotate.py --tri --spe a spe11a.vtk
#optionally some clean up
rm -iv spe11a.vtk
```

This should produce a file `spe11a_extruded.vtu` that can be inspected with [paraview](https://www.paraview.org/) or another
vtk enabled 3D reader. The cell data of porosity and permeability are tagged as _PORO_ and _PERM_ in this output. The facies labels
are registered under _attribute_.

An additional option exists for applying pore-volume modification for case spe11-b. Then the altered procedure, starting this time from
a structured mesh is:

```bash
python3 make_structured_mesh.py --variant B -nx 300 -ny 100
meshio convert spe11b_structured.msh spe11b_structured.vtu
python3 extrude_and_rotate.py --quad --poromult --spe b spe11b_structured.vtu
#optionally some clean up
rm -iv spe11b_structured.vtu spe11b_structured.msh
```
__Note__: In this last case we use [meshio](https://pypi.org/project/meshio/2.3.5/) to convert
to **vtu** format (not vtk).

__Note__: Though valid, the option `--poromult` will have no effect on spe11-a case.

The file `spe11b_structured_extruded.vtu` should be produced and can be inspected using [paraview](https://www.paraview.org/) or another
vtk enabled 3D reader.
