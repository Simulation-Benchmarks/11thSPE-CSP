# SPDX-FileCopyrightText: 2023 Dennis Gläser <dennis.glaeser@iws.uni-stuttgart.de>
#
# SPDX-License-Identifier: MIT
"""
Generate a structured mesh in gmsh file format version 2 (.msh) for one of the SPE11 variant.
This script is to be executed from the same folder that contains the .geo files of the variants.
"""
from __future__ import annotations

import os
import sys
import argparse
import textwrap
import typing
import subprocess
import itertools

import gmsh

from make_spe11c_geo import z_offset_at


PHYSICAL_INDEX_FACIES_7 = 7
PHYSICAL_INDEX_OUTSIDE_OF_DOMAIN = 1000
PHYSICAL_NAME_OUTSIDE_OF_DOMAIN = str(PHYSICAL_INDEX_OUTSIDE_OF_DOMAIN)


def _is_in_bbox(position, min, max) -> bool:
    return all(position[i] <= max[i] and position[i] >= min[i] for i in range(dim))


def _get_bounding_box(model, entity_dim=-1, entity_tag=-1) -> tuple:
    bbox = model.getBoundingBox(entity_dim, entity_tag)
    return tuple(bbox[:3]), tuple(bbox[3:])


def _get_variant_geo_file(variant: str) -> str:
    return f"spe11{variant.lower()}.geo"


def _get_domain_dimensions(variant: str) -> tuple:
    return {
        "A": ((0.0, 0.0, 0.0), (2.8, 1.2, 0.0)),
        "B": ((0.0, 0.0, 0.0), (8400, 1200.0, 0.0)),
        "C": ((0.0, 0.0, 0.0), (8400.0, 5000.0, 1354.0))
    }[variant]


class StructuredLattice:
    def __init__(self, min: tuple, max: tuple, num_cells: tuple) -> None:
        self._num_cells = num_cells
        self._dim = len(num_cells)
        self._dx = [(max[i] - min[i])/float(num_cells[i]) for i in range(self._dim)]
        self._points = [
            (
                min[0] + i*self._dx[0],
                min[1] + j*self._dx[1],
                min[2] + k*(self._dx[2] if self._dim == 3 else 0.0)
            )
            for k in range(num_cells[2] + 1 if self._dim == 3 else 1)
            for j in range(num_cells[1] + 1)
            for i in range(num_cells[0] + 1)
        ]

    @property
    def number_of_points(self) -> int:
        points = (self._num_cells[0]+1)*(self._num_cells[1]+1)
        if self._dim == 2:
            return points
        return points*(self._num_cells[2]+1)

    @property
    def number_of_cells(self) -> int:
        cells = self._num_cells[0]*self._num_cells[1]
        if self._dim == 2:
            return cells
        return cells*self._num_cells[2]

    @property
    def points(self) -> list:
        return self._points

    @property
    def cells(self) -> typing.Iterable:
        return (
            (i, j, k)
            for k in range(num_cells[2] if dim == 3 else 1)
            for j in range(num_cells[1])
            for i in range(num_cells[0])
        )

    def corners(self, cell: tuple) -> tuple:
        def _get_quad_corners(p0: int) -> tuple:
            return (
                p0,
                p0 + 1,
                p0 + 1 + num_cells[0] + 1,
                p0 + num_cells[0] + 1
            )

        if self._dim == 2:
            p0 = cell[1]*(num_cells[0] + 1) + cell[0]
            return _get_quad_corners(p0)

        x, y, z = cell
        nx, ny = num_cells[0], num_cells[1]
        z_layer_offset = (nx+1)*(ny+1)
        p0 = z*z_layer_offset + y*(nx + 1) + x
        return _get_quad_corners(p0) + _get_quad_corners(p0 + z_layer_offset)

    def center(self, cell: tuple) -> tuple:
        result = tuple([0.0 for _ in range(3)])
        corners = self.corners(cell)
        for pidx in self.corners(cell):
            result = tuple([
                result[i] + self._points[pidx][i]
                for i in range(self._dim)
            ])
        return tuple([(result[i]/len(corners) if i < self._dim else 0.0) for i in range(3)])


class PhysicalIndexMapper:
    def __init__(self, variant: str) -> None:
        gmsh.initialize()

        assert os.path.exists(_get_variant_geo_file(variant))
        gmsh.open(_get_variant_geo_file(variant))
        self._model_name = gmsh.model.getCurrent()
        self._variant = variant

        if variant == "C":
            gmsh.open(_get_variant_geo_file("B"))
            self._2d_model_name = gmsh.model.getCurrent()

        gmsh.model.setCurrent(self._model_name)
        self._physical_groups = self._with_model_for_physical_index_queries(
            lambda: self._read_physical_groups()
        )

    def physical_index(self, position: tuple) -> int:
        position = self._project_to_model_for_index_queries(position)
        index = self._with_model_for_physical_index_queries(
            lambda: self._get_physical_index(position)
        )
        return PHYSICAL_INDEX_OUTSIDE_OF_DOMAIN if index is None else index

    def physical_groups(self, dim: int) -> dict:
        return {
            self._get_physical_name(dim, tag): tag
            for _, tag in self._physical_groups
        }

    def _read_physical_groups(self) -> dict:
        return {
            (d, t): gmsh.model.getEntitiesForPhysicalGroup(dim=d, tag=t)
            for d, t in gmsh.model.getPhysicalGroups()
        }

    def _get_physical_name(self, dim: int, tag: int) -> str:
        name = gmsh.model.getPhysicalName(dim, tag)
        return name if name else str(tag)

    def _with_model_for_physical_index_queries(self, action):
        if self._variant == "C":
            gmsh.model.setCurrent(self._2d_model_name)
        result = action()
        if self._variant == "C":
            gmsh.model.setCurrent(self._model_name)
        return result


    def _project_to_model_for_index_queries(self, position: tuple) -> tuple:
        if self._variant == "C":
            return (
                position[0],
                position[2] - z_offset_at(position[1]),
                0.0
            )
        return position

    def _get_physical_index(self, position: tuple) -> typing.Optional[int]:
        for dim, tag in gmsh.model.getEntities(self._query_model_dimension()):
            min, max = _get_bounding_box(gmsh.model, dim, tag)
            if _is_in_bbox(position, min, max):
                if gmsh.model.isInside(dim, tag, position):
                    return self._get_entity_physical_index(dim, tag)
        return None

    def _get_entity_physical_index(self, dim: int, tag: int) -> typing.Optional[int]:
        for (physical_dim, physical_index), entity_tags in self._physical_groups.items():
            if physical_dim == dim and tag in entity_tags:
                return physical_index
        return None

    def _model_dimension(self) -> int:
        return 3 if self._variant == "C" else 2

    def _query_model_dimension(self) -> int:
        return 2


class FilteredLattice:
    def __init__(self,
                 lattice: StructuredLattice | FilteredLattice,
                 physical_cell_indices: list[int],
                 exclude_physical_index: int) -> None:
        self._lattice = lattice
        self._physical_cell_indices = physical_cell_indices
        self._exclude_physical_index = exclude_physical_index
        self._create_mappings()

    @property
    def number_of_cells(self) -> int:
        return len(self._cells)

    @property
    def number_of_points(self) -> int:
        return len(self._point_map_to_lattice)

    @property
    def cells(self) -> typing.Iterable:
        return self._cells

    @property
    def points(self) -> list:
        return self._points

    def corners(self, cell: tuple) -> typing.Iterable:
        return (self._point_map_from_lattice[p_idx] for p_idx in self._lattice.corners(cell))

    @property
    def physical_cell_indices(self) -> list:
        return list(
            self._physical_cell_indices[self._cell_index_map[i]]
            for i in range(self.number_of_cells)
        )

    def _create_mappings(self) -> None:
        self._cell_index_map = [i for i in range(self._lattice.number_of_cells) if self._is_included(i)]
        self._cells = [cell for i, cell in enumerate(self._lattice.cells) if self._is_included(i)]
        included_points = [False for _ in range(self._lattice.number_of_points)]
        for i, cell in enumerate(self._lattice.cells):
            if self._is_included(i):
                for corner_idx in self._lattice.corners(cell):
                    included_points[corner_idx] = True
        self._point_map_from_lattice = list(itertools.accumulate(included_points, initial=0))[:-1]
        self._point_map_to_lattice = [i for i in range(self._lattice.number_of_points) if included_points[i]]
        self._points = [
            self._lattice.points[self._point_map_to_lattice[i]]
            for i in range(self.number_of_points)
        ]

    def _is_included(self, cell_index: int) -> bool:
        return self._physical_cell_indices[cell_index] != self._exclude_physical_index



parser = argparse.ArgumentParser(description="Create a structured gmsh grid for one of the SPE11 variants")
parser.add_argument(
    "-v", "--variant",
    required=True,
    choices=["A", "B", "C"],
    help="Specify the variant for which to produce the grid"
)
parser.add_argument(
    "-r", "--remove-cells-in-facies-7",
    required=False,
    action="store_true",
    help="Remove all cells within the regions of facies 7"
)
parser.add_argument("-nx", "--number-of-cells-x", required=True, help="Desired number of cells in x-direction")
parser.add_argument("-ny", "--number-of-cells-y", required=True, help="Desired number of cells in y-direction")
parser.add_argument("-nz", "--number-of-cells-z", required=False, help="Desired number of cells in z-direction")
args = vars(parser.parse_args())

variant = args["variant"]
dim = 3 if variant == "C" else 2
if variant == "C":
    if args["number_of_cells_z"] is None:
        sys.exit("Variant C is three-dimensional, please provide the number of cells in z-direction")
    if not os.path.exists(_get_variant_geo_file(variant)):
        print("Could not find geometry of version C. Invoking the generator script to create it.")
        subprocess.run(["python3", "make_spe11c_geo.py", "-s", "100"], check=True)
        assert os.path.exists(_get_variant_geo_file(variant))

num_cells = tuple([int(args[f"number_of_cells_{['x', 'y', 'z'][i]}"]) for i in range(dim)])
gmsh_cell_type = (
    3
    if dim == 2  # quadrangle,
    else  5   # hexahedron
)

physical_index_mapper = PhysicalIndexMapper(variant)
lattice: StructuredLattice | FilteredLattice = StructuredLattice(
    *_get_domain_dimensions(variant),
    num_cells=num_cells
)
num_cells_total = lattice.number_of_cells

print("Determining physical groups for all cells")
physical_cell_indices = [0 for _ in range(num_cells_total)]
for cell_index, cell in enumerate(lattice.cells):
    print(f"Checking cell {cell_index} of {num_cells_total}", end="\r")
    physical_cell_indices[cell_index] = physical_index_mapper.physical_index(lattice.center(cell))

print("Removing all cells outside of the domain")
inside_lattice = FilteredLattice(lattice, physical_cell_indices, PHYSICAL_INDEX_OUTSIDE_OF_DOMAIN)
inside_physical_cell_indices = inside_lattice.physical_cell_indices

if args["remove_cells_in_facies_7"]:
    print("Removing all cells in the facies 7 layers")
    filtered_lattice = FilteredLattice(inside_lattice, inside_physical_cell_indices, PHYSICAL_INDEX_FACIES_7)
    filtered_physical_cell_indices = filtered_lattice.physical_cell_indices
else:
    filtered_lattice = inside_lattice
    filtered_physical_cell_indices = inside_physical_cell_indices

msh_file_name = os.path.splitext(_get_variant_geo_file(variant))[0] + "_structured.msh"
print(f"Writing mesh file '{msh_file_name}'")
with open(msh_file_name, "w") as msh_file:
    msh_file.write(textwrap.dedent("""
        $MeshFormat
        2.2 0 8
        $EndMeshFormat
        $PhysicalNames
        {}
    """.format(len(physical_index_mapper.physical_groups(dim))).lstrip("\n")))
    msh_file.write("{}".format(
        "\n".join(f'{dim} {tag} "{name}"'
        for name, tag in physical_index_mapper.physical_groups(dim).items())
    ))
    msh_file.write(textwrap.dedent(f"""
        $EndPhysicalNames
        $Nodes
        {filtered_lattice.number_of_points}
    """))
    for count, p in enumerate(filtered_lattice.points):
        msh_file.write(f"{count+1} {' '.join(str(c) for c in p)}\n".format())
    msh_file.write(textwrap.dedent(f"""
        $EndNodes
        $Elements
        {filtered_lattice.number_of_cells}
    """).lstrip("\n"))
    for cell_index, cell in enumerate(filtered_lattice.cells):
        phys_index = filtered_physical_cell_indices[cell_index]
        msh_file.write(f"{cell_index+1} {gmsh_cell_type} 2 {phys_index} {phys_index} ")
        msh_file.write(" ".join(str(i + 1) for i in filtered_lattice.corners(cell)))
        msh_file.write("\n")
    msh_file.write("\n$EndElements")
