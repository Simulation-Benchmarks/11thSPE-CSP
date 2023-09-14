# SPDX-FileCopyrightText: 2023 Jacques Franc <jacquesfrancdev@gmail.com>
#
# SPDX-License-Identifier: MIT
"""
Generate 1-cells extruded versions of quad and triangle meshes while painting porosity
and permeability on those. These celldata will be tagged respectively PORO and PERM.
The attribute celldata denotes facies label.

Porosity multipliers for *spe11b* can be triggered using --poromult option.
Simple painting of 2D meshes for *spe11a* or *spe11b* can be obtain using --paint option.
Simple painting for (already) 3D meshes for *spe11c* is also accessible through --paint.

"""

import vtk
import numpy as np
import argparse
from functools import partial


class spe11_preprocessing:

    def set_data(self, spe):
        # global values -- default is spe11-a
        self.E = 1.e-2  # dispersion -- unused for now
        self.ratio = {}
        self.modlist = []

        # region values
        self.perm = [4e-11, 5e-10, 1e-9, 2e-9, 4e-9, 1e-8, 1e-18]  # 7th is 0 by specs
        self.poro = [.44, .43, .44, .45, .43, .46, .0000001]  # ditto
        self.multipliers = [0., 0., 0., 0., 0., 0., 0.]
        self.Diw = 1e-9  # except on 7th -- unused for now
        self.Dig = 1.6e-5  # ditto
        self.geom_off = 0.019
        self.depth_off = -1.2
        self.xb = [0, 2.8]
        self.is_aniso = False

        if (spe[0] == "b" or spe[0] == "B") or (spe[0] == "c" or spe[0] == "C"):
            self.perm = [1e-16, 1e-13, 2e-13, 5e-13, 1e-12, 2e-12, 1e-18]  # 7th is 0 by specs
            self.poro = [.1, .2, .2, .2, .25, .35, .0000001]  # ditto
            self.multipliers = [0., 5e4, 0., 0., 5e4, 0., 0.]
            self.Dig = 2.e-8  # ditto
            self.geom_off = 1.0
            self.depth_off = -1200
            self.xb = [0, 8400]
            self.is_aniso = True

    def extrude(self, fname, callback):
        extension = fname.split('.')[-1]
        if extension == "vtk":
            r = vtk.vtkUnstructuredGridReader()
            r.SetFileName(fname)
            r.Update()
        elif extension == "vtu":
            r = vtk.vtkXMLUnstructuredGridReader()
            r.SetFileName(fname)
            r.Update()

        g = r.GetOutput()
        g3 = vtk.vtkUnstructuredGrid()

        self.precompute_fictive_volume_ration(g)
        (g, g3) = callback(g, g3, self.geom_off, self.depth_off)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(fname.split('.')[0] + '_extruded.vtu')
        writer.SetInputData(g3)
        writer.Write()

    # note do also swap y<->z
    def triangle_mesh(self, poromult, g, g3, geom_off, depth_off):
        off = g.GetNumberOfPoints()
        # set points
        points = vtk.vtkPoints()
        for i in range(g.GetNumberOfPoints()):
            pt = g.GetPoint(i)
            points.InsertPoint(i + off, (pt[0], pt[2], pt[1] + depth_off))
            points.InsertPoint(i, (pt[0], pt[2] + geom_off, pt[1] + depth_off))

        g3.SetPoints(points)

        # cells
        wedge = vtk.vtkWedge()
        cells = vtk.vtkCellArray()
        mask = [2, 1, 0]
        for i in range(g.GetNumberOfCells()):
            if g.GetCell(i).GetCellType() == vtk.VTK_TRIANGLE:  # triangles
                tri = g.GetCell(i)
                ids = tri.GetPointIds()
                swap = False
                pt = 3 * [None]
                for j in range(tri.GetNumberOfPoints()):
                    pt[j] = points.GetPoint(tri.GetPointId(j))

                n = np.empty(3)
                tri.ComputeNormal(pt[0], pt[1], pt[2], n)
                if n[1] < 0:
                    # swap indices
                    swap = True
                    # print("swapping ", i)
                for j in range(ids.GetNumberOfIds()):
                    id_ = tri.GetPointId(j)
                    if swap:
                        wedge.GetPointIds().SetId(mask[j], id_)
                        wedge.GetPointIds().SetId(mask[j] + 3, id_ + off)
                    else:
                        wedge.GetPointIds().SetId(j, id_)
                        wedge.GetPointIds().SetId(j + 3, id_ + off)
                cells.InsertNextCell(wedge)

        g3.SetCells(vtk.VTK_WEDGE, cells)

        former_tag = "CellEntityIds"
        self.paint_data(g, g3, poromult, former_tag=former_tag)

        return (g, g3)

    def precompute_fictive_volume_ration(self, g):

        eps = 1e-3
        for i in range(g.GetNumberOfCells()):

            for iface in range(g.GetCell(i).GetNumberOfEdges()):
                pe0 = [0, 0, 0]
                pe1 = [0, 0, 0]
                pts = g.GetCell(i).GetEdge(iface).GetPoints()
                pts.GetPoint(0, pe0)
                pts.GetPoint(1, pe1)

                if pe0[0] == pe1[0] and (
                        np.abs(self.xb[0] - pe0[0]) < eps or np.abs(self.xb[1] - pe0[0]) < eps or np.abs(
                    self.xb[0] - pe1[0]) < eps or np.abs(self.xb[1] - pe1[0]) < eps):
                    area = np.max(np.abs(np.asarray([pe1[1] - pe0[1],
                                                     pe1[2] - pe0[2]])))
                    if g.GetCell(i).GetCellType() == vtk.VTK_TRIANGLE:
                        p0 = [0, 0, 0]
                        p1 = [0, 0, 0]
                        p2 = [0, 0, 0]

                        g.GetCell(i).GetPoints().GetPoint(0, p0)
                        g.GetCell(i).GetPoints().GetPoint(1, p1)
                        g.GetCell(i).GetPoints().GetPoint(2, p2)
                        vc = vtk.vtkTriangle.TriangleArea(p0, p1, p2)
                    elif g.GetCell(i).GetCellType() == vtk.VTK_QUAD:
                        bounds = [0, 0, 0, 0, 0, 0]
                        g.GetCell(i).GetBounds(bounds)
                        lx = (1 if (bounds[0] == bounds[1]) else bounds[1] - bounds[0])
                        ly = (1 if (bounds[2] == bounds[3]) else bounds[3] - bounds[2])
                        lz = (1 if (bounds[4] == bounds[5]) else bounds[5] - bounds[4])
                        vc = lx * ly * lz

                    self.modlist.append(i)
                    self.ratio[i] = area / vc

    # note do also swap y<->z
    def quad_mesh(self, poromult, g, g3, geom_off, depth_off):
        off = g.GetNumberOfPoints()
        # set points
        points = vtk.vtkPoints()
        for i in range(g.GetNumberOfPoints()):
            pt = g.GetPoint(i)
            points.InsertPoint(i + off, (pt[0], pt[2], pt[1] + depth_off))
            points.InsertPoint(i, (pt[0], pt[2] + geom_off, pt[1] + depth_off))

        g3.SetPoints(points)

        # cells
        hexa = vtk.vtkHexahedron()
        cells = vtk.vtkCellArray()
        mask = [0, 1, 2, 3]
        for i in range(g.GetNumberOfCells()):
            if g.GetCell(i).GetCellType() == vtk.VTK_QUAD:  # filter QUAD
                quad = g.GetCell(i)
                ids = quad.GetPointIds()

                for j in range(ids.GetNumberOfIds()):
                    id_ = quad.GetPointId(j)

                    hexa.GetPointIds().SetId(mask[j], id_)
                    hexa.GetPointIds().SetId(mask[j] + 4, id_ + off)
                cells.InsertNextCell(hexa)
        g3.SetCells(vtk.VTK_HEXAHEDRON, cells)

        # transfer data
        self.paint_data(g, g3, poromult)

        return (g, g3)

    def paint_data(self, g, g3, poromult, former_tag="gmsh:physical"):
        for i in range(g.GetCellData().GetNumberOfArrays()):
            array = vtk.vtkTypeInt32Array()

            perm_array = vtk.vtkFloatArray()
            perm_array.SetName("PERM")
            perm_array.SetNumberOfComponents(3)

            poro_array = vtk.vtkFloatArray()
            poro_array.SetName("PORO")

            if g.GetCellData().GetArray(i).GetName() == former_tag:
                array.SetName("attribute")
                for j in range(g.GetCellData().GetArray(i).GetNumberOfTuples()):
                    if g.GetCell(j).GetCellType() in [vtk.VTK_QUAD, vtk.VTK_TRIANGLE, vtk.VTK_HEXAHEDRON,
                                                      vtk.VTK_TETRA]:
                        attribute = g.GetCellData().GetArray(i).GetValue(j)
                        array.InsertNextValue(attribute)
                        perm_array.InsertNextTuple3(self.perm[attribute - 1], self.perm[attribute - 1],
                                                    (0.1 if self.is_aniso else 1.) * self.perm[
                                                        attribute - 1])  # C-numbering
                        if j not in self.modlist or not poromult:
                            poro_array.InsertNextValue(self.poro[attribute - 1])  # C-numbering
                        else:
                            poro_array.InsertNextValue(self.poro[attribute - 1] * (
                                    1 + self.multipliers[attribute - 1] * self.ratio[j]))  # C-numbering
                    # gmsh generated triangle meshes still have line element, so have to skip them or flag them
                    elif g.GetCell(j).GetCellType() in [vtk.VTK_LINE]:
                        array.InsertNextValue(-1)
                        perm_array.InsertNextTuple3(-1., -1., -1.)
                        poro_array.InsertNextValue(-1.)

                g3.GetCellData().AddArray(array)
                g3.GetCellData().AddArray(perm_array)
                g3.GetCellData().AddArray(poro_array)

    def paint_only(self, fname, poromult, former_tag):
        extension = fname.split('.')[-1]

        if extension == "vtk":
            r = vtk.vtkUnstructuredGridReader()
            r.SetFileName(fname)
            r.Update()
        elif extension == "vtu":
            r = vtk.vtkXMLUnstructuredGridReader()
            r.SetFileName(fname)
            r.Update()

        g = r.GetOutput()
        self.paint_data(g, g, poromult, former_tag)

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetFileName(fname.split('.')[0] + '_painted.vtu')
        writer.SetInputData(g)
        writer.Write()

        return g


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("fname", help="input vtk filename")

    parser.add_argument("--tri", help="for use if triangle mesh to be extruded", action="store_true")
    parser.add_argument("--quad", help="for use if quad mesh to be extruded", action="store_true")
    parser.add_argument("--paint", help="paint only", action="store_true")

    parser.add_argument("--spe", help="paint data from spe11 a or b", nargs=1, default='a', required=True)
    #
    parser.add_argument("--poromult", help="apply pore mult", action="store_true")

    args = parser.parse_args()
    preproc = spe11_preprocessing()
    preproc.set_data(args.spe)

    if args.paint:
        preproc.paint_only(args.fname, args.poromult, former_tag="gmsh:physical")
    elif args.tri:
        preproc.extrude(args.fname, partial(preproc.triangle_mesh, args.poromult))
    elif args.quad:
        preproc.extrude(args.fname, partial(preproc.quad_mesh, args.poromult))
    else:
        raise NotImplementedError
