from math import pi
from argparse import ArgumentParser
import gmsh


def z_offset_at(y: float) -> float:
    """
    Compute the difference in z-coordinate between reference and physical space
    according to eq. (4.1) of the description for SPE11 Version C, given a
    y-coordinate in the reference space.
    """
    f = (y - 2500.0)/2500.0
    return 150.0*(1.0 - f*f)


if __name__ == "__main__":
    parser = ArgumentParser(description="Create the geometry for SPE11 version C")
    parser.add_argument(
        "-s",
        "--mesh-size",
        required=True,
        help="Set the characteristic mesh size in the domain. This will be defined in the"
            "last line of the created .geo file, where you can edit it manually afterwards"
    )
    args = vars(parser.parse_args())

    gmsh.initialize()
    gmsh.open("spe11b.geo")

    print("Reading entities & physical groups from the model")
    entities = {
        "points": [(tag, gmsh.model.getValue(0, tag, [])) for _, tag in gmsh.model.getEntities(0)],
        "curves": [(tag, gmsh.model.getAdjacencies(1, tag)[1]) for _, tag in gmsh.model.getEntities(1)],
        "surfaces": [(tag, gmsh.model.getAdjacencies(2, tag)[1]) for _, tag in gmsh.model.getEntities(2)]
    }
    surface_to_physical_properties = {}
    for dim, index in gmsh.model.getPhysicalGroups(dim=2):
        name = gmsh.model.getPhysicalName(dim, index)
        for tag in gmsh.model.getEntitiesForPhysicalGroup(dim, index):
            assert tag not in surface_to_physical_properties
            surface_to_physical_properties[tag] = {"name": name, "index": index}

    print("Defining new model using the OCC kernel")
    gmsh.model.add("occ_model")
    gmsh.model.setCurrent("occ_model")

    print("Adding points")
    for tag, coordinates in entities["points"]:
        gmsh.model.occ.addPoint(*coordinates, tag=tag)

    print("Adding curves")
    for tag, (p0, p1) in entities["curves"]:
        gmsh.model.occ.addLine(p0, p1, tag=tag)

    print("Adding surfaces")
    for tag, curve_tags in entities["surfaces"]:
        wire_tag = gmsh.model.occ.addCurveLoop(curve_tags, tag=tag)
        stag = gmsh.model.occ.addPlaneSurface([wire_tag], tag=tag)

    print("Rotating model such that the surfaces align with the z-axis")
    gmsh.model.occ.synchronize()
    bbox = gmsh.model.getBoundingBox(-1, -1)
    min, max = tuple(bbox[:3]), tuple(bbox[3:])
    gmsh.model.occ.rotate(
        dimTags=gmsh.model.getEntities(2),
        x=0.5*(max[0] - min[0]), y=min[1], z=0.0,
        ax=1.0, ay=0.0, az=0.0,
        angle=pi/2.0
    )

    print("Copying Surfaces to the back of the domain")
    size_domain_y = 5000.0
    gmsh.model.occ.synchronize()
    frontside_surface_tags = gmsh.model.getEntities(dim=2)
    backside_surface_tags = gmsh.model.occ.copy(frontside_surface_tags)
    gmsh.model.occ.translate(backside_surface_tags, dx=0.0, dy=size_domain_y, dz=0.0)

    def _make_connecting_spline(source_tag: int, target_tag: int) -> tuple[int, list]:
        source = gmsh.model.getValue(0, source_tag, [])
        target = gmsh.model.getValue(0, target_tag, [])
        assert abs(source[0] - target[0]) < 1e-6 \
            and abs(source[2] - target[2]) < 1e-6 \
            and "Expecting source and target to have only differing y coordinates"

        extrusion_length = 5000.0
        num_support_points = 10
        dy = extrusion_length/float(num_support_points+1)  # first/last point are source/target

        support_point_tags = [source_tag]
        for i in range(num_support_points):
            y_reference = source[1] + float(i + 1)*dy
            support_point_tags.append(
                gmsh.model.occ.addPoint(
                    x=source[0],
                    y=source[1] + y_reference,
                    z=source[2] + z_offset_at(y_reference)
                )
            )
        return gmsh.model.occ.addBSpline(support_point_tags + [target_tag]), support_point_tags[1:]

    print("Creating connecting surfaces and volumes")
    physical_volumes = {}
    front_point_to_connecting_spline_index = {}
    front_curve_to_connecting_surface_index = {}
    gmsh.model.occ.synchronize()
    for i, (front, back) in enumerate(zip(frontside_surface_tags, backside_surface_tags)):
        print(f"Creating volume {i+1} of {len(frontside_surface_tags)}", end="\r")
        front_boundary_curves = gmsh.model.getBoundary([front], recursive=False)
        back_boundary_curves = gmsh.model.getBoundary([back], recursive=False)
        assert len(front_boundary_curves) == len(back_boundary_curves)

        bounding_surface_tags = [front[1]]
        for front_curve, back_curve in zip(front_boundary_curves, back_boundary_curves):
            assert front_curve[1] > 0 and back_curve[1] > 0 \
                or front_curve[1] < 0 and back_curve[1] < 0

            abs_front_curve = abs(front_curve[1])
            if abs_front_curve in front_curve_to_connecting_surface_index:
                bounding_surface_tags.append(front_curve_to_connecting_surface_index.get(abs_front_curve))
            else:
                front_curve_points = gmsh.model.getBoundary([front_curve], recursive=False)
                back_curve_points = gmsh.model.getBoundary([back_curve], recursive=False)
                assert len(front_curve_points) == len(back_curve_points)
                assert len(front_curve_points) == 2

                pfront_0, pfront_1 = front_curve_points[0][1], front_curve_points[1][1]
                pback_0, pback_1 = back_curve_points[0][1], back_curve_points[1][1]
                spline1, points1 = front_point_to_connecting_spline_index.get(pfront_0), []
                spline2, points2 = front_point_to_connecting_spline_index.get(pfront_1), []

                if spline1 is None:
                    spline1, points1 = _make_connecting_spline(pfront_0, pback_0)
                    front_point_to_connecting_spline_index[pfront_0] = spline1
                if spline2 is None:
                    spline2, points2 = _make_connecting_spline(pfront_1, pback_1)
                    front_point_to_connecting_spline_index[pfront_1] = spline2

                curve_loop = [spline1, spline2]
                if front_curve[1] < 0:
                    curve_loop = list(reversed(curve_loop))
                wire_tag = gmsh.model.occ.addWire([front_curve[1], *curve_loop, -back_curve[1]])
                bounding_surface_tags.append(gmsh.model.occ.addSurfaceFilling(wire_tag, tag=wire_tag))
                front_curve_to_connecting_surface_index[abs_front_curve] = bounding_surface_tags[-1]

                # remove support points
                gmsh.model.occ.remove(dimTags=[(0, t) for t in points1], recursive=False)
                gmsh.model.occ.remove(dimTags=[(0, t) for t in points2], recursive=False)

        bounding_surface_tags.append(back[1])
        physical_props = surface_to_physical_properties[front[1]]
        physical_index = physical_props["index"]
        surf_loop = gmsh.model.occ.addSurfaceLoop(bounding_surface_tags)
        volume_tag = gmsh.model.occ.addVolume([surf_loop])
        if physical_index not in physical_volumes:
            physical_volumes[physical_index] = {"name": physical_props["name"], "volumes": []}
        physical_volumes[physical_index]["volumes"].append(volume_tag)

    print("Writing .brep geometry file")
    base_filename = "spe11c"
    brep_filename = f"{base_filename}.brep"
    gmsh.model.occ.synchronize()
    gmsh.write(brep_filename)
    print(f"Wrote '{brep_filename}'")

    print("Writing .geo file with physical volume definitions")
    geo_filename = f"{base_filename}.geo"
    with open(geo_filename, "w") as geo_file:
        geo_file.write(f'Merge "{brep_filename}";\n')
        for physical_index, properties in physical_volumes.items():
            volume_tags_str = ",".join(str(t) for t in properties["volumes"])
            geo_file.write(f"Physical Volume({physical_index}) = {{{volume_tags_str}}};\n")
        geo_file.write(f"Characteristic Length{{:}} = {args['mesh_size']};\n")
    print(f"Wrote '{geo_filename}'")
