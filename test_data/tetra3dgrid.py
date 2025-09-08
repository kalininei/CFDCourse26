from hybmeshpack import hmscript as hm
import xml.etree.ElementTree as ET

output_name = "tetra3dgrid.vtk"

step = 0.125
c1 = hm.add_rect_contour([0, 0], [1, 1])
c1 = hm.partition_contour(c1, "const", step)
g = hm.triangulate_domain(c1, fill="3")
hm.export_grid_vtk(g, output_name)
# edges->vertices
it = iter(list(hm.tab_grid2(g, "edge_vert")))
for ind, (n1, n2) in enumerate(zip(it, it)):
    print(ind, n1, n2)
quit()

# g = hm.extrude_grid(g, [0, 1.0])
# gf = hm.tetrahedral_fill(g)
# hm.export3d_grid_vtk(gf, output_name)

# vertices
vertices_indices = []
vertices_coo = []
it = iter(list(hm.tab_grid2(g, "vert")))
for x, y in zip(it, it):
    # xy z+
    vertices_coo.extend([x, y, 1])
    # xy z-
    vertices_coo.extend([x, 1 - y, 0])
    # yz x+
    vertices_coo.extend([1, y, 1 - x])
    # yz x-
    vertices_coo.extend([0, y, x])
    # xz y+
    vertices_coo.extend([x, 1, 1 - y])
    # xz y-
    vertices_coo.extend([x, 0, y])

# edges->vertices
edge_vert = []
it = iter(list(hm.tab_grid2(g, "edge_vert")))
for n1, n2 in zip(it, it):
    for i in range(6):
        edge_vert.append(6 * n1 + i)
        edge_vert.append(6 * n2 + i)

# faces->edges
face_edge = []
it = iter(list(hm.tab_grid2(g, "cell_edge")))
for n1, n2, n3 in zip(it, it, it):
    for i in range(6):
        face_edge.extend([3, 6 * n1 + i, 6 * n2 + i, 6 * n3 + i])


# create hmc
root = ET.Element("HybMeshData")
srf = ET.SubElement(root, "SURFACE3D")
srf.attrib["name"] = "Grid3D_1_surface"
ET.SubElement(srf, "N_VERTICES").text = str(len(vertices_coo) / 3)
ET.SubElement(srf, "N_EDGES").text = str(len(edge_vert) / 2)
ET.SubElement(srf, "N_FACES").text = str(len(face_edge) / 4)
vrt = ET.SubElement(srf, "VERTICES")
coo = ET.SubElement(vrt, "COORDS")
coo.text = " ".join([str(x) for x in vertices_coo])
coo.attrib["type"] = "double"
coo.attrib["format"] = "ascii"
edg = ET.SubElement(srf, "EDGES")
vc = ET.SubElement(edg, "VERT_CONNECT")
vc.text = " ".join([str(x) for x in edge_vert])
vc.attrib["type"] = "int"
vc.attrib["format"] = "ascii"
fc = ET.SubElement(srf, "FACES")
ef = ET.SubElement(fc, "EDGE_CONNECT")
ef.text = " ".join([str(x) for x in face_edge])
ef.attrib["type"] = "int"
ef.attrib["format"] = "ascii"
ef.attrib["dim"] = "variable"


et = ET.ElementTree(root)
sname = output_name + ".hmc"
et.write(sname, "utf-8", True)

loaded_srf = hm.import3d_surface_hmc(sname)
# hm.export3d_surface_vtk(loaded_srf, output_name)
print(loaded_srf)
print(hm.info_surface(loaded_srf))
gf = hm.tetrahedral_fill(loaded_srf)
hm.export3d_grid_vtk(gf, output_name)
