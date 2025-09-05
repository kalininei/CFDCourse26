from hybmeshpack import hmscript as hm

N = 5000
output_name = "cylgrid.vtk"

h = 22.2462 / (N**0.5034)

length = 15
height = 4
cyl_x = 3
cyl_y = -0.5

h0 = h / 3
hb1 = h0 / 4
hb2 = h0 / 1.2

L1 = cyl_x + 0.5
g1 = hm.add_unf_rect_grid([0, -height/2], [L1, height/2],
                          int(L1/h), int(height/h))

g2 = hm.add_unf_rect_grid([L1, -height/2], [length, height/2],
                          int((length-L1)/h0), int(height/h0))
bnd_h = 0.3
na = int(3.1415*(1 + 2*bnd_h)/h0)
c3 = hm.add_circ_contour([cyl_x, cyl_y], 0.5, na)
bnd = hm.partition_segment(0, bnd_h, hb1, hb2)
g3 = hm.build_boundary_grid1(c3, bnd, 'right')

g = hm.unite_grids1(g1, g2, 0.5, empty_holes=True, buffer_fill="4")
g = hm.unite_grids1(g, g3, 0.2, empty_holes=True, buffer_fill="4")
print(hm.info_grid(g))
hm.export_grid_vtk(g, output_name)
print("Grid was written into " + output_name)
