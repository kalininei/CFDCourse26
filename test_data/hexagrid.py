from hybmeshpack import hmscript as hm
import math

N = 1000
output_name = "hexagrid.vtk"

step = 2*math.sqrt(1.0/6/math.sqrt(3)/N)
c1 = hm.add_rect_contour([0, 0], [1, 1])
g = hm.add_unf_hex_grid([[0, 0], [1, 1]], step, strict=True)
g = hm.exclude_contours(g, c1, "outer")
print(hm.info_grid(g))
hm.export_grid_vtk(g, output_name)
