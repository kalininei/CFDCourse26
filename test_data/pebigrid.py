from hybmeshpack import hmscript as hm
import math

N = 1000
output_name = "pebigrid.vtk"

step = math.sqrt(1.4/N)
c1 = hm.add_rect_contour([0, 0], [1, 1])
c1 = hm.partition_contour(c1, 'const', step)
g = hm.pebi_fill(c1)
print(hm.info_grid(g))
hm.export_grid_vtk(g, output_name)
print("Grid was written into " + output_name)
