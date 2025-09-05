from hybmeshpack import hmscript as hm
import math

N = 1000
output_name = "trigrid.vtk"

step = math.pow(N/math.exp(0.87025576215809), -1.0/2.02352892745984)
c1 = hm.add_rect_contour([0, 0], [1, 1])
c1 = hm.partition_contour(c1, 'const', step)
g = hm.triangulate_domain(c1, fill='3')
print(hm.info_grid(g))
hm.export_grid_vtk(g, output_name)
print("Grid was written into " + output_name)
