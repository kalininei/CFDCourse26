from hybmeshpack import hmscript as hm
import math

N = 3_000
output_name = "cylgrid_vertical_3k.vtk"

h = math.exp((2.718385717 - math.log(N)) / 1.998375352)
step = 2 * h

na = int(3.1415 / h)
step = 2 * h

width = 3
height = 6
y_before = 1.5

c1 = hm.add_rect_contour([-width / 2, -y_before], [width / 2, height - y_before])
c1 = hm.partition_contour(c1, "const", step)
c2 = hm.add_circ_contour([0, 0], 0.5, na)
g = hm.triangulate_domain([c1, c2], fill="3")

print(hm.info_grid(g))
hm.export_grid_vtk(g, output_name)
print("Grid was written into " + output_name)
