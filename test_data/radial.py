from hybmeshpack import hmscript as hm
import math

N = 50000
output_name = f"radial_{N}.vtk"
r0 = 0.05
R = 1.0 + r0
coef = math.sqrt(N / 490)

na_outer = int(64 * coef)
c1 = hm.add_circ_contour([0, 0], R, na_outer)
g1 = hm.triangulate_domain(c1, fill="4")

na_inner = int(2 * r0 * na_outer / R)
h = 2 * math.pi * r0 / na_inner
g2 = hm.add_unf_ring_grid([0, 0], r0, r0 + h, na_inner, 1)
g = hm.unite_grids1(g1, g2, 3 * r0, True, buffer_fill="4")
hm.export_grid_vtk(g, output_name)
print(hm.info_grid(g))
print("Grid was written into " + output_name)
