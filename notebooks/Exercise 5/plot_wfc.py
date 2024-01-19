import pyvista as pv
import numpy as np

x, y, z = np.loadtxt("psi2_3d_z2.dat", unpack=True)
points = np.transpose(np.array([x, y, z]))

dist = np.sqrt(points[:, 0] ** 2 + points[:, 1] ** 2 + points[:, 2] ** 2)
point_cloud = pv.PolyData(points)

point_cloud["distance"] = dist
axes = pv.Axes(show_actor=True, actor_scale=2.0, line_width=5)
axes.origin = (0.0, 0.0, 0.0)
# Create a plotter
plotter = pv.Plotter(notebook=False, off_screen=True)
plotter.add_actor(axes.actor)

for i in range(60):
    
    rotate = point_cloud.rotate_z(6 * i, point=axes.origin, inplace=True)
    cloud = plotter.add_mesh(rotate, render_points_as_spheres=True,  style='points_gaussian', opacity=0.7, point_size=5,emissive=True)
    plotter.camera_position = "xz"
    plotter.screenshot(plotter.screenshot("./psi210/frame_{:03d}.png".format(i)))
    plotter.remove_actor(cloud)

import matplotlib.pyplot as plt
plt.hist(dist, 50)
plt.savefig("baka.png")