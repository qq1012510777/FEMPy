import numpy as np
import matplotlib.pyplot as plt
import h5py
import gmsh
import os
import sys

# Dimensions
L = 80
W = 60

gmsh.initialize()
gmsh.model.add("Structured_Mesh")

# Characteristic length
lc = 2

# Number of points in each dimension
if len(sys.argv) > 1:
    num_points_x = int(sys.argv[1])  # Number of points along the length
else:
    num_points_x = 11

if len(sys.argv) > 2:
    num_points_y = int(sys.argv[2])  # Number of points along the width
else:
    num_points_y = 11

# Create points
gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(L, 0, 0, lc, 2)
gmsh.model.geo.addPoint(L, W, 0, lc, 3)
gmsh.model.geo.addPoint(0, W, 0, lc, 4)

# Create lines
gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

# Create curve loop and plane surface
gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
gmsh.model.geo.addPlaneSurface([1], 1)

# Apply transfinite mesh constraints to the lines
gmsh.model.geo.mesh.setTransfiniteCurve(1, num_points_x)
gmsh.model.geo.mesh.setTransfiniteCurve(2, num_points_y)
gmsh.model.geo.mesh.setTransfiniteCurve(3, num_points_x)
gmsh.model.geo.mesh.setTransfiniteCurve(4, num_points_y)

# Apply transfinite surface constraints
gmsh.model.geo.mesh.setTransfiniteSurface(1)
gmsh.model.geo.mesh.setRecombine(2, 1)

# Synchronize and generate the mesh
gmsh.model.geo.synchronize()

# gmsh.model.mesh.setOrder(2)

gmsh.model.mesh.generate(2)

# Save the mesh to a .msh file
gmsh.write("Structured_Mesh.msh")

# Run the GUI (optional)
# gmsh.fltk.run()

# Save mesh data to HDF5 format
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
element_types, element_tags, node_tags_per_element = gmsh.model.mesh.getElements()
Points = node_coords.reshape(-1, 3)
Elements = node_tags_per_element[1].reshape(-1, 4)-1
# print(len(node_coords), len(node_tags_per_element))
# gmsh.model.mesh.generate(2)
# node_tags_h, node_coords_h, _ = gmsh.model.mesh.getNodes()
# element_types_h, element_tags_h, node_tags_per_element_h = gmsh.model.mesh.getElements()
# with h5py.File("Structured_Mesh.h5", "w") as f:
#     f.create_dataset("Points", data=node_coords.reshape(-1, 3))
#     f.create_dataset("Elements", data=node_tags_per_element[0].reshape(-1, 4))
#     f.create_dataset("Points_h", data=node_coords_h.reshape(-1, 3))
#     f.create_dataset("Elements_h", data=node_tags_per_element_h[0].reshape(-1, 8))
# print(len(node_coords_h), len(node_tags_per_element_h))
# Finalize Gmsh
# gmsh.fltk.run()

gmsh.finalize()

gmsh.initialize()
gmsh.model.add("Structured_Mesh_h")

gmsh.open("./Structured_Mesh.msh")

gmsh.model.mesh.setOrder(2)

gmsh.write("Structured_Mesh.msh")

node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
element_types, element_tags, node_tags_per_element = gmsh.model.mesh.getElements()

Points_h = node_coords.reshape(-1, 3)
Elements_h = node_tags_per_element[1].reshape(-1, 9)-1
Elements_h = np.array([Elements_h[:, 0], Elements_h[:, 4], Elements_h[:, 1], Elements_h[:, 7],
                      Elements_h[:, 8], Elements_h[:, 5], Elements_h[:, 3], Elements_h[:, 6], Elements_h[:, 2]])
Elements_h = np.transpose(Elements_h)
# print(Points_h, Elements_h)
gmsh.finalize()

with h5py.File("Structured_Mesh.h5", "w") as f:
    f.create_dataset("Points", data=Points)
    f.create_dataset("Elements", data=Elements)
    f.create_dataset("Points_h", data=Points_h)
    f.create_dataset("Elements_h", data=Elements_h)

# fig,ax = plt.subplots()
# ax.plot(Points[Elements[0, :], 0], Points[Elements[0, :], 1])
# for i in range(4):
#     ax.text(Points[Elements[0, i], 0], Points[Elements[0, i], 1], str(Elements[0, i]) + ', ' + str(i + 1))
# plt.show()

fig,ax=plt.subplots(1, 2)

for h in range(np.shape(Elements_h)[0]):
    ax[0].plot(Points_h[Elements_h[h, :], 0], Points_h[Elements_h[h, :], 1], 'o')
    for i in range(np.shape(Elements_h)[1]):
        ax[0].text(Points_h[Elements_h[h, i], 0], Points_h[Elements_h[h, i], 1], str(Elements_h[h, i]))# + ', ' + str(i + 1))

for h in range(np.shape(Elements)[0]):
    ax[1].plot(Points[Elements[h, :], 0], Points[Elements[h, :], 1], 'o')
    for i in range(np.shape(Elements)[1]):
        ax[1].text(Points[Elements[h, i], 0], Points[Elements[h, i], 1], str(Elements[h, i]))# + ', ' + str(i + 1))

plt.show()

# Specify the file path
file_path = "./Structured_Mesh.msh"

# Check if the file exists
if os.path.exists(file_path):
    # Delete the file
    os.remove(file_path)
    print(f"{file_path} has been deleted.")
else:
    print(f"{file_path} does not exist.")