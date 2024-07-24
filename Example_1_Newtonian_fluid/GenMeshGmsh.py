import numpy as np
import matplotlib.pyplot as plt
import h5py
import gmsh

# Dimensions
L = 80
W = 60

gmsh.initialize()
gmsh.model.add("Structured_Mesh")

# Characteristic length
lc = 2

# Number of points in each dimension
num_points_x = 11  # Number of points along the length
num_points_y = 11  # Number of points along the width

# Create points
gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(W, 0, 0, lc, 2)
gmsh.model.geo.addPoint(W, L, 0, lc, 3)
gmsh.model.geo.addPoint(0, L, 0, lc, 4)

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

gmsh.model.mesh.setOrder(2)

gmsh.model.mesh.generate(2)

# Save the mesh to a .msh file
##gmsh.write("Structured_Mesh.msh")

# Run the GUI (optional)
gmsh.fltk.run()

# Save mesh data to HDF5 format
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
element_types, element_tags, node_tags_per_element = gmsh.model.mesh.getElements()

print(len(node_coords), len(node_tags_per_element))

gmsh.model.mesh.generate(2)
node_tags_h, node_coords_h, _ = gmsh.model.mesh.getNodes()
element_types_h, element_tags_h, node_tags_per_element_h = gmsh.model.mesh.getElements()

with h5py.File("Structured_Mesh.h5", "w") as f:
    f.create_dataset("Points", data=node_coords.reshape(-1, 3))
    f.create_dataset("Elements", data=node_tags_per_element[0].reshape(-1, 4))
    f.create_dataset("Points_h", data=node_coords_h.reshape(-1, 3))
    f.create_dataset("Elements_h", data=node_tags_per_element_h[0].reshape(-1, 8))

print(len(node_coords_h), len(node_tags_per_element_h))
# Finalize Gmsh
#gmsh.fltk.run()
gmsh.finalize()
