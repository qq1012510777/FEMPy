import numpy as np
import matplotlib.pyplot as plt
import h5py
import gmsh
import os
import sys

gmsh.initialize()
gmsh.model.add("Mesh")

# Characteristic length
lc = 2

# Create points
gmsh.model.geo.addPoint(0, 0, 0, lc, 1)
gmsh.model.geo.addPoint(10, -5, 0, lc, 2)
gmsh.model.geo.addPoint(10, 15, 0, lc, 3)
gmsh.model.geo.addPoint(0, 10, 0, lc, 4)

# Create lines
gmsh.model.geo.addLine(1, 2, 1)
gmsh.model.geo.addLine(2, 3, 2)
gmsh.model.geo.addLine(3, 4, 3)
gmsh.model.geo.addLine(4, 1, 4)

# Create curve loop and plane surface
gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)
gmsh.model.geo.addPlaneSurface([1], 1)

# Synchronize and generate the mesh
gmsh.model.geo.synchronize()

# gmsh.model.mesh.setOrder(2)

gmsh.model.mesh.generate(2)

# Save the mesh to a .msh file
#gmsh.write("Structured_Mesh.msh")

# Run the GUI (optional)
# gmsh.fltk.run()

# Save mesh data to HDF5 format
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
element_types, element_tags, node_tags_per_element = gmsh.model.mesh.getElements()
Points = node_coords.reshape(-1, 3)
Elements = node_tags_per_element[1].reshape(-1, 3)
# print(len(node_coords), len(node_tags_per_element))
# gmsh.model.mesh.generate(2)
# node_tags_h, node_coords_h, _ = gmsh.model.mesh.getNodes()
# element_types_h, element_tags_h, node_tags_per_element_h = gmsh.model.mesh.getElements()

with h5py.File("Mesh.h5", "w") as f:
    f.create_dataset("Points", data=Points)
    f.create_dataset("Elements", data=Elements)

# print(len(node_coords_h), len(node_tags_per_element_h))
# Finalize Gmsh
gmsh.fltk.run()

gmsh.finalize()
