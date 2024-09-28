import numpy as np
import matplotlib.pyplot as plt
import h5py
import gmsh
import os
import sys
import math

# Dimensions
L = 100
W = 60

gmsh.initialize()
gmsh.model.add("Structured_Mesh")

# Characteristic length
lc = 2

D = 9
S = 10 ### S should be larger than D

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

Pnt_ID = 5
Arc_ID = 5
Cir_ID = 2
centerx_1 = S / 2.0
R = D / 2.0

for i in range(math.floor(L / S)):
    for j in range(math.floor(W / S)):

        CX = centerx_1 + i * S
        CY = centerx_1 + j * S

        gmsh.model.geo.addPoint(CX, CY, 0, lc, Pnt_ID) # CENTER

        gmsh.model.geo.addPoint(CX - R, CY, 0, lc, Pnt_ID + 1)
        gmsh.model.geo.addPoint(CX, CY + R, 0, lc, Pnt_ID + 2)
        gmsh.model.geo.addPoint(CX + R, CY, 0, lc, Pnt_ID + 3)
        gmsh.model.geo.addPoint(CX, CY - R, 0, lc, Pnt_ID + 4)

        gmsh.model.geo.addCircleArc(Pnt_ID + 1, Pnt_ID, Pnt_ID + 2, Arc_ID)
        gmsh.model.geo.addCircleArc(Pnt_ID + 2, Pnt_ID, Pnt_ID + 3, Arc_ID + 1)
        gmsh.model.geo.addCircleArc(Pnt_ID + 3, Pnt_ID, Pnt_ID + 4, Arc_ID + 2)
        gmsh.model.geo.addCircleArc(Pnt_ID + 4, Pnt_ID, Pnt_ID + 1, Arc_ID + 3)

        gmsh.model.geo.addCurveLoop([Arc_ID, Arc_ID + 1, Arc_ID + 2, Arc_ID + 3], Cir_ID)

        Pnt_ID = Pnt_ID + 5
        Arc_ID = Arc_ID + 4
        Cir_ID = Cir_ID + 1

gmsh.model.geo.addPlaneSurface(list(range(1, Cir_ID)), 1)

# Synchronize and generate the mesh
gmsh.model.geo.synchronize()
# gmsh.model.mesh.setOrder(2)

gmsh.model.mesh.generate(2)
gmsh.fltk.run()

# Save the mesh to a .msh file
#gmsh.write("Structured_Mesh.msh")

# Run the GUI (optional)
gmsh.finalize()
