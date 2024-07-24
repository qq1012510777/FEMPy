import numpy as np
import matplotlib.pyplot as plt
import h5py as h5py

# Dimensions
L = 80
W = 60

# Number of points in each dimension
num_points_x = 11  # num_points_x intervals along the length
num_points_y = 11  # num_points_y intervals along the width

# Generate mesh grid
x = np.linspace(0, L, num_points_x)
y = np.linspace(0, W, num_points_y)
X, Y = np.meshgrid(x, y)

# Flatten the grid points into a list of coordinates
points = np.vstack([X.ravel(), Y.ravel()]).T

# Initialize the list for elements (triangles) and edges
elements = []

# Generate triangles and edges
for i in range(num_points_y - 1):
    for j in range(num_points_x - 1):
        # Indices of the corners of the current cell
        p1 = i * num_points_x + j
        p2 = p1 + 1
        p3 = p1 + num_points_x
        p4 = p3 + 1

        # Two triangles for each cell
        elements.append([p1, p2, p4, p3])
       

# Convert lists to arrays
points = np.array(points)
elements = np.array(elements)

# Print the generated arrays
# print("Points:")
# print(points)
# print("\nElements (squares):")
# print(elements)

# Plotting
plt.figure(figsize=(10, 7))
#plt.plot(points[:, 0], points[:, 1], 'ko')  # plot points

# Plot elements (squares)
for square in elements:
    square_points = points[square]
    square_points = np.vstack([square_points, square_points[0]])  # close the square
    plt.plot(square_points[:, 0], square_points[:, 1], 'b-')

plt.grid(True)
plt.title('2D Structured Mesh')
plt.xlabel('Length (L)')
plt.ylabel('Width (W)')

# NumPnts = np.shape(points)[0]
# for i in range(0, NumPnts):
#     plt.text(points[i, 0], points[i, 1], str(i))

plt.show()

with h5py.File("Mesh.h5", "w") as f:
    f.create_dataset("Points", data=points)
    f.create_dataset("Elements", data=elements)

