import numpy as np
import matplotlib.pyplot as plt
import h5py as h5py

with h5py.File("Mesh.h5", "r") as f:
    Points = np.array(f["Points"][:])
    Elements = np.array(f["Elements"][:])

# ________ solver _____________
Triplet = []

NumElements = np.shape(Elements)[0]
w_i = np.array([0.17132, 0.36076, 0.46791, 0.17132, 0.36076, 0.46791])
xi_i = np.array([0.93247, 0.66121, 0.23862, -0.93247, -0.66121, -0.23862])

for e in range(0, NumElements):
    B1 = np.zeros((4, 9))
    B2 = np.zeros((4, 9))

    pnt_x = Points[Elements[e, :], 0]
    pnt_y = Points[Elements[e, :], 1]

    pnt_x1 = 1 / 2. * (pnt_x + pnt_x[[1,2,3,0]])
    pnt_y1 = 1 / 2. * (pnt_y + pnt_y[[1,2,3,0]])

    pnt_x = np.array([pnt_x[0], pnt_x1[0], pnt_x[1], pnt_x1[3], 1 / 2. * (pnt_x1[1] + pnt_x1[3]), pnt_x1[1], pnt_x[3], pnt_x1[2], pnt_x[2]])
    pnt_y = np.array([pnt_y[0], pnt_y1[0], pnt_y[1], pnt_y1[3], 1 / 2. * (pnt_y1[1] + pnt_y1[3]), pnt_y1[1], pnt_y[3], pnt_y1[2], pnt_y[2]])

    for i in range(0, 6):
        for j in range(0, 6):
            xi_e = xi_i[i]
            eta_e = xi_i[j]
            Psi = np.array(
                [1 / 4.0 * (1 - xi_e) * (1 - eta_e),
                 1 / 4.0 * (1 + xi_e) * (1 - eta_e),
                 1 / 4.0 * (1 + xi_e) * (1 + eta_e),
                 1 / 4.0 * (1 - xi_e) * (1 + eta_e)],
            )
            Phi_xi = np.array(
                [
                    1/4 * (2 * xi_e - 1) * (eta_e - 1) * eta_e,
                    -xi_e * (eta_e - 1) * eta_e,
                    1/4 * (2 * xi_e + 1) * (eta_e - 1) * eta_e,
                    -1/2 * (2 * xi_e - 1) * (eta_e ** 2 - 1),
                    2 * xi_e * (eta_e ** 2 - 1),
                    -1/2 * (2 * xi_e + 1) * (eta_e**2 - 1),
                    1/4 * (2 * xi_e - 1) * eta_e * (eta_e + 1),
                    -xi_e * eta_e * (eta_e + 1),
                    1/4 * (2 * xi_e + 1) * eta_e * (eta_e + 1)
                ]
            )
            Phi_eta = np.array(
                [
                    1/4 * (xi_e - 1) * xi_e * (2 * eta_e - 1),
                    -1/2 * (xi_e**2 - 1) * (2**eta_e - 1),
                    1/4 * xi_e * (xi_e + 1) * (2 * eta_e - 1),
                    -(xi_e - 1) * xi_e * eta_e,
                    2 * (xi_e**2 - 1) * eta_e,
                    -xi_e * (xi_e + 1) * eta_e,
                    1/4 * (xi_e - 1) * xi_e * (2 * eta_e + 1),
                    -1/2 * (xi_e**2 - 1) * (2 * eta_e + 1),
                    1/4 * (2 * xi_e + 1) * eta_e * (eta_e + 1)
                ]
            )
            x_xi = np.dot(Phi_xi, pnt_x)
            y_xi = np.dot(Phi_xi, pnt_y)
            x_eta = np.dot(Phi_eta, pnt_x)
            y_eta = np.dot(Phi_eta, pnt_x)
            Jacob = np.array(
                [[x_xi, y_xi],
                 [x_eta, y_eta]]
            )
            #print(np.shape(np.linalg.inv(Jacob)))
            #print(np.array([Phi_xi, Phi_eta]))
            AA = np.linalg.inv(Jacob) @ np.array([Phi_xi, Phi_eta])
            Phi_x = np.array([AA[0, :]])
            Phi_y = np.array([AA[1, :]])
            B1 = B1 + w_i[i] * w_i[j] * \
                (np.transpose([Psi]) @ Phi_x) * np.linalg.det(Jacob)
            B2 = B2 + w_i[i] * w_i[j] * \
                (np.transpose([Psi]) @ Phi_y) * np.linalg.det(Jacob)
