import numpy as np
import matplotlib.pyplot as plt
import h5py as h5py
import scipy.sparse as sp
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix, coo_matrix

with h5py.File("Mesh.h5", "r") as f:
    Points = np.array(f["Points"][:])
    Elements = np.array(f["Elements"][:])

# ____________ parameters ____________-
mu = 1.0e3

# ________ solver _____________
NumPnts = np.shape(Points)[0]

K = coo_matrix(shape=())
b = coo_matrix(shape=())

NumElements = np.shape(Elements)[0]
w_i = np.array([0.17132, 0.36076, 0.46791, 0.17132, 0.36076, 0.46791])
xi_i = np.array([0.93247, 0.66121, 0.23862, -0.93247, -0.66121, -0.23862])

for e in range(0, NumElements):
    B1 = np.zeros((4, 9))
    B2 = np.zeros((4, 9))

    D11 = np.zeros((9, 9))
    D12 = np.zeros((9, 9))
    D21 = np.zeros((9, 9))
    D22 = np.zeros((9, 9))

    C1 = np.zeros((9, 4))
    C2 = np.zeros((9, 4))

    F1 = np.zeros((9, 1))
    F2 = np.zeros((9, 1))

    pnt_x = Points[Elements[e, :], 0]
    pnt_y = Points[Elements[e, :], 1]
    
    pnt_x1 = 1 / 2. * (pnt_x + pnt_x[[1, 2, 3, 0]])
    pnt_y1 = 1 / 2. * (pnt_y + pnt_y[[1, 2, 3, 0]])

    pnt_x = np.array([pnt_x[0], pnt_x1[0], pnt_x[1], pnt_x1[3], 1 / 2. *
                     (pnt_x1[1] + pnt_x1[3]), pnt_x1[1], pnt_x[3], pnt_x1[2], pnt_x[2]])
    pnt_y = np.array([pnt_y[0], pnt_y1[0], pnt_y[1], pnt_y1[3], 1 / 2. *
                     (pnt_y1[1] + pnt_y1[3]), pnt_y1[1], pnt_y[3], pnt_y1[2], pnt_y[2]])

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
            # print(np.shape(np.linalg.inv(Jacob)))
            # print(np.array([Phi_xi, Phi_eta]))
            AA = np.linalg.inv(Jacob) @ np.array([Phi_xi, Phi_eta])
            Phi_x = np.array([AA[0, :]])
            Phi_y = np.array([AA[1, :]])

            deter_J = np.linalg.det(Jacob)

            B1 = B1 + w_i[i] * w_i[j] * \
                (np.transpose([Psi]) @ Phi_x) * deter_J
            B2 = B2 + w_i[i] * w_i[j] * \
                (np.transpose([Psi]) @ Phi_y) * deter_J

            D11 = D11 + mu * w_i[i] * w_i[j] * (2 * np.transpose(
                [Phi_x]) @ Phi_x + np.transpose([Phi_y]) @ Phi_y) * deter_J
            D12 = D12 + mu * w_i[i] * w_i[j] * (np.transpose(
                [Phi_x]) @ Phi_y) * deter_J
            D21 = D21 + mu * w_i[i] * w_i[j] * (np.transpose(
                [Phi_y]) @ Phi_x) * deter_J
            D22 = D22 + mu * w_i[i] * w_i[j] * (2 * np.transpose(
                [Phi_y]) @ Phi_y + np.transpose([Phi_x]) @ Phi_x) * deter_J

            C1 = C1 + mu * w_i[i] * w_i[j] * \
                (np.transpose([Phi_x]) @ [Psi]) * deter_J
            C2 = C2 + mu * w_i[i] * w_i[j] * \
                (np.transpose([Phi_y]) @ [Psi]) * deter_J

    # _________ boundary condition: pressure ______________
    pnt_x = Points[Elements[e, :], 0]
    pnt_y = Points[Elements[e, :], 1]
    for l in range(0, 4):
        if (pnt_x[l] == 0 and pnt_x[(l + 1) % 4] == 0) or (pnt_x[l] == 80. and pnt_x[(l + 1) % 4] == 80.):
            p_i = 1e3
            if (pnt_x[l] == 80. and pnt_x[(l + 1) % 4] == 80.):
                p_i = 0
            p_ie = np.zeros((4, 1))
            p_ie[l] = p_i
            p_ie[(l+1) % 4] = p_i
            lengthOfEdge = ((pnt_x[l] - pnt_x[(l + 1) % 4]) **
                            2 + (pnt_y[l] - pnt_y[(l + 1) % 4]) ** 2) ** 0.5
            normV = np.array([pnt_x[(l + 1) % 4] - pnt_x[l],
                             pnt_y[(l + 1) % 4] - pnt_y[l]])
            normV = normV / np.linalg.norm(normV)
            normV = np.array([normV[1], -normV[0]])
            for i in range(0, 6):
                xi_e = xi_i[i]
                eta_e = xi_i[i]

                if l == 0:
                    eta_e = -1
                elif l == 1:
                    xi_e = 1
                elif l == 2:
                    eta_e = 1
                elif l == 3:
                    xi_e = -1
                Psi = np.array(
                    [1 / 4.0 * (1 - xi_e) * (1 - eta_e),
                     1 / 4.0 * (1 + xi_e) * (1 - eta_e),
                     1 / 4.0 * (1 + xi_e) * (1 + eta_e),
                     1 / 4.0 * (1 - xi_e) * (1 + eta_e)],
                )
                Phi = np.array([
                    1 / 4. * xi_e * eta_e * (xi_e - 1) * (eta_e - 1),
                    1 / 2. * eta_e * (1 - xi_e ** 2) * (eta_e - 1),
                    1 / 4. * xi_e * eta_e * (xi_e + 1) * (eta_e - 1),
                    1 / 2. * xi_e * (1 - eta_e ** 2) * (xi_e - 1),
                    (1 - xi_e ** 2) * (1 - eta_e ** 2),
                    1 / 2. * xi_e * (1 - eta_e ** 2) * (xi_e + 1),
                    1 / 4. * xi_e * eta_e * (xi_e - 1) * (eta_e + 1),
                    1 / 2. * eta_e * (1 - xi_e ** 2) * (eta_e + 1),
                    1 / 4. * xi_e * eta_e * (xi_e + 1) * (eta_e + 1)
                ])
                F1 = F1 + w_i[i] * (np.transpose([Phi]) @
                                    [Psi] @ p_ie) * normV[0] * lengthOfEdge
                F2 = F2 + w_i[i] * (np.transpose([Phi]) @
                                    [Psi] @ p_ie) * normV[1] * lengthOfEdge

