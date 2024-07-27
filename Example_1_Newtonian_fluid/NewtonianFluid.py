import numpy as np
import matplotlib.pyplot as plt
import h5py as h5py
import scipy.sparse as sp
from scipy.sparse import csr_matrix, csc_matrix, lil_matrix, coo_matrix
import matplotlib.tri as tri

with h5py.File("Structured_Mesh.h5", "r") as f:
    Points = np.array(f["Points"][:, 0:2])
    Elements = np.array(f["Elements"][:])
    Points_h = np.array(f["Points_h"][:, 0:2])
    Elements_h = np.array(f["Elements_h"][:])

# ____________ parameters ____________-
mu = 1.0e3

# ________ solver _____________

row_indices_A = []
col_indices_A = []
values_A = []


def add_triplets_A(triplets):
    for (i, j, value) in triplets:
        row_indices_A.append(i)
        col_indices_A.append(j)
        values_A.append(value)


row_indices_b = []
col_indices_b = []
values_b = []


def add_triplets_b(triplets):
    for (i, j, value) in triplets:
        row_indices_b.append(i)
        col_indices_b.append(j)
        values_b.append(value)


NumPnts = np.shape(Points)[0]
NumPnts_h = np.shape(Points_h)[0]

Dims_K = NumPnts_h + NumPnts_h + NumPnts

rows = np.array([], dtype=int)
cols = np.array([], dtype=int)
data = np.array([])
K = coo_matrix((data, (rows, cols)), shape=(Dims_K, Dims_K))
b = coo_matrix((data, (rows, cols)), shape=(Dims_K, 1))

NumElements = np.shape(Elements)[0]
w_i = np.array([0.17132, 0.36076, 0.46791, 0.17132, 0.36076, 0.46791])
xi_i = np.array([0.93247, 0.66121, 0.23862, -0.93247, -0.66121, -0.23862])

Pnt_Top = []
Pnt_Bot = []

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

    pnt_x_h = Points_h[Elements_h[e, :], 0]
    pnt_y_h = Points_h[Elements_h[e, :], 1]

    triplet_k = []

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
            x_xi = np.dot(Phi_xi, pnt_x_h)
            y_xi = np.dot(Phi_xi, pnt_y_h)
            x_eta = np.dot(Phi_eta, pnt_x_h)
            y_eta = np.dot(Phi_eta, pnt_x_h)
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
                Phi_x) @ Phi_x + np.transpose(Phi_y) @ Phi_y) * deter_J
            D12 = D12 + mu * w_i[i] * w_i[j] * (np.transpose(
                Phi_x) @ Phi_y) * deter_J
            D21 = D21 + mu * w_i[i] * w_i[j] * (np.transpose(
                Phi_y) @ Phi_x) * deter_J
            D22 = D22 + mu * w_i[i] * w_i[j] * (2 * np.transpose(
                Phi_y) @ Phi_y + np.transpose(Phi_x) @ Phi_x) * deter_J

            C1 = C1 + mu * w_i[i] * w_i[j] * \
                (np.transpose(Phi_x) @ [Psi]) * deter_J
            C2 = C2 + mu * w_i[i] * w_i[j] * \
                (np.transpose(Phi_y) @ [Psi]) * deter_J

    Element_l_8 = np.array([Elements_h[e, 0], Elements_h[e, 1], Elements_h[e, 2],
                            Elements_h[e, 5], Elements_h[e, 8],
                            Elements_h[e, 7], Elements_h[e, 6], Elements_h[e, 3]])
    for l in range(0, 4):
        # _________ boundary condition: pressure ______________
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
                # print(np.shape(np.transpose([Phi])), np.shape([Psi]), np.shape(p_ie))
                F1 = F1 + w_i[i] * (np.transpose([Phi]) @
                                    [Psi] @ p_ie) * normV[0] * lengthOfEdge
                # print(np.shape(np.transpose([Phi])), np.shape([Psi]), np.shape(p_ie))
                F2 = F2 + w_i[i] * (np.transpose([Phi]) @
                                    [Psi] @ p_ie) * normV[1] * lengthOfEdge
        if (pnt_y[l] == 60 and pnt_y[(l + 1) % 4] == 60):
            Pnt_Top.append(Element_l_8[l * 2])
            Pnt_Top.append(Element_l_8[(l * 2 + 1)])
            Pnt_Top.append(Element_l_8[(l * 2 + 2) % 8])
        if (pnt_y[l] == 0 and pnt_y[(l + 1) % 4] == 0):
            Pnt_Bot.append(Element_l_8[l * 2])
            Pnt_Bot.append(Element_l_8[(l * 2 + 1)])
            Pnt_Bot.append(Element_l_8[(l * 2 + 2) % 8])

    # ________________assemble ______________
    for i in range(9):
        for j in range(9):
            triplet_1 = [(Elements_h[e, i], Elements_h[e, j], D11[i, j]),
                         (Elements_h[e, i], Elements_h[e, j] +
                          NumPnts_h, D12[i, j]),
                         (Elements_h[e, i] + NumPnts_h,
                          Elements_h[e, j], D21[i, j]),
                         (Elements_h[e, i] + NumPnts_h, Elements_h[e, j] + NumPnts_h, D22[i, j])]
            add_triplets_A(triplet_1)

    for i in range(9):
        for j in range(4):
            add_triplets_A(
                [(Elements_h[e, i], Elements[e, j] + 2 * NumPnts_h, -C1[i, j])])
            add_triplets_A(
                [(Elements_h[e, i] + NumPnts_h, Elements[e, j] + 2 * NumPnts_h, -C2[i, j])])

    for i in range(4):
        for j in range(9):
            add_triplets_A(
                [(Elements[e, i] + 2 * NumPnts_h, Elements_h[e, j], B1[i, j])])
            add_triplets_A(
                [(Elements[e, i] + 2 * NumPnts_h, Elements_h[e, j] + NumPnts_h, B2[i, j])])

    for i in range(9):
        add_triplets_b([(Elements_h[e, i], 0, -F1[i, 0])])
        add_triplets_b([(Elements_h[e, i] + NumPnts_h, 0, -F2[i, 0])])

    # for i in range(len(row_indices_A)):
        # if row_indices_A[i] >= Dims_K or col_indices_A[i] >= Dims_K:
        #    print(row_indices_A[i], col_indices_A[i], values_A[i])
       # print(row_indices_A[i], col_indices_A[i], values_A[i])
    K_s = coo_matrix((values_A, (row_indices_A, col_indices_A)),
                     shape=(Dims_K, Dims_K))
    b_s = coo_matrix(
        (values_b, (row_indices_b, col_indices_b)), shape=(Dims_K, 1))
    K = K + K_s
    b = b + b_s
    row_indices_A = []
    col_indices_A = []
    values_A = []
    row_indices_b = []
    col_indices_b = []
    values_b = []

K = K.tolil()
b = b.tolil()
# ___________veolocity boundary
for i in range(len(Pnt_Top)):
    PointID = int(Pnt_Top[i])
    b = b - K[:, PointID] * 1e-2
    K[PointID, :] = 0
    K[:, PointID] = 0
    K[PointID, PointID] = 1
    b[PointID, 0] = 1e-2

    b = b - K[:, PointID + NumPnts_h] * 0
    K[PointID + NumPnts_h, :] = 0
    K[:, PointID + NumPnts_h] = 0
    K[PointID + NumPnts_h, PointID + NumPnts_h] = 1
    b[NumPnts_h + PointID] = 0

for i in range(len(Pnt_Bot)):
    PointID = int(Pnt_Top[i])
    b = b - K[:, PointID] * 0
    K[PointID, :] = 0
    K[:, PointID] = 0
    K[PointID, PointID] = 1
    b[PointID, 0] = 0

    b = b - K[:, PointID + NumPnts_h] * 0
    K[PointID + NumPnts_h, :] = 0
    K[:, PointID + NumPnts_h] = 0
    K[PointID + NumPnts_h, PointID + NumPnts_h] = 1
    b[NumPnts_h + PointID] = 0

K = K.tocsr()
b = b.tocsr()

x = sp.linalg.spsolve(K, b)

# print(x)

fig,ax=plt.subplots(1, 2, figsize=(10,4))

def plot_fem_mesh(nodes_x, nodes_y, elements, f):
    for element in elements:
        x = [nodes_x[element[i]] for i in range(len(element))]
        y = [nodes_y[element[i]] for i in range(len(element))]
        f.fill(x, y, edgecolor='black', fill=False)
    
plot_fem_mesh(Points[:, 0], Points[:, 1], Elements,ax[0])
ax[0].quiver(Points_h[:, 0], Points_h[:, 1], x[0:NumPnts_h], x[NumPnts_h:2*NumPnts_h])

plot_fem_mesh(Points[:, 0], Points[:, 1], Elements,ax[1])
for i in range(NumPnts):
    ax[1].text(Points[i, 0], Points[i, 1], f"{x[2*NumPnts + i]:.2e}")

plt.show()