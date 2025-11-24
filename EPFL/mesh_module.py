import numpy as np
import matplotlib.pyplot as plt

def map_region(X, offset, parameters):
    rho = X[:, 0]
    theta = X[:, 1]

    a = parameters[0]
    b = 0
    d = (parameters[0] - parameters[1]) * (16) / (np.pi ** 3)
    c = -(3 / 4) * np.pi * d

    limit_rho = (
        a
        + b * (theta - np.pi / 4)
        + c * (theta - np.pi / 4) ** 2
        + d * (theta - np.pi / 4) ** 3
    )

    trho = (
        (1 / np.sin(theta)) * (rho)
        + (limit_rho) * (1 - rho)
    )

    tX = np.zeros(shape=X.shape)
    tX[:, 0] = trho * np.cos(theta + offset)
    tX[:, 1] = trho * np.sin(theta + offset)

    return tX

def create_subdomain_mesh(n=[10, 10], degree=1, p0=[0.0, 0.0], p1=[1.0, 1.0], parameters=np.array([0.5, 0.5, 0.5, 0.5])):
    conn = []

    for n_y in range(4 * n[1]):
        for n_x in range(n[0]):
            element_nodes = np.array([
                n_x * degree + n_y * degree * (n[0] * degree + 1),
                (n_x + 1) * degree + n_y * degree * (n[0] * degree + 1),
                n_x * degree + ((n_y + 1) % (4 * n[1])) * degree * (n[0] * degree + 1),
                (n_x + 1) * degree + ((n_y + 1) % (4 * n[1])) * degree * (n[0] * degree + 1)
            ])

            for m_x in range(1, degree):
                element_nodes = np.append(element_nodes, n_x * degree + m_x + n_y * degree * (n[0] * degree + 1))

            for m_y in range(1, degree):
                element_nodes = np.append(element_nodes, n_x * degree + (n_y * degree + m_y) * (n[0] * degree + 1))

            for m_y in range(1, degree):
                element_nodes = np.append(element_nodes, (n_x + 1) * degree + (n_y * degree + m_y) * (n[0] * degree + 1))

            for m_x in range(1, degree):
                element_nodes = np.append(element_nodes, n_x * degree + m_x + ((n_y + 1) % (4 * n[1])) * degree * (n[0] * degree + 1))

            for m_y in range(1, degree):
                for m_x in range(1, degree):
                    element_nodes = np.append(element_nodes, n_x * degree + m_x + (n_y * degree + m_y) * (n[0] * degree + 1))

            conn.append(element_nodes)

    x = np.linspace(0, 1, n[0] * degree + 1)
    y = np.linspace(np.pi / 4, 3 * np.pi / 4, n[1] * degree + 1)[:-1]
    X, Y = np.meshgrid(x, y)

    points = np.column_stack((X.ravel(), Y.ravel()))

    coords1 = map_region(points, 0, [parameters[3], parameters[2]])
    coords2 = map_region(points, 0.5 * np.pi, [parameters[2], parameters[0]])
    coords3 = map_region(points, np.pi, [parameters[0], parameters[1]])
    coords4 = map_region(points, 1.5 * np.pi, [parameters[1], parameters[3]])

    coords = np.vstack((coords1, coords2, coords3, coords4))

    coords[:, 0] = 0.5 * (p1[0] + p0[0]) + 0.5 * (p1[0] - p0[0]) * coords[:, 0]
    coords[:, 1] = 0.5 * (p1[1] + p0[1]) + 0.5 * (p1[1] - p0[1]) * coords[:, 1]

    return coords, conn

def plot_mesh(coords, conn):
    fig, ax = plt.subplots(figsize=(8, 8))
    coords = np.array(coords)
    for element in conn:
        element = np.array(element, dtype=int)
        poly = element[[0, 1, 3, 2, 0]]
        ax.plot(coords[poly, 0], coords[poly, 1], 'k-')
    ax.set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Mesh Visualization')
    plt.show()
