import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    psi = np.loadtxt("Computation\data")

    # Paramètres physiques
    radius = 0.5

    # Maillage espace
    sizeR = 50
    sizeTheta = 75

    # Maillage en coordonnée polaire
    lengthR = 19.0 * radius
    stepSizeR = lengthR / (sizeR - 1.0)
    r = np.linspace(radius, radius + lengthR, sizeR)

    lengthTheta = 2.0 * np.pi
    stepSizeTheta = lengthTheta / (sizeTheta - 1.0)
    theta = np.linspace(0, lengthTheta, sizeTheta)

    # Maillage en coordonnée cartésienne
    x = np.zeros([sizeR, sizeTheta])
    y = np.zeros([sizeR, sizeTheta])
    for i in range(sizeR):
        for j in range(sizeTheta):
            x[i, j] = r[i] * np.cos(theta[j])
            y[i, j] = r[i] * np.sin(theta[j])

    # APPROXIMATION DES VITESSES PAR DIFFERENCES FINIES
    # CALCUL DE UR
    ur = np.zeros([sizeR, sizeTheta])
    for i in range(0, sizeR):
        for j in range(1, sizeTheta-1):
            ur[i, j] = (psi[i, j+1] - psi[i, j-1]) / (2*r[i]*stepSizeTheta)
        ur[i, 0] = (psi[i, 1] - psi[i, 0]) / (r[i]*stepSizeTheta)
        ur[i, sizeTheta-1] = ur[i, 0]

    # CALCUL DE UT
    ut = np.zeros([sizeR, sizeTheta])
    for j in range(0, sizeTheta):
        for i in range(1, sizeR-1):
            ut[i, j] = - (psi[i+1, j] - psi[i-1, j]) / (2*stepSizeR)
        ut[0, j] = - (psi[1, j] - psi[0, j]) / stepSizeR
        ut[sizeR-1, j] = - (psi[sizeR-1, j] - psi[sizeR-2, j]) / stepSizeR

    # PROJECTION DANS LA BASE CARTESIENNE
    ux = ur * np.cos(theta) - ut * np.sin(theta)
    uy = ur * np.sin(theta) + ut * np.cos(theta)

    # CALCUL DE LA NORME DE LA VITESSE
    u = np.sqrt(ux**2 + uy**2)

    fig, ax = plt.subplots()
    ax.pcolormesh(x, y, u, shading="gouraud")
    ax.axis("equal")
    plt.show()
