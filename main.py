import numpy as np
import matplotlib.pyplot as plt

psi = np.loadtxt("Computation\data")

# Paramètres physiques
radius = 0.5
Re = 200.0
uInf = 1.0

# Maillage espace
sizeR = 40
sizeTheta = 40

# Maillage en coordonnée polaire
lengthR = 19.0 * radius
stepSizeR = lengthR / (sizeR - 1.0)
r = np.linspace(radius, radius + lengthR, sizeR)

lengthTheta = 2.0 * np.pi
stepSizeTheta = lengthTheta / (sizeTheta - 1.0)
theta = np.linspace(0, lengthTheta, sizeTheta)


x = np.zeros([sizeR,sizeTheta])
y = np.zeros([sizeR,sizeTheta])
for i in range(0, sizeR):
    for j in range(0, sizeTheta):
        x[i,j] = r[i] * np.cos(theta[j])
        y[i,j] = r[i] * np.sin(theta[j])

# APPROXIMATION DES VITESSES PAR DIFFERENCES FINIES
# CALCUL DE UR
ur = np.zeros([sizeR,sizeTheta])
for i in range(1, sizeR-2):
    for j in range(1, sizeTheta-2):
        ur[i, j] = (psi[i, j+1]-psi[i, j-1]) / (2*r[i]*stepSizeTheta)
    ur[i, 0] = (psi[i, 1] - psi[i, 0]) / (r[i]*stepSizeTheta)
    ur[i, sizeR-1] = (psi[i, sizeR-1] - psi[i, sizeR-2]) / (r[i]*stepSizeTheta)


# CALCUL DE UT
ut = np.zeros([sizeR,sizeTheta])
for j in range(1, sizeTheta-2):
    for i in range(1, sizeR-2):
        ut[i, j] = - (psi[i+1, j] - psi[i-1, j]) / (2*stepSizeR)
    ut[0, j] = - (psi[1, j] - psi[0, j]) / stepSizeR
    ut[sizeTheta-1, j] = - (psi[sizeTheta-1, j] - psi[sizeTheta-2, j]) / stepSizeR
    
# PROJECTION DANS LA BASE CARTESIENNE
ux = np.zeros([sizeR,sizeTheta])
uy = np.zeros([sizeR,sizeTheta])
for i in range(0, sizeR-1):
    for j in range(0, sizeTheta-1):
        ux[i,j] = ur[i,j] * np.cos(theta[j]) - ut[i,j] * np.sin(theta[j])
        uy[i,j] = ur[i,j] * np.sin(theta[j]) + ut[i,j] * np.cos(theta[j])

#CALCUL DE LA NORME DE LA VITESSE
u = np.zeros([sizeR,sizeTheta])
for i in range(0, sizeR-1):
    for j in range(0, sizeTheta-1):
        u[i,j] = ux[i,j]**2 + uy[i,j]**2

fig, ax = plt.subplots()
ax.pcolormesh(x,y,u, shading="nearest")
ax.axis("equal")
plt.show()