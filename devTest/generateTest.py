import sys

import numpy as np
import matplotlib.pyplot as plt

filename = "low_resolution.dat"

image = np.loadtxt(filename)

# Image is 2D array with format x y intensity
# x and y are the coordinates of the pixel (float)
# intensity is the intensity of the pixel (float)

# Determine resolution of image
# Unique x values
x  = image[:,0]
xu = np.unique(x)
Nx = len(xu)
# Unique y values
y  = image[:,1]
yu = np.unique(y)
Ny = len(yu)

dx = np.diff(xu)
# Check if dx is constant
if np.allclose(dx, dx[0]):
    dx = dx[0]
else:
    print("dx is not constant")
    sys.exit()

dy = np.diff(yu)
# Check if dy is constant
if np.allclose(dy, dy[0]):
    dy = dy[0]
else:
    print("dy is not constant")
    sys.exit()

Lx = dx*Nx
Ly = dy*Ny

print("Resolution of image: ", Nx, "x", Ny, "pixels")
print("dx = ", dx, "dy = ", dy, "Lx = ", Lx, "Ly = ", Ly)
print("Number of pixels: ", len(image),Nx*Ny)
print("x values: ", xu)
print("y values: ", yu)

# Convert image to 2D array, (Nx, Ny)
image2D_aux = image.copy()

# Convert x (y) values to indices in range 0 to Nx-1 (Ny-1)
xmin = np.min(xu)
xmax = np.max(xu)
ymin = np.min(yu)
ymax = np.max(yu)
for i in range(len(image2D_aux)):
    x = image2D_aux[i,0]
    y = image2D_aux[i,1]

    x_index = int((x - xmin)/dx)
    y_index = int((y - ymin)/dy)

    image2D_aux[i,0] = x_index
    image2D_aux[i,1] = y_index

image2D = np.zeros((Nx,Ny))

for i in range(len(image2D_aux)):
    x = int(image2D_aux[i,0])
    y = int(image2D_aux[i,1])
    image2D[x,y] = image2D_aux[i,2]

plt.imshow(image2D, origin='lower', extent=[-Lx/2, Lx/2, -Ly/2, Ly/2])
plt.show()

