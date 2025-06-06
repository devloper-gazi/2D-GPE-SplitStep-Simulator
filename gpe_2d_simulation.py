import numpy as np
import matplotlib.pyplot as plt

# --------------------------------------
# 2D Gross–Pitaevskii Simulation
# --------------------------------------

# 1) Parameters and Grid Definition
Nx = 256             # Number of grid points in x
Ny = 256             # Number of grid points in y
Lx = 10.0            # Physical length in x: [-Lx/2, Lx/2]
Ly = 10.0            # Physical length in y: [-Ly/2, Ly/2]
dx = Lx / Nx         # Spatial step in x
dy = Ly / Ny         # Spatial step in y

x = np.linspace(-Lx/2, Lx/2, Nx)
y = np.linspace(-Ly/2, Ly/2, Ny)
X, Y = np.meshgrid(x, y, indexing='ij')

dt = 0.0002          # Time step
g  = 1.0             # Interaction strength
V0 = 10.0            # Potential depth
a  = 2.0             # Potential period
steps = 2000         # Number of time steps

# 2) 2D Periodic Potential (Photonic-Crystalline)
#    V(x,y) = V0 * [cos^2(pi x / a) + cos^2(pi y / a)]
Vx = np.cos(np.pi * X / a)**2
Vy = np.cos(np.pi * Y / a)**2
V = V0 * (Vx + Vy)

# 3) Initial Wavefunction: 2D Gaussian Profile
sigma = 1.0
psi = np.exp(-(X**2 + Y**2) / (2 * sigma**2)).astype(np.complex128)
psi /= np.sqrt(np.sum(np.abs(psi)**2) * dx * dy)

# 4) Fourier-Space Wavenumbers
kx = 2 * np.pi * np.fft.fftfreq(Nx, d=dx)
ky = 2 * np.pi * np.fft.fftfreq(Ny, d=dy)
KX, KY = np.meshgrid(kx, ky, indexing='ij')
K2 = KX**2 + KY**2

# 5) Time Evolution: Split-Step Fourier (2D)
for _ in range(steps):
    # 5a) Potential + Nonlinear Half Step
    psi *= np.exp(-1j * (V + g * np.abs(psi)**2) * dt / 2)
    
    # 5b) Kinetic Half Step in Fourier Space
    psi_k = np.fft.fft2(psi)
    psi_k *= np.exp(-1j * K2 * dt / 2)
    psi = np.fft.ifft2(psi_k)
    
    # 5c) Potential + Nonlinear Half Step
    psi *= np.exp(-1j * (V + g * np.abs(psi)**2) * dt / 2)
    
    # 5d) Normalise
    norm = np.sqrt(np.sum(np.abs(psi)**2) * dx * dy)
    psi /= norm

# 6) Compute and Visualise Density
density = np.abs(psi)**2

plt.figure(figsize=(6, 5))
plt.imshow(
    density.T, 
    extent=[-Lx/2, Lx/2, -Ly/2, Ly/2], 
    origin='lower', 
    cmap='inferno'
)
plt.colorbar(label='Density $|\psi(x,y)|^2$')
plt.title("2D GPE: Density Modulation in a Periodic Potential")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.show()
