# 2D-GPE-Simulation 🚀

## Overview 🌟

This repository contains a Python implementation of the two-dimensional (2D) Gross–Pitaevskii Equation (GPE) solver using the Split-Step Fourier Method. It simulates the time evolution of a polariton condensate under a photonic-crystal-like periodic potential and visualises the resulting density modulation in two dimensions.

---

## Contents 📂

- **`gpe_2d_simulation.py`**  
  Main Python script that sets up the grid, initial wavefunction, periodic potential, and performs time evolution.

- **`README.md`**  
  This document, which explains the mathematical model, numerical method, usage instructions, and file structure.

- **`LICENSE`**  
  MIT License file.

---

## 1. Mathematical Model 📐

We solve the conservative 2D GPE (ℏ = 1, m* = 1) for a condensate wavefunction \(\psi(x,y,t)\):

<p>
  \( i \frac{\partial \psi(x,y,t)}{\partial t}
  = \left[ -\tfrac{1}{2}\bigl(\partial_x^2 + \partial_y^2\bigr)
  + V(x,y) + g\,|\psi(x,y,t)|^2 \right] \psi(x,y,t) \,.\)
</p>

\[
i \frac{\partial \psi(x,y,t)}{\partial t}
= \left[ -\tfrac{1}{2}\bigl(\partial_x^2 + \partial_y^2\bigr)
+ V(x,y) + g\,|\psi(x,y,t)|^2 \right] \psi(x,y,t)\,.
\]

- \(\psi(x,y,t)\): complex-valued condensate wavefunction  
- \(g\): interaction strength (nonlinear coefficient)  
- \(V(x,y)\): external periodic potential  

### 1.1 Periodic Potential 📏

We choose a simple separable 2D potential that mimics a photonic crystal:

\[
V(x,y) \;=\; V_0 \Bigl[\cos^2\!\bigl(\tfrac{\pi x}{a}\bigr) 
+ \cos^2\!\bigl(\tfrac{\pi y}{a}\bigr)\Bigr]\,,
\]

where:  
- \(V_0\) is the potential depth (e.g., 10)  
- \(a\) is the lattice period (e.g., 2)

---

## 2. Numerical Method: Split-Step Fourier ⚛️

We discretise the square domain \([-L_x/2,\,L_x/2]\times[-L_y/2,\,L_y/2]\) into \(N_x\times N_y\) grid points. Time step is \(\Delta t\). Each full step is performed by Strang splitting into three sub-steps:

1. **Potential + Nonlinear Half Step**  
   \[
   \psi(x,y,t) 
   \;\longleftarrow\; 
   \psi(x,y,t)\,\exp\!\Bigl[-\,i\bigl(V(x,y)+g\,|\psi|^2\bigr)\,\tfrac{\Delta t}{2}\Bigr].
   \]

2. **Kinetic Half Step (2D Fourier Space)**  
   - Compute \(\Psi(k_x,k_y) = \mathcal{F}_{2D}[\psi(x,y,t)]\).  
   - Multiply by \(\exp\bigl[-\,i\,(k_x^2 + k_y^2)\,\tfrac{\Delta t}{2}\bigr]\).  
   - Inverse 2D Fourier transform back to real space.  

   Mathematically:  
   \[
   \Psi(k_x,k_y) \;=\; \mathcal{F}_{2D}[\psi(x,y,t)],\quad
   \Psi(k_x,k_y)\;\longleftarrow\;\Psi(k_x,k_y)\,\exp\!\Bigl[-\,i\,(k_x^2 + k_y^2)\,\tfrac{\Delta t}{2}\Bigr],\quad
   \psi(x,y)\;=\;\mathcal{F}_{2D}^{-1}[\Psi(k_x,k_y)].
   \]

3. **Potential + Nonlinear Half Step (again)**  
   \[
   \psi(x,y,t + \Delta t) 
   \;\longleftarrow\; 
   \psi(x,y,t)\,\exp\!\Bigl[-\,i\bigl(V(x,y)+g\,|\psi|^2\bigr)\,\tfrac{\Delta t}{2}\Bigr].
   \]

4. **Normalisation**  
   \[
   \psi \;\longleftarrow\; 
   \frac{\psi}{\sqrt{\sum_{i,j} |\psi_{i,j}|^2\,\Delta x\,\Delta y}}\,,
   \]
   to preserve the total norm.

---

## 3. File Structure 🗂️

2D-GPE-Simulation/
│
├── gpe_2d_simulation.py # Main simulation script
├── README.md # This document
└── LICENSE # MIT License file

---

## 4. Usage Instructions 🛠️

1. **Clone the repository**  
   ```bash
   git clone https://github.com/<your-username>/2D-GPE-Simulation.git
   cd 2D-GPE-Simulation
2. (Optional) Create and activate a virtual environment ⚙️
This helps avoid conflicts if pip has issues installing globally.
```bash
● python3 -m venv venv
# On Windows:
● venv\Scripts\activate
# On macOS/Linux:
● source venv/bin/activate
```
3. Install dependencies 📦
Ensure you have Python 3 (≥ 3.7) and the required packages:
```bash
● pip install numpy matplotlib
```
⚠️ If you encounter permission errors, try:
```bash
● pip install --user numpy matplotlib
```
4. Run the simulation ✨
```bash
python gpe_2d_simulation.py
```
5. Deactivate the virtual environment (if used) ✓
```bash
deactivate
```
## 5. Code Excerpt 📋
Below is a minimal excerpt demonstrating the main loop structure. For the full code, see gpe_2d_simulation.py.
```bash
import numpy as np
import matplotlib.pyplot as plt

# 1) Parameters and Grid Definition
Nx, Ny = 256, 256
Lx, Ly = 10.0, 10.0
dx, dy = Lx/Nx, Ly/Ny
x = np.linspace(-Lx/2, Lx/2, Nx)
y = np.linspace(-Ly/2, Ly/2, Ny)
X, Y = np.meshgrid(x, y, indexing='ij')

dt = 0.0002
g = 1.0
V0 = 10.0
a = 2.0
steps = 2000

# 2) 2D Periodic Potential (Photonic Crystal–Like)
V = V0 * (np.cos(np.pi*X/a)**2 + np.cos(np.pi*Y/a)**2)

# 3) Initial Wavefunction (2D Gaussian) + Normalisation
sigma = 1.0
psi = np.exp(-(X**2 + Y**2)/(2*sigma**2)).astype(np.complex128)
psi /= np.sqrt(np.sum(np.abs(psi)**2) * dx * dy)

# 4) Fourier-Space Wavenumbers
kx = 2*np.pi * np.fft.fftfreq(Nx, d=dx)
ky = 2*np.pi * np.fft.fftfreq(Ny, d=dy)
KX, KY = np.meshgrid(kx, ky, indexing='ij')
K2 = KX**2 + KY**2

# 5) Time Evolution: Split-Step Fourier (2D)
for _ in range(steps):
    # 5a) Potential + Nonlinear Half Step
    psi *= np.exp(-1j * (V + g*np.abs(psi)**2) * dt/2)
    
    # 5b) Kinetic Half Step (2D Fourier Space)
    psi_k = np.fft.fft2(psi)
    psi_k *= np.exp(-1j * K2 * dt/2)
    psi = np.fft.ifft2(psi_k)
    
    # 5c) Potential + Nonlinear Half Step
    psi *= np.exp(-1j * (V + g*np.abs(psi)**2) * dt/2)
    
    # 5d) Normalisation
    norm = np.sqrt(np.sum(np.abs(psi)**2) * dx * dy)
    psi /= norm

# 6) Compute and Visualise Density
density = np.abs(psi)**2
plt.figure(figsize=(6,5))
plt.imshow(density.T, extent=[-Lx/2, Lx/2, -Ly/2, Ly/2], origin='lower', cmap='inferno')
plt.colorbar(label='Density $|\psi(x,y)|^2$')
plt.title("2D GPE: Density Modulation in a Periodic Potential")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.show()
```

## 6. Default Parameters 🔧

Below are the default values used throughout the simulation:

- **Spatial grid**  
  - Number of points: `N_x = N_y = 256`  
  - Physical extent: `L_x = L_y = 10.0` (so that \(\Delta x = \Delta y = 10.0 / 256\))  

- **Time step**  
  - \(\Delta t = 0.0002\)  

- **Interaction coefficient**  
  - \(g = 1.0\)  

- **Periodic potential**  
  - Depth: \(V_0 = 10.0\)  
  - Period: \(a = 2.0\)  

- **Number of time‐evolution steps**  
  - `steps = 2000`  

These values can be adjusted in `gpe_2d_simulation.py` to explore different regimes or to optimise performance and accuracy.  
