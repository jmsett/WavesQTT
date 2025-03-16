# WavesQTT
Simulation of wave dynamics in the QTT framework

This repository provides a simulation framework for the time evolution of the Nonlinear Schrödinger Equation (NLSE) using the Split-Step Fourier Method (SSFM) combined with Quantized Tensor Train (QTT) techniques. The implementation leverages FFT for efficient computation of spatial derivatives and QTT compression to reduce the computational cost of high-dimensional simulations.

## Nonlinear Schrödinger Equation (NLSE) Simulation

This repository simulates the time evolution of the Nonlinear Schrödinger Equation (NLSE) using the Split-Step Fourier Method (SSFM) combined with Quantized Tensor Train (QTT) compression. The NLSE is given by:

$$ \Huge i \frac{\partial \psi}{\partial t} = -\frac{\partial^2 \psi}{\partial x^2} -2 |\psi|^2 \psi,$$

where $\psi(x, t)$ is the complex wave function, $x$ is the spatial variable, and $t$ is time.

### Split-Step Fourier Method (SSFM)

The equation is split into two parts:

1. **Linear operator** $\large L(\psi) = -\frac{\partial^2 \psi}{\partial x^2}$, which accounts for dispersion and is solved efficiently in the Fourier domain.

3. **Nonlinear operator** $\large N(\psi) = -2|\psi|^2 \psi$, which accounts for the self-interaction of the wave function.

### Time Evolution

The time evolution over a small time step $\Delta t$ is computed as follows:

1. Apply the nonlinear evolution: $\large\psi(x, t_0 + \Delta t) = \psi(x, t_0) e^{-i N(\psi) \Delta t}.$

2. Apply the linear evolution using the Fourier Transform: $\large \hat{\psi}(k, t + \Delta t) = \hat{\psi}(k, t_0) e^{-i k^2 \Delta t}$, where $\large\hat{\psi}(k, t)$ is the Fourier transform of $\large\psi(x, t)$.

3. Transform back to real space via the inverse Fourier transform to get the result in real space.

$$\Huge\psi(x, t_0 + \Delta t) = \mathcal{F}^{-1}\left[e^{-i k^2 \Delta t}\mathcal{F}\left[e^{i N(\psi) \Delta t} \psi(x, t_0)\right]\right] $$

This method alternates between linear and nonlinear steps, efficiently solving the NLSE.

Note that for the linear step for a $n$-component discretized function is given by:

$$\Huge
e^{-i k^2 \Delta t} \hat{\psi}(k, t_0)  =
\left[e^{-i\Delta t\left(\frac{1\times2\pi}{n\Delta x}\right)^2}, e^{-i\Delta t\left(\frac{2\times2\pi}{n\Delta x}\right)^2}, \dots, e^{-i\Delta t\left(\frac{n\times2\pi}{n\Delta x}\right)^2}\right]
\cdot
\begin{bmatrix}
\hat{\psi}_1\\
\hat{\psi}_2\\
\vdots\\
\hat{\psi}_n
\end{bmatrix}
,$$

which is equivalent to a multiplication by a diagonal matrix (or diagonal MPO in the QTT context).

### Features

- **FFT-Based Implementation**: Efficient handling of spatial derivatives using the Fast Fourier Transform (FFT).
- **QTT Compression**: Reduces computational complexity for high-dimensional or large-scale problems.
- **Customizability**: Easily modify initial conditions, domain size, and simulation parameters.
