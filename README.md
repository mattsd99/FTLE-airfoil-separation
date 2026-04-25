# FTLE-airfoil-separation

The **Finite-Time Lyapunov Exponent (FTLE)** measures the maximum rate of stretching between neighbouring fluid particles over a finite integration time *T*. It quantifies how rapidly nearby trajectories diverge in a flow, revealing hidden material surfaces called **Lagrangian Coherent Structures (LCS)** that act as transport barriers in the flow.

## Data
The full dataset (31 velocity field snapshots, NACA 0012 α=12°, Re=5×10⁵):

[![Google Drive](https://img.shields.io/badge/Google%20Drive-Download%20Dataset-blue?logo=googledrive)](https://drive.google.com/drive/folders/1mZT73kUixwulYAlUfqTtegU_ue_HAqiE?usp=sharing)

Extract all CSV files into the same folder as the MATLAB scripts and set 'DATA_DIR' in 'run_ftle_backwards.m' to that folder path.


### Formula

The FTLE at a point **x₀** is defined as:

σ(x₀, t₀, T) = (1 / 2|T|) · ln( λ_max(Δ) )

where **Δ** is the **right Cauchy-Green deformation tensor**:

Δ = (dϕ/dx₀)ᵀ · (dϕ/dx₀)

and **ϕ** is the **flow map**, the final position of a particle seeded at **x₀** at time *t₀* after being advected for time *T*.

**λ_max** is the largest eigenvalue of **Δ**. The flow map Jacobian **dϕ/dx₀** is the 2×2 matrix:

F = | ∂xf/∂x₀   ∂xf/∂y₀ |
    | ∂yf/∂x₀   ∂yf/∂y₀ |

Results
Backward FTLE — 6 time instances
The 2×3 comparison shows the attracting LCS evolving through one shedding cycle. The bright ridges mark the boundary of the separated shear layer, with clear time dependency indication of shear layer rollup location and occurrence with time.
Animated GIF
The GIF cycles through t₀ = 0 to 1.0 s, showing the evolution of flow separation and inception of the leading edge vortex with respect to time and location.


Dependencies

MATLAB R2017a or later
No additional toolboxes required


Attribution
FTLE computation framework adapted from:

Shadden, S.C., Lekien, F., Marsden, J.E. (2005). Definition and properties of Lagrangian coherent structures from finite-time Lyapunov exponents in two-dimensional aperiodic flows. Physica D, 212(3-4), 271–304.

Extended and generalised for unstructured CSV velocity data, analytical NACA airfoil masking, phase-resolved visualisation, and animated output.

Author
Matt Duran — PhD, Mechanical Engineering, University of Central Florida
