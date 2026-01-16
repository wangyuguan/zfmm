# Laplace and Helmholtz Fast Multipole Method with Complex Coordinates

This repository contains research code for the paper:

**Fast Multipole Method with Complex Coordinates**  
arXiv:2509.05458: https://arxiv.org/abs/2509.05458

**Authors:** Tristan Goodwill, Leslie Greengard, Jeremy Hoskins, Manas Rachh, and Yuguan Wang

**Third-party dependency:** mwrap  
https://github.com/zgimbutas/mwrap

---

## Overview

This package provides Fast Multipole Method (FMM) implementations for Laplace and Helmholtz kernels in 2D and 3D, including **complex-coordinate** variants where source/target points may be complex-valued (e.g., contour deformations and complex shifts).

### Available MATLAB-callable functions

1. `zhfmm2d` — 2D Helmholtz kernel (complex coordinates)
2. `zlfmm2d` — 2D Laplace kernel (complex coordinates)
3. `zhfmm3d` — 3D Helmholtz kernel (complex coordinates)
4. `zlfmm3d` — 3D Laplace kernel (complex coordinates)

---

## Build (MATLAB MEX)

To build the Fortran code into MATLAB MEX functions, run:

- `matlab_linux/runmatlab.sh` on Linux
- `matlab_mac/runmatlab.sh` on macOS

It is recommended to place `mwrap` at the same directory level as this repository so that the relative path below works:

../../mwrap/mwrap

### Configure paths

Edit `runmatlab.sh` to point to your local MATLAB `mex` and your `gfortran` runtime library path, for example:

```bash
MATLAB_BIN="/usr/local/MATLAB/R2023a/bin/mex"
GFORTRAN_LIB_PATH="/usr/lib/gcc/x86_64-linux-gnu/"
```


## Usage

### zhfmm2d — 2D complex Helmholtz FMM

This routine evaluates the 2D Helmholtz potential with optional dipole contributions using a fast multipole method with complex coordinates.

For target points `x_i`, the computed field is
$$
u(x_i)
=
\sum_{j=1}^{n_s} q_j \, \frac{i}{4} \, H_0\!\left( z_k \, R(x_i, y_j) \right)
\;-\;
\sum_{j=1}^{n_s} d_j \, \frac{i}{4}
\left\langle v_j,\,
\nabla_y H_0\!\left( z_k \, R(x_i, y_j) \right)
\right\rangle .
$$

where

$$
R(x,y) = \sqrt{(x_1 - y_1)^2 + (x_2 - y_2)^2}.
$$

---

#### MATLAB call

```matlab
[pot, grad] = zhfmm2d(eps, zk, ns, zsrc, ifcharge, charge, ...
                      ifdipole, dipstr, dipvec, nt, ztarg, ifpgh, isep);
```

#### Input arguments

- **eps** *(real scalar)*  
  Requested FMM accuracy.

- **zk** *(complex scalar)*  
  Helmholtz wavenumber.

- **ns** *(integer)*  
  Number of source points.

- **zsrc** *(2 × ns complex array)*  
  Source locations. Complex coordinates are supported.

- **ifcharge** *(0 or 1)*  
  Flag indicating whether charge sources are included.

- **charge** *(ns complex array)*  
  Charge strengths.

- **ifdipole** *(0 or 1)*  
  Flag indicating whether dipole sources are included.

- **dipstr** *(ns complex array)*  
  Dipole strengths.

- **dipvec** *(2 × ns complex array)*  
  Dipole orientation vectors.

- **nt** *(integer)*  
  Number of target points.

- **ztarg** *(2 × nt complex array)*  
  Target locations. Complex coordinates are supported.

- **ifpgh** *(integer)*  
  Output selector:  
  - `1` — compute potential only  
  - `2` — compute potential and gradient

- **isep** *(integer)*  
  Well-separation parameter used in FMM tree construction.

---

## Output arguments

- **pot** *(nt complex array)*  
  Potential evaluated at target points.

- **grad** *(2 × nt complex array)*  
  Gradient of the potential at target points *(only valid if `ifpgh = 2`)*.
