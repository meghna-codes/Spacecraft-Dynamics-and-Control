# Spacecraft Dynamics and Control

## Overview

This repository contains simulation frameworks, control strategies, and analysis tools for spacecraft motion in complex dynamical environments. The primary focus is on trajectory design, station-keeping, and feedback control in nonlinear astrodynamical systems, with emphasis on the Circular Restricted Three-Body Problem (CR3BP).

The codebase supports research-oriented experimentation and is suitable for applications in mission design, orbital stability analysis, and advanced control implementation.

---

## Key Features

* Nonlinear spacecraft dynamics modeling (CR3BP framework)
* Linearized dynamics and variational equations
* Periodic orbit generation (e.g., Halo orbits)
* Control strategies:

  * Target Point Approach (TPA)
  * Time-Invariant LQR (TI-LQR)
  * Time-Varying LQR (TV-LQR)
* Station-keeping simulations and performance evaluation
* Numerical integration and trajectory propagation tools

---

## Getting Started

### Prerequisites

* MATLAB (recommended) or Python (NumPy, SciPy, Matplotlib)
* Basic understanding of orbital mechanics and control systems

### Running Simulations

1. Clone the repository:

   ```
   git clone https://github.com/your-username/your-repo-name.git
   ```

2. Navigate to the project directory:

   ```
   cd your-repo-name
   ```

3. Run main scripts (example in MATLAB):

   ```
   main_simulation.m
   ```

---

## Control Strategies Implemented

### Target Point Approach (TPA)

* Impulsive correction-based method
* Suitable for discrete station-keeping
* Lower computational cost

### Time-Invariant LQR (TI-LQR)

* Constant gain feedback control
* Based on linearized dynamics at equilibrium
* Moderate performance and simplicity

### Time-Varying LQR (TV-LQR)

* Gain scheduling along reference trajectory
* High-precision tracking capability
* Computationally intensive but optimal

---

## Results and Analysis

The repository includes tools to evaluate:

* Position and velocity tracking error
* RMS and maximum deviation
* Control effort (ΔV consumption)
* Stability characteristics of closed-loop systems

Comparative studies between TPA, TI-LQR, and TV-LQR are provided for station-keeping performance.

---

## Applications

* Lagrange point mission design (L1/L2)
* Halo orbit maintenance
* Autonomous spacecraft navigation
* Optimal control in nonlinear dynamical systems

---

## Future Work

* Robust and adaptive control under uncertainty
* Model predictive control (MPC)
* Machine learning-assisted control policies
* Extension to multi-body and perturbed environments

---

## License

This project is licensed under the MIT License. See the LICENSE file for details.

---

## Citation

If you use this work in your research, please cite:

```
@misc{MeghnaSantra2026,
  title={Spacecraft Dynamics and Control},
  author={Meghna Santra},
  year={2026},
  note={GitHub repository}
}
```

---

## Contact

For questions, collaborations, or suggestions, feel free to reach out.

---
