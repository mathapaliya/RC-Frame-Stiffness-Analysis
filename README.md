Dynamic Modal Analysis of 6-Storey RC Frame
This repository contains a Python-based implementation for performing dynamic modal analysis on a 6-storey reinforced concrete (RC) frame using shear-type behavior assumptions. The script computes natural frequencies, mode shapes, modal participation, modal displacements, inter-storey drifts, and equivalent static forces using both SRSS and CQC methods. All results are visualized using matplotlib.

🧮 Model Summary
Number of storeys: 6

Floor height (uniform): 3.2 m

System: Shear-type multi-storey RC frame (no rotational DOFs)

Degrees of freedom: 6 translational DOFs (1 per floor level)

Modal analysis type: Linear undamped free vibration (eigenvalue problem)

Mode combination methods: SRSS and CQC

🧾 Mass and Stiffness Input
The mass matrix (M) and stiffness matrix (K) have been pre-computed externally (based on building geometry, material properties, tributary areas, and member dimensions).

These matrices are directly inserted into the script.

Units:

Mass: kg (implicitly; via kN/g conversion if applicable)

Stiffness: N/m

🧠 Features of the Analysis
Solves generalized eigenvalue problem:
Kϕ = λMϕ

Extracts:

Natural circular frequencies ωₙ

Time periods Tₙ

Mode shapes φₙ

Modal participation factors τₙ

Effective modal mass ratios (EMMR)

Calculates:

Modal displacements (uₙ)

Inter-storey drifts

Equivalent static force distributions

Storey shears

Bending moments per frame

Combines modal results using:

SRSS (Square Root of Sum of Squares)

CQC (Complete Quadratic Combination)

Visualizes:

Mode shape displacements

Inter-storey drifts

Storey shear force diagrams

Combined displacement, drift, and force diagrams (SRSS vs. CQC)

🖥️ Dependencies
This code requires the following Python libraries:

numpy

scipy

matplotlib



📊 Output Summary
Natural frequencies (ωₙ) and time periods (Tₙ)

Mode shapes normalized to max displacement

Modal displacements (in mm)

Inter-storey drift diagrams

Modal participation factors (τₙ)

Equivalent static force per mode and per frame

Combined results using SRSS and CQC for:

Displacement

Drift

Shear

Force per frame

Bending moments

📘 Notes
The script assumes no torsional modes (pure shear type lateral displacement only).

The height of each storey is 3.2 meters.

You may edit the mass or stiffness matrices to model other building configurations.

Visualization blocks (e.g., for plotting) are organized by response type.
