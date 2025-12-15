# Wall Hindrance Mobility Analysis (BEM)

This project simulates the hydrodynamic mobility of a rigid particle near a plane wall under Stokes flow conditions using a Fortran Boundary Element Method (BEM) solver (`prtcl_3d_mob`) and Python scripts for automation and visualization.

The experiment systematically measures the particle's translational and angular velocities as a function of the **Surface Gap ($\delta$)** between the particle surface and a fixed, no-slip wall.

## 1. Prerequisites

### A. Software Requirements

You must have the following software installed:

1.  **Fortran Compiler:** `gfortran` (or similar) to compile the BEM solver.
2.  **Python 3:** (3.6 or newer).
3.  **Python Libraries:** `numpy`, `pandas`, and `matplotlib`.

```bash
# Install Python libraries
pip install numpy pandas matplotlib
````

### B. File Structure

Ensure all necessary files are in the main project directory:

```
├── prtcl_3d_mob.f      # Fortran Source Code (BEM Solver)
├── prtcl_3d_mob.dat    # Input Data File (Particle/Flow Parameters)
├── run_and_extract.py  # Script 1: Runs Fortran, extracts U/Omega, and writes CSV
├── plot_traction_vectors.py # Script 2: Plots 3D Traction Field and Velocities on PNG
├── run_wall_experiment.py # Script 3: MASTER SCRIPT - Automates the entire loop
├── plot_mobility_summary.py # Script 4: Generates final 2D summary plots
└── README.md
```

-----

## 2\. Setup and Compilation

### A. Compile the Fortran Solver

The BEM solver must be compiled first to create the executable:

```bash
# Compile the Fortran source code
gfortran prtcl_3d_mob.f -o prtcl_3d_mob
```

This command generates the executable file: `./prtcl_3d_mob`.

### B. Verify Input File

Ensure the primary input file (`prtcl_3d_mob.dat`) is configured for the **Wall Flow** case and includes the initial geometry you wish to test (e.g., Sphere or Ellipsoid).

-----

## 3\. Running the Experiment

The entire experiment is controlled by the master script, `run_wall_experiment.py`.

### A. Experiment Execution

Execute the script from your terminal:

```bash
python3 run_wall_experiment.py
```

### B. What the Script Does (The Workflow Loop)

The master script automates the following steps for every wall position defined in its `WALL_DISTANCES` array:

1.  **Update `prtcl_3d_mob.dat`:** Sets the current `wall` position and ensures the particle is driven by a consistent parallel force ($\mathbf{F}=(1.0, 0, 0)$) and zero torque ($\mathbf{T}=0$).
2.  **Run Simulation:** Executes `run_and_extract.py`.
      * `run_and_extract.py` calls `./prtcl_3d_mob`, feeding it the necessary interactive inputs for Wall Flow (Type 2).
      * The BEM solver runs and outputs the solution to `simulation_output.log` and the temporary traction data to `solution_vector.csv`.
3.  **Extract Mobility Data:** `run_and_extract.py` parses the log file, extracts the linear ($\mathbf{U}$) and angular ($\mathbf{\Omega}$) velocities, and **appends them to `mobility_summary.csv`**.
4.  **Visualize Traction:** Executes `plot_traction_vectors.py`, which reads `solution_vector.csv`, `mobility_summary.csv`, and generates the 3D traction plot (`traction_visualization.png`) annotated with the current velocity and gap distance.
5.  **Archive Results:** Moves the generated `.log` and `.png` files into the `Wall_Hindrance_Results/` directory, naming them with the current gap distance (e.g., `Gap_0.10_Wall_0.90_Traction.png`).

-----

## 4\. Analysis and Visualization

After the `run_wall_experiment.py` script finishes, you will have two primary outputs:

### A. Summary Data

The file `mobility_summary.csv` contains the raw results for every run.

### B. Generating Mobility Plots

The final step is to generate the 2D plots summarizing the hindrance effect.

```bash
python3 plot_mobility_summary.py
```

This script reads `mobility_summary.csv` and generates a multi-panel plot saved as **`Wall_Hindrance_Results/wall_hindrance_mobility_summary.png`**, illustrating $U_x$ vs. Surface Gap ($\delta$).

-----

## 5\. Experiment Customization (Ellipsoid Testing)

To test different scenarios (e.g., the ellipsoid orientations), modify the `run_wall_experiment.py` file before running the master script.

| Test Case | Aspect Ratios (`b/a c/a`) | Rotation Angles ($\phi_1, \phi_2, \phi_3$) | Effect |
| :--- | :--- | :--- | :--- |
| **Sphere (Baseline)** | `1.0 1.0` | `0.0 0.0 0.0` | Standard Hindrance |
| **Ellipsoid (Perpendicular)** | `5.0 1.0` | `0.0 0.0 0.0` | **Maximizes** Hindrance (Long axis perpendicular to wall) |
| **Ellipsoid (Parallel)** | `5.0 1.0` | `0.0 0.0 0.5` | **Minimizes** Hindrance (Long axis parallel to wall) |

The master script will automatically delete prior results, ensuring a clean run for the new configuration.
