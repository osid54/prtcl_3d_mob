import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

# --- Configuration ---
SUMMARY_FILE = "mobility_summary.csv"
PLOT_FILENAME = "wall_hindrance_mobility_summary.png"
OUTPUT_DIR = "Wall_Hindrance_Results" # Where the final plot will be saved

def plot_mobility_data():
    """Reads mobility data and generates 2D plots of velocity vs. surface gap."""
    
    # Check if the summary file exists
    if not os.path.exists(SUMMARY_FILE):
        print(f"ERROR: Mobility summary file not found: {SUMMARY_FILE}")
        print("Please ensure run_wall_experiment.py has run successfully.")
        return

    # 1. Load Data
    df = pd.read_csv(SUMMARY_FILE)
    
    # Sort data by Surface_Gap to ensure lines plot correctly
    df = df.sort_values(by='Surface_Gap', ascending=True)

    # 2. Create Figure and Subplots
    fig, axes = plt.subplots(3, 1, figsize=(8, 12), sharex=True)
    fig.suptitle('Particle Mobility vs. Surface Gap ($\delta$)', fontsize=14)
    
    # Define Gap variable for plotting
    gap = df['Surface_Gap']

    # --- Subplot 1: Primary Translation (Ux vs. Gap) ---
    # This plot shows the main hindrance effect.
    axes[0].plot(gap, df['Ux'], marker='o', linestyle='-', color='b')
    axes[0].set_ylabel('Translational Velocity $U_x$')
    axes[0].set_title('A) Parallel Mobility (Hindrance)')
    axes[0].grid(True, linestyle='--', alpha=0.6)
    
    # --- Subplot 2: Secondary Translation (Uy vs. Gap) ---
    # This plot should show values very close to zero if force is only in X.
    axes[1].plot(gap, df['Uy'], marker='s', linestyle='--', color='g')
    axes[1].set_ylabel('Translational Velocity $U_y$')
    axes[1].set_title('B) Perpendicular Mobility')
    axes[1].grid(True, linestyle='--', alpha=0.6)
    
    # --- Subplot 3: Angular Velocity (Oz vs. Gap) ---
    # This plot should show values very close to zero if torque/shear is zero.
    axes[2].plot(gap, df['Oz'], marker='^', linestyle=':', color='r')
    axes[2].set_ylabel('Angular Velocity $\Omega_z$')
    axes[2].set_title('C) Angular Mobility')
    axes[2].set_xlabel('Surface Gap $\delta = c_y - R - y_{wall}$')
    axes[2].grid(True, linestyle='--', alpha=0.6)
    
    # Finalize Layout
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    
    # 3. Save Plot
    output_path = os.path.join(OUTPUT_DIR, PLOT_FILENAME)
    plt.savefig(output_path)
    print(f"\nSuccessfully generated mobility plot: {output_path}")

# --- Execution ---
if __name__ == "__main__":
    plot_mobility_data()