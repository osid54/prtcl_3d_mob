import os
import shutil
import re
import numpy as np

# --- Configuration ---
DATA_FILE = "prtcl_3d_mob.dat"
EXTRACT_SCRIPT = "run_and_extract.py"
PLOT_SCRIPT = "plot_traction_vectors.py"
OUTPUT_DIR = "Wall_Hindrance_Results"
SIMULATION_LOG = "simulation_output.log" # Used by run_and_extract.py

# Define the experiment range (Wall Positions)
# This controls the gap distance (Distance = Center_Y - Wall_Y)
# Assuming particle center (cy) is fixed at 1.0 (from your .dat file)
# Distance from center to surface is R=1.0 (assuming aspect ratios are 1.0 1.0)
WALL_DISTANCES = np.arange(1.0, -3.01, -.05)

# Define the fixed inputs needed for the Fortran executable via run_and_extract.py
# Iflow=2 (Wall Flow), Ienrd=1 (Read .dat file), Icompute=0 (Skip inverse)
INTERACTIVE_INPUT = '2\n1\n0\n'

SUMMARY_FILE = "mobility_summary.csv"

def prepare_environment():
    """Sets up the results directory and cleans up old data files."""
    if os.path.exists(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)
    print(f"Created output directory: {OUTPUT_DIR}")
    
    # --- NEW CLEANUP LOGIC ---
    # Delete the cumulative summary file to ensure only current experiment data is logged.
    if os.path.exists(SUMMARY_FILE):
        os.remove(SUMMARY_FILE)
    print(f"Cleaned up previous mobility summary file ({SUMMARY_FILE}).")

def update_dat_file(wall_position):
    """
    Edits prtcl_3d_mob.dat to set the new wall position and ensures 
    other parameters are clean for the wall hindrance test.
    """
    with open(DATA_FILE, 'r') as f:
        lines = f.readlines()

    updated_lines = []
    
    # Flags to check for the correct line positions
    is_wall_line = False
    
    # Pattern matching the data lines (check the space alignment carefully)
    # This is a fixed-width format, so careful matching is necessary
    center_pattern = r"^0\.00\s+[\d\.\s]+\s+cx, cy, cz"
    force_pattern = r"^0\.0\s+[\d\.\s]+\s+force"
    torque_pattern = r"^0\.0\s+[\d\.\s]+\s+torque"
    rotation_pattern = r"^0\.\s+0\.\s+0\.0\s+[\d\.\s]+"
    wall_pattern = r"^\s*[\-\.\d]+\s*wall:" 

    for line in lines:
        if "aspect ratios: b/a c/a" in line:
            # Enforce SPHERE (aspect ratios b/a=1.0, c/a=1.0)
            updated_lines.append(f"5.00 1.00        aspect ratios: b/a c/a\n") # <--- FIXED TO 1.0 1.0
        
        elif re.search(center_pattern, line):
            updated_lines.append("0.00 1.00 0.00   cx, cy, cz (coordinates of the particle center)\n")
        
        elif re.search(force_pattern, line):
            # Enforce Fx=1.0, Fy=Fz=0.0 (Force parallel to wall)
            updated_lines.append("1.0 0.0 0.0   force\n") # <--- FIXED TO PARALLEL FORCE
        
        elif re.search(torque_pattern, line):
            # Enforce zero torque
            updated_lines.append("0.0 0.0 0.0   torque\n")
        
        elif re.search(rotation_pattern, line):
            # Enforce zero rotation
            updated_lines.append("0. 0. 0.5 0. 0. 0.   rotation angles: phix, phiy,phiz\n")

        elif 'shear rate' in line:
            # Enforce zero shear rate
            updated_lines.append("0.0 0.0    shear rate\n")

        elif re.search(wall_pattern, line) or is_wall_line:
            # This is the line where the wall position is set
            updated_lines.append(f"{wall_position:.2f}    wall:  position of wall at y=wall (for Iflow=2 and 7)\n")
            is_wall_line = False # Reset flag

        else:
            updated_lines.append(line)

    with open(DATA_FILE, 'w') as f:
        f.writelines(updated_lines)
        
        # --- NEW CODE: ENFORCE FILE SYNC ---
        # This is the critical step to prevent the Fortran program from reading old, cached data.
        try:
            os.fsync(f.fileno())
            print(f"   -> Updated {DATA_FILE} with wall position {wall_position:.2f} and synchronized to disk.")
        except AttributeError:
            # fsync might not be available in some environments, but we try it.
            pass
        # --- END NEW CODE ---

def run_step(wall_position):
    """Executes the extraction and plotting scripts for one wall distance."""
    
    # 1. Update the .dat file
    update_dat_file(wall_position)
    
    # Calculate the gap for the filename
    gap_distance = 1.0 - wall_position # cy=1.0 is fixed
    file_prefix = f"Gap_{gap_distance:.2f}_Wall_{wall_position:.2f}"
    
    print(f"\n--- Running: {file_prefix} (Surface Gap: {gap_distance:.2f}) ---")
    
    # 2. Run the extraction script (Fortran simulation)
    # Execute the extraction script, passing the wall position as an argument
    try:
        # Note: We execute the main script (which calls run_simulation and then extract_data)
        # We pass the wall position as the second argument (sys.argv[1] in the target script)
        # This execution is now handled directly by the updated main in run_and_extract.py 
        # which you will make next.
        subprocess.run(["python3", EXTRACT_SCRIPT, str(wall_position)], check=True, capture_output=False) 
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Extraction failed for {file_prefix}. Check {SIMULATION_LOG}.")
        return
    """
    # 3. Rename the log file and plot output, and move them to the results folder
    
    # Plotting script output
    plot_filename = 'traction_visualization.png'
    new_plot_path = os.path.join(OUTPUT_DIR, f"{file_prefix}_Traction.png")
    
    # Execution of plotting script
    try:
        # The plotting script uses the generated CSV and OUT files
        subprocess.run(["python3", PLOT_SCRIPT, str(wall_position)], check=True, capture_output=True)
        shutil.move(plot_filename, new_plot_path)
        print(f"Image saved to {new_plot_path}")

        # Rename and move the log file for completeness
        new_log_path = os.path.join(OUTPUT_DIR, f"{file_prefix}_Log.log")
        shutil.move(SIMULATION_LOG, new_log_path)
        
    except FileNotFoundError:
        print(f"ERROR: Could not find visualization output file '{plot_filename}'. Plotting script failed.")
    except Exception as e:
        print(f"An unexpected error occurred during plotting or file moving: {e}")
    """


def main_experiment():
    """Main execution block for the experiment."""
    prepare_environment()
    
    # CRITICAL: Verify the INTERACTIVE_INPUT in run_and_extract.py is '2\n1\n0\n'
    print("WARNING: Ensure run_and_extract.py is hardcoded to use Flow Type '2' (Wall Flow).")
    
    for wall_pos in WALL_DISTANCES:
        run_step(wall_pos)

    print("\n=================================================")
    print("Experiment Complete. Results are in the 'Wall_Hindrance_Results' folder.")
    print("=================================================")

if __name__ == "__main__":
    import subprocess
    main_experiment()