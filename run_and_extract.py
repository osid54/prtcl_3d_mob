import subprocess
import re
import os

# --- Configuration ---
EXECUTABLE = "./prtcl_3d_mob"
LOG_FILE = "simulation_output.log"
CSV_FILE = "solution_vector.csv"

# --- Main function to run simulation and extract data ---

def run_simulation():
    """Runs the Fortran executable and saves all output to a log file."""
    print(f"1. Running simulation: {EXECUTABLE}")
    try:
        # Run the executable and capture output
        # We need to simulate the interactive input for the matrix inverse prompt
        # Send '0' to the program to skip computing the matrix inverse
        process = subprocess.run(
            [EXECUTABLE],
            capture_output=True,
            text=True,
            input='2\n1\n0\n',
            check=True  # Raise error if non-zero exit code
        )
        
        # Save the full output to the log file
        with open(LOG_FILE, 'w') as f:
            f.write(process.stdout)
            f.write(process.stderr)

        print(f"   -> Simulation complete. Output saved to {LOG_FILE}")
        return process.stdout

    except FileNotFoundError:
        print(f"\nERROR: Executable '{EXECUTABLE}' not found. Did you run 'make' successfully?")
        return None
    except subprocess.CalledProcessError as e:
        print(f"\nERROR: Simulation failed to run. Check the compiler setup.")
        print("--- Fortran Error Output ---")
        print(e.stdout)
        print(e.stderr)
        return None
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")
        return None

def extract_data(log_content, wall_position):
    """Parses the log content to extract key results and the solution vector."""
    print("2. Extracting key data and solution vector...")
    
    # 1. Extract System Size (Mdim)
    size_match = re.search(r"system size:\s*(\d+)", log_content, re.IGNORECASE)
    if not size_match:
        print("   -> ERROR: Could not find system size. Execution terminated prematurely.")
        return
        
    Mdim = int(size_match.group(1))
    M = int((Mdim - 6) / 3)
    
    print(f"   -> Detected System Size (Mdim): {Mdim}")
    print(f"   -> Deduced Number of Elements (M): {M}")

    # 2. Extract U and Omega (Linear and Angular Velocity)
    velocity_section = re.search(r"Linear and Angular Velocity\s*([\s\S]+)", log_content, re.IGNORECASE)
    
    if velocity_section:
        all_numbers = re.findall(r"[-]?\d+\.\d{1,}", log_content)

        # We expect Mdim numbers for the full solution vector (3M tractions + 6 velocities)
        num_found = len(all_numbers)
        if num_found >= Mdim:

            # Take the last Mdim numbers as the solution vector (sln)
            solution_vector = all_numbers[-Mdim:] 
            
            # Velocities are the last 6 entries
            Ux, Uy, Uz, Ox, Oy, Oz = [float(x) for x in solution_vector[-6:]]

            # -------------------------------------------------------------
            # --- NEW CODE: Append to MOBILITY SUMMARY CSV ---
            # -------------------------------------------------------------
            
            SUMMARY_FILE = "mobility_summary.csv"
            
            # Determine if we need to write the header (check if file is new/empty)
            write_header = not os.path.exists(SUMMARY_FILE) or os.path.getsize(SUMMARY_FILE) == 0
            
            with open(SUMMARY_FILE, 'a') as f:
                if write_header:
                    f.write("Wall_Position,Surface_Gap,Ux,Uy,Uz,Ox,Oy,Oz\n")
                
                # Calculate surface gap (assuming cy=1.0)
                gap = 1.0 - wall_position 
                
                f.write(f"{wall_position:.2f},{gap:.2f},{Ux:.10f},{Uy:.10f},{Uz:.10f},{Ox:.10f},{Oy:.10f},{Oz:.10f}\n")
            
            print("\n--- Summary of Results ---")
            print(f"Wall Position: {wall_position:.2f}, Gap: {gap:.2f}")
            print(f"Linear Velocity (U): {Ux:.10f}, {Uy:.10f}, {Uz:.10f}")
            print(f"Angular Velocity (Omega): {Ox:.10f}, {Oy:.10f}, {Oz:.10f}")
            print(f"Mobility data appended to {SUMMARY_FILE}")
            
            # -------------------------------------------------------------
            # --- TRACTION VECTOR CSV (Rest of your existing code) ---
            # -------------------------------------------------------------
            
            # Tractions are the first 3*M entries
            traction_vector = solution_vector[:3*M]
            
            with open(CSV_FILE, 'w') as f:
                f.write("Element_Index,Traction_X,Traction_Y,Traction_Z\n")
                
                # Write M rows of organized data
                # Tractions are stored sequentially by component (Fx1, Fx2..., Fy1, Fy2..., Fz1, Fz2...)
                for i in range(M):
                    fx = traction_vector[i]          # Fx component for element i
                    fy = traction_vector[i + M]      # Fy component for element i
                    fz = traction_vector[i + 2 * M]  # Fz component for element i
                    f.write(f"{i+1},{fx},{fy},{fz}\n")

        else:
            print(f"   -> ERROR: Found {num_found} numbers, expected {Mdim}. Data incomplete.")

    else:
        print("   -> ERROR: Could not find 'Linear and Angular Velocity' section.")
        
# NOTE: The main block of run_and_extract.py also needs modification
# to pass 'wall_position' to 'extract_data' (see run_wall_experiment.py step 2).

# --- Execution (at the end of run_and_extract.py) ---
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("ERROR: Wall position argument missing. Please run via run_wall_experiment.py.")
        exit(1)
        
    try:
        current_wall_pos = float(sys.argv[1])
    except ValueError:
        print("ERROR: Wall position argument is not a valid number.")
        exit(1)

    if not os.path.exists(EXECUTABLE):
        print(f"ERROR: Executable '{EXECUTABLE}' not found. Please ensure you ran 'make'.")
    else:
        log_data = run_simulation()
        if log_data:
            # Pass the wall position to the extraction function
            extract_data(log_data, current_wall_pos) 
        print("\n*** Workflow complete. ***")