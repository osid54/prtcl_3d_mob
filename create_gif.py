import os
import glob
from PIL import Image

# --- Configuration ---
OUTPUT_DIR = "Wall_Hindrance_Results"
GIF_FILENAME = "wall_hindrance_animation.gif"

def create_animation(directory, output_filename):
    """
    Creates an animated GIF from all PNG files found in the specified directory.
    The frames are ordered alphabetically by filename.
    """
    print(f"1. Looking for PNG files in: {directory}")
    
    # 1. Get the list of all PNG files
    # The glob pattern ensures they are loaded in alphabetical/numeric order
    # (e.g., Gap_0.05...png comes before Gap_0.50...png)
    png_files = sorted(glob.glob(os.path.join(directory, "*.png")))
    
    if not png_files:
        print("ERROR: No PNG files found in the directory. Did the simulations run?")
        return

    # 2. Open all images using Pillow
    images = []
    
    # Open the first image separately to get the base for saving
    first_image = Image.open(png_files[0])
    
    # Load the rest of the images
    for filename in png_files[1:]:
        images.append(Image.open(filename))

    print(f"2. Found {len(png_files)} frames. Generating GIF...")
    
    # 3. Save the animated GIF
    # duration controls the speed in milliseconds (e.g., 500ms = 0.5s per frame)
    first_image.save(
        output_filename,
        save_all=True,
        append_images=images,
        duration=1000,  # Milliseconds per frame (Adjust speed here!)
        loop=0          # 0 means loop forever
    )
    
    print(f"3. Successfully created GIF: {output_filename}")
    print("   -> Check your project directory in Windows File Explorer to view.")


# --- Execution ---
if __name__ == "__main__":
    output_path = os.path.join(OUTPUT_DIR, GIF_FILENAME)
    
    # Run the function, but look inside the output directory for the frames
    create_animation(OUTPUT_DIR, output_path)