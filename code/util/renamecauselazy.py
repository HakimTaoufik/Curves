import yaml
import os

def get_file_paths(folder, all_files=None):
    if all_files is None:
        all_files = []
    
    # Iterate through all items in the folder
    for item in os.listdir(folder):
        full_path = os.path.join(folder, item)
        
        # If it's a file, add it to the list
        if os.path.isfile(full_path) and full_path.endswith('.yaml'):
            all_files.append(full_path)
        # If it's a directory, recursively search it
        elif os.path.isdir(full_path):
            get_file_paths(full_path, all_files)
            
    return all_files

# Example usage
# Specify your desired output directory
output_directory = "generated/text"  # Change this to your target directory
# Ensure the output directory exists
os.makedirs(output_directory, exist_ok=True)

# Source directory
source_folder = "generated/data"
all_file_paths = get_file_paths(source_folder)

for file_path in all_file_paths:
    # Load the YAML file
    with open(file_path, 'r') as f:
        data = yaml.safe_load(f)
    
    # Extract family name from the path - assuming it's the first subdirectory after source_folder
    relative_path = os.path.relpath(file_path, source_folder)
    path_parts = relative_path.split(os.sep)
    
    if len(path_parts) > 0:
        family_name = path_parts[0]  # Get the first directory as family name
    else:
        family_name = "unknown"  # Fallback if no subdirectory exists
    
    # Create a family-specific subfolder in the output directory
    family_output_dir = os.path.join(output_directory, family_name)
    os.makedirs(family_output_dir, exist_ok=True)
    
    # Set output filename to include family name
    file_name = os.path.basename(file_path).replace(".yaml", ".txt")
    output_file_path = os.path.join(family_output_dir, f"{family_name}_{file_name}")
    
    # Keep only unique (x,y) pairs
    unique_points = set()
    for point in data:
        # Convert to tuple to make it hashable for the set
        point_tuple = (point['x'], point['y'])
        unique_points.add(point_tuple)
    
    # Sort points by x-coordinate for consistency
    sorted_unique_points = sorted(unique_points)
    
    # Open the output text file
    with open(output_file_path, 'w') as f:
        # Write the header
        f.write("# CURVE 2D, normalized by CHORD\n")
        f.write("# COLUMN#1 = X-coordinate, COLUMN#2 = Y-coordinate\n")
        f.write("#-------------------------------------------------------------------------------\n")
        
        # Write each unique point in the format "x y"
        for x, y in sorted_unique_points:
            f.write(f"{x} {y}\n")