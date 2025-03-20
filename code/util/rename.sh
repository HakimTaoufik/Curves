#!/bin/bash

# Set the directory to search
TEXT_DIR="generated/text"

# Find all files with the pattern X_points_Y in the text directory and its subdirectories
find "$TEXT_DIR" -type f -name "*_points_*" | while read file; do
    # Extract the directory path and filename
    dir=$(dirname "$file")
    filename=$(basename "$file")
    
    # Check if the filename matches the pattern X_points_Y
    if [[ "$filename" =~ ^(.+)_points_(.+)$ ]]; then
        # Extract the Y part
        y_part="${BASH_REMATCH[2]}"
        
        # Create the new filename: points_Y
        new_filename="points_$y_part"
        
        # Create the full path for the new file
        new_file="$dir/$new_filename"
        
        # Rename the file
        mv "$file" "$new_file"
        echo "Renamed: $filename -> $new_filename"
    fi
done

echo "File renaming complete."