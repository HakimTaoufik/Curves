#!/bin/bash

# Create base destination directory if it doesn't exist
mkdir -p ~/PSC/tst

# Navigate to PSC directory
# Note: Replace /path/to/PSC with the actual path if needed
cd ~/BDD-PSC || { echo "Cannot find PSC directory"; exit 1; }

# Find all json files in the PSC/A*/B/C/ structure
find . -path "./Forme_*/Formes_samplees/json/*.json" -type f | while read -r file; do
    # Extract just the filename from the path
    filename=$(basename "$file")
    
    # Extract family name (everything before the second underscore)
    # This handles Z_1_23.json where Z_1 is the family name
    if [[ "$filename" =~ ^([^_]+_[^_]+)_.*\.json$ ]]; then
        family_name="${BASH_REMATCH[1]}"
    else
        # Fallback for filenames with only one underscore
        family_name=${filename%%_*}
    fi
    
    # Create family directory if it doesn't exist
    mkdir -p ~/json/"$family_name"
    
    # Copy the file to the appropriate family directory
    cp "$file" ~/json/"$family_name"/
done

echo "All JSON files have been organized by family name in ~/json/"