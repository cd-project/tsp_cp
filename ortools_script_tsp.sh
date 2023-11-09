#!/bin/bash

# Set the directory where your files are located
directory="/home/cuong/CLionProjects/TSP_CP/tsp_formatted_instance"

# Iterate over files in the directory
for file in "$directory"/*; do
    if [ -f "$file" ]; then
        # Check if the item is a regular file (not a directory)
        # Call your program with the file as an argument
        ./TSP_CP "$file" makespan
    fi
done