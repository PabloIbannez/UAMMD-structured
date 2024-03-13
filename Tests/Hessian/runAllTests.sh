#!/bin/bash

# Program to check
program_name="UAMMDlauncherDouble"

# Directory to skip
skip_dir="scripts"

# Check if UAMMDlauncherDouble is in PATH
if ! command -v "$program_name" &> /dev/null
then
    echo "$program_name could not be found in your PATH. Please ensure it is installed and accessible."
    exit 1
fi

# Loop through each item in the current directory
for dir in */ ; do
    # Remove the trailing '/'
    dir=${dir%*/}

    # Skip the directory if it matches the skip_dir
    if [ "$dir" = "$skip_dir" ]; then
        continue
    fi

    # Check if it's a directory and not a file
    if [ -d "$dir" ]; then
        (cd "$dir" && ./run.sh)
    fi
done
