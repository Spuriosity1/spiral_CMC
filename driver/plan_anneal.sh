#!/bin/bash

# Check if input file is provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <input_file> <L> [program_args...]"
    echo "Example: $0 numbers.txt 10"
    exit 1
fi

INPUT_FILE="$1"
L="$2"
shift 2  # Remove first two arguments, leaving any additional program args

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: File '$INPUT_FILE' not found"
    exit 1
fi

PROGRAM=build/anneal

# Check if program exists and is executable
if [ ! -x "$PROGRAM" ]; then
    echo "Error: Program '$PROGRAM' not found or not executable"
    exit 1
fi

Tcold=0.0001
SEED=0xa0ff19102

# Read file line by line
while IFS= read -r line || [ -n "$line" ]; do
    # Skip empty lines and comments
    [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
    
    # Read two numbers from the line
    read -r j2 j3 <<< "$line"
    
    # Check if both numbers were found
    if [ -z "$j2" ] || [ -z "$j3" ]; then
        echo "Warning: Skipping invalid line: $line"
        continue
    fi
    
    # Pass the numbers to the program
    echo "$PROGRAM" "$L" --J1 1 --J2 $j2 --J3 $j3 "$@" --T_cold $Tcold -o ../tmp -s $SEED
    echo "$PROGRAM" "$L" --J1 -1 --J2 $j2 --J3 $j3 "$@" --T_cold $Tcold -o ../tmp -s $SEED
    
done < "$INPUT_FILE"

#echo "Processing complete"
