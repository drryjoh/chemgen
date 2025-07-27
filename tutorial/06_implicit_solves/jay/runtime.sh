#!/bin/bash

total=0
count=0

for i in $(seq 1 50); do
    # Run the executable and capture its output
    output=$(./bin/bin/chemgen)

    # Extract the Backward Euler time line
    time_line=$(echo "$output" | grep "\[Backward Euler\] Time elapsed")

    # Extract the numeric time value using grep
    time_val=$(echo "$time_line" | grep -oE '[0-9]+\.[0-9]+')

    # If a time value was found, add to total
    if [[ -n "$time_val" ]]; then
        total=$(echo "$total + $time_val" | bc)
        count=$((count + 1))
    fi
done

# Print the average time
if [[ $count -gt 0 ]]; then
    average=$(echo "scale=6; $total / $count" | bc)
    echo "avg be time: $average seconds"
    echo "runs: $count"
else
    echo "No .exe"
fi
