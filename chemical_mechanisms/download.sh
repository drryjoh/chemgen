#!/bin/bash

# URL of the file to be downloaded
URL="https://web.stanford.edu/group/haiwanglab/FFCM2/assets/data/optmodel/FFCM2_model.yaml"

# Output file name
OUTPUT="FFCM2_model.yaml"

# Download the file using curl
echo "Downloading $URL..."

curl -fLo "$OUTPUT" "$URL"

# Check if the download was successful
if [ $? -eq 0 ]; then
    echo "File downloaded successfully and saved as $OUTPUT"
else
    echo "Failed to download the file"
    exit 1
fi

