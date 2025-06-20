#!/bin/bash
cg
# Array of mechanism names
NAMES=("OConnaire" "burke" "gri30" "FFCM2_model" "sandiego")

# Loop over each name and run chemgen
for NAME in "${NAMES[@]}"
do
    python3 ~/chemgen/bin/chemgen.py "./$NAME.yaml" . \
        --compile \
        --custom-test custom_test.py \
        --n-points-test 10000
done
