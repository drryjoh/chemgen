#!/bin/bash
set -e
source ~/python_environments/chemgen/bin/activate

# Array of mechanism names
NAMES=("ffcm2_h2" "OConnaire" "burke" "gri30" "sandiego" "FFCM2_model")
OPTIONS=("--ignore-temp-dependence" "--temperature-equation" "")
RUNS=("Ignoring Temperature Derivative" "Temperature Equation" "Internal Energy")

# Loop over each mechanism
for NAME in "${NAMES[@]}"
do
	for i in "${!OPTIONS[@]}"
	do
		OPTION="${OPTIONS[$i]}"
		DESC="${RUNS[$i]}"
		echo ">>> Running $NAME with setting: $DESC"
		python3 ~/chemgen/bin/chemgen.py "./$NAME.yaml" . \
			--compile \
			--custom-test custom_test.py \
			$OPTION
	done
done
