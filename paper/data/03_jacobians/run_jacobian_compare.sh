#!/bin/bash
set -e
source ~/python_environments/chemgen/bin/activate

# Array of mechanism names
NAMES=("ffcm2_h2" )
#"OConnaire" "burke" "gri30" "sandiego" "FFCM2_model")
OPTIONS=("--ignore-temp-dependence" "")
RUNS=("Ignoring Temperature Derivative" "Internal Energy")

# Loop over each mechanism
for NAME in "${NAMES[@]}"
do
	./create_random_data.py "$NAME.yaml" 1000
	for i in "${!OPTIONS[@]}"
	do
		OPTION="${OPTIONS[$i]}"
		DESC="${RUNS[$i]}"
		echo ">>> Running $NAME with setting: $DESC"
		python3 ~/chemgen/bin/chemgen.py "./$NAME.yaml" . \
			--pybind \
			$OPTION
		python3 ./src/setup_chemgen.py build_ext --inplace
		./run_random_data.py "$NAME.yaml" 1000
	done
done
