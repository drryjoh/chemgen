#!/usr/bin/bash
#SBATCH --job-name=test_job
#SBATCH --output=chemgen_h2.out
#SBATCH --error=chemgen.err
#SBATCH --time=50:00
#SBATCH -p haiwang 
#SBATCH --cpus-per-task=10
#SBATCH --mem=8GB

. ~/.ryan_modules
echo "Nodes allocated for this job: $SLURM_NODELIST"
export TBB_NUM_THREADS=10
NUMBER=100

# Loop to double the number 10 times
for ((i=1; i<=10; i++))
do
    echo "Number of points $NUMBER"
    python3 ../chemgen/bin/chemgen.py ffcm2_h2 . --custom-source ../chemgen/test/test_thread_points_source/write_source_threaded.py --custom-test ../chemgen/test/test_thread_points_source/test_thread_points_source.py --cmake --n-points-test $NUMBER
    cd build
    make
    ./bin/chemgen
    cd ../
    # Double the value of NUMBER
    NUMBER=$((NUMBER * 2))
done
echo "running!"
