#!/usr/bin/bash
#SBATCH --job-name=test_job
#SBATCH --output=test_job.%j.out
#SBATCH --error=test_job.%j.err
#SBATCH --time=10:00
#SBATCH -p haiwang 
#SBATCH -c 10
#SBATCH --mem=8GB

. ~/.ryan_modules

NUMBER=500

# Loop to double the number 10 times
for ((i=1; i<=10; i++))
do
    echo "Number of points $NUMBER"
    python3 ../chemgen/bin/chemgen.py FFCM2_model . --custom-source ../chemgen/test/test_thread_points_source/write_source_threaded.py --custom-test ../chemgen/test/test_thread_points_source/test_thread_points_source.py --cmake --n-points-test $NUMBER
    cd build
    make
    ./bin/chemgen
    cd ../
    # Double the value of NUMBER
    NUMBER=$((NUMBER * 2))
done
echo "running!"
