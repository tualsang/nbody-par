#!/bin/bash
#SBATCH --partition=Centaurus
#SBATCH --time=01:00:00
#SBATCH --mem=64G

echo "Simulating solar system with dt = 200 and 5000000 steps..."
start_time=$(date +%s)
./nbody planet 200 5000000 10000 > solar.out
end_time=$(date +%s)
echo "Simulation finished in $(($end_time - $start_time)) seconds"

echo "Simulating 100 particles with dt = 1 and 10000 steps..."
start_time=$(date +%s)
./nbody 100 1 10000 100 > random_100.out
end_time=$(date +%s)
echo "Simulation finished in $(($end_time - $start_time)) seconds"

echo "Simulating 1000 particles with dt = 1 and 10000 steps..."
start_time=$(date +%s)
./nbody 1000 1 10000 100 > random_1000.out
end_time=$(date +%s)
echo "Simulation finished in $(($end_time - $start_time)) seconds"
