#!/bin/bash
{
gcc makeMnist.c -o makeMnist.out

mpicc ./MpiMnist.c -lm -o mpiMnist.out

echo "--------------------------------------------------------16384 points with 784 diamensions of MNIST"

./makeMnist.out 16384

mpirun -np 2 ./mpiMnist.out 16384 784 

mpirun -np 4 ./mpiMnist.out 16384 784 

mpirun -np 8 ./mpiMnist.out 16384 784 

mpirun -np 16 ./mpiMnist.out 16384 784 

mpirun -np 32 ./mpiMnist.out 16384 784 

mpirun -np 64 ./mpiMnist.out 16384 784

echo "--------------------------------------------------------32768 points with 784 diamensions of MNIST"

./makeMnist.out 32768 

mpirun -np 2 ./mpiMnist.out 32768 784 

mpirun -np 4 ./mpiMnist.out 32768 784 

mpirun -np 8 ./mpiMnist.out 32768 784 

mpirun -np 16 ./mpiMnist.out 32768 784

echo "--------------------------------------------------------59200 points with 784 diamensions of MNIST"

./makeMnist.out 59200 

mpirun -np 2 ./mpiMnist.out 59200 784 

mpirun -np 4 ./mpiMnist.out 59200 784 

mpirun -np 8 ./mpiMnist.out 59200 784 

mpirun -np 16 ./mpiMnist.out 59200 784
} 2>&1 | tee resultsMnist.txt
