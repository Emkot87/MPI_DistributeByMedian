#!/bin/bash
{
gcc makeMnist.c -o makeMnist.out

gcc makeRand.c -o makeRand.out

mpicc ./MpiMnist.c -lm -o mpiMnist.out

mpicc ./Mpi.c -lm -o mpiRand.out

echo "--------------------------------------------------------32768 points with 386 diamensions of randomly generated values 0-100"


./makeRand.out 32768 386

mpirun -np 2 ./mpiRand.out 32768 386 

mpirun -np 4 ./mpiRand.out 32768 386 

mpirun -np 8 ./mpiRand.out 32768 386 

mpirun -np 16 ./mpiRand.out 32768 386 

mpirun -np 32 ./mpiRand.out 32768 386 


echo "--------------------------------------------------------262144 points with 256 diamensions of randomly generated values 0-100"


./makeRand.out 262144 256 

mpirun -np 2 ./mpiRand.out 262144 256 

mpirun -np 4 ./mpiRand.out 262144 256 

mpirun -np 8 ./mpiRand.out 262144 256 


echo "--------------------------------------------------------524288 points with 256 diamensions of randomly generated values 0-100"


./makeRand.out 524288 256 

mpirun -np 2 ./mpiRand.out 524288 256 

mpirun -np 4 ./mpiRand.out 524288 256 

mpirun -np 8 ./mpiRand.out 524288 256 

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
} 2>&1 | tee results.txt




