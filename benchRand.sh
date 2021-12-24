#!/bin/bash
{

gcc makeRand.c -o makeRand.out

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


} 2>&1 | tee resultsRand.txt
