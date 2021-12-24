# Parallel-2
First unzip the data folder  
Running benchScript.sh compiles all the files and runs benchmarks with both random values and mnist data set  

Running benchRand.sh compiles the rand variants and runs benchmarks with random values  

Running benchMnist.sh compiles the Mnist variants and runs benchmarks with the Mnist dataset  

if you want to run it yourself run benchScript to compile the files and run ./makeRand.out "number of points" "number of diamensions"
and then mpirun -np "number of processes" ./mpiRand.out "number of points" "number of diamensions"   


if you want to run Mnist yourself run benchScript to compile the files and run ./makeMnist.out "number of points" up to 60000
and then mpirun -np "number of processes" ./mpiRand.out "number of points" 784    
