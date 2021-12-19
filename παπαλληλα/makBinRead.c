#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


float **alloc_2d_init(int rows, int cols) {
    float *data = (float *)malloc(rows*cols*sizeof(float));
    float **array= (float **)malloc(rows*sizeof(float*));
    for (int i=0; i<rows; i++)
        array[i] = &(data[cols*i]);

    return array;
}

int main(int argc, char* argv[]){
    int N, d;
    N = (int) strtol(argv[1],NULL,10);
	d = (int) strtol(argv[2],NULL,10);

    float **pointsForFile;

    pointsForFile = alloc_2d_init(N,d);
    for(int i = 0; i < N ; i++){
        for(int j = 0; j < d ; j++){
            pointsForFile[i][j] = (float)rand()/((float)RAND_MAX/100);
        }
    }
    
    FILE *fptr;
    
    fptr = fopen("program2.bin","wb+");
    
    fwrite(pointsForFile, N*d,sizeof(float),fptr);
    fclose(fptr);

    float **vals;
    vals = alloc_2d_init(N,d);
    FILE *sad;

    sad = fopen("program2.bin","rb");

    fread(vals,sizeof(float),N*d,sad);
    fclose(sad);

    for(int i = 0 ; i<N ;i++){
        for(int j = 0 ; j < d ; j++){
            printf("%f\n",vals[i][j]);
        }
    }
	
}