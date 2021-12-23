#include "mnist.h"

int main(int argc, char* argv[])
{
    // call to store mnist in array
    load_mnist();
    
    
    int N =(int) strtol(argv[1],NULL,10);
    int d = 784;

    float* valuesForBin = (float*)malloc(sizeof(float)*N*d);

    for(int j = 0 ; j < N ; j++){
        for(int k = 0 ; k < d ; k++){
            float temp = train_image[j][k];
            valuesForBin[j*d + k] = temp;
        }
    }

    FILE *fptr;
    
    fptr = fopen("programMnist.bin","wb+");
    
    fwrite(valuesForBin, N*d,sizeof(float),fptr);
    fclose(fptr);


    return 0;
}