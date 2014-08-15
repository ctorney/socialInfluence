
/*test.cc*/
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <curand_kernel.h>
#include <iostream>
#include <algorithm>

#include <time.h>
#include "cuda_call.h"
#include "switcherKernel.h"
#include "params.h"
#include "cnpy.h"
#include "plotGrid.h"

int main() 
{

    // number of reps
    int numBlocks = 1;
    plotGrid* pg = new plotGrid;

    // length of grid
    int Nx = 64;
    int N = Nx * Nx;
    int N2 = 0.5 * N;
    int N4 = 0.5 * N2;
    int N_ALL = N * numBlocks;
    float* d_OU;
    cudaMalloc((void**)&d_OU, sizeof(float) *  (N_ALL) );
    CUDA_CALL(cudaMemset (d_OU, 0, sizeof(float) * (N_ALL)));

    float wg = 0.1f;


    dim3 threadGrid(Nx, Nx);
    curandState *devRands;
    CUDA_CALL(cudaMalloc((void **)&devRands, N_ALL * sizeof(curandState)));
    srand (time(NULL));
    initRands(threadGrid, numBlocks, devRands, rand());


    float* h_OU = new float [N_ALL];
    float* d_wg;
    CUDA_CALL(cudaMalloc((void**)&d_wg, sizeof(float) *  (N_ALL) ));
    float* h_wg = new float [N_ALL];
    int* d_states;
    CUDA_CALL(cudaMalloc((void**)&d_states, sizeof(int) * N_ALL));

    int* h_states = new int[N_ALL];

    int* d_blockTotals;
    CUDA_CALL(cudaMalloc((void**)&d_blockTotals, sizeof(int) * numBlocks));

    int* h_blockTotals = new int[numBlocks];
    int* h_blockTimes = new int[numBlocks];
    int infCount = N;
    //int infNums [8] = {1024,512,256,128,64,32,16,8};
    int infNums[8] = {256,128,64,32,16,8,4,2};
    const unsigned int shape[] = {1,N,infCount};

    float* results = new float[N*infCount];
    for (int i=0;i<N*infCount;i++)
       results[i]=0.0f;


    int net = 0;


    char fileName[18];
    sprintf(fileName, "cost2q%d.npy", net);


    CUDA_CALL(cudaMemset (d_states, 0, sizeof(int) * (N_ALL)));
    CUDA_CALL(cudaMemset (d_blockTotals, 0, sizeof(int) * (numBlocks)));


    for (int i=0;i<N_ALL;i++)
        h_wg[i]=1.0;
    CUDA_CALL(cudaMemcpy(d_wg, h_wg, (N_ALL) * sizeof(float), cudaMemcpyHostToDevice));


    for (int b=0;b<numBlocks;b++)
        h_blockTimes[b] = -1;
    int maxTime = 10000;
    int checkTime = 100000000;
    float sw = 1.0f - (float)net*0.2;;

    for (int t=0;t<maxTime;t++)
    {
        cout<<t<<endl;
        advanceTimestep(threadGrid, numBlocks, devRands, d_OU, d_wg, d_states, Nx, sw);
        CUDA_CALL(cudaMemcpy(h_states, d_states, (N_ALL) * sizeof(int), cudaMemcpyDeviceToHost));
        pg->draw(Nx, h_states);
        if (t%checkTime == -1 ) 
        {
            countStates(N, numBlocks, d_states, d_blockTotals, N_ALL);

            CUDA_CALL(cudaMemcpy(h_blockTotals, d_blockTotals, (numBlocks) * sizeof(int), cudaMemcpyDeviceToHost));
            bool allDone = true;
            for (int b=0;b<numBlocks;b++)
                if (h_blockTotals[b]>N2)
                {
                    if (h_blockTimes[b]<0)
                        h_blockTimes[b]=t;
                }
                else
                    allDone = false;
            if (allDone)
                break;
        }

    }

    /*
       float avTime = 0.0f;
       int count=0;
       for (int b=0;b<numBlocks;b++)
       if (h_blockTimes[b]>0)
       {
       avTime += (float)h_blockTimes[b];
       count++;
       }
       if (count>0)
       results[exnum * infCount + infnum] = avTime/(float)count;
       else
       results[exnum * infCount + infnum] = maxTime;
       if (avTime/(float)count > 100*checkTime)
       checkTime = checkTime * 10;
       if (checkTime > 10000)
       checkTime = 10000;

       cnpy::npy_save(fileName,results,shape,3,"w");
       }
     */


return 0;
}
