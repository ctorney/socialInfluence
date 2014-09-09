
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

#define VISUAL 0
#if VISUAL
#include "plotGrid.h"
#endif

using namespace std;
int main() 
{

#if VISUAL
    plotGrid* pg = new plotGrid;
#endif

    // number of reps
    int numBlocks = 256;
    // length of grid
    int Nx = 8;
    int N = Nx * Nx;
    int N2 = 0.5 * N;
    int N4 = 0.5 * N2;
    int N_ALL = N * numBlocks;



    dim3 threadGrid(Nx, Nx);
    curandState *devRands;
    CUDA_CALL(cudaMalloc((void **)&devRands, N_ALL * sizeof(curandState)));

    srand (time(NULL));
    initRands(threadGrid, numBlocks, devRands, rand());


    float* d_wg;
    CUDA_CALL(cudaMalloc((void**)&d_wg, sizeof(float) *  (N_ALL) ));
    int* d_states;
    CUDA_CALL(cudaMalloc((void**)&d_states, sizeof(int) * N_ALL));
    int* d_states2;
    CUDA_CALL(cudaMalloc((void**)&d_states2, sizeof(int) * N_ALL));

    float* d_up;
    CUDA_CALL(cudaMalloc((void**)&d_up, sizeof(float) *  (N + 1) ));

    float* h_up = new float [N+1];

    float* d_down;
    CUDA_CALL(cudaMalloc((void**)&d_down, sizeof(float) *  (N + 1) ));

    float* h_down = new float [N+1];

    int* d_upcount;
    CUDA_CALL(cudaMalloc((void**)&d_upcount, sizeof(int) *  (N + 1) ));

    int* h_upcount = new int [N+1];

    int* d_downcount;
    CUDA_CALL(cudaMalloc((void**)&d_downcount, sizeof(int) *  (N + 1) ));

    int* h_downcount = new int [N+1];



    int* d_blockTotals;
    CUDA_CALL(cudaMalloc((void**)&d_blockTotals, sizeof(int) * numBlocks));

    float* h_wg = new float [N_ALL];
    int* h_states = new int[N_ALL];
    int* h_blockTotals = new int[numBlocks];
    int* h_blockTimes = new int[numBlocks];
    int wgCount = 1;

    const unsigned int shape[] = {N+1,2};

    float* results = new float[(N+1)*2];
    for (int i=0;i<(N+1)*2;i++)
       results[i]=0.0f;



    for (int G=0;G<wgCount;G++)
    {
        float wg = 0.3;//5 + 0.2 * float(G);
        for (int i=0;i<N_ALL;i++)
        {
            h_wg[i]=wg;
            float unum =rand()/double(RAND_MAX);
   //         cout<<unum<<endl;
 //           if (unum<0.2345)
 //               h_wg[i]+=0.833;
        }
        CUDA_CALL(cudaMemcpy(d_wg, h_wg, (N_ALL) * sizeof(float), cudaMemcpyHostToDevice));


        for (int b=0;b<numBlocks;b++)
            h_blockTimes[b] = -1;
        int maxTime = 1000;
        int checkTime = 1000;
        float sw = 20.0f;

        char fileName[30];
        sprintf(fileName, "potential%d-%d.npy", int(10*wg),int(100.0*sw));
//        cout<<fileName<<endl;

        CUDA_CALL(cudaMemset (d_states, 0, sizeof(int) * (N_ALL)));
        CUDA_CALL(cudaMemset (d_blockTotals, 0, sizeof(int) * (numBlocks)));
        CUDA_CALL(cudaMemset (d_up, 0, sizeof(float) * (N + 1)));
        CUDA_CALL(cudaMemset (d_down, 0, sizeof(float) * (N + 1)));
        CUDA_CALL(cudaMemset (d_upcount, 0, sizeof(int) * (N + 1)));
        CUDA_CALL(cudaMemset (d_downcount, 0, sizeof(int) * (N + 1)));





        for (int t=0;t<maxTime;t++)
        {
            
            advanceTimestep(threadGrid, numBlocks, devRands, d_wg, d_states, Nx, sw, t);
            recordData(threadGrid, numBlocks, d_states, d_states2, Nx, d_up, d_down, d_upcount, d_downcount, t);
            /*
            CUDA_CALL(cudaMemcpy(h_states, d_states, (N_ALL) * sizeof(int), cudaMemcpyDeviceToHost));
            int countUp = 0;
            for (int i=0;i<N_ALL;i++)
                if (h_states[i]>0)
                    countUp++;
            cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl<<countUp<<endl;
//            */
#if VISUAL
            CUDA_CALL(cudaMemcpy(h_states, d_states, (N_ALL) * sizeof(int), cudaMemcpyDeviceToHost));
            pg->draw(Nx, h_states);
#endif
            if (t%checkTime == 0 ) 
            {
               // cout<<t<<endl;
                countStates(N, numBlocks, d_states, d_blockTotals, N_ALL);

                CUDA_CALL(cudaMemcpy(h_blockTotals, d_blockTotals, (numBlocks) * sizeof(int), cudaMemcpyDeviceToHost));
                bool allDone = true;
                for (int b=0;b<numBlocks;b++)
                {
                    if (h_blockTotals[b]>0.75*N)
                    {
    //                    cout<<"block total : "<<h_blockTotals[b]<<endl;
                        if (h_blockTimes[b]<0)
                            h_blockTimes[b]=t;
                    }
                    else
                    {
             //           cout<<b<<" block done"<<endl;
                        allDone = false;
                    }
                }
                if (allDone)
                {
                    break;
                    /*
                    for (int b=0;b<numBlocks;b++)
                        h_blockTimes[b] = -1;
                    CUDA_CALL(cudaMemset (d_states, 0, sizeof(int) * (N_ALL)));
                    */
                } 
            }

        }
        CUDA_CALL(cudaMemcpy(h_up, d_up, (N + 1) * sizeof(float), cudaMemcpyDeviceToHost));
        CUDA_CALL(cudaMemcpy(h_down, d_down, (N + 1) * sizeof(float), cudaMemcpyDeviceToHost));
        CUDA_CALL(cudaMemcpy(h_upcount, d_upcount, (N + 1) * sizeof(int), cudaMemcpyDeviceToHost));

        for (int i=0;i<N+1;i++)
        {
            results[2*i]=h_up[i];
            results[2*i+1]=h_down[i];
            cout<<i/float(N)<<" : "<<h_up[i]<<" : "<<h_down[i]<<" : "<<h_upcount[i]<<endl;
        }

        cnpy::npy_save(fileName,results,shape,2,"w");
    }
    return 0;
}
