

/*
 * Convention model on a dynamic small world network 
*/

#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include "main.h"
#include "cuda_call.h"
#include <curand_kernel.h>
#include "switcherKernel.h"
#include "params.h"
#include <algorithm>
#include "mkl.h"
#include "cnpy.h"
#include <time.h>

using namespace std;

int main(int argc, char** argv) 
{
    int inNetwork=atoi(argv[1]);

    // number of reps
    int numBlocks = 256;

    // length of grid
    int N_x = 8;
    int N = N_x * N_x;
    int N2 = 0.5 * N;
    int N4 = 0.5 * N2;
    int N_ALL = N * numBlocks;
    float* d_OU;
    cudaMalloc((void**)&d_OU, sizeof(float) *  (N_ALL) );
    CUDA_CALL(cudaMemset (d_OU, 0, sizeof(float) * (N_ALL)));

    float wg = 0.1f;


    dim3 threadGrid(N_x, N_x);
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
    int* d_states2;
    CUDA_CALL(cudaMalloc((void**)&d_states2, sizeof(int) * N_ALL));

    int* h_states = new int[N_ALL];
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

    CUDA_CALL(cudaMemset (d_up, 0, sizeof(float) * (N + 1)));
    CUDA_CALL(cudaMemset (d_down, 0, sizeof(float) * (N + 1)));
    CUDA_CALL(cudaMemset (d_upcount, 0, sizeof(int) * (N + 1)));
    CUDA_CALL(cudaMemset (d_downcount, 0, sizeof(int) * (N + 1)));




    int* d_blockTotals;
    CUDA_CALL(cudaMalloc((void**)&d_blockTotals, sizeof(int) * numBlocks));

    int* h_blockTotals = new int[numBlocks];
    int* h_blockTimes = new int[numBlocks];
    int infCount = N;
    //int infNums [8] = {1024,512,256,128,64,32,16,8};
    int infNums[8] = {256,128,64,32,16,8,4,2};
    const unsigned int shape[] = {N+1,2};

    float* results = new float[N*infCount];
    for (int i=0;i<N*2+2;i++)
       results[i]=0.0f;


 //   for (int net=0;net<10;net++)
        int net = inNetwork;

    {
        char fileName[20];
//        for (int exnum=N-1;exnum>=0;exnum--)
        int exnum =50;
        {
            sprintf(fileName, "potential%d-%d.npy", net,exnum);
 //           for (int infnum=0;infnum<infCount;infnum++)
            {
                int infnum = 14;
                cout<<exnum<<endl;
                float sw = 1.0f - (float)net*0.2;;
                CUDA_CALL(cudaMemset (d_states, 0, sizeof(int) * (N_ALL)));
                CUDA_CALL(cudaMemset (d_blockTotals, 0, sizeof(int) * (numBlocks)));

                //setInformation(h_wg, wg, exnum, infnum + 1, N, numBlocks);
                setInformation(h_wg, wg, exnum, N, N, numBlocks);

                CUDA_CALL(cudaMemcpy(d_wg, h_wg, (N_ALL) * sizeof(float), cudaMemcpyHostToDevice));


                for (int b=0;b<numBlocks;b++)
                    h_blockTimes[b] = -1;
                int maxTime = 100000;
                int checkTime = 100;

                for (int t=0;t<maxTime;t++)
                {
                    advanceTimestep(threadGrid, numBlocks, devRands, d_OU, d_wg, d_states, N_x, sw);
                    recordData(threadGrid, numBlocks, d_states, d_states2, N_x, d_up, d_down, d_upcount, d_downcount, t);
                    if (t%checkTime == 0 ) 
                    {
                        countStates(N, numBlocks, d_states, d_blockTotals, N_ALL);

                        CUDA_CALL(cudaMemcpy(h_blockTotals, d_blockTotals, (numBlocks) * sizeof(int), cudaMemcpyDeviceToHost));
                        bool allDone = true;
                        for (int b=0;b<numBlocks;b++)
                        {
                            if (h_blockTotals[b]>N2)
                            {
                                if (h_blockTimes[b]<0)
                                    h_blockTimes[b]=t;
                            }
                            else
                                allDone = false;
                        }
                        if (allDone)
                            break;
                    }

                }
                CUDA_CALL(cudaMemcpy(h_up, d_up, (N + 1) * sizeof(float), cudaMemcpyDeviceToHost));
                CUDA_CALL(cudaMemcpy(h_down, d_down, (N + 1) * sizeof(float), cudaMemcpyDeviceToHost));
                CUDA_CALL(cudaMemcpy(h_upcount, d_upcount, (N + 1) * sizeof(int), cudaMemcpyDeviceToHost));

                for (int i=0;i<N+1;i++)
                {
                    results[i*2]=h_up[i];
                    results[i*2+1]=h_down[i];
                    cout<<h_up[i]<<":"<<h_down[i]<<":"<<h_upcount[i]<<endl;
                }

            }
        }
        cnpy::npy_save(fileName,results,shape,2,"w");
    }
    /*
       CUDA_CALL(cudaMemcpy(h_states, d_states, (N_ALL) * sizeof(float), cudaMemcpyDeviceToHost));
       for (int b=0;b<numBlocks;b++)
       for (int i=0;i<N_x;i++)
       {
       for (int j=0;j<N_x;j++)
       cout<<h_states[(b*N) + (j*N_x) + i]<<" ";
       cout<<endl;
       }
       */
}


void setInformation(float* wg, float base, int ex_num, int split_num, int N, int numBlocks)
{
    // ex_num is the expected number of individuals that will be made to switch
    // split_num is the number of individuals that will be influenced, e.g. if
    // ex_num is 1 and split_num is N, then each individual will be given an extra 
    // 1/N probability of being correct.



    //calculate the baseline probability to switch for the uninfluenced population
    float base_prob = 0.5 - 0.5*erf(sqrt(base)*( (0.25*8.0*log(ws/(1.0-ws))/base)  - 1.0f));
    // calculate the boost given to each of the split_num individuals
    float p = (float)ex_num/(float)split_num;
    // now we need to know the value of wg that gives the probability boost
    float* in = new float[1];
    float* out = new float[1];
    float wginc=0.0;
    if ((base_prob+p)<1.0f)
    {
        in[0] =1.0-2.0*(base_prob+p);
        vsErfInv(1, in ,out);
        float A = 1;
        float B = out[0];
        float C = -2.0*log(ws/(1.0-ws));
        // this is the increment that needs to be added
        wginc = powf((-B + sqrt(B*B - 4.0*A*C))/(2.0*A),2) - base;
    }
    else
    {
        wginc = 99.0f;
    }

    float totalP = 16.0f * ex_num / (float)N;
    wginc = totalP/(float)split_num;

    // first apply the baseline wg
    for (int b=0;b<numBlocks;b++)
        for (int i=0;i<N;i++)
            wg[b * N + i] = base;

    // next make a random selection and boost those individuals
    int* listN = new int[N];
    for (int i=0;i<N;i++)
        listN[i]=i;

    int* randN = new int[N];
    srand(0);

    for (int b=0;b<numBlocks;b++)
    {
        memcpy(randN,listN,N*sizeof(int));

        random_shuffle(randN, randN+N); 
        for (int i=0;i<split_num;i++)
            wg[b * N + randN[i]] += wginc;



    }
 //   cout<<ex_num<<":"<<split_num<<":"<<wginc<<endl;
    delete [] randN;
    delete [] listN;
    delete [] in;
    delete [] out;
}
