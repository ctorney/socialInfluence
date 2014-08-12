
#include <curand_kernel.h>
#include <stdio.h>
#include "params.h"

__device__ int getIndex(int t_x, int t_y)
{
    // calculate full index from a grid position 
    int indx = __mul24(t_y,blockDim.x) + t_x;
    return __mul24(blockDim.y, __mul24(blockIdx.x, blockDim.x)) + indx;

}
        

__global__ void d_initRands(curandState *state, int seed)
{
    int id = getIndex(threadIdx.x, threadIdx.y);

    /* Each thread gets same seed, a different sequence 
     *        number, no offset */
    curand_init(seed, id, 0, &state[id]);
}

__global__ void d_generateOU(curandState *state, float* ou_process, float* wg)
{

    int id = getIndex(threadIdx.x, threadIdx.y);

    float last = ou_process[id];
    float nm = curand_normal(&state[id]);
    
    float dt = 1.0f;
    float wgi = wg[id];
    float lambda = 1.0f;
    float sigma = 1.0f;

    ou_process[id] = (last*exp(-wgi * dt)) + (lambda * (1 - exp(-wgi*dt))) + ( nm * sigma * sqrtf((1.0f-exp(-2.0f*wgi*dt))/(2.0f*wgi)) );

}
__global__ void d_updateStates(int* states, float* ou_process, float* wg, int N_x, curandState* d_rands, float sw)
{
    int id = getIndex(threadIdx.x, threadIdx.y);
    int edges=8;
    int neigh[8][2] = { { 1, 1 }, { 1, 0 }, { 1, -1 } , { 0, 1 }, { 0, -1 }, { -1, -1 } , { -1, 0 }, { -1, 1 } };
    int deltan = 0;
   
    for (int e=0;e<edges;e++)
    {
        if (curand_uniform(&d_rands[id])<sw)
        {
            int n2 = curand_uniform(&d_rands[id])*8;
            if (n2==8) n2 = 0;

            int x_n = (((threadIdx.x + neigh[n2][0]) % N_x) + N_x) % N_x;
            int y_n = (((threadIdx.y + neigh[n2][1]) % N_x) + N_x) % N_x;

            int n2_id = getIndex(x_n, y_n);
            if (states[n2_id]>0.5)
                deltan++;

        }
        else
        {
            int x_n = curand_uniform(&d_rands[id])*N_x;
            int y_n = curand_uniform(&d_rands[id])*N_x;

            if (x_n==N_x) x_n = 0;
            if (y_n==N_x) y_n = 0;

            int n2_id = getIndex(x_n, y_n);
            int n2 = curand_uniform(&d_rands[id])*N_x*N_x;
 //           n2 += (N_x*N_x) * blockIdx.x;
            if (states[n2_id]>0.5)
                deltan++;

        }


    }

    // deltan is N+ right now but we want (N+ - N-)
    deltan*=2;
    deltan-=edges;


    float pup = exp(-4.0f*wg[id]*ou_process[id]);
    float pall = pup*powf((1.0f - ws)/ws,deltan);
    int newState;
    if (pall<1.0f)
        newState = 1;
    else
        newState = 0;

    __syncthreads();
 //   printf("%d %d %d %d %d %d %d %d %d\n", threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, blockDim.x, blockDim.y, gridDim.x, gridDim.y, id);

    states[id] = newState;
 //   if (curand_uniform(&d_rands[id])<0.01)
 //       states[id] = 1;
}
__global__ void d_recordData(int* states, int* states2, int N_x, float* d_up, float* d_down, int* d_upcount, int* d_downcount, int t)
{

    int group_id = threadIdx.y * N_x + threadIdx.x;

    int N = N_x*N_x;

    if ((group_id==0)&&(blockIdx.x==0))
        for (int b=0;b<gridDim.x;b++)
        {
            if (t==0)
                for (int i=0;i<N;i++)
                    states2[b * N + i] = states[b * N + i];
            else
            {
                int totalUp = 0;
                for (int i=0;i<N;i++)
                    if (states2[b * N + i] > 0.5)
                        totalUp++;


                int nowDown = 0;
                for (int i=0;i<N;i++)
                    if ((states2[b * N + i] > 0.5)&&(states[b * N + i] < 0.5))
                        nowDown++;

                int nowUp = 0;
                for (int i=0;i<N;i++)
                    if ((states2[b * N + i] < 0.5)&&(states[b * N + i] > 0.5))
                        nowUp++;


                d_upcount[totalUp]+=1;
                int c = d_upcount[totalUp];
                //           printf("%d %d %d %d\n",t, totalUp,nowDown, nowUp);
                d_down[totalUp] = (nowDown/(float)N)/(float)c + (c-1)*d_down[totalUp]/(float)c;
                d_up[totalUp] = (nowUp/(float)N)/(float)c + (c-1)*d_up[totalUp]/(float)c;



                //         res[blockIdx.y] = counter/float(t+1) + t*res[blockIdx.y]/float(t+1);

            for (int i=0;i<N;i++)
                states2[b * N + i] = states[b * N + i];

            }

    
        //res[t * gridDim.y + blockIdx.y] = counter;
      //  if (t==0)
  //          res[blockIdx.y] = counter;
     //   else
   //         res[blockIdx.y] = counter/float(t+1) + t*res[blockIdx.y]/float(t+1);
        }
}

__global__ void block_sum(const int *input, int *per_block_results, const size_t n)
{
    extern __shared__ int sdata[];

    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

    // load input into __shared__ memory
    int x = 0;
    if(i < n)
    {
        x = input[i];
    }
    sdata[threadIdx.x] = x;
    __syncthreads();

    // contiguous range pattern
    for(int offset = blockDim.x / 2;
            offset > 0;
            offset >>= 1)
    {
        if(threadIdx.x < offset)
        {
            // add a partial sum upstream to our own
            sdata[threadIdx.x] += sdata[threadIdx.x + offset];
        }

        // wait until all threads in the block hav
        // updated their partial sums
        __syncthreads();
    }

    // thread 0 writes the final result
    if(threadIdx.x == 0)
    {
 //       printf("%d %d\n", blockIdx.x, sdata[0]);
        per_block_results[blockIdx.x] = sdata[0];
    }
}

void initRands(dim3 threadGrid, int numBlocks, curandState *state, int seed) 
{
    d_initRands<<< numBlocks, threadGrid >>>(state, seed);
    if (cudaSuccess != cudaGetLastError()) printf( "cuda error!\n" );
}
void advanceTimestep(dim3 threadGrid, int numBlocks, curandState *rands, float* OU, float* wg, int* states, int N_x, float sw)
{
    d_generateOU<<< numBlocks, threadGrid >>>(rands, OU, wg);
    if (cudaSuccess != cudaGetLastError()) printf( "cuda error!\n" );

    d_updateStates<<< numBlocks, threadGrid >>>(states, OU, wg, N_x, rands, sw);
    if (cudaSuccess != cudaGetLastError()) printf( "cuda error!\n" );

}
void countStates(int numThreads, int numBlocks, int* states, int* blockTotals, int N_ALL)
{
    block_sum<<< numBlocks, numThreads, numThreads * sizeof(int) >>>(states, blockTotals, N_ALL);
    if (cudaSuccess != cudaGetLastError()) printf( "cuda error!\n" );

}
void recordData(dim3 threadGrid, int numBlocks, int* states, int* states2, int N_x, float* d_up, float* d_down, int* d_upcount, int* d_downcount, int t)
{
     d_recordData<<< numBlocks, threadGrid >>>(states, states2, N_x, d_up, d_down, d_upcount, d_downcount, t);
     if (cudaSuccess != cudaGetLastError()) printf( "cuda error!\n" );
}
