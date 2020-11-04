#include <iostream>
#include <chrono>
#include "graph.h"
#include "MyStack.h"
#include <stdio.h>
#include <unistd.h>
using namespace std;
using namespace std::chrono;
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define White 0
#define Gray 1
#define Black 2
#define INF (1 << 31) - 1
#define bitMapSize 10000000
#define StackSize 1000

#define INF 0x3f3f3f3f

void printlist(int *a, int size)
{
    // int size = sizeof(&a)/sizeof(a[0]); only use for array() not for int*
    // std::cout << "sizeof a :" << size << std::endl;
    for (int i = 2; i < size; i++)
    {
        std::cout << a[i] << " ";
    }
    std::cout << std::endl;
}

void printAdj(int **Adjacency, int numofVertex)
{
    // int size = sizeof(&a)/sizeof(a[0]); only use for array() not for int*
    std::cout << "numofVertex a :" << numofVertex << std::endl;

    for (int i = 0; i < numofVertex; i++)
    {
        cout<<i<<" -> ";
        int size = Adjacency[i][1];
        printlist(Adjacency[i], size + 2);
    }
    std::cout << std::endl;
}
void CreateAdjacency(int **Adjacency, int numofVertex, graph<long, long, /*int*/ long, long, long, /* char*/ long> *ginst)
{
    Adjacency[3] = new int[5];
    std::cout << ginst->vert_count << " size CreateAdjacency...\n";
    for (int i = 0; i < ginst->vert_count; i++)
    {
        int beg = ginst->beg_pos[i];
        int end = ginst->beg_pos[i + 1];
        int numofneighbor = end - beg;
        int size = numofneighbor + 2;
        // std::cout << i << " size= " << size<< "'s neighor list: \n";
        Adjacency[i] = new int[size];
        // std::cout<<" i "<<i<<"  ";
        Adjacency[i][0] = 0;
        Adjacency[i][1] = numofneighbor; //default color is 0

        if (numofneighbor > 0)
        {
            for (int j = beg; j < end; j++)
            {
                Adjacency[i][j - beg + 2] = ginst->csr[j];
                // std::cout << "j= " << j << " " << ginst->csr[j] << " ";
            }
        }
        // std::cout << "\n";
        // printlist(Adjacency[i], size);
    }
    std::cout << " \nCreateAdjacency finish \n";
}

void initAdj(int **&Adjacency, int numofVertex)
{
    std::cout << "initAdj... numofVertex :" << numofVertex << std::endl;

    for (int i = 0; i < numofVertex; i++)
    {
        Adjacency[i][0] = 0;
        // printlist(Adjacency[i], size+2);
    }
    std::cout << std::endl;
}
//-----------------------------------------------------------------------------
// MurmurHash3 was written by Austin Appleby, and is placed in the public
// domain. The author hereby disclaims copyright to this source code.
// Note - The x86 and x64 versions do _not_ produce the same results, as the
// algorithms are optimized for their respective platforms. You can still
// compile and run any of them on any platform, but your performance with the
// non-native version will be less than optimal.

__forceinline__  __host__ __device__ uint32_t rotl32( uint32_t x, int8_t r ) {
    return (x << r) | (x >> (32 - r));
  }
  
  __forceinline__ __host__ __device__ uint32_t fmix32( uint32_t h ) {
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
  }
  
  __forceinline__  __host__ __device__ uint32_t hash_murmur(const int64_t& key) {
  
    constexpr int len = sizeof(int64_t);
    const uint8_t * const data = (const uint8_t*)&key;
    constexpr int nblocks = len / 4;
    uint32_t h1 = 0;
    constexpr uint32_t c1 = 0xcc9e2d51;
    constexpr uint32_t c2 = 0x1b873593;
    //----------
  
    // body
    const uint32_t * const blocks = (const uint32_t *)(data + nblocks*4);
    for(int i = -nblocks; i; i++)
    {
      uint32_t k1 = blocks[i];
      k1 *= c1;
      k1 = rotl32(k1,15);
      k1 *= c2;
      h1 ^= k1;
      h1 = rotl32(h1,13); 
      h1 = h1*5+0xe6546b64;
    }
    //----------
    // tail
    const uint8_t * tail = (const uint8_t*)(data + nblocks*4);
    uint32_t k1 = 0;
    switch(len & 3)
    {
      case 3: k1 ^= tail[2] << 16;
      case 2: k1 ^= tail[1] << 8;
      case 1: k1 ^= tail[0];
              k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1;
    };
    //----------
    // finalization
    h1 ^= len;
    h1 = fmix32(h1);
    return h1;
  }
//-----------------------------------------------------------------------------
//end of hash

__global__
void addinHash(int** Adjacency_partitioned, long* rowcompress, 
    long* dst, int* d_hash, int d_hashSize,
    int numofVertex,int numofEdges){
        
        for (int i = blockIdx.x * blockDim.x + threadIdx.x; 
            i < numofVertex; 
            i += blockDim.x * gridDim.x) 
         {

            if(i%2 !=0){
                int64_t pair =((int64_t)(i)<<32L);
                uint32_t hval = (hash_murmur(pair)%d_hashSize);
                d_hash[hval]=1;
                // printf(" %d in hash %d, ",i, hval);
            }
         }
}

__global__
void merge(int** Adjacency_partitioned, long* rowcompress, 
    long* dst, int* d_hash, int d_hashSize,
    int numofVertex,int numofEdges, long* d_rowcompress_new, long*d_dst_new){

    // if( threadIdx.x==0) {
        // for (int i = 0; i < numofVertex+1; i++)
        // printf(" %d, ",rowcompress[i]);
        // printf("\n");
        // for (int i = 0; i < numofEdges; i++)
        // printf(" %d, ",dst[i]);
        // printf("\n");
        // for (int i = 0; i < numofEdges; i++)
        // printf(" %d, ",d_dst_new[i]);
        // printf("\n");
        // }

        for (int i = blockIdx.x * blockDim.x + threadIdx.x; 
            i < numofEdges; 
            i += blockDim.x * gridDim.x) 
         {
            int64_t pair =((int64_t)(i)<<32L);
            uint32_t hval = (hash_murmur(pair)%d_hashSize);
            if(d_hash[hval]==0){
                d_dst_new[i]=dst[i];
               
            }
            // else{
            //     printf("%d merged ",i);
            // }
         }
}


int main(int args, char **argv)
{
    std::cout << "Input: ./exe beg csr weight source(optional) report(optional)\n";
    if (args >= 7)
    {
        std::cout << "Wrong input\n";
        return -1;
    }

    int source = 0;
    bool report_nb_vertices_visited = false;

    const char *beg_file = argv[1];
    const char *csr_file = argv[2];
    const char *weight_file = argv[3];

    if (args >= 5)
        source = atoi(argv[4]);
    if (args >= 6)
        report_nb_vertices_visited = atoi(argv[5]);

    auto start = high_resolution_clock::now();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    //template <file_vertex_t, file_index_t, file_weight_t
    //new_vertex_t, new_index_t, new_weight_t>
    graph<long, long, /*int*/ long, long, long, /* char*/ long>
        *ginst = new graph<long, long, /*int*/ long, long, long, /*char*/ long>
        (beg_file, csr_file, weight_file);

    int numofVertex = ginst->vert_count;
    int numofEdges = ginst->edge_count;

    // for (int i = 0; i < numofVertex+1; i++)
    //     cout<<ginst->beg_pos[i]<<",";
    // cout<<"\n";

    // for (int i = 0; i < numofEdges; i++)
    //     cout<<ginst->csr[i]<<",";
    // cout<<"\n";
    
    int **Adjacency = new int *[numofVertex];
    CreateAdjacency(Adjacency, numofVertex, ginst);

    // printAdj(Adjacency, numofVertex);


    int **Adjacency_partitioned = new int *[numofVertex];
    long* rowcompress=ginst->beg_pos;
    long* dst=ginst->csr;



    int* d_hash;
    const int HASH_TABLE_SIZE = 10000000;  int d_hashSize=HASH_TABLE_SIZE;
    printf("start malloc\n");

    cudaMalloc(&d_hash, HASH_TABLE_SIZE*sizeof(int));
    cudaMemset(d_hash,0,HASH_TABLE_SIZE*sizeof(int));


    int thread=256;
    int block = numofVertex/256 +1;
    cout<<block<<" "<<thread<<endl;
    addinHash<<<block,thread>>> (Adjacency_partitioned,rowcompress,dst,
                                d_hash,d_hashSize, numofVertex,numofEdges);

    
    ////-------------------------test addinHash success
    // int count=0; 
    // int* h_hash=(int*)malloc(HASH_TABLE_SIZE * sizeof(int));
    // cudaMemcpy(h_hash,d_hash,HASH_TABLE_SIZE*sizeof(int),cudaMemcpyDeviceToHost);
    // for(int i=0; i<d_hashSize; i++){
    //     if (h_hash[i]==1)
    //     count++;
    // }
    // cout<<count;
    ////-------------------------test addinHash success

    long* d_rowcompress; long*  d_dst;
    
    cudaMalloc(&d_rowcompress, (numofVertex+1)*sizeof(long));
    cudaMalloc(&d_dst, (numofEdges)*sizeof(long));
    cudaMemcpy(d_rowcompress,rowcompress,(numofVertex+1)*sizeof(long),cudaMemcpyHostToDevice);
    cudaMemcpy(d_dst,dst,numofEdges*sizeof(long),cudaMemcpyHostToDevice);

    /// output
    long* d_rowcompress_new; long*  d_dst_new;
    cudaMalloc(&d_rowcompress_new, (numofVertex+1)*sizeof(long));
    cudaMalloc(&d_dst_new, (numofEdges)*sizeof(long));

    cudaMemset(d_rowcompress_new, INF ,(numofVertex+1)*sizeof(long));
    cudaMemset(d_dst_new, INF , (numofEdges)*sizeof(long) );



    block = numofEdges/256 +1; //numofEdges
    cout<<block<<" "<<thread<<endl;


    start = high_resolution_clock::now();

    merge<<<block,thread>>> (Adjacency_partitioned,d_rowcompress, d_dst,
        d_hash,d_hashSize, numofVertex,numofEdges,d_rowcompress_new,d_dst_new);
        cudaDeviceSynchronize();

    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start);
    std::cout << "\n   merge<<<block,thread>>>  " << duration.count() << "\n ";

    auto time1=duration.count();

    //-------------------------test merge success
    int count=0; 
    int* h_dst_new=(int*)malloc(numofEdges*sizeof(long));
    cudaMemcpy(h_dst_new,d_dst_new,numofEdges*sizeof(long),cudaMemcpyDeviceToHost);
    for(int i=0; i<numofEdges; i++){
        
        if (h_dst_new[i]!=1061109567 && h_dst_new[i]>0  ){
            // printf("%d ",h_dst_new[i]);
            count++;
        }
    }

    for(int i=0; i<(numofEdges<100?numofEdges:100); i++)
        cout<<" "<<h_dst_new[i]<<" ";
        cout<<"\n";

    cout<<"\n"<<numofEdges<<endl;
    cout<<count<<endl;
    cout<<time1;
    //-------------------------test merge success

// cudaMalloc(&d_src, HASH_TABLE_SIZE*sizeof(vid_t));
// cudaMalloc(&d_dest, HASH_TABLE_SIZE*sizeof(vid_t));
// cudaMalloc(&d_hash, HASH_TABLE_SIZE*sizeof(int));














































    return 0;
}




//     int **Adjacency = new int *[numofVertex];
//     CreateAdjacency(Adjacency, numofVertex, ginst);
//     Adjacency[0][0] = 1;
//     Adjacency[1][0] = 1;
//     Adjacency[2][0] = 1;
//     Adjacency[3][0] = 1;
//     // printAdj(Adjacency, numofVertex);
//     initAdj(Adjacency, numofVertex);
//     // printAdj(Adjacency, numofVertex);

//     start = high_resolution_clock::now();
//     ///////////
//     // int *visit = dfs_by_vertexid_array(Adjacency, 0, numofVertex);

//     int *visit = dfs_by_vertexid_array(Adjacency, source, numofVertex, report_nb_vertices_visited);
//     // printlist(visit, numofVertex);
//     ///////////
//     stop = high_resolution_clock::now();
//     duration = duration_cast<milliseconds>(stop - start);
//     std::cout << "\nAverage time SC15 dfs_by_vertexid_array time spends (ms): \n"
//               << duration.count() << endl;
// float time1=duration.count(), time2;

//     initAdj(Adjacency, numofVertex);
//     start = high_resolution_clock::now();
//     ///////////
//     // DFS_optimize5_removeGraphStruct(Adjacency, 0, numofVertex); DFS_optimize5_removeGraphStruct
//     DFS_optimize5_removeGraphStruct(Adjacency, source, numofVertex, report_nb_vertices_visited);
//     ///////////
//     stop = high_resolution_clock::now();
//     duration = duration_cast<milliseconds>(stop - start);
//     std::cout << "\n DFS_optimize5_removeGraphStruct time spends (ms): \n" 
//     << duration.count() << endl;
// time2=duration.count();
//     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DFS_visit_removeGraph(int **Adjacency, int start, int &NumDiscover)
{ //   //7/15/2020 removeGraph from DFS_visit_continue_removeSTLstack
    Stack<int> stack(1000);

    if (Adjacency[start][0] == White)
    {
        NumDiscover++;
    }
    if (Adjacency[start][0] == Black)
    {
        return;
    }

    Adjacency[start][0] = Gray;
    stack.push(start);
    while (!stack.isEmpty())
    {
        int currValue = 0;
        stack.pop(&currValue);
        int currColor = Adjacency[start][0];

        ///////////////////////////////////////// ///////////////////////////////////////// TODO pop and push has problem
        //cout << currValue<< "poped, stack.size() " << stack.size() << endl;
        //cout << " at: " << currValue;
        if (Adjacency[currValue][1] == 0)
        {
            continue;
        }
        for (int i = 2; i < Adjacency[currValue][1] + 2; i++)
        {
            //Vertex u = g.vertexs[neighbors[i]];
            int Uvalue = Adjacency[currValue][i];
            int Ucolor = Adjacency[Uvalue][0];
            //int Uvalue = g.vertexs[g.vertexs[currValue].neighbors[i]].value;
            //if (i == neighbors.size() - 1 && u.color == Gray)                ///////////////////////////////////////// never  called on toy graph
            if (i == Adjacency[currValue][1] - 1 && Ucolor == Gray)
            {
                //g.vertexs[currValue].color = Black;
                Adjacency[currValue][0] = Black;
                //cout << " Black: " << currValue;              /////////////////////////////////////////
            }
            if (Ucolor == White)
            {
                //cout << " discover: " << Uvalue;                  /////////////////////////////////////////
                //Ucolor = Gray;
                //set_color(g, u, 1);
                Adjacency[Uvalue][0] = Gray;
                stack.push(Uvalue);
                NumDiscover++;
                continue; // if continue search all neighbors, became  DFS break
            }
        }
        //cout << "\n";                                  /////////////////////////////////////////
    }
}

void DFS_optimize5_removeGraphStruct(int **Adjacency, int start, int numofVertex, bool report_nb_vertices_visited)
{                        //7/15/2020 removeGraphStruct, other same as optimize4
    int NumDiscover = 0; //////////use int* in mutil thread or use
    for (int i = start; i < numofVertex; i++)
    {
        if (Adjacency[i][0] == White)
        {
            DFS_visit_removeGraph(Adjacency, i, NumDiscover);
        }
        if (NumDiscover >= numofVertex)
        {
            if (report_nb_vertices_visited)
                cout << "DFS_optimize5_removeGraphStruct NumDiscover: " << NumDiscover << endl;
            return;
        }
    } //cout << "g.NumDiscover: " << g.NumDiscover << endl;
    if (report_nb_vertices_visited)
        cout << "DFS_optimize5_removeGraphStruct NumDiscover: " << NumDiscover << endl;
}

template <class Value>
Value *mynew_array(size_t nb)
{
    Value *res = (Value *)malloc(size_t(sizeof(Value)) * nb);
    if (res == NULL)
        cout << "mynew_array returned NULL";
    return res;
}
template <class Number, class Size>
void fill_array_seq(Number *array, Size sz, Number val)
{
    memset(array, val, sz * sizeof(Number));
    // for (Size i = Size(0); i < sz; i++)
    //   array[i] = val;
}
static inline void myfree(void *p)
{
    free(p);
}

int *dfs_by_vertexid_array(int **Adjacency, int start, int numofVertex, bool report_nb_vertices_visited)
{
    long nb_vertices_visited = 1;
    typedef int vtxid_type;
    int *visited;
    vtxid_type nb_vertices = numofVertex;
    auto source = start;
    visited = mynew_array<int>(nb_vertices);
    fill_array_seq(visited, nb_vertices, 0); // init all color as 0, means unvisited
    cout << "finish init! ALGO_PHASE:\n";
    vtxid_type *frontier = mynew_array<vtxid_type>(nb_vertices);
    vtxid_type frontier_size = 0;
    frontier[frontier_size++] = source;
    visited[source] = 1;

    while (frontier_size > 0)
    {
        vtxid_type vertex = frontier[--frontier_size];
        vtxid_type degree = Adjacency[vertex][1];
        // vtxid_type *neighbors = graph.adjlists[vertex].get_out_neighbors();
        // cout<<"\nat: "<<vertex<<" degree"<<degree<<" visit: ";

        for (vtxid_type edge = 2; edge < degree + 2; edge++)
        {
            // cout<<Adjacency[vertex][edge]<<" ";
            vtxid_type other = Adjacency[vertex][edge];
            if (visited[other])
                continue;
            if (report_nb_vertices_visited)
                (nb_vertices_visited)++;

            visited[other] = 1;
            frontier[frontier_size++] = other;
        }
    }
    if (report_nb_vertices_visited)
        cout << "dfs_by_vertexid_array nb_vertices_visited： " << nb_vertices_visited << endl;
    myfree(frontier);
    return visited;
}
