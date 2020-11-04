#include <iostream>
#include <chrono>
#include "graph.h"
#include "MyStack.h"
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
    auto duration = duration_cast<milliseconds>(stop - start);

    //template <file_vertex_t, file_index_t, file_weight_t
    //new_vertex_t, new_index_t, new_weight_t>
    graph<long, long, /*int*/ long, long, long, /* char*/ long>
        *ginst = new graph<long, long, /*int*/ long, long, long, /*char*/ long>
        (beg_file, csr_file, weight_file);

    int numofVertex = ginst->vert_count;
    int numofEdges = ginst->edge_count;

    for (int i = 0; i < numofVertex; i++)
        cout<<ginst->beg_pos[i]<<",";
    cout<<"\n";

    for (int i = 0; i < numofEdges; i++)
        cout<<ginst->csr[i]<<",";
    cout<<"\n";
    
    int **Adjacency = new int *[numofVertex];
    CreateAdjacency(Adjacency, numofVertex, ginst);
    printAdj(Adjacency, numofVertex);
// std::cout<<time1<<"\n"<<time2;


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
