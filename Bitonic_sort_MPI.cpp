/* MPI Program Template */

#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono> 
#include "mpi.h"
using namespace std;
using namespace std::chrono; 

typedef long long ll;

vector<ll> vec,ar,br,rec,tmp;

void bitonic_merge(int start, int end, int dir)
{
    int no_of_elements = end-start+1;
    for(int j = no_of_elements/2; j > 0 ; j = j / 2)
        for(int i = start; i + j <= end ; i++)
            if(dir == (ar[i] > ar[i+j])) swap( ar[i], ar[i+j]);
}

void bitonic_sort(int start, int end)
{
    int no_of_elements = end-start+1;
    for(int j = 2; j <= no_of_elements ; j = j * 2)
    {
            for(int i = 0; i < no_of_elements ; i = i + j)
            {
                if( ( i / j ) % 2 == 0) bitonic_merge( i , i + j - 1 , 1);
                else bitonic_merge( i , i + j - 1 , 0);
            }
    }
}

void pad_and_sort()
{
    int m = (int) ar.size();
    int l = m;
    while(m&(m-1))
    {
        ar.push_back(0);
        ++m;
    }
    bitonic_sort(0,m-1);
    for(int i=m-l;i<m;++i)
        br.push_back(ar[i]);
    ar.clear();
    ar = br;
    br.clear();
}

void trace(vector<ll> td)
{
    cout << (int) td.size() << "\n";
    for(int i=0;i<(int) td.size();++i)
        cout << td[i] << " ";
    cout << "\n";
}

int main( int argc, char **argv ) {
    int rank, numprocs;

    /* start up MPI */
    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &numprocs );

    /* user input from file */
    int n, G = 1e4;
    /*  partially partitioning the data here
        keeping a small chunk of data separately in the rank 0 process and sorting it
        later and merging it with the parallely sorted significant array, 
        the small chunk of data is of size n%numprocs and cut off from the end of array
    */
    if(!(rank))
    {
        n = stoi(argv[1]);
        for(int i=0;i<n;++i)
            vec.push_back(rand()%G + 1);
        // for(auto v:vec)
        //     cout << v << "\n";
        n = vec.size();
        n -= n%numprocs;
    }
    
    /*synchronize all processes*/
    MPI_Barrier( MPI_COMM_WORLD );
    double tbeg = MPI_Wtime();

    /* write your code here */
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    int m = n/numprocs; // m is the size of chunk each process receive
    ar.resize(m);
    MPI_Scatter(vec.data(),m,MPI_LONG_LONG_INT,ar.data(),m,MPI_LONG_LONG_INT,0,MPI_COMM_WORLD);
    pad_and_sort();

    MPI_Status stat;

    int i = 1;
    while(i<numprocs)
    {
        int check = i<<1;
        if(rank%check)
        {
            int sz = (int) ar.size();
            int rec_proc = rank-i;
            MPI_Send(ar.data(),sz,MPI_LONG_LONG_INT,rec_proc,0,MPI_COMM_WORLD);
            break;
        }
        if(rank+i < numprocs)
        {
            int sz;
            if(n >= (rank+check)*m) sz = m*i;
            else sz = n-(rank+i)*m;
            rec.resize(sz);
            int send_proc = rank+i;
            MPI_Recv(rec.data(),sz,MPI_LONG_LONG_INT,send_proc,0,MPI_COMM_WORLD,&stat);
            tmp.resize((int) ar.size() + (int) rec.size());
            merge(ar.begin(), ar.end(), rec.begin(), rec.end(), tmp.begin());
            ar = tmp;
            tmp.clear();
            rec.clear(); 
        }
        i = check;
    }

    if(!(rank))
    {
        int sz = (int) vec.size();
        if(sz%numprocs)
        {
            sort(vec.begin()+n,vec.end());
            tmp.resize(sz);
            merge(ar.begin(), ar.end(), vec.begin() + n, vec.end(), tmp.begin());
            ar = tmp;
            tmp.clear();
        }
        
        // for(int i=0;i<(int) ar.size();++i)
        //     cout << ar[i] << " ";
        // cout << "\n";
    }

    MPI_Barrier( MPI_COMM_WORLD );
    double elapsedTime = MPI_Wtime() - tbeg;
    double maxTime;
    MPI_Reduce( &elapsedTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    if ( rank == 0 ) {
        printf( "duration for parallel code to run : %f\n", maxTime);
        
        tmp = ar;
        ar = vec;
        vec = tmp;
        
        auto start = high_resolution_clock::now();
        pad_and_sort();
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        cout << "duration for serial code to run : " << (duration.count())*(1e-6) << endl;
        
        for(int i=0;i<(int)vec.size();++i)
        {
            if(ar[i] != vec[i])
            {
                cout << "Outputs of serial and parallel did not match." << endl;
                /* shut down MPI */
                MPI_Finalize();
                return 0;
            }
        }
        cout << "Outputs of serial and parallel matched successfully." << endl;
    }

    /* shut down MPI */
    MPI_Finalize();
    return 0;
}