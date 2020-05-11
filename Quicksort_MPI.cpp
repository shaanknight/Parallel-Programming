/* MPI Program Template */

#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "mpi.h"
using namespace std;

typedef long long ll;

vector<ll> vec,ar,rec,tmp;
string input,output;

void insertion_sort(int low,int high)
{
    for(int i = low+1; i<=high; ++i)
    {
        ll val = ar[i];
        int j = i-1;

        while(j >= low && val < ar[j])
        {
            ar[j+1] = ar[j];
            j--;
        }

        ar[j+1] = val;
    }
    return;
}

ll choose_pivot(int l,int h)
{
    int p = rand()%(h-l+1)+l;
    int q = rand()%(h-l+1)+l;
    int r = rand()%(h-l+1)+l;

    if(ar[q] < ar[p])
        swap(ar[p],ar[q]);
    if(ar[r] < ar[p])
        swap(ar[p],ar[r]);
    if(ar[q] < ar[r])
        swap(ar[q],ar[r]);

    return ar[r];
}

pair<int,int> partition(int low,int high,ll pivot)
{
    int j = low;
    while(j<=high)
    {
        if(ar[j] < pivot)
        {
            swap(ar[j],ar[low]);
            ++low;
            ++j;
        }
        else if(ar[j] == pivot)
            ++j;
        else
        {
            swap(ar[j],ar[high]);
            high--;
        }
    }

    return {low,high};
}

void quicksort(int low,int high)
{
    if(high <= low)
        return;

    if(high-low <= 20)
    {
        insertion_sort( low, high);
        return;
    }

    ll pivot = choose_pivot( low, high);
    
    pair<int,int> lh = partition( low, high, pivot);
    
    quicksort( low, lh.first - 1);
    quicksort( lh.second + 1, high);
    
    return;
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
    int n;
    /*  partially partitioning the data here
        keeping a small chunk of data separately in the rank 0 process and sorting it
        later and merging it with the parallely sorted significant array, 
        the small chunk of data is of size n%numprocs and cut off from the end of array
    */
    if(!(rank))
    {
        input = argv[1];
        output = argv[2];
        ifstream fin (input);
        ll t;
        while(fin >> t)
            vec.push_back(t);
        n = vec.size();
        n -= n%numprocs;
        fin.close();
    }
    
    /*synchronize all processes*/
    MPI_Barrier( MPI_COMM_WORLD );
    double tbeg = MPI_Wtime();

    /* write your code here */
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
    int m = n/numprocs; // m is the size of chunk each process receive
    ar.resize(m);
    MPI_Scatter(vec.data(),m,MPI_LONG_LONG_INT,ar.data(),m,MPI_LONG_LONG_INT,0,MPI_COMM_WORLD);
    quicksort(0,m-1);

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
            
            // trace(ar);
            // trace(rec);

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
            ar.insert(ar.end(),vec.begin()+n,vec.end());
            quicksort(n,(int) ar.size() - 1);
            
            tmp.resize(sz);
            merge(ar.begin(), ar.begin() + n, ar.begin() + n, ar.end(), tmp.begin());
            ar = tmp;
            tmp.clear();
        }
        ofstream fout;
        fout.open (output);
        for(int i=0;i<(int) ar.size();++i)
            fout << ar[i] << " ";
        fout << "\n";
        fout.close();
    }

    MPI_Barrier( MPI_COMM_WORLD );
    double elapsedTime = MPI_Wtime() - tbeg;
    double maxTime;
    MPI_Reduce( &elapsedTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    if ( rank == 0 ) {
        printf( "Total time (s): %f\n", maxTime );
    }

    /* shut down MPI */
    MPI_Finalize();
    return 0;
}