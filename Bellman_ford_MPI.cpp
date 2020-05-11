/* MPI Program Template */

#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits.h>
#include "mpi.h"
using namespace std;

typedef long long ll;

vector<ll> vec,graph,dist,ar;
string input,output;

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

    ll n,m,source,i,j,div,req;
    
    if(!(rank))
    {
        input = argv[1];
        output = argv[2];
        ifstream fin (input);
        ll t;
        while(fin >> t)
            vec.push_back(t);
        fin.close();
        n = vec[0];
        m = vec[1];
        source = vec.back();
        for(i=2;i<(int) vec.size() - 1;++i)
            graph.push_back(vec[i]);
        div = (m+numprocs-1)/numprocs;
        req = div*numprocs;
        while(m<req)
        {
            ++m;
            graph.push_back(0);
            graph.push_back(0);
            graph.push_back(0);
        }
    }
    
    /*synchronize all processes*/
    MPI_Barrier( MPI_COMM_WORLD );
    double tbeg = MPI_Wtime();

    /* write your code here */
    MPI_Bcast(&n,1,MPI_LONG_LONG_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&m,1,MPI_LONG_LONG_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&source,1,MPI_LONG_LONG_INT,0,MPI_COMM_WORLD);

    ll h = (3*m)/numprocs; // h is the size of chunk each process receive
    ar.resize(h);
    MPI_Scatter(graph.data(),h,MPI_LONG_LONG_INT,ar.data(),h,MPI_LONG_LONG_INT,0,MPI_COMM_WORLD);

    dist.resize(n+1);
    fill(dist.begin(),dist.end(),LONG_MAX);

    dist[source] = 0;
    dist[0] = 0;

    bool chg = 1;

    for(i=1;i<n && chg == 1;++i)
    {
        chg = 0;
        for(j=0;j<(int) ar.size();j+=3)
        {
            if(dist[ar[j]]-ar[j+2] > dist[ar[j+1]])
            {
                chg = 1;
                dist[ar[j]] = dist[ar[j+1]]+ar[j+2];
            }
            if(dist[ar[j+1]]-ar[j+2] > dist[ar[j]])
            {
                chg = 1;
                dist[ar[j+1]] = dist[ar[j]]+ar[j+2];
            }
        }

        MPI_Allreduce(MPI_IN_PLACE,&chg,1,MPI_CXX_BOOL,MPI_LOR,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,dist.data(),n+1,MPI_LONG_LONG_INT,MPI_MIN,MPI_COMM_WORLD);
    }

    if(!(rank))
    { 
        ofstream fout;
        fout.open (output);
        for(int i=1;i<(int) dist.size();++i)
        {
            if(dist[i] == LONG_MAX)
                fout << i << " " << -1 << "\n";
            else
                fout << i << " " << dist[i] << "\n";
        }
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