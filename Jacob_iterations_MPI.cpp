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

vector<int> vec,gr,dist,ar;
string input,output;

void trace(vector<ll> td)
{
    cout << (int) td.size() << "\n";
    for(int i=0;i<(int) td.size();++i)
        cout << td[i] << " ";
    cout << "\n";
}

float prv[105][105],h[105][105];

double simulate(int s,int e,int n)
{
    double err = 0;
    for(int i=s;i<=e;++i)
    {
        if(i == 1 || i == n)
            continue;
        for(int j=2;j<n;++j)
        {
            h[i][j] = 0.25*(prv[i-1][j] + prv[i+1][j] + prv[i][j-1] + prv[i][j+1]);
            err = max(err,(double) abs(prv[i][j]-h[i][j]));
        }
    }
    for(int i=s;i<=e;++i)
    {
        if(i == 1 || i == n)
            continue;
        for(int j=2;j<n;++j)
        {
            prv[i][j] = h[i][j];
        }
    }
    return err;
}

int main( int argc, char **argv ) {
    int rank, numprocs;

    /* start up MPI */
    MPI_Init( &argc, &argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &numprocs );

    /* user input from file */

    int n,m,i,j,div,req;
    
    if(!(rank))
    {
        n = 100;
        m = n*n;
        for(i=1;i<=n;++i)
        {
            for(j=1;j<=n;++j)
            {
                if(i == 1 && j>30 && j<=70)
                    gr.push_back(100);
                else
                    gr.push_back(20); 
            }
        }
        if(n%numprocs)
            div = numprocs-n%numprocs;
        else
            div = 0;
        req = m + n*div;
        while(m<req)
        {
            ++m;
            gr.push_back(0);
        }
    }
    
    /*synchronize all processes*/
    MPI_Barrier( MPI_COMM_WORLD );
    double tbeg = MPI_Wtime();

    /* write your code here */
    MPI_Bcast(&n,1,MPI_LONG_LONG_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&m,1,MPI_LONG_LONG_INT,0,MPI_COMM_WORLD);

    ll ht = m/numprocs; // ht is the size of chunk each process receive
    ar.resize(ht);
    MPI_Scatter(gr.data(),ht,MPI_INT,ar.data(),ht,MPI_INT,0,MPI_COMM_WORLD);

    int p = 0;
    int v = (n+numprocs-1)/numprocs;
    int e = (rank+1)*v;
    int s = e-v+1;
    e = min(e,n);

    for(i=s;i<=e;++i)
    {
        for(j=1;j<=n;++j)
        {
            prv[i][j] = ar[p];
            ++p;
        }
    }

    float buffer[n];
    float buffer2[n];
    float buffer3[n];
    float buffer4[n];

    double err = 1.0;
    int iters = 0;
    
    while(1)
    {
        if(numprocs >= 2)
        {
            MPI_Request request, request2, request3, request4;
            MPI_Status status1,status2;
            
            int right = rank + 1;
            int left = rank - 1;

            for(i=1;i<=n;++i)
                buffer2[i-1] = prv[s][i];
            
            for(i=1;i<=n;++i)
                buffer4[i-1] = prv[e][i];
            
            if (left >= 0 && right < numprocs)
            {
                MPI_Send(buffer2, n, MPI_FLOAT, left, 1, MPI_COMM_WORLD);
                MPI_Recv(buffer, n, MPI_FLOAT, right, 1, MPI_COMM_WORLD, &status1);
                MPI_Send(buffer4, n, MPI_FLOAT, right, 2, MPI_COMM_WORLD);
                MPI_Recv(buffer3, n, MPI_FLOAT, left, 2, MPI_COMM_WORLD, &status2);

                for(i=1;i<=n;++i)
                {
                    prv[e+1][i] = buffer[i-1];
                    prv[s-1][i] = buffer3[i-1];
                }
            }
            else if(left < 0)
            {
                MPI_Recv(buffer, n, MPI_FLOAT, right, 1, MPI_COMM_WORLD, &status1);
                MPI_Send(buffer4, n, MPI_FLOAT, right, 2, MPI_COMM_WORLD);

                for(i=1;i<=n;++i)
                    prv[e+1][i] = buffer[i-1];
            }
            else
            {
                MPI_Send(buffer2, n, MPI_FLOAT, left, 1, MPI_COMM_WORLD);
                MPI_Recv(buffer3, n, MPI_FLOAT, left, 2, MPI_COMM_WORLD, &status1);

                for(i=1;i<=n;++i)
                    prv[s-1][i] = buffer3[i-1];
            }
        }

        ++iters;
        
        err = 0; 
        err = simulate(s,e,n);
        MPI_Allreduce(MPI_IN_PLACE,&err,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        if(err < 0.01)
            break;
    }

    if(!(rank))
    { 
        cout << "mpi simulation stops in : " << iters << "\n";
    }

    MPI_Barrier( MPI_COMM_WORLD );
    double elapsedTime = MPI_Wtime() - tbeg;
    double maxTime;
    MPI_Reduce( &elapsedTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );
    if ( rank == 0 ) {
        printf( "duration for mpi code to run : %f seconds\n", maxTime );
    }

    /* shut down MPI */
    MPI_Finalize();
    return 0;
}