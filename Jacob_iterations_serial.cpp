#include <iostream>
#include <chrono> 
#include <algorithm>
#include <omp.h> 
#include <math.h>
#include <vector>
using namespace std;
using namespace std::chrono; 

#define ll long long int
#define ld long double
#define ff first
#define ss second
#define pb push_back
#define pi pair<ll,ll>
#define all(X) X.begin(),X.end()

const int M = (1<<20)+5;
const int md = 1e9+7;
const int sz = 205;

double h[sz][sz],prv[sz][sz],org[sz][sz];

int serial_simulation()
{
	int iters = 0,i,j;

	for(i=1;i<=100;++i)
		for(j=1;j<=100;++j)
			prv[i][j] = h[i][j];

	while(1)
	{
		++iters;
		double err = 0;
		for(i=2;i<=99;++i)
		{
			for(j=2;j<=99;++j)
			{
				h[i][j] = 0.25*(prv[i-1][j] + prv[i+1][j] + prv[i][j-1] + prv[i][j+1]);
				err = max(err,abs(prv[i][j]-h[i][j]));
			}
		}

		if(err < 0.01)
			return iters;

		for(i=1;i<=100;++i)
			for(j=1;j<=100;++j)
				prv[i][j] = h[i][j];
	}
}

int main()
{
	int i,j;

	for(i=1;i<=100;++i)
	{
		for(j=1;j<=100;++j)
		{
			if(i == 1 && j>30 && j<=70)
				h[i][j] = 100;
			else
				h[i][j] = 20;
		}
	}

	for(i=1;i<=100;++i)
	{
		for(j=1;j<=100;++j)
		{
			org[i][j] = h[i][j];
		}
	}

	auto start = high_resolution_clock::now();
	cout << "serial simulation stops in : " << serial_simulation() << "\n";
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << "duration for serial code to run : " << (duration.count()*1e-6) << " seconds" << endl;


	return 0;
}