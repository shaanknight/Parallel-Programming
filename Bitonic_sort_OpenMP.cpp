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
const int sz = 2005;

vector<int> ar,tmp,res;

void bitonic_merge(int start, int end, int dir)
{
	int no_of_elements = end-start+1;
	for(int j = no_of_elements/2; j > 0 ; j = j / 2)
		for(int i = start; i + j <= end ; i++)
			if(dir == (ar[i] > ar[i+j])) swap( ar[i], ar[i+j]);
}

// ind : 1 -> parallelize the sort
void bitonic_sort(int start, int end,int ind)
{
	int no_of_elements = end-start+1;
	for(int j = 2; j <= no_of_elements ; j = j * 2)
	{
		if(ind)
		{
			#pragma omp parallel for
			for(int i = 0; i < no_of_elements ; i = i + j)
			{
				if( ( i / j ) % 2 == 0) bitonic_merge( i , i + j - 1 , 1);
				else bitonic_merge( i , i + j - 1 , 0);
			}
		}
		else
		{
			for(int i = 0; i < no_of_elements ; i = i + j)
			{
				if( ( i / j ) % 2 == 0) bitonic_merge( i , i + j - 1 , 1);
				else bitonic_merge( i , i + j - 1 , 0);
			}
		}
	}
}

void trace_array(int start,int end)
{
	for(int i=start;i<end;++i)
		cout << ar[i] << " ";
	cout << "\n";
}

int main()
{
	int n,nt = 6,m,i;
	int G = 1e4;
	srand(0);
	cin >> n >> nt;

	for(i=0;i<n;++i)
		ar.pb(rand()%G + 1);

	m = n;
	while(m&(m-1))
	{
		ar.pb(0);
		++m;
	}

	tmp = ar;
	auto start = high_resolution_clock::now();
	bitonic_sort(0,m-1,0);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "duration for serial code to run : " << (duration.count()*1e-6) << endl;

	res = ar;
	ar = tmp;
	start = high_resolution_clock::now();
	omp_set_num_threads(min(nt,n));
	bitonic_sort(0,m-1,1);
	stop = high_resolution_clock::now();
	duration = duration_cast<microseconds>(stop - start);

	cout << "duration for parallel code to run : " << (duration.count()*1e-6) << endl;

	for(i=m-n;i<m;++i)
	{
		if(ar[i] != res[i])
		{
			cout << "Outputs of serial and parallel did not match." << endl;
			return 0;
		}
	}

	cout << "Outputs of serial and parallel matched successfully." << endl;

	return 0;
}


