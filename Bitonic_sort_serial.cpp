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

vector<int> ar;

void merge(int start, int m, int ord) 
{ 
	if(m == 1)
		return;

    for (int i=start; i<start+(m>>1); ++i) 
        if(ord == (ar[i] > ar[i+(m>>1)])) swap( ar[i], ar[i+(m>>1)]); 
       
  	merge(start, m>>1, ord); 
    merge(start+(m>>1), m>>1, ord); 
} 

void bitonic_sort(int start, int m, int ord) 
{
	if(m == 1)
		return;  
  
    bitonic_sort(start, m>>1, 1); 
    bitonic_sort(start+(m>>1), m>>1, 0);  
    merge(start, m, ord); 
}

void trace_array(int start,int end)
{
	for(int i=start;i<end;++i)
		cout << ar[i] << " ";
	cout << "\n";
}

int main()
{
	int n,m,i;
	int G = 1e4;
	srand(0);
	cin >> n;

	for(i=0;i<n;++i)
		ar.pb(rand()%G + 1);
	auto start = high_resolution_clock::now();
	m = n;
	// trace_array(m-n,m);
	while(m&(m-1))
	{
		ar.pb(0);
		++m;
	}
	// trace_array(m-n,m);
	bitonic_sort(0,m,1);
	
	auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    cout << "duration for serial code to run : " << (duration.count())*(1e-6) << endl;

	return 0;
}


