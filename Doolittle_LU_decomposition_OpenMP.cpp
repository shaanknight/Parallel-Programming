#include <iostream>
#include <chrono> 
#include <algorithm>
#include <omp.h> 
#include <math.h>
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

// utility function to print a matrix
void trace_matrix(int n,double mat[sz][sz])
{
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<n;++j)
			cout << mat[i][j] << " ";
		cout << "\n";
	}
}

// utility function to multiply 2 n*n matrices
void multiply(int n,int m,int l,double mat1[sz][sz],double mat2[sz][sz],double prd[sz][sz])
{
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<l;++j)
		{
			prd[i][j] = 0;
			for(int k=0;k<m;++k)
			{
				prd[i][j] += mat1[i][k]*mat2[k][j];
			}
		}
	}
}

void par_multiply(int n,int m,int l,double mat1[sz][sz],double mat2[sz][sz],double prd[sz][sz])
{
	int i,j,k;
	#pragma omp parallel shared (mat1, mat2, prd, n) private(i,j,k) 
	{
		#pragma omp for schedule(static)
		for(i=0;i<n;++i)
		{
			for(j=0;j<l;++j)
			{
				prd[i][j] = 0;
				for(k=0;k<m;++k)
				{
					prd[i][j] += mat1[i][k]*mat2[k][j];
				}
			}
		}
	}
}

double ltm[sz][sz], utm[sz][sz];

void mylu(int n,double mat[sz][sz])
{
	for(int i=0;i<n;++i)
		for(int j=0;j<n;++j)
			ltm[i][j] = utm[i][j] = 0;

	for(int i=0;i<n;++i)
	{
		// filling ith column of lower triangular matrix
		for(int j=i;j<n;++j)
		{
			ltm[j][i] = mat[j][i];
			for(int k=0;k<i;++k)
			{
				ltm[j][i] -= ltm[j][k]*utm[k][i];
			}
		}

		// filling ith row of upper triangular matrix

		utm[i][i] = 1;
		for(int j=i+1;j<n;++j)
		{
			utm[i][j] = mat[i][j];
			for(int k=0;k<i;++k)
			{
				utm[i][j] -= ltm[i][k]*utm[k][j];
			}
			utm[i][j] /= ltm[i][i];
		}
	}

	// overwriting the matrix A
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<=i;++j)
			mat[i][j] = ltm[i][j];
		for(int j=i+1;j<n;++j)
			mat[i][j] = utm[i][j];
	}
}

void par_lu(int n,double mat[sz][sz])
{
	int i,j,k;
	#pragma omp parallel shared (mat, ltm, utm, n) private(i,j,k) 
	{
		#pragma omp for schedule(static)
		for(i=0;i<n;++i)
		{
			for(int j=0;j<n;++j)
				ltm[i][j] = utm[i][j] = 0;
		}

		for(i=0;i<n;++i)
		{
			// filling ith column of lower triangular matrix
			#pragma omp for schedule(static)
			for(j=i;j<n;++j)
			{
				ltm[j][i] = mat[j][i];
				for(k=0;k<i;++k)
				{
					ltm[j][i] -= ltm[j][k]*utm[k][i];
				}
			}

			// filling ith row of upper triangular matrix

			utm[i][i] = 1;
			#pragma omp for schedule(static)
			for(j=i+1;j<n;++j)
			{
				utm[i][j] = mat[i][j];
				for(k=0;k<i;++k)
				{
					utm[i][j] -= ltm[i][k]*utm[k][j];
				}
				utm[i][j] /= ltm[i][i];
			}
		}

		// overwriting the matrix A
		#pragma omp for schedule(static)
		for(i=0;i<n;++i)
		{
			for(j=0;j<=i;++j)
				mat[i][j] = ltm[i][j];
			for(j=i+1;j<n;++j)
				mat[i][j] = utm[i][j];
		}
	}
}

// forward sweep
void forward(int n,double mat[sz][sz],double b[sz][sz],double y[sz][sz])
{
	for(int i=0;i<n;++i)
	{
		double val = b[i][0];
		for(int j=0;j<i;++j)
			val -= mat[i][j]*y[j][0];
		y[i][0] = val/mat[i][i];
	}
}

// backward sweep
void backward(int n,double mat[sz][sz],double y[sz][sz],double x[sz][sz])
{
	for(int i=n-1;i>=0;i--)
	{
		double val = y[i][0];
		for(int j=n-1;j>i;j--)
			val -= mat[i][j]*x[j][0];
		x[i][0] = val;
	}
}

int main()
{
	ios::sync_with_stdio(false);
	cin.tie(0);
	cout.tie(0);

	int n,m,i,j;
	int G = 1e4;
	static double mat[sz][sz],tmp[sz][sz],b[sz][sz],y[sz][sz],x[sz][sz],prd[sz][sz];

	srand(6*rand());
	cin >> n >> m;
	for(i=0;i<n;++i)
	{
		for(j=0;j<n;++j)
		{
			mat[i][j] = rand()%G;
			tmp[i][j] = mat[i][j];
		}
	}

	for(i=0;i<n;++i)
		for(j=0;j<n;++j)
			b[i][j] = y[i][j] = x[i][j] = 0;

	for(i=0;i<n;++i)
		b[i][0] = 1ll*rand()*rand()%(G*G);

	// serial implemenation

	auto start = high_resolution_clock::now();
	
	mylu(n,mat);
	forward(n,mat,b,y);
	backward(n,mat,y,x);
	
	multiply(n,n,1,tmp,x,prd);

	double err = 0;
	for(i=0;i<n;++i)
		err += (prd[i][0]-b[i][0])*(prd[i][0]-b[i][0]);
	
	cout << "solution found using serial implemenatation with error = " << sqrt(err) << endl;

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);

	cout << "duration for serial code to run : " << (duration.count()*1e-6) << " seconds" << endl;

	// parallel implementation

	for(i=0;i<n;++i)
		for(j=0;j<n;++j)
			mat[i][j] = tmp[i][j];

	for(i=0;i<n;++i)
		for(j=0;j<n;++j)
			y[i][j] = x[i][j] = prd[i][j]= 0;

	omp_set_num_threads(m);

	start = high_resolution_clock::now();
	
	par_lu(n,mat);
	forward(n,mat,b,y);
	backward(n,mat,y,x);
	
	par_multiply(n,n,1,tmp,x,prd);

	err = 0;
	for(i=0;i<n;++i)
		err += (prd[i][0]-b[i][0])*(prd[i][0]-b[i][0]);
	
	cout << "solution found using parallel implemenatation with error = " << sqrt(err) << endl;

	stop = high_resolution_clock::now();
	duration = duration_cast<microseconds>(stop - start);

	cout << "duration for parallel code to run : " << (duration.count()*1e-6) << " seconds" << endl;

	return 0;
}
	