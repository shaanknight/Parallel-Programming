#include<bits/stdc++.h>
using namespace std;

#define ll long long int
#define ld long double
#define ff first
#define ss second
#define pb push_back
#define pi pair<ll,ll>
#define all(X) X.begin(),X.end()

const int M = (1<<20)+5;
const int md = 1e9+7;
const int sz = 200;

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

void hessenberg_LU_decompose(int n,double mat[sz][sz])
{
	double ltm[sz][sz];
	double utm[sz][sz];

	for(int i=0;i<n;++i)
		for(int j=0;j<n;++j)
			ltm[i][j] = utm[i][j] = 0;

	for(int i=0;i<n;++i)
	{
		ltm[i][i] = 1;
		utm[0][i] = mat[0][i];
	}

	for(int i=1;i<n;++i)
	{
		ltm[i][i-1] = mat[i][i-1]/utm[i-1][i-1];
		for(int j=i;j<n;++j)
			utm[i][j] = mat[i][j]-ltm[i][i-1]*utm[i-1][j];
	}

	cout << "printing lower triangular matrix" << "\n";
	trace_matrix(n,ltm);
	cout << "printing upper triangular matrix" << "\n";
	trace_matrix(n,utm);

	double prd[sz][sz];
	multiply(n,n,n,ltm,utm,prd);
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<n;++j)
		{
			if(abs(prd[i][j]-mat[i][j]) > 1e-3)
			{
				cout << "error in LU decomposition" << "\n";
				return;
			}
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


int main()
{
	ios::sync_with_stdio(false);
	cin.tie(0);
	cout.tie(0);

	int n,i,j;
	double mat[sz][sz],tmp[sz][sz];

	cout << "enter size of Hessenberg matrix A : " << endl;
	cin >> n;

	cout << "enter the Hessenberg matrix A : " << endl;
	for(i=0;i<n;++i)
	{
		for(j=0;j<n;++j)
		{
			cin >> mat[i][j];
			tmp[i][j] = mat[i][j];
		}
	}

	hessenberg_LU_decompose(n,mat);
	return 0;
}
	