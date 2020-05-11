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

void mylu(int n,double mat[sz][sz])
{
	double ltm[sz][sz];
	double utm[sz][sz];

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

	 
	// cout << "printing lower triangular matrix" << "\n";
	// trace_matrix(n,ltm);
	// cout << "printing upper triangular matrix" << "\n";
	// trace_matrix(n,utm);

	// double prd[sz][sz];
	// multiply(n,n,n,ltm,utm,prd);
	// for(int i=0;i<n;++i)
	// {
	// 	for(int j=0;j<n;++j)
	// 	{
	// 		if(abs(prd[i][j]-mat[i][j]) > 1e-3)
	// 		{
	// 			cout << "error in LU decomposition" << "\n";
	// 			return;
	// 		}
	// 	}
	// }
	

	// overwriting the matrix A
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<=i;++j)
			mat[i][j] = ltm[i][j];
		for(int j=i+1;j<n;++j)
			mat[i][j] = utm[i][j];
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

	int n,i,j;
	double mat[sz][sz],tmp[sz][sz];

	cout << "enter size of matrix A : " << endl;
	cin >> n;

	cout << "enter the matrix A : " << endl;
	for(i=0;i<n;++i)
	{
		for(j=0;j<n;++j)
		{
			cin >> mat[i][j];
			tmp[i][j] = mat[i][j];
		}
	}

	mylu(n,mat);
	// cout << "Done with LU decomposition" << endl;

	// solving Ax = b

	double b[sz][sz],y[sz][sz],x[sz][sz],prd[sz][sz];
	
	for(i=0;i<n;++i)
		for(j=0;j<n;++j)
			b[i][j] = y[i][j] = x[i][j] = 0;

	cout << "enter the matrix b : " << endl;

	for(i=0;i<n;++i)
		cin >> b[i][0];

	forward(n,mat,b,y);
	backward(n,mat,y,x);

	cout << "corresponding solution to the system of equations : " << endl;
	for(i=0;i<n;++i)
		cout << x[i][0] << " ";
	cout << "\n";

	multiply(n,n,1,tmp,x,prd);
	double err = 0;
	for(i=0;i<n;++i)
		err += (prd[i][0]-b[i][0])*(prd[i][0]-b[i][0]);
	
	cout << "solution found with error = " << err << endl;
	return 0;
}
	