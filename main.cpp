#include "matrlib.h"


using namespace std;

int main()
{	
	double pi = 3.1415926535;
	matrix Result(11,1);
	for(int j = 2; j<10; j++)
	{	
		int N = pow(2, j);
		double h = 1/double(N-1);
		N = N + 1;
		matrix A(N,1);
	 	matrix B(N,1);
		matrix C(N,1);
		matrix X(N,1);
		matrix Y(N,1);
		matrix Res_SP(N,N);
		A.filmat_h_first(N,0,0,h);
	 	B.filmat_diag(N,1/(h*h),h);
	 	C.filmat_h_first(N,1,0,h);
	 	for(int m = 1; m < N-1; m++)
	 	{
	 		for(int k = 0; k < N; k++)
	 			{
					X.SetMij(k, 0, sin(pi * double(m * (2*k-1))/double(2*(N-2)))); // N - 2 т.к. N -> N+1 в начале 
				}
			smart_mult_mat(A,B,C,X,Y);
			for(int k = 0; k < N; k++)
			{
				Res_SP.SetMij(k,m,X.GetMij(k,0));
			}
			double L = 4 * sin(pi*double(m)/double(2*(N-2))) * sin(pi*double(m)/double(2*(N-2)))/(h*h) ;
	 		mat_mult_scalar(L, X);
			for(int k = 1; k < N-1; k++)
			{
				double res = abs(X.GetMij(k,0) - Y.GetMij(k,0))/abs(L);
				if(Result.GetMij(j-2,0) < res)
				{
					Result.SetMij(j-2,0,res);
				}
			}
	 	}
		double Sum = 0;
		for(int i = 0; i < N; i++)
		{
			for(int l = i + 1; l < N; l++)
			{
				double S = 0;
				for(int k = 0; k < N; k++)
				{
					if(k==0 || k ==1 || k==N-2 || k == N-1)
						S += 0.5*Res_SP.GetMij(k,i)*Res_SP.GetMij(k,l);
					else 
						S += Res_SP.GetMij(k,i)*Res_SP.GetMij(k,l);
				}
				if(S > Sum)
					Sum = S; 
			}
		}
		cout << "N = 2^" << j << "   scal. prod.  = " << Sum << endl;
		cout << "precision" << Result.GetMij(j-2,0) << endl;  
	}
	return 0;
}
