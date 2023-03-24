#ifndef matrlib
#define matrlib
#include <iostream>
#include <vector>
#include <fstream>
#include <typeinfo>
#include <ctime>
#include<math.h>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>

using namespace std;


class matrix
{
private:
	int n; // строк 
	int m; // столбцов
	double* mat;
public:
	matrix()
	{
		n = 0;
		m = 0;
		mat = nullptr;
	}
	matrix(int m_rows, int m_cols)
	{
		n = m_rows;
		m = m_cols;
		mat = new double[n * m];
		for(int i=0;i<n*m;i++)
			mat[i] = 0;
	}
	matrix(int m_rows)
	{
		n = m_rows;
		m = m_rows;
		mat = new double[n * m];
	}
	int Dim_rows()
	{
		return n;
	}
	int Dim_cols()
	{
		return m;
	}
	double GetMij(int i, int j)
	{
		if ((m >= 0) && (n >= 0))
			return mat[i * m + j];
		else
			return 0;
	}
	void SetMij(int i, int j, double value)
	{
		if ((i < 0) || (i >= n))
			return;
		if ((j < 0) || (j >= m))
			return;
		mat[i * m + j] = value;
	}
	void print_mat()
	{
		for(int i=0; i<m;i++)
		{
			for(int j=0; j<n;j++)
			{
				cout << mat[i*n+j] << ' ';
			}
			cout << '\n';
		}
		cout << endl << endl;
	}
	// zero_ind = индекс где надо ноль
	void filmat_h_first(int n, int zero_ind, double spec_h, double h)
	{
		mat[(n-1)*zero_ind] = spec_h;
		for(int i = zero_ind;i<n+zero_ind-1;i++)
		{
			mat[i]=-1/(h*h);
		}
	}
	void filmat_diag(int n, double spec_h, double h)
	{
		mat[0]= spec_h;
		mat[n-1] = spec_h;
		for(int i = 0;i<n;i++)
		{
			mat[i]= 2 /(h * h);
		}
	}
	void operator=(const matrix& M)
	{
		if ((n > 0) && (m > 0))
		{
			delete[] mat;
		}

		m = M.m;
		n = M.n;
		mat = new double[n * m];

		for (int i = 0; i < m; i++)
			for (int j = 0; j < n; j++)
				mat[i * n + j] = M.mat[i * n + j];
		//return *this;
	}

	void out_mat()
	{

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				//cout << setprecision(3);
				//cout << mat[i * m + j] << ' ';
				printf("%10.7lf ", mat[i * m + j]);
			}
			printf("\n");
		}
		cout << endl << endl;
	}
	int fil_mat_file(string filename)
	{
		ifstream fin;
		fin.open(filename);
		if (!fin.is_open())
		{
			cout << "File was not opened" << endl;
			return -1;
		}
		for (int i = 0; i < n * m ; i++)
		{
			if (!fin.eof())
			{
				fin >> mat[i];
			}
			/*if ((fin.eof()) && (i != (n * m - 1)))
			{
				cout << "Filling error" << endl;
				return -1;
			}
			if ((!fin.eof()) && (i == (n * m - 1)))
			{
				cout << "Filling error" << endl;
				return -1;
			}*/
		}
		fin.close();
		return 0;
	}
	void fil_mat_by_x(const matrix& x)
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				mat[i * m + j] = pow(x.mat[i], j);
	}
	void fil_mat_by_x_and_h(const matrix& x)
	{
		//cout << "here" << endl;
		for(int i = 0;i<n;i++)
		{
			for(int j=1;j<m;j++)
			{
				//cout << i << ' ' << j << ' ' << x.mat[i] << j << endl;
				mat[i*m+j] = pow(x.mat[i], j-1);
			}
		}
		for(int i=0;i<n;i++)
		{
			mat[i*m] = (-1+2*(i%2));
		}
	}
	void fil_mat_by_2d_fun(double f(double,double), int N_x, int N_y, int x0, int y0, double h_x, double h_y)
	{
		int i = 0;
		int j = 0;
		int I = 0;
		int J = 0;
		x0=0;
		y0=0;
		/*for(i = 0; i < N_x; i++)
		{
			mat[i*n] = 0;
		}
		for(i = 0; i < N_x; i++)
		{
			mat[i*n + N_y-1] = 0;
		}
		for(j = 0; j < N_y; j++)
		{
			mat[j] = 0;
		}
		for(j = 0; j < N_y; j++)
		{
			mat[(N_x-1)*n+j] = 0;
		}*/
		for(i = 0; i < N_x; i++)
		{
			for (j = 0; j < N_y ; j++)
			{
				I = double(i);
				J = double(j);
				mat[i*n+j] = f(x0+I*h_x,y0+ J*h_y);
			}
		}
	}
	void fil_mat_by_exp(const matrix& x)
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				mat[i * m + j] = exp(x.mat[i]);
	}
	void fil_mat_by_const(int n, double c)
	{
		for(int i=0;i<n;i++)
		{
			mat[i] = c;
		}
	}
	void fil_mat_by_first(int N, const matrix& x)
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < m; j++)
				mat[i * m + j] = x.mat[i];
	}
	void fil_mat_by_seq()
	{
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				mat[i * m + j] = i;
	}
	void friend up_triangle(matrix& A, matrix& b)
	{
		double d;
		int n = A.Dim_rows();
		int m = A.Dim_cols();
		for (int k = 0; k < n; k++) // ������ ���
		{
			if (abs(A.mat[k * m + k]) > 1e-14)
			{
				for (int i = k + 1; i < m; i++)
				{
					if (abs(A.mat[i * m + k]) > 1e-14)
					{
						d = A.mat[i * m + k] / A.mat[k * m + k];
						for (int j = k; j < n; j++)
						{
							A.mat[i * m + j] = A.mat[i * m + j] - d * A.mat[k * m + j];
						}
						b.mat[i] = b.mat[i] - d * b.mat[k];
					}
				}
			}
			else
			{
				double buf = 0;
				for (int i = k + 1; i < m; i++)
				{
					if (abs(A.mat[i * m + k]) > 1e-14)
					{
						for (int j = k; j < m; j++)
						{
							buf = A.mat[k * m + j];
							A.mat[k * m + j] = A.mat[i * m + j];
							A.mat[i * m + j] = buf;
						}
					}
				}
				for (int i = k + 1; i < m+1; i++)
				{
					if (i == m + 1)
					{
						cout << "System can not be solved";
						break;
					}
					if (abs(A.mat[i * m + k]) > 1e-14)
					{
						d = A.mat[i * m + k] / A.mat[k * m + k];
						for (int j = k; j < n; j++)
						{
							A.mat[i * m + j] = A.mat[i * m + j] - d * A.mat[k * m + j];
						}
						b.mat[i] = b.mat[i] - d * b.mat[k];
					}
				}
			}
		}
	}
	int friend mat_mult(matrix& A, matrix& B)
	{
		int n_A = A.Dim_rows();
		int m_A = A.Dim_cols();
		int n_B = B.Dim_rows();
		int m_B = B.Dim_cols();
		if(m_A != m_B)
		{
			return -1;
		}
		for(int i = 0; i < n_A; i++)
		{
			for(int j = 0; j < n_B; j++)
			{
				for(int k = 0; k < m_B; k++)
				{
				}
			}
		}
	}
	//только для векторов работает
	void friend smart_mult_mat(matrix &A, matrix &B, matrix &C, matrix& X, matrix& Y)
	{
		int N_plus_1 = X.Dim_rows()-1;
		for(int i = 1; i < N_plus_1; i++)
		{
			Y.SetMij(i,0, X.mat[i-1]*A.mat[i] + X.mat[i]*B.mat[i]+X.mat[i+1] * C.mat[i]);
		}
	}
	void friend mat_mult_scalar(double L, matrix& B) //need add to lib
	{
		int n_B = B.Dim_rows();
		int m_B = B.Dim_cols();
		//cout << n_B << ' ' << m_B << endl;
		//L = L*2;
		for(int i = 0; i < n_B; i++)
		{
			for(int j = 0; j < m_B; j++)
			{
				B.mat[i] = L*B.mat[i] ;
			}
		}
	}
	void friend gauss_back(matrix& A, matrix& b)
	{
		double d;
		int n = A.Dim_rows();
		int m = A.Dim_cols();
		for (int i = 0; i < n; i++)
		{
			d = A.mat[i * m + i];
			for (int j = i; j < n; j++)
			{
				A.mat[i * m + j] /= d;
			}
			b.mat[i] /= d;
		}
		for (int j = n - 1; j >= 0; j--)
		{
			for (int i = 0; i < j; i++)
			{
				b.mat[i] = b.mat[i] - b.mat[j] * A.mat[i * m + j];
				A.mat[i * m + j] = 0;
			}
		}
		//x = b;
	}
	~matrix()
	{
		delete[] mat;
	}

};
double stand_pol_no_first(double x, matrix& coef);
	void filmat_diag(int n, double spec_h, double h);
double solution_mnrp(int n, int N, matrix& x, matrix& y, matrix& p, matrix& f, matrix& ind, matrix& vec);
int change_points(double n, double N, double index, matrix& x, matrix& y, matrix& p, matrix& f,matrix& vec, matrix& ind);
double H_at_point(double index, matrix& p, matrix& f, matrix& vec);
double H_max(int N, double& ind, matrix& x, matrix& y, matrix& vec);
	void filmat_h_first(int n, int zero_ind, double spec_h, double h);
void fil_mat_by_x_and_h(const matrix& x);
void fil_mat_by_first(int N, const matrix& x);
void fil_mat_by_seq();
void cheb_points(double a, double b, double N, matrix& p);
double stand_pol(double x, matrix& coef);
double lagr_pol(double x, matrix& p, matrix& f); // �������
void uniform_points(double a, double b, double N, matrix& p);
void c_rand_points(double a, double b, double N, matrix& p);
void generation(string filename, double a, double b, double N, int p_type, double(*f)(double));
double abs_x(double x);
double runge_x(double x); void fill_mat(string filename, matrix& x, matrix& f);
void val_filler(string filename, double a, double b, double(*f)(double));
void val_filler(string filename, double a, double b, double(*f)(double, matrix&), matrix& p);
void val_filler(string filename, double a, double b, double(*f)(double, matrix&, matrix&), matrix& p, matrix& val);
double pol_1(double x);
void val_filler(int N, string filename, double a, double b, double(*f)(double));
void val_filler(int N, string filename, double a, double b, double(*f)(double, matrix&), matrix& p);
void val_filler(int N, string filename, double a, double b, double(*f)(double, matrix&, matrix&), matrix& p, matrix& val);
double step_h(int n, matrix& x, matrix& y);
	void print_mat();
	int solution_progonka(matrix& x, matrix& a, matrix& b, matrix& c, matrix& d);
	void mult_three_diag_and_vec(int n, matrix& a,matrix& b, matrix& c, matrix& x, matrix& d);
	double residual(int n, matrix& A, matrix& B);
	double fourier_2ds(matrix &u, matrix &d, matrix &c, int N_x, int N_y, double x0, double y0, double h_x, double h_y);
double sin_2d(double x,double y);
double mat_dif(matrix &a, matrix &b, int N_x, int N_y);
void out_mat_to_file_fun(matrix &u, string filename, int N_x, int N_y, double h_x, double h_y, double x0, double y0);
#endif
