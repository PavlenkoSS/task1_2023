#include "matrlib.h"

using namespace std;
# define M_PI           3.14159265358979323846  /* pi */

double pol_1(double x)
{
	return x;
}
double abs_x(double x)
{
	return abs(x);
}

double runge_x(double x)
{
	return 1 / (25 * x * x + 1);
}
void fill_mat(string filename, matrix& x, matrix& f)
{
	ifstream fin;
	fin.open(filename);
	int i = 0,n=0;
	n = x.Dim_rows();
	fin >> n;
	double buf;
	for (i = 0; i < n; i++)
	{
		fin >> buf;
		x.SetMij(i, 0, buf);
	}
	for (i = 0; i < n; i++)
	{
		fin >> buf;
		f.SetMij(i, 0, buf);
	}
}
void uniform_points(double a, double b, double N, matrix& p)
{
	N = N - 1;
	double I = 0;
	for (int i = 0; i < N + 1; i++)
	{
		p.SetMij(i, 0, a + (b - a) * I / N);
		I++;
	}
}

void cheb_points(double a, double b, double N, matrix& p)
{
	double PI = 3.141592653589793;
	N = N + 1;
	double I = N;
	for (int i = N-1; i >= 0; i = i - 1)
	{
		p.SetMij(i, 0, (a + b) / 2 + (b - a) * cos(PI * (2 * I - 1) / (2 * N)) / 2);
		I--;
	}
}

// ����������� �����
void c_rand_points(double a, double b, double N, matrix& p)
{
	double I = 0;
	for (int i = 0; i < N; i++)
	{
		p.SetMij(i, 0, (double)(rand()) / RAND_MAX * (b - a) + a);
		I++;
	}
}	

void generation(string filename, double a, double b, double N, int p_type, double(*f)(double))
{
	ofstream fout;
	fout.open(filename);
	fout << N << endl;
	matrix p(N,1);
	switch (p_type)
	{
		case(2):
			cheb_points(a, b, N, p);
			break;
		case(3):
			c_rand_points(a, b, N, p);
			break;
		default:
			uniform_points(a, b, N, p);
			break;
	}
	for (int i = 0; i < N; i++)
	{
		fout << p.GetMij(i,0) << ' ';
	}
	fout << endl;
	for (int i = 0; i < N; i++)
	{
		fout << f(p.GetMij(i, 0)) << ' ';
	}
	fout.close();

}

double lagr_pol(double x, matrix& p, matrix& f) // �������
{
	double S = 0;
	double s = 1;
	int n = p.Dim_rows();
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i != j)
			{
				s = s * (x - p.GetMij(j, 0)) / (p.GetMij(i, 0) - p.GetMij(j, 0));
			}	
		}
		S = S+f.GetMij(i, 0) * s;
		s = 1;
	}
	return S;
}

void mult_three_diag_and_vec(int n, matrix& a,matrix& b, matrix& c, matrix& x, matrix& d)
{
	d.SetMij(0,0, b.GetMij(0,0)*x.GetMij(0,0)+c.GetMij(0,0)*x.GetMij(1,0));
	d.SetMij(n-1,0, a.GetMij(n-1,0)*x.GetMij(n-2,0)+b.GetMij(n-1,0)*x.GetMij(n-1,0));
	for(int i=1;i<n-1;i++)
	{
		d.SetMij(i,0, a.GetMij(i,0)*x.GetMij(i-1,0)+b.GetMij(i,0)*x.GetMij(i,0)+c.GetMij(i,0)*x.GetMij(i+1,0));
	}
}

double stand_pol(double x, matrix& coef) // 
{
	double s = 0;
	int n = coef.Dim_rows();
	double X = 1;
	for (int i = 0; i < n+1; i++) // n+1 -> n 
	{
		s = s + X * coef.GetMij(i, 0);
		X = X * x;
		//cout << i << ' ' << coef.GetMij(i, 1) << ' ' << X << ' ' << s <<'\n';
	}

	return s;
}
double stand_pol_no_first(double x, matrix& coef) // 
{
	double s = 0;
	int n = coef.Dim_rows();
	double X = 1;
	for (int i = 1; i < n; i++) // n+1 -> n 
	{
		s = s + X * coef.GetMij(i, 0);
		X = X * x;
		//cout << i << ' ' << coef.GetMij(i, 1) << ' ' << X << ' ' << s <<'\n';
	}

	return s;
}
void val_filler(string filename, double a, double b, double(*f)(double))
{
	ofstream fout;
	fout.open(filename);
	for (double i = a; i < b; i += 1e-3)
	{
		fout << i << ' ' << f(i) << endl;
	} 
	fout.close();
}

void val_filler(string filename, double a, double b, double(*f)(double, matrix&), matrix& p)
{
	ofstream fout;
	fout.open(filename);
	for (double i = a; i < b; i += 1e-3)
	{
		fout << i << ' ' << f(i, p) << endl;
	}
	fout.close();
}

void val_filler(string filename, double a, double b, double(*f)(double, matrix&, matrix&), matrix& p, matrix& val)
{
	ofstream fout;
	fout.open(filename);
	for (double i = a; i < b; i += 1e-3)
	{
		fout << i << ' ' << f(i, p, val) << endl;
	}
	fout.close();
}	


void val_filler(int N, string filename, double a, double b, double(*f)(double))
{
	N = N - 1;
	ofstream fout;
	fout.open(filename);
	int I=0;
	for (double i = a; i < b; i +=  a + (b - a) * I / N)
	{
		fout << i << ' ' << f(i) << endl;
		I++;
	} 
	fout.close();
}

void val_filler(int N, string filename, double a, double b, double(*f)(double, matrix&), matrix& p)
{
	N = N - 1;
	ofstream fout;
	fout.open(filename);
	int I = 0;
	for (double i = a; i < b; i += a + (b - a) * I / N)
	{
		I++;
		fout << i << ' ' << f(i, p) << endl;
	}
	fout.close();
}

void val_filler(int N, string filename, double a, double b, double(*f)(double, matrix&, matrix&), matrix& p, matrix& val)
{
	N = N - 1;
	ofstream fout;
	fout.open(filename);
	int I = 0;
	for (double i = a; i < b; i +=  a + (b - a) * I / N)
	{
		I++;
		fout << i << ' ' << f(i, p, val) << endl;
	}
	fout.close();
}	

double step_h(int n, matrix& x, matrix& y)
{
	double h = 0;
	double denum = 0;
	double mult_part = 1;
	for(int j = 0; j < n; j++)
	{
		for(int k = 0; k < n;k++)
		{
			if(k != j)
				mult_part = mult_part*(x.GetMij(j,0)-x.GetMij(k,0));
		}
		denum = denum + (1-2*(j%2))/mult_part;
		mult_part = 1;
	}
	for(int i = 0; i < n ; i++)
	{
		for(int k = 0; k < n;k++)
		{
			if(k != i)
				mult_part = mult_part*(x.GetMij(i,0)-x.GetMij(k,0));
		}
		h += (1-2*(i%2))*mult_part*y.GetMij(i,0)/denum;
	}
	return 0;
}

// h(\sigma) в конспекте
double h_max(int n, matrix& x, matrix& y)
{
	double max = 0;
	double buf;
	for(int i=0;i<n;i++)
	{
		buf = abs(y.GetMij(i, 0) - lagr_pol(x.GetMij(i, 0), x,y));
		if ((buf) > max)
			max = buf;
	}
	return max;
}
double phi_max_in_basis(int N, matrix& x, matrix& y, matrix& p, matrix& f)
{
	double max = 0;
	double buf;
	for(int i=0;i<N;i++)
	{
		buf = abs(f.GetMij(i, 0) - lagr_pol(p.GetMij(i, 0), x,y));
		if (abs(buf) > max)
			max = buf;
	}
	return max;
}

double delta_x(double x_i, double y_i, matrix& x, matrix& y)
{
	return y_i - lagr_pol(x_i, x, y);
}
// y -> f 
double H_max(int N, double& ind, matrix& x, matrix& y, matrix& vec)
{
	double max=0;
	double buf;
	for(int i=0;i<N;i++)
	{
		buf = abs(y.GetMij(i,0)-stand_pol_no_first(x.GetMij(i,0),vec));
		//cout << buf << ' '; 
		if (buf>max)
		{
			max = buf;
			ind = i;
			//cout << "\n max  = " << max << "  ind = " << ind << '\n';
		}
	}
	//cout << '\n' << '\n';
	return max;
}
double H_at_point(double index, matrix& p, matrix& f, matrix& vec)
{
	return f.GetMij(index,0)-stand_pol_no_first(p.GetMij(index,0),vec);
}

int change_points(double n, double N, double index, matrix& x, matrix& y, matrix& p, matrix& f,matrix& vec, matrix& ind)
{
	N = N+1;
	if(index<ind.GetMij(0,0))
	{
		if(H_at_point(index,p,f,vec)*H_at_point(ind.GetMij(0,0),p,f,vec)>0)
		{
			x.SetMij(0,0, p.GetMij(index,0));
			y.SetMij(0,0, f.GetMij(index,0));
			ind.SetMij(0,0, index);
		}else
		{
			for(int j=n+1;j>0;j--)
			{
				x.SetMij(j,0, x.GetMij(j-1,0));
				y.SetMij(j,0, y.GetMij(j-1,0));
				ind.SetMij(j,0, ind.GetMij(j-1,0));
			}
				x.SetMij(0,0, p.GetMij(index,0));
				y.SetMij(0,0, f.GetMij(index,0));
				ind.SetMij(0,0, index);
		}
		//cin >> N;
		return 0;
	}
	if(index > ind.GetMij(n+1,0))
	{
		if(H_at_point(index,p,f,vec)*H_at_point(ind.GetMij(n+1,0),p,f,vec)>0)
		{
			x.SetMij(n+1,0, p.GetMij(index,0));
			y.SetMij(n+1,0, f.GetMij(index,0));
			ind.SetMij(n+1,0, index);

		}else
		{
			for(int j=0;j<n+1;j++)
			{
				x.SetMij(j,0, x.GetMij(j+1,0));
				y.SetMij(j,0, y.GetMij(j+1,0));
				ind.SetMij(j, 0, ind.GetMij(j+1,0));

			}
				x.SetMij(n+1,0, p.GetMij(index,0));
				y.SetMij(n+1,0, f.GetMij(index,0));
				ind.SetMij(n+1,0, index);
		}
		return 1;
	}
	int i;
	bool ind_find = false;
	for(i = 0; i<n+2;i++)
	{
		if(((ind.GetMij(i,0)<index)&&(index < ind.GetMij(i+1,0))))
		{
			ind_find = true;
			break;
		}
	}
	if(!ind_find)
		return -1;
	if(H_at_point(index,p,f,vec)*H_at_point(ind.GetMij(i,0),p,f,vec)>0)
	{
		x.SetMij(i,0, p.GetMij(index,0));
		y.SetMij(i,0, f.GetMij(index,0));
		ind.SetMij(i,0, index);
	}else
	{
		i++;
		x.SetMij(i,0, p.GetMij(index,0));
		y.SetMij(i,0, f.GetMij(index,0));
		ind.SetMij(i,0, index);
	}

	return 2;
}

double solution_mnrp(int n, int N, matrix& x, matrix& y, matrix& p, matrix& f, matrix& ind, matrix& vec)
{
	//начальная инициализация
	int i;
	int j = 0;
	int stepi = (N-N%(n+2))/n+2-2;
	stepi = max(1,stepi);
	for( i=0; i<n+2; i++)
	{
		j+=stepi;
		x.SetMij(i, 0, p.GetMij(stepi,0));
		y.SetMij(i, 0, f.GetMij(stepi,0));
		ind.SetMij(i, 0, stepi);
	}
	matrix A(n+2);
	A.fil_mat_by_x_and_h(x);
	vec = y;
	up_triangle(A, vec);
	gauss_back(A, vec);
	double h = vec.GetMij(0,0);
	double index = 0;
	double H = H_max(N, index, p, f, vec);
	//cout << H << " 388 \n";
	int stop=0;
	if(H<h)
		return 0;
	else
	{
		while((H>=h)&&(H>1e-14))
		{
			stop = change_points(n,N,index,x,y,p,f,vec,ind);
			if(stop==-1)
				return -1;
			A.fil_mat_by_x_and_h(x);
			vec = y;
			up_triangle(A, vec);
			gauss_back(A, vec);
			h = abs(vec.GetMij(0,0));
			index = 0;
			H = H_max(N, index, p, f, vec); // итерирует index 
		}
	}
	return 0;
}
bool y_is_zero(double y)
{
	if(abs(y)<1e-15)
		return true;
	return false;
}
int solution_progonka(matrix& x, matrix& a, matrix& b, matrix& c, matrix& d)
{
	int n = a.Dim_rows();
	double y = b.GetMij(0,0);
	matrix alpha(n,1);
	matrix beta(n,1);
	if(y_is_zero(y))
		return -1;
	alpha.SetMij(0,0,-c.GetMij(0,0)/y);
	beta.SetMij(0,0,d.GetMij(0,0)/y);

	for(int i = 1;i<n-1;i++)
	{
		y = b.GetMij(i,0)+a.GetMij(i,0)*alpha.GetMij(i-1,0);
		if(y_is_zero(y))
			return -1;
		alpha.SetMij(i,0,-c.GetMij(i,0)/y);
		beta.SetMij(i,0,(d.GetMij(i,0)-a.GetMij(i,0)*beta.GetMij(i-1,0))/y);
	}
	y = b.GetMij(n-1,0)+a.GetMij(n-1,0)*alpha.GetMij(n-2,0);
	if(y_is_zero(y))
		return -1;
	beta.SetMij(n-1,0,(d.GetMij(n-1,0)-a.GetMij(n-1,0)*beta.GetMij(n-2,0))/y);
	x.SetMij(n-1,0,beta.GetMij(n-1,0));
	for(int i = n-2; i>-1;i--)
	{
		x.SetMij(i,0,alpha.GetMij(i,0)*x.GetMij(i+1,0)+beta.GetMij(i,0));
	}
	return 0;
}

double residual(int n, matrix& A, matrix& B)
{
	double s =0;
	for(int i =0;i<n;i++)
	{
		s+=abs(A.GetMij(i,0)-B.GetMij(i,0));
	}
	return s;
}
double fourier_1d_coef(matrix &f,  int N_x, double x0, double h_x, int j, double n)
{
	double sum = 0;
	x0 = 0;
	for(int i = 1; i < N_x ;i++)
	{
		sum += f.GetMij(j,i)*sin(x0+ M_PI*n*double(i)*h_x);
	}
	sum = sqrt(2)*sum/N_x;
	return sum;
}
double fourier_2ds(matrix &u, matrix &d, matrix &c, int N_x, int N_y, double x0, double y0, double h_x, double h_y)
{
	for(int j = 1; j< N_y; j++) // вниз по матрице u, вправо по d
	{
		for(int n = 1; n<N_x;n++) // вправо по матрице u, вниз по d
		{
			d.SetMij(n,j,fourier_1d_coef(u,N_x, x0, h_x,j,double(n))); // она транспонированная, чтобы потом сразу подсунуть во второй фурье
		}
	}
	for(int n =1;n<N_x;n++) // вниз по матрице c
	{
		for(int m = 1; m < N_y;m++) // вправо по матрице c 
		{
			c.SetMij(n,m, fourier_1d_coef(d,N_y,y0,h_y,n,double(m)));
		}
	}
	return 0;
}
double sin_2d(double x,double y)
{
	return 2*sin(M_PI*2*x)*sin(M_PI*3*y);
}
void out_mat_to_file_fun(matrix &u, string filename, int N_x, int N_y, double h_x, double h_y, double x0, double y0)
{
	ofstream fout; 
	fout.open(filename);
	for(int i =0; i<N_x;i++)
	{
		for(int j =0; j<N_y;j++)
			fout << x0 + double(i)*h_x << ' ' << y0+ double(j)*h_y << ' ' << u.GetMij(i,j)<< endl;
	}
	fout.close();
}
double mat_dif(matrix &a, matrix &b, int N_x, int N_y)
{
	double sum = 0; 
	double S = 0;
	for(int i = 0; i<N_x;i++)
	{
		sum = 0;
		for(int j = 0; j<N_y;j++)
		{
			sum += abs(a.GetMij(i,j)-b.GetMij(i,j));
		}
		if(sum>S)
			S = sum;
	}
	return sum;
}
/*double (matrix &c, int i, int j, double h_x,double h_y)
{
}
*/
/*void fourier_like(matrix &u,matrix &c,int N_x,int N_y,double x0,double y0,double h_x,double h_y)
{
	for(int i = 1;i<N_x-1;i++)
	{
		for(int j  = 1 ; j < N_y - 1; j++)
		{
			//u.SetMij(i,j,)
		}
	}
}
*/