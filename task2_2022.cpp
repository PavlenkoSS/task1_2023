// task2_2022.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include "matrlib.h"


using namespace std;

int main()
{	
	// Генерация точек
	int m = 5;
	double a = 0;
	double b = 1;
	double (*f)(double);
	/* ТОЧКИ: 1 - равномерно 2 - чебышев 3 - случайные
	* ФУНКЦИИ: exp, runge_x, abs_x
	*/
	f = exp;
	string filename = "data_points.txt";
	generation(filename, a, b, m, 1, f);

	// заполнение массивов точками и значениями
	int n = 1;
	ifstream fin;
	fin.open(filename);
	fin >> n;
	fin.close();
	matrix x_points(n, 1);
	matrix f_points(n, 1);
	fill_mat(filename, x_points, f_points);
	x_points.out_mat();
	f_points.out_mat();

	// решение СЛУ для коэф-ов 

	matrix A(n);
	A.fil_mat_by_x(x_points);
	matrix vec;
	vec = f_points;
	up_triangle(A,vec);
	gauss_back(A, vec);

	// заполнение файлов

	val_filler("dat1.txt", a, b, exp);
	val_filler("dat2.txt", a, b, stand_pol, vec);
	val_filler("dat3.txt", a, b, lagr_pol, x_points, f_points);

	return 0;
}
