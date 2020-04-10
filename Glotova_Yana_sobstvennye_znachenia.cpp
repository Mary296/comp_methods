/* студент: Глотова Яна Сергеевна
 группа: 16111
 задание 3 (найти минимальное и максимальное собственные значения А степенным методом)
*/ 
#include <iostream> 
#include <math.h> 
#include <stdio.h> 
#include <iomanip> 
#include <fstream> 
#include <conio.h> 

using namespace std;

void Operator_Act(int N, double** U, double** U1, double a, double b, double h) // действие матрицы А
{
	for (int i = 1; i < N; i++)
	{
		for (int j = 1; j < N; j++)
		{
			U1[i][j] = -a * (U[i - 1][j] - 2 * U[i][j] + U[i + 1][j]) - b * (U[i][j - 1] - 2 * U[i][j] + U[i][j + 1]);
			//cout << U1[i][j] << "  ";
		}
		//cout << endl;
	}
}

void Operator_Bact(int N, double** U2, double** U3, double a, double b, double h, double lambda_plus) // действие вспомогательной матрицы, необходимой для нахождения минимального собственного значения
{
	for (int i = 1; i < N; i++)
	{
		for (int j = 1; j < N; j++)
		{
	        U3[i][j] = lambda_plus * U2[i][j] * (h*h) + a * (U2[i - 1][j] - 2 * U2[i][j] + U2[i + 1][j]) + b * (U2[i][j - 1] - 2 * U2[i][j] + U2[i][j + 1]);
			//cout << U1[i][j] << "  ";
		}
		//cout << endl;
	}
}

double Norma(int N, double** U) // норма l2
{
	double norm = 0;
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < N + 1; j++)
		{
			norm += U[i][j] * U[i][j];
		}
	}
	return sqrt(norm);
}

double Scalar_Prod(int N, double** U, double** U1) // скалярное произведение
{
	double sc_pr = 0;
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < N + 1; j++)
		{
			sc_pr += U[i][j] * U1[i][j];
		}
	}
	return sc_pr;
}

int main()
{
	setlocale(LC_ALL, "rus");
	
	double a, b, h;
	a = 0.9;   b = 1.1;
	
	int N;
	cout << "Введите N: ";
	cin >> N;
	h = 1. / (double)N;
	cout << "шаг = " << h << endl;

/////////////////////////////////////////////////////////////////////////////////	
	double** U = new double*[N + 1];
	for (int i = 0; i <= N; i++)
	{
		U[i] = new double[N + 1];
	}
	
	// массив (начальный вектор из метода) строится по области:
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < N + 1; j++)
		{
			U[i][j] = 0;
		}
	}
	for (int i = 1; i < N / 2; i++)
	{
		for (int j = 1 + N / 2 - i; j < N / 2 + i; j++)
		{
			U[i][j] = 1;
		}
	}
	for (int i = N / 2; i < N; i++)
	{
		for (int j = 1; j < N; j++)
		{
			U[i][j] = 1;
		}
	}
	/*
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < N + 1; j++)
		{
			cout << U[i][j] << "  ";
		}
		cout << endl;
	}
	cout << "\n";
	*/
	double** U1 = new double*[N + 1];
	for (int i = 0; i <= N; i++)
	{
		U1[i] = new double[N + 1];
	}
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < N + 1; j++)
		{
			U1[i][j] = 0;
			//cout << U1[i][j] << "  ";
		}
		//cout << endl;
	}
/////////////////////////////////////////////////////////////////////////////////
	double lambda, lambda2, eps, e1;
	lambda2 = 0;
	eps = 1e-12;
	cout << "остановка при изменении на итерации <= " << setprecision(10) <<eps << endl << endl;
	
	int iter1=0;
	
	// цикл для нахождения максимального с.зн.:
	do
	{
		iter1++;
		lambda = lambda2;
		double nor = Norma(N, U);
		for (int i = 0; i < N + 1; i++) // нормировка массива
		{
			for (int j = 0; j < N + 1; j++)
			{
				U[i][j] = U[i][j] / nor;
			}
		}

		Operator_Act(N, U, U1, a, b, h);
		lambda2 = Scalar_Prod(N, U, U1)/(h*h);
		
		for (int i = 1; i < N; i++)
		{
			for (int j = 1; j < N; j++)
			{
				U[i][j] = U1[i][j]; // / Norma(N, U);
			}
		}
		e1 = fabs((lambda - lambda2) / lambda2);
		
//		cout << "итерация: " << iter1 << "       макс. с.з = " << setprecision(20) << lambda2 << endl;
	} while (e1 > eps);

	// вывод макс. собственного значения:
	cout << "максимальное с.з.: " << setprecision(10) << lambda2 << endl;
	cout << "количество итераций n: " << iter1 << endl;
	cout << "изменение на итерации |l(n)-l(n-1)|: " << setprecision(10) << e1 << endl << endl;


/////////////////////////////////////////////////////////////////////////////////
	// вычисление минимального с.з.:
	double lam, lam_2, e2;
	lam_2 = 0;
	int iter2=0;
	
	double** U2 = new double*[N + 1];
	for (int i = 0; i <= N; i++)
	{
		U2[i] = new double[N + 1];
	}
	
	double** U3 = new double*[N + 1];
	for (int i = 0; i <= N; i++)
	{
		U3[i] = new double[N + 1];
	}
	
	// массив (начальный вектор из метода) строится по области:
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < N + 1; j++)
		{
			U2[i][j] = 0;
		}
	}
	for (int i = 1; i < N / 2; i++)
	{
		for (int j = 1 + N / 2 - i; j < N / 2 + i; j++)
		{
			U2[i][j] = 1;
		}
	}
	for (int i = N / 2; i < N; i++)
	{
		for (int j = 1; j < N; j++)
		{
			U2[i][j] = 1;
		}
	}
	/*
	for (int i = 0; i < N + 1; i++)
	{
		for (int j = 0; j < N + 1; j++)
		{
			cout << U2[i][j] << "  ";
		}
		cout << endl;
	}
	cout << "\n";
	*/
	for (int i = 0; i < N + 1; i++)
	{
		U3[i][0] = 0;
		U3[i][N] = 0;
		U3[0][i] = 0;
		U3[N][i] = 0;
	}
	// цикл для нахождения минимального с.зн.:
	double lambda_plus=lambda2+1;
	do
	{
		iter2++;
		lam = lam_2;

		double nor1 = Norma(N, U2);
		for (int i = 0; i < N + 1; i++)
		{
			for (int j = 0; j < N + 1; j++)
			{
				U2[i][j] = U2[i][j] / nor1;
			}
		}

		Operator_Bact(N, U2, U3, a, b, h, lambda_plus);
		lam_2 = Scalar_Prod(N, U2, U3) / (h*h);

		for (int i = 1; i < N; i++)
		{
			for (int j = 1; j < N; j++)
			{
				U2[i][j] = U3[i][j];// / Norma(N, U);
				//cout << U2[i][j] << "  ";
			}
			//cout << endl;
		}
		e2 = fabs((lam - lam_2) / lam_2);
		
//		cout << lam <<  "  " << lam_2 << endl;
//		cout << "итерация: " << iter2 << "       мин. с.з = " << setprecision(20) << lambda_plus - lam_2 << endl;
	} while (e2 > eps);

	
	double z = lambda_plus - lam_2;
	// вывод мин. собственного значения:
	cout << "минимальное с.з.: " << setprecision(10) << z << endl;
	cout << "количество итераций n: " << iter2 << endl;
	cout << "изменение на итерации |l(n)-l(n-1)|: " << setprecision(10) << e2 << endl << endl;
	
	// освобождение памяти:
	delete[] U; delete[] U1; delete[] U2; delete[] U3;
	
	_getch();
	return 0;
}
