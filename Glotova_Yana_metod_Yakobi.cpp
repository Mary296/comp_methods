#include <iostream>
#include <math.h>
#include <iomanip> 
#include <fstream>

using namespace std;

int main()
{
	setlocale(LC_ALL, "rus");
	double a, b; //коэф-ты
	int n;
	cout << "Введите n: ";
	cin >> n;
	a = -0.9;
	b = -1.1;
	
	double *x = new double[n + 1];
	double *y = new double[n + 1];
	double **u;
	u = new double*[n + 1];
	for(int i = 0; i < n + 1; i++)
		u[i] = new double[n + 1];
	x[0] = 0;
	x[n] = 1;
	y[0] = 0;
	y[n] = 1;

	double h = 1 / double(n); // шаг

	for (int i = 1; i <= n - 1; i++)
	{
		x[i] = x[0] + i*h;
	}
	for (int i = 1; i <= n - 1; i++)
	{
		y[i] = y[0] + i*h;
	}

	double max = -1;
	double e = 1e-20; // ������� ��������� ������ �����
	double **f;
	f = new double*[n + 1];
	for(int i = 0; i < n + 1; i++)
		f[i] = new double[n + 1];

	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			u[i][j] = 0;
		}
	}

	double **v; // ������� ��������� �������
	v = new double*[n + 1];
	for (int i = 0; i <= n; i++)
		v[i] = new double[n + 1];
	for(int i = 0; i <= n; i++)
		for(int j = 0; j <= n; j++)
			v[i][j] = 0;
	for(int i = 0; i < n / 2; i++)
	{
		for(int j = 0; j < n / 2; j++)
		{
			if(i + j >= n / 2)
				v[i][j] = 1 / (1 + x[i]*x[i] + y[j]*y[j]);
		}
		for(int j = n / 2; j < n; j++)
		{
			if(j - i <= n / 2)
				v[i][j] = 1 / (1 + x[i]*x[i] + y[j]*y[j]);
		}
	}
	for (int i = n / 2; i <= n; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			v[i][j] = 1 / (1 + x[i]*x[i] + y[j]*y[j]);
		}
	}
	//������� �� �������:
	for(int i = n / 2; i <= n; i++)
	{
		u[i][0] = v[i][0];
		u[i][n] = v[i][n];
	}
	for(int i = 0; i <= n; i++)
	{
		u[n][i] = v[n][i];
	}
	for(int i = 0; i <= n / 2; i++)
	{
		for(int j = 0; j <= n / 2; j++)
		{
			if(i + j == n / 2)
				u[i][j] = v[i][j];
		}
		for(int j = n / 2; j <= n; j++)
		{
			if(j - i == n / 2)
				u[i][j] = v[i][j];
		}
	}
	/*
	for (int i = 0; i < n + 1; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			cout << u[i][j] << "  ";
		}
		cout << endl;
	}
	cout << "\n";
	*/
	for (int i = 0; i <= n; i++){
		for (int j = 0; j <= n; j++){
			f[i][j] = (4 - 3.2*x[i]*x[i] - 4.8*y[j]*y[j]) / pow(1 + x[i]*x[i] + y[j]*y[j], 3); // ������� �� ������ �����
		}
	}

	double **u1;
	u1 = new double*[n + 1];
	for(int i = 0; i < n + 1; i++)
		u1[i] = new double[n + 1];
	int q = 0; // ���������� ��������

	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			u1[i][j] = 0;
		}
	}
	for(int i = n / 2; i <= n; i++)
	{
		u1[i][0] = v[i][0];
		u1[i][n] = v[i][n];
	}
	for(int i = 0; i <= n; i++)
	{
		u1[n][i] = v[n][i];
	}
	for(int i = 0; i <= n / 2; i++)
	{
		for(int j = 0; j <= n / 2; j++)
		{
			if(i + j == n / 2)
				u1[i][j] = v[i][j];
		}
		for(int j = n / 2; j <= n; j++)
		{
			if(j - i == n / 2)
				u1[i][j] = v[i][j];
		}
	}
	//----------����----------------�
	do
	{
		for(int i = 0; i < n / 2; i++)
		{
			for(int j = 0; j < n / 2; j++)
			{
				if(i + j > n / 2)
					u1[i][j] = (a*(u[i - 1][j] + u[i + 1][j]) + b*(u[i][j - 1] + u[i][j + 1]) - f[i][j] * h*h) / (2 * a + 2 * b);
			}
			for(int j = n / 2; j < n; j++)
			{
				if(j - i < n / 2)
					u1[i][j] = (a*(u[i - 1][j] + u[i + 1][j]) + b*(u[i][j - 1] + u[i][j + 1]) - f[i][j] * h*h) / (2 * a + 2 * b);
			}
		}
		for (int i = n / 2; i <= n - 1; i++)
		{
			for (int j = 1; j <= n - 1; j++)
			{
				u1[i][j] = (a*(u[i - 1][j] + u[i + 1][j]) + b*(u[i][j - 1] + u[i][j + 1]) - f[i][j] * h*h) / (2 * a + 2 * b); // ���������� ������������� �������
			}
		}
	
		q++; // �������
		// ��� ������� ���������:
		max = -1;
		for (int i = 0; i <= n; i++){
			for (int j = 0; j <= n; j++){
				if ((fabs(u1[i][j] - u[i][j])) > max)
				{
					max = fabs(u1[i][j] - u[i][j]);
				}
			}
		}
		
		double p;
		p = -1;
		for (int i = 0; i <= n; i++){
			for (int j = 0; j <= n; j++)
			{
				if ((fabs(u[i][j] - v[i][j])) > p)
				{
					p = fabs(u[i][j] - v[i][j]);
				}
			}
		}
		/*
		// ����� ������������� ����, ����� �������� � �����������:
		cout << "������������ ���: " << q << endl;
		cout << "����� �������� �� ���� ����: " << max << endl;
		cout << "����������� (��������� ���������� � ��������� �������) �� ���� ����: " << setprecision(10) << p << endl;
		cout << "\n";
		*/
		for (int i = 0; i <= n; i++)
		{
			for (int j = 0; j <= n; j++)
			{
				u[i][j] = u1[i][j];
			}
		}

	} while (max > e);

	cout << "���������� ��������: " << q << endl;
	cout << "��� �����: " << h << endl;
	cout << "max: " << max << endl;
	
	double p;
	p = 0;
	for (int i = 0; i <= n; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			if ((fabs(u[i][j] - v[i][j])) > p)
			{
				p = fabs(u[i][j] - v[i][j]);
			}
		}
	}
	cout << "����������� (��������� ���������� � ��������� �������): " << p << endl;
	return 0;
}
