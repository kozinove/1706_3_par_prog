// integral.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>
#include <math.h>
#include "mpi.h"

using namespace std;

// evaluate integral of this function :
double fun_1(double x) {
	return sin(x);
}

double fun_2(double x) {
	return pow(x, 3);
}

double fun_3(double x) {
	return exp(x);
}


double calc_integral(double a, double b, int n, int id) {
	double dx = (b - a) / n;
	double sum = 0;
	double f;
	for (int i = 0; i < n; i++) {
		double x = a + i * dx + dx * .5; // central
		switch (id) {
			case 1: f = fun_1(x); break;
			case 2: f = fun_2(x); break;
			case 3: f = fun_3(x); break;
			default:;
		}
		sum += f * dx;
	}
	return sum;
}

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	int myrank, mysize;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &mysize);
	cout << " # rank: " << myrank << " size: " << mysize << endl;

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Status status;

	double a0, b0;
	int id, n;
	if (myrank == 0) { // master process

		cout << " ------------------------- " << endl ;

		// select function to integrate
		do {
			cout << "select function to calc integral: " << endl;
			cout << " 1: sin(x)" << endl;
			cout << " 2: x^3" << endl;
			cout << " 3: e^x" << endl;
			cin >> id;
		} while (id <= 0 && id > 3);

		double a, b;

		do {
			cout << "enter integrity bounds [a, b] :" << endl;
			cin >> a;
			cin >> b;
		} while (a >= b);

		do {
			cout << "enter integrity precision(1000 ~100000) :" << endl;
			cin >> n;
		} while (n <= 0);

		double dx = (b - a) / mysize;

		for (int i = 1; i < mysize; i++) {
			a0 = a + dx * i;
			b0 = a0 + dx;
			// send data to all slave processes (2-points exchange)
			MPI_Send(&a0, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD);
			MPI_Send(&b0, 1, MPI_DOUBLE, i, 99, MPI_COMM_WORLD);
			MPI_Send(&n, 1, MPI_INTEGER, i, 99, MPI_COMM_WORLD);
			MPI_Send(&id, 1, MPI_INTEGER, i, 99, MPI_COMM_WORLD);

		}

		// calc bounds on master process
		a0 = a;
		b0 = a0 + dx;

	}

	// slave processes
	else {
		// get data on slave processes

		MPI_Recv(&a0, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
		MPI_Recv(&b0, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
		MPI_Recv(&n, 1, MPI_INTEGER, 0, 99, MPI_COMM_WORLD, &status);
		MPI_Recv(&id, 1, MPI_INTEGER, 0, 99, MPI_COMM_WORLD, &status);

	}

	// all processes

	cout << "#" << myrank << ": [a0, bO]: (" << a0 << ", " << b0 << ") ; n = " << n << endl;

	double sum = calc_integral(a0, b0, n, id);

	MPI_Barrier(MPI_COMM_WORLD);

	cout << "#" << myrank << " local sum = " << sum << endl;

	// get sum on master process
	double all_sum;
	MPI_Reduce(&sum, &all_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	if (myrank == 0) { // master process
		cout << " ------------------------- "<< endl << endl;
		cout << " Integral = " << all_sum << endl;

		// check

		double a = a0;

		double b = a + (b0 - a0) * mysize;
		int n2 = n * mysize;

		double check_sum = calc_integral(a, b, n2, id);
		cout << " check = " << check_sum << endl;
	}

	MPI_Finalize();
	return 0;
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
