// integral.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>
#include "mpi.h"

using namespace std;


int My_Scatter(void*sendbuf, int sendcount, MPI_Datatype sendtype, void*recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
	int myrank, mysize;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &mysize);
	int tag = 99;
	if (myrank == root) { // корень рассылает
		int lb, stype_ext;
		MPI_Type_get_extent(sendtype, &lb, &stype_ext);
		for (int i = 0; i < mysize; i++) {
			if (i != root) { // отправляем всем, кроме себя
				MPI_Send((char*)sendbuf + i * sendcount*stype_ext, sendcount, sendtype, i, tag, comm);
			}
			else { // себе копируем
				for (int j = 0; j < sendcount*stype_ext; j++) {
					*((char*)recvbuf + j) = *((char*)sendbuf + i * sendcount*stype_ext + j);
				}
			}
		}
	}
	else { // некорневые принимают
		MPI_Status status;
		MPI_Recv(recvbuf, recvcount, recvtype, root, tag, comm, &status);
	}
	return 0;
}


int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	int myrank, mysize;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &mysize);

	int *sbuf = 0, *rbuf, *rbuf2;
	int n = 10; 

	// создаем и инициализируем отправной буфер (на корневом процессоре)
	if (myrank == 0) {
		sbuf = new int[n*mysize];
		for (int i = 0; i < n*mysize; i++) {
			sbuf[i] = i;
		}
	}

	// выделяем память под принимающий буфер (на всех процессорах)
	rbuf = new int[n];
	rbuf2 = new int[n];

	// так работает оригинальный MPI_Scatter
	MPI_Scatter(sbuf, n, MPI_INTEGER, rbuf, n, MPI_INTEGER, 0, MPI_COMM_WORLD);
	cout << " # proc(" << myrank << "/" << mysize << ") MPI_Scatter: ";
	for (int i = 0; i < n; i++) {
		cout << rbuf[i] << " ";
	}
	cout << endl;

	// а так - собственная реализация
	My_Scatter(sbuf, n, MPI_INTEGER, rbuf2, n, MPI_INTEGER, 0, MPI_COMM_WORLD);
	cout << " # proc(" << myrank << "/" << mysize << ") My_Scatter: ";
	for (int i = 0; i < n; i++) {
		cout << rbuf2[i] << " ";
	}
	cout << endl;

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
