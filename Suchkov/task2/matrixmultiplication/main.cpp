#include <mpi.h>
#include <iostream>
#include <ctime>
#include <random>
#include <cstdlib>
#include <time.h>
#include <Windows.h>
#include <queue>
#include <iostream>

using namespace std;
int main(int argc, char** argv)
{
	srand(time(0));
	int ProcNum, ProcRank;
	MPI_Init(&argc, &argv);
	MPI_Status status;

	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	//Начало
		//Чтение данных		
		//Матрицы M x N (M - строк, N - столбцов)
	int L = 3;
	int N = 3;
	int M = 3;

	int min = 1;
	int max = 9;

	if (argc > 1)
		L = atoi(argv[1]);
	if (argc > 2)
		M = atoi(argv[2]);
	if (argc > 3)
		N = atoi(argv[3]);

	if (argc > 4)
		min = atoi(argv[4]);
	if (argc > 5)
		max = atoi(argv[5]);

	int sizeA = L * M;
	int sizeB = M * N;
	int sizeC = L * N;
	float* ArrB = new float[sizeB];

	MPI_Group UsingProcces;
	MPI_Comm_group(MPI_COMM_WORLD, &UsingProcces);

	int usingProcces;
	if (ProcNum > sizeC)
		usingProcces = sizeC;
	else usingProcces = ProcNum;
	int* ArrUsingProcces = new int[usingProcces];
	for (int i = 0; i < usingProcces; i++)
	{
		ArrUsingProcces[i] = i;
	}
	MPI_Group_incl(UsingProcces, usingProcces, ArrUsingProcces, &UsingProcces);
	delete[] ArrUsingProcces;


	MPI_Comm Using;
	MPI_Comm_create(MPI_COMM_WORLD, UsingProcces, &Using);
	MPI_Group_free(&UsingProcces);

	ProcNum = usingProcces;

	//нагругка на каждый поток
	int* loadProcces = new int[ProcNum];
	for (int rank = 0; rank < ProcNum; rank++)
	{
		loadProcces[rank] = sizeC / ProcNum;
	}
	for (int i = 0; i < sizeC % ProcNum; i++)
	{
		loadProcces[i] ++;
	}


	int* displs = new int[ProcNum];
	displs[0] = 0;
	for (int i = 1; i < ProcNum; i++)
	{
		displs[i] = loadProcces[i - 1] + displs[i - 1];
	}


	int* sizeScat = new int[ProcNum];
	for (int i = 0; i < ProcNum; i++)
	{
		int sizeL = 1;
		for (int j = 0; j < loadProcces[i]; j++)
		{
			if ((((displs[i] + j) % N) == (N - 1)) && (j + 1 < loadProcces[i]))
				sizeL++;
		}
		sizeScat[i] = sizeL * M;
	}

	int* displsScat = new int[ProcNum];
	displsScat[0] = 0;
	for (int i = 1; i < ProcNum; i++)
	{
		displsScat[i] = displs[i] / N * M;
	}

	float* ArrBuff = new float[sizeScat[ProcRank]];
	float* ArrResult = new float[loadProcces[ProcRank]];

	int* ArrIndexB = new int[loadProcces[ProcRank]];
	for (int i = 0; i < loadProcces[ProcRank]; i++)
	{
		ArrIndexB[i] = (displs[ProcRank] + i) % N;
	}


	//НУЛЕВОЙ
	if (ProcRank == 0)
	{
		//Инициализация

		float* ArrA = new float[sizeA];
		float* ArrC = new float[sizeC];
		float* ArrCC = new float[sizeC];

		for (int i = 0; i < sizeA; i++)
		{
			ArrA[i] = rand() % (max - min + 1) + min;
		}
		for (int i = 0; i < sizeB; i++)
		{
			ArrB[i] = rand() % (max - min + 1) + min;
		}
		cout << "size matrixA " << sizeA << " (" << L << "x" << M << ") " << endl;
		cout << "matrix A" << endl;
		for (int i = 0; i < L; i++)
		{
			for (int j = 0; j < M; j++)
			{
				cout << ArrA[i * M + j] << " ";
			}
			cout << endl;
		}
		cout << endl;
		cout << "size matrixB " << sizeB << " (" << M << "x" << N << ") " << endl;
		cout << "matrix B" << endl;
		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < N; j++)
			{
				cout << ArrB[i * N + j] << " ";
			}
			cout << endl;
		}

		//рассылка матрицы B
		MPI_Bcast(ArrB, sizeB, MPI_FLOAT, 0, Using);

		//расчёт нагрузки на процессы
		// максимум sizeC операций и рассылкой N раз каждой из L строк матрицы A 


		cout << endl;

		cout << "rankC" << endl;
		int rank = 0;
		int counter = 0;

		for (int j = 0; j < ProcNum; j++)
		{
			for (int i = 0; i < loadProcces[rank]; i++)
			{
				cout << rank << " ";
				counter++;

				if ((counter) % N == 0)
					cout << endl;
			}
			rank++;
		}
		cout << endl;

		//рассылка строк
		MPI_Scatterv(ArrA, sizeScat, displsScat, MPI_FLOAT, ArrBuff, sizeScat[ProcRank], MPI_FLOAT, 0, Using);


		int sizeStr = 0;
		for (int i = 0; i < loadProcces[ProcRank]; i++)
		{
			ArrResult[i] = 0;
			for (int j = 0; j < M; j++)
			{
				ArrResult[i] += ArrBuff[j + M * sizeStr] * ArrB[j * N + ArrIndexB[i]];    //Cij = summ(Aim*Bmj);
				//cout << " ( " << j + M * sizeStr << " : " << j * N + ArrIndexB[i] << " ) ";
			}
			if ((((displs[ProcRank] + i) % N) == (N - 1)) && (i + 1 < loadProcces[ProcRank]))
				sizeStr++;
		}
		//cout << endl;

		// сборка
		MPI_Gatherv(ArrResult, loadProcces[ProcRank], MPI_FLOAT, ArrC, loadProcces, displs, MPI_FLOAT, 0, Using);

		//печать С


		cout << "size matrixC " << sizeC << " (" << L << "x" << N << ") " << endl;
		cout << "matrix C" << endl;
		for (int i = 0; i < L; i++)
		{
			for (int j = 0; j < N; j++)
			{
				cout << ArrC[i * N + j] << " ";
			}
			cout << endl;
		}

		for (int i = 0; i < L; i++)
		{
			for (int j = 0; j < N; j++)
			{
				ArrCC[i * N + j] = 0;
				for (int k = 0; k < M; k++)
				{
					ArrCC[i * N + j] += ArrA[i * M + k] * ArrB[k * N + j];
				}
			}
		}

		cout << "matrix CC" << endl;
		for (int i = 0; i < L; i++)
		{
			for (int j = 0; j < N; j++)
			{
				cout << ArrCC[i * N + j] << " ";
			}
			cout << endl;
		}

		// clear
		delete[] ArrA;
		delete[] ArrC;

		delete[] ArrCC;
	}
	//работа других процессов
	else
	{
		// получение матрицы B
		MPI_Bcast(ArrB, sizeB, MPI_FLOAT, 0, Using);


		//принятие строк
		MPI_Scatterv(nullptr, sizeScat, displsScat, MPI_FLOAT, ArrBuff, sizeScat[ProcRank], MPI_FLOAT, 0, Using);


		//cout << "Procces " << ProcRank <<" "<< sizeScat[ProcRank] << endl;

		/*cout << "J " << endl;
		for (int i = 0; i < loadProcces[ProcRank]; i++)
		{
			cout << " " << ArrIndexB[i];
		}
		cout << endl;
*/

		int sizeStr = 0;
		for (int i = 0; i < loadProcces[ProcRank]; i++)
		{
			ArrResult[i] = 0;
			for (int j = 0; j < M; j++)
			{
				ArrResult[i] += ArrBuff[j + M * sizeStr] * ArrB[j * N + ArrIndexB[i]];    //Cij = summ(Aim*Bmj);
			}
			if ((((displs[ProcRank] + i) % N) == (N - 1)) && (i + 1 < loadProcces[ProcRank]))
				sizeStr++;
		}

		MPI_Gatherv(ArrResult, loadProcces[ProcRank], MPI_FLOAT, nullptr, loadProcces, displs, MPI_FLOAT, 0, Using);

	}

	//общая часть
	delete[] sizeScat;
	delete[] displsScat;

	delete[] ArrIndexB;
	delete[] ArrResult;

	delete[] ArrB;
	delete[] ArrBuff;

	delete[]loadProcces;
	delete[]displs;


	MPI_Comm_free(&Using);
	MPI_Finalize();

}