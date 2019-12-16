#include "mpi.h"
#include <cstdlib>
//#include <time.h>
#include <Windows.h>
#include <queue>
#include <iostream>
#include <vector>


#define INFINITI 5000
using namespace std;

void printD(int* d, int size)
{
	for (int i = 0; i < size; i++)
	{
		cout << d[i] << " ";
	}
	cout << endl;
}
int* initG(int countEdge, int countVertex)
{
	cout << "Matrix" << endl;
	int rank = 0, procNum = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int* oneG = new int[countVertex * countVertex];
	int** G = new int* [countVertex];
	int ver = 25, veroyatn;
	for (int i = 0; i < countVertex; i++)
	{
		G[i] = new int[countVertex];
	}
	//srand(time(NULL));
	int t = 0;
	for (int i = 0; i < countVertex; i++) {
		cout << i << ": ";
		for (int j = 0; j < countVertex; j++) {
			if (i == j)
			{
				G[i][j] = 0;

			}
			else
			{
				G[i][j] = rand() % 5;

				veroyatn = rand() % 100 + 1;
				//if (veroyatn <= ver)
				if (G[i][j] == 0)
				{
					G[i][j] = INFINITI;
				}
			}
			cout << G[i][j] << " ";
		}
		cout << endl;
	}
	for (int i = 0, t = 0; i < countVertex; i++)
	{
		for (int j = 0; j < countVertex; j++)
		{
			oneG[t] = G[i][j];
			t++;
		}
	}
	//printD(oneG, countVertex * countVertex);
	//cout << endl;
	return oneG;
}


int* initD(int size)
{
	int* d = new int[size];
	for (int i = 0; i < size; i++) {
		d[i] = INFINITI;
	}
	return d;
}

int* stepByStep(int* G, int start, int countVertex) {
	double t1 = 0.0, t2 = 0.0, tTotal = 0.0;
	t1 = MPI_Wtime();
	
	int* d = initD(countVertex);//все infinity кроме стартовой
	d[start] = 0;
	queue<pair<int, int>>  q;
	q.push(make_pair(0, start));
	t2 = MPI_Wtime();
	tTotal += t2 - t1;
		
	while (!q.empty())
	{
		t1 = MPI_Wtime();
		int v = q.front().second, cur_d = q.front().first;
		q.pop();
		if (cur_d > d[v])  continue;
		t2 = MPI_Wtime();
		tTotal += t2 - t1;
		for (int j = 0; j < countVertex; ++j)
		{
			int to = j,
				len = G[v * countVertex + j];//находим вес
			if (d[v] + len < d[to])
			{
				d[to] = d[v] + len;
				q.push(make_pair(d[to], to));
			}
		}
	}
	//cout << "time sending step by step = " << tTotal << endl;
	return d;
}
int* parallel(int* G, int start, int countVertex)
{
	int* d = nullptr;
	int rank = 0, procNum = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int partSize = countVertex / procNum;
	   	
	int* Graph = new int[countVertex * countVertex];
	if (rank == 0)
	{
		for (int i = 0; i < countVertex * countVertex; i++)
			Graph[i] = G[i];
	}
	double vremya1 = MPI_Wtime();
	MPI_Bcast(Graph, countVertex * countVertex, MPI_INT, 0, MPI_COMM_WORLD);
	double vremya2 = MPI_Wtime();
	d = initD(countVertex);
	d[start] = 0;
	queue<pair<int, int>>  q;
	q.push(make_pair(0, start));

	while (!q.empty())
	{
		int v = q.front().second,//вершина
			dest = q.front().first;//расстояние
		q.pop();

		for (int j = 0; j < countVertex; j++)//диапазон вершины
		{
			int to = j;
			int	len = Graph[v * countVertex + j];

			if (dest + len < d[to])
			{
				
				d[to] = dest + len;
								

				if (to >= rank * partSize && to < (rank + 1) * partSize)//если to входит в нужный диапазон
				{
					q.push(make_pair(d[to], to));
					//cout << rank << "   " << to << endl;
				}

				else if (to >= rank * partSize && to < countVertex)//если to входит в нужный диапазон//если не взаимно простые, то зайдет в диапазон
				{
					q.push(make_pair(d[to], to));
					//cout << rank << "   " << to << endl;
				}

			}
		}

	}
	   	
	int* final_d = new int[countVertex];
	MPI_Reduce(d, final_d, countVertex, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
		//if (rank == 0)
			//printD(final_d, countVertex);

	if (rank == 0)
	{
		return final_d;
	}
}
void isCorrectImplementation(int* d1, int* d2, int size) {
	int countMistakes = 0;
	for (int i = 0; i < size; i++) {
		if (d1[i] != d2[i]) {
			cout << "something isn`t correct" << endl;
			countMistakes++;
			break;
		}
	}
	if (countMistakes == 0)
		cout << "well done, all is correct" << endl;
}
int main(int argc, char* argv[])
{
	int procNum = 0, rank = 0, countVertex = 0, countEdge = 0;
	MPI_Init(&argc, &argv);
	countVertex = atoi(argv[1]);////////////
	//countVertex = 8;
	//countVertex = 5;
	//srand(time(NULL));
	countEdge = (countVertex - 1) + rand() % ((countVertex * (countVertex - 1)) / 2);//слева минимум, справа
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int* dStep = nullptr;
	int* dParallel = nullptr;
	int* G = nullptr;
	double t1Step = 0.0, t2Step = 0.0, t1Parallel = 0.0, t2Parallel = 0.0;
	int start = rand() % (countVertex - 1);
	if (rank == 0)
		G = initG(countEdge, countVertex);
	//start = 0;
	//G = Route(countEdge, countVertex);
	if (procNum == 1) {
		//int* GStep = initG(countEdge, countVertex);
		t1Step = MPI_Wtime();
		//dStep = stepByStep(GStep, start, countVertex);
		dStep = stepByStep(G, start, countVertex);
		t2Step = MPI_Wtime();
		cout << "step by step implementation: " << endl;
		for (int i = 0; i < countVertex; i++)
			cout << dStep[i] << " ";
		cout << endl;
		cout << "count vertexs = " << countVertex << endl;
		//cout << "count edges = " << countEdge << endl;
		cout << "time of job = " << t2Step - t1Step << endl;
	}
	else
	{
		double par1, par2;
		if (!rank)//если не нулевой процесс
		{
			//cout << "start = " << start << endl;
			t1Parallel = MPI_Wtime();
		}
		//par1 = MPI_Wtime();
		dParallel = parallel(G, start, countVertex);
		//par2 = MPI_Wtime();
		if (rank == 0)
		{
			t2Parallel = MPI_Wtime();
			//cout << "parallel implementation: " << endl;
			cout << "count vertexs = " << countVertex << endl;
			cout << "start = " << start << endl;
			cout << "time of job Parallel = " << t2Parallel - t1Parallel << endl;
			//int** GStep = initGStep(countEdge, countVertex);
			t1Step = MPI_Wtime();
			dStep = stepByStep(G, start, countVertex);
			t2Step = MPI_Wtime();
			//cout << "step by step implementation: " << endl;
			cout << "time of job Lin = " << t2Step - t1Step << endl;
			isCorrectImplementation(dStep, dParallel, countVertex);
			cout << "Dest S: ";
			for (int i = 0; i < countVertex; i++)
				cout << dStep[i] << " ";
			cout << endl;
			
		}
	}
	MPI_Finalize();
}