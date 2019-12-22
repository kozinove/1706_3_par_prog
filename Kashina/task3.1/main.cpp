#include "mpi.h"
#include <cstdlib>
//#include <time.h>
#include <Windows.h>
#include <queue>
#include <iostream>
#include <vector>


#define INFINITI 5000
using namespace std;


void PrintD(int* d, int size)//печатает массив
{
	for (int i = 0; i < size; i++)
	{
		cout << d[i] << " ";
	}
	cout << endl;
}
void PrintG(int** G, int size)//печатает матрицу
{
	cout << "Matrix: " << endl;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			cout << G[i][j] << " ";
		cout << endl;
	}
}

int* InitOneLineGraph(int** G, int countEdge, int countVertex)//создает из матрицы массив
{
	int* oneG = new int[countVertex * countVertex];

	for (int i = 0, t = 0; i < countVertex; i++)
	{
		for (int j = 0; j < countVertex; j++)
		{
			oneG[t] = G[j][i];
			t++;
		}
	}
	return oneG;
}

int** IniGraph(int countEdge, int countVertex)//создает матрицу
{
	
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
	int t = 0;
	for (int i = 0; i < countVertex; i++) {
		for (int j = 0; j < countVertex; j++) {
			if (i == j)
			{
				G[i][j] = 0;

			}
			else
			{
				G[i][j] = rand() % 5;

				veroyatn = rand() % 100 + 1;
				if (G[i][j] == 0)
				{
					G[i][j] = INFINITI;
				}
			}
		}
	}

	return G;
}

int* ManuallyGraph(int countVertex)//задаем вручную
{
	cout << "Matrix" << endl;
	int rank = 0, procNum = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int* oneG = new int[countVertex * countVertex];
	int** G = new int* [countVertex];
	for (int i = 0; i < countVertex; i++)
	{
		G[i] = new int[countVertex];
	}
	//int t = 0;
	for (int i = 0; i < countVertex; i++)
	{
		cout << i << ": ";
		for (int j = 0; j < countVertex; j++)
		{
			cin >> G[i][j];
			if (i == j)
				G[i][j] = 0;
			cout << G[i][j] << " ";
		}
		cout << endl;
	}
	for (int i = 0, t = 0; i < countVertex; i++)
	{
		for (int j = 0; j < countVertex; j++)
		{
			oneG[t] = G[j][i];
			t++;
		}
	}
	for (int i = 0; i < countVertex; i++)
	{
		cout << i << ": ";
		for (int j = 0; j < countVertex; j++)
		{
			cout << G[i][j] << " ";
		}
		cout << endl;
	}
	PrintD(oneG, countVertex * countVertex);
	//cout << endl;
	return oneG;
}

int* Route(int countEdge, int countVertex)//граф путь (над диагональю)
{
	int rank = 0, procNum = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int* oneG = new int[countVertex * countVertex];
	int** G = new int* [countVertex];
	for (int i = 0; i < countVertex; i++)
	{
		G[i] = new int[countVertex];
	}
	//srand(time(NULL));
	int t = 0;
	for (int i = 0; i < countVertex; i++) {
		for (int j = 0; j < countVertex; j++) {
			if (i < j && i >= j - 1)
				G[j][i] = 1;
			else if (i == j)
				G[j][i] = 0;
			else
				G[j][i] = INFINITI;
			cout << G[j][i] << " ";
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
	PrintD(oneG, countVertex * countVertex);
	return oneG;
}

int* Round(int countEdge, int countVertex)//граф цикл (над диагональю и в начале последнего)
{
	int rank = 0, procNum = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int* oneG = new int[countVertex * countVertex];
	int** G = new int* [countVertex];
	for (int i = 0; i < countVertex; i++)
	{
		G[i] = new int[countVertex];
	}
	//srand(time(NULL));
	int t = 0;
	for (int i = 0; i < countVertex; i++) {
		for (int j = 0; j < countVertex; j++) {
			if (i < j && i >= j - 1)
				G[j][i] = 1;
			else if (i == j)
				G[j][i] = 0;
			else
				G[j][i] = INFINITI;
			G[0][countVertex - 1] = 1;
			cout << G[j][i] << " ";
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
	PrintD(oneG, countVertex * countVertex);
	return oneG;
}

int* Star(int countEdge, int countVertex)//граф звездочка (нулевой столбец)
{
	int rank = 0, procNum = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int* oneG = new int[countVertex * countVertex];
	int** G = new int* [countVertex];
	for (int i = 0; i < countVertex; i++)
	{
		G[i] = new int[countVertex];
	}
	//srand(time(NULL));
	int t = 0;
	for (int i = 0; i < countVertex; i++) {
		for (int j = 0; j < countVertex; j++) {
			if (i == 0 && i != j)
				G[j][i] = 1;
			else if (i == j)
				G[j][i] = 0;
			else
				G[j][i] = INFINITI;
			cout << G[j][i] << " ";
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
	PrintD(oneG, countVertex * countVertex);
	return oneG;
}

int* ReverseStar(int countEdge, int countVertex)//граф обратная звездочка(нулевой столбец и диангональ)
{
	int rank = 0, procNum = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int* oneG = new int[countVertex * countVertex];
	int** G = new int* [countVertex];
	for (int i = 0; i < countVertex; i++)
	{
		G[i] = new int[countVertex];
	}
	//srand(time(NULL));
	int t = 0;
	for (int i = 0; i < countVertex; i++) {
		for (int j = 0; j < countVertex; j++) {
			if (j == 0 && i != j)
				G[j][i] = 1;
			else if (i < j && i >= j - 1)
				G[j][i] = 1;
			else if (i == j)
				G[j][i] = 0;
			else
				G[j][i] = INFINITI;
			cout << G[j][i] << " ";
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
	PrintD(oneG, countVertex * countVertex);
	return oneG;
}

int* InitD(int size)
{
	int* d = new int[size];
	for (int i = 0; i < size; i++) {
		d[i] = INFINITI;
	}
	return d;
}

int* LinMoor(int* G, int start, int countVertex, vector<vector<int>>& vec) {
	double t1 = 0.0, t2 = 0.0, tTotal = 0.0;
	t1 = MPI_Wtime();

	vec.resize(countVertex);//сколько строчек в векторе векторов
	
	int* d = InitD(countVertex);//все infinity кроме стартовой
	d[start] = 0;
	queue<pair<int, int>>  q;
	q.push(make_pair(0, start));
	t2 = MPI_Wtime();
	tTotal += t2 - t1;

	for (int i = 0; i < vec.size(); i++)
	{
		vec[i].push_back(start);
	}

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
				len = G[j * countVertex + v];//находим вес
			if (d[v] + len < d[to])
			{
				d[to] = d[v] + len;
				vec[j] = vec[v];
				vec[j].push_back(j);
				q.push(make_pair(d[to], to));
			}
		}
	}

	return d;
}

int* ParallelMoor(int* G, int start, int countVertex, vector<vector<int>>& vec)
{
	int* d = nullptr;
	int* rowGraf = nullptr;
	int rank = 0, procNum = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int* distD = new int[procNum];//смещение
	int* distG = new int[procNum];//смещение
	int* countElementD = new int[procNum];//количество вершин каждому процессу
	int* countElementG = new int[procNum];//количество
	int partSize = countVertex / (procNum - 1);//часть для каждого процесса(колич вершин)
	countElementD[0] = 0;
	distD[0] = 0;
	countElementG[0] = 0;
	distG[0] = 0;
	for (int i = 1; i < procNum - 1; i++)
	{
		countElementD[i] = partSize;
		distD[i] = (i - 1) * partSize;

		countElementG[i] = countVertex * partSize;
		distG[i] = (i - 1) * countVertex * partSize;
	}
	////////
	countElementG[procNum - 1] = countVertex * (partSize + (countVertex % (procNum - 1)));//последнему добавляем то, что осталось(остаток)
	countElementD[procNum - 1] = (partSize + (countVertex % (procNum - 1)));//добавляем остаток по вершинам
	if (countElementG[procNum - 1] > countElementG[1])
	{
		int razn = (countElementG[procNum - 1] - countElementG[1]) / countVertex; //вынесли из цикла, т.к. countElG[0] меняется в цикле
																				//(количество вершин, на сколько больше в посл, чем в первом
		for (int i = 0; i < razn; i++)//раздаем из последнего другим
		{
			countElementG[i + 1] += countVertex;//количество весов из матрицы смежности
			countElementD[i + 1] ++;
			countElementD[procNum - 1] --;
			countElementG[procNum - 1] -= countVertex;
		}
		
	}
	for (int i = 1; i < procNum; i++)
	{
		distG[i] = distG[i - 1] + countElementG[i - 1];//смещение предыдущего и количество, которое достанется этому процессу
		distD[i] = distD[i - 1] + countElementD[i - 1];
	}

	int* partG = new int[countElementG[rank]];
	int* partD = new int[countElementD[rank]];
	int flag = 1, flagFinishFindCurrMinDistInfRank = 1;
	pair<int, int> currentVertex;
	pair<int, int> currtVertexWithCurrMinDest;
	pair<int, int> currentFlag;
	MPI_Status st;
	double t1SendingG = 0.0, t2SendingG = 0.0, t1SendingGTotal = 0.0, t1SendingD = 0.0, t2SendingD = 0.0, t1SendingDTotal = 0.0;
	double t1RecvingG = 0.0, t2RecvingG = 0.0, t1RecvingGTotal = 0.0, t1RecvingD = 0.0, t2RecvingD = 0.0, t1RecvingDTotal = 0.0;
	if (rank == 0)
	{			
		vec.resize(countVertex);//сколько строчек в векторе векторов
		for (int i = 0; i < vec.size(); i++)//ставим на начальную позицию стартовую вершину
		{
			vec[i].push_back(start);
		}		
		t1SendingG = MPI_Wtime();
		MPI_Scatterv(G, countElementG, distG, MPI_INT, partG, countElementG[rank], MPI_INT, 0, MPI_COMM_WORLD);//раздает процессам части матрицы
		t2SendingG = MPI_Wtime();
		t1SendingGTotal = t2SendingG - t1SendingG;
		d = InitD(countVertex);//заполнение d INF
		d[start] = 0;
		queue<pair<int, int>>  q;
		q.push(make_pair(0, start));//вставляет в очередь стартовую вершину
		int flagMaster = 1;
		while (!q.empty())
		{
			int v = q.front().second,//берет номер вершины
				dest = q.front().first;//берет вес
			q.pop();
			if (dest > d[v])
				continue;
			t1SendingD = MPI_Wtime();
			MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
			currentVertex = make_pair(dest, v);//создаем пару и отправляем ее
			MPI_Bcast(&currentVertex, 1, MPI_2INT, 0, MPI_COMM_WORLD);
			MPI_Scatterv(d, countElementD, distD, MPI_INT, partD, countElementD[rank], MPI_INT, 0, MPI_COMM_WORLD);
			t2SendingD = MPI_Wtime();
			t1SendingDTotal += t2SendingD - t1SendingD;
			for (int i = 1; i < procNum; i++) {
				flagMaster = 1;				
				while (flagMaster != 0) {//если == 0 значит !0 процесс завершил рассмотрение
					//cout << "rank = 0" << endl;
					t1RecvingD = MPI_Wtime();
					MPI_Recv(&flagFinishFindCurrMinDistInfRank, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &st);
					t2RecvingD = MPI_Wtime();
					t1RecvingDTotal += t2RecvingD - t1RecvingD;
					flagMaster = flagFinishFindCurrMinDistInfRank;
					if (flagMaster == 0) {

						continue;
					}
					t1RecvingD = MPI_Wtime();
					MPI_Recv(&currtVertexWithCurrMinDest, 1, MPI_2INT, i, 0, MPI_COMM_WORLD, &st);
					t2RecvingD = MPI_Wtime();
					t1RecvingDTotal += t2RecvingD - t1RecvingD;
					d[currtVertexWithCurrMinDest.second] = currtVertexWithCurrMinDest.first;//d[to] = d[v] + len;					
					q.push(make_pair(currtVertexWithCurrMinDest.first, currtVertexWithCurrMinDest.second));

					vec[currtVertexWithCurrMinDest.second] = vec[v];
					vec[currtVertexWithCurrMinDest.second].push_back(currtVertexWithCurrMinDest.second);
				}				
			}
			flagMaster = 1;
		}		
		flag = 0;//очередь закончилась, сообщаем другим
		MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		//cout << "time sending = " << t1SendingDTotal + t1SendingGTotal << endl;
		//cout << "time recving = " << t1RecvingDTotal + t1RecvingGTotal << endl;
		//PrintD(d, countVertex);

	}
	else {
		MPI_Scatterv(G, countElementG, distG, MPI_INT, partG, countElementG[rank], MPI_INT, 0, MPI_COMM_WORLD);
		
		while (flag) {
			MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if (flag == 0) {//когда очередь пустая
				continue;
			}
			MPI_Bcast(&currentVertex, 1, MPI_2INT, 0, MPI_COMM_WORLD);//принимает 
			MPI_Scatterv(d, countElementD, distD, MPI_INT, partD, countElementD[rank], MPI_INT, 0, MPI_COMM_WORLD);//принимает часть d
						
			for (int j = 0; j < countElementD[rank]; j++)
			{				
				int to = j;//вершина, до которой будем рассматривать
				int v = currentVertex.second;
				
				int len = partG[j * countVertex + v];//берет расстояние до вершины, умножаем, чтобы перейти на нужную строчку

				if (currentVertex.first + len < partD[to])
				{
					flagFinishFindCurrMinDistInfRank = 1;
					MPI_Send(&flagFinishFindCurrMinDistInfRank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
					currtVertexWithCurrMinDest.first = currentVertex.first + len;
					currtVertexWithCurrMinDest.second = to + distD[rank];			
					MPI_Send(&currtVertexWithCurrMinDest, 1, MPI_2INT, 0, 0, MPI_COMM_WORLD);
				}
			}			
			flagFinishFindCurrMinDistInfRank = 0;
			MPI_Send(&flagFinishFindCurrMinDistInfRank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}
	return d;
}

void IsCorrectImplementation(int* d1, int* d2, int size) {
	int countMistakes = 0;
	for (int i = 0; i < size; i++) {
		if (d1[i] != d2[i]) {
			cout << "dest isn`t correct" << endl;
			countMistakes++;
			break;
		}
	}
	if (countMistakes == 0)
		cout << "well done, dest is correct" << endl;
}

void IsCorrectPath(vector<vector<int>> vec1, vector <vector<int>> vec2, int size) {

	int countMistakes = 0;

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < vec1[i].size(); j++)
		{
			if (vec1[i][j] != vec2[i][j])
			{
				cout << "path isn`t correct" << endl;
				countMistakes++;
				break;
			}
		}
	}
	if (countMistakes == 0)
		cout << "well done, path is correct" << endl;
}
int main(int argc, char* argv[])
{
	int procNum = 0, rank = 0, countVertex = 0, countEdge = 0;
	MPI_Init(&argc, &argv);
	countVertex = atoi(argv[1]);
	//countVertex = 8;
	//countVertex = 5;
	//srand(time(NULL));
	int** Matrix = new int* [countVertex];//матрица
	Matrix = IniGraph(countEdge, countVertex);
	vector<vector<int>> vec1, vec2;
	countEdge = (countVertex - 1) + rand() % ((countVertex * (countVertex - 1)) / 2);//слева минимум, справа
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int* dStep = nullptr;
	int* dParallel = nullptr;
	int* G = nullptr;
	double t1Step = 0.0, t2Step = 0.0, t1Parallel = 0.0, t2Parallel = 0.0;
	int start = rand() % (countVertex - 1);
	start = atoi(argv[2]);
	if (rank == 0)
	{
		cout << "1.Auto\n2.Manually" << endl;
		int input;
		cin >> input;
		if (input == 2)
		{
			cout << "Start: " << endl;
			cin >> start;
			G = ManuallyGraph(countVertex);
		}
		else
		{
			cout << "1.Random\n2.Route\n3.Round\n4.Star\n5.ReverseStar" << endl;
			cin >> input;
			if (input == 1)
				G = InitOneLineGraph(Matrix, countEdge, countVertex);
			else if (input == 2)
				G = Route(countEdge, countVertex);
			else if (input == 3)
				G = Round(countEdge, countVertex);
			else if (input == 4)
				G = Star(countEdge, countVertex);
			else if (input == 5)
				G = ReverseStar(countEdge, countVertex);
		}
		//G = initG(countEdge, countVertex);
		if (countVertex < 11 && input == 1)
			PrintG(Matrix, countVertex);
	}
	//start = 0;
	//G = Route(countEdge, countVertex);
	if (procNum == 1) {
		//int* GStep = initG(countEdge, countVertex);
		t1Step = MPI_Wtime();
		//dStep = LinMoor(GStep, start, countVertex);
		dStep = LinMoor(G, start, countVertex, vec1);
		t2Step = MPI_Wtime();
		cout << "step by step implementation: " << endl;
		for (int i = 0; i < countVertex; i++)
			cout << dStep[i] << " ";
		cout << endl;
		cout << "count vertexs = " << countVertex << endl;	
		cout << "time of job = " << t2Step - t1Step << endl;
	}
	else
	{
		double par1, par2;
		if (!rank)//если не нулевой процесс
		{			
			t1Parallel = MPI_Wtime();
		}
		//par1 = MPI_Wtime();
		dParallel = ParallelMoor(G, start, countVertex, vec2);
		//par2 = MPI_Wtime();
		if (rank == 0)
		{
			t2Parallel = MPI_Wtime();	
			cout << "count vertexs = " << countVertex << endl;
			cout << "start = " << start << endl;
			cout << "time of job Parallel = " << t2Parallel - t1Parallel << endl;
			t1Step = MPI_Wtime();
			dStep = LinMoor(G, start, countVertex, vec1);
			t2Step = MPI_Wtime();
			cout << "time of job Lin = " << t2Step - t1Step << endl;
			cout << "Difference = " << (t2Parallel - t1Parallel) - (t2Step - t1Step) << endl;
			IsCorrectImplementation(dStep, dParallel, countVertex);
			if (countVertex < 20)
			{
				cout << "Dest S: ";
				for (int i = 0; i < countVertex; i++)
					cout << dStep[i] << " ";
				cout << endl;
				cout << "Dest P: ";
				for (int i = 0; i < countVertex; i++)
					cout << dParallel[i] << " ";
				cout << endl;
			}
			IsCorrectPath(vec1, vec2, countVertex);

			if (countVertex < 11)
				for (int i = 0; i < vec1.size(); i++)//линейный
				{				
					cout << "vertex " << i << ": ";
					for (int j = 0; j < vec1[i].size(); j++)
					{
						cout << vec1[i][j] << " ";
					}
					cout << endl;
				}

			//for (int i = 0; i < vec2.size(); i++)//линейный
			//{
			//	//vec[i].push_back(i);
			//	cout << "vertex " << i << ": ";
			//	for (int j = 0; j < vec2[i].size(); j++)
			//	{
			//		cout << vec2[i][j] << " ";
			//	}
			//	cout << endl;
			//}

		}
	}
	MPI_Finalize();
}