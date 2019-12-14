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

#define INFINITI 5000

void printD(int* d, int size) {
	for (int i = 0; i < size; i++) {
		cout << d[i] << " ";
	}
	cout << endl;
}
int* initG(int countEdge, int countVertex) {
	int rank = 0, procNum = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int* oneG = new int[countVertex * countVertex];
	int** G = new int* [countVertex];
	for (int i = 0; i < countVertex; i++) {
		G[i] = new int[countVertex];
	}
	srand(time(NULL));
	int t = 0;
	for (int i = 0; i < countVertex; i++) {
		for (int j = 0; j < countVertex; j++) {
			if (i == j) {
				G[i][j] = 0;

			}
			else {
				G[i][j] = rand() % 3;
				if (G[i][j] == 0) {
					G[i][j] = INFINITI;
				}
			}
			//cout << G[i][j] << " ";
		}
		//cout << endl;
	}
	for (int i = 0, t = 0; i < countVertex; i++) {
		for (int j = 0; j < countVertex; j++) {
			oneG[t] = G[j][i];
			t++;
		}
	}
	//printD(oneG, countVertex * countVertex);
	return oneG;
}
int* initD(int size) {
	int* d = new int[size];
	for (int i = 0; i < size; i++) {
		d[i] = INFINITI;
	}
	return d;
}
int* stepbystep(int* G, int start, int countVertex) {
	double t1 = 0.0, t2 = 0.0, tTotal = 0.0;
	t1 = MPI_Wtime();
	int* d = initD(countVertex);
	d[start] = 0;
	queue<pair<int, int>>  q;
	q.push(make_pair(0, start));
	t2 = MPI_Wtime();
	tTotal += t2 - t1;
	while (!q.empty()) {
		t1 = MPI_Wtime();
		int v = q.front().second, cur_d = -q.front().first;
		q.pop();
		if (cur_d > d[v])  continue;
		t2 = MPI_Wtime();
		tTotal += t2 - t1;
		for (int j = 0; j < countVertex; ++j) {
			int to = j,
				len = G[j * countVertex + v];
			if (d[v] + len < d[to]) {
				d[to] = d[v] + len;
				q.push(make_pair(-d[to], to));
			}
		}
	}
	cout << "time sending step by step = " << tTotal << endl;

	return d;
}
int* parallel(int* G, int start, int countVertex) {
	int* d = nullptr;
	int* rowGraf = nullptr;
	int rank = 0, procNum = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int* distD = new int[procNum];
	int* distG = new int[procNum];
	int* countElementD = new int[procNum];
	int* countElementG = new int[procNum];
	int partSize = countVertex / (procNum - 1);
	//int sizeCurrentPartD = countVertex / procNum ;
	countElementD[0] = 0;
	distD[0] = 0;
	countElementG[0] = 0;
	distG[0] = 0;
	for (int i = 1; i < procNum - 1; i++) {
		countElementD[i] = partSize;
		distD[i] = (i - 1) * partSize;

		countElementG[i] = countVertex * partSize;
		distG[i] = (i - 1) * countVertex * partSize;
	}
	distD[procNum - 1] = (procNum - 2) * partSize;
	countElementD[procNum - 1] = countVertex - (procNum - 2) * partSize;
	distG[procNum - 1] = (procNum - 2) * (countVertex * partSize);
	countElementG[procNum - 1] = countVertex * countVertex - (procNum - 2) * (countVertex * partSize);

	int* partG = new int[countElementG[rank]];
	int* partD = new int[countElementD[rank]];
	int flag = 1, flag_finish_find_current_min_distation_in_sigment_of_rank = 1;
	pair<int, int> currentVertex;
	pair<int, int> currentVertexWithCurrentMinDestation;
	pair<int, int> currentFlag;
	MPI_Status st;
	double t1SendingG = 0.0, t2SendingG = 0.0, t1SendingGTotal = 0.0, t1SendingD = 0.0, t2SendingD = 0.0, t1SendingDTotal = 0.0;
	double t1RecvingG = 0.0, t2RecvingG = 0.0, t1RecvingGTotal = 0.0, t1RecvingD = 0.0, t2RecvingD = 0.0, t1RecvingDTotal = 0.0;
	if (rank == 0) {
		//printD(G, countVertex*countVertex);
		t1SendingG = MPI_Wtime();
		MPI_Scatterv(G, countElementG, distG, MPI_INT, partG, countElementG[rank], MPI_INT, 0, MPI_COMM_WORLD);
		t2SendingG = MPI_Wtime();
		t1SendingGTotal = t2SendingG - t1SendingG;
		d = initD(countVertex);
		d[start] = 0;
		queue<pair<int, int>>  q;
		q.push(make_pair(0, start));
		//int flafRank0 = 1;
		int flafRank0 = 1;
		while (!q.empty()) {
			//flafRank0 = 1;
			int v = q.front().second,
				dest = -q.front().first;
			q.pop();
			if (dest > d[v])
				continue;
			t1SendingD = MPI_Wtime();
			MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
			currentVertex = make_pair(dest, v);
			MPI_Bcast(&currentVertex, 1, MPI_2INT, 0, MPI_COMM_WORLD);
			MPI_Scatterv(d, countElementD, distD, MPI_INT, partD, countElementD[rank], MPI_INT, 0, MPI_COMM_WORLD);
			t2SendingD = MPI_Wtime();
			t1SendingDTotal += t2SendingD - t1SendingD;
			for (int i = 1; i < procNum; i++) {
				//cout << "rank = 0" << endl;
				//flagFinishFindCurrentMinDistationInSigmentOfRank = 1;
				flafRank0 = 1;
				//cout << "rank = " << i << endl;
				while (flafRank0 != 0) {
					//cout << "rank = 0" << endl;
					t1RecvingD = MPI_Wtime();
					MPI_Recv(&flag_finish_find_current_min_distation_in_sigment_of_rank, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &st);
					t2RecvingD = MPI_Wtime();
					t1RecvingDTotal += t2RecvingD - t1RecvingD;
					flafRank0 = flag_finish_find_current_min_distation_in_sigment_of_rank;
					if (flafRank0 == 0) {

						continue;
					}
					t1RecvingD = MPI_Wtime();
					MPI_Recv(&currentVertexWithCurrentMinDestation, 1, MPI_2INT, i, 0, MPI_COMM_WORLD, &st);
					t2RecvingD = MPI_Wtime();
					t1RecvingDTotal += t2RecvingD - t1RecvingD;
					d[currentVertexWithCurrentMinDestation.second] = currentVertexWithCurrentMinDestation.first;
					//cout << " currntVetex d = " << currentVertex.first << " currentVertex = " << currentVertex.second << endl;
					q.push(make_pair(-currentVertexWithCurrentMinDestation.first, currentVertexWithCurrentMinDestation.second));
				}
				//currentVertex.second = 1;
			}
			flafRank0 = 1;
		}
		//cout << "end" << endl;
		flag = 0;
		MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		cout << "time sending = " << t1SendingDTotal + t1SendingGTotal << endl;
		cout << "time recving = " << t1RecvingDTotal + t1RecvingGTotal << endl;
		//printD(d, countVertex);
	}
	else {
		MPI_Scatterv(G, countElementG, distG, MPI_INT, partG, countElementG[rank], MPI_INT, 0, MPI_COMM_WORLD);
		//printD(partG, countElementG[rank]);
		while (flag) {
			MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if (flag == 0) {
				continue;
			}
			MPI_Bcast(&currentVertex, 1, MPI_2INT, 0, MPI_COMM_WORLD);
			MPI_Scatterv(d, countElementD, distD, MPI_INT, partD, countElementD[rank], MPI_INT, 0, MPI_COMM_WORLD);

			//printD(partD, countElementD[rank]);
			//printD(partG, countElementG[rank]);
			//cout << "rank = " << rank << "index = "<< 1*countVertex + currentVertex.second << " count Elements = "<<countElementG[rank] << endl;
			for (int j = 0; j < countElementD[rank]; j++) {
				//cout << " j = " << j << endl;
				int to = j,
					len = partG[j * countVertex + currentVertex.second];
				//cout <<"index = " << j * countVertex + currentVertex.second << "len = " << partG[j * countVertex + currentVertex.second] << "d[" << j << "] = " << partD[j] << endl;;
				//cout << "index = in circle" << partG[j*countVertex + currentVertex.second] << " currntVetex d = "<<currentVertex.first << " currentVertex = " << currentVertex.second <<endl;
				if (currentVertex.first + len < partD[to]) {
					flag_finish_find_current_min_distation_in_sigment_of_rank = 1;
					MPI_Send(&flag_finish_find_current_min_distation_in_sigment_of_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
					currentVertexWithCurrentMinDestation.first = currentVertex.first + len;
					currentVertexWithCurrentMinDestation.second = to + distD[rank];
					//cout << " to + distD = " << currentVertex.first;
					MPI_Send(&currentVertexWithCurrentMinDestation, 1, MPI_2INT, 0, 0, MPI_COMM_WORLD);
				}
			}
			//cout << endl;
			//cout <<"rank = "<<rank << endl;
			//currentVertex.second = -1;
			flag_finish_find_current_min_distation_in_sigment_of_rank = 0;
			MPI_Send(&flag_finish_find_current_min_distation_in_sigment_of_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}
	return d;
}
void isCorrectImplementation(int* d1, int* d2, int size) {
	int countMistakes = 0;
	for (int i = 0; i < size; i++) {
		if (d1[i] != d2[i]) {
			cout << " sorry, something isn`t correct" << endl;
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
	try 
	{
		countVertex = atoi(argv[1]);
	}
	catch (...)
	{
		cout << "ERROR INPUT DATA\n\n";
	}
	srand(time(NULL));
	countEdge = (countVertex - 1) + rand() % ((countVertex * (countVertex - 1)) / 2);
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int* dStep = nullptr;
	int* dParallel = nullptr;
	int* G = nullptr;
	double t1Step = 0.0, t2Step = 0.0, t1Parallel = 0.0, t2Parallel = 0.0;
	int start = rand() % (countVertex - 1);
	if (rank == 0)
		G = initG(countEdge, countVertex);
	if (procNum == 1) {
		int* GStep = initG(countEdge, countVertex);
		t1Step = MPI_Wtime();
		dStep = stepbystep(GStep, start, countVertex);
		t2Step = MPI_Wtime();
		cout << "step by step implementation: " << endl;
		cout << "count vertexs = " << countVertex << endl;
		//cout << "count edges = " << countEdge << endl;
		cout << "time of job = " << t2Step - t1Step << endl;
	}
	else {
		if (!rank) {
			//cout << "start = " << start << endl;
			t1Parallel = MPI_Wtime();
		}
		dParallel = parallel(G, start, countVertex);
		if (rank == 0) {
			t2Parallel = MPI_Wtime();
			cout << "parallel implementation: " << endl;
			cout << "count vertexs = " << countVertex << endl;
			cout << "start = " << start << endl;
			cout << "time of job = " << t2Parallel - t1Parallel << endl;
			//int** GStep = initGStep(countEdge, countVertex);
			t1Step = MPI_Wtime();
			dStep = stepbystep(G, start, countVertex);
			t2Step = MPI_Wtime();
			cout << "step by step implementation: " << endl;
			cout << "time of job = " << t2Step - t1Step << endl;
			isCorrectImplementation(dStep, dParallel, countVertex);
		}
	}
	return MPI_Finalize();
}
