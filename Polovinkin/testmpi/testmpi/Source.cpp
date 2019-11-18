#include <iostream>
#include "mpi.h"
#include <queue>
#include <windows.h>
using namespace std;

const int READER_OPEN = 0;
const int READER_CLOSE = 1;
const int WRITER_OPEN = 2;
const int WRITER_CLOSE = 3;

const int RESPONSE_READER_YES = 4;
const int RESPONSE_READER_WAIT = 5;
const int RESPONSE_WRITER_YES = 6;
const int RESPONSE_WRITER_WAIT = 7;



class Request
{
public:
	int rank;
	int step;
};

void RunReader(int rank, int time = 1000)
{
	
	Request request;
	int response;
	request.rank = rank;
	request.step = READER_OPEN;
    MPI_Status status;
	MPI_Send(&request, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);
	MPI_Recv(&response, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	if (response == RESPONSE_READER_WAIT)
	{
		MPI_Recv(&response, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		cout << "Reader" << rank << " wait " << endl;
	}
	cout << "Reader" << rank << " ready" << endl;
	Sleep(time);
	request.step = READER_CLOSE;
	MPI_Send(&request, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);
	MPI_Recv(&response, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	cout << "Reader" << rank << " finish" << endl;
}

void RunWriter(int rank, int time = 1000)
{
	
	Request request;
	int response;
	request.rank = rank;
	request.step = WRITER_OPEN;
    MPI_Status status;
	MPI_Send(&request, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);
	MPI_Recv(&response, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	if (response == RESPONSE_WRITER_WAIT)
	{
		MPI_Recv(&response, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		cout << "Writer" << rank << " wait" << endl;
	}
	cout << "Writer" << rank << " ready" << endl;
	Sleep(time);
	request.step = WRITER_CLOSE;
	MPI_Send(&request, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);
	MPI_Recv(&response, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	cout << "Writer" << rank << " finish" << endl;
}
int main(int argc, char**argv)
{
	int rank;
	int size;
	const int k = 2;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;
	if (rank == 0)
	{
		int ReaderCount = 0;
		int WriterCount = 0;
		queue <int> q_readers, q_writers;
		for (int i = 0; i < (size-1)*k * 2; i++)
		{
			MPI_Status status;
			Request request;
			MPI_Recv(&request, 2, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
			switch (request.step)
			{
			case READER_OPEN:
			{
				if (WriterCount == 0)
				{
					ReaderCount++;
					cout << "Reader" << request.rank << " open - yes\n" << endl;
					MPI_Send(&RESPONSE_READER_YES, 1, MPI_INT, request.rank, 0, MPI_COMM_WORLD);
				}
				else
				{
					cout << "Reader" << request.rank << " open - no\n" << endl;
					q_readers.push(request.rank);
					MPI_Send(&RESPONSE_READER_WAIT, 1, MPI_INT, request.rank, 0, MPI_COMM_WORLD);
				}
				break;
			}
			case READER_CLOSE:
			{
				ReaderCount--;
				cout << "Reader" << request.rank << " close - yes\n" << endl;
				MPI_Send(&RESPONSE_READER_YES, 1, MPI_INT, request.rank, 0, MPI_COMM_WORLD);
				if (q_writers.size() > 0 && ReaderCount == 0)
				{
					WriterCount++;
					MPI_Send(&RESPONSE_WRITER_YES, 1, MPI_INT, q_writers.front(), 0, MPI_COMM_WORLD);
					q_writers.pop();

				}
				break;
			}
			case WRITER_OPEN:
			{
				if (ReaderCount == 0 && WriterCount == 0)
				{
					WriterCount++;
					cout << "Writer" << request.rank << " open - yes\n" << endl;
					MPI_Send(&RESPONSE_WRITER_YES, 1, MPI_INT, request.rank, 0, MPI_COMM_WORLD);
				}
				else
				{
					cout << "Writer" << request.rank << " open - no\n" << endl;
					q_writers.push(request.rank);
					MPI_Send(&RESPONSE_WRITER_WAIT, 1, MPI_INT, request.rank, 0, MPI_COMM_WORLD);
				}
				break;
			}
			case WRITER_CLOSE:
			{
				WriterCount--;
				cout << "Writer" << request.rank << " close - yes\n" << endl;
				MPI_Send(&RESPONSE_WRITER_YES, 1, MPI_INT, request.rank, 0, MPI_COMM_WORLD);

				if (q_writers.size() > 0 && ReaderCount == 0)
				{
					WriterCount++;
					MPI_Send(&RESPONSE_WRITER_YES, 1, MPI_INT, q_writers.front(), 0, MPI_COMM_WORLD);
					q_writers.pop();

				}
				else if (q_readers.size() > 0)
				{
					int q_size = q_readers.size();
					for (int i = 0; i < q_size; i++)
					{
						ReaderCount++;
						MPI_Send(&RESPONSE_READER_YES, 1, MPI_INT, q_readers.front(), 0, MPI_COMM_WORLD);
						q_readers.pop();
					}
				}
				break;
			}
			}
		}
	}
	
	if ( rank % 2 == 0 && rank != 0)
	{
		for (int i = 0; i <k; i++)
		{
			RunReader(rank);
		}
	}
	if ( rank % 2 == 1)
	{
		for (int i = 0; i < k; i++)
		{
			RunWriter(rank);
		}
	}
	MPI_Finalize();
	return 0;
}