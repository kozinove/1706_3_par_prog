#include <iostream>
#include <vector>

#include <mpi.h>

int ROOT = 0;

template <class T>
T DoOperation(T** array, int id, int size, MPI_Op op)
{
	T var;

	switch (op)
	{
	case MPI_MAX:
	{
		var = array[0][id];
		for (int i = 0; i < size; i++)
		{
			if (array[i][id] > var)
				var = array[i][id];
		}
		break;
	}
	case MPI_MIN:
	{
		var = array[0][id];
		for (int i = 0; i < size; i++)
		{
			if (array[i][id] < var)
				var = array[i][id];
		}
		break;
	}
	case MPI_SUM:
	{
		var = 0;
		for (int i = 0; i < size; i++)
		{
			var += array[i][id];
		}
		break;
	}
	default:
		break;
	}

	return var;
}


template <class T>
void DoAllreduce(
	T* sendbuf,
	T* recvbuf,
	int count,
	_In_ MPI_Datatype datatype,
	_In_ MPI_Op op,
	_In_ MPI_Comm comm
)
{
	int numtasks, taskid;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	T** array;

	if (numtasks == 1)
	{
		if (taskid == ROOT)
		{
			for (int i = 0; i < count; i++)
			{
				recvbuf[i] = sendbuf[i];
			}
		}
		return;
	}

	if (count == 1)
	{
		if (taskid != ROOT)
		{
			MPI_Send(&(*sendbuf), count, datatype, ROOT, 1, comm);
			MPI_Recv(&(*recvbuf), count, datatype, ROOT, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
		}
		else
		{
			array = new T * [numtasks];
			for (int i = 0; i < numtasks; i++)
				array[i] = new T[count];
			array[ROOT][0] = *(sendbuf);

			for (int i = 0; i < numtasks; i++)
				if (i != ROOT) MPI_Recv(&(array[i][0]), count, datatype, i, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);

			*(recvbuf) = DoOperation<T>(array, 0, numtasks, op);

			for (int i = 0; i < numtasks; i++)
				if (i != ROOT) MPI_Send(&(*recvbuf), count, datatype, i, 0, comm);

			for (int i = 0; i < numtasks; i++)
				delete[] array[i];
			delete[] array;
		}
		return;
	}

	if (taskid != ROOT)
	{
		MPI_Send(sendbuf, count, datatype, ROOT, 1, comm);
		MPI_Recv(recvbuf, count, datatype, ROOT, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);
	}
	else
	{
		array = new T * [numtasks];
		for (int i = 0; i < numtasks; i++)
			array[i] = new T[count];
		array[ROOT][0] = *(sendbuf);

		for (int i = 0; i < count; i++)
			array[ROOT][i] = sendbuf[i];

		for (int i = 0; i < numtasks; i++)
			if (i != ROOT) MPI_Recv(array[i], count, datatype, i, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);

		for (int i = 0; i < count; i++)
			recvbuf[i] = DoOperation<T>(array, i, numtasks, op);

		for (int i = 0; i < numtasks; i++)
			if (i != ROOT) MPI_Send(recvbuf, count, datatype, i, 0, comm);

		for (int i = 0; i < numtasks; i++)
			delete[] array[i];
		delete[] array;
	}
}



MPI_METHOD MPI_Allreduce(
	_In_range_(!= , recvbuf) _In_opt_ const void* sendbuf,
	_Out_opt_ void* recvbuf,
	_In_range_(>= , 0) int count,
	_In_ MPI_Datatype datatype,
	_In_ MPI_Op op,
	_In_ MPI_Comm comm
)
{
	if (datatype == MPI_INT)
		DoAllreduce<int>((int*)sendbuf, (int*)recvbuf, count, datatype, op, comm);
	else if (datatype == MPI_DOUBLE)
		DoAllreduce<double>((double*)sendbuf, (double*)recvbuf, count, datatype, op, comm);
	else if (datatype == MPI_FLOAT)
		DoAllreduce<float>((float*)sendbuf, (float*)recvbuf, count, datatype, op, comm);
	else return MPI_ERR_TYPE;


	return MPI_SUCCESS;
}


int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	int numtasks, taskid;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
	
	bool status;
	if (taskid == 0)
	{
		std::cout << "Set ROOT process: ";
		std::cin >> ROOT;

		for (int i = 1; i < numtasks; i++)
			MPI_Send(&ROOT, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Recv(&ROOT, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	if (ROOT >= numtasks || ROOT < 0) return -1;

	int count = 5;
	int* sendbuf,* recvbuf;

	sendbuf = new int[count];
	recvbuf = new int[count];

	for (int i = 0; i < count; i++)
	{
		sendbuf[i] = i;
		recvbuf[i] = 0;
	}

	MPI_Allreduce(sendbuf, recvbuf, count, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	for (int i = 0; i < count; i++)
		std::cout << "id: " << taskid << " recvbuf[" << i << "] = " << recvbuf[i] << std::endl;


	delete[] recvbuf;
	delete[] sendbuf;

	MPI_Finalize();
}