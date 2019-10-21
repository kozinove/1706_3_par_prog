#include <iostream>
#include <vector>

#include <mpi.h>

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);

	int* vector = nullptr;
	int vector_size = 1000;

	int numtasks, taskid;

	int global_min	=	INT_MAX, 
		local_min	=	INT_MAX, 
		real_min	=	INT_MAX;

	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

	if (taskid == 0)
	{
		vector = new int[vector_size];
		for (size_t i = 0; i < vector_size; i++)
		{
			vector[i] = rand();
			if (vector[i] < real_min)
			{
				real_min = vector[i];
			}
		}
	}

	int receive_vector_size = vector_size / numtasks;
	int remainder = vector_size % numtasks;
	if (remainder) receive_vector_size += 1;
	int* receive_vector = new int[receive_vector_size];

	MPI_Scatter(vector, receive_vector_size, MPI_INT, receive_vector, receive_vector_size, MPI_INT, 0, MPI_COMM_WORLD);

	if (receive_vector_size * taskid + 1 <= vector_size)
	{
		if (vector_size - receive_vector_size * (taskid) < receive_vector_size)
		{
			receive_vector_size = vector_size - (taskid) * (receive_vector_size);
		}

		std::cout << "taskid " << taskid << ": receive_vector_size = " << receive_vector_size << std::endl;

		local_min = receive_vector[0];
		for (size_t i = 0; i < receive_vector_size; i++)
		{
			if (local_min > receive_vector[i]) local_min = receive_vector[i];
		}
	}

	MPI_Reduce(&local_min, &global_min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

	std::cout << "taskid " << taskid << ": local_min = " << local_min << std::endl;
	if (taskid == 0)
	{
		std::cout << "real_min = " << real_min <<
		std::endl << "global_min = " <<	global_min;

		delete[] vector;
	}

	delete[] receive_vector;

	MPI_Finalize();
}
