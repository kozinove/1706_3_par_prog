#include <iostream>
#include <queue>
#include <mpi.h>
#include <time.h>

using namespace std;

int* merge(int *a, int *b, int n, int m)
{
	int* c;
	int size = n + m;
	c = new int[size];

	int i = 0, j = 0, k = 0;

	while (i < n && j < m)
	{
		if (a[i] <= b[j])
		{
			c[k] = a[i];
			i++;
		}
		else
		{
			c[k] = b[j];
			j++;
		}
		k++;
	}
	while (i < n)
	{
		c[k] = a[i];
		k++;
		i++;
	}
	while (j < m)
	{
		c[k] = b[j];
		k++;
		j++;
	}
	return c;
}

void swap(int* a, int* b)
{
	int t = *a;
	*a = *b;
	*b = t;
}

int partition(int arr[], int low, int high)
{
	int pivot = arr[high];
	int i = (low - 1);

	for (int j = low; j <= high - 1; j++)
	{
		if (arr[j] < pivot)
		{
			i++;
			swap(&arr[i], &arr[j]);
		}
	}
	swap(&arr[i + 1], &arr[high]);
	return (i + 1);
}

void quickSort(int arr[], int low, int high)
{
	if (low < high)
	{
		int pi = partition(arr, low, high);

		quickSort(arr, low, pi - 1);
		quickSort(arr, pi + 1, high);
	}
}


int main(int argc, char* argv[])
{
	int* start_array = nullptr;
	int* sorted_start_arry = nullptr;
	int* array = nullptr;
	int* remainder_array = nullptr;
	clock_t start;
	clock_t end;
	clock_t start_single;
	clock_t end_single;

	int size = 1000000;
	int remainder = 0;
	
	MPI_Init(&argc, &argv);
	int numtasks, taskid;
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &taskid);


	if (taskid == 0)
	{
		cout << "Enter array size: ";
		std::cin >> size;

		for (size_t i = 1; i < numtasks; i++)
			MPI_Send(&size, 1, MPI_INT, i, 1, MPI_COMM_WORLD);

		start_array = new int[size];
		sorted_start_arry = new int[size];

		for (size_t i = 0; i < size; i++)
			start_array[i] = rand() % 1000;
		for (size_t i = 0; i < size; i++)
			sorted_start_arry[i] = start_array[i];

		start_single = clock();
		quickSort(sorted_start_arry, 0, size - 1);
		end_single = clock();
	}
	else
	{
		MPI_Recv(&size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
	}

	remainder = size % numtasks;
	size = size - remainder;

	if (taskid == 0)
	{
		array = new int[size];

		if (remainder)
			remainder_array = new int[remainder];

		int i = 0;

		for (; i < size; i++)
			array[i] = start_array[i];

		for (; i < remainder; i++)
			remainder_array[i] = start_array[i];

		if (remainder)
			quickSort(remainder_array, 0, remainder - 1);

		cout << endl;
	}

	int part_size = size / numtasks;
	int recieve_array_size = size / numtasks;
	int* recieve_array = new int[recieve_array_size];

	if (taskid == 0) start = clock();
	MPI_Scatter(array, recieve_array_size, MPI_INT, recieve_array, recieve_array_size, MPI_INT, 0, MPI_COMM_WORLD);

	quickSort(recieve_array, 0, recieve_array_size - 1);

	int* recieve_sorted_array = new int[recieve_array_size];
	int recieve_sorted_array_size = recieve_array_size;

	int* sorted_array = nullptr;

	int step = 2;
	int post_receive = 0;
	std::queue<int> post_receive_sizes;
	std::queue<int> post_receive_tasks;

	if (taskid % 2) MPI_Send(recieve_array, recieve_array_size, MPI_INT, taskid - 1, 0, MPI_COMM_WORLD);
	else if (!(taskid == numtasks - 1 && !(numtasks % 2)))
	{
		MPI_Recv(recieve_sorted_array, recieve_array_size, MPI_INT, taskid + 1, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		sorted_array = merge(recieve_array, recieve_sorted_array, recieve_array_size, recieve_array_size);

		recieve_array_size *= 2;
		while (step*2 <= numtasks)
		{
			delete[] recieve_sorted_array;
			delete[] recieve_array;

			recieve_array = sorted_array;

			if ((numtasks / step) % 2)
			{
				if (taskid == 0)
				{
					post_receive++;
					post_receive_sizes.push(recieve_array_size);
					post_receive_tasks.push(numtasks - step);
				}
				else if (taskid == numtasks - step)
				{
					MPI_Send(recieve_array, recieve_array_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
					break;
				}
			}

			if ((taskid / step)%2)
			{
				MPI_Send(recieve_array, recieve_array_size, MPI_INT, taskid - step, 0, MPI_COMM_WORLD);
			}
			else
			{
				recieve_sorted_array = new int[recieve_array_size];
				MPI_Recv(recieve_sorted_array, recieve_array_size, MPI_INT, taskid + step, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
				sorted_array = merge(recieve_array, recieve_sorted_array, recieve_array_size, recieve_array_size);
				recieve_array_size *= 2;
			}

			step *= 2;
			if (taskid != 0 && (taskid % step)) break;
		}
	}
	else
	{
		MPI_Send(recieve_array, recieve_array_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	if (taskid == 0 && !numtasks % 2)
	{
		delete[] recieve_sorted_array;
		delete[] recieve_array;
		recieve_array = sorted_array;

		recieve_sorted_array = new int[part_size];
		MPI_Recv(recieve_sorted_array, part_size, MPI_INT, numtasks - 1, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		sorted_array = merge(recieve_array, recieve_sorted_array, recieve_array_size, part_size);
	}
	while (post_receive)
	{
		post_receive--;
		delete[] recieve_sorted_array;
		delete[] recieve_array;
		recieve_array = sorted_array;

		recieve_sorted_array = new int[post_receive_sizes.front()];
		MPI_Recv(recieve_sorted_array, post_receive_sizes.front(), MPI_INT, post_receive_tasks.front(), 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		sorted_array = merge(recieve_array, recieve_sorted_array, recieve_array_size, post_receive_sizes.front());
		post_receive_tasks.pop();
		post_receive_sizes.pop();
	}
	
	if (taskid == 0)
	{
		if (remainder) {
			delete[] recieve_array;
			recieve_array = sorted_array;
			sorted_array = merge(recieve_array, remainder_array, size, remainder);
			delete[] remainder_array;
		}
		if (taskid == 0) end = clock();

		bool succes = true;
		for (int i = 0; i < size + remainder; i++)
		{
			if (sorted_array[i] != sorted_start_arry[i])
			{
				succes = false;
				break;
			}
		}
		std::cout << "Arrays check pass: ";
		if (succes) std::cout << "true!" << std::endl;
		else std::cout << "false!" << std::endl;

		double seconds_single = (double)(end_single - start_single) / CLOCKS_PER_SEC;
		std::cout << "The time of single process quick sort " << seconds_single << " seconds." << std::endl;
		double seconds = (double)(end - start) / CLOCKS_PER_SEC;
		std::cout << "The time of" << numtasks << "x process quick sort " << seconds << " seconds." << std::endl;

		std::cout << std::endl << "Result: " << seconds_single/seconds << "x acceleration" << std::endl;

		delete[] array;
	}

	MPI_Finalize();
	return 0;
}
