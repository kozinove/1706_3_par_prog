#include <iostream>
#include <time.h>
#include <thread>
#include <vector>

using namespace std;

struct part
{
	int* array;
	unsigned size;
};

int* merge(int* a, int* b, int n, int m)
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

void ref_merge(int*& c, int* a, int* b, int n, int m)
{
	c = merge(a, b, n, m);
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

int* parallel_merge(int numtasks, part* parts);


int main(int argc, char* argv[])
{
	vector<thread> threads;
	int* start_array = nullptr;
	int* sorted_start_array = nullptr;
	int* array = nullptr;
	int* remainder_array = nullptr;
	int* result_array = nullptr;
	clock_t start;
	clock_t end;
	clock_t start_single;
	clock_t end_single;
	int numtasks;
	int size;
	int remainder = 0;

	int** parts;

	cout << "Enter number of threads: ";
	std::cin >> numtasks;

	cout << "Enter array size: ";
	std::cin >> size;

	start_array = new int[size];
	sorted_start_array = new int[size];

	for (size_t i = 0; i < size; i++)
		start_array[i] = rand() % 1000;
	for (size_t i = 0; i < size; i++)
		sorted_start_array[i] = start_array[i];

	start_single = clock();
	quickSort(sorted_start_array, 0, size - 1);
	end_single = clock();

	remainder = size % numtasks;
	size = size - remainder;


	array = new int[size];

	if (remainder)
		remainder_array = new int[remainder];

	int i = 0;

	for (; i < size; i++)
		array[i] = start_array[i];

	for (int j = 0; j < remainder; j++, i++)
		remainder_array[j] = start_array[i];

	if (remainder > 1)
		quickSort(remainder_array, 0, remainder - 1);

	delete[] start_array;

	int part_size = size / numtasks;
	parts = new int* [numtasks];
	for (size_t i = 0; i < numtasks; i++)
	{
		parts[i] = new int[part_size];
		int shift = i * part_size;
		for (int j = 0; j < part_size; j++)
		{
			parts[i][j] = array[j + shift];
		}
	}

	start = clock();

	for (size_t i = 0; i < numtasks; i++)
	{
		threads.push_back(thread(quickSort, ref(parts[i]), 0, part_size - 1));
	}

	for (size_t i = 0; i < numtasks; i++)
	{
		threads.at(i).join();
	}

	threads.erase(threads.begin(), threads.begin() + numtasks);

	if (numtasks == 2)
	{
		if (remainder)
		{
			result_array = merge(remainder_array, merge(parts[0], parts[1], part_size, part_size), remainder, size);
			size += remainder;
		}
		else result_array = merge(parts[0], parts[1], part_size, part_size);
	}
	else if (numtasks == 3)
	{
		if (remainder)
		{
			result_array = merge(remainder_array, merge(parts[2], merge(parts[0], parts[1], part_size, part_size),
				part_size, part_size + part_size), remainder, size);
			size += remainder;
		}
		else result_array = merge(parts[2], merge(parts[0], parts[1], part_size, part_size), part_size, part_size + part_size);
	}
	else if (numtasks < 2)
	{
		result_array = parts[0];
	}
	else
	{
		int current_parts = numtasks / 2 + numtasks % 2;
		part* sorted_parts = new part[current_parts];

		int _numtasks = numtasks;
		int _current_parts = current_parts;

		for (size_t i = 0; i < numtasks / 2; i++)
		{
			sorted_parts[i].array = new int[part_size * 2];
			sorted_parts[i].size = part_size * 2;
		}
		if (numtasks % 2)
		{
			sorted_parts[current_parts - 1].size = part_size;
			sorted_parts[current_parts - 1].array = new int[part_size];
			for (size_t i = 0; i < part_size; i++)
			{
				sorted_parts[current_parts - 1].array[i] = parts[numtasks - 1][i];
			}

			_numtasks--;
			_current_parts--;
		}

		for (size_t i = 0; i < _current_parts; i++)
		{
			int part_id = ((i + 1) * 2) - 1;
			threads.push_back(thread(ref_merge, ref(sorted_parts[i].array), ref(parts[part_id]), ref(parts[part_id - 1]), part_size, part_size));
		}

		for (size_t i = 0; i < _current_parts; i++)
		{
			threads.at(i).join();
		}

		threads.erase(threads.begin(), threads.begin() + _current_parts);

		for (size_t i = 0; i < numtasks; i++)
		{
			delete[] parts[i];
		}
		delete[] parts;
		if (remainder)
		{
			result_array = merge(remainder_array, parallel_merge(current_parts, sorted_parts), remainder, size);
			size += remainder;
		}
		else result_array = parallel_merge(current_parts, sorted_parts);
	}

	end = clock();

	bool same = true;
	for (size_t i = 0; i < size; i++)
	{
		if (sorted_start_array[i] != result_array[i])
		{
			same = false;
			break;
		}

	}

	cout << "Arrays check: ";
	if (same) cout << "PASS!" << endl;
	else cout << "ERROR!" << endl;

	double seconds_single = (double)(end_single - start_single) / CLOCKS_PER_SEC;
	std::cout << "The time of single thread quick sort " << seconds_single << " seconds." << std::endl;
	double seconds = (double)(end - start) / CLOCKS_PER_SEC;
	std::cout << "The time of " << numtasks << "x threads quick sort " << seconds << " seconds." << std::endl;
	std::cout << std::endl << "Result: " << seconds_single / seconds << "x acceleration" << std::endl;


	delete[] sorted_start_array;
	delete[] array;
	delete[] remainder_array;
	delete[] result_array;

	return 0;
}


int* parallel_merge(int numtasks, part* parts)
{
	if (numtasks == 2)
	{
		return merge(parts[0].array, parts[1].array, parts[0].size, parts[1].size);
	}
	else if (numtasks == 3)
	{
		return merge(parts[2].array, merge(parts[0].array, parts[1].array, parts[0].size, parts[1].size),
			parts[2].size, parts[0].size + parts[1].size);
	}
	else if (numtasks < 2) return nullptr;

	int current_parts = numtasks / 2 + numtasks % 2;
	part* sorted_parts = new part[current_parts];

	int _numtasks = numtasks;
	int _current_parts = current_parts;
	int part_size = parts[0].size;

	for (size_t i = 0; i < numtasks / 2; i++)
	{
		sorted_parts[i].array = new int[part_size * 2];
		sorted_parts[i].size = part_size * 2;
	}
	if (numtasks % 2)
	{
		sorted_parts[current_parts - 1].size = parts[numtasks - 1].size;
		sorted_parts[current_parts - 1].array = new int[sorted_parts[current_parts - 1].size];
		for (size_t i = 0; i < sorted_parts[current_parts - 1].size; i++)
		{
			sorted_parts[current_parts - 1].array[i] = parts[numtasks - 1].array[i];
		}

		_numtasks--;
		_current_parts--;
	}


	vector<thread> threads;

	for (size_t i = 0; i < _current_parts; i++)
	{
		int part_id = ((i + 1) * 2) - 1;
		threads.push_back(thread(ref_merge, ref(sorted_parts[i].array), ref(parts[part_id].array), ref(parts[part_id - 1].array), parts[part_id].size, parts[part_id - 1].size));
	}

	
	for (size_t i = 0; i < _current_parts; i++)
	{
		threads.at(i).join();
	}

	threads.erase(threads.begin(), threads.begin() + _current_parts);


	for (size_t i = 0; i < numtasks; i++)
	{
		delete[] parts[i].array;
	}

	return parallel_merge(current_parts, sorted_parts);
}
