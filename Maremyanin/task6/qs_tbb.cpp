#include "tbb\task_scheduler_init.h"
#include "tbb\parallel_for.h"
#include "tbb\blocked_range.h"
#include "tbb\tbb.h"
#include <time.h>
#include <iostream>
#include <cstdio>


using namespace tbb;

void swap(int* a, int* b)
{
	int temp = *a;
	*a = *b;
	*b = temp;
}
void quickSort(int arr[], int first, int last)
{
	int f = first;
	int l = last;
	int pivot = arr[(f + l) / 2];
	while (f <= l)
	{
		while (arr[f] < pivot)
			f++;
		while (arr[l] > pivot)
			l--;
		if (f <= l)
			swap(&arr[f++], &arr[l--]);
	}
	if (first < l)
		quickSort(arr, first, l);
	if (last > f)
		quickSort(arr, f, last);

}

class SeparationEven : public task
{
private:
	int* arr;
	int size1, size2;

public:
	SeparationEven(int* _arr, int _size1, int _size2) : arr(_arr), size1(_size1), size2(_size2) {
	}

	task* execute()
	{
		int* arr2 = arr + size1;

		int num = (size1 + size2 + 1) / 2;
		int* tmp = new int[num];

		int a = 0, b = 0, i = 0;
		while (a < size1 && b < size2) {
			if (arr[a] <= arr2[b])
			{
				tmp[i] = arr[a];
				a += 2;
			}
			else
			{
				tmp[i] = arr2[b];
				b += 2;
			}
			i++;
		}

		if (a >= size1)
			for (int j = b; j < size2; j += 2, i++)
				tmp[i] = arr2[j];
		else
			for (int j = a; j < size1; j += 2, i++)
				tmp[i] = arr[j];

		
		for (int j = 0; j < num; ++j)
			arr[j * 2] = tmp[j];
		return NULL;
	}
};

class SeparationOdd : public task
{
private:
	int* arr;
	int size1, size2;

public:
	SeparationOdd(int* _arr, int _size1, int _size2) : arr(_arr), size1(_size1), size2(_size2) {
	}

	task* execute()
	{
		int* arr2 = arr + size1;

		int num = (size1 + size2) - (size1 + size2 + 1) / 2;
		int* tmp = new int[num];

		int a = 1, b = 1, i = 0;
		while (a < size1 && b < size2) {
			if (arr[a] <= arr2[b])
			{
				tmp[i] = arr[a];
				a += 2;
			}
			else
			{
				tmp[i] = arr2[b];
				b += 2;
			}
			i++;
		}

		if (a >= size1)
			for (int j = b; j < size2; j += 2, i++)
				tmp[i] = arr2[j];
		else
			for (int j = a; j < size1; j += 2, i++)
				tmp[i] = arr[j];

		
		for (int j = 0; j < num; ++j)
			arr[j * 2 + 1] = tmp[j];
		return NULL;
	}
};


class Checker
{
private:
	int* arr;
public:
	Checker(int* _arr) : arr(_arr) {}

	void operator()(const blocked_range<int>& r) const
	{
		int begin = r.begin(), end = r.end();
		for (int i = begin; i < end; i++)
			if (arr[i - 1] > arr[i])
			{
				int tmp = arr[i - 1];
				arr[i - 1] = arr[i];
				arr[i] = tmp;
			}
	}
};

class Sorting : public task
{
private:
	int* arr;
	int size;
	int parts;

public:
	Sorting	(int* _arr, int _size, int _parts) : arr(_arr), size(_size), parts(_parts) {}

	task* execute()
	{
		if (size <= parts)
		{
			quickSort(arr,0,size-1);
		}
		else
		{
			int s = size / 2 + (size / 2) % 2;
			Sorting& sorting1 = *new (allocate_child()) Sorting(arr, s, parts);
			Sorting& sorting2 = *new (allocate_child()) Sorting(arr + s, size - s, parts);
			set_ref_count(3);
			spawn(sorting1);
			spawn_and_wait_for_all(sorting2);
			SeparationEven& separationEven = *new (allocate_child()) SeparationEven(arr, s, size - s);
			SeparationOdd& separationOdd = *new (allocate_child()) SeparationOdd(arr, s, size - s);
			set_ref_count(3);
			spawn(separationEven);
			spawn_and_wait_for_all(separationOdd);

			parallel_for(blocked_range<int>(1, size), Checker(arr));
		}
		return NULL;
	}
};

void quickSortTBB(int* arr, int size, int threads)
{
	int parts = size / threads + 1;
	Sorting& sorting = *new (task::allocate_root())
		Sorting(arr, size, parts);
	task::spawn_root_and_wait(sorting);
}


int main(int argc, char* argv[])
{
	int* arr = nullptr;
	int pivot = 0;
	int size = 0;
	bool check = false;
	int threads;
	std :: cout << "Enter array size: ";
	std :: cin >> size;
	std::cout << "Enter number of processes: ";
	std::cin >> threads;
	if (size < 1) return 0;
	if (threads < 1) return 0;
	arr = new int[size];
	for (int i = 0; i < size; i++)
		arr[i] = rand() % 1000;
	if (size <= 100)
	{
		for (int i = 0; i < size; i++)
			std:: cout << arr[i] << " ";
		std:: cout << std:: endl;
	}
	double start = clock();
	quickSortTBB(arr, size, threads);
	double end = clock();
	double time = end - start;
	std::cout << "Time tbb: " << time / 1000 << std::endl;
	if (size <= 100)
	{
		for (int i = 0; i < size; i++)
			std::cout << arr[i] << " ";
		std::cout << std::endl;
	}
	for (int i = 0; i < size; i++) // проверка на правильность
	{
		if (arr[i] >= arr[i - 1])
			continue;
		else
		{
			check = true;
			break;
		}
	}
	if (check == true)
		std::cout << "Array isnt sorted";
	else
		std::cout << "Array is sorted";
	delete arr;
	return 0;
}
