#include <iostream>
#include <omp.h>
#include <vector>
using namespace std;

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
void SeparationOddEven(int chetnost, int* first_arr, int size1, int* second_arr, int size2,vector<int>& parts) // выбирает четные или нечетные элементы и записывает в вектор 
{
	int i, j;
	if (chetnost == 0)
	{
		i = 0, j = 0;
	}
	else
	{
		i = 1, j = 1;
	}
	parts.reserve(size1 + size2);
	while (i < size1 && j < size2) 
	{
		if (first_arr[i] <= second_arr[j])
		{
			parts.push_back(first_arr[i]);
			i += 2;
		}
		else
		{
			parts.push_back(second_arr[j]);
			j += 2;
		}
	}
	if (i >= size1)
		while (j < size2) 
		{
			parts.push_back(second_arr[j]);
			j += 2;
		}
	else
		while (i < size1) 
		{
			parts.push_back(first_arr[i]);
			i += 2;
		}
}
void Merger(vector<int> first_vec,vector<int> second_vec, int* result) // сливает вектора нечетные и четных значений
{
	int i = 0, j = 0;
	int size1 = first_vec.size(), size2 = second_vec.size();

	while (i < size1 && j < size2) 
	{
		result[i + j] = first_vec[i];
		result[i + j + 1] = second_vec[j];
		++i; ++j;
	}

	while (i < size1) 
	{
		result[size2 + i] = first_vec[i];
		i++;
	}
	while (j < size2) 
	{
		result[size1 + j] = second_vec[j];
		j++;
	}

	i = 1; 
	while (i < size1 + size2 - 1) 
	{
		if (result[i] > result[i + 1])
		{
			j = result[i];
			result[i] = result[i + 1];
			result[i + 1] = j;
		}
		++i;
	}
}


int main(int argc, char* argv[])
{
	int* arr = nullptr;
	int pivot = 0;
	int size = 0;
	bool check = false;
	int numtasks;
	cout << "Enter array size: ";
	cin >> size;
	cout << "Enter number of processes: ";
	cin >> numtasks;


	if (size < 1) return 0;
	if (numtasks < 1) return 0;
	arr = new int[size];

	for (int i = 0; i < size; i++)
		arr[i] = rand() % 1000;
	if (size <= 100)
	{
		for (int i = 0; i < size; i++)
			cout << arr[i] << " ";
		cout << endl;
	}
	omp_set_num_threads(numtasks);
	double time = omp_get_wtime();
	int step;
	vector<int>* temparr = new vector<int>[numtasks];
	int* sdvig = new int[numtasks];
	int* fragments_size = new int[numtasks];
#pragma omp parallel shared(arr, step, sdvig, fragments_size, temparr) num_threads(numtasks) 
	{	
		int task_id;
		int thread_pair;
		task_id = omp_get_thread_num();
		
		sdvig[task_id] = task_id * (size / numtasks);
		fragments_size[task_id] = (task_id == numtasks - 1) ? size - task_id * (size / numtasks) : size / numtasks;
		quickSort(arr + sdvig[task_id], 0, fragments_size[task_id] - 1);
#pragma omp barrier

		step = 1;
		while (step < numtasks)
		{
			thread_pair = (int)pow(2, step - 1);

			if (task_id % (thread_pair * 2) == 0)
			{
				SeparationOddEven(0, arr + sdvig[task_id], fragments_size[task_id], arr + sdvig[task_id + thread_pair], fragments_size[task_id + thread_pair], temparr[task_id]);
			}
			else if (task_id % thread_pair == 0)
			{
				SeparationOddEven(1, arr + sdvig[task_id], fragments_size[task_id], arr + sdvig[task_id - thread_pair], fragments_size[task_id - thread_pair], temparr[task_id]);
			}

#pragma omp barrier

			if (task_id % (thread_pair * 2) == 0)
			{
				Merger(temparr[task_id], temparr[task_id + thread_pair], arr + sdvig[task_id]);
				fragments_size[task_id] += fragments_size[task_id + thread_pair];
				temparr[task_id].clear();
				temparr[task_id].shrink_to_fit();
				temparr[task_id + thread_pair].clear();
				temparr[task_id + thread_pair].shrink_to_fit();
			}
#pragma omp single 
			{
				step *= 2;
			}
#pragma omp barrier
		}
	}

	delete[] temparr;
	delete[] fragments_size;
	delete[] sdvig;

	time = omp_get_wtime() - time;
	cout << "Parallel time:" << time << endl;
	if (size <= 100)
	{
		for (int i = 0; i < size; i++)
			cout << arr[i] << " ";
		cout << endl;
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
		cout << "Array isnt sorted";
	else
		cout << "Array is sorted";
	delete arr;

	return 0;
}