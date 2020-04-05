#include <iostream>

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
	int pivot = arr[(f+l)/2];
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
	if(last>f)
		quickSort(arr, f, last);

}

int main(int argc, char* argv[])
{
	int* start_array = nullptr;
	int pivot = 0;
	int size = 0;
	bool check = false;

	cout << "Enter array size: ";
	cin >> size;

	if (size < 1) return 0;

	start_array = new int[size];

	for (int i = 0; i < size; i++)
		start_array[i] = rand() % 1000;

	for (int i = 0; i < size; i++)
		cout << start_array[i] << " ";
	cout << endl;
	quickSort(start_array, 0, size - 1);

	for (int i = 0; i < size; i++)
		cout << start_array[i] << " ";
	cout << endl;
	//проверка на правильность
		for (int i = 0; i < size; i++)
		{
			if (start_array[i] >= start_array[i - 1])
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
	delete start_array;

	return 0;
}