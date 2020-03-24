#include <iostream>
#include <time.h>
using namespace std;

void Shell(int A[], int n) 
{
	int i, j, d;
	int count = 0;
	d = n;
	d = d / 2;
	while (d > 0) {
		for (i = 0; i < n - d; i++) {
			j = i;
			while (j >= 0 && A[j] > A[j + d]) {
				count = A[j];
				A[j] = A[j + d];
				A[j + d] = count;
				j--;
			}
		}
		d = d / 2;
	}
}
bool correctCheck(int A[], int size){

	for (int i = 0; i < size-1; i++) 
		if (A[i] > A[i + 1])
			return false;

	return true;
}

int main(int argc, char* argv[])
{
	int* massive = nullptr;
	int size;	 
	
	cout << "Enter array size: ";
	std::cin >> size;

	massive = new int[size]; 

	for (size_t i = 0; i < size; i++) {
		massive[i] = rand() % 1000; 
		cout << massive[i]<<" ";
	}
	cout << endl << endl;

	Shell(massive, size);

	for (size_t i = 0; i < size; i++) {
		cout << massive[i] << " ";
	}
	cout << endl << endl;

	cout << "Check correct: ";
	if (correctCheck(massive, size))
	 cout <<"True"  << endl;
	else cout << "False" << endl;

	delete[] massive;


	return 0;
}

