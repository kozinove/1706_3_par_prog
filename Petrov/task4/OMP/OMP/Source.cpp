#include <iostream>
#include <time.h>
#include <omp.h>
#include <vector>
using namespace std;

vector<int>c;

void ShellSort(vector<int>& A, int n)
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

bool correctCheck(vector<int> A, vector<int> B) {

    for (int i = 0; i < A.size(); i++)
        if (A[i] != B[i])
            return false;
    return true;
}

vector<int> merge(vector<int> a, vector<int> b, vector<int> c, int n, int m)
{

    int size = n + m;
    c.resize(size);

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

vector<int> recursiveMerge(int numOfThreads, vector<vector<int>>& part) {

    if (numOfThreads == 2)
    {
        return merge(part[0], part[1], c, part[0].size(), part[1].size());
    }
    else if (numOfThreads < 2) return part[0];

    int parts = numOfThreads / 2;
    vector<vector<int>> sortedParts(parts);
    int partSize = part[0].size();

    for (int i = 0; i < numOfThreads / 2; i++) {
        sortedParts[i].resize(part[0].size() * 2);
    }

    omp_set_num_threads(parts);
#pragma omp parallel shared(parts)
    {
        int threadNum = omp_get_thread_num();
        int partNum = threadNum + 1;
        sortedParts[threadNum] = merge(part[(partNum * 2) - 1], part[(partNum * 2) - 2], c, part[(partNum * 2) - 1].size(), part[(partNum * 2) - 2].size());
    }

    return recursiveMerge(parts, sortedParts);
}



int main(int argc, char* argv[])
{

    int k = 0;
    double startParallel;
    double endParallel;
    double parallelTime;
    double startLinear;
    double endLinear;
    double linearTime;
    int numOfThreads;
    int size;
    
    cout << "Enter array size:";
    cin >> size;
    cout << "Enter number of threads:";
    cin >> numOfThreads;
    vector<int> mainVector(size);
    vector<int> mainFinalVector;
    omp_set_num_threads(numOfThreads);
    int remain = size % numOfThreads;
    vector<int>remainVector(remain);

    for (size_t i = 0; i < size; i++) {
        mainVector[i] = rand() % 100;
    }

    vector<vector<int>> finalVector(numOfThreads / 2);
    vector<vector<int>> threadVector(numOfThreads);

    size -= remain;
    int partSize = size / numOfThreads;

    for (int i = 0; i < numOfThreads / 2; i++) {
        finalVector[i].resize(partSize);
    }
    for (int i = 0; i < numOfThreads; i++) {
        vector<int> part(partSize);
        for (int j = 0; j < partSize; k++, j++) {
            part[j] = mainVector[k];
        }
        threadVector[i] = part;
    }
    for (int i = 0; i < remain; i++) {
        remainVector[i] = mainVector[size + i];
    }
    cout << endl;
    startLinear = omp_get_wtime();
    ShellSort(mainVector, size + remain);
    endLinear = omp_get_wtime();
    linearTime = endLinear - startLinear;
    cout << "Linear time: " << linearTime << endl << endl;
    startParallel = omp_get_wtime();



#pragma omp parallel for shared(threadVector)
    for (int i = 0; i < numOfThreads; i++)
        ShellSort(threadVector[i], partSize);

    if (remain)
        ShellSort(remainVector, remain);


    mainFinalVector = recursiveMerge(numOfThreads, threadVector);
    if (remain) mainFinalVector = merge(mainFinalVector, remainVector, c, mainFinalVector.size(), remain);
    endParallel = omp_get_wtime();
    parallelTime = endParallel - startParallel;
    cout << "Parallel time: " << parallelTime << endl << endl;

    double acc;
    acc = linearTime / parallelTime;
    cout << "Aceleration increased by " << acc << " times." << endl;
    cout << endl << endl;
    cout << "Check correct: ";
    if (correctCheck(mainFinalVector, mainVector))
        cout << "True" << endl;
    else cout << "False" << endl;
    system("pause");
    return 0;
}

