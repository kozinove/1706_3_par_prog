#pragma once
#include <iostream>
#include <queue>
#include <limits>       

#include <algorithm>
#include <functional>
#include <vector>

#include <omp.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>

using namespace tbb;

#define dataType int
#define NO_EBGE 0
#define NUM_THREAD 2
#define CHUNK 100


using std::priority_queue;
using std::pair;

using std::cout;
using std::numeric_limits;
using std::endl;

using std::vector;
using std::greater;

using std::make_heap;
using std::pop_heap;


//последовательный алгоритм Дейкстры на основе кучи
int dikstra(dataType** Gragh, dataType* Rez, int size);

// алгоритм Дейкстры на основе кучи c использованием openMP
int dikstraMP(dataType** Gragh, dataType* Rez, int size);

// алгоритм Дейкстры на основе кучи c использованием TBB
int dikstraTBB(dataType** Gragh, dataType* Rez, int size);