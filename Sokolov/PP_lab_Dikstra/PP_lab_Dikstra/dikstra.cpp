
#include "pch.h"
#include "dikstra.h"


//последовательный алгоритм Дейкстры на основе кучи
int dikstra(dataType** Gragh, dataType* Rez, int size)
{
	vector< pair<dataType, int>> Marks(size);

	Marks[0] = pair<dataType, int>(0, 0);
	for (int i = 1; i < size; i++)
	{
		Marks[i] = pair<dataType, int>(numeric_limits<dataType>::max(), i);
	}

	make_heap(Marks.begin(), Marks.end(), std::greater<>{});

	/*for (auto i = Marks.begin(); i != Marks.end(); ++i)
	{
		cout << (*i).first << " " << (*i).second << endl;
	}*/



	while (!Marks.empty())
	{
		pop_heap(Marks.begin(), Marks.end(), std::greater<>{});
		pair<dataType, int> vertex = Marks.back();
		Marks.pop_back();
		Rez[vertex.second] = vertex.first;

		int size_vec = Marks.size();
		for (int change_vertex_ined = 0; change_vertex_ined < size_vec; ++change_vertex_ined)
		{
			if (Gragh[vertex.second][Marks[change_vertex_ined].second] != NO_EBGE)
			{
				dataType delta = Marks[change_vertex_ined].first - (vertex.first + Gragh[vertex.second][Marks[change_vertex_ined].second]);
				if (delta > 0)
				{
					Marks[change_vertex_ined].first -= delta;
				}
			}
		}
	}
	return 0;
}

// алгоритм Дейкстры на основе кучи c использованием openMP
int dikstraMP(dataType** Gragh, dataType* Rez, int size)
{
	omp_set_num_threads(NUM_THREAD);
	vector< pair<dataType, int>> Marks(size);


	Marks[0] = pair<dataType, int>(0, 0);
	int i;
	int chunk = CHUNK;

#pragma omp parallel shared(Marks) private(i) 
	{
	//	if (omp_in_parallel())
	//		cout << omp_get_num_threads() << endl;
#pragma omp for schedule(static)
		for (i = 1; i < size; i++)                                                                          //here
		{
			Marks[i] = pair<dataType, int>(numeric_limits<dataType>::max(), i);
		}
	}

	make_heap(Marks.begin(), Marks.end(), std::greater<>{});

	/*for (auto i = Marks.begin(); i != Marks.end(); ++i)
	{
		cout << (*i).first << " " << (*i).second << endl;
	}*/

	while (!Marks.empty())
	{
		pop_heap(Marks.begin(), Marks.end(), std::greater<>{});
		pair<dataType, int> vertex = Marks.back();
		Marks.pop_back();
		Rez[vertex.second] = vertex.first;

		int size_vec = Marks.size();
		int change_vertex_ined;

#pragma omp parallel shared(Marks, vertex, size_vec, Gragh) private(change_vertex_ined) 
		{
		//	if (omp_in_parallel())
		//		cout << omp_get_num_threads() << endl;
#pragma omp for schedule(static)
			for (change_vertex_ined = 0; change_vertex_ined < size_vec; ++change_vertex_ined)                             //here
			{
				if (Gragh[vertex.second][Marks[change_vertex_ined].second] != NO_EBGE)
				{
					dataType delta = Marks[change_vertex_ined].first - (vertex.first + Gragh[vertex.second][Marks[change_vertex_ined].second]);
					if (delta > 0)
					{
						Marks[change_vertex_ined].first -= delta;
					}
				}
			}
		}
	}
	return 0;
}


class SetDistance //Функтор
{
	vector< pair<dataType, int>>* vecMarks;
public:
	SetDistance(vector< pair<dataType, int>> *tVecMarks) : vecMarks(tVecMarks)
	{}
	void operator()(const blocked_range<int>& r) const
	{
		int begin = r.begin(), end = r.end();
		for (int i = begin; i != end; i++)
			(*vecMarks)[i] = pair<dataType, int>(numeric_limits<dataType>::max(), i);
	}
};

class CalculateChangeDistanse //Функтор
{
	dataType** Gragh;

	vector< pair<dataType, int>>* vecMarks; 
	pair<dataType, int> *vertex;

public:
	CalculateChangeDistanse(dataType** tGragh, vector< pair<dataType, int>> *tMarks, pair<dataType, int> *tvertex) : Gragh(tGragh), vecMarks(tMarks), vertex(tvertex)
	{}
	
	void operator()(const blocked_range<int>& r) const
	{
		int begin = r.begin(), end = r.end();
		for (int change_vertex_ined = begin; change_vertex_ined != end; ++change_vertex_ined)
		{
			if (Gragh[(*vertex).second][(*vecMarks)[change_vertex_ined].second] != NO_EBGE)
			{
				dataType delta = (*vecMarks)[change_vertex_ined].first - ((*vertex).first + Gragh[(*vertex).second][(*vecMarks)[change_vertex_ined].second]);
				if (delta > 0)
				{
					(*vecMarks)[change_vertex_ined].first -= delta;
				}
			}
		}
	}
};

// алгоритм Дейкстры на основе кучи c использованием TBB
int dikstraTBB(dataType** Gragh, dataType* Rez, int size)
{

	vector< pair<dataType, int>> Marks(size);

	Marks[0] = pair<dataType, int>(0, 0);
	int i;
	

	task_scheduler_init init;
	parallel_for(blocked_range<int>(1, size), SetDistance(&Marks)); // у блок ранге есть 3 аргумент - размер порции

	init.terminate();

	make_heap(Marks.begin(), Marks.end(), std::greater<>{});
	
	/*for (auto i = Marks.begin(); i != Marks.end(); ++i)
	{
		cout << (*i).first << " " << (*i).second << endl;
	}*/

	init.initialize();

	while (!Marks.empty())
	{
		pop_heap(Marks.begin(), Marks.end(), std::greater<>{});
		pair<dataType, int> vertex = Marks.back();
		Marks.pop_back();
		Rez[vertex.second] = vertex.first;

		int size_vec = Marks.size();
		int change_vertex_ined;

		parallel_for(blocked_range<int>(0, size_vec), CalculateChangeDistanse(Gragh, &Marks, &vertex));// у блок ранге есть 3 аргумент - размер порции

		
		//	for (change_vertex_ined = 0; change_vertex_ined < size_vec; ++change_vertex_ined)//here
		//{
		//	if (Gragh[vertex.second][Marks[change_vertex_ined].second] != NO_EBGE)
		//	{
		//		dataType delta = Marks[change_vertex_ined].first - (vertex.first + Gragh[vertex.second][Marks[change_vertex_ined].second]);
		//		if (delta > 0)
		//		{
		//			Marks[change_vertex_ined].first -= delta;
		//		}
		//	}
		//}
	}

	return 0;
}






void InitFunction(vector< pair<dataType, int>> &Marks, int pos, int length) {
	for (int i = pos; i < (length + pos); i++)
	{
		Marks[i] = pair<dataType, int>(numeric_limits<dataType>::max(), i);
	}
}

std::mutex Thread_mutex[NUM_THREAD];

void JoinFunction(dataType** Gragh, vector< pair<dataType, int>> &Marks, pair<dataType, int> &vertex, int index, std::vector<bool> &Thread_barier_work, bool& exit) {
	while (!exit)
	{
		if (Thread_barier_work[index] == true)
		{
			int pos[NUM_THREAD];
			int length[NUM_THREAD];

			int _t_size = Marks.size() / NUM_THREAD;

			for (int i = 0; i < NUM_THREAD; i++)
			{
				length[i] = _t_size;
			}
			for (int i = 0; i < (Marks.size() % NUM_THREAD); i++)
			{
				length[i] += 1;
			}


			pos[0] = 0;
			for (int i = 1; i < NUM_THREAD; i++)
			{
				pos[i] = pos[i - 1] + length[i - 1];
			}

			for (int change_vertex_ined = pos[index]; change_vertex_ined < (length[index] + pos[index]); ++change_vertex_ined)
			{
				if (Gragh[vertex.second][Marks[change_vertex_ined].second] != NO_EBGE)
				{
					dataType delta = Marks[change_vertex_ined].first - (vertex.first + Gragh[vertex.second][Marks[change_vertex_ined].second]);
					if (delta > 0)
					{
						Marks[change_vertex_ined].first -= delta;
					}
				}
			}

			Thread_barier_work[index] = false;
		}
	}
}

// алгоритм Дейкстры на основе кучи c использованием Std Thread
int dikstraStdThread(dataType** Gragh, dataType* Rez, int size)
{
	int N = std::thread::hardware_concurrency();
	
	vector< pair<dataType, int>> Marks(size);

	int length[NUM_THREAD];
	int pos[NUM_THREAD];
	for (int i = 0; i < NUM_THREAD; i++)
	{
		length[i] = (size-1) / NUM_THREAD;
	}
	for (int i = 0; i < ((size-1) % NUM_THREAD); i++)
	{
		length[i] += 1;
	}
	pos[0] = 1;
	for (int i = 1; i < NUM_THREAD; i++)
	{
		pos[i] = pos[i - 1] + length[i - 1];
	}

	Marks[0] = pair<dataType, int>(0, 0);

	std::vector<std::thread> vecThread;
	for (int i = 1; i < NUM_THREAD; i++)
	{
		vecThread.push_back(std::thread(InitFunction, std::ref(Marks), std::ref(pos[i]), std::ref(length[i])));
	}

	for (int i = 1; i < (length[0] + 1); i++)
	{
		Marks[i] = pair<dataType, int>(numeric_limits<dataType>::max(), i);
	}

	for (auto& t: vecThread)
	{
		t.join();
	}
	vecThread.clear();
	
	make_heap(Marks.begin(), Marks.end(), std::greater<>{});

	//for (auto i = Marks.begin(); i != Marks.end(); ++i)
	//{
	//	cout << (*i).first << " " << (*i).second << endl;
	//}

	std::vector<bool> Thread_barier_work(NUM_THREAD);
	for (bool flag : Thread_barier_work)
	{
		flag = false;
	}

	bool exit = false;

	pair<dataType, int> vertex;

	for (int i = 1; i < NUM_THREAD; i++)
	{
		vecThread.push_back(std::thread(JoinFunction, std::ref(Gragh), std::ref(Marks), std::ref(vertex), i, std::ref(Thread_barier_work), std::ref(exit)));
	}

	while (!Marks.empty())
	{
		pop_heap(Marks.begin(), Marks.end(), std::greater<>{});
		vertex = Marks.back();
		Marks.pop_back();
		Rez[vertex.second] = vertex.first;
		
		//start thread
		for (int i = 1; i < NUM_THREAD; i++)
		{
			Thread_barier_work[i] = true;
		}

		int length = Marks.size() / NUM_THREAD;
		if (Marks.size() % NUM_THREAD != 0)
			length += 1;
		
		for (int change_vertex_ined = 0; change_vertex_ined < (length + 0); ++change_vertex_ined)
		{
			if (Gragh[vertex.second][Marks[change_vertex_ined].second] != NO_EBGE)
			{
				dataType delta = Marks[change_vertex_ined].first - (vertex.first + Gragh[vertex.second][Marks[change_vertex_ined].second]);
				if (delta > 0)
				{
					Marks[change_vertex_ined].first -= delta;
				}
			}
		}

		//wait end work thread
		for (int i = 1; i < NUM_THREAD; i++)
		{
			bool flag = true;
			while (flag)
			{
				if (Thread_barier_work[i] == false)
				{
					flag = false;
				}
			}
		}
	}
	exit = true;
	for (auto& t : vecThread)
	{
		t.join();
	}
	return 0;
}
