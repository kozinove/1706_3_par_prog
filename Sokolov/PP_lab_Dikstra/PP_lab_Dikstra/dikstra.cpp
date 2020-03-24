
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
	//omp_set_num_threads(NUM_THREAD);
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