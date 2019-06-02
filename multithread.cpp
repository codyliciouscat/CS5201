#include <iostream>
#include <thread>
#include <deque>
#include <mutex>
#include <fstream>
#include <vector>
using namespace std;

mutex count_mutex;
mutex file_mutex;

/*
=======================================================
--------------------COUNT_OCCURENCE--------------------
======================================================= */
int count_occurence(const string file, const string input)
{
  int count = 0;
  string tmp;
  ifstream fin;

  fin.open(file.c_str());
  // runs through all words in file
  while(fin >> tmp)
  {
    // increments count if a match for input is found
    if(tmp == input)
      count++;
  }
  fin.close();

  return count;
}

/*
==============================================
--------------------MAPPER--------------------
============================================== */
void mapper(deque<string> & file_q, deque<int> & count_q, const string input)
{
  string file;
  int count;

  while(!file_q.empty())
  {
    /* CRITICAL SECTION */
    file_mutex.lock();
    file = file_q.front(); // set file equal to first element of file_q
    file_q.pop_front();
    file_mutex.unlock();

    count = count_occurence(file, input); // count number of times input appears in file

    /* CRITICAL SECTION */
    count_mutex.lock();
    count_q.push_back(count); // push count onto count_q
    count_mutex.unlock();
  }
}

/*
===============================================
--------------------REDUCER--------------------
=============================================== */
void reducer(deque<int> & count_q, deque<string> & file_q)
{
  while(count_q.size() != 1)
  {
      int a, b, sum;

      /* CRITICAL SECTION */
      count_mutex.lock();
      a = count_q.front(); // set a equal to first element of count_q
      count_q.pop_front();
      b = count_q.front(); // set b equal to second element of count_q
      count_q.pop_front();
      count_mutex.unlock();

      sum = a + b;

      /* CRITICAL SECTION */
      count_mutex.lock();
      count_q.push_back(sum); // push a + b onto count_q
      count_mutex.unlock();
  }
}

/*
============================================
--------------------MAIN--------------------
============================================ */
int main()
{
  string input;
  ifstream fin;
  const int num_mappers = 4, num_reducers = 1;
  int num_files;
  thread mappers[num_mappers], reducers[num_reducers];
  deque<string> file_q;
  deque<int> count_q;

  fin.open("files.dat");
  // reads in all files listed in files.dat and adds them to file_q
  while(fin >> input)
    file_q.push_back(input);
  fin.close();
  num_files = file_q.size();

  cout << endl << "Please enter word you wish to be counted: ";
  cin >> input;
  cout << endl;

  // launch mapper threads
  for(int i = 0; i < num_mappers; i++)
    mappers[i] = thread(mapper, ref(file_q), ref(count_q), input);
  // pauses until mappers finish
  for(int i = 0; i < num_mappers; i++)
    mappers[i].join();
  // launch reducer threads
  for(int i = 0; i < num_reducers; i++)
    reducers[i] = thread(reducer, ref(count_q), ref(file_q));
  // pauses until reducers finish
  for(int i = 0; i < num_reducers; i++)
    reducers[i].join();

  cout << "\"" << input << "\" appears " << count_q.front()
       << " times in " << num_files << " files" << endl << endl;

  return 0;
}
