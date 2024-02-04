#include <iostream>
#include <utility>
#include <vector>
#include <chrono>
#include <algorithm>
#include <random>
#include <thread>
#include <queue>
#include <semaphore>

using namespace std;
using namespace chrono;

typedef pair<int, int>  ii;
#define SEED 31

binary_semaphore semaphoreQueue {1};
binary_semaphore semaphoreFinished {1};

/*
This function generates all the intervals for merge sort iteratively, given the 
range of indices to sort. Algorithm runs in O(n).

Parameters:
start : int - start of range
end : int - end of range (inclusive)

Returns a list of integer pairs indicating the ranges for merge sort.
*/
vector<ii> generate_intervals (int start, int end);

/*
This function performs the merge operation of merge sort.

Parameters:
array : vector<int> - array to sort
s     : int         - start index of merge
e     : int         - end index (inclusive) of merge
*/
void merge (vector<int> &array, int s, int e);

/*
This function is the sort of "main()" of the threads. Keeps looping while there are jobs unfinished. Pulls a job if unemployed.

Parameters:
intervalsQueue : queue<ii>   - queue of ready jobs
intervals      : vector<ii>  - list of unreleased intervals
numbers        : vector<int> - array to sort
finishedPairs  : vecotr<ii>  - list of intervals that are finished
*/
void threadMoment (queue<ii> &intervalsQueue, vector<ii> &intervals, vector<int> &numbers, vector<ii> &finishedPairs);

/*
This function checks if an interval is ready to be processed i.e., their required halves are finished merging.

Parameters:
interval      : ii         - interval to be checked
finishedPairs : vector<ii> - list of intervals that are finished.
*/
bool check_if_doable (ii interval, vector<ii> &finishedPairs);
/*
This function performs the sanity operation of my sanity
i.e., checking if the array was sorted correctly (ascending)

Parameters:
array : vector<int> - array that is allegedly sorted
size  : int         - size of the array
*/
bool sanity_check (vector<int> &array, int size);

int main () {
    int random, nRand;
    int nSize, nPower, nThreadCount;
    int i;
    high_resolution_clock::time_point timeStart, timeEnd;
    milliseconds execTime;

    // TODO: Seed your randomizer
    srand (SEED);
    nRand  = rand ();


    // TODO: Get array size and thread count from user
    cout << "Input array size (power of 2): ";
    cin >> nPower;

    cout << "Number of threads to use: ";
    cin >> nThreadCount;

    //Get the array size: 2^nPower
    nSize = 2;
    for (i = 1; i < nPower; i++)
        nSize *= 2;

    vector<thread> threadPool;
    threadPool.reserve (nSize);

    // TODO: Generate a random array of given size
    vector<int> numbersV;
    numbersV.reserve (nSize);

    for (i = 0; i < nSize; i++)
        numbersV.push_back (i + 1);

    shuffle (begin (numbersV), end (numbersV),
                    default_random_engine (nRand));
    
    // TODO: Call the generate_intervals method to generate the merge sequence
    vector<ii> intervals = generate_intervals(0, nSize - 1);
    queue<ii> intervalQueue;    //Queue to store intervals that are ready to be processed

    for (i = 0; i < nSize; i++) {
        intervalQueue.push (intervals.at (0));
        intervals.erase (intervals.begin ());
    }

    //Vector to store finished pairs
    vector<ii> finishedPairs;

    cout << endl << "Beginning sort" << endl << endl;
    
    //Get time started
    timeStart = high_resolution_clock::now ();
    
    //Instantiate threads
    for (i = 0; i < nThreadCount; i++)
        threadPool.emplace_back (threadMoment, ref (intervalQueue),
                                               ref (intervals),
                                               ref (numbersV),
                                               ref (finishedPairs));

    while (!intervals.empty ()) {
        for (i = 0; i < intervals.size (); i++) {
            if (check_if_doable (intervals.at (i), ref (finishedPairs))) {
                semaphoreQueue.acquire ();

                intervalQueue.push (intervals.at (i));
                intervals.erase (intervals.begin () + i);

                semaphoreQueue.release ();
            }
        }
    }

    for (auto& thread : threadPool) {
        thread.join ();
    }

    timeEnd = high_resolution_clock::now ();
    execTime = duration_cast <milliseconds> (timeEnd - timeStart);

    if (sanity_check (numbersV, numbersV.size ()))
        cout << "Array is correctly sorted" << endl;
    else cout << "NOT SORTED INSANE TIME" << endl;

    cout << "Time it took: " << execTime.count () << endl;

    system ("pause");
    return 0;
}

vector<ii> generate_intervals (int start, int end) {
    vector<ii> frontier;
    frontier.push_back (ii (start, end));
    int i = 0;

    while (i < (int)frontier.size ()) {
        int s = frontier[i].first;
        int e = frontier[i].second;

        i++;

        if (s == e)
            continue;

        int m = s + (e - s) / 2;

        frontier.push_back (ii (m + 1, e));
        frontier.push_back (ii (s, m));
    }

    vector<ii> retval;

    for (i = (int)frontier.size () - 1; i >= 0; i--) {
        retval.push_back (frontier[i]);
    }

    return retval;
}

void merge (vector<int> &array, int s, int e) {
    int m = s + (e - s) / 2;
    vector<int> left, right;
    int i, l_ptr, r_ptr;

    for (i = s; i <= e; i++) {
        if (i <= m) {
            left.push_back (array[i]);
        } else {
            right.push_back (array[i]);
        }
    }

    l_ptr = r_ptr = 0;

    for (i = s; i <= e; i++) {
        if (l_ptr == (int)left.size ()) {
            array[i] = right[r_ptr];
            r_ptr++;
        } else if (r_ptr == (int)right.size () || left[l_ptr] <= right[r_ptr]) {
            array[i] = left[l_ptr];
            l_ptr++;
        } else {
            array[i] = right[r_ptr];
            r_ptr++;
        }
    }
}

void threadMoment (queue<ii> &intervalsQueue, vector<ii> &intervals, vector<int> &numbersV, vector<ii> &finishedPairs) {
    bool die = false;
    
    while (!die) {
        semaphoreQueue.acquire ();

        if (intervalsQueue.empty ()) {
            if (intervals.empty ())
                die = true;
            semaphoreQueue.release ();
            continue;
        }

        ii interval = intervalsQueue.front ();
        intervalsQueue.pop ();
        semaphoreQueue.release ();

        merge (numbersV, interval.first, interval.second);

        semaphoreFinished.acquire ();
        finishedPairs.push_back (interval);
        semaphoreFinished.release ();
    }
}

bool sanity_check (vector<int> &array, int size) {
    int i;
    for (i = 0; i < size - 1; i++)
        if (array[i] > array[i + 1])
            return false;

    return true;
}

bool check_if_doable (ii interval, vector<ii> &finishedPairs) {
    int i;

    semaphoreFinished.acquire ();
    if (finishedPairs.empty ()) {
        semaphoreFinished.release ();
        return false;
    }

    int finishedSize = finishedPairs.size ();
    semaphoreFinished.release ();

    ii preLeft, preRight;
    int gap = (interval.second - interval.first) / 2;
    bool foundLeft, foundRight;
    foundLeft = foundRight = false;

    preLeft = {interval.first, interval.first + gap};
    preRight = {interval.first + gap + 1, interval.second};

    for (i = 0; i < finishedSize; i++) {
        if (finishedPairs.at (i).first == preLeft.first && finishedPairs.at (i).second == preLeft.second)
            foundLeft = true;
        else if (finishedPairs.at (i).first == preRight.first && finishedPairs.at (i).second == preRight.second)
            foundRight = true;

        if (foundLeft && foundRight) return true;
    }

    return false;
}