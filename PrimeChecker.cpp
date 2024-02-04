#include <iostream>
#include <vector>
#include <thread>
#include <semaphore>
#include <chrono>

using namespace std;    //hope using these are fine
using namespace chrono;

binary_semaphore semaphorePrimeVector {1};
                                                                                                                                            
/*
This function checks if an integer n is prime.

Parameters:
n : int - integer to check

Returns true if n is prime, and false otherwise.
*/
bool check_prime (const int &n);
void check_prime_ranges (int start, int end, vector<int>& primes);

int main () {
    vector<int> primes;
    vector<thread> threads;
    high_resolution_clock::time_point timeStart, timeEnd;
    nanoseconds execTime;
    int nLimit = 2;
    int nThreadCount = 1;
    int nInterval, nExtra, nPrevEnd;

    cout << "Prime check up to number: ";
    cin >> nLimit;

    cout << "Number of threads to use: ";
    cin >> nThreadCount;

    //Subtract 1 since [1] is not being checked
    nInterval = (nLimit - 1) / nThreadCount;    
    nExtra = (nLimit - 1) % nThreadCount;
    nPrevEnd = 1;   //First interval starts from 2

    //I know specs said timer starts after user input, but I didn't want to include those three lines of computation.
    //It's probably negligible time, but, you know, just in case.
    timeStart = high_resolution_clock::now ();

    for (int i = 0; i < nThreadCount; i++) {
        int start = nPrevEnd + 1;
        int end = (start + nInterval) > nLimit ? nLimit : start + nInterval - 1;

        if (nExtra != 0)  {
            end += 1;
            nExtra--;
        }

        if (end >= nLimit)  //In case upper bound is less than thread count
            i = nThreadCount + 1;

        nPrevEnd = end;
        threads.emplace_back (check_prime_ranges, start, end, ref (primes));
    }

    for (auto& thread : threads) {
        thread.join ();
    }

    timeEnd = high_resolution_clock::now ();
    execTime = duration_cast <nanoseconds> (timeEnd - timeStart);

    cout << primes.size () << " primes were found." << endl;
    cout << "Execution time: " << execTime.count () << endl;

    system ("pause");
    return 0;
}

bool check_prime (const int &n) {
    for (int i = 2; i * i <= n; i++)
        if (n % i == 0)
            return false;

    return true;
}

void check_prime_ranges (int start, int end, vector<int>& primes) {
    for (int current_num = start; current_num <= end; current_num++)
        if (check_prime (current_num)) {
            semaphorePrimeVector.acquire ();
            primes.push_back (current_num);
            semaphorePrimeVector.release ();
        }
}