/////////////////////////////////////////////////////////////////////
// Methods for sorting arrays of data.                             //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 11/06/2019 - created - Brennan Young                            //
// 12/30/2019 - modified - Brennan Young                           //
// - objects sorted in quickSort now only need the < logical       //
//   operator (instead of <=).                                     //
// 04/22/2020 - modified - Brennan Young                           //
// - added macro guard.
/////////////////////////////////////////////////////////////////////

#ifndef BYSTDLIB_SORT_20200422
#define BYSTDLIB_SORT_20200422

#include <vector>

namespace bystd { // Brennan Young standard namespace

// Swap the values of two variables.
// -- Arguments --
// a : a variable.
// b : a variable.
// -- Returns --
// Nothing. Swaps the values of the given variables.
template <class T>
void swap ( T * a, T * b )
{
    T c = *a;
    *a = *b;
    *b = c;
}

// Quick sort.
// -- Arguments --
// A    : input array.
// low  : first index to be sorted.
// high : last index to be sorted.
// -- Returns --
// Nothing. Sorts the given array.
template <typename T>
int quickSort_partition ( std::vector<T> * A, int low, int high )
{
    T pivot = (*A)[high];
    int i = -1;
    for ( int j = 0; j < high; ++j ) {
        if ( !(pivot < (*A)[j]) )
            swap(&(*A)[++i], &(*A)[j]);
    }
    swap(&(*A)[i+1], &(*A)[high]);
    return i+1;
}

template <typename T>
void quickSort ( std::vector<T> * A, int low, int high )
{
    if ( low < high ) {
        int i = quickSort_partition(A, low, high);
        quickSort(A, low, i-1);
        quickSort(A, i+1, high);
    }
}

template <typename T>
void quickSort ( std::vector<T> * A )
{
    quickSort(A, 0, A->size()-1);
}

template <typename T>
int quickSort_partition ( T * A, int low, int high )
{
    T pivot = A[high];
    int i = -1;
    for ( int j = 0; j < high; ++j ) {
        if ( !(pivot < A[j]) )
            swap(&A[++i], &A[j]);
    }
    swap(&A[i+1], &A[high]);
    return i+1;
}

template <typename T>
void quickSort ( T * A, int low, int high )
{
    if ( low < high ) {
        int i = quickSort_partition(A, low, high);
        quickSort(A, low, i-1);
        quickSort(A, i+1, high);
    }
}

template <typename T>
void quickSort ( T * A, int n )
{
    quickSort(A, 0, n-1);
}

// Bubble sort
template <typename T>
void bubbleSort ( std::vector<T> * P )
{
    int i, j, k;
    for ( i = 0; i < P->size() - 1; ++i ) {
        for ( j = 0; j < P->size()-i-1; ++j )
            if ( (*P)[j+1] < (*P)[j] ) swap(&(*P)[j+1], &(*P)[j]);
    }
}

template <typename T>
void bubbleSort ( T * P, int n )
{
    int i, j, k;
    for ( i = 0; i < n - 1; ++i ) {
        for ( j = 0; j < n-i-1; ++j )
            if ( P[j+1] < P[j] ) swap(&P[j+1], &P[j]);
    }
}

} // namespace bystd

#endif // BYSTDLIB_SORT_20200422