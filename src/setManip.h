/////////////////////////////////////////////////////////////////////
// Functions for manipulating sets.                                //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/10/2020 - Created - Brennan Young                            //
// 04/14/2020 - Modified - Brennan Young                           //
// - added setDifference.                                          //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_STDLIB_SETMANIP
#define YOUNG_STDLIB_SETMANIP

#include <set>

namespace bystd {

//Return a set of items found in A but not in B.
template<typename T>
std::set<T> setDifference ( const std::set<T>& A,
    const std::set<T>& B )
{
    std::set<T> C;
    typename std::set<T>::iterator Ai = A.begin();
    typename std::set<T>::iterator Bi = B.begin();
    for ( ; Ai != A.end() && Bi != B.end(); ) {
        if ( *Ai < *Bi ) {
            C.insert(*Ai);
            ++Ai;
        }
        else if ( *Bi < *Ai ) ++Bi;
        else {
            ++Ai;
            ++Bi;
        }
    }
    C.insert(Ai, A.end());
    return C;
}

// Return a set of all items that are in A and B.
template<typename T>
std::set<T> setIntersect ( const std::set<T>& A,
    const std::set<T>& B )
{
    std::set<T> C;
    typename std::set<T>::iterator Ai = A.begin();
    typename std::set<T>::iterator Bi = B.begin();
    for ( ; Ai != A.end() && Bi != B.end(); ) {
        if ( *Ai < *Bi ) ++Ai;
        else if ( *Bi < *Ai ) ++Bi;
        else {
            C.insert(*Ai);
            ++Ai;
            ++Bi;
        }
    }
    return C;
}

// Return a set of all items that are in A or B.
template<typename T>
std::set<T> setUnion ( const std::set<T>& A, const std::set<T>& B )
{
    std::set<T> C = A;
    C.insert(B.begin(), B.end());
    return C;
}

}; // namespace bystd

#endif // YOUNG_STDLIB_SETMANIP