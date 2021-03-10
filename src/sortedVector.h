// created 09/20/2017 by Brennan Young
// modified through 10/12/2017 by Brennan Young
// modified 01/11/2018 by Brennan Young
// - included stddef.h library for size_t

// Functions for treating a std::vector as a sorted vector. Only uses < to
// compare elements

#ifndef YOUNG_STDLIB_SORTEDVECTOR_20170920
#define YOUNG_STDLIB_SORTEDVECTOR_20170920

#include <stddef.h> // size_t
#include <vector>   // std::vector

namespace bystd { // Brennan Young standard namespace

// Returns the index of an element, or the index where the element would be
// inserted.
template <class T> size_t sortedVector_index(const std::vector<T> & v, T x){
    if ( v.size() == 0 ){ return v.size(); }
    
    size_t i = v.size() / 2; // index being checked
    size_t ii;
    size_t li = 0;             // lower bound
    size_t ui = v.size();      // upper bound
    
    while ( v[i] < x || x < v[i] ){
        if ( v[i] < x ){
            if ( i == v.size() - 1 ){ return v.size(); }
            li = i + 1;
        }
        if ( x < v[i] ) ui = i;
        ii = li + (ui - li) / 2;
        if ( i == ii ) return i;
        i = ii;
    }
    
    return i;
}

// Returns the index of an element, or std::vector::size() if the element is
// not in the vector.
template <class T> size_t sortedVector_find(const std::vector<T> & v, T x){
    size_t i = sortedVector_index(v, x);
    return ( v[i] < x || x < x[i] ) ? v.size() : i;
}

// Returns the index of the inserted element.
template <class T> size_t sortedVector_insert(std::vector<T> & v, T x){
    size_t i = sortedVector_index(v, x);
    v.insert(v.begin() + i, x);
    return i;
}

// Erases the element from the vector, if it exists.
template <class T> void sortedVector_remove(std::vector<T> & v, T x){
    size_t i = sortedVector_find(v, x);
    if ( i == v.size() ) return;
    v.erase(v.begin() + i, x);
}

} // namespace bystd

#endif // YOUNG_STDLIB_SORTEDVECTOR_20170920