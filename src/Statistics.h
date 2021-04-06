/////////////////////////////////////////////////////////////////////
// Statistics                                                      //
// Object for computing and storing basic 1D statistics.           //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- History ---------------------------------------------------- //
// 09/22/2020 - Brennan Young                                      //
// - created.                                                      //
// 12/09/2020 - Brennan Young                                      //
// - added macro guard.                                            //
// 04/02/2021 - Brennan Young                                      //
// - added to bystd namespace.                                     //
// - added median, mode members.                                   //
// - added 'middle' method to compute the median and mode.         //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_STATISTICS_20200922
#define YOUNG_STATISTICS_20200922

#include <cstddef>
#include <cmath>
#include <map>


namespace bystd {


class Statistics {
public:
    size_t n;       // count of not-ignored elements
    double min;     // minimum value
    double max;     // maximum value
    double sum;     // sum
    double sum2;    // sum of squares
    double mean;    // mean value
    double var;     // variance
    double stdev;   // standard deviation
    double median;  // middle value
    double mode;    // most common value
    
    // constructors, destructors
    Statistics();
    template<typename T, class Iterator> Statistics(
        const Iterator&, const Iterator&, const T&);
    template<typename T, class Container> Statistics(
        const Container&, const size_t, const size_t, const size_t,
        const T&);
    Statistics(const Statistics&);
    ~Statistics();
    
    // operators
    Statistics& operator=(const Statistics&);
    
    // operations
    template<typename T, class Iterator> Statistics& middle(
        const Iterator&, const Iterator&, const T&);
    template<typename T, class Container> Statistics& middle(
        const Container&, const size_t, const size_t, const size_t,
        const T&);
};


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

// Default constructor.
Statistics::Statistics () : n(0), min(0.0), max(0.0), sum(0.0),
    sum2(0.0), mean(0.0), stdev(0.0), median(0.0), mode(0.0) {}

// Construct from container, from begin to and not including end,
// ignoring all values with the indicated 'ignore' value. Does not
// compute median and mode.
template <typename T, class Iterator>
Statistics::Statistics ( const Iterator& begin, const Iterator& end,
    const T& ignore )
    : n(0), min(0.0), max(0.0), sum(0.0), sum2(0.0), mean(0.0),
    stdev(0.0), median(0.0), mode(0.0)
{
    Iterator it = begin;
    while ( it != end && *it == ignore ) ++it;
    if ( it == end ) return;
    
    min = max = *it;
    
    for ( ; it != end; ++it ) {
        if ( *it == ignore ) continue;
        ++n;
        if ( *it < min ) min = *it;
        if ( max < *it ) max = *it;
        sum += *it;
        sum2 += ((*it) * (*it));
    }
    
    mean = sum / n;
    var = sum2 / n - mean * mean;
    stdev = sqrt(var);
}

// Construct from random access container, from index being to and
// not including index end, ignoring all values with the indicated
// 'ignore' value, skipping every step elements. Does not compute
// median and mode.
template <typename T, class Container>
Statistics::Statistics ( const Container& A, const size_t begin,
    const size_t end, const size_t step, const T& ignore )
    : n(0), min(0.0), max(0.0), sum(0.0), sum2(0.0), mean(0.0),
    stdev(0.0)
{
    size_t i = begin;
    while ( i != end && A[i] == ignore ) i += step;
    if ( i == end ) return;
    
    min = max = A[i];
    
    for ( ; i != end; i += step ) {
        if ( A[i] == ignore ) continue;
        ++n;
        if ( A[i] < min ) min = A[i];
        if ( max < A[i] ) max = A[i];
        sum += A[i];
        sum2 += (A[i] * A[i]);
    }
    
    mean = sum / n;
    var = sum2 / n - mean * mean;
    stdev = sqrt(var);
}

// Copy constructor.
Statistics::Statistics ( const Statistics& stat )
    : n(stat.n), min(stat.min), max(stat.max), sum(stat.sum),
    sum2(stat.sum2), mean(stat.mean), stdev(stat.stdev),
    median(stat.median), mode(stat.mode)
{}

// Destructor.
Statistics::~Statistics () {}


// OPERATORS ////////////////////////////////////////////////////////

// Assignment to copy another Statistics object.
Statistics& Statistics::operator= ( const Statistics& stat )
{
    if ( this == &stat) return *this;
    
    n      = stat.n;
    min    = stat.min;
    max    = stat.max;
    sum    = stat.sum;
    sum2   = stat.sum2;
    mean   = stat.mean;
    stdev  = stat.stdev;
    median = stat.median;
    mode   = stat.mode;
    
    return *this;
}


// OPERATIONS ///////////////////////////////////////////////////////

// Compute the median and mode in a container, from begin to and not
// including end, ignoring all values with the indicated 'ignore'
// value.
template <typename T, class Iterator>
Statistics& Statistics::middle ( const Iterator& begin,
    const Iterator& end, const T& ignore )
{
    Iterator it = begin;
    while ( it != end && *it == ignore ) ++it;
    
    std::map<T, int> values;
    T valMaxCount = 0;
    values[valMaxCount] = 0;
    
    for ( ; it != end; ++it ) {
        if ( *it == ignore ) continue;
        
        if ( values.find(*it) == values.end() ) values[*it] = 1;
        else values[*it] += 1;
        
        if ( values[*it] > values[valMaxCount] ) valMaxCount = *it;
    }
    
    mode = valMaxCount;
    
    typename std::map<T, int>::iterator item = values.begin();
    int n = values.size() / 2;
    for ( ; n > 0 && item != values.end(); --n ) ++item;
    median = (item == values.end() ? 0.0 : item->first);
    
    return *this;
}

// Construct from random access container, from index being to and
// not including index end, ignoring all values with the indicated
// 'ignore' value, skipping every step elements. Does not compute
// median and mode.
template <typename T, class Container>
Statistics& Statistics::middle ( const Container& A,
    const size_t begin, const size_t end, const size_t step,
    const T& ignore )
{
    size_t i = begin;
    while ( i < A.size() && A[i] == ignore ) i += step;
    
    std::map<T, int> values;
    T valMaxCount = 0;
    values[valMaxCount] = 0;
    
    for ( ; i < A.size(); i += step ) {
        if ( A[i] == ignore ) continue;
        
        if ( values.find(A[i]) == values.end() ) values[A[i]] = 1;
        else values[A[i]] += 1;
        
        if ( values[A[i]] > values[valMaxCount] )
            valMaxCount = A[i];
    }
    
    mode = valMaxCount;
    
    typename std::map<T, int>::iterator item = values.begin();
    int n = values.size() / 2;
    for ( ; n > 0 && item != values.end(); --n ) ++item;
    median = (item == values.end() ? 0.0 : item->first);
    
    return *this;
}


}; // namespace bystd

#endif // YOUNG_STATISTICS_20200922