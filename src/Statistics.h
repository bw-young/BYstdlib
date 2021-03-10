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
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_STATISTICS_20200922
#define YOUNG_STATISTICS_20200922

#include <cstddef>
#include <cmath>

class Statistics {
public:
    size_t n;       // number of not-ignored elements
    double min;     // minimum value
    double max;     // maximum value
    double sum;     // sum
    double sum2;    // sum of squares
    double mean;    // mean value
    double var;     // variance
    double stdev;   // standard deviation
    
    // constructors, destructors
    Statistics();
    template<typename T, class Iterator>
        Statistics(const Iterator&, const Iterator&, const T&);
    template<typename T, class Container>
        Statistics(const Container&, const size_t, const size_t,
        const size_t, const T&);
    Statistics(const Statistics&);
    ~Statistics();
    
    // operators
    Statistics& operator=(const Statistics&);
};


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

// Default constructor.
Statistics::Statistics () : n(0), min(0.0), max(0.0), sum(0.0),
    sum2(0.0), mean(0.0), stdev(0.0) {}

// Construct from container, from begin to and not including end,
// ignoring all values with the indicated 'ignore' value.
template <typename T, class Iterator>
Statistics::Statistics ( const Iterator& begin, const Iterator& end,
    const T& ignore ) : n(0), min(0.0), max(0.0), sum(0.0),
    sum2(0.0), mean(0.0), stdev(0.0)
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
// 'ignore' value, skipping every step elements.
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
Statistics::Statistics ( const Statistics& stat ) : n(stat.n),
    min(stat.min), max(stat.max), sum(stat.sum), sum2(stat.sum2),
    mean(stat.mean), stdev(stat.stdev)
{}

// Destructor.
Statistics::~Statistics () {}


// OPERATORS ////////////////////////////////////////////////////////

Statistics& Statistics::operator= ( const Statistics& stat )
{
    if ( this == &stat) return *this;
    
    n = stat.n;
    min = stat.min;
    max = stat.max;
    sum = stat.sum;
    sum2 = stat.sum2;
    mean = stat.mean;
    stdev = stat.stdev;
    
    return *this;
}

#endif // YOUNG_STATISTICS_20200922