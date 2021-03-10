/////////////////////////////////////////////////////////////////////
// Integration functions.                                          //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// 11/21/2019 - Created - Brennan Young                            //
// 05/13/2020 - Modified - Brennan Young                           //
// - moved under bystd namespace and BYstdlib.                     //
// - added macro guard.                                            //
/////////////////////////////////////////////////////////////////////

#ifndef BYOUNG_BYSTDLIB_20200513
#define BYOUNG_BYSTDLIB_20200513

namespace bystd { // Brennan Young standard namespace

// Integrate the function y(x) using the trapezoidal rule.
template<typename T>
double integrateTrapezoid ( std::vector<T>& x, std::vector<T>& y )
{
    double sum = 0.0;
    for ( size_t i = 1; i < x.size(); ++i )
        sum += 0.5 * (y[i-1] + y[i]) * (x[i] - x[i-1]);
    return sum;
}

template<typename T>
double integrateTrapezoid ( T* x, T* y, int n )
{
    double sum = 0.0;
    for ( int i = 1; i < n; ++i )
        sum += 0.5 * (y[i-1] + y[i]) * (x[i] - x[i-1]);
    return sum;
}

} // bystd

#endif // BYOUNG_BYSTDLIB_20200513