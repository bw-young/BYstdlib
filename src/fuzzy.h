/////////////////////////////////////////////////////////////////////
// Fuzzy transformation functions.                                 //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/11/2020 - Brennan Young                                      //
// - modified.                                                     //
// 03/26/2021 - Brennan Young                                      //
// - added linearMembership.                                       //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_STDLIB_FUZZY_H
#define YOUNG_STDLIB_FUZZY_H

#include <cmath>

// Linear fuzzy membership.
// -- Arguments --
// x : value to transform to membership.
// a : x-value for membership = 0.
// b : x-value for membership = 1.
double linearMembership ( double x, double a, double b )
{
    if ( a < b ) {
        if ( x < a ) return 0;
        if ( x > b ) return 1;
    }
    else {
        if ( x < b ) return 0;
        if ( x > a ) return 1;
    }
    return (x - a) / (b - a);
}

// Triangular fuzzy function.
// -- Arguments --
// x : value to transform to membership.
// a : left x-value where membership is 0.
// b : middle x-vaue where membership is 1.
// c : right x-value where membership is 0.
double triangularMembership ( double x, double a, double b,
    double c )
{
    if ( x <= a ) return 0.0;
    else if ( x < b ) return (x-a) / (b-a);
    else if ( x < c ) return (c-x) / (c-b);
    else return 0.0;
}

// Power fuzzy function.

// Log fuzzy function.

// Binary sigmoid fuzzy function, or logistic function.
// -- Arguments --
// x  : value to transform to membership.
// x0 : x-value at sigmoid's midpoint.
// L  : curve's maximum value.
// k  : logistic growth rate (curve steepness).
// -- Notes --
// Standard logistic function = logisticSigmoid(x,0,1,1);
// As used in neural networks = logisticSigmoid(x,x0,1,2*beta);
double logisticSigmoid ( double x, double x0, double L, double k )
{
    return L / (1.0 + exp(-k * (x - x0)));
}

#endif // YOUNG_STDLIB_FUZZY_H