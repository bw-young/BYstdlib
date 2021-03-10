/////////////////////////////////////////////////////////////////////
// Functions for regression equations.                             //
/////////////////////////////////////////////////////////////////////

#ifndef REGRESSION_H
#define REGRESSION_H

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>


// LINEAR ///////////////////////////////////////////////////////////


// -- Source --
// created by Bishop, M. P., Geography, Texas A&M University
// received 03/2018 by Brennan Young
// modified 03/30/2018 by Brennan Young
// - remove approximate probability.
// - add constructor.
// modified 06/11/2018 by Brennan Young
// - overloaded PearsonsCorrelation for float pointers.
// - corrected error in PearsonsCorrelation where sumxy was not initialized to
//   0.
struct Regression {
    double mean;
    double df;               // Degrees of freedom
    double r;                // Pearsons product moment correlation coefficient
    double r2;               // Coefficient of determination
    double t;                // Student's t statistic
    double z;                // Fisher's z transformation
    double prob;            // Student's t probability
    double slope;            // Linear regression coefficient
    double intercept;        // Y-intercept of regression line
    
    Regression() : mean(0.0), df(0.0), r(0.0), r2(0.0), t(0.0), z(0.0),
        prob(0.0), slope(0.0), intercept(0.0) {}
};

// -- Source --
// created by Bishop, M. P., Geography, Texas A&M University
// modified 03/30/2018 by Brennan Young
// - some reorganization and optimization.
double gammln (double xx){
    register int i;
    
    double x = xx;
    double y = xx;
    static double cof[6] = { 76.18009172947146, -86.50532032941677,
                             24.01409824083091,  -1.231739572450155,
                             0.1208650973866179e-2, -0.5395239384953e-5 };
    double tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    double ser = 1.000000000190015;
    
    for ( i = 0; i < 6; ++i ) ser += cof[i] / ++y;
    
    return -tmp + log(2.5066282746310005 * ser / x);
}

// -- Source --
// created by Bishop, M. P., Geography, Texas A&M University
double betacf (double a, double b, double x){
    const int MAXIT = 100;
    const double EPS = 3.0e-7;
    const double FPMIN = 1.0e-30;
    int m, m2;
    double aa, c, d, del, h, qab, qam, qap;
    
    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;
    c = 1.0;
    d = 1.0 - qab * x / qap;
    if ( fabs(d) < FPMIN ) d = FPMIN;
    d = 1.0 / d;
    h = d;
    
    for ( m = 1; m <= MAXIT; ++m ) {
        m2 = 2.0 * m;
        aa = m*(b-m)*x / ((qam+m2)*(a+m2));
        d = 1.0+aa*d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 +aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN; //c + FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del-1.0) < EPS) break;
    }
    if ( m > MAXIT ) std::cout << "ERR: a or b too big, or MAXIT too small in betacf\n";
    return h;
}

// -- Source --
// created by Bishop, M. P., Geography, Texas A&M University
// modified 03/30/2018 by Brennan Young
// - some minor rearranging.
double betai (double a, double b, double x){
    double bt;
    
    if ( x < 0.0 || x > 1.0 ) {
        std::cout << "ERR: Bad x in function betai\n";
        exit (1);
    }
    if ( x == 0.0 || x == 1.0 ) bt = 0.0;
    else bt = exp(gammln(a+b) - gammln(a) - gammln(b)
        + a * log(x) + b * log(1.0-x));
    
    if ( x < (a+1.0) / (a+b+2.0) ) return bt * betacf (a, b, x) / a;
    else return 1.0 - bt * betacf (b, a, 1.0-x) / b;
}

// -- Source --
// created by Bishop, M. P., Geography, Texas A&M University
// modified 03/30/2018 by Brennan Young
// - returns a Regression object instead of taking a pointer to one as an
//   argument.
// - overloaded to accept x and y vectors instead of pointer arrays.
Regression PearsonsCorrelation (double * x, double * y, int n){
    register int i;
    
    const double TINY = 1.0e-20;
    double sumx, sumy, sumx2, sumy2, sumxy;
    Regression r;
    
    r.df = n - 2.0;
    sumx = sumy = sumx2 = sumy2 = sumxy = 0.0;
    
    for ( i = 0; i < n; ++i ) {
        sumx += x[i];
        sumy += y[i];
        sumx2 += x[i]*x[i];
        sumy2 += y[i]*y[i];
        sumxy += x[i]*y[i];
    }
    
    r.r = (n*sumxy - sumx*sumy)
        / sqrt((n*sumx2 - sumx*sumx)*(n*sumy2 - sumy*sumy));
    r.r2 = r.r * r.r;
    if ( r.r < -1.0 || r.r > 1.0 ) {
        r.r = 0.0;
        r.r2 = 0.0;
    }
    else {
        r.t = r.r * sqrt(r.df / (1.0 - r.r2));
        r.z = 0.5 * log ((1.0 + r.r + TINY) / (1.0 - r.r + TINY));
        r.prob = betai ((0.5 * r.df), 0.5, (r.df / (r.df + r.t * r.t)));
        r.slope = (n * sumxy - sumx * sumy) / (n * sumx2 - sumx * sumx);
        r.intercept = sumy / n - (r.slope * sumx / n);
    }
    
    return r;
}

Regression PearsonsCorrelation (float * x, float * y, int n){
    register int i;
    
    const double TINY = 1.0e-20;
    double sumx, sumy, sumx2, sumy2, sumxy;
    Regression r;
    
    r.df = n - 2.0;
    sumx = sumy = sumx2 = sumy2 = sumxy = 0.0;
    
    for ( i = 0; i < n; ++i ) {
        sumx += x[i];
        sumy += y[i];
        sumx2 += x[i]*x[i];
        sumy2 += y[i]*y[i];
        sumxy += x[i]*y[i];
    }
    
    r.r = (n*sumxy - sumx*sumy)
        / sqrt((n*sumx2 - sumx*sumx)*(n*sumy2 - sumy*sumy));
    r.r2 = r.r * r.r;
    if ( r.r < -1.0 || r.r > 1.0 ) {
        r.r = 0.0;
        r.r2 = 0.0;
    }
    else {
        r.t = r.r * sqrt(r.df / (1.0 - r.r2));
        r.z = 0.5 * log ((1.0 + r.r + TINY) / (1.0 - r.r + TINY));
        r.prob = betai ((0.5 * r.df), 0.5, (r.df / (r.df + r.t * r.t)));
        r.slope = (n * sumxy - sumx * sumy) / (n * sumx2 - sumx * sumx);
        r.intercept = sumy / n - (r.slope * sumx / n);
    }
    
    return r;
}

template <class T>
Regression PearsonsCorrelation (const std::vector<T> & x,
        const std::vector<T> & y) {
    
    register int i;
    
    int n = x.size() < y.size() ? x.size() : y.size();
    const double TINY = 1.0e-20;
    double sumx, sumy, sumx2, sumy2, sumxy;
    Regression r;
    
    r.df = n - 2.0;
    sumx = sumy = sumx2 = sumy2 = sumxy = 0.0;
    
    for ( i = 0; i < n; ++i ) {
        sumx += x[i];
        sumy += y[i];
        sumx2 += x[i] * x[i];
        sumy2 += y[i] * y[i];
        sumxy += x[i] * y[i];
    }
    r.mean = sumy / n;
    
    r.r = (n*sumxy - sumx*sumy)
        / sqrt((n*sumx2 - sumx*sumx) * (n*sumy2 - sumy*sumy));
    r.r2 = r.r * r.r;
    
    if ( r.r < -1.0 || r.r > 1.0 ) {
        r.r = 0.0;
        r.r2 = 0.0;
    }
    else {
        r.t = r.r * sqrt(r.df / (1.0 - r.r2));
        r.z = 0.5 * log ((1.0 + r.r + TINY) / (1.0 - r.r + TINY));
        r.prob = betai ((0.5 * r.df), 0.5, (r.df / (r.df + r.t * r.t)));
        r.slope = (n * sumxy - sumx * sumy) / (n * sumx2 - sumx * sumx);
        r.intercept = sumy / n - (r.slope * sumx / n);
    }
    
    return r;
}


// NON-LINEAR ///////////////////////////////////////////////////////


// Power-law regression for the equation y = a * x^b.
// -- Source --
// Created 05/08/2020 by Brennan W. Young
// -- Arguments --
// x : independent variable array.
// y : dependent variable array.
// n : length of x and y arrays.
// a : (output) coefficient.
// b : (output) power coefficient.
template<typename T, typename U>
void powerLawRegression ( T* x, T* y, int n, U* a, U* b )
{
    double lnx, lny;
    double sumx = 0.0;
    double sumy = 0.0;
    double sumx2 = 0.0;
    double sumxy = 0.0;
    
    for ( int i = 0; i < n; ++i ) {
        lnx = log(x[i]);
        lny = log(y[i]);
        
        sumx += lnx;
        sumy += lny;
        sumx2 += lnx * lnx;
        sumxy += lnx * lny;
    }
    
    *b = (n * sumxy - sumx * sumy) / (n * sumx2 - sumx * sumx);
    *a = exp((sumy - *b * sumx) / n);
}

#endif // REGRESSION_H