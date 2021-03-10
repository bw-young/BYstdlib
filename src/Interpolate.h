/////////////////////////////////////////////////////////////////////
// Functions for interpolating datasets.                           //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// ??/??/2019 - Created - Brennan Young                            //
// 05/13/2020 - Modified - Brennan Young                           //
// - Added bilinearInterpolate and barycentricInterpolate.         //
/////////////////////////////////////////////////////////////////////

#include <cmath>
#include "Matrix.h"


// 1D Interpolation /////////////////////////////////////////////////


// Linear interpolation.
// -- Source --
// "Interpolation Methods" by Paul Bourke, Dec 1999.
// -- Arguments --
// y1 : value at point 1.
// y2 : value at point 2.
// mu : where to interpolate [0=position of point 1, 1=position of
//      point 2].
// -- Returns --
// Value interpolated at mu.
double linearInterpolate ( double y1, double y2, double mu )
{
    return y1 * (1.0 - mu) + y2 * mu;
}

// Cosine interpolation.
// -- Source --
// "Interpolation Methods" by Paul Bourke, Dec 1999.
// -- Arguments --
// y1 : value at point 1.
// y2 : value at point 2.
// mu : where to interpolate [0=position of point 1, 1=position of
//      point 2].
// -- Returns --
// Value interpolated at mu.
double cosineInterpolate ( double y1, double y2, double mu )
{
    static const double pi = 4.0 * atan(1.0);
    double mu2 = 0.5 * (1.0 - cos(mu * pi));
    return y1 * (1.0 - mu2) + y2 * mu2;
}

// Cubic spline interpolation. Does not address issue for how to
// interpolate between first and last segments; a common solution is
// to invent extra points at the start and end such that they have
// the same slope as the start or end segment.
// -- Source --
// "Interpolation Methods" by Paul Bourke, Dec 1999.
// -- Arguments --
// y0 : value at point before point 1.
// y1 : value at point 1.
// y2 : value at point 2.
// y3 : value at point after point 2.
// mu : where to interpolate [0=position of point 1, 1=position of
//      point 2].
// -- Returns --
// Value interpolated at mu.
double cubicInterpolate ( double y0, double y1, double y2,
    double y3, double mu )
{
    double mu2 = mu * mu;
    double a0 = y3 -y2 - y0 + y1;
    double a1 = y0 - y1 - a0;
    double a2 = y2 - y0;
    double a3 = y1;
    return (a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3);
}

// Catmull-Rom spline. Does not address issue for how to interpolate
// between first and last segments; see notes for cubicInterpolate.
// -- Source --
// "Interpolation Methods" by Paul Bourke, Dec 1999.
double catmullRomInterpolate ( double y0, double y1, double y2,
    double y3, double mu )
{
    double mu2 = mu * mu;
    double a0 = -0.5 * y0 + 1.5 * y1 - 1.5 * y2 + 0.5 * y3;
    double a1 = y0 - 2.5 * y1 + 2.0 * y2 - 0.5 * y3;
    double a2 = -0.5 * y0 + 0.5 * y2;
    double a3 = y1;
    return (a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3);
}

// Cubic hermite spline. Does not address issue for how to
// interpolate between first and last segments; see notes for
// cubicInterpolate.
// -- Source --
// "Interpolation Methods" by Paul Bourke, Dec 1999.
// -- Arguments --
// y0      : value at point before point 1.
// y1      : value at point 1.
// y2      : value at point 2.
// y3      : value at point after point 2.
// mu      : where to interpolate [0=position of point 1, 1=position
//           of point 2].
// tension : tightness of curvature at known points [1=high,
//           0=normal, -1=low].
// bias    : twistiness of curve about known points [0=even,
//           +=towards first segment, -=towards other segment].
// -- Returns --
// Value interpolated at mu.
double hermiteInterpolate ( double y0, double y1, double y2,
    double y3, double mu, double tension, double bias )
{
    double mu2 = mu * mu;
    double mu3 = mu2 * mu;
    double t = 0.5 * (1.0 - tension);
    double m0 = (y1 - y0) * (1.0 + bias) * t
              + (y2 - y1) * (1.0 - bias) * t;
    double m1 = (y2 - y1) * (1.0 + bias) * t
              + (y3 - y2) * (1.0 - bias) * t;
    double a0 =  2.0 * mu3 - 3.0 * mu2 + 1.0;
    double a1 =        mu3 - 2.0 * mu2 + mu;
    double a2 =        mu3 -       mu2;
    double a3 = -2.0 * mu3 + 3.0 * mu2;
    return (a0 * y1 + a1 * m0 + a2 * m1 + a3 * y2);
}


// 2D Interpolation /////////////////////////////////////////////////


// Interpolate using barycentric coordinates of a triangle.
// -- Source --
// https://codeplea.com/triangular-interpolation
// -- Arguments --
// x, y : x and y coordinates to be interpolated.
// X, Y : X- and Y-coordinate vectors, n=3.
// Z : the value vector at (X, Y), n=3.
double barycentricInterpolate ( double x, double y, double * X, double * Y,
    double * Z )
{
    double denom = (Y[1] - Y[2])*(X[0] - X[2]) + (X[2] - X[1])*(Y[0] - Y[2]);
    double W[3];
    W[0] = ((Y[1] - Y[2]) * (x - X[2]) + (X[2] - X[1]) * (y - Y[2])) / denom;
    W[1] = ((Y[2] - Y[0]) * (x - X[2]) + (X[0] - X[2]) * (y - Y[2])) / denom;
    W[2] = 1.0 - W[0] - W[1];
    return Z[0] * W[0] + Z[1] * W[1] + Z[2] * W[2];
}

// Interpolate using bilinear interpolation of points in an irregular grid.
// -- Source --
// http://www.ahinson.com/algorithms_general/Sections/InterpolationRegression/InterpolationIrregularBilinear.pdf
// -- Arguments --
// x, y : x and y coordinates to be interpolated.
// X, Y : X- and Y-coordinate vectors, n=4.
// Z : the value vector at (X, Y), n=4.
double bilinearInterpolate ( double x, double y, double * X, double * Y,
    double * Z )
{
    typedef struct point {
        double x;
        double y;
        double z;
        point(double xx=0.0, double yy=0.0, double zz=0.0) :
            x(xx), y(yy), z(zz) {}
    } Point;
    
    int     i, j, k;
    Point   P[4]; // points at UL, UR, BL, BR
    Point   A, B, C, D;
    Point   tmp;
    bool    vParallel, hParallel;
    double  a, b, c, d, s, t;
    
    P[0] = Point(X[0], Y[0], Z[0]);
    P[1] = Point(X[1], Y[1], Z[1]);
    P[2] = Point(X[2], Y[2], Z[2]);
    P[3] = Point(X[3], Y[3], Z[3]);
    
    // sort by y to get top points in P[0] and P[1], bottom in P[2] and P[3]
    for ( i = 0; i < 4; ++i ) {
        k = i;
        for ( j = i + 1; j < 4; ++j ) {
            if ( P[k].y < P[j].y ) k = j;
        }
        if ( k != i ) {
            tmp = P[i];
            P[i] = P[k];
            P[k] = tmp;
        }
    }
    
    // sort left and right points
    if ( P[1].x < P[0].x ) {
        tmp = P[0];
        P[0] = P[1];
        P[1] = tmp;
    }
    if ( P[3].x < P[2].x ) {
        tmp = P[2];
        P[2] = P[3];
        P[3] = tmp;
    }
    
    // determine edge parallelism using the cross product of edges
    vParallel = (P[2].x - P[0].x) * (P[3].y - P[1].y)
        - (P[3].x - P[1].x) * (P[2].y - P[0].y) == 0.0;
    hParallel = (P[1].x - P[0].x) * (P[3].y - P[2].y)
        - (P[3].x - P[2].x) * (P[1].y - P[0].y) == 0.0;
    
    if ( !vParallel ) {         // method 1
        std::cout << "method1\n";
        double x21 = P[1].x - P[0].x;
        double x31 = P[2].x - P[0].x;
        double x42 = P[3].x - P[1].x;
        double y21 = P[1].y - P[0].y;
        double y31 = P[2].y - P[0].y;
        double y42 = P[3].y - P[1].y;
        
        a = x31 * y42 - y31 * x42;
        b = y * (x42 - x31) - x * (y42 - y31)
            + x31 * P[1].y - y31 * P[1].x + P[0].x * y42 - P[0].y * x42;
        c = y * x21 - x * y21 + P[0].x * P[1].y - P[1].x * P[0].y;
        
        // solve for t by quadratic equation
        d = sqrt(b*b - 4.0 * a * c);
        std::cout << "a=" << a << " b=" << b << " c=" << c << " d=" << d << "\n";
        t = (-b + d) / (2.0 * a);
        if ( t < 0.0 || t > 1.0 ) t = (-b - d) / (2.0 * a);
        
        // solve for s
        s = (y - P[0].y - y31 * t)
            / (P[1].y + y42 * t - P[0].y - y31 * t);
        
        std::cout << "t=" << t << " s=" << s
                  << " ?=" << (P[1].y + y42 * t - P[0].y - y31 * t)
                  << " <- " << P[1].y
                  << " + " << (y42 * t)
                  << " - " << P[0].y 
                  << " - " << (y31 * t)
                  << "\n";
    }
    else if ( !hParallel ) {    // method 2
        std::cout << "method 2\n";
        double x21 = P[1].x - P[0].x;
        double x31 = P[2].x - P[0].x;
        double x43 = P[3].x - P[2].x;
        double y21 = P[1].y - P[0].y;
        double y31 = P[2].y - P[0].y;
        double y43 = P[3].y - P[2].y;
        
        a = x21 * y43 - y21 * x43;
        b = y * (x43 - x21) - x * (y43 - y21) + P[0].x * y43
            - P[0].y * x43 + x21 * P[2].y - y21 * P[2].x;
        c = y * x31 - x * y31 + P[0].x * P[2].y - P[2].x * P[0].y;
        
        // solve for s by quadratic equation
        d = sqrt(b*b - 4.0 * a * c);
        s = (-b + d) / (2.0 * a);
        if ( s < 0.0 || s > 1.0 ) s = (-b - d) / (2.0 * a);
        
        // solve for t
        t = (y - P[0].y - y21 * s)
            / P[2].y + y43 * s - P[0].y - y21 * s;
    }
    else {                      // method 3
        double x21 = P[1].x - P[0].x;
        double x31 = P[2].x - P[0].x;
        double y21 = P[1].y - P[0].y;
        double y31 = P[2].y - P[0].y;
        
        Matrix<double> mA (2, 2);
        mA[0] = x21;    mA[1] = x31;
        mA[2] = y21;    mA[3] = y31;
        
        Matrix<double> mB (2, 1); // s, t
        
        Matrix<double> mC (2, 1);
        mC[0] = x - P[0].x;
        mC[1] = y - P[0].y;
        
        Matrix<double> mAi = mA.inverse();
        mB = mA.inverse() * mC;
        s = mB[0];
        t = mB[1];
    }
    
    // solve interpolation problem
    return P[0].z * (1.0 - s) * (1.0 - t)
         + P[1].z * s * (1.0 - t)
         + P[2].z * (1.0 - s) * t
         + P[3].z * s * t;
}