/////////////////////////////////////////////////////////////////////
// General-use structures and functions for manipulating geometric //
// primitives.                                                     //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// created 09/14/2017 by Brennan Young                             //
// modified 12/06/2017 by Brennan Young                            //
// - added avgAzimuth.                                             //
// modified 12/07/2017 by Brennan Young                            //
// - added lineSlope and lineYint to compute the equation of a     //
//   line.                                                         //
// - added lineSlopeOrth to compute the perpendicular to a line.   //
// - added intercept to compute the intersection of two lines.     //
// - overloaded dist to compute the distance between a point and a //
//   line.                                                         //
// - added circularDiff.                                           //
// modified 12/08/2017 by Brennan Young                            //
// - azimuth and avgAzimuth now correctly detects 'no solution'.   //
// - circularDiff now returns a value.                             //
// modified 01/11/2018 by Brennan Young                            //
// - added inLine and inPolygon functions.                         //
// - added Jarvis, chain, and Chan convex hull functions, with     //
//   helper functions less, equal, and sortedXY_index.             //
// - created table of contents                                     //
// modified 01/24/2018 by Brennan Young                            //
// - created visualCenter.                                         //
// modified 04/01/2019 by Brennan Young                            //
// - made avgAzimuth a template function.                          //
// modified 06/06/2019 by Brennan Young                            //
// - made correction to azimuth.                                   //
// modified 06/13/2019 by Brennan Young                            //
// - made correction to azimuth.                                   //
// modified 09/19/2019 by Brennan Young                            //
// - made the logic more specific in inPolygon.                    //
// modified 10/17/2019 by Brennan Young                            //
// - updated lineSlopeOrth to handle 0 slopes.                     //
// - added nearestXY to calculate the nearest XY to a line         //
//   segment.                                                      //
// modified 11/06/2019 by Brennan Young                            //
// - dist(...) for finding the distance between a point and a line //
//   now returns the correct distance.                             //
// modified 03/13/2020 by Brennan Young                            //
// - simplified azimuth -- no need to check for near-zero values   //
//   with atan2.                                                   //
// modified 04/24/2020 by Brennan Young                            //
// - fixed some compiler warnings.                                 //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// 
// CONTENTS
// 
// GLOBAL CONSTANTS
// 
// double pi = 3.14...
// 
// FUNCTIONS
// - T is template for primitive data types (int, float, double, etc.)
// 
// T max (T a, T b)                         larger of two values
// T min (T a, T b)                         smaller of two values
// double azimuth (double dx, double dy)    bearing from dx, dy
// double avgAzimuth (vector<double> az)    average azimuth bearing
// double avgAzimuth (double * az, int n)   average azimuth bearing
// double circularDiff (double a, double b, shortest distance around circle
//                      double m)
// double lineSlope (double dx, double dy)  compute m in y=mx+b
// double lineYint (double x, double y,     compute b in y=mx+b
//                  double m)
// double lineSlopeOrth (double m)          compute perpendicular slope
// void intercept (double m1, double b1,    intercept of two lines
//                 double xx1, double m2,
//                 double b2, double xx2,
//                 double * x, double * y)
// double dist (double dx, double dy)       distance between two points
// double dist (double xx, double yy,       distance between point and line
//              double m, double b,
//              double x)
// 
// TOPOLOGICAL OPERATIONS
// - P is template for points, which are structures with members x and y
// 
// bool less (P a, P b)                     a.x<b.x or (a.x==b.x and a.y<b.y)
// bool equal (P a, P b)                    x and y are the same
// size_t sortedXY_index (vector<P>, P)     index for vector sorted by x then y
// vector<P> convexHull_Jarvis              convex hull Jarvis march O(nh)
// vector<P> convexHull_chain               c. hull Anderson method O(n log n)
// vector<P> convexHull_Chan                convex hull Chan method O(n log h)
//                                            using Jarvis and chain methods
// bool inLine(P p1, P p2, P p3)            if p3 is on line p1p2
// bool inPolygon(vector<P> polygon, P p)   if p is in polygon
// 
///////////////////////////////////////////////////////////////////////////////

#ifndef YOUNG_STDLIB_GEOMETRY_20170914
#define YOUNG_STDLIB_GEOMETRY_20170914

#include <cmath>    // sqrt, fabs, fmod, tan, tan2
#include <stddef.h> // size_t
#include <vector>   // std::vector

namespace bystd { // Brennan Young standard namespace

// GLOBAL CONSTANTS ///////////////////////////////////////////////////////////

const double pi = 4.0 * atan(1.0);

// FUNCTIONS //////////////////////////////////////////////////////////////////

template <class T> T max(const T & a, const T & b){ return a < b ? b : a; }
template <class T> T min(const T & a, const T & b){ return a < b ? a : b; }

// Compute the azimuth of a vector (in radians) given the magnitude of the x
// and y component vectors (i.e., the direction that a is from b). In
// geographic context, +x and +y are understood to be east and north,
// respectively. For a - b, +dx = eastward (a is east of b) and +dy = northward
// (a is north of b).
double azimuth ( double dx, double dy )
{
    static const double pi2 = 2.0 * pi;
    double az = 0.5 * pi - atan2(dy, dx);
    while ( az < 0.0 ) az += pi2;
    while ( az >= pi2 ) az -= pi2;
    return az;
}

// Compute the average azimuth of a dataset of azimuth bearings (in radians).
// Average is understood to mean to 'middle' azimuth bearing in the same
// semicircle as the majority of the given azimuth bearings.
template <typename T>
double avgAzimuth(std::vector<T> az) {
    double sumCosAz = 0.0;
    double sumSinAz = 0.0;
    for ( unsigned int i = 0; i < az.size(); ++i ){
        sumCosAz += cos(az[i]);
        sumSinAz += sin(az[i]);
    }
    if ( fabs(sumCosAz) < 0.000001 && fabs(sumSinAz) < 0.000001 ){
        return -9999.0; // no solution
    }
    const double pi2 = 2.0 * pi;
    double result = atan2(sumSinAz, sumCosAz);
    if ( result < 0.0 ) result += pi2;
    return result < pi2 ? result : result - pi2;
}

template <typename T>
double avgAzimuth(T * az, int n) {
    double sumCosAz = 0.0;
    double sumSinAz = 0.0;
    for ( int i = 0; i < n; ++i ){
        sumCosAz += cos(az[i]);
        sumSinAz += sin(az[i]);
    }
    if ( fabs(sumCosAz) < 0.000001 && fabs(sumSinAz) < 0.000001 ){
        return -9999.0; // no solution
    }
    const double pi2 = 2.0 * pi;
    double result = atan2(sumSinAz, sumCosAz);
    if ( result < 0.0 ) result += pi2;
    return result < pi2 ? result : result - pi2;
}

// Compute the minimum distance about a circular scale between two values,
// given each value (a, b) and the circular scale (m).
double circularDiff(double a, double b, double m){
    while ( a >= m ) a -= m;
    while ( b >= m ) b -= m;
    double c = a - b;
    if ( c < 0.0 ) c *= -1.0;
    if ( c > m / 2.0 ) c = m - c;
    return c;
}

// Computes the slope (m) of the line defined by dx and dy, of the form y=mx+b.
// If the slope is undefined (vertical), returns -9999.0.
double lineSlope(double dx, double dy){
    if ( fabs(dx) < 0.000001 ) return -9999.0;
    return dy / dx;
}

// Computes the y-intercept (b) of the line defined by x, y, and m, of the form
// y=mx+b. If the y-intercept is undefined (slope is vertical = -9999.0),
// returns -9999.0.
double lineYint(double x, double y, double m){
    if ( m == -9999.0 ) return -9999.0;
    return y - m * x;
}

// Computes the slope of the line that is perpendicular to a line
// with slope (m), of the form y=mx+b. If the slope is undefined
// (vertical=-9999.0), returns 0. If the slope if 0, returns
// undefined (vertical=-9999.0).
double lineSlopeOrth(double m){
    if ( m == -9999.0 ) return 0.0;
    else if ( m == 0.0 ) return -9999.0;
    return -1.0 / m;
}

// Computes the intercept of two lines of the form y=mx+b and stores the point
// of interception in the given variables x and y. If they do not intersect
// (parallel), sets x and y to -9999.0. xx1 and xx2 are x coordinates that
// lines 1 and 2 pass through, respectively.
void intercept(double m1, double b1, double xx1, double m2, double b2,
        double xx2, double * x, double * y){
    double dm = m1 - m2;
    if ( dm < 0.0 ) dm *= -1.0;
    if ( dm < 0.000001 ){ // lines are parallel
        *x = -9999.0;
        *y = -9999.0;
    }
    else if ( m1 == -9999.0 ){ // line 1 is vertical
        *x = xx1;
        *y = m2 * (*x) + b2;
    }
    else if ( m2 == -9999.0 ){ // line 2 is vertical
        *x = xx2;
        *y = m1 * (*x) + b1;
    }
    else {
        *x = (b2 - b1) / (m1 - m2);
        *y = m1 * (*x) + b1;
    }
}

// Compute the magnitude of a 2D vector given the magnitude of the x and y
// component vectors.
double dist (double dx, double dy){ return sqrt(dx*dx + dy*dy); }

// Compute the shortest 2D distance between a point (xx, yy) and a line of the
// form y=mx+b that passes through a point with x-coordinate x.
double dist (double xx, double yy, double m, double b, double x){
    double mm = lineSlopeOrth(m);
    double bb = lineYint(xx, yy, mm);
    double xxx, yyy;
    intercept(m, b, x, mm, bb, xx, &xxx, &yyy);
    return dist(xx - xxx, yy - yyy);
}

// Compute the shortest 2D distance between a point (xx,yy) and the
// line segment between (x0,y0) and (x1,y1).
// (solving for the height of the triangle h = 2A/b, where the
// vertices are the given points.)
double dist ( double xx, double yy, double x0, double y0,
    double x1, double y1 )
{
    double dx = x1 - x0;
    double dy = y1 - y0;
    return fabs(dy * xx - dx * yy + x1 * y0 - y1 * x0)
        / sqrt(dy*dy + dx*dx);
}

// Computes the nearest x and y coordinates between a given point at
// (xx,yy) and the line defined by (x0,y0) and (x1,y1)
void nearestXY ( double xx, double yy, double x0, double y0,
    double x1, double y1, double * x, double * y )
{
    double m = lineSlope(x0-y0, x1-y1);
    double b = lineYint(x0, y0, m);
    double mm = lineSlopeOrth(m);
    double bb = lineYint(xx, yy, mm);
    intercept(m, b, x0, mm, bb, xx, x, y);
}

// TOPOLOGICAL OPERATIONS /////////////////////////////////////////////////////

// Comparing template points with x and y coordinates
template <class P> bool less(const P & a, const P & b){
    if ( a.x < b.x ) return true;
    if ( b.x < a.x ) return false;
    return a.y < b.y;
}
template <class P> bool equal(const P & a, const P & b){
    float d = a.x - b.x;
    if ( d < 0.0 ) d *= -1.0;
    if ( d > 0.0000001 ) return false;
    d = a.y - b.y;
    if ( d < 0.0 ) d *= -1.0;
    return d < 0.0000001;
}

// Algorithm for sorting a vector based on a point's x, then y values.
// Source: created 01/11/2018 by Brennan Young
template <class P> size_t sortedXY_index(const std::vector<P> & v, P p){
    if ( v.size() == 0 ){ return v.size(); }
    
    size_t i = v.size() / 2; // index being checked
    size_t ii;
    size_t li = 0;             // lower bound
    size_t ui = v.size();      // upper bound
    
    while ( less(v[i], p) || less(p, v[i]) ){
        if ( less(v[i], p) ){
            if ( i == v.size() - 1 ){ return v.size(); }
            li = i + 1;
        }
        if ( less(p, v[i]) ) ui = i;
        ii = li + (ui - li) / 2;
        if ( i == ii ) return i;
        i = ii;
    }
    
    return i;
}

// 2D cross product
template <class P>
double cross (const P & a, const P & b, const P & c){
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

// Jarvis' march ("gift-wrapping") algorithm.
// mode = 0 if S is unsorted by x amd then y.
template <class P>
std::vector<P> convexHull_Jarvis(const std::vector<P> & S, unsigned char mode){
    
    if ( S.size() <= 3 ) return S;
    
    int i, j;
    
    P pointOnHull = S[0]; // leftmost point
    
    if ( mode == 0 ){ // unsorted
        // find point with minimum x
        for ( i = 0; i < S.size(); ++i ){
            if ( S[i].x < pointOnHull.x ) pointOnHull = S[i];
        }
    }
    
    std::vector<P> H; // hull
    P endpoint = S[0];
    
    i = 0;
    do {
        H.push_back(pointOnHull);
        for ( j = 0; j < S.size(); ++j ){
            if ( equal(endpoint, pointOnHull)
                    || cross(H[i], endpoint, S[j]) < 0.0 )
                endpoint = S[j]; // greater left turn, update endpoint
        }
        ++i;
        pointOnHull = endpoint;
    } while ( !equal(endpoint, H[0]) ); // wrapped around to first point on hull
    
    return H;
}
template <class P> std::vector<P> convexHull_Jarvis(const std::vector<P> & S){
    return convexHull_Jarvis(S, 0);
}

// Andrew's algorithm for determining convex hull ("monotone chain").
// mode = 0 if S is unsorted by x and then y.
// Returns points counterclockwise, beginning with the point with minimum x (or
// minimum x and y, if tied for x).
// Source: Andrew, A. M. 1979. Another efficient algorithm for convex hulls in
//   in two dimensions. Information Processing Letters 9(5), 216-219.
// O (n log n)
template <class P>
std::vector<P> convexHull_chain(const std::vector<P> & S, unsigned char mode){
    
    if ( S.size() <= 3 ) return S;
    
    int i, k, t;
    std::vector<P> Q;
    
    if ( mode == 0 ){ // unsorted
        // sort all points by x then by y
        for ( i = 0; i < S.size(); ++i ){
            k = sortedXY_index(Q, S[i]);
            Q.insert(Q.begin() + k, S[i]);
        }
    }
    
    // build lower hull
    std::vector<P> H;
    H.resize(Q.size() * 2);
    k = 0;
    for ( int i = 0; i < Q.size(); ++i ){
        while ( k >= 2 && cross(H[k-2], H[k-1], Q[i]) <= 0 ) --k;
        H[k++] = Q[i];
    }
    
    // build upper hull
    for ( int i = Q.size() - 2, t = k+1; i >= 0; --i ){
        while ( k >= t && cross(H[k-2], H[k-1], Q[i]) <= 0 ) --k;
        H[k++] = Q[i];
    }
    
    H.resize(k - 1);
    return H;
}
template <class P>
std::vector<P> convexHull_chain(const std::vector<P> & S){
    return convexHull_chain(S, 0);
}

// Returns the convex hull of the given container of points using Chan's
// optimal output-sensitive algorithm.
// Source: Chan, T. M. 1996. Optimal output-sensitive convex hull algorithms in
//   two and three dimensions. Discrete and Computational Geometry 16(4), 361-
//   368.
// O(n log h)
template <class P>
std::vector<P> convexHull_Chan(const std::vector<P> & S){
    
    if ( S.size() <= 3 ) return S;
    
    int i, j, idx;
    
    int n = S.size();
    int m = S.size() / 4;
    if ( m == 0 ) m = S.size();
    
    std::vector< std::vector<P> > Q; // subsets
    Q.resize(1 + n / m);
    std::vector<P> H; // hull
    
    // divide points into subsets, sorted by x then y
    j = 0;
    for ( i = 0; i < n; ++i ){
        idx = sortedXY_index(Q[j], S[i]);
        Q[j].insert(Q[j].begin() + idx, S[i]);
        if ( Q[j].size() >= m ) ++j;
    }
    
    // compute convext hull by Andrew's monotone chain algorithm (stores
    // vertices in counterclockwise order)
    for ( i = 0; i < Q.size(); ++i ){
        Q[i] = convexHull_chain(Q[i], 1);
        H.insert(H.end(), Q[i].begin(), Q[i].end());
    }
    
    // Jarvis' march ("gift wrapping algorithm") on pre-computed convex hulls
    H = convexHull_Jarvis(H);
    
    return H;
}

template <class P>
std::vector<P> convexHull(const std::vector<P> & S){
    return convexHull_Chan(S);
}

// Returns true if the point is on the line.
// 
// Source: https://stackoverflow.com/questions/26849632/see-if-a-point-lies-on-a-linevector
// - (author) Mark Schlosser, edited by lennon310.
// - "An efficient way to solve this problem is to use the signed area of a
//   triangle. When the signed area of the triangle created by points {x1,y1},
//   {x2,y2}, and {x,y} is near-zero, you can consider {x,y} to be on the
//   line."
// - Modified 01/11/2018 by Brennan Young
//   - changed function name to inLine and variable name p3 to q.
//   - made function receive a template point structure, with MUST have members
//     x and y. This required different treatment of the variables va and vb,
//     which are removed, their components moved into computation of 'area'.
//   - replaced fabs with an if statement and multiplication, to reduce
//     dependencies.
//   - specified type of 'area' as double.
//   - added a check to see if the point lies in the bounds of the line segment
template <class P>
bool inLine(const P & p1, const P & p2, const P & q){
    // determine if the point is in the bounds of the line segment
    if ( !(q.x <= max(p1.x, p2.x) && q.x >= min(p1.x, p2.x)
            && q.y <= max(p1.y, p2.y) && q.y >= min(p1.y, p2.y)) )
        return false;
    
    // calculate signed area of triangle made by the three points
    double area = (p1.x - p2.x) * (q.y - p2.y)
        - (p1.y - p2.y) * (q.x - p2.x);
    if ( area < 0.0 ) area *= -1;
    
    return area < 0.0000001;
}

// Returns true if the point is inside the polygon. Note that division by zero
// is avoided because the division is protected by the "if" clause which
// surrounds it.
// 
// Source: http://alienryderflex.com/polygon/
// - (c) Darel Rex Finley
// - "This complete article, unmodified, may be freely distributed for
//   educational purposes."
// - Modified 01/11/2018 by Brennan Young: vector of points and a point as
//   arguments. Points are template structures which MUST have members x and y.
// - Modified 01/11/2018 by Brennan Young: Ambiguous result at the edges and
//   vertices remediated by checking if the point lies on a vertex or on an
//   edge first - if so, returns true.
template <class P>
bool inPolygon(const std::vector<P> & polygon, const P & p){
    int  i, j = polygon.size() - 1;
    bool oddNodes = false;
    
    for (i = 0; i < polygon.size(); ++i) {
        // special case: lies on edge or vertex
        if ( inLine(polygon[i], polygon[j], p) ) return true;
        
        // general case
        if ( ((polygon[i].y < p.y && polygon[j].y >= p.y)
                || (polygon[j].y < p.y && polygon[i].y >= p.y))
                && (polygon[i].x <= p.x || polygon[j].x <= p.x) ) {
            oddNodes ^=
                (polygon[i].x + (p.y - polygon[i].y) /
                (polygon[j].y - polygon[i].y) *
                (polygon[j].x - polygon[i].x) < p.x);
        }
    j = i;
    }
    
    return oddNodes;
}

// Returns the point at the visual center of the polygon in the xy plane.
// 
// Source: https://github.com/mapbox/polylabel
// - No code was taken from the source, only the concept for the algorithm.
// - (c) 2016 Mapbox (for the code, which was not used)
// - "Permission to use, copy, modify, and/or distribute this software for any
//   purpose with or without fee is hereby granted, provided that the above
//   copyright notice and this permission notice appear in all copies."
// <!> needs to be able to deal with holes:
//   std::vector<std::vector<P> > & polygon
template <class P>
P visualCenter(const std::vector<P> & polygon){
    //struct extent {
    //    double x0; double x1; double y0; double y1; double d; double p;
    //};
    
    // Container of extent objects (each also needs "distance" and "potential"
    // numbers)
    
    //std::vector<extent> extents;
    
    // 1. Initial extent as the polygon extent
    // 2. Compute centroid and its distance from the polygon edge = this is the
    //   initial "best distance"
    // 3. Divide each extent into four extents
    // 4. Compute distance of center to edge of polygon for each extent
    //   (negative if outside polygon)
    // 5. Compute potential for extent (center distance from polygon edge
    //   + distance from center to corner)
    // 6. If potential < best distance, discard
    // 7. Sort by potential and simultaneously find the least (best) distance
    // 8. Discard all extents where potential < best distance
    // 9. Repeat 3-8 until a solution is found at the desired precision
    
    P a;
    return a;
}

} // namespace bystd

#endif // YOUNG_STDLIB_GEOMETRY_20170914