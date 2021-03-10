/////////////////////////////////////////////////////////////////////
// RTree class to spatially index points, areas, volumes, or       //
// hypervolumes. Inserted entities must have ID >= 0. The maximum  //
// number of dimensions that can be represented in the RTree is    //
// 255, and the maximum number of children per node is 255.        //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// HISTORY ------------------------------------------------------- //
// 06/10/2019 - Created - Brennan Young                            //
// 06/18/2019 - Modified - Brennan Young                           //
// - find functions are now all const.                             //
// - added find within radius.                                     //
// - added find between two lines.                                 //
// 06/21/2019 - Modified - Brennan Young                           //
// - replaced find between two lines with find within an arc and   //
//   find within a polygon.                                        //
// 07/12/2019 - Modified - Brennan Young                           //
// - moved member function definitions outside of class            //
//   definition, for nested classes Box and Page.                  //
// 07/15/2019 - Modified - Brennan Young                           //
// - moved member function definitions outside of class            //
//   definition.                                                   //
// 08/01/2019 - Modified - Brennan Young                           //
// - now includes the cmath library.                               //
// 08/02/2019 - Modified - Brennan Young                           //
// - added size method.                                            //
// 10/09/2019 - Modified - Brennan Young                           //
// - added findExtent method.                                      //
// 10/28/2019 - Modified - Brennan Young                           //
// - worked on memory leak associated with Page.                   //
// 01/27/2020 - Declared Unstable - Brennan Young                  //
// - need to figure out memory issues, probably need to completely //
//   restructure the RTree.                                        //
// - changed page children to a std::vector... because dynamically //
//   allocated memory is giving me a headache. Fixed the problem.  //
// 01/31/2020 - Modified - Brennan Young                           //
// - added RTree::Box::dist and sqdist to compute the minimum      //
//   distance between boxes.                                       //
// - made nearestNeighbors work.                                   //
// 03/13/2020 - Modified - Brennan Young                           //
// - removed RTree::minsqdist, maxsqdist, and sqdist.              //
// - fixed find within radius.                                     //
// 03/19/2020 - Modified - Brennan Young                           //
// - RTree::Box::overlap now returns the amount of overlap, rather //
//   than binary true/false overlap (returns 1 if point-overlap).  //
// - removed RTree::Page::volChange, in favor of using             //
//   RTree::Box::overlap.                                          //
// - implemented R* tree-like insertion criteria.                  //
// 03/20/2020 - Modified - Brennan Young                           //
// - fixed a memory leak: Box wasn't deallocating with assignment  //
//   operator!                                                     //
// - fixed insertion criteria to be more like R* tree.             //
// 04/07/2020 - Modified - Brennan Young                           //
// - intersect now correctly always returns a value.               //
// - corrected compiler warnings.                                  //
// 04/08/2020 - Modified - Brennan Young                           //
// - corrected some minor errors and style.                        //
// 04/22/2020 - Modified - Brennan Young                           //
// - enabled duplicate IDs, to permit representation of sub-       //
//   divided objects.                                              //
// 04/23/2020 - Modified - Brennan Young                           //
// - modified to manage tree structure with pointers instead of a  //
//   vector.                                                       //
// - altered insert and remove operations to internally trade and  //
//   modify pointers.                                              //
// - removed findExtent - object extent should be tracked external //
//   to this data structure.                                       //
// 07/31/2020 - Modified - Brennan Young                           //
// - corrected some minor compiler warnings.                       //
// 10/01/2020 - Modified - Brennan Young                           //
// - added findLine.                                               //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- NOTES ------------------------------------------------------ //
// - Need to test node removal.                                    //
// - Need more complete R* tree insertion and removal              //
//   implementation.                                               //
/////////////////////////////////////////////////////////////////////

#ifndef RTREE_H
#define RTREE_H

#include <cmath>
#include <set>
#include <vector>

class RTree {
public:
    // bounding box structure
    class Box {
    public:
        unsigned char n;    // number of dimensions
        double* a;         // minimum bound
        double* b;         // maximum bound
        
        // constructors / destructor
        Box(unsigned char);
        Box(double, double); // 1D
        Box(double, double, double, double); // 2D
        Box(double, double, double, double, double, double); // 3D
        Box(const Box&);
        ~Box();
        
        // operators
        Box& operator= (const Box&);    // assignment
        
        // operations
        double volume() const;              // compute hypervolume
        double sqdist(const Box&) const;    // sq dist btwn boxes
        double dist(const Box&) const;      // dist between boxes
        double overlap(const Box&) const;   // intersection volume
        void expand(const Box&);            // include new box
    }; // Box
    
private:
    // node structure
    class Page {
    public:
        int id;                 // -1 unless a leaf node
        Box box;                // mininimum bounding volume
        unsigned char n;        // current number of children
        unsigned char nC;       // maximum number of children
        int size;               // number descendents + self
        Page** children;        // child nodes
        
        // constructors / destructor
        Page(unsigned char, int, const RTree::Box&);
        Page(const Page&);
        ~Page();
        
        // operators
        Page& operator=(const Page&);
        bool operator<(const Page&) const;
        
        // operations
        void update();                  // update based on children
    }; // Page
    
    // structure to help organize pages by some value
    struct PageVal {
        RTree::Page * page;
        double val;
        
        PageVal(RTree::Page* p=NULL, double v=0.0)
            : page(p), val(v) {}
        bool operator<(const PageVal& v) const {
            if ( val < v.val ) return true;
            if ( val > v.val ) return false;
            return page < v.page;
        }
    }; // Pageval
    
private:
    unsigned int nD;    // number of dimensions represented in tree
    unsigned int nC;    // maximum number of children in each page
    Page * root;        // root node
    
    // helper functions
    bool intersection(double*, double*, // get x,y where lines cross
        double, double, double, double,
        double, double) const;
    bool rayIntersect(const Box&,       // get ray intersects box
        const Box&, std::vector<double>*) const;
    bool polygonIntersect(double*,      // determine if polys overlap
        double*, int, double*, double*, int) const;
    
    // query
    void find (std::set<int>*,          // IDs overlapping box
        const Box&,
        Page const * const) const;
    void find(std::set<int>*,           // IDs within radius
        const Box&, double,
        Page const * const) const;
    void findInPolygon(std::set<int>*,  // IDs overlapping polygon
        double*, double*, int, int,
        int, Page const * const) const;
    void findRay(std::vector<int>*,     // IDs intersecting ray
        std::vector<double>*,
        std::vector<std::vector<double> >*,
        const Box&, Page const * const ) const;
    void findLine(std::vector<int>*,    // IDs intersecting line
        std::vector<double>*,
        std::vector<std::vector<double> >*,
        const Box&, Page const * const ) const;
    
    // insertion / removal
    unsigned int insertChooseSubtree(   // choose where to insert
        const Box&, Page const * const) const;
    Page* insert(int, const Box&, Page*);   // add to tree
    Page* remove(int, const Box&, Page*,    // remove from tree
        std::set<Page>*);
    
    // I/O
    void printTree (Page const * const, unsigned int) const;
public:
    // constructors / destructor
    RTree (unsigned char d=1, unsigned char c=3);
    RTree (const RTree&);
    ~RTree ();
    
    //operators
    RTree & operator= (const RTree&);
    
    // getters
    int size () const;                  // get number of objects
    int dimensions () const;            // get nD
    int pageSize () const;              // get nC
    Box extent () const;                // get root box
    
    // retrieval
    std::set<int> find (const Box&)     // IDs overlapping volume
        const;
    std::set<int> find (const Box&,     // IDs in radius
        double r ) const;
    std::set<int> find (int, double,    // IDs in range in dimension
        double) const;
    std::set<int> findInPolygon (       // IDs overlapping polygon
        double*, double*, int, int,
        int) const;
    std::vector<int> findRay (          // IDs intersecting ray
        const Box&,
        std::vector<std::vector<double> >*) const;
    std::vector<int> findLine (         // IDs intersecting line
        const Box&,
        std::vector<std::vector<double> >*) const;
    std::set<int> findInArc (double,    // IDs in arc in dimension
        double, double, double, int,
        int) const;
    std::vector<int> nearestNeighbors ( // nearest k IDs
        const Box&, unsigned int) const;
    
    // insertion / removal
    void insert(int, const Box&);
    void remove(int, const Box&);
    void update(int, const Box&, const Box&);
    
    // I/O
    void printTree() const;
}; // RTree


/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// Box                                                             //
/////////////////////////////////////////////////////////////////////

// Constructors / Destructor ////////////////////////////////////////

// Default constructor.
RTree::Box::Box ( unsigned char nD=1 ) : n(nD)
{
    a  = new double [n];
    b  = new double [n];
    for ( int i = 0; i < n; ++i ) {
        a[i] = 0.0;
        b[i] = 0.0;
    }
}

// 1D constructor.
RTree::Box::Box ( double x0, double x1 ) : n(1)
{
    a = new double [n];
    b = new double [n];
    a[0] = x0;
    b[0] = x1;
}

// 2D constructor.
RTree::Box::Box ( double x0, double x1, double y0, double y1 ) : n(2)
{
    a = new double [n];
    b = new double [n];
    a[0] = x0;
    b[0] = x1;
    a[1] = y0;
    b[1] = y1;
}

// 3D constructor.
RTree::Box::Box ( double x0, double x1, double y0, double y1,
    double z0, double z1 ) : n(3)
{
    a = new double [n];
    b = new double [n];
    a[0] = x0;
    b[0] = x1;
    a[1] = y0;
    b[1] = y1;
    a[2] = z0;
    b[2] = z1;
}

// Copy constructor.
RTree::Box::Box ( const RTree::Box & B ) : n(B.n)
{
    a = new double[n];
    b = new double[n];
    for ( int i = 0; i < n; ++i ) {
        a[i] = B.a[i];
        b[i] = B.b[i];
    }
}

// Destructor.
RTree::Box::~Box ()
{
    delete [] a;
    delete [] b;
}

// Operators ////////////////////////////////////////////////////////

RTree::Box & RTree::Box::operator= ( const RTree::Box & B )
{
    delete [] a;
    delete [] b;
    
    n = B.n;
    a = new double [n];
    b = new double [n];
    for ( int i = 0; i < n; ++i ) {
        a[i] = B.a[i];
        b[i] = B.b[i];
    }
    
    return *this;
}

// Operations ///////////////////////////////////////////////////////

// Compute the box's volume (or hypervolume).
double RTree::Box::volume () const
{
    double v = b[0] - a[0];
    for ( int i = 1; i < n; ++i ) v *= (b[i] - a[i]);
    return v;
}

// Compute the square of the distance between this and the given box.
double RTree::Box::sqdist ( const RTree::Box & B ) const
{
    double d, sum = 0.0;
    for ( int i = 0; i < n && i < B.n; ++i ) {
        if ( b[i] < B.a[i] ) d = B.a[i] - b[i];
        else if ( B.b[i] < a[i] ) d = a[i] - B.b[i];
        else d = 0.0;
        sum += d*d;
    }
    return sum;
}

// Compute distance between this and the given box.
double RTree::Box::dist ( const RTree::Box & B ) const
{
    return sqrt(this->sqdist(B));
}

// Compute the overlap hypervolume at the intersection between boxes.
// Returns a 1 if either box is a point within the extent of the
// other.
double RTree::Box::overlap ( const RTree::Box & B ) const
{
    double amax, bmin, ab, sum = 1.0;
    for ( int i = 0; i < n; ++i ) {
        // check if any overlap
        if ( b[i] < B.a[i] ) return 0.0;
        if ( B.b[i] < a[i] ) return 0.0;
        
        // get overlap in dimension i
        amax = a[i] > B.a[i] ? a[i] : B.a[i];
        bmin = b[i] < B.b[i] ? b[i] : B.b[i];
        ab = bmin - amax;
        if ( ab <= 0.0 ) continue; // is a point
        sum *= ab;
    }
    return sum;
}

// Expand to include new box.
void RTree::Box::expand ( const RTree::Box & B )
{
    for ( int i = 0; i < n; ++i ) {
        if ( B.a[i] < a[i] ) a[i] = B.a[i];
        if ( B.b[i] > b[i] ) b[i] = B.b[i];
    }
}


/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// Page                                                            //
/////////////////////////////////////////////////////////////////////

// Constructors / Destructor ////////////////////////////////////////

// Default constructor
RTree::Page::Page ( unsigned char n=3, int i=-1,
    const RTree::Box & b=RTree::Box() )
: id(i), box(b), n(0), nC(n), size(1)
{
    if ( nC < 2 ) nC = 2;
    children = new Page* [nC];
    for ( unsigned int i = 0; i < nC; ++i ) children[i] = NULL;
}

// Copy constructor.
RTree::Page::Page ( const RTree::Page & p )
: id(p.id), box(p.box), n(p.n), nC(p.nC), size(p.size)
{
    children = new Page* [nC];
    unsigned int i = 0;
    for ( ; i < n; ++i ) children[i] = new Page(*p.children[i]);
    for ( ; i < nC; ++i ) children[i] = NULL;
}

// Destructor
RTree::Page::~Page ()
{
    for ( unsigned int i = 0; i < n; ++i ) delete children[i];
    delete [] children;
}

// Operators ////////////////////////////////////////////////////////

// Assignment operator.
RTree::Page & RTree::Page::operator= ( const RTree::Page & p )
{
    if ( &p == this ) return *this;
    
    for ( unsigned int i = 0; i < n; ++i ) delete children[i];
    delete [] children;
    
    id = p.id;
    box = p.box;
    n = p.n;
    nC = p.nC;
    size = p.size;
    
    children = new Page* [nC];
    unsigned int i = 0;
    for ( ; i < n; ++i ) children[i] = new Page(*p.children[i]);
    for ( ; i < nC; ++i ) children[i] = NULL;
    
    return *this;
}

// Less-than comparison operator.
bool RTree::Page::operator< ( const RTree::Page & p ) const
{
    return id < p.id;
}

// Operations ///////////////////////////////////////////////////////

// Determine bounding box based on children.
void RTree::Page::update ()
{
    unsigned char i;
    unsigned char j;
    double a, b;
    
    size = id >= 0 ? 1 : 0;
    
    if ( n == 0 ) return;
    
    // bounding box
    for ( i = 0; i < box.n; ++i ) {
        a = children[0]->box.a[i];
        b = children[0]->box.b[i];
        for ( j = 1; j < n; ++j ) {
            if ( children[j]->box.a[i] < a )
                a = children[j]->box.a[i];
            if ( children[j]->box.b[i] > b )
                b = children[j]->box.b[i];
        }
        box.a[i] = a;
        box.b[i] = b;
    }
    
    // size of sub-tree
    for ( i = 0; i < n; ++i ) size += children[i]->size;
}

/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// RTree                                                           //
/////////////////////////////////////////////////////////////////////

// Constructor / Destructor /////////////////////////////////////////

// Default constructor.
RTree::RTree ( unsigned char d, unsigned char c )
    : nD(d), nC(c), root(NULL) {}

// Copy constructor.
RTree::RTree ( const RTree & r ) : nD(r.nD), nC(r.nC), root(NULL)
{
    if ( r.root != NULL ) root = new Page (*r.root);
}

// Destructor.
RTree::~RTree () { if ( root != NULL ) delete root; }

// Operators ////////////////////////////////////////////////////////

// Assignment operator.
RTree & RTree::operator= ( const RTree & r )
{
    nD = r.nD;
    nC = r.nC;
    if ( root != NULL ) delete root;
    if ( r.root != NULL ) root = new Page (*r.root);
    
    return *this;
}

// Helper Functions /////////////////////////////////////////////////

// Get the point of intersection between two lines. Use m < -9000
// to represent vertical lines. Returns false if no intersection.
bool RTree::intersection ( double * x, double * y, double m1,
    double x1, double b1, double m2, double x2, double b2 ) const
{
    if ( (m1 < -9000 && m2 < -9000) || m1 == m2 ) return false;
    if ( m1 < -9000 ) {
        *x = x1;
        *y = m2 * (*x) + b2;
    }
    else if ( m2 < -9000 ) {
        *x = x2;
        *y = m1 * (*x) + b1;
    }
    else {
        *x = (b2 - b1) / (m1 - m2);
        *y = m1 * (*x) + b1;
    }
    return true;
}

// Returns if the ray intersects the box. The ray is defined as a
// box, where the a coordinates are treated as the ray's origin and
// the b coordinates as the ray's direction. Stores the coordinates
// of the intersection point in the hit vector, of the same length
// as the number of dimensions in the given box.
bool RTree::rayIntersect ( const RTree::Box & ray,
    const RTree::Box & box, std::vector<double>* hit ) const
{
    static const unsigned char LEFT = 0;
    static const unsigned char MIDDLE = 1;
    static const unsigned char RIGHT = 2;
    
    bool inside = true;
    unsigned char * quadrant = new unsigned char [box.n];
    int i;
    int whichPlane;
    double * maxT = new double [box.n];
    double * candidatePlane = new double [box.n];
    
    hit->resize(box.n);
    for ( i = 0; i < box.n; ++i ) (*hit)[i] = 0.0;
    
    // find candidate planes
    for ( i = 0; i < box.n; ++i ) {
        if ( ray.a[i] < box.a[i] ) {
            quadrant[i] = LEFT;
            candidatePlane[i] = box.a[i];
            inside = false;
        }
        else if ( ray.a[i] > box.b[i] ) {
            quadrant[i] = RIGHT;
            candidatePlane[i] = box.b[i];
            inside = false;
        }
        else quadrant[i] = MIDDLE;
    }
    
    // ray origin inside bounding box
    if ( inside ) {
        for ( i = 0; i < ray.n; ++i ) (*hit)[i] = ray.a[i];
        return true;
    }
    
    // calculate T distances to candidate planes
    for ( i = 0; i < box.n; ++i ) {
        if ( quadrant[i] != MIDDLE && ray.b[i] != 0.0 )
            maxT[i] = (candidatePlane[i] - ray.a[i]) / ray.b[i];
        else maxT[i] = -1.0;
    }
    
    // get largest of the maxT's for final choice of intersection
    whichPlane = 0;
    for ( i = 1; i < box.n; ++i )
        if ( maxT[whichPlane] < maxT[i] ) whichPlane = i;
    
    // check that the final candidate is actually inside the box
    if ( maxT[whichPlane] < 0.0 ) return false;
    for ( i = 0; i < box.n; ++i ) {
        if ( whichPlane != i ) {
            (*hit)[i] = ray.a[i] + maxT[whichPlane] * ray.b[i];
            if ( (*hit)[i] < box.a[i] || (*hit)[i] > box.b[i] )
                return false;
        }
        else (*hit)[i] = candidatePlane[i];
    }
    
    return true;
}

// Return true if the two polygons intersect using the separating
// axis theorem. Polygons are assumed to be convex and vertices
// ordered counter-clockwise.
bool RTree::polygonIntersect ( double * xa, double * ya, int na,
    double * xb, double * yb, int nb ) const
{
    int i, ii, j;
    double dx1, dy1, dx2, dy2, mag, dxn, dyn, len = 0.0;
    
    // check from faces of a
    for ( i = 0; i < na; ++i ) {
        // get vector for axis parallel to edge i, i+1
        ii = (i + 1) % na;
        dx1 = xa[ii] - xa[i];
        dy1 = ya[ii] - ya[i];
        
        // get normal vector (right)
        dxn = dy1;
        dyn = -dx1;
        mag = sqrt(dxn*dxn + dyn*dyn);
        
        if ( mag == 0.0 ) continue;
        
        // normalize normal vector
        dxn = dxn / mag;
        dyn = dyn / mag;
        
        // project vertices of b onto normal vector
        for ( j = 0; j < nb; ++j ) {
            dx2 = xb[j] - xa[i];
            dy2 = yb[j] - ya[i];
            len = dx2 * dxn + dy2 * dyn;
            
            if ( len <= 0.0 ) break;    // overlap
        }
        if ( len > 0.0 ) return false;  // no overlap
    }
    
    // check from faces of b
    for ( i = 0; i < nb; ++i ) {
        // get vector for axis parallel to edge i, i+1
        ii = (i + 1) % na;
        dx1 = xb[ii] - xb[i];
        dy1 = yb[ii] - yb[i];
        
        // get normal vector (right)
        dxn = dy1;
        dyn = -dx1;
        mag = sqrt(dxn*dxn + dyn*dyn);
        
        if ( mag == 0.0 ) continue;
        
        // normalize normal vector
        dxn = dxn / mag;
        dyn = dyn / mag;
        
        // project vertices of b onto normal vector
        for ( j = 0; j < na; ++j ) {
            dx2 = xa[j] - xb[i];
            dy2 = ya[j] - yb[i];
            len = dx2 * dxn + dy2 * dyn;
            
            if ( len <= 0.0 ) break;    // overlap
        }
        if ( len > 0.0 ) return false;  // no overlap
    }
    
    return true;
}

// Getters //////////////////////////////////////////////////////////

// Return the number of objects in the tree.
int RTree::size () const
{
    if ( root == NULL ) return 0;
    return root->size;
}

// Return the number of dimensions represented in the tree.
int RTree::dimensions () const
{
    return nD;
}

// Return the maximum number of children in each node.
int RTree::pageSize () const
{
    return nC;
}

// Return the box that contains the index.
RTree::Box RTree::extent () const
{
    if ( root == NULL ) return Box(nD);
    return root->box;
}

// Retrieval ////////////////////////////////////////////////////////

// Get all IDs that overlap the given volume.
void RTree::find ( std::set<int> * out, const RTree::Box & box,
    RTree::Page const * const p ) const
{
    if ( p == NULL || box.overlap(p->box) == 0.0 ) return;
    if ( p->id >= 0 ) out->insert(p->id);
    for ( unsigned int i = 0; i < p->n; ++i )
        find(out, box, p->children[i]);
}

// Get all IDs that overlap within the given radius.
void RTree::find ( std::set<int> * out, const RTree::Box & box,
    double r, RTree::Page const * const p ) const
{
    // check if the distance between boxes is within radius
    if ( p->box.sqdist(box) > r*r ) return;
    
    if ( p->id >= 0 ) out->insert(p->id);
    for ( unsigned int i = 0; i < p->n; ++i )
        find(out, box, r, p->children[i]);
}

// Get all IDs that lie within a given n-vertex polygon (counter-
// clockwise order).
void RTree::findInPolygon ( std::set<int> * out, double * x,
    double * y, int n, int Dx, int Dy, RTree::Page const * const p )
    const
{
    if ( p == NULL ) return;
    
    double boxX [4];
    double boxY [4];
    boxX[0] = p->box.a[Dx]; boxY[0] = p->box.a[Dy];
    boxX[1] = p->box.b[Dx]; boxY[1] = p->box.a[Dy];
    boxX[2] = p->box.b[Dx]; boxY[2] = p->box.b[Dy];
    boxX[3] = p->box.a[Dx]; boxY[3] = p->box.b[Dy];
    
    if ( !polygonIntersect(boxX, boxY, 4, x, y, n) ) return;
    
    if ( p->id >= 0 ) out->insert(p->id);
    for ( unsigned int i = 0; i < p->n; ++i )
        findInPolygon(out, x, y, n, Dx, Dy, p->children[i]);
}

// Get all IDs that intersect the given ray (represented by a box,
// where a is the ray's origin and b is the ray's direction), in
// order of distance from the ray's origin.
void RTree::findRay ( std::vector<int>* out,
    std::vector<double>* dist,
    std::vector<std::vector<double> >* hit,
    const RTree::Box & ray,
    RTree::Page const * const p ) const
{
    if ( p == NULL ) return;
    
    std::vector<double> temphit;
    if ( !rayIntersect(ray, p->box, &temphit) ) return;
    
    if ( p->id >= 0 ) {
        // get distance to intersection point
        double di, d = 0.0;
        for ( unsigned int i = 0; i < nD; ++i ) {
            di = ray.a[i] - temphit[i];
            d += (di * di); // square distance (avoid sqrt)
        }
        
        // insert into distance-sorted vector
        int k = 0;
        if ( out->size() > 0 ) {
            int u = dist->size();
            int b = 0;
            int p = -1;
            k = u / 2;
            
            while ( k != p ) {
                if ( d > (*dist)[k] ) {
                    if ( k == (int) dist->size() - 1 ) {
                        k = dist->size();
                        break;
                    }
                    b = k + 1;
                }
                else if ( d < (*dist)[k] ) u = k;
                else break;
                
                p = k;
                k = (b + u) / 2;
            }
        }
        out->insert(out->begin() + k, p->id);
        dist->insert(dist->begin() + k, d);
        hit->insert(hit->begin() + k, temphit);
    }
    for ( unsigned int i = 0; i < p->n; ++i )
        findRay(out, dist, hit, ray, p->children[i]);
}

// Get all IDs that intersect the given line (represented by a box,
// where a is the line's origin and b is the line's destination), in
// order of distance from the line's origin.
void RTree::findLine ( std::vector<int>* out,
    std::vector<double>* dist,
    std::vector<std::vector<double> >* hit,
    const RTree::Box & line,
    RTree::Page const * const p ) const
{
    if ( p == NULL ) return;
    
    std::vector<double> temphit;
    if ( !rayIntersect(line, p->box, &temphit) ) return;
    
    // sq length of line
    double len2 = 0.0;
    for ( unsigned int i = 0; i < nD; ++i )
        len2 += (line.b[i] - line.a[i]) * (line.b[i] - line.a[i]);
    
    if ( p->id >= 0 ) {
        // get distance to intersection point
        double di, d = 0.0;
        for ( unsigned int i = 0; i < nD; ++i ) {
            di = line.a[i] - temphit[i];
            d += (di * di); // square distance (avoid sqrt)
        }
        
        // if doesn't intersect line, return
        if ( d > len2 ) return;
        
        // insert into distance-sorted vector
        int k = 0;
        if ( out->size() > 0 ) {
            int u = dist->size();
            int b = 0;
            int p = -1;
            k = u / 2;
            
            while ( k != p ) {
                if ( d > (*dist)[k] ) {
                    if ( k == (int) dist->size() - 1 ) {
                        k = dist->size();
                        break;
                    }
                    b = k + 1;
                }
                else if ( d < (*dist)[k] ) u = k;
                else break;
                
                p = k;
                k = (b + u) / 2;
            }
        }
        out->insert(out->begin() + k, p->id);
        dist->insert(dist->begin() + k, d);
        hit->insert(hit->begin() + k, temphit);
    }
    for ( unsigned int i = 0; i < p->n; ++i )
        findLine(out, dist, hit, line, p->children[i]);
}

// Return a set of ids found within a space.
std::set<int> RTree::find ( const RTree::Box & box ) const
{
    std::set<int> out;
    find(&out, box, root);
    return out;
}

// Return a set of ids found within radius of given point.
std::set<int> RTree::find ( const RTree::Box & box, double r ) const
{
    std::set<int> out;
    find(&out, box, r, root);
    return out;
}

// Return a set of ids found between a and b (a<b) in dimension d.
std::set<int> RTree::find ( int d, double a, double b ) const
{
    std::set<int> out;
    if ( root == NULL ) return out;
    
    RTree::Box box (nD);
    for ( int i = 0; i < d; ++i ) {
        box.a[i] = root->box.a[i];
        box.b[i] = root->box.b[i];
    }
    box.a[d] = a;
    box.b[d] = b;
    for ( int i = d+1; (unsigned int) i < nD; ++i ) {
        box.a[i] = root->box.a[i];
        box.b[i] = root->box.b[i];
    }
    
    find(&out, box, root);
    return out;
}

// Return a set of nodes found within the given n-vertex polygon
// in dimensions Dx and Dy. Assumes vertices are in counter-
// clockwise order.
std::set<int> RTree::findInPolygon ( double * x, double * y, int n,
    int Dx, int Dy ) const
{
    std::set<int> out;
    findInPolygon(&out, x, y, n, Dx, Dy, root);
    return out;
}

// Return a vector of nodes found intersecting a ray, in order of
// nearest-to-furthest. The ray is represented with a box, which has
// its origin in a and its direction in b. Fills the hit vector with
// the location of the intersection.
std::vector<int> RTree::findRay ( const RTree::Box& ray,
    std::vector<std::vector<double> >* hit ) const
{
    std::vector<int> out;
    std::vector<double> dist;
    findRay(&out, &dist, hit, ray, root);
    return out;
}

// Return a vector of nodes found intersecting a line, in order of
// nearest-to-furthest. The line is represented with a box, which has
// its origin in a and its destination in b. Fills the hit vector
// with the location of the intersection.
std::vector<int> RTree::findLine ( const RTree::Box& line,
    std::vector<std::vector<double> >* hit ) const
{
    std::vector<int> out;
    std::vector<double> dist;
    findLine(&out, &dist, hit, line, root);
    return out;
}

// Return a set of object IDs found within the given arc defined
// by a point (xc, yc) and clockwise of dir1 and counter-
// clockwise of dir2 [radians] (0 = +y, increasing clockwise),
// within the plane defined by dimensions Dx and Dy.
std::set<int> RTree::findInArc ( double xc, double yc, double dir1,
    double dir2, int Dx, int Dy ) const
{
    std::set<int> out;
    if ( root == NULL ) return out;
    
    double x [6];
    double y [6];
    int n = 0;
    
    x[0] = xc;
    y[0] = yc;
    
    // get x, y components and slopes of each direction
    double Ydir1 = cos(dir1);
    double Xdir1 = sin(dir1);
    double m1, b1;
    if ( Xdir1 == 0.0 ) {
        m1 = -9999.0;
        b1 = 0.0;
    }
    else {
        m1 = Ydir1 / Xdir1;
        b1 = yc - m1 * xc;
    }
    
    double Ydir2 = cos(dir2);
    double Xdir2 = sin(dir2);
    double m2, b2;
    if ( Xdir2 == 0.0 ) {
        m2 = -9999.0;
        b2 = 0.0;
    }
    else {
        m2 = Ydir2 / Xdir2;
        b2 = yc - m2 * xc;
    }
    
    // get intersection with bounding rectangle
    if ( Ydir2 < 0.0 )
        intersection(&x[1], &y[1],
            m2, xc, b2, 0.0, 0.0, root->box.a[Dy]);
    else if ( Ydir2 > 0.0 )
        intersection(&x[1], &y[1],
            m2, xc, b2, 0.0, 0.0, root->box.b[Dy]);
    else if ( Xdir2 < 0.0 )
        intersection(&x[1], &y[1],
            m2, xc, b2, -9999.0, root->box.a[Dx], 0.0);
    else intersection(&x[1], &y[1],
            m2, xc, b2, -9999.0, root->box.b[Dx], 0.0);
    
    if ( Ydir1 < 0.0 )
        intersection(&x[2], &y[2],
            m1, xc, b1, 0.0, 0.0, root->box.a[Dy]);
    else if ( Ydir1 > 0.0 )
        intersection(&x[2], &y[2],
            m1, xc, b1, 0.0, 0.0, root->box.b[Dy]);
    else if ( Xdir1 < 0.0 )
        intersection(&x[2], &y[2],
            m1, xc, b1, -9999.0, root->box.a[Dx], 0.0);
    else intersection(&x[2], &y[2],
            m1, xc, b1, -9999.0, root->box.b[Dx], 0.0);
    
    // get last vertices to fill out arc
    n = 3;
    double ex1, ey1, ex2, ey2;
    for ( int i = 1; i < 4; ++i ) {
        // relative positions on edges
        ex1 = (x[i] - root->box.a[Dx])
            / (root->box.b[Dx] - root->box.a[Dx]);
        ey1 = (y[i] - root->box.a[Dy])
            / (root->box.b[Dy] - root->box.a[Dy]);
        ex2 = (x[i+1] - root->box.a[Dx])
            / (root->box.b[Dx] - root->box.a[Dx]);
        ey2 = (y[i+1] - root->box.a[Dy])
            / (root->box.b[Dy] - root->box.a[Dy]);
        
        // vertex i on (-y) edge
        if ( ey1 <= 0.0 ) {
            // vertex i+1 ahead on same edge
            if ( ey2 <= 0.0 && ex2 >= ex1 ) break;
            
            // vertex i at (+x, -y) corner
            if ( ex1 >= 1.0 ) {
                // vertex i+1 on (+x) edge
                if ( ex2 >= 1.0 ) break;
                else {
                    // push last intersection vertex back
                    ++n;
                    x[i+2] = x[i+1];
                    y[i+2] = y[i+1];
                    
                    // add new vertex at (+x, +y) corner
                    x[i+1] = root->box.b[Dx];
                    y[i+1] = root->box.b[Dy];
                }
            }
            else {
                // push last intersection vertex back
                ++n;
                x[i+2] = x[i+1];
                y[i+2] = y[i+1];
                
                // add new vertex at (+x,-y) corner
                x[i+1] = root->box.b[Dx];
                y[i+1] = root->box.a[Dy];
            }
        }
        
        // vertex i on (+x) edge
        else if ( ex1 >= 1.0 ) {
            // vertex i+1 ahead on same edge
            if ( ex2 >= 0 && ey2 >= ey1 ) break;
            
            // vertex i at (+x, +y) corner
            if ( ey1 >= 1.0 ) {
                // vertex i+1 on (+y) edge
                if ( ey2 >= 1.0 ) break;
                else {
                    // push last intersection vertex back
                    ++n;
                    x[i+2] = x[i+1];
                    y[i+2] = y[i+1];
                    
                    // add new vertex at (-x, +y) corner
                    x[i+1] = root->box.a[Dx];
                    y[i+1] = root->box.b[Dy];
                }
            }
            else {
                // push last intersection vertex back
                ++n;
                x[i+2] = x[i+1];
                y[i+2] = y[i+1];
                
                // add new vertex at (+x, +y) corner
                x[i+1] = root->box.b[Dx];
                y[i+1] = root->box.b[Dy];
            }
        }
        
        // vertex i on (+y) edge
        else if ( ey1 >= 1.0 ) {
            // vertex i+1 ahead on same edge
            if ( ey2 >= 1.0 && ex2 <= ex1 ) break;
            
            // vertex i at (-x, +y corner)
            if ( ex1 <= 0.0 ) {
                // vertex i+1 on (-x) edge
                if ( ex2 <= 0.0 ) break;
                else {
                    // push last intersection vertex back
                    ++n;
                    x[i+2] = x[i+1];
                    y[i+2] = y[i+1];
                    
                    // add new vertex at (-x, -y) corner
                    x[i+1] = root->box.a[Dx];
                    y[i+1] = root->box.a[Dy];
                }
            }
            else {
                // push last intersection vertex back
                ++n;
                x[i+2] = x[i+1];
                y[i+2] = y[i+1];
                
                // add new vertex at (-x, +y) corner
                x[i+1] = root->box.a[Dx];
                y[i+1] = root->box.b[Dy];
            }
        }
        
        // vertex on (-x) edge
        else {
            // vertex i+1 ahead on same edge
            if ( ex2 <= 0.0 && ey2 <= ey1 ) break;
            
            // vertex i at (-x, -y) corner
            if ( ey1 <= 0.0 ) {
                // vertex i+1 on (-y) edge
                if ( ey2 <= 0.0 ) break;
                else {
                    // push last intersection vertex back
                    ++n;
                    x[i+2] = x[i+1];
                    y[i+2] = y[i+1];
                    
                    // add new vertex at (+x, -y) corner
                    x[i+1] = root->box.b[Dx];
                    y[i+1] = root->box.a[Dy];
                }
            }
            else {
                // push last intersection vertex back
                ++n;
                x[i+2] = x[i+1];
                y[i+2] = y[i+1];
                
                // add new vertex at (-x, -y) corner
                x[i+1] = root->box.a[Dx];
                y[i+1] = root->box.a[Dy];
            }
        }
    }
    
    findInPolygon(&out, x, y, n, Dx, Dy, root);
    
    return out;
}

// Return the nearest k neighbors of the given volume, in the order
// of proximity (nearest first).
std::vector<int> RTree::nearestNeighbors ( const RTree::Box & box,
    unsigned int k ) const
{
    unsigned int i;
    std::vector<int> out;
    if ( root == NULL ) return out;
    
    out.reserve(k);
    
    PageVal p;
    std::set<PageVal> pages;
    pages.insert(PageVal(root, root->box.sqdist(box)));
    
    while ( pages.size() > 0 && out.size() < k ) {
        p = *pages.begin();
        pages.erase(pages.begin());
        if ( p.page->id >= 0 ) out.push_back(p.page->id);
        else {
            for ( i = 0; i < p.page->n; ++i ) {
                pages.insert(
                    PageVal(p.page->children[i],
                    p.page->children[i]->box.sqdist(box)));
            }
        }
    }
    
    return out;
}

// Insertion / Removal //////////////////////////////////////////////

// choose which subtree to insert into.
unsigned int RTree::insertChooseSubtree ( const RTree::Box& box,
    RTree::Page const * const p ) const
{
    if ( p->n < p->nC ) return p->n;
    
    unsigned int k = 0;
    
    Box temp = p->box;
    temp.expand(box);
    
    double pvol = p->box.volume();      // volume of p's box
    double lap;                         // overlap
    double vol;                         // volume change
    int    dn = p->size;                // optimal subtree size
    double dlap = 0.0;                  // optimal overlap
    double dvol = temp.volume() - pvol; // optimal volume change
    
    for ( unsigned int i = 0; i < p->n; ++i ) {
        // prioritize smallest volume change (minimize overlap
        //   between other boxes)
        // in tie, favor maximum overlap (minimize overlap between
        //   other boxes)
        // in tie, favor smallest subtree size (favor balanced
        //   tree)
        temp = p->children[i]->box;
        temp.expand(box);
        
        lap = p->children[i]->box.overlap(box);
        vol = temp.volume() - p->box.volume();
        
        if ( vol < dvol ) {
            k = i;
            dvol = vol;
            dlap = lap;
            dn = p->children[i]->size;
        }
        else if ( vol == dvol ) {
            if ( lap > dlap ) {
                k = i;
                dlap = lap;
                dn = p->children[i]->size;
            }
            else if ( lap == dlap && p->children[i]->size < dn ) {
                k = i;
                dn = p->children[i]->size;
            }
        }
    }
    
    return k;
}

// Insert volume; return 1 if page's box must be updated.
RTree::Page* RTree::insert ( int id, const RTree::Box& box,
    RTree::Page* p )
{
    // new leaf node
    if ( p == NULL ) p = new Page(nC, id, box);
    
    // leaf node
    else if ( p->n == 0 ) {
        p->children[0] = new Page(nC, p->id, p->box);
        p->children[1] = insert(id, box, p->children[1]);
        p->id = -1;
        p->n = 2;
    }
    
    // subtree
    else {
        unsigned int i = insertChooseSubtree(box, p);
        p->children[i] = insert(id, box, p->children[i]);
        if ( i == p->n ) ++p->n;
    }
    
    p->update();
    return p;
}

// Insert volume with id.
void RTree::insert ( int id, const RTree::Box & box )
{
    root = insert(id, box, root);
}

RTree::Page* RTree::remove ( int id, const RTree::Box& box,
    RTree::Page* p, std::set<RTree::Page>* orphaned )
{
    // no overlap
    if ( p == NULL || p->box.overlap(box) == 0.0 )
        return p;
    
    // remove leaf
    if ( p->id == id ) {
        delete p;
        return NULL;
    }
    
    // remove from children
    bool flag = false;
    unsigned int i, j;
    for ( i = 0; i < p->n; ) {
        p->children[i] = remove(id, box, p->children[i], orphaned);
        if ( p->children[i] == NULL ) {
            for ( j = i+1; j < p->n; ++j )
                p->children[j-1] = p->children[j];
            p->children[p->n - 1] = NULL;
            --p->n;
            flag = true;
        }
        else ++i;
    }
    
    // remove if now empty
    if ( flag && p->n == 0 ) {
        delete p;
        return NULL;
    }
    
    p->update();
    return p;
}

// Remove intersecting volumes with id.
void RTree::remove ( int id, const RTree::Box & box )
{
    std::set<Page> orphaned;
    root = remove(id, box, root, &orphaned);
    
    if ( root == NULL ) return;
    
    // reinsert nodes from dissolved subtree
    while ( orphaned.size() > 0 ) {
        insert((*orphaned.begin()).id, (*orphaned.begin()).box);
        orphaned.erase(orphaned.begin());
    }
}

// Updates the volume with id.
void RTree::update ( int id, const RTree::Box & oldBox,
    const RTree::Box & newBox )
{
    remove(id, oldBox);
    insert(id, newBox);
}

// I/O //////////////////////////////////////////////////////////////

// Print out the tree.
void RTree::printTree ( RTree::Page const * const p, unsigned int L )
    const
{
    if ( p == NULL ) return;
    unsigned int i;
    for ( i = 0; i < L; ++i ) std::cout << "  ";
    if ( p->id >= 0 ) std::cout << p->id << " ";
    else std::cout << ". ";
    std::cout << "[" << p->box.a[0] << " " << p->box.b[0];
    for ( i = 1; i < nD; ++i )
        std::cout << "; " << p->box.a[i] << " " << p->box.b[i];
    std::cout << "] N=" << p->size << "\n";
    for ( i = 0; i < p->n; ++i ) printTree(p->children[i], L+1);
}

// print out the tree
void RTree::printTree () const
{
    printTree(root, 0);
}

#endif // RTREE_H