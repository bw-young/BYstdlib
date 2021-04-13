/////////////////////////////////////////////////////////////////////
// Constraint                                                      //
// Object for managing numerical constraints.                      //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- History ---------------------------------------------------- //
// 04/02/2021 - Brennan Young                                      //
// - created.                                                      //
// 04/07/2021 - Brennan Young                                      //
// - added allData static method.                                  //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_CONSTRAINT_20210402
#define YOUNG_CONSTRAINT_20210402

#include <list>


namespace bystd {


class Constraint {
public:
    double min;     // minimum of range
    double max;     // maximum of range
    bool inclusive; // true to include max
    bool inverse;   // true for NOT in range
    
    // constructors, destructors
    Constraint(double, double, bool, bool);
    Constraint(const Constraint&);
    ~Constraint();
    
    // operators
    Constraint& operator=(const Constraint&); // copy
    bool operator()(double) const; // check value against constraint
    
    // basic constraints
    static Constraint allData(){ return Constraint(0,0,false,true); }
};


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

// Default constructor.
Constraint::Constraint ( double a=0.0, double b=1.0,
    bool inc=false, bool inv=false )
: min(a), max(b), inclusive(inc), inverse(inv)
{}

// Copy constructor.
Constraint::Constraint ( const Constraint& C )
: min(C.min), max(C.max), inclusive(C.inclusive), inverse(C.inverse)
{}

// Destructor.
Constraint::~Constraint () {}


// OPERATORS ////////////////////////////////////////////////////////

// Assignment.
Constraint& Constraint::operator= ( const Constraint& C )
{
    min = C.min;
    max = C.max;
    inclusive = C.inclusive;
    inverse = C.inverse;
    
    return *this;
}

// Check value against constraint - return true if within constraint.
bool Constraint::operator() ( double value ) const
{
    bool x = value >= min
        && ((value < max) || (inclusive && value <= max));
    return (!inverse && x) || (inverse && !x);
}


}; // namespace bystd

#endif // YOUNG_CONSTRAINT_20210402