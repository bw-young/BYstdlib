/////////////////////////////////////////////////////////////////////
// BYstdlib                                                        //
// CardinalDir                                                     //
//                                                                 //
// Object for managing working with the cardinal directions (N, E, //
// S, W) and their diagonals (NE, SE, SW, NW).                     //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/15/2020                                                      //
// - created in GridObjPolyOps.h.                                  //
// 03/08/2021                                                      //
// - migrated to this file from GridObjPolyOps.h.                  //
// - added to bystd namespace.                                     //
/////////////////////////////////////////////////////////////////////

namespace bystd {

// Cardinal direction object, for helping to track 8-mode direction.
class CardinalDir {
private:
    static const unsigned char N  = 0;
    static const unsigned char NE = 1;
    static const unsigned char E  = 2;
    static const unsigned char SE = 3;
    static const unsigned char S  = 4;
    static const unsigned char SW = 5;
    static const unsigned char W  = 6;
    static const unsigned char NW = 7;
    
    unsigned char dir;
    char dx;
    char dy;
    
    void setdxdy ()
    {
        if      ( dir == NE ) { dx= 1; dy= 1; }
        else if ( dir == E  ) { dx= 1; dy= 0; }
        else if ( dir == SE ) { dx= 1; dy=-1; }
        else if ( dir == S  ) { dx= 0; dy=-1; }
        else if ( dir == SW ) { dx=-1; dy=-1; }
        else if ( dir == W  ) { dx=-1; dy= 0; }
        else if ( dir == NW ) { dx=-1; dy= 1; }
        else                  { dx= 0; dy= 1; }
    }
public:
    // default constructor
    CardinalDir () : dir (0) { setdxdy(); }
    
    // construct with string indicating the direction
    CardinalDir ( const std::string& s ) : dir(0)
    {
        if      ( s == "N"  ) dir = N;
        else if ( s == "NE" ) dir = NE;
        else if ( s == "E"  ) dir = E;
        else if ( s == "SE" ) dir = SE;
        else if ( s == "S"  ) dir = S;
        else if ( s == "SW" ) dir = SW;
        else if ( s == "W"  ) dir = W;
        else if ( s == "NW" ) dir = NW;
        setdxdy();
    }
    
    // construct with difference in x and y indicating the direction
    CardinalDir ( float x, float y )
    {
        if ( y < 0.0f ) {
            if      ( x < 0.0f ) dir = SW;
            else if ( x > 0.0f ) dir = SE;
            else                 dir = S;
        }
        else if ( y > 0.0f ) {
            if      ( x < 0.0f ) dir = NW;
            else if ( x > 0.0f ) dir = NE;
            else                 dir = N;
        }
        else {
            if ( x < 0.0f )      dir = W;
            else if ( x > 0.0f ) dir = E;
            else                 dir = N; // no direction
        }
        setdxdy();
    }
    
    // copy constructor
    CardinalDir ( const CardinalDir& d ) :
        dir(d.dir), dx(d.dx), dy(d.dy) {}
    
    // destructor
    ~CardinalDir () {}
    
    // assignment
    CardinalDir& operator= ( const CardinalDir& d )
    {
        dir = d.dir;
        dx  = d.dx;
        dy  = d.dy;
        return *this;
    }
    
    // tell if same as another direction
    bool operator== ( const CardinalDir& d ) const
    {
        return dir == d.dir;
    }
    bool operator!= ( const CardinalDir& d ) const
    {
        return dir != d.dir;
    }
    
    // turn clockwise
    CardinalDir& operator++ ()
    {
        if ( dir == 7 ) dir = 0;
        else  ++dir;
        setdxdy();
        return *this;
    }
    
    // turn counter-clockwise
    CardinalDir& operator-- ()
    {
        if ( dir == 0 ) dir = 7;
        else --dir;
        setdxdy();
        return *this;
    }
    
    // turn clockwise 45 degrees x times.
    CardinalDir operator+ ( int x ) const
    {
        CardinalDir out;
        int d = ((int) dir + x) % 8;
        if ( d < 0 ) d += 8;
        out.dir = (unsigned char) d;
        out.setdxdy();
        return out;
    }
    
    // turn counter-clockwise 45 degrees x times
    CardinalDir operator- ( int x ) const
    {
        CardinalDir out;
        int d = ((int) dir - x) % 8;
        if ( d < 0 ) d += 8;
        out.dir = (unsigned char) d;
        out.setdxdy();
        return out;
    }
    
    // get CCW difference to turn from this to d; negative if CW.
    int operator- ( const CardinalDir& d ) const
    {
        CardinalDir out;
        int x = (int) dir - (int) d.dir;
        if ( x < -3 ) return x+8;
        if ( x > 4 ) return x-8;
        return x;
    }
    
    int x () const { return (int) dx; }
    int y () const { return (int) dy; }
}; // CardinalDir

}; // namespace bystd