/////////////////////////////////////////////////////////////////////
// Vector for a field in a data table.                             //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// 06/29/2020 - Brennan Young                                      //
// - created.                                                      //
// 07/02/2020 - Brennan Young                                      //
// - added units data member.                                      //
// 11/03/2020 - Brennan Young                                      //
// - moved size() to 'getters.                                     //
// - added number of elements to constructor.                      //
// 11/10/2020 - Brennan Young                                      //
// - copy constructor and assignment operator now copy full length //
//   of data array (replaced v.size() with n*dlen).                //
// 12/09/2020 - Brennan Young                                      //
// - update to use new Statistics object.                          //
// - added calculateStatistics.                                    //
// - added operator[], as an alternative to getf().                //
// - added nul member.                                             //
// - added setNull.                                                //
// - added isNull.                                                 //
// - removed changeType.                                           //
// 02/22/2021 - Brennan Young                                      //
// - constructor now sets nul.                                     //
// - overloaded setNull to assume no change to current data.       //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_DATAVECTOR_20200629
#define YOUNG_DATAVECTOR_20200629

#include <cstdlib>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <string>

#include "Statistics.h"


class DataVector {
public:
    // constants
    static const unsigned char T_INT   = 0; // int
    static const unsigned char T_UINT  = 1; // unsigned int
    static const unsigned char T_FLOAT = 2; // float
    static const unsigned char T_STR   = 3; // string
private:
    unsigned char dtyp; // type of data
    unsigned char dlen; // element byte length
    size_t n;           // number of elements
    char* data;         // data array
    double nul;         // no-data value
    Statistics statistics;
    
    void toBytes(char*, char*, int) const;
    void fromBytes(char*, char*, int) const;
public:
    std::string name;
    std::string units;
    
    // constructors / destructor
    DataVector(unsigned char, unsigned char, size_t,
        const std::string&, const std::string&);
    DataVector(const DataVector&);
    ~DataVector();
    
    // operators
    DataVector& operator=(const DataVector&);
    double operator[](size_t) const;
    
    // getters
    const Statistics& stats() const;
    size_t size() const;
    unsigned char getType() const;
    unsigned char getLen() const;
    template<class T> void get(size_t, T*) const;
    int geti(size_t) const;
    double getf(size_t) const;
    std::string getStr(size_t) const;
    double getNull() const;
    
    // setters
    template<class T> void set(size_t, const T&);
    void setStr(size_t, const std::string&);
    void setNull(double, bool);
    void setNull(double);
    
    // operations
    void resize(size_t);
    void zero();
    bool isNull(size_t) const;
    bool isNull(double) const;
    void calculateStatistics();
}; // DataVector


// CONSTRUCTORS / DESTRUCTORS ///////////////////////////////////////

// Default constructor.
DataVector::DataVector ( unsigned char t=T_INT, unsigned char s=4,
    size_t nr=1, const std::string& nm="field",
    const std::string& un="" )
: dtyp(t), dlen(s), n(nr), nul(-9999.0), name(nm), units(un)
{
    if ( dtyp != T_INT && dtyp != T_UINT && dtyp != T_FLOAT
            && dtyp != T_STR )
        dtyp = T_INT;
    
    try {
        data = new char [n*dlen];
    }
    catch ( ... ) {
        std::cout << "Memory allocation error.\n";
        exit(0);
    }
}

// Copy constructor.
DataVector::DataVector ( const DataVector& v )
: dtyp(v.dtyp), dlen(v.dlen), n(v.n), name(v.name), units(v.units)
{
    try {
        data = new char [n*dlen];
    }
    catch ( ... ) {
        std::cout << "Memory allocation error.\n";
        exit(0);
    }
    for ( size_t i = 0; i < n*dlen; ++i ) data[i] = v.data[i];
    
    nul = v.nul;
    statistics = v.statistics;
}

// Destructor.
DataVector::~DataVector ()
{
    delete [] data;
}


// OPERATORS ////////////////////////////////////////////////////////

// Assignment.
DataVector& DataVector::operator= ( const DataVector& v )
{
    if ( this == &v ) return *this;
    
    dtyp = v.dtyp;
    dlen = v.dlen;
    n = v.n;
    
    delete [] data;
    try {
        data = new char [n*dlen];
    }
    catch ( ... ) {
        std::cout << "Memory allocation error.\n";
        exit(0);
    }
    for ( size_t i = 0; i < n*dlen; ++i ) data[i] = v.data[i];
    
    nul = v.nul;
    statistics = v.statistics;
    name = v.name;
    units = v.units;
    
    return *this;
}

// Get element as floating-point number.
double DataVector::operator[] ( size_t i ) const
{
    return getf(i);
}


// HELPERS //////////////////////////////////////////////////////////

void DataVector::toBytes ( char* bytes, char* x, int len ) const
{
    for ( int i = 0; i < len; ++i ) bytes[i] = x[i];
}

void DataVector::fromBytes ( char* bytes, char* x, int len ) const
{
    for ( int i = 0; i < len; ++i ) x[i] = bytes[i];
}


// GETTERS //////////////////////////////////////////////////////////

// Return the object's statistics object.
const Statistics& DataVector::stats () const
{
    return statistics;
}

// Returns the number of elements in the vector.
size_t DataVector::size () const { return n; }

// Returns the data type of the vector.
unsigned char DataVector::getType () const { return dtyp; }

// Returns the byte length of the vector's data type.
unsigned char DataVector::getLen () const { return dlen; }

// Interpret the value of an element based on the vector's type.
template<class T> void DataVector::get ( size_t i, T* x ) const
{
    if ( dtyp == T_INT ) {
        if ( dlen == sizeof(char) ) {
            char a;
            fromBytes(&data[i*dlen], &a, dlen);
            *x = (T) a;
        }
        else if ( dlen == sizeof(short) ) {
            short a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            *x = (T) a;
        }
        else if ( dlen == sizeof(int) ) {
            int a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            *x = (T) a;
        }
        else if ( dlen == sizeof(long) ) {
            long a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            *x = (T) a;
        }
        else if ( dlen == sizeof(long long) ) {
            long long a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            *x = (T) a;
        }
        else {
            int L = sizeof(T) < dlen ? sizeof(T) : dlen;
            fromBytes(&data[i*dlen], (char*) x, L);
        }
    }
    else if ( dtyp == T_UINT ) {
        if ( dlen == sizeof(unsigned char) ) {
            unsigned char a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            *x = (T) a;
        }
        else if ( dlen == sizeof(unsigned short) ) {
            unsigned short a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            *x = (T) a;
        }
        else if ( dlen == sizeof(unsigned int) ) {
            unsigned int a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            *x = (T) a;
        }
        else if ( dlen == sizeof(unsigned long) ) {
            unsigned long a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            *x = (T) a;
        }
        else if ( dlen == sizeof(unsigned long long) ) {
            unsigned long long a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            *x = (T) a;
        }
        else {
            int L = sizeof(T) < dlen ? sizeof(T) : dlen;
            fromBytes(&data[i*dlen], (char*) x, L);
        }
    }
    else if ( dtyp == T_FLOAT ) {
        if ( dlen == sizeof(float) ) {
            float a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            *x = a;
        }
        else if ( dlen == sizeof(double) ) {
            double a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            *x = a;
        }
        else {
            int L = sizeof(T) < dlen ? sizeof(T) : dlen;
            fromBytes(&data[i*dlen], (char*) x, L);
        }
    }
    else *x = 0;
}

// Return the value of the element, interpreted as a signed integer.
int DataVector::geti ( size_t i ) const
{
    int x;
    get(i, &x);
    return x;
}

// Return the value of the element, interpreted as a double.
double DataVector::getf ( size_t i ) const
{
    double x;
    get(i, &x);
    return x;
}

// Return the value of an element, interpreted as a string. If the
// vector's data type is a number, returns a string containing that
// number.
std::string DataVector::getStr ( size_t i ) const
{
    std::stringstream ss;
    
    if ( dtyp == T_INT ) {
        if ( dlen == sizeof(char) ) {
            char a;
            fromBytes(&data[i*dlen], &a, dlen);
            ss << a;
        }
        else if ( dlen == sizeof(short) ) {
            short a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            ss << a;
        }
        else if ( dlen == sizeof(int) ) {
            int a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            ss << a;
        }
        else if ( dlen == sizeof(long) ) {
            long a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            ss << a;
        }
        else if ( dlen == sizeof(long long) ) {
            long long a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            ss << a;
        }
    }
    else if ( dtyp == T_UINT ) {
        if ( dlen == sizeof(unsigned char) ) {
            unsigned char a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            ss << a;
        }
        else if ( dlen == sizeof(unsigned short) ) {
            unsigned short a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            ss << a;
        }
        else if ( dlen == sizeof(unsigned int) ) {
            unsigned int a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            ss << a;
        }
        else if ( dlen == sizeof(unsigned long) ) {
            unsigned long a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            ss << a;
        }
        else if ( dlen == sizeof(unsigned long long) ) {
            unsigned long long a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            ss << a;
        }
    }
    else if ( dtyp == T_FLOAT ) {
        if ( dlen == sizeof(float) ) {
            float a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            ss << a;
        }
        else if ( dlen == sizeof(double) ) {
            double a;
            fromBytes(&data[i*dlen], (char*) &a, dlen);
            ss << a;
        }
    }
    else if ( dtyp == T_STR ) {
        for ( size_t j = 0; j < dlen; ++j ) {
            if ( data[i*dlen+j]== '\0' ) break;
            ss << data[i*dlen+j];
        }
    }
    
    return ss.str();
}


// SETTERS //////////////////////////////////////////////////////////

template<class T> void DataVector::set ( size_t i, const T& x )
{
    if ( dtyp == T_INT ) {
        if ( dlen == sizeof(char) ) {
            char a = (char) x;
            toBytes(&data[i*dlen], &a, dlen);
        }
        else if ( dlen == sizeof(short) ) {
            short a = (short) x;
            toBytes(&data[i*dlen], (char*) &a, dlen);
        }
        else if ( dlen == sizeof(int) ) {
            int a = (int) x;
            toBytes(&data[i*dlen], (char*) &a, dlen);
        }
        else if ( dlen == sizeof(long) ) {
            long a = (long) x;
            toBytes(&data[i*dlen], (char*) &a, dlen);
        }
        else if ( dlen == sizeof(long long) ) {
            long long a = (long long) x;
            toBytes(&data[i*dlen], (char*) &a, dlen);
        }
        else {
            T a = x;
            int L = sizeof(T) < dlen ? sizeof(T) : dlen;
            toBytes(&data[i*dlen], (char*) &a, L);
        }
    }
    else if ( dtyp == T_UINT ) {
        if ( dlen == sizeof(unsigned char) ) {
            unsigned char a = (unsigned char) x;
            toBytes(&data[i*dlen], (char*) &a, dlen);
        }
        else if ( dlen == sizeof(unsigned short) ) {
            unsigned short a = (unsigned short) x;
            toBytes(&data[i*dlen], (char*) &a, dlen);
        }
        else if ( dlen == sizeof(unsigned int) ) {
            unsigned int a = (unsigned int) x;
            toBytes(&data[i*dlen], (char*) &a, dlen);
        }
        else if ( dlen == sizeof(unsigned long) ) {
            unsigned long a = (unsigned long) x;
            toBytes(&data[i*dlen], (char*) &a, dlen);
        }
        else if ( dlen == sizeof(unsigned long long) ) {
            unsigned long long a = (unsigned long long) x;
            toBytes(&data[i*dlen], (char*) &a, dlen);
        }
        else {
            T a = x;
            int L = sizeof(T) < dlen ? sizeof(T) : dlen;
            toBytes(&data[i*dlen], (char*) &a, L);
        }
    }
    else if ( dtyp == T_FLOAT ) {
        if ( dlen == sizeof(float) ) {
            float a = x;
            toBytes(&data[i*dlen], (char*) &a, dlen);
        }
        else if ( dlen == sizeof(double) ) {
            double a = x;
            toBytes(&data[i*dlen], (char*) &a, dlen);
        }
        else {
            T a = x;
            int L = sizeof(T) < dlen ? sizeof(T) : dlen;
            toBytes(&data[i*dlen], (char*) &a, L);
        }
    }
}

void DataVector::setStr ( size_t i, const std::string& x )
{
    size_t j = 0;
    for ( ; j < dlen && j < x.size(); ++j )
        data[i*dlen+j] = x[j];
    if ( j < dlen ) data[i*dlen+j] = '\0';
}

// Set the no-data value, and change all current no-data values to
// the new one if change is true.
void DataVector::setNull ( double x, bool change )
{
    if ( x == nul ) return;
    if ( change ) {
        for ( size_t i = 0; i < n; ++i )
            if ( getf(i) == nul ) set(i, x);
    }
    nul = x;
}
void DataVector::setNull ( double x )
{
    setNull(x, false);
}


// OPERATIONS ///////////////////////////////////////////////////////

// Resize the array.
void DataVector::resize ( size_t nr ) {
    if ( n == nr ) return;
    char* x = new char [nr*dlen];
    for ( size_t i = 0; i < n*dlen && i < nr*dlen; ++ i )
        x[i] = data[i];
    delete [] data;
    data = x;
    n = nr;
}

// Set every element of the vector to zero. If string, sets it to
// the null character.
void DataVector::zero ()
{
    if ( dtyp == T_STR )
        for ( size_t i = 0; i < n*dlen; ++i ) data[i] = '\0';
    else
        for ( size_t i = 0; i < n*dlen; ++i ) data[i] = 0;
}

// Get the no-data value
double DataVector::getNull () const
{
    return nul;
}

// Determine if the given element is the same as the no-data value.
bool DataVector::isNull ( size_t i ) const
{
    return getf(i) == nul;
}

bool DataVector::isNull ( double x ) const
{
    return x == nul;
}

// Compute the vector's statistics.
void DataVector::calculateStatistics ()
{
    statistics = Statistics (*this, 0, this->size(), 1, nul);
}


#endif // YOUNG_DATAVECTOR_20200629