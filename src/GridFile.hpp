/////////////////////////////////////////////////////////////////////
// Grid.                                                           //
// Class for creating objects to directly query and manipulate     //
// raster data on the disk (not held in memory).                   //
/////////////////////////////////////////////////////////////////////

#include <cmath>    // nan
#include <cstdlib>  // exit
#include <iostream> // std::cout
#include <fstream>  // std::fstream
#include <string>   // std::string

class Grid {
public:
    // data type constants
    const unsigned char INT = 0;
    const unsigned char FLOAT = 1;
    
    // interleave format constants
    const unsigned char BSQ = 0;
    const unsigned char BIL = 1;
    const unsigned char BIP = 2;
    
    // helper function
    void byteswap(void*, int) const;
private:
    unsigned char dataType; // how to interpret file binary
    unsigned char dataBytes; // byte-length for each element
    unsigned char format; // interleave format
    std::string filename; // file name, with full path and extension
    std::fstream file; // the active file to read/write data
    bool swap; // true if byte-swap is necessary
    int nr; // number of rows
    int nc; // number of columns
    int nb; // number of bands
public:
    // constructors, destructor
    Grid(const std::string&);
    Grid(unsigned char, unsigned char, int, int, int);
    Grid(int, int, int);
    Grid(const Grid&);
    ~Grid();
    
    // operators
    Grid& operator=(const Grid&);
    
    // file access
    void open(const std::string&);
    void close();
    
    // structure getters
    int nrows() const;
    int ncols() const;
    int nbands() const;
    int size() const;
    int volume() const;
    int index(int, int) const;
    int index(int, int, int) const;
    
    // structure setters
    void setByteSwap(bool);
    void setDimensions(int, int, int);
    void setInterleave(unsigned char);
    void setDataType(unsigned char, unsigned char);
    
    // element access
    double get(int);
    double get(int, int, int);
    double get(int, int);
    template <typename T> void set(T, int);
    template <typename T> void set(T, int, int, int);
    template <typename T> void set(T, int, int);
}; // Grid


// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

Grid::Grid ( const std::string& fn="" )
: dataType(FLOAT), dataBytes(4), format(BSQ),
    swap(false), nr(0), nc(0), nb(0)
{}

Grid::Grid ( unsigned char newDataType, unsigned char byteLength,
    int rows, int cols, int bands )
: dataType(FLOAT), dataBytes(4), format(BSQ),
    swap(false), nr(rows), nc(cols), nb(bands)
{
    setDataType(newDataType, byteLength);
}

Grid::Grid ( int rows, int cols, int bands )
: dataType(FLOAT), dataBytes(4), format(BSQ),
    swap(false), nr(rows), nc(cols), nb(bands)
{}

Grid::Grid ( const Grid& g )
: dataType(g.dataType), dataBytes(g.dataBytes), format(g.format),
    swap(false), nr(g.nr), nc(g.nc), nb(g.nb)
{
    if ( g.file.is_open() ) open(g.filename);
    else filename = g.filename;
}

Grid::~Grid () { if ( file.is_open() ) file.close(); }


// OPERATORS ////////////////////////////////////////////////////////

Grid& Grid::operator= ( const Grid& g )
{
    if ( &g == this ) return *this;
    
    if ( file.is_open() ) file.close();
    
    dataType = g.dataType;
    dataBytes = g.dataBytes;
    format = g.format;
    nr = g.nr;
    nc = g.nc;
    nb = g.nb;
    if ( g.file.is_open() ) open(g.filename);
    else filename = g.filename;
    
    return *this;
}


// FILE ACCESS //////////////////////////////////////////////////////

void Grid::open ( const std::string& fn )
{
    if ( file.is_open() ) file.close();
    
    filename = fn;
    
    try {
        file = std::fstream(filename,
            std::ios::in | std::ios::out | std::ios::binary);
        if ( !file.is_open() ) {
            // try to create the file
            file = std::fstream(filename,
                std::ios::out | std::ios::binary);
            if ( !file.is_open() ) throw 1;
            else {
                int n = nr * nc * nb;
                double buff = nan("");
                for ( int i = 0; i < n; ++i ) {
                    file.write((char*)& buff, dataBytes);
                }
                file.close();
                file = std::fstream(filename,
                    std::ios::in | std::ios::out | std::ios::binary);
            }
        }
    }
    catch (...) {
        std::cout << "Error: could not open file '"
                  << filename << "'\n";
        exit(0);
    }
}

void Grid::close () { if ( file.is_open() ) file.close(); }


// GETTERS FOR NAVIGATING FILE STRUCTURE ////////////////////////////

int Grid::nrows () const { return nr; }
int Grid::ncols () const { return nc; }
int Grid::nbands () const { return nb; }
int Grid::size () const { return nr * nc; }
int Grid::volume () const { return nr * nc * nb; }
int Grid::index ( int i, int j, int k ) const
{
    if ( format == BIL ) return i*(nc*nb) + k*nb + j;
    if ( format == BIP ) return i*(nc*nb) + j*nb + k;
    return k*(nr*nc) + i*nc + j;
}
int Grid::index ( int i, int j ) const { return index(i,j,0); }


// SETTERS FOR NAVIGATING FILE STRUCTURE ////////////////////////////

void Grid::setByteSwap ( bool byteSwap )
{
    swap = byteSwap;
}

void Grid::setDimensions ( int rows, int cols, int bands )
{
    nr = rows;
    nc = cols;
    nb = bands;
}

void Grid::setInterleave ( unsigned char newFormat )
{
    if ( newFormat == BSQ || newFormat == BIL || newFormat == BIP )
        format = newFormat;
}

void Grid::setDataType ( unsigned char newDataType,
    unsigned char newByteLength )
{
    if ( newDataType == INT | newDataType == FLOAT )
    {
        dataType = newDataType;
        dataBytes = newByteLength;
    }
}


// ELEMENT ACCESS ///////////////////////////////////////////////////

void Grid::byteswap ( void* data, int n ) const
{
    unsigned char* start = (unsigned char*) data;
    unsigned char* end = start + n - 1;
    unsigned char swap;
    for ( ; start < end; ++start, --end ){
        swap = *start;
        *start = *end;
        *end = swap;
    }
}

double Grid::get ( int i )
{
    if ( !file.is_open() || i < 0 || i > size() ) return nan("");
    file.seekg(i * dataBytes);
    
    try {
        if ( dataType == INT ) {
            if ( dataBytes == 1 ) {
                char buff;
                file.read(&buff, dataBytes);
                if ( swap ) byteswap(&buff, dataBytes);
                return buff;
            }
            else if ( dataBytes == 2 ) {
                short buff;
                file.read((char*) &buff, dataBytes);
                if ( swap ) byteswap(&buff, dataBytes);
                return buff;
            }
            else if ( dataBytes == 4 ) {
                int buff;
                file.read((char*) &buff, dataBytes);
                if ( swap ) byteswap(&buff, dataBytes);
                return buff;
            }
            else if ( dataBytes == 8 ) {
                long long buff;
                file.read((char*) &buff, dataBytes);
                if ( swap ) byteswap(&buff, dataBytes);
                return buff;
            }
        }
        else {
            if ( dataBytes == 4 ) {
                float buff;
                file.read((char*) &buff, dataBytes);
                if ( swap ) byteswap(&buff, dataBytes);
                return buff;
            }
            else if ( dataBytes == 8 ) {
                double buff;
                file.read((char*) &buff, dataBytes);
                if ( swap ) byteswap(&buff, dataBytes);
                return buff;
            }
        }
    }
    catch (...) {
        std::cout << "Error: unable to read from '" << filename << "'";
        exit(0);
    }
    
    return nan("");
}

double Grid::get ( int i, int j, int k )
{
    if ( i < 0 || i >= nr || j < 0 || j >= nc || k < 0 || k >= nb )
        return nan("");
    return get(index(i,j,k));
}

double Grid::get ( int i, int j ) { return get(i,j,0); }

template <typename T> void Grid::set ( T x, int i )
{
    if ( !file.is_open() || i < 0 || i > size() ) return;
    file.seekp(i * dataBytes);
    
    try {
        if ( dataType == INT ) {
            if ( dataBytes == 1 ) {
                char buff = x;
                if ( swap ) byteswap(&buff, dataBytes);
                file.write(&buff, dataBytes);
            }
            else if ( dataBytes == 2 ) {
                short buff = x;
                if ( swap ) byteswap(&buff, dataBytes);
                file.write((char*) &buff, dataBytes);
            }
            else if ( dataBytes == 4 ) {
                int buff = x;
                if ( swap ) byteswap(&buff, dataBytes);
                file.write((char*) &buff, dataBytes);
            }
            else if ( dataBytes == 8 ) {
                long long buff = x;
                if ( swap ) byteswap(&buff, dataBytes);
                file.write((char*) &buff, dataBytes);
            }
        }
        else {
            if ( dataBytes == 4 ) {
                float buff = x;
                if ( swap ) byteswap(&buff, dataBytes);
                file.write((char*) &buff, dataBytes);
            }
            else if ( dataBytes == 8 ) {
                double buff = x;
                if ( swap ) byteswap(&buff, dataBytes);
                file.write((char*) &buff, dataBytes);
            }
        }
    }
    catch (...) {
        std::cout << "Error: unable to write to '" << filename << "'";
        exit(0);
    }
}

template <typename T> void Grid::set ( T x, int i, int j, int k )
{
    if ( i < 0 || i >= nr || j < 0 || j >= nc || k < 0 || k >= nb )
        return;
    set(x, index(i,j,k));
}

template <typename T> void Grid::set ( T x, int i, int j )
{
    set(x, i,j,0);
}