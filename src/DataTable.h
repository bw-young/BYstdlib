/////////////////////////////////////////////////////////////////////
// Data table for managing multiple fields of diverse data types.  //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 06/29/2020 - Brennan Young                                      //
// - created.                                                      //
// 07/02/2020 - Brennan Young                                      //
// - overloaded addField to permit omition of insertion position.  //
// 08/17/2020 - Brennan Young                                      //
// - added findField.                                              //
// 11/03/2020 - Brennan Young                                      //
// - pass number of rows when constructing new DataVector objects. //
//   in addField.                                                  //
// 12/09/2020 - Brennan Young                                      //
// - update to use Statistics object.                              //
// 02/22/2021 - Brennan Young                                      //
// - added setNull, getNull, and isNull methods.                   //
// 03/10/2021 - Brennan Young                                      //
// - overloaded geti, getf, and getStr to take the field name.     //
/////////////////////////////////////////////////////////////////////

#ifndef YOUNG_DATATABLE_20200629
#define YOUNG_DATATABLE_20200629

#include <vector>
#include "DataVector.h"


class DataTable {
private:
    std::vector<DataVector> fields;
    size_t n; // number of records
public:
    // constructors / destructor
    DataTable();
    DataTable(const DataTable&);
    ~DataTable();
    
    // operators
    DataTable& operator=(const DataTable&);
    
    // getters
    const Statistics& stats(size_t) const;
    size_t numFields() const;
    size_t numRecords() const;
    size_t findField(const std::string&) const;
    std::string& fieldName(size_t);
    std::string fieldName(size_t) const;
    unsigned char fieldType(size_t) const;
    unsigned char fieldLength(size_t) const;
    template<typename T> void get(size_t, size_t, T*) const;
    int geti(size_t, size_t) const;
    int geti(size_t, const std::string&) const;
    double getf(size_t, size_t) const;
    double getf(size_t, const std::string&) const;
    std::string getStr(size_t, size_t) const;
    std::string getStr(size_t, const std::string&) const;
    double getNull(size_t) const;
    
    // setters
    template<typename T> void set(size_t, size_t, const T&);
    void setStr(size_t, size_t, const std::string&);
    void setNull(size_t, double);
    
    // operations
    void addField(size_t, unsigned char, unsigned char,
        const std::string&, const std::string&);
    void addField(size_t, unsigned char, unsigned char,
        const std::string&);
    void addField( unsigned char, unsigned char, const std::string&,
        const std::string&);
    void addField( unsigned char, unsigned char, const std::string&);
    void removeField(size_t);
    void removeField(const std::string&);
    void setNumRecords(size_t);
    void zero();
    bool isNull(size_t, size_t) const;
    bool isNull(double, size_t) const;
}; // DataTable


// CONSTRUCTORS / DESTRUCTORS ///////////////////////////////////////

// Default constructor.
DataTable::DataTable () : n(0) {}

// Copy constructor.
DataTable::DataTable ( const DataTable& t ) : n(t.n)
{
    fields = t.fields;
}

// Destructor.
DataTable::~DataTable () {}


// OPERATORS ////////////////////////////////////////////////////////

// Assignment.
DataTable& DataTable::operator= ( const DataTable& t )
{
    if ( this == &t ) return *this;
    
    n = t.n;
    fields = t.fields;
    
    return *this;
}


// GETTERS //////////////////////////////////////////////////////////

// Get the statistics for field j.
const Statistics& DataTable::stats ( size_t j ) const
{
    return fields[j].stats();
}

// Get the number of fields in the table.
size_t DataTable::numFields() const { return fields.size(); }

// Get the number of records in the table (number of elements in
// each field).
size_t DataTable::numRecords() const { return n; }

// Get the index of field with name. If not found, returns the size
// of the fields container.
size_t DataTable::findField ( const std::string& name ) const
{
    for ( size_t i = 0; i < fields.size(); ++i )
        if ( fields[i].name == name ) return i;
    return fields.size();
}

// Get the name of field j.
std::string& DataTable::fieldName ( size_t j )
{
    return fields[j].name;
}
std::string DataTable::fieldName ( size_t j ) const
{
    return fields[j].name;
}

// Get the data type of field j.
unsigned char DataTable::fieldType ( size_t j ) const
{
    return fields[j].getType();
}

// Get the byte length of elements of field j.
unsigned char DataTable::fieldLength ( size_t j ) const
{
    return fields[j].getLen();
}

// Get the value of element i in field j.
template<typename T> void DataTable::get ( size_t i, size_t j, T* x )
    const
{
    fields[j].get(i, x);
}

// Get the value of element i in field j, interpreted as an integer.
int DataTable::geti ( size_t i, size_t j ) const
{
    return fields[j].geti(i);
}
int DataTable::geti ( size_t i, const std::string& fld ) const
{
    size_t j = findField(fld);
    if ( j == fields.size() ) return -9999;
    return fields[j].getf(i);
}

// Get the value of element i in field j, interpreted as a double.
double DataTable::getf ( size_t i, size_t j ) const
{
    return fields[j].getf(i);
}
double DataTable::getf ( size_t i, const std::string& fld ) const
{
    size_t j = findField(fld);
    if ( j == fields.size() ) return -9999.0;
    return fields[j].getf(i);
}

// Get the value of element i in field j, interpreted as a string.
std::string DataTable::getStr ( size_t i, size_t j ) const
{
    return fields[j].getStr(i);
}
std::string DataTable::getStr ( size_t i, const std::string& fld )
    const
{
    size_t j = findField(fld);
    if ( j == fields.size() ) return "";
    return fields[j].getStr(i);
}

// Get the null value for field j.
double DataTable::getNull ( size_t j ) const
{
    return fields[j].getNull();
}


// SETTERS //////////////////////////////////////////////////////////

// Set the value of element i in numeric field j.
template<typename T> void DataTable::set ( size_t i, size_t j,
    const T& x )
{
    fields[j].set(i, x);
}

// Set the value of element i in string field j.
void DataTable::setStr ( size_t i, size_t j,
    const std::string& x )
{
    fields[j].setStr(i, x);
}

// Set the null value for field j.
void DataTable::setNull ( size_t j, double x )
{
    fields[j].setNull(x);
}


// OPERATIONS ///////////////////////////////////////////////////////

// Add a field at position j of type dtyp of byte size dlen. If j
// is greater than the current number of fields, appends the field
// to the end of the fields vector.
void DataTable::addField ( size_t j, unsigned char dtyp,
    unsigned char dlen, const std::string& fldname,
    const std::string& units )
{
    if ( j > fields.size() ) j = fields.size();
    fields.insert(fields.begin() + j,
        DataVector(dtyp, dlen, n, fldname, units));
    fields[j].zero();
}
void DataTable::addField ( size_t j, unsigned char dtyp,
    unsigned char dlen, const std::string& fldname )
{
    addField(j, dtyp, dlen, fldname, "");
}
void DataTable::addField ( unsigned char dtyp, unsigned char dlen,
    const std::string& fldname, const std::string& units )
{
    addField(fields.size(), dtyp, dlen, fldname, units);
}
void DataTable::addField ( unsigned char dtyp, unsigned char dlen,
    const std::string& fldname )
{
    addField(fields.size(), dtyp, dlen, fldname, "");
}

// Remove the field at position j, if j is within bounds.
void DataTable::removeField ( size_t j )
{
    if ( j < fields.size() ) fields.erase(fields.begin() + j);
}

// Remove all fields with the given name.
void DataTable::removeField ( const std::string& fldname )
{
    std::vector<DataVector> newfields;
    newfields.reserve(fields.size());
    for ( size_t j = 0; j < fields.size(); ++j ) {
        if ( fields[j].name == fldname ) continue;
        newfields.push_back(fields[j]);
    }
    fields = newfields;
}

// Set the number of records in each field. If this enlarges the
// field, zeroes the new elements.
void DataTable::setNumRecords ( size_t nr )
{
    size_t old_n = n;
    n = nr;
    if ( fields.size() == 0 || n == old_n ) return;
    for ( size_t j = 0; j < fields.size(); ++j ) {
        fields[j].resize(nr);
        if ( fields[j].getType() == DataVector::T_STR ) {
            for ( size_t i = old_n; i < n; ++i )
                fields[j].set(i, 0);
        }
        else {
            for ( size_t i = old_n; i < n; ++i )
                fields[j].set(i, '\0');
        }
    }
}

// Set all records to 0. String records are set to the null
// character.
void DataTable::zero ()
{
    for ( size_t j = 0; j < fields.size(); ++j ) fields[j].zero();
}

// Report true if the given element is the same value as the no-data
// value in field j, given either element i or value x.
bool DataTable::isNull ( size_t i, size_t j ) const
{
    return fields[j].isNull(i);
}
bool DataTable::isNull ( double x, size_t j ) const
{
    return fields[j].isNull(x);
}


#endif // YOUNG_DATATABLE_20200629