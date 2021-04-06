# BYstdlib
Brennan Young's Standard Library

<!> This is an ongoing work-in-progress. It is also meant to be a personal educational tool, where I code up things from scratch to learn their ins-and-outs.

<!> And yes, I know that "good style" or convention of whatever gives C++ header files the .hpp extension instead of .h. It's been a long journey for me, and the .h is basically left over from the olden days when I was somehow dumber than I am now.

| Table | Type | Status | Description |
| --- | --- | --- | --- |
| [BinManip](#BinManip) | Functions | stable, incomplete | Functions for reading/writing binary. |
| [BSTree](#BSTree) | Object | stable | AVL tree, or balanced binary search tree |
| [CardinalDir](#CardinalDir) | Object | stable | Object for abstracting the concept of directionality in a square-tesselated framework. |
| [DataVector](#DataVector) | Object | stable | Object for storing a data type without needing for a template to define type. |
| [DataTable](#DataTable) | Object | stable | Object for managing a consistent set of DataVectors representing various data types. |
| [DLL](#DLL) | Object | stable | Doubly-linked list. |
| [fuzzy](#fuzzy) | Functions | stable | Fuzzy membership functions. |
| [Geometry](#Geometry) | Functions | stable | Functions for managing various miscellaneous aspects of geometry. This is slowly getting phased out. |
| [Grid](#Grid) | Object | stable | Object for interacting with a grid stored in a binary data file. |
| [integration](#integration) | Functions | stable | Functions for calculating the integral of a function. |
| [interpolate](#interpolate) | Functions | stable | Functions for interpolating between values. |
| [Matrix](#Matrix) | Object | stable | Matrix. |
| [ProgCounter](#ProgCounter) | Object | stable | Object for simplifying progress-tracking through various levels of a program. |
| [Regression](#Regression) | Functions | stable | Regression functions. |
| [RTree](#RTree) | Object | stable-ish | R* tree for spatial indexing. |
| [setManip](#setManip) | Functions | stable | Functions for manipulating std::set objects. |
| [sort](#sort) | Functions | stable | Sorting functions. |
| [sortedVector](#sortedVector) | Functions | stable | Functions for creating and maintaining a sorted vector. |
| [Statistics](#Statistics) | Object | stable | Object for storing and computing basic statistics. |

## BinManip

**readStrBinary** - read a string from an open binary file.
```cpp
std::string readStrBinary ( std::ifstream& f );
```

**writeStrBinary** - write a string to an open binary file.
```cpp
void writeStrBinary ( std::ofstream& f, const std::string& s );
```

## BSTree

The BSTree object is my own stab at implementing what is essentially the std::map, but requires the key to be an integer and is probably more messy and less efficient.

```cpp
bystd::BSTree<T> tree;
bystd::BSTree<T>::iterator it = tree.begin();
bystd::BSTree<T>::const_iterator cit = tree.begin(); // if tree is const
```

| Method | Arguments | Description |
| --- | --- | --- |
| ***Iteration*** | | |
| begin | | {iterator OR const_iterator} Get the element with the lowest-valued key. |
| last | | {iterator OR const_iterator} Get the element with the highest-valued key. |
| end | | {iterator OR const_iterator} Get the invalid element demarking the end of the tree. |
| at | int i | {iterator OR const_iterator} Get the element with key i. |
| find | int i | {iterator OR const_iterator} Get the element with key i
| ***Getters*** | | |
| size | | {size_t} Get the total number of elements stored in the tree. |
| height | | {int} Get the maximum height of the tree. |
| getKey | iterator | {int} Get the key of the given iterator. |
| ***Setters*** | | |
| insert | int i, T x | Add element x with key i to the tree. |
| remove | int i | Remove element with key i to the tree. |
| clear | | Remove all elements from the tree. |
| printKeys | | Prints the tree's keys using std::cout |

## CardinalDir

An object for representing cardinal directions and their diagonals (N, NE, E, SE, S, SW, W, NW). I created this object because I found keeping track of a "forward" direction and "turning" left or right in a grid was just annoying without the level of abstraction this object provides.

```cpp
CardinalDir d; // initialize to N
CardinalDir d("N"); // initialize to N
CardinalDir d(0,1); // initialize to direction from difference in (x,y) -- this example is N
```

| Method | Arguments | Description |
| --- | --- | --- |
| operator++ | | {CardinalDir} turn clockwise; return this. |
| operator-- | | {CardinalDir} counterclockwise; return this. |
| operator+ | int x | {CardinalDir} Get a CardinalDir turned 45 degrees to the right x times. |
| operator- | int x | {CardinalDir} Get a CardinalDir turned 45 degrees to the left x times. |
| operator- | CardinalDir d | {int} Get the number of 45-degree turns to turn clockwise toward d; negative if counter-clockwise. |
| x | | {int} Get 1 if direction is toward +x, -1 if toward -x, or 0 otherwise. |
| y | | {int} Get 1 if direction is toward +y, -1 if toward -y, or 0 otherwise. |

## DataVector

I wanted a data table handling any kind of data to get around the strong typing of C++, so I made the data vector to store values in binary format.

```cpp
DataVector v (DataVector::T_INT, 4, 100, "field name", "units");
```

| Constant | Description |
| --- | --- |
| T_INT | Signed integer-type. |
| T_UINT | Unsigned integer-type. |
| T_FLOAT | Floating-point or real type. |
| T_STR | Character string type. |

| Method | Arguments | Description |
| --- | --- | --- |
| ***Operators*** | | |
| operator[] | size_t i | {double} Get the value of element i interpreted as a floating-point number. |
| ***Getters*** | | |
| stats | | {Statistics} Get the vector's statistics (assuming they've already been computed using the calculateStatistics method). |
| size | | {size_t} Get the length of the vector. |
| getType | | {unsigned char} Get the constant for the vector's represented data type.
| getLen | | {unsigned char} Get the byte length of each of the vector's elements.
| get | size_t i, T* x | Get the value of element i and cast it and store it in x. |
| geti | size_t i | {int} Get the value of element i cast as an integer. |
| getf | size_t i | {double} Get the value of element j cast as a floating-point number. |
| getStr | size_t i | {std::string} Get the value of element j interpreted as a string. |
| getNull | | {double} Get the vector's no-data value. |
| ***Setters*** | | |
| set | size_t i, T x | Cast x to the vector's data type and store it in element i. |
| setStr | size_t i, std::string s | Store s in element i. |
| setNull | double x, bool change=false | Change the no-data value to x. If change is true, changes all existing no-data values to the new no-data value. |
| ***Operations*** | | |
| resize | size_t n | Resize the vector to n elements. Retains up to n elements that already existed. |
| zero | | Assigns all elements to zero. |
| isNull | size_t i | {bool} Get true if element i is no-data. |
| isNull | double x | {bool} Get true if x is no-data. |
| calculateStatistics | | Calculate the vector's statistics. |

## DataTable

Basically a collection of DataVector objects.

```cpp
DataTable table ();
```

| Method | Arguments | Description |
| --- | --- | --- |
| ***Getters*** | | |
| stats | size_t i | {Statistics} Get the statistics of vector i. |
| numFields | | {size_t} Get the number of DataVectors in the data table. |
| numRecords | | {size_t} Get the length of all DataVectors in the data table. |
| findField | std::string s | {size_t} Get the index of the DataVector with the name s. |
| fieldName | size_t i | {std::string} Get the name of DataVector i.
| fieldType | size_t i | {unsigned char} Get the constant for the type of DataVector i. |
| fieldLength | size_t i | {unsigned char} Get the byte-length of each element in DataVector i. |
| get | size_t i, size_t j, T* x | Get the value of element i of DataVector j and cast it and store it in x. |
| geti | size_t i, size_t j | {int} Get the value of element i of DataVector j interpreted as an integer. |
| geti | size_t i, std::string s | {int} Get the value of element i in the DataVector named s, interpreted as an integer. |
| getf | size_t i, size_t j | {double} Get the value of element i of DataVector j interpreted as a floating-point number. |
| getf | size_t i, std::string s | {double} Get the value of element i in the DataVector named s, interpreted as a float-point number. |
| getStr | size_t i, size_t j | {std::string} Get the string stored in element i of DataVector j. |
| getStr | size_t i, std::string s | {std::string} Get the string stored in element i of the DataVector named s. |
| getNull | size_t i | {double} Get the no-data value of DataVector i. |
| ***Setters*** | | |
| set | size_t i, size_t j, T x | Cast x to the vector's data type and store it in element i of DataVector j. |
| setStr | size_t i, size_t j, std::string s | Store s in element i of DataVector j. |
| setNull | size_t i, double x | Change the no-data value of DataVector i to x. |
| ***Operations*** | | |
| addField | size_t i, unsigned char t, unsigned chan n, std::string s, std::string u | Add a new DataVector before DataVector i of type t, byte-length n, named s and with values in units u. |
| addField | size_t i, unsigned char t, unsigned char n, std::string s | Add a new unitless DataVector before DataVector i of type t, byte-length n, named s. |
| addField | unsigned char, unsigned char, std::string, std::string | Add a new DataVector at the end of the list of DataVectors, of type t, byte-length n, named s and with values in units u. |
| addField | unsigned char, unsigned char, std::string | Add a new unitless DataVector at the end of the list of DataVectors, of type t, byte-length n, named s. |
| removeField | size_t i | Removed DataVector i. |
| removeField | std::string | Remove DataVector named i. |
| setNumRecords | size_t n | Resize all DataVectors to accomodate n elements. |
| zero | | Assign all elements to 0. |
| isNUll | size_t i, size_t j | {bool} Get true if element i in DataVector j is no-data. |
| isNull | double x, size_t j | {bool} Get true if x is no-data according to DataVector j. |

## DLL

Double-linked list. Basically, std::list but with shorthand for getting the iterator to element i, using at(i).

## fuzzy

Functions for various fuzzy functions.

**triangularMembership** - Get the membership \[0,1] of x in the triangle that is 0 at a, peaks at b, and returns to 0 at c.
```cpp
double triangularMembership (double x, double a, double b, double c);
```

***logisticSigmoid*** - Transform x with a logistic sigmoid \[0,1] centered on x0, with a curve peaking with a value of L, with a logistic growth rate of k.
```cpp
double logisticSigmoid (double x, double x0, double L, double k);
```

## Geometry

Various geometric functions. Not to be confused with my Geometry library.

## Grid

Grid object for interacting directly with data stored in binary on the disk -- the grid's data is not loaded into memory. This is a useful base class for other grid types.

```cpp
Grid grid ();                             // create an empty grid handle
Grid grid (1000, 1000, 1);                // create a grid handle for 1000 rows and columns and 1 band
Grid grid (Grid::INT, 4, 1000, 1000, 1);  // create a grid handle for 1000 rows and columns and 1 band

Grid newgrid (grid);                      // create a copy of the grid. If it has  a file open, tries to open that file.
Grid newgrid = grid;                      // create a copy of the grid. If it has a file open, tries to open that file.
```

| Constant | Description |
| --- | --- |
| INT | Signed integer-type. |
| UINT | Unsigned integer-type. |
| FLOAT | Floating-point or real type. |
| BSQ | Band-sequential file format. |
| BIL | Band interleave-by-line format. |
| BIP | Band interleave-by-pixel format. |

| Method | Arguments | Description |
| --- | --- | --- |
| ***File Access*** | | |
| open | std::string s | Opens the file with full path and filename s. |
| close | | Closes the file, if open. |
| ***Structure Getters*** | | |
| nrows | | {int} Get the number of rows represented in the grid file. |
| ncols | | {int} Get the number of columns represented in the grid file. |
| nbands | | {int} Get the number of bands represented in the grid file. |
| size | | {int} Get the number of elements in one band (nrows * ncols). |
| volume | | {int} Get the number of elements in the grid file (size * nbands). |
| index | int i, int j | Get the index for accessing the element at row i, column j. |
| index | int i, int j, int k | {int} Get the index for accessing the element at row i, column j, band k. |
| ***Structure Setters*** | | |
| setByteSwap | bool b | Set whether values should be byte-swapped when read from/written to the file. |
| setDimensions | int rows, int cols, int bands | Sets the expected dimensions of the grid. |
| setInterleave | unsigned char c | Sets the interleave to BSQ, BIL, or BIP, or no change if c is not recognized. |
| setDataType | unsigned char t, unsigned char s | Sets the data type to one of INT, UINT, or FLOAT with byte-length s for each element, or no change if t is not recognized. |
| ***Element Access*** | | |
| get | int i | {double} Get the value of element i, or nan if out of bounds. |
| get | int i, int j | {double} Get the value of the element at row i, column j, or nan if out of bounds. |
| get | int i, int j, int k | {double} Get the value of the element at row i, column j, band k, or nan if out of bounds. |
| set | T x, int i | Set element i to x. |
| set | T x, int i, int j | Set the element at row i, column j to x. |
| set | T x, int i, int j, int k | Set the element at row i, column j, band k to x. |

## Integration

Functions for performing integrations.

## Interpolate

Functions for interpolating.

## Matrix

My Matrix object.

## ProgCounter

My progress-counter object.

## Regression

## RTree

Basically an n-dimensional ABAB tree for spatial indexing. Designed to be like an R* tree.

## setManip

Functions for manupulating std::set objects.

## sort

Sorting functions.

## sortedVector

Like a vector, but keeps things in order.

## Statistics

My Statistics object.

| Members | Description |
| --- | --- |
| n | Number of valid or not-ignored samples in population. |
| min | Minimum value in population. |
| max | Maximum value in population. |
| sum | Sum of values in population. |
| sum2 | Sum of squares. |
| mean | Mean value of population. |
| var | Variance of population. |
| stdev | Standard deviation of population. |
