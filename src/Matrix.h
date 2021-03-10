/////////////////////////////////////////////////////////////////////
// Matrix                                                          //
// Class object for matrix operations                              //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 01/18/2019  - Brennan Young                                     //
// - created.                                                      //
// 03/15/2019  - Brennan Young                                     //
// - included library: cstdlib.                                    //
// - added operator() for two-dimensional element access.          //
// 02/14/2020 - Brennan Young                                      //
// - operator= now returns by reference.                           //
// - default constructor now takes optional arguments, removed     //
//   redundant constructors.                                       //
// 02/17/2020 - Brennan Young                                      //
// - moved default arguments to constructor declaration, as needed //
//   for template classes.                                         //
// 05/08/2020 - Brennan Young                                      //
// - removed some redundant variables.                             //
// 05/13/2020 - Brennan Young                                      //
// - removed some unused variables.                                //
// 06/16/2020 - Brennan Young                                      //
// - added operator/ scalar division.                              //
// 06/17/2020 - Brennan Young                                      //
// - removed redundant get methods.                                //
// - changed element access operators' return type to const T&.    //
// - added implicit column vector support: operator(...) and       //
//   indexColVec(...).                                             //
// - added addition, subtraction, and negation operators.          //
// - moved matrix algebra operations outside of the class          //
//   definition.                                                   //
// - made matrix algebra operations permit operations with other   //
//   data types (added U to template).                             //
// - added non-member functions: abs, sqrt, cos, sin, and          //
//   skew.                                                         //
// - added insertSubMatrix().                                      //
// - added diag().                                                 //
// 06/18/2020 - Brennan Young                                      //
// - added sub().                                                  //
// - added reshape().                                              //
// - added addrow().                                               //
// - added addcol().                                               //
// - added mag().                                                  //
// 06/22/2020 - Brennan Young                                      //
// - corrected typo of Matrix addition.                            //
// 06/30/2020 - Brennan Young                                      //
// - wrapped dynamic memory allocation in try/catch statement.     //
// 07/01/2020 - Brennan Young                                      //
// - added pow (scalar power) function.                            //
// - modified +,-,* operators between matrices now treat one-      //
//   element matrices as scalars.                                  //
// 03/10/2021 - Brennan Young                                      //
// - added elemMult and elemDivide.                                //
/////////////////////////////////////////////////////////////////////

#ifndef MATRIX_20190118
#define MATRIX_20190118

#include <cstdlib>
#include <cmath>
#include <iostream>

template <typename T>
class Matrix {
private:
    T * data;
    int nr, nc, n;
public:
    // constructors, destructor
    Matrix(int r=1, int c=1);
    Matrix(const Matrix<T>&);
    ~Matrix();
    
    // operators
    Matrix<T> & operator=(const Matrix<T>&);
    const T& operator[](int) const;
    T& operator[](int);
    const T& operator()(int, int) const;
    T& operator()(int, int);
    const T& operator()(int) const;
    T& operator()(int);
    Matrix<T> operator-() const;
    
    // getters
    int index(int, int) const;
    int indexColVec(int) const;
    int size() const;
    int nrows() const;
    int ncols() const;
    Matrix<T> diag() const;
    Matrix<T> implicitColVec() const;
    Matrix<T> implicitColVec(int, int) const;
    
    // operations
    Matrix<T>& identity();
    Matrix<T> identity() const;
    Matrix<T> transpose() const;
    Matrix<T> reshape(int, int) const;
    Matrix<T>& addrow(int);
    Matrix<T>& addcol(int);
    
    Matrix<T>& insertSubMatrix(int, int, const Matrix<T>&);
    Matrix<T> sub(int, int, int, int) const;
    
    double mag() const;
    double determinant() const;
    Matrix<T> inverse(double) const;
    Matrix<T> inverse() const;
    
    Matrix<T> leastSquares(const Matrix<T> &) const;
    Matrix<T> gaussianElimination(const Matrix<T>&);
    void LUDecomp(Matrix<T>*, Matrix<T>*);
    Matrix<T> LUDbacksub(const Matrix<T>&);
    Matrix<T> Crout(const Matrix<T>&);
}; // Matrix

// CONSTRUCTORS / DESTRUCTOR ////////////////////////////////////////

// Default constructor. If given just one argument, constructs a
// column vector.
template<typename T> Matrix<T>::Matrix ( int r, int c )
: nr(r), nc(c), n(r*c)
{
    try {
        data = new T[n];
    }
    catch (...) {
        std::cout << "Matrix: memory allocation error.\n";
        exit(0);
    }
    for ( int i = 0; i < n; ++i ) data[i] = 0;
}

// Copy constructor.
template<typename T> Matrix<T>::Matrix ( const Matrix<T>& A )
{
    nr = A.nr;
    nc = A.nc;
    n = A.n;
    try {
        data = new T[n];
    }
    catch (...) {
        std::cout << "Matrix: memory allocation error.\n";
        exit(0);
    }
    for ( int i = 0; i < n; ++i ) data[i] = A.data[i];
}

// Destructor.
template<typename T> Matrix<T>::~Matrix ()
{
    delete [] data;
}

// OPERATORS ////////////////////////////////////////////////////////

// Assignment.
template<typename T>
Matrix<T>& Matrix<T>::operator= ( const Matrix<T>& A )
{
    if ( this == &A ) return *this;
    
    nr = A.nr;
    nc = A.nc;
    n = A.n;
    
    if (data) delete [] data;
    try {
        data = new T[n];
    }
    catch (...) {
        std::cout << "Matrix: memory allocation error.\n";
        exit(0);
    }
    
    for ( int i = 0; i < n; ++i ) data[i] = A.data[i];
    
    return *this;
}

// Element access.
template<typename T> const T& Matrix<T>::operator[] ( int i ) const
{
    return data[i];
}

// Element access.
template<typename T> T& Matrix<T>::operator[] ( int i )
{
    return data[i];
}

// Element access.
template<typename T>
const T& Matrix<T>::operator() ( int i, int j ) const
{
    return data[index(i, j)];
}

// Element access.
template <typename T> T& Matrix<T>::operator() ( int i , int j )
{
    return data[index(i, j)];
}

// Element access in the matrix's implicit column vector.
template <typename T> const T& Matrix<T>::operator() ( int i ) const
{
    return data[indexColVec(i)];
}

// Element access in the matrix's implicit column vector.
template <typename T> T& Matrix<T>::operator() ( int i )
{
    return data[indexColVec(i)];
}

// Matrix negation operator.
template <typename T>
Matrix<T> Matrix<T>::operator- () const
{
    Matrix<T> A (nr, nc);
    for ( int i = 0; i < n; ++i ) A[i] = -data[i];
    return A;
}

// GETTERS //////////////////////////////////////////////////////////

// Get index in matrix data array for row i, column j.
template<typename T> int Matrix<T>::index ( int i, int j ) const
{
    return i * nc + j;
}

// Get index in matrix data array for the matrix's implicit column
// vector.
template<typename T> int Matrix<T>::indexColVec ( int k ) const
{
    int i = k % nr;
    int j = (k - i) / nr;
    return index(i, j);
}

template<typename T> int Matrix<T>::size () const
{
    return n;
}

template<typename T> int Matrix<T>::nrows () const
{
    return nr;
}

template<typename T> int Matrix<T>::ncols () const
{
    return nc;
}

// Get the implicit column vector of the given matrix.
template <typename T>
Matrix<T> Matrix<T>::implicitColVec () const
{
    Matrix<T> A (n);
    
    // jump to column-major position
    int i, k = 0;
    for ( int j = 0; j < nc; ++j ) {
        for ( i = 0; i < nr; ++i )
            A[k++] = (*this)(i, j);
    }
    return A;
}

// Get the a column vector filled with the diagonal elements.
template <typename T>
Matrix<T> Matrix<T>::diag () const
{
    Matrix<T> A(nr < nc ? nr : nc);
    for ( int i = 0; i < A.size(); ++i ) A[i] = (*this)(i,i);
    return A;
}

// Get the implicit column vector of the given matrix for the given
// index range [k,k1) in that vector (not including k1).
template <typename T>
Matrix<T> Matrix<T>::implicitColVec ( int k, int k1 ) const
{
    Matrix<T> A (k1 - k);
    
    // jump to column-major position
    int i = k % nr;
    int j = (k - i) / nr;
    int kk = 0;
    
    for ( ; j < nc && k < k1; ++j ) {
        for ( ; i < nr && k < k1; ++i, ++k )
            A[kk++] = (*this)(i, j);
        i = 0;
    }
    return A;
}

// MATRIX OPERATIONS ////////////////////////////////////////////////

// Sets the matrix to its identity.
template<typename T> Matrix<T>& Matrix<T>::identity ()
{
    int i, j, k = 0;
    for ( i = 0; i < nr; ++i )
        for ( j = 0; j < nc; ++j )
            data[k++] = (i == j ? 1 : 0);
    return *this;
}

// Returns this matrix's identity matrix.
template<typename T> Matrix<T> Matrix<T>::identity () const
{
    Matrix<T> I (nr, nc);
    int i, j, k = 0;
    for ( i = 0; i < nr; ++i )
        for ( j = 0; j < nc; ++j )
            I[k++] = (i == j ? 1 : 0);
    return I;
}

// Returns a transpose of the matrix.
template<typename T> Matrix<T> Matrix<T>::transpose () const
{
    int i, j, k = 0;
    Matrix<T> A(nc, nr);
    for ( i = 0; i < nr; ++i )
        for ( j = 0; j < nc; ++j )
            A(j, i) = data[k++];
    return A;
}

// Returns a matrix reshaped to r rows and c columns. Preserves
// coumn-major ordering. If the new matrix is too small, truncates.
// If the new matrix is too large, fills with zeroes.
template<typename T>
Matrix<T> Matrix<T>::reshape ( int r , int c ) const
{
    Matrix<T> A (r,c);
    for ( int i = 0; i < n && i < A.size(); ++i ) A(i) = (*this)(i);
    return A;
}

// Adds x new blank rows the matrix.
template<typename T>
Matrix<T>& Matrix<T>::addrow( int x )
{
    Matrix<T> A(nr+x, nc);
    for ( int i = 0; i < n; ++i ) A[i] = data[i];
    *this = A;
    return *this;
}

// Adds x new blank columns to the matrix.
template<typename T>
Matrix<T>& Matrix<T>::addcol( int x )
{
    Matrix<T> A(nr, nc+x);
    for ( int i = 0; i < nr; ++i ) {
        for ( int j = 0; j < nr; ++j )
            A(i,j) = (*this)(i,j);
    }
    *this = A;
    return *this;
}

// Inserts a sub-matrix into the matrix, or as much of the sub-matrix
// as would fit, with the sub-matrix's (0,0) being located at (r,c).
template<typename T>
Matrix<T>& Matrix<T>::insertSubMatrix ( int r, int c ,
    const Matrix<T>& A )
{
    int j;
    for ( int i = 0; i < A.nr && i+r < nr; ++i ) {
        for ( j = 0; j < A.nc && j+c < nc; ++j )
            (*this)(i+r, j+c) = A(i,j);
    }
    return *this;
}

// Returns a sub matrix for rows [i0,i1) and columns [j0,j1), not
// including row i1 or column j1.
template<typename T>
Matrix<T> Matrix<T>::sub ( int i0, int i1, int j0, int j1 ) const
{
    int r = i1-i0;
    int c = j1-j0;
    Matrix<T> A (r, c);
    int j;
    for ( int i = 0; i < r; ++i ) {
        for ( j = 0; j < c; ++j )
            A(i,j) = (*this)(i+i0,j+j0);
    }
    return A;
}

// Returns the magnitude of the matrix. If the matrix has more than
// one row and more than one column, this is the Frobenius norm of
// the matrix.
template<typename T> double Matrix<T>::mag () const
{
    double sum = 0.0;
    for ( int i = 0; i < n; ++i ) sum += (data[i] * data[i]);
    return sqrt(sum);
}

// Returns the determinant of the matrix. The matrix should be
// square.
template<typename T> double Matrix<T>::determinant () const
{
    if ( n == 1 ) return data[0];
    
    int     i, j, m, d, c;
    double  det = 0.0;
    double  s = 1.0;
    
    Matrix<T> B (nr-1, nc-1);
    
    for ( c = 0; c < nc; ++c ) {
        // minor expansion
        m = 0;
        for ( i = 1; i < nr; ++i ) {
            d = 0;
            for ( j = 0; j < nc; ++j ) {
                if ( j == c ) continue;
                B(m, d++) = (*this)(i, j);
            }
            ++m;
        }
        det += s * (*this)(0, c) * B.determinant();
        s *= -1;
    }
    
    return det;
}

// Returns the inverse of the matrix. The matrix should be square.
// -- Arguments --
// det : the determinant of the matrix. This is to force
//       determination of whether or not the matrix has an inverse
//       prior to calling this function.
template<typename T>
Matrix<T> Matrix<T>::inverse ( double det ) const
{
    Matrix<T> B (nr-1, nc-1);
    Matrix<T> I (nr, nc);
    double cofactor;
    int p, q, m, d, i, j;
    for ( q = 0; q < nc; q++ ) {
        for ( p = 0; p < nc; p++ ) {
            // minor expansion
            m = 0;
            for ( i = 0; i < nr; ++i ) {
                if ( i == q ) continue;
                d = 0;
                for ( j = 0; j < nc; ++j ) {
                    if ( j == p ) continue;
                    B(m, d++) = (*this)(i, j);
                }
                ++m;
            }
            
            cofactor = pow(-1, p+q) * B.determinant();
            I(p, q) = cofactor / det;
        }
    }
    
    return I;
}
template<typename T> Matrix<T> Matrix<T>::inverse() const
{
    return this->inverse(this->determinant());
}

// Solves for the coefficient vector b in the system of linear
// equations Mb=z.
// -- Source --
// Created 5/19/2017 by Brennan Young for BishopMatrix.
// -- Arguments --
// Z : vector of solutions. Z should be a column vector.
// -- Returns --
// A column vector of coefficients.
template<typename T>
Matrix<T> Matrix<T>::leastSquares ( const Matrix<T> & Z ) const
{
    Matrix Dt = this->transpose();
    return (Dt * (*this)).inverse() * (Dt * Z);
}

// Returns the augmented U matrix in the LU decomposition of the
// matrix, by Gaussian elimination.
// -- Source --
// Created 11/28/2018 by Brennan Young for BishopMatrix.
// -- Arguments --
// Z : Solution vector for the system of linear equations represented
//     by the matrix. Z should be a column vector.
// -- Returns --
// Augmented matrix U.
template<typename T>
Matrix<T> Matrix<T>::gaussianElimination ( const Matrix<T> & Z )
{
    if ( nr != nc ) {
        std::cout << "Must be a square matrix\n";
        exit(0);
    }
    
    int i, j, m, n;
    Matrix U (nrows, ncols + 1);
    Matrix I (nrows, ncols);
    
    // initialize augmented matrix
    for ( i = 0; i < nrows; ++i ) {
        for ( j = 0; j < ncols; ++j )
            U(i, j) = (*this)(i, j);
        U(i, j) = Z[i];
    }
    
    // initialize I as identity matrix
    I.identity();
    
    // Gaussian elimination
    for ( j = 0; j < ncols; ++j ) {
        // set up elementary row operator
        for ( m = 0; m < nrows; ++m ) {
            for ( n = 0; n < ncols; ++n ) {
                if ( m == n ) break; // next row: upper half of I
                                     //   remains as identity
                else if ( n == j )
                    I.SetValue(
                        m, n, -U.GetValue(m, n) / U.GetValue(n, n));
                else I.SetValue(m, n, 0.0);
            }
        }
        
        // perform the row operation
        U = I * U;
    }
    
    return U;
}

// Decomposes the matrix into upper and lower triangular matrices by
// Gaussian elimination. The matrix must be square.
// -- Source --
// Created 11/28/2018 by Brennan Young for BishopMatrix.
// -- Arguments --
// L, U : pointers to the lower and upper triangular matrices.
// -- Returns --
// Nothing -- resizes and fills L and U.
template<typename T>
void Matrix<T>::LUDecomp ( Matrix<T> * L, Matrix<T> * U )
{
    if ( nr != nc ) {
        std::cout << "Must be a square matrix\n";
        exit(0);
    }
    
    int j, m, d;
    Matrix I (nr, nc);
    *L = Matrix<T>(nr, nc);
    *U = *this;
    
    // initialize I and L as identity matrix
    I.identity();
    L->identity();
    
    // Gaussian elimination
    for ( j = 0; j < nc; ++j ) {
        // set up elementary row operator
        for ( m = 0; m < nr; ++m ) {
            for ( d = 0; d < nc; ++d ) {
                if ( m == d ) break; // next row: upper half of I
                                     //   remains as identity
                else if ( d == j ) {
                    I(m, d) = -(*U)(m, d) / (*U)(d, d);
                    (*L)(m, d) = -I(m, d);
                }
                else I(m, d) = 0;
            }
        }
        
        // perform the row operation
        *U = I * *U;
    }
}

// Returns the coefficients to the linear equation by LU
// decomposition, forward-substitution to solve Ly=b, and then
// back-substitution.
// -- Source --
// Created 11/28/2018 by Brennan Young for BishopMatrix.
// -- Arguments --
// z : solution vector for the system of linear equations represented
//     by the matrix. z should be a column vector.
// -- Returns --
// Vector of coefficients c. c is a column vector.
template<typename T>
Matrix<T> Matrix<T>::LUDbacksub ( const Matrix<T> & z )
{
    // U matrix from Gaussian elimination
    Matrix U = this->gaussianElimination(z);
    
    // back-substitution
    int i, j;
    Matrix<T> c (nr);
    for ( i = nc - 1; i >= 0; --i ) {
        c[i] = U(i, nc); // solution after decomposition
        // subtract previously solved coefficients
        for ( j = nc - 1; j > i; --j ) c[i] -= U(i, j) * c[j];
        c[i] /= U(i, i);
    }
    
    return c;
}

// Solves for the coefficients of a linear system of equations Ac=z
// by Crout's method of LU decomposition.
// -- Source --
// Created 11/28/2018 by Brennan Young for BishopMatrix.
// -- Arguments --
// z : vector of solutions to the linear equation. z should be a row
//     vector.
// -- Returns --
// Vector of coefficients.
template<typename T>
Matrix<T> Matrix<T>::Crout ( const Matrix<T> & z )
{
    int         i, j;
    T           sum;
    Matrix      L (nr, nc);
    Matrix      U (nr, nc);
    Matrix<T>   y (nr);
    Matrix<T>   c (nr);
    
    // construct the lower and upper triangular matrices
    this->LUDecomp(&L, &U);
    
    // forward substitution to solve Ly=z
    y[0] = z[0] / L(0, 0);
    for ( i = 1; i < nr; ++i ) {
        sum = 0.0;
        for ( j = 0; j < i; ++j )
            sum += L(i, j) * y[j];
        y[i] = (z[i] - sum) / L(i, i);
    }
    
    // backward substitution to solve Uc=y
    i = nr - 1;
    c[i] = y[i] / U(i, i);
    for ( i = nr - 2; i >= 0; --i ) {
        sum = 0.0;
        for ( j = i + 1; j < nr; ++j )
            sum += U(i, j) * c[j];
        c[i] = (y[i] - sum) / U(i, i);
    }
    
    return c;
}


/////////////////////////////////////////////////////////////////////
// Non-member Matrix operations /////////////////////////////////////
/////////////////////////////////////////////////////////////////////

// Matrix addition with scalar.
template <typename T, typename U>
Matrix<T> operator+ ( const Matrix<T>& M, const U& x )
{
    Matrix<T> A (M.nrows(), M.ncols());
    for ( int i = 0; i < M.size(); ++i ) A[i] = M[i] + x;
    return A;
}

template <typename T, typename U>
Matrix<T> operator+ ( const U& x, const Matrix<T>& M )
{
    return M+x;
}

// Matrix entrywise sum. Limits the output matrix to the smallest
// dimensions between the two input matrices.
template <typename T, typename U>
Matrix<T> operator+ ( const Matrix<T>& A, const Matrix<U>& B )
{
    if ( A.size() == 1 ) return A[0] + B;
    else if ( B.size() == 1 ) return A + B[0];
    
    int nr = A.nrows() < B.nrows() ? A.nrows() : B.nrows();
    int nc = A.ncols() < B.ncols() ? A.ncols() : B.ncols();
    Matrix<T> C (nr, nc);
    int j;
    for ( int i = 0; i < nr; ++i ) {
        for ( j = 0; j < nc; ++j )
            C(i, j) = A(i, j) + B(i, j);
    }
    return C;
}

// Matrix subtraction with scalar.
template <typename T, typename U>
Matrix<T> operator- ( const Matrix<T>& M, const U& x )
{
    Matrix<T> A (M.nrows(), M.ncols());
    for ( int i = 0; i < M.size(); ++i ) A[i] = M[i] - x;
    return A;
}

template <typename T, typename U>
Matrix<T> operator- ( const U& x, const Matrix<T>& M )
{
    Matrix<T> A (M.nrows(), M.ncols());
    for ( int i = 0; i < M.size(); ++i ) A[i] = x - M[i];
    return A;
}

// Matrix subtraction with matrix. Limits the output matrix to the
// smallest dimensions between the two input matrices.
template <typename T, typename U>
Matrix<T> operator- ( const Matrix<T>& A, const Matrix<U>& B )
{
    if ( A.size() == 1 ) return A[0] - B;
    else if ( B.size() == 1 ) return A - B[0];
    
    int nr = A.nrows() < B.nrows() ? A.nrows() : B.nrows();
    int nc = A.ncols() < B.ncols() ? A.ncols() : B.ncols();
    
    Matrix<T> C (nr, nc);
    int j;
    for ( int i = 0; i < nr; ++i ) {
        for ( j = 0; j < nc; ++j )
            C(i, j) = A(i, j) - B(i, j);
    }
    return C;
}

// Matrix multiplication with scalar.
template<typename T, typename U>
Matrix<T> operator* ( const Matrix<T>& M, const U& x )
{
    Matrix<T> A (M.nrows(), M.ncols());
    for ( int i = 0; i < M.size(); ++i ) A[i] = M[i] * x;
    return A;
}

template <typename T, typename U>
Matrix<T> operator* ( const U& x, const Matrix<T>& M )
{
    return M*x;
}

// Matrix multiplication with matrix. The number of columns in the
// left matrix must be the same as the number of rows in the right
// matrix.
template<typename T, typename U>
Matrix<T> operator* ( const Matrix<T>& A, const Matrix<U>& B )
{
    if ( A.size() == 1 ) return A[0] * B;
    else if ( B.size() == 1 ) return A * B[0];
    
    int i, j, k;
    Matrix<T> C (A.nrows(), B.ncols());
    for ( i = 0; i < C.nrows(); ++i ) {
        for ( j = 0; j < C.ncols(); ++j ) {
            C(i,j) = 0;
            for ( k = 0; k < A.ncols(); ++k )
                C(i,j) += A(i, k) * B(k, j);
        }
    }
    
    return C;
}

// Elementwise multiplication with matrix. Both matrices must have
// the same dimensions, or the the product will be incomplete.
template<typename T, typename U>
Matrix<T> elemMult ( const Matrix<T>& A, const Matrix<U>& B )
{
    int i, j, k;
    int nr = A.nrows() < B.nrows() ? A.nrows() : B.nrows();
    int nc = A.ncols() < B.ncols() ? A.ncols() : B.ncols();
    Matrix<T> C (nr, nc);
    for ( i = 0; i < nr; ++i ) {
        for ( j = 0; j < nc; ++j )
            C(i,j) = A(i,j) * B(i,j);
    }
    
    return C;
}

// Matrix division with scalar.
template <typename T, typename U>
Matrix<T> operator/ ( const Matrix<T>& M, const U& x )
{
    Matrix<T> A (M.nrows(), M.ncols());
    for ( int i = 0; i < M.size(); ++i ) A[i] = M[i] / x;
    return A;
}

template <typename T, typename U>
Matrix<T> operator/ ( const U& x, const Matrix<T>& M )
{
    Matrix<T> A (M.nrows(), M.ncols());
    for ( int i = 0; i < M.size(); ++i ) A[i] = x / M[i];
    return A;
}

// Elementwise division with matrix. Both matrices must have
// the same dimensions, or the the result will be incomplete.
template<typename T, typename U>
Matrix<T> elemDivide ( const Matrix<T>& A, const Matrix<U>& B )
{
    int i, j, k;
    int nr = A.nrows() < B.nrows() ? A.nrows() : B.nrows();
    int nc = A.ncols() < B.ncols() ? A.ncols() : B.ncols();
    Matrix<T> C (nr, nc);
    for ( i = 0; i < nr; ++i ) {
        for ( j = 0; j < nc; ++j )
            C(i,j) = A(i,j) / B(i,j);
    }
    
    return C;
}

// Get the absolute value for each element in the matrix.
template <typename T>
Matrix<T> abs ( Matrix<T> M )
{
    for ( int i = 0; i < M.size(); ++i ) M[i] = fabs(M[i]);
    return M;
}

// Raise each element of the matrix to the given power.
template <typename T>
Matrix<T> pow ( Matrix<T> M, double x )
{
    for ( int i = 0; i < M.size(); ++i ) M[i] = pow(M[i], x);
    return M;
}

// Get the square root of each element in the matrix.
template <typename T>
Matrix<T> sqrt ( Matrix<T> M )
{
    for ( int i = 0; i < M.size(); ++i ) M[i] = sqrt(M[i]);
    return M;
}

// Get the cosine for each element in the matrix.
Matrix<double> cos ( Matrix<double> M )
{
    for ( int i = 0; i < M.size(); ++i ) M[i] = cos(M[i]);
    return M;
}

// Get the sine for each element in the matrix.
Matrix<double> sin ( Matrix<double> M )
{
    for ( int i = 0; i < M.size(); ++i ) M[i] = sin(M[i]);
    return M;
}

// Create a skew-symmetric matrix from 3-vector S. The matrix is
// equivalent to taking the cross product of the vector:
//   cross(a,b) == skew(a)*b
// See also: UNSKEW
// -- Source --
// Andrew Simon, TAMU, 2017.
Matrix<double> skew ( const Matrix<double>& s )
{
    Matrix<double> S(3,3);
    S(0,0) = 0.0;   S(0,1) = -s[2]; S(0,2) = s[1];
    S(1,0) = s[2];  S(1,1) = 0.0;   S(1,2) = -s[0];
    S(2,0) = -s[1]; S(2,1) = s[0];  S(2,2) = 0.0;
    return S;
}


#endif // MATRIX_20190118