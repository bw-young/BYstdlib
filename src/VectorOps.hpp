/////////////////////////////////////////////////////////////////////
// Vector operations.                                              //
// i.e., for linear algebra.                                       //
//                                                                 //
// This should work on any container with operator[] element       //
// access and the size() method that reports the valid length of,  //
// or number of elements in the container.                         //
/////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////
// -- HISTORY ---------------------------------------------------- //
// 04/08/2021 - Brennan Young                                      //
// - created.                                                      //
/////////////////////////////////////////////////////////////////////

template <class VecA, class VecB>
double dot ( const VecA& A, const VecB& B )
{
    double sum = 0.0;
    for ( unsigned int i = 0; i < (unsigned int) A.size()
            && i < (unsigned int) B.size(); ++i )
        sum += A[i] * B[i];
    return sum;
}