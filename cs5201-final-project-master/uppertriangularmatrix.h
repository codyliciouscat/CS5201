/**
 * @file uppertriangularmatrix.h
 * @author Caleb Berg (cabr29@mst.edu)
 * @brief Header file for UpperTriangularMatrix class
 * @date 2019-04-03
 */

#ifndef UPPER_TRIANGULAR_MATRIX_H
#define UPPER_TRIANGULAR_MATRIX_H

template< typename T >
class DenseMatrix;

template< typename T >
class LowerTriangularMatrix;

#include <algorithm>
#include "abstractmatrix.h"
#include "densematrix.h"
#include "lowertriangularmatrix.h"

using namespace std;

/**
 * @brief An upper triangular matrix
 * 
 * @pre Type T must have =, +=, -=, unary -, binary * defined
 */
template< typename T >
class UpperTriangularMatrix : public AbstractMatrix< T >
{
    private:
        Vector< T > * data;
        int n;

    public:
        /**
         * @brief Default constructor
         * 
         * @pre None
         * @post Matrix is constructed to be size 0x0
         */
        UpperTriangularMatrix();

        /**
         * @brief Parameterized constructor
         * 
         * @param n Dimensions of square matrix
         * 
         * @pre Dimensions must be nonnegative
         * @post Matrix is constructed to be size n x n
         * @throws Error< PARAMETER_ERROR > if precondition isn't met
         */
        UpperTriangularMatrix( const int n );

        /**
         * @brief Copy constructor
         * 
         * @param source Matrix to copy
         * 
         * @pre None
         * @post Matrix is constructed as a copy of source
         */
        UpperTriangularMatrix( const UpperTriangularMatrix< T > & source );

        /**
         * @brief Conversion constructor
         * 
         * @param source Matrix to copy
         * 
         * @pre None
         * @post Matrix is constructed as a conversion of source
         */
        UpperTriangularMatrix( const DiagonalMatrix< T > & source );

        /**
         * @brief Destructor
         * 
         * @pre None
         * @post Matrix is destroyed
         */
        ~UpperTriangularMatrix();

        /**
         * @brief Assignment operator
         * 
         * @param source Matrix to copy
         * 
         * @pre Vector< T > must have operator=() defined
         * @post Old matrix data is replaced by a copy of source
         * @return Reference to calling matrix
         */
        UpperTriangularMatrix< T > & operator=( const UpperTriangularMatrix< T > & source );
        
        /**
         * @brief Accessor
         * 
         * @param row Row in matrix
         * @param col Column in matrix
         * 
         * @pre 0 <= row < getRows() and 0 <= col < getColumns()
         * @return Element at matrix location
         * @throws Error< INDEX_ERROR > if precondition not met
         */
        T get( const int row, const int col ) const;

        /**
         * @brief Mutator
         * 
         * @param row Row in matrix
         * @param col Column in matrix
         * @param val Value to store in matrix
         * 
         * @pre 0 <= row < getRows() and 0 <= col < getColumns()
         * @post Value is stored in matrix at location (row, col)
         * @throws Error< INDEX_ERROR > if precondition not met
         */
        void set( const int row, const int col, const T & val );
        
        /**
         * @brief Vector multiplication
         * 
         * @param rhs Vector to multiply by
         * 
         * @pre getColumns() == rhs.size()
         * @return Resulting vector from multiplication
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition not met
         */
        Vector< T > operator*( const Vector< T > & rhs ) const;
        
        /**
         * @brief Get number of columns in matrix
         * 
         * @pre None
         * @return Number of columns in matrix
         */
        int getColumns() const;

        /**
         * @brief Get number of rows in matrix
         * 
         * @pre None
         * @return Number of rows in matrix 
         */
        int getRows() const;

        /**
         * @brief Check if matrix is lower triangular
         * 
         * @pre None
         * @return true if matrix is lower triangular
         */
        bool isLowerTriangular() const;

        /**
         * @brief Check if matrix is upper triangular
         * 
         * @pre None
         * @return true if matrix is upper triangular
         */
        bool isUpperTriangular() const;

        /**
         * @brief Check if matrix is diagonal
         * 
         * @pre None
         * @return true if matrix is diagonal
         */
        bool isDiagonal() const;

        /**
         * @brief Check if matrix is symmetric
         * 
         * @pre None
         * @return true if matrix is symmetric
         */
        bool isSymmetric() const;

        /**
         * @brief Check if matrix is tridiagonal
         * 
         * @pre None
         * @return true if matrix is tridiagonal
         */
        bool isTridiagonal() const;

        /**
         * @brief Matrix addition
         * 
         * @param rhs Matrix to add by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by adding corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator+( const DenseMatrix< T > & rhs ) const;
        
        /**
         * @brief Matrix addition
         * 
         * @param rhs Matrix to add by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by adding corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator+( const LowerTriangularMatrix< T > & rhs ) const;
        
        /**
         * @brief Matrix addition
         * 
         * @param rhs Matrix to add by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by adding corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        UpperTriangularMatrix< T > operator+( const UpperTriangularMatrix< T > & rhs ) const;

        /**
         * @brief Matrix addition
         * 
         * @param rhs Matrix to add by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by adding corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        UpperTriangularMatrix< T > operator+( const DiagonalMatrix< T > & rhs ) const;

        /**
         * @brief Matrix addition
         * 
         * @param rhs Matrix to add by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by adding corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator+( const SymmetricMatrix< T > & rhs ) const;

        /**
         * @brief Matrix addition
         * 
         * @param rhs Matrix to add by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by adding corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator+( const TridiagonalMatrix< T > & rhs ) const;
        
        /**
         * @brief Matrix addition
         * 
         * @param rhs Matrix to add by
         * 
         * @pre Matrices must have the same number of rows and columns
         *      Type T must have operator+=() defined
         * @post Matrix elements are summed with corresponding elements in other matrix
         * @return Reference to calling matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        UpperTriangularMatrix< T > & operator+=( const DiagonalMatrix< T > & rhs );

        /**
         * @brief Matrix addition
         * 
         * @param rhs Matrix to add by
         * 
         * @pre Matrices must have the same number of rows and columns
         *      Type T must have operator+=() defined
         * @post Matrix elements are summed with corresponding elements in other matrix
         * @return Reference to calling matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        UpperTriangularMatrix< T > & operator+=( const UpperTriangularMatrix< T > & rhs );
        
        /**
         * @brief Matrix subtraction
         * 
         * @param rhs Matrix to subtract by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by subtracting corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator-( const DenseMatrix< T > & rhs ) const;
        
        /**
         * @brief Matrix subtraction
         * 
         * @param rhs Matrix to subtract by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by subtracting corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator-( const LowerTriangularMatrix< T > & rhs ) const;
        
        /**
         * @brief Matrix subtraction
         * 
         * @param rhs Matrix to subtract by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by subtracting corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        UpperTriangularMatrix< T > operator-( const UpperTriangularMatrix< T > & rhs ) const;
        
        /**
         * @brief Matrix subtraction
         * 
         * @param rhs Matrix to subtract by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by subtracting corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        UpperTriangularMatrix< T > operator-( const DiagonalMatrix< T > & rhs ) const;

        /**
         * @brief Matrix subtraction
         * 
         * @param rhs Matrix to subtract by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by subtracting corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator-( const SymmetricMatrix< T > & rhs ) const;

        /**
         * @brief Matrix subtraction
         * 
         * @param rhs Matrix to subtract by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by subtracting corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator-( const TridiagonalMatrix< T > & rhs ) const;

        /**
         * @brief Matrix subtraction
         * 
         * @param rhs Matrix to subtract by
         * 
         * @pre Matrices must have the same number of rows and columns
         *      Type T must have operator-=() defined
         * @post Matrix elements are subtracted by corresponding elements in other matrix
         * @return Reference to calling matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        UpperTriangularMatrix< T > & operator-=( const DiagonalMatrix< T > & rhs );

        /**
         * @brief Matrix subtraction
         * 
         * @param rhs Matrix to subtract by
         * 
         * @pre Matrices must have the same number of rows and columns
         *      Type T must have operator-=() defined
         * @post Matrix elements are subtracted by corresponding elements in other matrix
         * @return Reference to calling matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        UpperTriangularMatrix< T > & operator-=( const UpperTriangularMatrix< T > & rhs );
        
        /**
         * @brief Negate matrix
         * 
         * @pre None
         * @return Copy of calling matrix but with all elements negated
         */
        UpperTriangularMatrix< T > operator-() const;
        
        /**
         * @brief Matrix multiplication
         * 
         * @param rhs Matrix to multiply by
         * 
         * @pre Calling matrix must have equal number of columns to the passed matrix's number of rows
         * @return Resulting matrix with the number of rows of the calling matrix and the number of columns of the passed matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator*( const DenseMatrix< T > & rhs ) const;
        
        /**
         * @brief Matrix multiplication
         * 
         * @param rhs Matrix to multiply by
         * 
         * @pre Calling matrix must have equal number of columns to the passed matrix's number of rows
         * @return Resulting matrix with the number of rows of the calling matrix and the number of columns of the passed matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator*( const LowerTriangularMatrix< T > & rhs ) const;
        
        /**
         * @brief Matrix multiplication
         * 
         * @param rhs Matrix to multiply by
         * 
         * @pre Calling matrix must have equal number of columns to the passed matrix's number of rows
         * @return Resulting matrix with the number of rows of the calling matrix and the number of columns of the passed matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        UpperTriangularMatrix< T > operator*( const UpperTriangularMatrix< T > & rhs ) const;
        
        /**
         * @brief Matrix multiplication
         * 
         * @param rhs Matrix to multiply by
         * 
         * @pre Calling matrix must have equal number of columns to the passed matrix's number of rows
         * @return Resulting matrix with the number of rows of the calling matrix and the number of columns of the passed matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        UpperTriangularMatrix< T > operator*( const DiagonalMatrix< T > & rhs ) const;

        /**
         * @brief Matrix multiplication
         * 
         * @param rhs Matrix to multiply by
         * 
         * @pre Calling matrix must have equal number of columns to the passed matrix's number of rows
         * @return Resulting matrix with the number of rows of the calling matrix and the number of columns of the passed matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator*( const SymmetricMatrix< T > & rhs ) const;

        /**
         * @brief Matrix multiplication
         * 
         * @param rhs Matrix to multiply by
         * 
         * @pre Calling matrix must have equal number of columns to the passed matrix's number of rows
         * @return Resulting matrix with the number of rows of the calling matrix and the number of columns of the passed matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator*( const TridiagonalMatrix< T > & rhs ) const;

        /**
         * @brief Scalar multiplication
         * 
         * @param scalar Scalar to multiply by
         * 
         * @pre None
         * @return Resulting matrix of multiplying each element by the scalar
         */
        UpperTriangularMatrix< T > operator*( const double scalar ) const;
        
        /**
         * @brief Scalar multiplication
         * 
         * @param scalar Scalar to multiply by
         * 
         * @pre None
         * @post Calling matrix elements are multiplied by the scalar
         * @return Reference to calling matrix
         */
        UpperTriangularMatrix< T > & operator*=( const double scalar );
        
        /**
         * @brief Tranpose matrix
         * 
         * @pre None
         * @return Resulting matrix created by swapping all elements' rows and columns
         */
        LowerTriangularMatrix< T > transpose() const;

};  /* class UpperTriangularMatrix */

#include "uppertriangularmatrix.hpp"

#endif  /* UPPER_TRIANGULAR_MATRIX_H */