/**
 * @file abstractmatrix.h
 * @author Caleb Berg (cabr29@mst.edu)
 * @brief Header file for AbstractMatrix class
 * @date 2019-04-03
 */

#ifndef ABSTRACT_MATRIX_H
#define ABSTRACT_MATRIX_H

#include "vector.h"
#include "error.h"

template< typename T >
class DenseMatrix;

/**
 * @brief An abstract matrix class to derive from
 */
template< typename T >
class AbstractMatrix
{
    public:
        /**
         * @brief Accessor
         * 
         * @param row Row in matrix
         * @param col Column in matrix
         * 
         * @pre 0 <= row < getRows() and 0 <= col < getColumns()
         * @return Element at matrix location
         */
        virtual T get( const int row, const int col ) const = 0;

        /**
         * @brief Mutator
         * 
         * @param row Row in matrix
         * @param col Column in matrix
         * @param val Value to store in matrix
         * 
         * @pre 0 <= row < getRows() and 0 <= col < getColumns()
         * @post Value is stored in matrix at location (row, col)
         */
        virtual void set( const int row, const int col, const T & val ) = 0;

        /**
         * @brief Vector multiplication
         * 
         * @param rhs Vector to multiply by
         * 
         * @pre getColumns() == rhs.size()
         * @return Resulting vector from multiplication
         */
        virtual Vector< T > operator*( const Vector< T > & rhs ) const = 0;

        /**
         * @brief Get number of columns in matrix
         * 
         * @pre None
         * @return Number of columns in matrix
         */
        virtual int getColumns() const = 0;

        /**
         * @brief Get number of rows in matrix
         * 
         * @pre None
         * @return Number of rows in matrix 
         */
        virtual int getRows() const = 0;

        /**
         * @brief Check if matrix is lower triangular
         * 
         * @pre None
         * @return true if matrix is lower triangular
         */
        virtual bool isLowerTriangular() const = 0;

        /**
         * @brief Check if matrix is upper triangular
         * 
         * @pre None
         * @return true if matrix is upper triangular
         */
        virtual bool isUpperTriangular() const = 0;

        /**
         * @brief Check if matrix is diagonal
         * 
         * @pre None
         * @return true if matrix is diagonal
         */
        virtual bool isDiagonal() const = 0;

        /**
         * @brief Check if matrix is symmetric
         * 
         * @pre None
         * @return true if matrix is symmetric
         */
        virtual bool isSymmetric() const = 0;

        /**
         * @brief Check if matrix is tridiagonal
         * 
         * @pre None
         * @return true if matrix is tridiagonal
         */
        virtual bool isTridiagonal() const = 0;

        /**
         * @brief Matrix addition
         * 
         * @param rhs Matrix to add by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by adding corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator+( const AbstractMatrix< T > & rhs ) const;

        /**
         * @brief Matrix subtraction
         * 
         * @param rhs Matrix to subtract by
         * 
         * @pre Matrices must have the same number of rows and columns
         * @return Result matrix created by subtracting corresponding elements of each matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator-( const AbstractMatrix< T > & rhs ) const;

        /**
         * @brief Matrix multiplication
         * 
         * @param rhs Matrix to multiply by
         * 
         * @pre Calling matrix must have equal number of columns to the passed matrix's number of rows
         * @return Resulting matrix with the number of rows of the calling matrix and the number of columns of the passed matrix
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        DenseMatrix< T > operator*( const AbstractMatrix< T > & rhs ) const;

        /**
         * @brief Matrix comparison
         * 
         * @param rhs Matrix to compare with
         * 
         * @pre None
         * @return true if matrices are the same size with the same elements
         */
        bool operator==( const AbstractMatrix< T > & rhs ) const;

};  /* class AbstractMatrix */

/**
 * @brief Stream insertion
 * 
 * @param os Stream to output into
 * @param matrix Matrix to output
 * 
 * @pre Stream is open and in good health
 * @post Matrix data is outputted into the stream
 * @return Reference to stream
 */
template< typename T >
ostream & operator<<( ostream & os, const AbstractMatrix< T > & matrix )
{
    for( int i = 0; i < matrix.getRows(); ++i )
    {
        for( int j = 0; j < matrix.getColumns(); ++j )
        {
            os << ( j != 0 ? "\t" : "" ) << matrix.get( i, j );
        }
        os << endl;
    }
    return( os );

}   /* operator<<() */

/**
 * @brief Stream extraction
 *
 * @param is Stream to extract from
 * @param matrix Matrix to fill with data
 * 
 * @pre Stream is open and in good health
 * @post Matrix is updated with data from the stream
 * @return Reference to stream
 */
template< typename T >
istream & operator>>( istream & is, AbstractMatrix< T > & matrix )
{
    for( int i = 0; i < matrix.getRows(); ++i )
    {
        for( int j = 0; j < matrix.getColumns(); ++j )
        {
            T val;
            is >> val;
            matrix.set( i, j , val );
        }
    }
    return( is );

}   /* operator>>() */

#include "abstractmatrix.hpp"

#endif  /* MATRIX_INTERFACE_H */