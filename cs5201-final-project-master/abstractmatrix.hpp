/**
 * @file abstractmatrix.hpp
 * @author Caleb Berg (cabr29@mst.edu)
 * @brief Implementation header file for AbstractMatrix class
 * @date 2019-04-03
 */

#ifndef ABSTRACT_MATRIX_HPP
#define ABSTRACT_MATRIX_HPP

#include "abstractmatrix.h"
#include "densematrix.h"

template< typename T >
DenseMatrix< T > AbstractMatrix< T >::operator+( const AbstractMatrix< T > & rhs ) const
{
    // Error check
    if( ( this->getRows() != rhs.getRows() ) || ( this->getColumns() != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be added together" );
    }

    DenseMatrix< T > result( *this );

    // Add corresponding rows
    for( int i = 0; i < this->getRows(); ++i )
    {
        for( int j = i; j < this->getColumns(); ++j )
        {
            result.set( i, j, ( result.get( i, j ) + rhs.get( i, j ) ) );
        }
    }

    return( *this );

}   /* AbstractMatrix::operator+() */

template< typename T >
DenseMatrix< T > AbstractMatrix< T >::operator-( const AbstractMatrix< T > & rhs ) const
{
    // Error check
    if( ( this->getRows() != rhs.getRows() ) || ( this->getColumns() != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be added together" );
    }

    DenseMatrix< T > result( *this );

    // Subtract corresponding rows
    for( int i = 0; i < this->getRows(); ++i )
    {
        for( int j = i; j < this->getColumns(); ++j )
        {
            result.set( i, j, ( result.get( i, j ) - rhs.get( i, j ) ) );
        }
    }

    return( *this );

}   /* AbstractMatrix::operator-() */

template< typename T >
DenseMatrix< T > AbstractMatrix< T >::operator*( const AbstractMatrix< T > & rhs ) const
{
    // Error check
    if( this->getColumns() != rhs.getRows() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be multiplied together" );
    }

    DenseMatrix< T > result( this->getRows(), rhs.getColumns() );

    // Multiply matrices
    for( int i = 0; i < result.getRows(); ++i )
    {
        for( int j = 0; j < result.getColumns(); ++j )
        {
            T sum = 0;
            for( int k = 0; k < this->getRows(); ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );

}   /* AbstractMatrix::operator*() */

template< typename T >
bool AbstractMatrix< T >::operator==( const AbstractMatrix< T > & rhs ) const
{
    bool isEqual = true;

    if( ( this->getRows() != rhs.getRows() ) || ( this->getColumns() != rhs.getColumns() ) )
    {
        isEqual = false;
    }

    int row = 0;
    int col = 0;
    int rows = this->getRows();
    int cols = this->getColumns();
    while( isEqual && ( row < rows ) )
    {
        while( isEqual && ( col < cols ) )
        {
            if( this->get( row, col ) != rhs.get( row, col ) )
            {
                isEqual = false;
            }
            ++col;
        }
        ++row;
    }
    return( isEqual );

}   /* AbstractMatrix< T >::operator==() */

#endif /* ABSTRACT_MATRIX_HPP */