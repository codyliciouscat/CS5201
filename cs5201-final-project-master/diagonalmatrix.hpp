/**
 * @file diagonalmatrix.hpp
 * @author Caleb Berg (cabr29@mst.edu)
 * @brief Implementation header file for DiagonalMatrix class
 * @date 2019-04-03
 */

#ifndef DIAGONAL_MATRIX_HPP
#define DIAGONAL_MATRIX_HPP

#include "diagonalmatrix.h"

template< typename T >
DiagonalMatrix< T >::DiagonalMatrix() :
    data( nullptr ),
    n( 0 )
{

}	/* DiagonalMatrix< T >::DiagonalMatrix()*/

template< typename T >
DiagonalMatrix< T >::DiagonalMatrix( const int n ) :
    n( n )
{
    // Error check
    if( this->n < 0 )
    {
        throw Error< PARAMETER_ERROR >( "Invalid matrix dimensions" );
    }

    // Construct matrix
    this->data = new Vector< T >[ this->n ]();
    for( int i = 0; i < this->n ; ++i )
    {
        // Diagonal matrix has 1 element per row
        this->data[ i ] = Vector< T >( 1 );
    }

}	/* DiagonalMatrix< T >::DiagonalMatrix() */

template< typename T >
DiagonalMatrix< T >::DiagonalMatrix( const DiagonalMatrix< T > & source ) :
    data( nullptr ),
    n( source.n )
{
    // Check if any deep copying needs to be done (dimensions > 0)
    if( this->n > 0 )
    {
        // Copy data
        this->data = new Vector< T >[ this->n ]();
        for( int i = 0; i < this->n; ++i )
        {
            data[ i ] = source.data[ i ];
        }
    }

}	/* DiagonalMatrix< T >::DiagonalMatrix() */

template< typename T >
DiagonalMatrix< T >::~DiagonalMatrix()
{
    // Check if destruction is necessary
    if( this->data != nullptr )
    {
        delete[] this->data;
        this->n = 0;
    }

}   /* DiagonalMatrix< T >::~DiagonalMatrix() */

template< typename T >
DiagonalMatrix< T > & DiagonalMatrix< T >::operator=( const DiagonalMatrix< T > & source )
{
    // Copy Swap Idiom
    DiagonalMatrix< T > sourceCopy( source );
    swap( this->data, sourceCopy.data );
    swap( this->n, sourceCopy.n );
    return( *this );

}   /* DiagonalMatrix< T >::operator=() */

template< typename T >
T DiagonalMatrix< T >::get( const int row, const int col ) const
{
    // Error check
    if( ( row < 0 ) || ( row >= this->n ) || ( col < 0 ) || ( col >= this->n ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix indices" );
    }
    
    // Check if not diagonal
    return( col != row ? 0 : this->data[ row ][ 0 ] );

}	/* DiagonalMatrix< T >::get() */

template< typename T >
void DiagonalMatrix< T >::set( const int row, const int col, const T & val )
{
    // Error check
    if( ( row < 0 ) || ( row >= this->n ) || ( col < 0 ) || ( col >= this->n ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix indices" );
    }
    else if( col != row )
    {
        // throw Error< INDEX_ERROR >( "Invalid index to set in diagonal matrix" );
    }
    else
    {
        this->data[ row ][ 0 ] = val;
    }

}	/* DiagonalMatrix< T >::set() */

template< typename T >
Vector< T > DiagonalMatrix< T >::operator*( const Vector< T > & rhs ) const
{
    // Error check
    if( this->n != rhs.size() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrix is unable to be multiplied by vector" );
    }

    Vector< T > result( this->n );

    // Multiply matrix and vector
    for( int i = 0; i < result.size(); ++i )
    {
        result[ i ] = ( this->get( i, i ) * rhs[ i ] );
    }

    return( result );

}	/* DiagonalMatrix< T >::operator*() */

template< typename T >
int DiagonalMatrix< T >::getColumns() const
{
    return( this->n );

}   /* DiagonalMatrix< T >::getColumns() */

template< typename T >
int DiagonalMatrix< T >::getRows() const
{
    return( this->n );

}   /* DiagonalMatrix< T >::getRows() */

template< typename T >
bool DiagonalMatrix< T >::isLowerTriangular() const
{
    return( true );

}   /* DiagonalMatrix< T >::isLowerTriangular() */

template< typename T >
bool DiagonalMatrix< T >::isUpperTriangular() const
{
    return( true );

}   /* DiagonalMatrix< T >::isUpperTriangular() */

template< typename T >
bool DiagonalMatrix< T >::isDiagonal() const
{
    return( true );

}   /* DiagonalMatrix< T >::isDiagonal() */

template< typename T >
bool DiagonalMatrix< T >::isSymmetric() const
{
    return( true );

}   /* DiagonalMatrix< T >::isSymmetric() */

template< typename T >
bool DiagonalMatrix< T >::isTridiagonal() const
{
    return( true );

}   /* DiagonalMatrix< T >::isTridiagonal() */

template< typename T >
DenseMatrix< T > DiagonalMatrix< T >::operator+( const DenseMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* DiagonalMatrix< T >::operator+() */

template< typename T >
LowerTriangularMatrix< T > DiagonalMatrix< T >::operator+( const LowerTriangularMatrix< T > & rhs ) const
{
    LowerTriangularMatrix< T > result( *this );
    return( result += rhs );

}   /* DiagonalMatrix< T >::operator+() */

template< typename T >
UpperTriangularMatrix< T > DiagonalMatrix< T >::operator+( const UpperTriangularMatrix< T > & rhs ) const
{
    UpperTriangularMatrix< T > result( *this );
    return( result += rhs );

}   /* DiagonalMatrix< T >::operator+() */

template< typename T >
DiagonalMatrix< T > DiagonalMatrix< T >::operator+( const DiagonalMatrix< T > & rhs ) const
{
    DiagonalMatrix< T > result( *this );
    return( result += rhs );

}   /* DiagonalMatrix< T >::operator+() */

template< typename T >
SymmetricMatrix< T > DiagonalMatrix< T >::operator+( const SymmetricMatrix< T > & rhs ) const
{
    SymmetricMatrix< T > result( *this );
    return( result += rhs );

}   /* DiagonalMatrix< T >::operator+() */

template< typename T >
TridiagonalMatrix< T > DiagonalMatrix< T >::operator+( const TridiagonalMatrix< T > & rhs ) const
{
    TridiagonalMatrix< T > result( *this );
    return( result += rhs );

}   /* DiagonalMatrix< T >::operator+() */

template< typename T >
DiagonalMatrix< T > & DiagonalMatrix< T >::operator+=( const DiagonalMatrix< T > & rhs )
{
    // Error check
    if( ( this->n != rhs.getRows() ) || ( this->n != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be added together" );
    }

    // Add corresponding rows
    for( int i = 0; i < this->n; ++i )
    {
        this->data[ i ] += rhs.data[ i ];
    }

    return( *this );

}   /* DiagonalMatrix< T >::operator+=() */

template< typename T >
DenseMatrix< T > DiagonalMatrix< T >::operator-( const DenseMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* DiagonalMatrix< T >::operator-() */

template< typename T >
LowerTriangularMatrix< T > DiagonalMatrix< T >::operator-( const LowerTriangularMatrix< T > & rhs ) const
{
    LowerTriangularMatrix< T > result( *this );
    return( result -= rhs );

}   /* DiagonalMatrix< T >::operator-() */

template< typename T >
UpperTriangularMatrix< T > DiagonalMatrix< T >::operator-( const UpperTriangularMatrix< T > & rhs ) const
{
    UpperTriangularMatrix< T > result( *this );
    return( result -= rhs );

}   /* DiagonalMatrix< T >::operator-() */

template< typename T >
DiagonalMatrix< T > DiagonalMatrix< T >::operator-( const DiagonalMatrix< T > & rhs ) const
{
    DiagonalMatrix< T > result( *this );
    return( result -= rhs );

}   /* DiagonalMatrix< T >::operator-() */

template< typename T >
SymmetricMatrix< T > DiagonalMatrix< T >::operator-( const SymmetricMatrix< T > & rhs ) const
{
    SymmetricMatrix< T > result( *this );
    return( result -= rhs );

}   /* DiagonalMatrix< T >::operator-() */

template< typename T >
TridiagonalMatrix< T > DiagonalMatrix< T >::operator-( const TridiagonalMatrix< T > & rhs ) const
{
    TridiagonalMatrix< T > result( *this );
    return( result -= rhs );

}   /* DiagonalMatrix< T >::operator-() */

template< typename T >
DiagonalMatrix< T > & DiagonalMatrix< T >::operator-=( const DiagonalMatrix< T > & rhs )
{
    // Error check
    if( ( this->n != rhs.getRows() ) || ( this->n != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be subtracted together" );
    }

    // Subtract corresponding rows
    for( int i = 0; i < this->n; ++i )
    {
        this->data[ i ] -= rhs.data[ i ];
    }

    return( *this );

}   /* DiagonalMatrix< T >::operator-=() */

template< typename T >
DiagonalMatrix< T > DiagonalMatrix< T >::operator-() const
{
    DiagonalMatrix< T > result( *this );
    
    // Negate rows
    for( int i = 0; i < result.getRows(); ++i )
    {
        result.data[ i ] = -( result.data[ i ] );
    }
    return( result );

}   /* DiagonalMatrix< T >::operator-=() */

template< typename T >
DenseMatrix< T > DiagonalMatrix< T >::operator*( const DenseMatrix< T > & rhs ) const
{
    // Error check
    if( this->n != rhs.getRows() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be multiplied together" );
    }

    DenseMatrix< T > result( this->n, rhs.getColumns() );

    // Multiply matrices
    for( int i = 0; i < result.getRows(); ++i )
    {
        for( int j = 0; j < result.getColumns(); ++j )
        {
            result.set( i, j, this->get( i, i ) * rhs.get( i, j ) );
        }
    }
    return( result );
}   /* DiagonalMatrix< T >::operator*() */

template< typename T >
LowerTriangularMatrix< T > DiagonalMatrix< T >::operator*( const LowerTriangularMatrix< T > & rhs ) const
{
    // Error check
    if( this->n != rhs.getRows() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be multiplied together" );
    }

    LowerTriangularMatrix< T > result( this->n );

    // Multiply matrices
    for( int i = 0; i < result.getRows(); ++i )
    {
        for( int j = 0; j < result.getColumns(); ++j )
        {
            result.set( i, j, this->get( i, i ) * rhs.get( i, j ) );
        }
    }
    return( result );

}   /* DiagonalMatrix< T >::operator*() */

template< typename T >
UpperTriangularMatrix< T > DiagonalMatrix< T >::operator*( const UpperTriangularMatrix< T > & rhs ) const
{
    // Error check
    if( this->n != rhs.getRows() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be multiplied together" );
    }

    UpperTriangularMatrix< T > result( this->n );

    // Multiply matrices
    for( int i = 0; i < result.getRows(); ++i )
    {
        for( int j = 0; j < result.getColumns(); ++j )
        {
            result.set( i, j, this->get( i, i ) * rhs.get( i, j ) );
        }
    }
    return( result );

}   /* DiagonalMatrix< T >::operator*() */

template< typename T >
DiagonalMatrix< T > DiagonalMatrix< T >::operator*( const DiagonalMatrix< T > & rhs ) const
{
    // Error check
    if( this->n != rhs.getRows() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be multiplied together" );
    }

    DiagonalMatrix< T > result( this->n );

    // Multiply matrices
    for( int i = 0; i < result.getRows(); ++i )
    {
        result.set( i, i, this->get( i, i ) * rhs.get( i, i ) );
    }
    return( result );

}   /* DiagonalMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > DiagonalMatrix< T >::operator*( const SymmetricMatrix< T > & rhs ) const
{
    // Error check
    if( this->n != rhs.getRows() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be multiplied together" );
    }

    DenseMatrix< T > result( this->n, rhs.getColumns() );

    // Multiply matrices
    for( int i = 0; i < result.getRows(); ++i )
    {
        for( int j = 0; j < result.getColumns(); ++j )
        {
            result.set( i, j, this->get( i, i ) * rhs.get( i, j ) );
        }
    }
    return( result );

}   /* DiagonalMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > DiagonalMatrix< T >::operator*( const TridiagonalMatrix< T > & rhs ) const
{
    // Error check
    if( this->n != rhs.getRows() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be multiplied together" );
    }

    DenseMatrix< T > result( this->n, rhs.getColumns() );

    // Multiply matrices
    for( int i = 0; i < result.getRows(); ++i )
    {
        for( int j = 0; j < result.getColumns(); ++j )
        {
            result.set( i, j, this->get( i, i ) * rhs.get( i, j ) );
        }
    }
    return( result );
    
}   /* DiagonalMatrix< T >::operator*() */

template< typename T >
DiagonalMatrix< T > DiagonalMatrix< T >::operator*( const double scalar ) const
{
    DiagonalMatrix< T > result( *this );
    return( result *= scalar );

}   /* DiagonalMatrix< T >::operator*() */

template< typename T >
DiagonalMatrix< T > & DiagonalMatrix< T >::operator*=( const double scalar )
{
    for( int i = 0; i < this->n; ++i )
    {
        this->data[ i ] *= scalar;
    }
    return( *this );

}   /* DiagonalMatrix< T >::operator*=() */

template< typename T >
DiagonalMatrix< T > DiagonalMatrix< T >::transpose() const
{
    DiagonalMatrix< T > result( *this );
    return( result );

}   /* DiagonalMatrix< T >::transpose() */


#endif  /* DIAGONAL_MATRIX_HPP */