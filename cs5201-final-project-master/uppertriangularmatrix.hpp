/**
 * @file uppertriangularmatrix.hpp
 * @author Caleb Berg (cabr29@mst.edu)
 * @brief Implementation header file for UpperTriangularMatrix class
 * @date 2019-04-03
 */

#ifndef UPPER_TRIANGULAR_MATRIX_HPP
#define UPPER_TRIANGULAR_MATRIX_HPP

#include "uppertriangularmatrix.h"

template< typename T >
UpperTriangularMatrix< T >::UpperTriangularMatrix() :
    data( nullptr ),
    n( 0 )
{

}	/* UpperTriangularMatrix< T >::UpperTriangularMatrix()*/

template< typename T >
UpperTriangularMatrix< T >::UpperTriangularMatrix( const int n ) :
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
        // Upper triangular matrix needs no storage for entries below main diagonal
        this->data[ i ] = Vector< T >( this->n - i );
    }

}	/* UpperTriangularMatrix< T >::UpperTriangularMatrix() */

template< typename T >
UpperTriangularMatrix< T >::UpperTriangularMatrix( const UpperTriangularMatrix< T > & source ) :
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

}	/* UpperTriangularMatrix< T >::UpperTriangularMatrix() */

template< typename T >
UpperTriangularMatrix< T >::UpperTriangularMatrix( const DiagonalMatrix< T > & source ) :
    data( nullptr ),
    n( source.getRows() )
{
    // Check if any deep copying needs to be done (dimensions > 0)
    if( this->n > 0 )
    {
        // Copy data
        this->data = new Vector< T >[ this->n ]();
        for( int i = 0; i < this->n; ++i )
        {
            this->data[ i ] = Vector< T >( this->n - i );
            this->set( i, i, source.get( i, i ) );
        }
    }

}   /* UpperTriangularMatrix< T >::UpperTriangularMatrix() */

template< typename T >
UpperTriangularMatrix< T >::~UpperTriangularMatrix()
{
    // Check if destruction is necessary
    if( this->data != nullptr )
    {
        delete[] this->data;
        this->n = 0;
    }

}   /* UpperTriangularMatrix< T >::~UpperTriangularMatrix() */

template< typename T >
UpperTriangularMatrix< T > & UpperTriangularMatrix< T >::operator=( const UpperTriangularMatrix< T > & source )
{
    // Copy Swap Idiom
    UpperTriangularMatrix< T > sourceCopy( source );
    swap( this->data, sourceCopy.data );
    swap( this->n, sourceCopy.n );
    return( *this );

}   /* UpperTriangularMatrix< T >::operator=() */

template< typename T >
T UpperTriangularMatrix< T >::get( const int row, const int col ) const
{
    // Error check
    if( ( row < 0 ) || ( row >= this->n ) || ( col < 0 ) || ( col >= this->n ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix indices" );
    }
    else
    
    // Check if below main diagonal
    return( col < row ? 0 : this->data[ row ][ col - row ] );

}	/* UpperTriangularMatrix< T >::get() */

template< typename T >
void UpperTriangularMatrix< T >::set( const int row, const int col, const T & val )
{
    // Error check
    if( ( row < 0 ) || ( row >= this->n ) || ( col < 0 ) || ( col >= this->n ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix indices" );
    }
    else if( col < row )
    {
        // throw Error< INDEX_ERROR >( "Invalid index to set in upper triangular matrix" );
    }
    else
    {
        this->data[ row ][ col - row ] = val;
    }

}	/* UpperTriangularMatrix< T >::set() */

template< typename T >
Vector< T > UpperTriangularMatrix< T >::operator*( const Vector< T > & rhs ) const
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
        T sum = 0;
        for( int j = i; j < this->n ; ++j )
        {
            sum += ( this->get( i, j ) * rhs[ j ] );
        }
        result[ i ] = sum;
    }

    return( result );

}	/* UpperTriangularMatrix< T >::operator*() */

template< typename T >
int UpperTriangularMatrix< T >::getColumns() const
{
    return( this->n );

}   /* UpperTriangularMatrix< T >::getColumns() */

template< typename T >
int UpperTriangularMatrix< T >::getRows() const
{
    return( this->n );

}   /* UpperTriangularMatrix< T >::getRows() */

template< typename T >
bool UpperTriangularMatrix< T >::isLowerTriangular() const
{
    bool isLower = true;
    int row = 0;
    while( isLower && ( row < this->getRows() ) )
    {
        int col = ( row + 1 );
        while( isLower && ( col < this->getColumns() ) )
        {
            if( this->get( row, col ) != 0 )
            {
                isLower = false;
            }
            ++col;
        }
        ++row;
    }
    return( isLower );

}   /* UpperTriangularMatrix< T >::isLowerTriangular() */

template< typename T >
bool UpperTriangularMatrix< T >::isUpperTriangular() const
{
    return( true );

}   /* UpperTriangularMatrix< T >::isUpperTriangular() */

template< typename T >
bool UpperTriangularMatrix< T >::isDiagonal() const
{
    bool isDiagonal = true;
    int row = 0;
    while( isDiagonal && ( row < this->getRows() ) )
    {
        int col = 0;
        while( isDiagonal && ( col < this->getColumns() ) )
        {
            if( ( row != col ) && ( this->get( row, col ) != 0 ) )
            {
                isDiagonal = false;
            }
            ++col;
        }
        ++row;
    }
    return( isDiagonal );

}   /* UpperTriangularMatrix< T >::isDiagonal() */

template< typename T >
bool UpperTriangularMatrix< T >::isSymmetric() const
{
    bool isSymmetric = true;
    int row = 0;
    while( isSymmetric && ( row < this->getRows() ) )
    {
        int col = ( row + 1 );
        while( isSymmetric && ( col < this->getColumns() ) )
        {
            if( this->get( row, col ) != this->get( col, row ) )
            {
                isSymmetric = false;
            }
            ++col;
        }
        ++row;
    }
    return( isSymmetric );

}   /* UpperTriangularMatrix< T >::isSymmetric() */

template< typename T >
bool UpperTriangularMatrix< T >::isTridiagonal() const
{
    bool isTridiagonal = true;
    int row = 0;
    while( isTridiagonal && ( row < this->getRows() ) )
    {
        int col = 0;
        while( isTridiagonal && ( col < this->getColumns() ) )
        {
            if( ( ( row > col + 1 ) || ( col > row + 1 ) ) && ( this->get( row, col ) != 0 ) )
            {
                isTridiagonal = false;
            }
            ++col;
        }
        ++row;
    }
    return( isTridiagonal );

}   /* UpperTriangularMatrix< T >::isTridiagonal() */

template< typename T >
DenseMatrix< T > UpperTriangularMatrix< T >::operator+( const DenseMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* UpperTriangularMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > UpperTriangularMatrix< T >::operator+( const LowerTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* UpperTriangularMatrix< T >::operator+() */

template< typename T >
UpperTriangularMatrix< T > UpperTriangularMatrix< T >::operator+( const UpperTriangularMatrix< T > & rhs ) const
{
    UpperTriangularMatrix< T > result( *this );
    return( result += rhs );

}   /* UpperTriangularMatrix< T >::operator+() */

template< typename T >
UpperTriangularMatrix< T > UpperTriangularMatrix< T >::operator+( const DiagonalMatrix< T > & rhs ) const
{
    UpperTriangularMatrix< T > result( *this );
    return( result += rhs );

}   /* UpperTriangularMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > UpperTriangularMatrix< T >::operator+( const SymmetricMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* UpperTriangularMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > UpperTriangularMatrix< T >::operator+( const TridiagonalMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* UpperTriangularMatrix< T >::operator+() */

template< typename T >
UpperTriangularMatrix< T > & UpperTriangularMatrix< T >::operator+=( const DiagonalMatrix< T > & rhs )
{
    // Error check
    if( ( this->n != rhs.getRows() ) || ( this->n != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be added together" );
    }

    // Add corresponding elements
    for( int i = 0; i < this->n; ++i )
    {
        this->set( i, i, this->get( i, i ) + rhs.get( i, i ) );
    }

    return( *this );

}   /* UpperTriangularMatrix< T >::operator+=() */

template< typename T >
UpperTriangularMatrix< T > & UpperTriangularMatrix< T >::operator+=( const UpperTriangularMatrix< T > & rhs )
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

}   /* UpperTriangularMatrix< T >::operator+=() */

template< typename T >
DenseMatrix< T > UpperTriangularMatrix< T >::operator-( const DenseMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* UpperTriangularMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > UpperTriangularMatrix< T >::operator-( const LowerTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* UpperTriangularMatrix< T >::operator-() */

template< typename T >
UpperTriangularMatrix< T > UpperTriangularMatrix< T >::operator-( const UpperTriangularMatrix< T > & rhs ) const
{
    UpperTriangularMatrix< T > result( *this );
    return( result -= rhs );

}   /* UpperTriangularMatrix< T >::operator-() */

template< typename T >
UpperTriangularMatrix< T > UpperTriangularMatrix< T >::operator-( const DiagonalMatrix< T > & rhs ) const
{
    UpperTriangularMatrix< T > result( *this );
    return( result -= rhs );

}   /* UpperTriangularMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > UpperTriangularMatrix< T >::operator-( const SymmetricMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* UpperTriangularMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > UpperTriangularMatrix< T >::operator-( const TridiagonalMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* UpperTriangularMatrix< T >::operator-() */

template< typename T >
UpperTriangularMatrix< T > & UpperTriangularMatrix< T >::operator-=( const DiagonalMatrix< T > & rhs )
{
    // Error check
    if( ( this->n != rhs.getRows() ) || ( this->n != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be subtracted together" );
    }

    // Subtract corresponding elements
    for( int i = 0; i < this->n; ++i )
    {
        this->set( i, i, this->get( i, i ) - rhs.get( i, i ) );
    }

    return( *this );

}   /* UpperTriangularMatrix< T >::operator-=() */

template< typename T >
UpperTriangularMatrix< T > & UpperTriangularMatrix< T >::operator-=( const UpperTriangularMatrix< T > & rhs )
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

}   /* UpperTriangularMatrix< T >::operator-=() */

template< typename T >
UpperTriangularMatrix< T > UpperTriangularMatrix< T >::operator-() const
{
    UpperTriangularMatrix< T > result( *this );
    
    // Negate rows
    for( int i = 0; i < result.getRows(); ++i )
    {
        result.data[ i ] = -( result.data[ i ] );
    }
    return( result );

}   /* UpperTriangularMatrix< T >::operator-=() */

template< typename T >
DenseMatrix< T > UpperTriangularMatrix< T >::operator*( const DenseMatrix< T > & rhs ) const
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
            T sum = 0;
            for( int k = i; k < this->getColumns(); ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );

}   /* UpperTriangularMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > UpperTriangularMatrix< T >::operator*( const LowerTriangularMatrix< T > & rhs ) const
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
            T sum = 0;
            for( int k = i; k < rhs.getRows(); ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );

}   /* UpperTriangularMatrix< T >::operator*() */

template< typename T >
UpperTriangularMatrix< T > UpperTriangularMatrix< T >::operator*( const UpperTriangularMatrix< T > & rhs ) const
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
            T sum = 0;
            for( int k = 0; k <= j; ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );

}   /* UpperTriangularMatrix< T >::operator*() */

template< typename T >
UpperTriangularMatrix< T > UpperTriangularMatrix< T >::operator*( const DiagonalMatrix< T > & rhs ) const
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
        for( int j = i; j < result.getColumns(); ++j )
        {
            result.set( i, j, this->get( i, j ) * rhs.get( j, j ) );
        }
    }
    return( result );

}   /* UpperTriangularMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > UpperTriangularMatrix< T >::operator*( const SymmetricMatrix< T > & rhs ) const
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
            T sum = 0;
            for( int k = i; k < this->getColumns(); ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );

}   /* UpperTriangularMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > UpperTriangularMatrix< T >::operator*( const TridiagonalMatrix< T > & rhs ) const
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
            T sum = 0;
            for( int k = i; k < this->getColumns(); ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );

}   /* UpperTriangularMatrix< T >::operator*() */

template< typename T >
UpperTriangularMatrix< T > UpperTriangularMatrix< T >::operator*( const double scalar ) const
{
    UpperTriangularMatrix< T > result( *this );
    return( result *= scalar );

}   /* UpperTriangularMatrix< T >::operator*() */

template< typename T >
UpperTriangularMatrix< T > & UpperTriangularMatrix< T >::operator*=( const double scalar )
{
    for( int i = 0; i < this->n; ++i )
    {
        this->data[ i ] *= scalar;
    }
    return( *this );

}   /* UpperTriangularMatrix< T >::operator*=() */

template< typename T >
LowerTriangularMatrix< T > UpperTriangularMatrix< T >::transpose() const
{
    LowerTriangularMatrix< T > result( this->n );
    for( int i = 0; i < result.getRows(); ++i )
    {
        for( int j = 0; j <= i; ++j )
        {
            result.set( i, j, this->get( j, i ) );
        }
    }
    return( result );

}   /* UpperTriangularMatrix< T >::transpose() */

#endif  /* UPPER_TRIANGULAR_MATRIX_HPP */