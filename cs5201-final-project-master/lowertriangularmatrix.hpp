/**
 * @file lowertriangularmatrix.hpp
 * @author Caleb Berg (cabr29@mst.edu)
 * @brief Implementation header file for LowerTriangularMatrix class
 * @date 2019-04-03
 */

#ifndef LOWER_TRIANGULAR_MATRIX_HPP
#define LOWER_TRIANGULAR_MATRIX_HPP

#include "lowertriangularmatrix.h"

template< typename T >
LowerTriangularMatrix< T >::LowerTriangularMatrix() :
    data( nullptr ),
    n( 0 )
{

}	/* LowerTriangularMatrix< T >::LowerTriangularMatrix()*/

template< typename T >
LowerTriangularMatrix< T >::LowerTriangularMatrix( const int n ) :
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
        // Lower triangular matrix needs no storage for entries above main diagonal
        this->data[ i ] = Vector< T >( i + 1 );
    }

}	/* LowerTriangularMatrix< T >::LowerTriangularMatrix() */

template< typename T >
LowerTriangularMatrix< T >::LowerTriangularMatrix( const LowerTriangularMatrix< T > & source ) :
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

}	/* LowerTriangularMatrix< T >::LowerTriangularMatrix() */

template< typename T >
LowerTriangularMatrix< T >::LowerTriangularMatrix( const DiagonalMatrix< T > & source ) :
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
            this->data[ i ] = Vector< T >( i + 1 );
            this->data[ i ][ i ] = source.get( i, i );
        }
    }

}   /* LowerTriangularMatrix< T >::LowerTriangularMatrix() */

template< typename T >
LowerTriangularMatrix< T >::~LowerTriangularMatrix()
{
    // Check if destruction is necessary
    if( this->data != nullptr )
    {
        delete[] this->data;
        this->n = 0;
    }

}   /* LowerTriangularMatrix< T >::~LowerTriangularMatrix() */

template< typename T >
LowerTriangularMatrix< T > & LowerTriangularMatrix< T >::operator=( const LowerTriangularMatrix< T > & source )
{
    // Copy Swap Idiom
    LowerTriangularMatrix< T > sourceCopy( source );
    swap( this->data, sourceCopy.data );
    swap( this->n, sourceCopy.n );
    return( *this );

}   /* LowerTriangularMatrix< T >::operator=() */

template< typename T >
T LowerTriangularMatrix< T >::get( const int row, const int col ) const
{
    // Error check
    if( ( row < 0 ) || ( row >= this->n ) || ( col < 0 ) || ( col >= this->n ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix indices" );
    }
    
    // Check if above main diagonal
    return( col > row ? 0 : this->data[ row ][ col ] );

}	/* LowerTriangularMatrix< T >::get() */

template< typename T >
void LowerTriangularMatrix< T >::set( const int row, const int col, const T & val )
{
    // Error check
    if( ( row < 0 ) || ( row >= this->n ) || ( col < 0 ) || ( col >= this->n ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix indices" );
    }
    else if( col > row )
    {
        // throw Error< INDEX_ERROR >( "Invalid index to set in lower triangular matrix" );
    }
    else
    {
        this->data[ row ][ col ] = val;
    }

}	/* LowerTriangularMatrix< T >::set() */

template< typename T >
Vector< T > LowerTriangularMatrix< T >::operator*( const Vector< T > & rhs ) const
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
        for( int j = 0; j <= i ; ++j )
        {
            sum += ( this->data[ i ][ j ] * rhs[ j ] );
        }
        result[ i ] = sum;
    }

    return( result );

}	/* LowerTriangularMatrix< T >::operator*() */

template< typename T >
int LowerTriangularMatrix< T >::getColumns() const
{
    return( this->n );

}   /* LowerTriangularMatrix< T >::getColumns() */

template< typename T >
int LowerTriangularMatrix< T >::getRows() const
{
    return( this->n );

}   /* LowerTriangularMatrix< T >::getRows() */

template< typename T >
bool LowerTriangularMatrix< T >::isLowerTriangular() const
{
    return( true );

}   /* LowerTriangularMatrix< T >::isLowerTriangular() */

template< typename T >
bool LowerTriangularMatrix< T >::isUpperTriangular() const
{
    bool isUpper = true;
    int row = 0;
    while( isUpper && ( row < this->getRows() ) )
    {
        int col = 0;
        while( isUpper && ( col < row ) )
        {
            if( this->get( row, col ) != 0 )
            {
                isUpper = false;
            }
            ++col;
        }
        ++row;
    }
    return( isUpper );

}   /* LowerTriangularMatrix< T >::isUpperTriangular() */

template< typename T >
bool LowerTriangularMatrix< T >::isDiagonal() const
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

}   /* LowerTriangularMatrix< T >::isDiagonal() */

template< typename T >
bool LowerTriangularMatrix< T >::isSymmetric() const
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

}   /* LowerTriangularMatrix< T >::isSymmetric() */

template< typename T >
bool LowerTriangularMatrix< T >::isTridiagonal() const
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

}   /* LowerTriangularMatrix< T >::isTridiagonal() */

template< typename T >
DenseMatrix< T > LowerTriangularMatrix< T >::operator+( const DenseMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* LowerTriangularMatrix< T >::operator+() */

template< typename T >
LowerTriangularMatrix< T > LowerTriangularMatrix< T >::operator+( const LowerTriangularMatrix< T > & rhs ) const
{
    LowerTriangularMatrix< T > result( *this );
    return( result += rhs );

}   /* LowerTriangularMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > LowerTriangularMatrix< T >::operator+( const UpperTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* LowerTriangularMatrix< T >::operator+() */

template< typename T >
LowerTriangularMatrix< T > LowerTriangularMatrix< T >::operator+( const DiagonalMatrix< T > & rhs ) const
{
    LowerTriangularMatrix< T > result( *this );
    return( result += rhs );

}   /* LowerTriangularMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > LowerTriangularMatrix< T >::operator+( const SymmetricMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* LowerTriangularMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > LowerTriangularMatrix< T >::operator+( const TridiagonalMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* LowerTriangularMatrix< T >::operator+() */

template< typename T >
LowerTriangularMatrix< T > & LowerTriangularMatrix< T >::operator+=( const DiagonalMatrix< T > & rhs )
{
    // Error check
    if( ( this->n != rhs.getRows() ) || ( this->n != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be added together" );
    }

    // Add corresponding elements
    for( int i = 0; i < this->n; ++i )
    {
        this->data[ i ][ i ] += rhs.get( i, i );
    }

    return( *this );

}   /* LowerTriangularMatrix< T >::operator+=() */

template< typename T >
LowerTriangularMatrix< T > & LowerTriangularMatrix< T >::operator+=( const LowerTriangularMatrix< T > & rhs )
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

}   /* LowerTriangularMatrix< T >::operator+=() */

template< typename T >
DenseMatrix< T > LowerTriangularMatrix< T >::operator-( const DenseMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* LowerTriangularMatrix< T >::operator-() */

template< typename T >
LowerTriangularMatrix< T > LowerTriangularMatrix< T >::operator-( const LowerTriangularMatrix< T > & rhs ) const
{
    LowerTriangularMatrix< T > result( *this );
    return( result -= rhs );

}   /* LowerTriangularMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > LowerTriangularMatrix< T >::operator-( const UpperTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* LowerTriangularMatrix< T >::operator-() */

template< typename T >
LowerTriangularMatrix< T > LowerTriangularMatrix< T >::operator-( const DiagonalMatrix< T > & rhs ) const
{
    LowerTriangularMatrix< T > result( *this );
    return( result -= rhs );

}   /* LowerTriangularMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > LowerTriangularMatrix< T >::operator-( const SymmetricMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* LowerTriangularMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > LowerTriangularMatrix< T >::operator-( const TridiagonalMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* LowerTriangularMatrix< T >::operator-() */

template< typename T >
LowerTriangularMatrix< T > & LowerTriangularMatrix< T >::operator-=( const DiagonalMatrix< T > & rhs )
{
    // Error check
    if( ( this->n != rhs.getRows() ) || ( this->n != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be subtracted together" );
    }

    // Subtract corresponding elements
    for( int i = 0; i < this->n; ++i )
    {
        this->data[ i ][ i ] -= rhs.get( i, i );
    }

    return( *this );

}   /* LowerTriangularMatrix< T >::operator-=() */

template< typename T >
LowerTriangularMatrix< T > & LowerTriangularMatrix< T >::operator-=( const LowerTriangularMatrix< T > & rhs )
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

}   /* LowerTriangularMatrix< T >::operator-=() */

template< typename T >
LowerTriangularMatrix< T > LowerTriangularMatrix< T >::operator-() const
{
    LowerTriangularMatrix< T > result( *this );
    
    // Negate rows
    for( int i = 0; i < result.getRows(); ++i )
    {
        result.data[ i ] = -( result.data[ i ] );
    }
    return( result );

}   /* LowerTriangularMatrix< T >::operator-=() */

template< typename T >
DenseMatrix< T > LowerTriangularMatrix< T >::operator*( const DenseMatrix< T > & rhs ) const
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
            for( int k = 0; k <= i; ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );

}   /* LowerTriangularMatrix< T >::operator*() */

template< typename T >
LowerTriangularMatrix< T > LowerTriangularMatrix< T >::operator*( const LowerTriangularMatrix< T > & rhs ) const
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
        for( int j = 0; j <= i; ++j )
        {
            T sum = 0;
            for( int k = j; k < rhs.getRows(); ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );

}   /* LowerTriangularMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > LowerTriangularMatrix< T >::operator*( const UpperTriangularMatrix< T > & rhs ) const
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
            for( int k = 0; k <= j; ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );

}   /* LowerTriangularMatrix< T >::operator*() */

template< typename T >
LowerTriangularMatrix< T > LowerTriangularMatrix< T >::operator*( const DiagonalMatrix< T > & rhs ) const
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
        for( int j = 0; j <= i; ++j )
        {
            result.set( i, j, this->get( i, j ) * rhs.get( j, j ) );
        }
    }
    return( result );
}   /* LowerTriangularMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > LowerTriangularMatrix< T >::operator*( const SymmetricMatrix< T > & rhs ) const
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
            for( int k = 0; k <= j; ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );
}   /* LowerTriangularMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > LowerTriangularMatrix< T >::operator*( const TridiagonalMatrix< T > & rhs ) const
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
            for( int k = 0; k <= j; ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );
}   /* LowerTriangularMatrix< T >::operator*() */

template< typename T >
LowerTriangularMatrix< T > LowerTriangularMatrix< T >::operator*( const double scalar ) const
{
    LowerTriangularMatrix< T > result( *this );
    return( result *= scalar );

}   /* LowerTriangularMatrix< T >::operator*() */

template< typename T >
LowerTriangularMatrix< T > & LowerTriangularMatrix< T >::operator*=( const double scalar )
{
    for( int i = 0; i < this->n; ++i )
    {
        this->data[ i ] *= scalar;
    }
    return( *this );

}   /* LowerTriangularMatrix< T >::operator*=() */

template< typename T >
UpperTriangularMatrix< T > LowerTriangularMatrix< T >::transpose() const
{
    UpperTriangularMatrix< T > result( this->n );
    for( int i = 0; i < result.getRows(); ++i )
    {
        for( int j = i; j < result.getColumns(); ++j )
        {
            result.set( i, j, this->data[ j ][ i ] );
        }
    }
    return( result );

}   /* LowerTriangularMatrix< T >::transpose() */


#endif  /* LOWER_TRIANGULAR_MATRIX_HPP */