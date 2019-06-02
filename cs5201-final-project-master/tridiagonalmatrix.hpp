/**
 * @file tridiagonalmatrix.hpp
 * @author Caleb Berg (cabr29@mst.edu)
 * @brief Implementation header file for TridiagonalMatrix class
 * @date 2019-04-03
 */

#ifndef TRIDIAGONAL_MATRIX_HPP
#define TRIDIAGONAL_MATRIX_HPP

#include "tridiagonalmatrix.h"

template< typename T >
TridiagonalMatrix< T >::TridiagonalMatrix() :
    data( nullptr ),
    n( 0 )
{

}	/* TridiagonalMatrix< T >::TridiagonalMatrix()*/

template< typename T >
TridiagonalMatrix< T >::TridiagonalMatrix( const int n ) :
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
        // Tridiagonal matrix has 2 elements on first and last rows and 3 on the rest
        if( this->n > 2 )
        {
            this->data[ i ] = Vector< T >( ( ( i == 0 ) || ( i == ( this->n - 1 ) ) ) ? 2 : 3 );
        }
        else
        {
            this->data[ i ] = Vector< T >( this->n );
        }
    }

}	/* TridiagonalMatrix< T >::TridiagonalMatrix() */

template< typename T >
TridiagonalMatrix< T >::TridiagonalMatrix( const TridiagonalMatrix< T > & source ) :
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

}	/* TridiagonalMatrix< T >::TridiagonalMatrix() */

template< typename T >
TridiagonalMatrix< T >::TridiagonalMatrix( const DiagonalMatrix< T > & source ) :
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
            // Tridiagonal matrix has 2 elements on first and last rows and 3 on the rest
            if( this->n > 2 )
            {
                this->data[ i ] = Vector< T >( ( ( i == 0 ) || ( i == ( this->n - 1 ) ) ) ? 2 : 3 );
            }
            else
            {
                this->data[ i ] = Vector< T >( this->n );
            }
            this->set( i, i, source.get( i, i ) );
        }
    }

}   /* TridiagonalMatrix< T >::TridiagonalMatrix() */

template< typename T >
TridiagonalMatrix< T >::~TridiagonalMatrix()
{
    // Check if destruction is necessary
    if( this->data != nullptr )
    {
        delete[] this->data;
        this->n = 0;
    }

}   /* TridiagonalMatrix< T >::~TridiagonalMatrix() */

template< typename T >
TridiagonalMatrix< T > & TridiagonalMatrix< T >::operator=( const TridiagonalMatrix< T > & source )
{
    // Copy Swap Idiom
    TridiagonalMatrix< T > sourceCopy( source );
    swap( this->data, sourceCopy.data );
    swap( this->n, sourceCopy.n );
    return( *this );

}   /* TridiagonalMatrix< T >::operator=() */

template< typename T >
T TridiagonalMatrix< T >::get( const int row, const int col ) const
{
    // Error check
    if( ( row < 0 ) || ( row >= this->n ) || ( col < 0 ) || ( col >= this->n ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix indices" );
    }
    
    return( ( ( row > col + 1 ) || ( col > row + 1 ) ) ? 0 : this->data[ row ][ col - ( row - 1 ) - ( row == 0 ) ] );

}	/* TridiagonalMatrix< T >::get() */

template< typename T >
void TridiagonalMatrix< T >::set( const int row, const int col, const T & val )
{
    // Error check
    if( ( row < 0 ) || ( row >= this->n ) || ( col < 0 ) || ( col >= this->n ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix indices" );
    }
    else if( ( row > col + 1 ) || ( col > row + 1 ) )
    {
        // throw Error< INDEX_ERROR >( "Invalid index to set in diagonal matrix" );
    }
    else
    {
        this->data[ row ][ col - ( row - 1 ) - ( row == 0 ) ] = val;
    }

}	/* TridiagonalMatrix< T >::set() */

template< typename T >
Vector< T > TridiagonalMatrix< T >::operator*( const Vector< T > & rhs ) const
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
        for( int j = 0; j < this->n ; ++j )
        {
            sum += ( this->get( i, j ) * rhs[ j ] );
        }
        result[ i ] = sum;
    }

    return( result );

}	/* TridiagonalMatrix< T >::operator*() */

template< typename T >
int TridiagonalMatrix< T >::getColumns() const
{
    return( this->n );

}   /* TridiagonalMatrix< T >::getColumns() */

template< typename T >
int TridiagonalMatrix< T >::getRows() const
{
    return( this->n );

}   /* TridiagonalMatrix< T >::getRows() */

template< typename T >
bool TridiagonalMatrix< T >::isLowerTriangular() const
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

}   /* TridiagonalMatrix< T >::isLowerTriangular() */

template< typename T >
bool TridiagonalMatrix< T >::isUpperTriangular() const
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

}   /* TridiagonalMatrix< T >::isUpperTriangular() */

template< typename T >
bool TridiagonalMatrix< T >::isDiagonal() const
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

}   /* TridiagonalMatrix< T >::isDiagonal() */

template< typename T >
bool TridiagonalMatrix< T >::isSymmetric() const
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

}   /* TridiagonalMatrix< T >::isSymmetric() */

template< typename T >
bool TridiagonalMatrix< T >::isTridiagonal() const
{
    return( true );

}   /* TridiagonalMatrix< T >::isTridiagonal() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator+( const DenseMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* TridiagonalMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator+( const LowerTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* TridiagonalMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator+( const UpperTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* TridiagonalMatrix< T >::operator+() */

template< typename T >
TridiagonalMatrix< T > TridiagonalMatrix< T >::operator+( const DiagonalMatrix< T > & rhs ) const
{
    TridiagonalMatrix< T > result( *this );
    return( result += rhs );

}   /* TridiagonalMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator+( const SymmetricMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* TridiagonalMatrix< T >::operator+() */


template< typename T >
TridiagonalMatrix< T > TridiagonalMatrix< T >::operator+( const TridiagonalMatrix< T > & rhs ) const
{
    TridiagonalMatrix< T > result( *this );
    return( result += rhs );

}   /* TridiagonalMatrix< T >::operator+() */

template< typename T >
TridiagonalMatrix< T > & TridiagonalMatrix< T >::operator+=( const DiagonalMatrix< T > & rhs )
{
    // Error check
    if( ( this->n != rhs.getRows() ) || ( this->n != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be added together" );
    }

    // Add corresponding elements on diagonal
    for( int i = 0; i < this->n; ++i )
    {
        this->set( i, i, this->get( i, i ) + rhs.get( i, i ) );
    }

    return( *this );

}   /* TridiagonalMatrix< T >::operator+=() */

template< typename T >
TridiagonalMatrix< T > & TridiagonalMatrix< T >::operator+=( const TridiagonalMatrix< T > & rhs )
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

}   /* TridiagonalMatrix< T >::operator+=() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator-( const DenseMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* TridiagonalMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator-( const LowerTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* TridiagonalMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator-( const UpperTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* TridiagonalMatrix< T >::operator-() */

template< typename T >
TridiagonalMatrix< T > TridiagonalMatrix< T >::operator-( const DiagonalMatrix< T > & rhs ) const
{
    TridiagonalMatrix< T > result( *this );
    return( result -= rhs );

}   /* TridiagonalMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator-( const SymmetricMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* TridiagonalMatrix< T >::operator-() */

template< typename T >
TridiagonalMatrix< T > TridiagonalMatrix< T >::operator-( const TridiagonalMatrix< T > & rhs ) const
{
    TridiagonalMatrix< T > result( *this );
    return( result -= rhs );

}   /* TridiagonalMatrix< T >::operator-() */

template< typename T >
TridiagonalMatrix< T > & TridiagonalMatrix< T >::operator-=( const DiagonalMatrix< T > & rhs )
{
    // Error check
    if( ( this->n != rhs.getRows() ) || ( this->n != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be subtracted together" );
    }

    // Subtract corresponding elements on diagonal
    for( int i = 0; i < this->n; ++i )
    {
        this->set( i, i, this->get( i, i ) - rhs.get( i, i ) );
    }

    return( *this );

}   /* TridiagonalMatrix< T >::operator-=() */

template< typename T >
TridiagonalMatrix< T > & TridiagonalMatrix< T >::operator-=( const TridiagonalMatrix< T > & rhs )
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

}   /* TridiagonalMatrix< T >::operator-=() */

template< typename T >
TridiagonalMatrix< T > TridiagonalMatrix< T >::operator-() const
{
    TridiagonalMatrix< T > result( *this );
    
    // Negate rows
    for( int i = 0; i < result.getRows(); ++i )
    {
        result.data[ i ] = -( result.data[ i ] );
    }
    return( result );

}   /* TridiagonalMatrix< T >::operator-=() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator*( const DenseMatrix< T > & rhs ) const
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
            for( int k = 0; k < this->getRows(); ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );

}   /* TridiagonalMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator*( const LowerTriangularMatrix< T > & rhs ) const
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
            for( int k = j; k < rhs.getRows(); ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );

}   /* TridiagonalMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator*( const UpperTriangularMatrix< T > & rhs ) const
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

}   /* TridiagonalMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator*( const DiagonalMatrix< T > & rhs ) const
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
            result.set( i, j, this->get( i, j ) * rhs.get( j, j ) );
        }
    }
    return( result );

}   /* TridiagonalMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator*( const SymmetricMatrix< T > & rhs ) const
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
            for( int k = 0; k < this->getRows(); ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );

}   /* TridiagonalMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > TridiagonalMatrix< T >::operator*( const TridiagonalMatrix< T > & rhs ) const
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
            for( int k = 0; k < this->getRows(); ++k )
            {
                sum += ( this->get( i, k ) * rhs.get( k, j ) );
            }
            result.set( i, j, sum );
        }
    }
    return( result );

}   /* TridiagonalMatrix< T >::operator*() */

template< typename T >
TridiagonalMatrix< T > TridiagonalMatrix< T >::operator*( const double scalar ) const
{
    TridiagonalMatrix< T > result( *this );
    return( result *= scalar );

}   /* TridiagonalMatrix< T >::operator*() */

template< typename T >
TridiagonalMatrix< T > & TridiagonalMatrix< T >::operator*=( const double scalar )
{
    for( int i = 0; i < this->n; ++i )
    {
        this->data[ i ] *= scalar;
    }
    return( *this );

}   /* TridiagonalMatrix< T >::operator*=() */

template< typename T >
TridiagonalMatrix< T > TridiagonalMatrix< T >::transpose() const
{
    TridiagonalMatrix< T > result( *this );
    for( int i = 1; i < result.getRows(); ++i )
    {
        result.set( i, i - 1, this->get( i - 1, i ) );
    }
    for( int i = 0; i < result.getRows() - 1; ++i )
    {
        result.set( i, i + 1, this->get( i + 1, i ) );
    }
    return( result );

}   /* TridiagonalMatrix< T >::transpose() */


#endif  /* TRIDIAGONAL_MATRIX_HPP */