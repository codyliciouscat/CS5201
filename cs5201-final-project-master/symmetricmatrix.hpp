/**
 * @file symmetricmatrix.hpp
 * @author Caleb Berg (cabr29@mst.edu)
 * @brief Implementation header file for SymmetricMatrix class
 * @date 2019-04-03
 */

#ifndef SYMMETRIC_MATRIX_HPP
#define SYMMETRIC_MATRIX_HPP

#include "symmetricmatrix.h"

template< typename T >
SymmetricMatrix< T >::SymmetricMatrix() :
    data( nullptr ),
    n( 0 )
{

}	/* SymmetricMatrix< T >::SymmetricMatrix()*/

template< typename T >
SymmetricMatrix< T >::SymmetricMatrix( const int n ) :
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
        this->data[ i ] = Vector< T >( i + 1 );
    }

}	/* SymmetricMatrix< T >::SymmetricMatrix() */

template< typename T >
SymmetricMatrix< T >::SymmetricMatrix( const SymmetricMatrix< T > & source ) :
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

}	/* SymmetricMatrix< T >::SymmetricMatrix() */

template< typename T >
SymmetricMatrix< T >::SymmetricMatrix( const DiagonalMatrix< T > & source ) :
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

}   /* SymmetricMatrix< T >::SymmetricMatrix() */

template< typename T >
SymmetricMatrix< T >::~SymmetricMatrix()
{
    // Check if destruction is necessary
    if( this->data != nullptr )
    {
        delete[] this->data;
        this->n = 0;
    }

}   /* SymmetricMatrix< T >::~SymmetricMatrix() */

template< typename T >
SymmetricMatrix< T > & SymmetricMatrix< T >::operator=( const SymmetricMatrix< T > & source )
{
    // Copy Swap Idiom
    SymmetricMatrix< T > sourceCopy( source );
    swap( this->data, sourceCopy.data );
    swap( this->n, sourceCopy.n );
    return( *this );

}   /* SymmetricMatrix< T >::operator=() */

template< typename T >
T SymmetricMatrix< T >::get( const int row, const int col ) const
{
    // Error check
    if( ( row < 0 ) || ( row >= this->n ) || ( col < 0 ) || ( col >= this->n ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix indices" );
    }
    
    // Check if not diagonal
    return( col > row ? this->data[ col ][ row ] : this->data[ row ][ col ] );

}	/* SymmetricMatrix< T >::get() */

template< typename T >
void SymmetricMatrix< T >::set( const int row, const int col, const T & val )
{
    // Error check
    if( ( row < 0 ) || ( row >= this->n ) || ( col < 0 ) || ( col >= this->n ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix indices" );
    }
    else if( col > row )
    {
        this->data[ col ][ row ] = val;
    }
    else
    {
        this->data[ row ][ col ] = val;
    }

}	/* SymmetricMatrix< T >::set() */

template< typename T >
Vector< T > SymmetricMatrix< T >::operator*( const Vector< T > & rhs ) const
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

}	/* SymmetricMatrix< T >::operator*() */

template< typename T >
int SymmetricMatrix< T >::getColumns() const
{
    return( this->n );

}   /* SymmetricMatrix< T >::getColumns() */

template< typename T >
int SymmetricMatrix< T >::getRows() const
{
    return( this->n );

}   /* SymmetricMatrix< T >::getRows() */

template< typename T >
bool SymmetricMatrix< T >::isLowerTriangular() const
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

}   /* SymmetricMatrix< T >::isLowerTriangular() */

template< typename T >
bool SymmetricMatrix< T >::isUpperTriangular() const
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

}   /* SymmetricMatrix< T >::isUpperTriangular() */

template< typename T >
bool SymmetricMatrix< T >::isDiagonal() const
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

}   /* SymmetricMatrix< T >::isDiagonal() */

template< typename T >
bool SymmetricMatrix< T >::isSymmetric() const
{
    return( true );

}   /* SymmetricMatrix< T >::isSymmetric() */

template< typename T >
bool SymmetricMatrix< T >::isTridiagonal() const
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

}   /* SymmetricMatrix< T >::isTridiagonal() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator+( const DenseMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* SymmetricMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator+( const LowerTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* SymmetricMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator+( const UpperTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* SymmetricMatrix< T >::operator+() */

template< typename T >
SymmetricMatrix< T > SymmetricMatrix< T >::operator+( const DiagonalMatrix< T > & rhs ) const
{
    SymmetricMatrix< T > result( *this );
    return( result += rhs );

}   /* SymmetricMatrix< T >::operator+() */

template< typename T >
SymmetricMatrix< T > SymmetricMatrix< T >::operator+( const SymmetricMatrix< T > & rhs ) const
{
    SymmetricMatrix< T > result( *this );
    return( result += rhs );

}   /* SymmetricMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator+( const TridiagonalMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* SymmetricMatrix< T >::operator+() */

template< typename T >
SymmetricMatrix< T > & SymmetricMatrix< T >::operator+=( const DiagonalMatrix< T > & rhs )
{
    // Error check
    if( ( this->n != rhs.getRows() ) || ( this->n != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be added together" );
    }

    // Add corresponding elements in diagonal
    for( int i = 0; i < this->n; ++i )
    {
        this->set( i, i, this->get( i, i ) + rhs.get( i, i ) );
    }

    return( *this );

}   /* SymmetricMatrix< T >::operator+=() */

template< typename T >
SymmetricMatrix< T > & SymmetricMatrix< T >::operator+=( const SymmetricMatrix< T > & rhs )
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

}   /* SymmetricMatrix< T >::operator+=() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator-( const DenseMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* SymmetricMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator-( const LowerTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* SymmetricMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator-( const UpperTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* SymmetricMatrix< T >::operator-() */

template< typename T >
SymmetricMatrix< T > SymmetricMatrix< T >::operator-( const DiagonalMatrix< T > & rhs ) const
{
    SymmetricMatrix< T > result( *this );
    return( result -= rhs );

}   /* SymmetricMatrix< T >::operator-() */

template< typename T >
SymmetricMatrix< T > SymmetricMatrix< T >::operator-( const SymmetricMatrix< T > & rhs ) const
{
    SymmetricMatrix< T > result( *this );
    return( result -= rhs );

}   /* SymmetricMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator-( const TridiagonalMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* SymmetricMatrix< T >::operator-() */

template< typename T >
SymmetricMatrix< T > & SymmetricMatrix< T >::operator-=( const DiagonalMatrix< T > & rhs )
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

}   /* SymmetricMatrix< T >::operator-=() */

template< typename T >
SymmetricMatrix< T > & SymmetricMatrix< T >::operator-=( const SymmetricMatrix< T > & rhs )
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

}   /* SymmetricMatrix< T >::operator-=() */

template< typename T >
SymmetricMatrix< T > SymmetricMatrix< T >::operator-() const
{
    SymmetricMatrix< T > result( *this );
    
    // Negate rows
    for( int i = 0; i < result.getRows(); ++i )
    {
        result.data[ i ] = -( result.data[ i ] );
    }
    return( result );

}   /* SymmetricMatrix< T >::operator-=() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator*( const DenseMatrix< T > & rhs ) const
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

}   /* SymmetricMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator*( const LowerTriangularMatrix< T > & rhs ) const
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

}   /* SymmetricMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator*( const UpperTriangularMatrix< T > & rhs ) const
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

}   /* SymmetricMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator*( const DiagonalMatrix< T > & rhs ) const
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

}   /* SymmetricMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator*( const SymmetricMatrix< T > & rhs ) const
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

}   /* SymmetricMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > SymmetricMatrix< T >::operator*( const TridiagonalMatrix< T > & rhs ) const
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

}   /* SymmetricMatrix< T >::operator*() */

template< typename T >
SymmetricMatrix< T > SymmetricMatrix< T >::operator*( const double scalar ) const
{
    SymmetricMatrix< T > result( *this );
    return( result *= scalar );

}   /* SymmetricMatrix< T >::operator*() */

template< typename T >
SymmetricMatrix< T > & SymmetricMatrix< T >::operator*=( const double scalar )
{
    for( int i = 0; i < this->n; ++i )
    {
        this->data[ i ] *= scalar;
    }
    return( *this );

}   /* SymmetricMatrix< T >::operator*=() */

template< typename T >
SymmetricMatrix< T > SymmetricMatrix< T >::transpose() const
{
    SymmetricMatrix< T > result( *this );
    return( result );

}   /* SymmetricMatrix< T >::transpose() */


#endif  /* SYMMETRIC_MATRIX_HPP */