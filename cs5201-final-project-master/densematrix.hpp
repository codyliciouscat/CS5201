/**
 * @file densematrix.hpp
 * @author Caleb Berg (cabr29@mst.edu)
 * @brief Implementation header file for DenseMatrix class
 * @date 2019-04-03
 */

#ifndef DENSE_MATRIX_HPP
#define DENSE_MATRIX_HPP

#include "densematrix.h"

template< typename T >
DenseMatrix< T >::DenseMatrix() :
    data( nullptr ),
    rows( 0 ),
    cols( 0 )
{

}   /* DenseMatrix< T >::DenseMatrix()*/

template< typename T >
DenseMatrix< T >::DenseMatrix( const int rows, const int cols ) :
    rows( rows ),
    cols( cols )
{
    // Error check
    if( ( ( this->rows < 0 ) || ( this->cols < 0 ) ) ||
        ( ( this->rows > 0 ) && ( this->cols == 0 ) ) ||
        ( ( this->rows == 0 ) && ( this->cols > 0 ) ) )
    {
        throw Error< PARAMETER_ERROR >( "Invalid matrix dimensions" );
    }

    // Construct matrix
    this->data = new Vector< T >[ this->rows ]();
    for( int i = 0; i < this->rows ; ++i )
    {
        this->data[ i ] = Vector< T >( this->cols );
    }

}   /* DenseMatrix< T >::DenseMatrix() */

template< typename T >
DenseMatrix< T >::DenseMatrix( const DenseMatrix< T > & source ) :
    data( nullptr ),
    rows( source.rows ),
    cols( source.cols )
{
    // Check if any deep copying needs to be done (dimensions > 0)
    if( this->rows > 0 )
    {
        // Copy data
        this->data = new Vector< T >[ this->rows ]();
        for( int i = 0; i < this->rows; ++i )
        {
            this->data[ i ] = source.data[ i ];
        }
    }

}   /* DenseMatrix< T >::DenseMatrix() */

template< typename T >
DenseMatrix< T >::DenseMatrix( const LowerTriangularMatrix< T > & source ) :
    data( nullptr ),
    rows( source.getRows() ),
    cols( source.getColumns() )
{
    // Check if any deep copying needs to be done (dimensions > 0)
    if( this->rows > 0 )
    {
        // Copy data
        this->data = new Vector< T >[ this->rows ]();
        for( int i = 0; i < this->rows; ++i )
        {
            this->data[ i ] = Vector< T >( this->cols );
            for( int j = 0; j <= i; ++j )
            {
                this->data[ i ][ j ] = source.get( i, j );
            }
        }
    }

}   /* DenseMatrix< T >::DenseMatrix() */

template< typename T >
DenseMatrix< T >::DenseMatrix( const UpperTriangularMatrix< T > & source ) :
    data( nullptr ),
    rows( source.getRows() ),
    cols( source.getColumns() )
{
    // Check if any deep copying needs to be done (dimensions > 0)
    if( this->rows > 0 )
    {
        // Copy data
        this->data = new Vector< T >[ this->rows ]();
        for( int i = 0; i < this->rows; ++i )
        {
            this->data[ i ] = Vector< T >( this->cols );
            for( int j = i; j < this->cols; ++j )
            {
                this->data[ i ][ j ] = source.get( i, j );
            }
        }
    }

}   /* DenseMatrix< T >::DenseMatrix() */

template< typename T >
DenseMatrix< T >::DenseMatrix( const DiagonalMatrix< T > & source ) :
    data( nullptr ),
    rows( source.getRows() ),
    cols( source.getColumns() )
{
    // Check if any deep copying needs to be done (dimensions > 0)
    if( this->rows > 0 )
    {
        // Copy data
        this->data = new Vector< T >[ this->rows ]();
        for( int i = 0; i < this->rows; ++i )
        {
            this->data[ i ] = Vector< T >( this->cols );
            this->data[ i ][ i ] = source.get( i, i );
        }
    }

}   /* DenseMatrix< T >::DenseMatrix() */

template< typename T >
DenseMatrix< T >::DenseMatrix( const SymmetricMatrix< T > & source ) :
    data( nullptr ),
    rows( source.getRows() ),
    cols( source.getColumns() )
{
    // Check if any deep copying needs to be done (dimensions > 0)
    if( this->rows > 0 )
    {
        // Copy data
        this->data = new Vector< T >[ this->rows ]();
        for( int i = 0; i < this->rows; ++i )
        {
            this->data[ i ] = Vector< T >( this->cols );
            for( int j = 0; j <= i; ++j )
            {
                this->data[ i ][ j ] = source.get( i, j );
                this->data[ j ][ i ] = source.get( i, j );
            }
        }
    }

}   /* DenseMatrix< T >::DenseMatrix() */

template< typename T >
DenseMatrix< T >::DenseMatrix( const TridiagonalMatrix< T > & source ) :
    data( nullptr ),
    rows( source.getRows() ),
    cols( source.getColumns() )
{
    // Check if any deep copying needs to be done (dimensions > 0)
    if( this->rows > 0 )
    {
        // Copy data
        this->data = new Vector< T >[ this->rows ]();
        for( int i = 0; i < this->rows; ++i )
        {
            this->data[ i ] = Vector< T >( this->cols );
            for( int j = 0; j < this->cols; ++j )
            {
                this->data[ i ][ j ] = source.get( i, j );
            }
        }
    }

}   /* DenseMatrix< T >::DenseMatrix() */

template< typename T >
DenseMatrix< T >::DenseMatrix( const AbstractMatrix< T > & source ) :
    data( nullptr ),
    rows( source.getRows() ),
    cols( source.getColumns() )
{
    // Check if any deep copying needs to be done (dimensions > 0)
    if( this->rows > 0 )
    {
        // Copy data
        this->data = new Vector< T >[ this->rows ]();
        for( int i = 0; i < this->rows; ++i )
        {
            this->data[ i ] = Vector< T >( this->cols );
            for( int j = 0; j < this->cols; ++j )
            {
                this->data[ i ][ j ] = source.get( i, j );
            }
        }
    }

}   /* DenseMatrix< T >::DenseMatrix() */

template< typename T >
DenseMatrix< T >::~DenseMatrix()
{
    // Check if destruction is necessary
    if( this->data != nullptr )
    {
        delete[] this->data;
        this->rows = 0;
        this->cols = 0;
    }

}   /* DenseMatrix< T >::~DenseMatrix() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator=( const DenseMatrix< T > & source )
{
    // Copy Swap Idiom
    DenseMatrix< T > sourceCopy( source );
    swap( this->data, sourceCopy.data );
    swap( this->rows, sourceCopy.rows );
    swap( this->cols, sourceCopy.cols );
    return( *this );

}   /* DenseMatrix< T >::operator=() */

template< typename T >
T DenseMatrix< T >::get( const int row, const int col ) const
{
    // Error check
    if( ( row < 0 ) || ( row >= this->rows ) || ( col < 0 ) || ( col >= this->cols ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix indices" );
    }
    
    return( this->data[ row ][ col ] );

}   /* DenseMatrix< T >::get() */

template< typename T >
void DenseMatrix< T >::set( const int row, const int col, const T & val )
{
    // Error check
    if( ( row < 0 ) || ( row >= this->rows ) || ( col < 0 ) || ( col >= this->cols ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix indices" );
    }

    this->data[ row ][ col ] = val;

}   /* DenseMatrix< T >::set() */

template< typename T >
Vector< T > DenseMatrix< T >::operator*( const Vector< T > & rhs ) const
{
    // Error check
    if( this->cols != rhs.size() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrix is unable to be multiplied by vector" );
    }

    Vector< T > result( this->rows );

    // Multiply matrix and vector
    for( int i = 0; i < result.size(); ++i )
    {
        result[ i ] = ( this->data[ i ] * rhs );
    }

    return( result );

}   /* DenseMatrix< T >::operator*() */

template< typename T >
int DenseMatrix< T >::getColumns() const
{
    return( this->cols );

}   /* DenseMatrix< T >::getColumns() */

template< typename T >
int DenseMatrix< T >::getRows() const
{
    return( this->rows );

}   /* DenseMatrix< T >::getRows() */

template< typename T >
bool DenseMatrix< T >::isLowerTriangular() const
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

}   /* DenseMatrix< T >::isLowerTriangular() */

template< typename T >
bool DenseMatrix< T >::isUpperTriangular() const
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

}   /* DenseMatrix< T >::isUpperTriangular() */

template< typename T >
bool DenseMatrix< T >::isDiagonal() const
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

}   /* DenseMatrix< T >::isDiagonal() */

template< typename T >
bool DenseMatrix< T >::isSymmetric() const
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

}   /* DenseMatrix< T >::isSymmetric() */

template< typename T >
bool DenseMatrix< T >::isTridiagonal() const
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

}   /* DenseMatrix< T >::isTridiagonal() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator+( const DenseMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* DenseMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator+( const LowerTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* DenseMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator+( const UpperTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* DenseMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator+( const DiagonalMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* DenseMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator+( const SymmetricMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* DenseMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator+( const TridiagonalMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result += rhs );

}   /* DenseMatrix< T >::operator+() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator+=( const DenseMatrix< T > & rhs )
{
    // Error check
    if( ( this->rows != rhs.getRows() ) || ( this->cols != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be added together" );
    }

    // Add corresponding rows
    for( int i = 0; i < this->rows; ++i )
    {
        this->data[ i ] += rhs.data[ i ];
    }

    return( *this );

}   /* DenseMatrix< T >::operator+=() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator+=( const LowerTriangularMatrix< T > & rhs )
{
    // Error check
    if( ( this->rows != rhs.getRows() ) || ( this->cols != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be added together" );
    }

    // Add corresponding rows
    for( int i = 0; i < this->rows; ++i )
    {
        for( int j = 0; j <= i; ++j )
        {
            this->data[ i ][ j ] += rhs.get( i, j );
        }
    }

    return( *this );

}   /* DenseMatrix< T >::operator+=() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator+=( const UpperTriangularMatrix< T > & rhs )
{
    // Error check
    if( ( this->rows != rhs.getRows() ) || ( this->cols != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be added together" );
    }

    // Add corresponding rows
    for( int i = 0; i < this->rows; ++i )
    {
        for( int j = i; j < this->cols; ++j )
        {
            this->data[ i ][ j ] += rhs.get( i, j );
        }
    }

    return( *this );

}   /* DenseMatrix< T >::operator+=() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator+=( const DiagonalMatrix< T > & rhs )
{
    // Error check
    if( ( this->rows != rhs.getRows() ) || ( this->cols != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be added together" );
    }

    // Add corresponding rows
    for( int i = 0; i < this->rows; ++i )
    {
        this->data[ i ][ i ] += rhs.get( i, i );
    }

    return( *this );

}   /* DenseMatrix< T >::operator+=() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator+=( const SymmetricMatrix< T > & rhs )
{
    // Error check
    if( ( this->rows != rhs.getRows() ) || ( this->cols != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be added together" );
    }

    // Add corresponding rows
    for( int i = 0; i < this->rows; ++i )
    {
        for( int j = 0; j < this->cols; ++j )
        {
            this->data[ i ][ j ] += rhs.get( i, j );
        }
    }

    return( *this );

}   /* DenseMatrix< T >::operator+=() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator+=( const TridiagonalMatrix< T > & rhs )
{
    // Error check
    if( ( this->rows != rhs.getRows() ) || ( this->cols != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be added together" );
    }

    // Add corresponding rows
    for( int i = 0; i < this->rows; ++i )
    {
        for( int j = 0; j < this->cols; ++j )
        {
            this->data[ i ][ j ] += rhs.get( i, j );
        }
    }

    return( *this );

}   /* DenseMatrix< T >::operator+=() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator-( const DenseMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* DenseMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator-( const LowerTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* DenseMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator-( const UpperTriangularMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* DenseMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator-( const DiagonalMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* DenseMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator-( const SymmetricMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* DenseMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator-( const TridiagonalMatrix< T > & rhs ) const
{
    DenseMatrix< T > result( *this );
    return( result -= rhs );

}   /* DenseMatrix< T >::operator-() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator-=( const DenseMatrix< T > & rhs )
{
    // Error check
    if( ( this->rows != rhs.getRows() ) || ( this->cols != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be subtracted together" );
    }

    // Subtract corresponding rows
    for( int i = 0; i < this->rows; ++i )
    {
        this->data[ i ] -= rhs.data[ i ];
    }

    return( *this );

}   /* DenseMatrix< T >::operator-=() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator-=( const LowerTriangularMatrix< T > & rhs )
{
    // Error check
    if( ( this->rows != rhs.getRows() ) || ( this->cols != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be subtracted together" );
    }

    // Subtract corresponding elements
    for( int i = 0; i < this->rows; ++i )
    {
        for( int j = 0; j <= i; ++j )
        {
            this->data[ i ][ j ] -= rhs.get( i, j );
        }
    }

    return( *this );

}   /* DenseMatrix< T >::operator-=() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator-=( const UpperTriangularMatrix< T > & rhs )
{
    // Error check
    if( ( this->rows != rhs.getRows() ) || ( this->cols != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be subtracted together" );
    }

    // Subtract corresponding elements
    for( int i = 0; i < this->rows; ++i )
    {
        for( int j = i; j < this->cols; ++j )
        {
            this->data[ i ][ j ] -= rhs.get( i, j );
        }
    }

    return( *this );

}   /* DenseMatrix< T >::operator-=() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator-=( const DiagonalMatrix< T > & rhs )
{
    // Error check
    if( ( this->rows != rhs.getRows() ) || ( this->cols != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be subtracted together" );
    }

    // Subtract corresponding elements
    for( int i = 0; i < this->rows; ++i )
    {
        this->data[ i ][ i ] -= rhs.get( i, i );
    }

    return( *this );

}   /* DenseMatrix< T >::operator-=() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator-=( const SymmetricMatrix< T > & rhs )
{
    // Error check
    if( ( this->rows != rhs.getRows() ) || ( this->cols != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be subtracted together" );
    }

    // Subtract corresponding elements
    for( int i = 0; i < this->rows; ++i )
    {
        for( int j = 0; j < this->cols; ++j )
        {
            this->data[ i ][ j ] -= rhs.get( i, j );
        }
    }

    return( *this );

}   /* DenseMatrix< T >::operator-=() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator-=( const TridiagonalMatrix< T > & rhs )
{
    // Error check
    if( ( this->rows != rhs.getRows() ) || ( this->cols != rhs.getColumns() ) )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be subtracted together" );
    }

    // Subtract corresponding elements
    for( int i = 0; i < this->rows; ++i )
    {
        for( int j = 0; j < this->cols; ++j )
        {
            this->data[ i ][ j ] -= rhs.get( i, j );
        }
    }

    return( *this );

}   /* DenseMatrix< T >::operator-=() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator-() const
{
    DenseMatrix< T > result( *this );
    
    // Negate rows
    for( int i = 0; i < result.getRows(); ++i )
    {
        result.data[ i ] = -( result.data[ i ] );
    }
    return( result );

}   /* DenseMatrix< T >::operator-=() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator*( const DenseMatrix< T > & rhs ) const
{
    // Error check
    if( this->cols != rhs.getRows() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be multiplied together" );
    }

    DenseMatrix< T > result( this->rows, rhs.getColumns() );

    // Multiply matrices
    for( int i = 0; i < result.rows; ++i )
    {
        for( int j = 0; j < result.cols; ++j )
        {
            T sum = 0;
            for( int k = 0; k < rhs.getRows(); ++k )
            {
                sum += ( this->data[ i ][ k ] * rhs.get( k, j ) );
            }
            result.data[ i ][ j ] = sum;
        }
    }
    return( result );

}   /* DenseMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator*( const LowerTriangularMatrix< T > & rhs ) const
{
    // Error check
    if( this->cols != rhs.getRows() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be multiplied together" );
    }

    DenseMatrix< T > result( this->rows, rhs.getColumns() );

    // Multiply matrices
    for( int i = 0; i < result.rows; ++i )
    {
        for( int j = 0; j < result.cols; ++j )
        {
            T sum = 0;
            for( int k = j; k < rhs.getRows(); ++k )
            {
                sum += ( this->data[ i ][ k ] * rhs.get( k, j ) );
            }
            result.data[ i ][ j ] = sum;
        }
    }
    return( result );

}   /* DenseMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator*( const UpperTriangularMatrix< T > & rhs ) const
{
    // Error check
    if( this->cols != rhs.getRows() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be multiplied together" );
    }

    DenseMatrix< T > result( this->rows, rhs.getColumns() );

    // Multiply matrices
    for( int i = 0; i < result.rows; ++i )
    {
        for( int j = 0; j < result.cols; ++j )
        {
            T sum = 0;
            for( int k = 0; k <= j; ++k )
            {
                sum += ( this->data[ i ][ k ] * rhs.get( k, j ) );
            }
            result.data[ i ][ j ] = sum;
        }
    }
    return( result );

}   /* DenseMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator*( const DiagonalMatrix< T > & rhs ) const
{
    // Error check
    if( this->cols != rhs.getRows() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be multiplied together" );
    }

    DenseMatrix< T > result( this->rows, rhs.getColumns() );

    // Multiply matrices
    for( int i = 0; i < result.rows; ++i )
    {
        for( int j = 0; j < result.cols; ++j )
        {
            result.set( i, j, this->get( i, j ) * rhs.get( j, j ) );
        }
    }
    return( result );

}   /* DenseMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator*( const SymmetricMatrix< T > & rhs ) const
{
    // Error check
    if( this->cols != rhs.getRows() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be multiplied together" );
    }

    DenseMatrix< T > result( this->rows, rhs.getColumns() );

    // Multiply matrices
    for( int i = 0; i < result.rows; ++i )
    {
        for( int j = 0; j < result.cols; ++j )
        {
            T sum = 0;
            for( int k = 0; k < rhs.getRows(); ++k )
            {
                sum += ( this->data[ i ][ k ] * rhs.get( k, j ) );
            }
            result.data[ i ][ j ] = sum;
        }
    }
    return( result );

}   /* DenseMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator*( const TridiagonalMatrix< T > & rhs ) const
{
    // Error check
    if( this->cols != rhs.getRows() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are unable to be multiplied together" );
    }

    DenseMatrix< T > result( this->rows, rhs.getColumns() );

    // Multiply matrices
    for( int i = 0; i < result.rows; ++i )
    {
        for( int j = 0; j < result.cols; ++j )
        {
            T sum = 0;
            for( int k = 0; k < rhs.getRows(); ++k )
            {
                sum += ( this->data[ i ][ k ] * rhs.get( k, j ) );
            }
            result.data[ i ][ j ] = sum;
        }
    }
    return( result );

}   /* DenseMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::operator*( const double scalar ) const
{
    DenseMatrix< T > result( *this );
    return( result *= scalar );

}   /* DenseMatrix< T >::operator*() */

template< typename T >
DenseMatrix< T > & DenseMatrix< T >::operator*=( const double scalar )
{
    for( int i = 0; i < this->rows; ++i )
    {
        this->data[ i ] *= scalar;
    }
    return( *this );

}   /* DenseMatrix< T >::operator*=() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::transpose() const
{
    DenseMatrix< T > result( this->cols, this->rows );
    for( int i = 0; i < result.rows; ++i )
    {
        for( int j = 0; j < result.cols; ++j )
        {
            result.data[ i ][ j ] = this->data[ j ][ i ];
        }
    }
    return( result );

}   /* DenseMatrix< T >::transpose() */

template< typename T >
DenseMatrix< T > DenseMatrix< T >::augment( const DenseMatrix< T > & denseMatrix ) const
{
    // Error check
    if( this->rows != denseMatrix.rows )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrices are incompatible sizes for augmentation" );
    }

    DenseMatrix< T > result( this->rows, this->cols + denseMatrix.cols );
    for( int i = 0; i < result.rows; ++i )
    {
        for( int j = 0; j < this->cols; ++j )
        {
            result.data[ i ][ j ] = this->data[ i ][ j ];
        }
        for( int j = 0; j < denseMatrix.cols; ++j )
        {
            result.data[ i ][ this->cols + j ] = denseMatrix.data[ i ][ j ];
        }
    }
    return( result );

}   /* DenseMatrix< T >::augment() */

template< typename T >
Vector< T > & DenseMatrix< T >::operator[]( const int index )
{
    // Error check
    if( ( index < 0 ) || ( index >= this->rows ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix index" );
    }
    return( this->data[ index ] );

}   /* DenseMatrix< T >::operator[]() */

template< typename T >
const Vector< T > & DenseMatrix< T >::operator[]( const int index ) const
{
    // Error check
    if( ( index < 0 ) || ( index >= this->rows ) )
    {
        throw Error< INDEX_ERROR >( "Invalid matrix index" );
    }
    return( this->data[ index ] );

}   /* DenseMatrix< T >::operator[]() */

template< typename T >
DenseMatrix< T >::DenseMatrix( const Vector< T > & source ) :
    data( new Vector< T >[ source.size() ]() ),
    rows( source.size() ),
    cols( 1 )
{
    // Copy elements
    for( int i = 0; i < this->rows ; ++i )
    {
        data[ i ] = Vector< T >( 1 );
    }
    for( int i = 0; i < this->rows; ++i )
    {
        data[ i ][ 0 ] = source[ i ];
    }

}   /* DenseMatrix< T >::DenseMatrix() */

#endif  /* DENSE_MATRIX_HPP */