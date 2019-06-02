/**
 * @file vector.hpp
 * @author Caleb Berg (cabr29@mst.edu)
 * @brief Implementation of templated vector class
 * @date 2019-02-26
 */

#ifndef VECTOR_HPP
#define VECTOR_HPP

#include "vector.h"

template< typename T >
void Vector< T >::clear()
{
    // Check if vector isn't already cleared
    if( this->data != nullptr )
    {
        delete[] this->data;
        this->data = nullptr;
        this->numElements = 0;
    }

}   /* Vector< T >::clear() */

template< typename T >
Vector< T >::Vector()  : data( nullptr ), numElements( 0 )
{

}   /* Vector< T >::Vector() */

template< typename T >
Vector< T >::Vector( const int size )
{
    // Error check
    if( size < 0 )
    {
        throw( Error< PARAMETER_ERROR >( "Invalid size parameter" ) );
    }
    this->data = new T[ size ];
    this->numElements = size;

    // Initialize all values to 0 for the type
    for( int i = 0; i < this->numElements; ++i )
    {
        this->data[ i ] = 0;
    }

}   /* Vector< T >::Vector() */

template< typename T >
Vector< T >::Vector( const Vector< T > & source )
{
    this->numElements = source.numElements;
    this->data = new T[ this->numElements ];
    for( int i = 0; i < source.numElements; ++i )
    {
        this->data[ i ] = source.data[ i ];
    }

}   /* Vector< T >::Vector() */

template< typename T >
Vector< T >::Vector( Vector< T > && source ) : data( source.data ), numElements( source.numElements )
{
    source.data = nullptr;
    source.numElements = 0;

}   /* Vector< T >::Vector() */

template< typename T >
Vector< T >::~Vector()
{
    this->clear();

}   /* Vector< T >::~Vector() */

template< typename T >
Vector< T > & Vector< T >::operator=( const Vector< T > & source )
{
    // Check if self-assign
    if( this != &source )
    {
        // Copy and swap
        Vector< T > temp( source );
        swap( *this, temp );
    }
    return( *this );

}   /* Vector< T >::operator=() */

template< typename T >
Vector< T > & Vector< T >::operator=( Vector< T > && source )
{
    // Check if self-assign
    if( this != &source )
    {
        delete[] this->data;
        this->data = source.data;
        this->numElements = source.numElements;
        source.data = nullptr;
        source.numElements = 0;
    }
    return( *this );

}   /* Vector< T >::operator=() */

template< typename T >
T & Vector< T >::operator[]( const int index )
{
    // Error check
    if( ( index < 0 ) || ( index >= this->numElements ) )
    {
        throw( Error< INDEX_ERROR >( "Invalid element index" ) );
    }
    return( this->data[ index ] );

}   /* Vector< T >::operator[]() */

template< typename T >
const T & Vector< T >::operator[]( const int index ) const
{
    // Error check
    if( ( index < 0 ) || ( index >= this->numElements ) )
    {
        throw( Error< INDEX_ERROR >( "Invalid element index" ) );
    }
    return( *( this->data + index ) );

}   /* Vector< T >::operator[]() */

template< typename T >
bool Vector< T >::operator==( const Vector< T > & vector ) const
{
    const double EPSILON = 0.05;
    bool isEqual = true;
    if( this->size() != vector.size() )
    {
        isEqual = false;
    }
    int index = 0;
    while( isEqual && ( index < vector.size() ) )
    {
        if( abs( this->data[ index ] - vector.data[ index ] ) > EPSILON )
        {
            isEqual = false;
        }
        ++index;
    }
    return( isEqual );

}   /* Vector< T >::operator==() */

template< typename T >
bool Vector< T >::operator!=( const Vector< T > & vector ) const
{
    return( !( ( *this ) == vector ) );

}   /* Vector< T >::operator==() */

template< typename T >
Vector< T > Vector< T >::operator+( const Vector< T > & vector ) const
{
    Vector< T > result( *this );
    return( result += vector );

}   /* Vector< T >::operator+() */

template< typename T >
Vector< T > & Vector< T >::operator+=( const Vector< T > & vector )
{
    // Error check
    if( this->numElements != vector.numElements )
    {
        throw( Error< SIZE_MISMATCH_ERROR >( "Vectors are different sizes" ) );
    }

    // Sum corresponding elements
    for( int i = 0; i < vector.size(); ++i )
    {
        this->data[ i ] += vector.data[ i ];
    }
    return( *this );

}   /* Vector< T >::operator+=() */

template< typename T >
Vector< T > Vector< T >::operator-( const Vector< T > & vector ) const
{
    Vector< T > result( *this );
    return( result -= vector );

}   /* Vector< T >::operator-() */

template< typename T >
Vector< T > & Vector< T >::operator-=( const Vector< T > & vector )
{
    // Error check
    if( this->numElements != vector.numElements )
    {
        throw( Error< SIZE_MISMATCH_ERROR >( "Vectors are different sizes" ) );
    }

    // Get difference of corresponding elements
    for( int i = 0; i < vector.size(); ++i )
    {
        this->data[ i ] -= vector.data[ i ];
    }
    return( *this );

}   /* Vector< T >::operator-=() */

template< typename T >
Vector< T > Vector< T >::operator-() const
{
    // Negate all elements
    Vector< T > result( *this );
    for( int i = 0; i < this->numElements; ++i )
    {
        result.data[ i ] = -result.data[ i ];
    }
    return( result );

}   /* Vector< T >::operator-() */

template< typename T >
T Vector< T >::operator*( const Vector< T > & vector ) const
{
    T result = 0;

    // Error check
    if( this->numElements != vector.numElements )
    {
        throw( Error< SIZE_MISMATCH_ERROR >( "Vectors are different sizes" ) );
    }

    // Compute dot product (sum of products)
    for( int i = 0; i < vector.size(); ++i )
    {
        result += ( this->data[ i ] * vector.data[ i ] );
    }
    return( result );

}   /* Vector< T >::operator*() */

template< typename T >
Vector< T > Vector< T >::operator*( const T & scalar ) const
{
    Vector< T > result( *this );
    return( result *= scalar );

}   /* Vector< T >::operator*() */

template< typename T >
Vector< T > & Vector< T >::operator*=( const T & scalar )
{
    for( int i = 0; i < this->numElements; ++i )
    {
        this->data[ i ] *= scalar;
    }
    return( *this );

}   /* Vector< T >::operator*() */

template< typename T >
Vector< T > Vector< T >::operator/( const T & scalar ) const
{
    Vector< T > result( *this );
    return( result /= scalar );

}   /* Vector< T >::operator/() */

template< typename T >
Vector< T > & Vector< T >::operator/=( const T & scalar )
{
    for( int i = 0; i < this->numElements; ++i )
    {
        this->data[ i ] /= scalar;
    }
    return( *this );

}   /* Vector< T >::operator/=() */

template< typename T >
int Vector< T >::size() const
{
    return( this->numElements );

}   /* Vector< T >::size() */

/**
 * @brief Swap the data of two vectors
 * 
 * @param vector1 First vector
 * @param vector2 Second vector
 * 
 * @pre Type U must have operator=() defined
 * @post {vector1} and {vector2} have their data swapped
 */
template< typename U >
void swap( Vector< U > & vector1, Vector< U > & vector2 )
{
    std::swap( vector1.numElements, vector2.numElements );
    std::swap( vector1.data, vector2.data );

}   /* swap() */

template< typename T >
std::ostream & operator<<( std::ostream & os, const Vector< T > & vector )
{
    for( int i = 0; i < vector.size(); ++i )
    {
        os << ( ( i == 0 ) ? "" : " " ) << vector[ i ];
    }
    return( os );

}   /* operator<<() */

template< typename T >
std::istream & operator>>( std::istream & is, Vector< T > & vector )
{
    for( int i = 0; i < vector.size(); ++i )
    {
        is >> vector[ i ];
        if( is.fail() )
        {
            throw( Error< PARSE_ERROR >( "Invalid input" ) );
        }
    }
    return( is );

}   /* operator>>() */

#endif