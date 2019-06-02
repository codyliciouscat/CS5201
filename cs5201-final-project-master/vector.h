/**
 * @file vector.h
 * @author Caleb Berg (cabr29@mst.edu)
 * @brief Declaration of templated vector class
 * @date 2019-02-26
 */

#ifndef VECTOR_H
#define VECTOR_H

#include <iostream>
#include <cmath>
#include "error.h"

/**
 * @brief A vector (collection of objects)
 * 
 * @pre T must have =, +=, -=, unary -, binary *, != defined,
 */
template< typename T >
class Vector
{
    private:
        T * data;
        int numElements;

    private:
        /**
         * @brief Destroys the vector and makes it of 0 size
         * 
         * @pre None
         * @post Vector is destroyed and made to be size 0
         */
        void clear();

    public:
        /**
         * @brief Default constructor
         *
         * @pre None
         * @post Vector is constructed with 0 elements
         */
        Vector();

        /**
         * @brief Parameterized constructor
         * 
         * @param size the number of elements the vector should have
         * 
         * @pre Type T must have operator=() defined and {size} must be nonnegative
         * @post Vector is constructed with {size} elements 
         * @throws Error< PARAMETER_ERROR > if {size} is negative
         */
        Vector( const int size );

        /**
         * @brief Copy constructor
         * 
         * @param source the vector to be copied
         * 
         * @pre Type T must have operator=() defined
         * @post Vector is constructed as a copy of {source}
         */
        Vector( const Vector< T > & source );

        /**
         * @brief Move constructor
         * 
         * @param source the vector to be move copied
         * 
         * @pre None
         * @post Vector is constructed as a moved copy of {source}
         */
        Vector( Vector< T > && source );

        /**
         * @brief Destructor
         * 
         * @pre None
         * @post Vector is destructed and memory is deallocated
         */
        ~Vector();

        /**
         * @brief Assignment operator
         * 
         * @param source the vector to be copied
         * 
         * @pre Type T must have operator=() defined
         * @post Vector is made to be a copy of {source}
         * @return Reference to calling vector
         */
        Vector< T > & operator=( const Vector< T > & source );

        /**
         * @brief Assignment operator
         * 
         * @param source the vector to be copied
         * 
         * @pre None
         * @post Vector is made to be a moved copy of {source} 
         * @return Reference to calling vector
         */
        Vector< T > & operator=( Vector< T > && source );

        /**
         * @brief Accesses element of vector
         * 
         * @param index of element to access
         * 
         * @pre {index} is within range of [0, size())
         * @return Reference to element at index
         * @throws Error< INDEX_ERROR > if index is out of range
         */
        T & operator[]( const int index );

        /**
         * @brief Accesses element of vector
         * 
         * @param index of element to access
         * 
         * @pre {index} is within range of [0, size())
         * @return Const reference to element at index
         * @throws Error< INDEX_ERROR > if index is out of range
         */
        const T & operator[]( const int index ) const;

        /**
         * @brief Vector comparison
         * 
         * @param vector vector to compare to
         * 
         * @pre None
         * @return true if vectors are the same, else false
         */
        bool operator==( const Vector< T > & vector ) const;

        /**
         * @brief Vector comparison
         * 
         * @param vector vector to compare to
         * 
         * @pre None
         * @return false if vectors are the same, else true
         */
        bool operator!=( const Vector< T > & vector ) const;

        /**
         * @brief Sum corresponding elements of vectors
         * 
         * @param vector Vector to sum with
         * 
         * @pre {vector}.size() must equal this->size()
         * @return Vector with elements as sums of corresponding elements in operands
         * @throws Error< SIZE_MISMATCH_ERROR > if vector sizes differ
         */
        Vector< T > operator+( const Vector< T > & vector ) const;

        /**
         * @brief Sum corresponding elements of vector to calling vector
         * 
         * @param vector Vector to sum with
         * 
         * @pre {vector}.size() must equal this->size()
         *      Type T must have operator+=() defined
         * @post Calling vector's elements are updated as sums
         * @return Reference to calling vector
         * @throws Error< SIZE_MISMATCH_ERROR > if vector sizes differ
         */
        Vector< T > & operator+=( const Vector< T > & vector );

        /**
         * @brief Substract corresponding elements of vectors
         * 
         * @param vector Vector to subtract with
         * 
         * @pre Vectors must be equal in size
         * @return Vector with elements as difference of corresponding elements in operands
         * @throws Error< SIZE_MISMATCH_ERROR > if vector sizes differ
         */
        Vector< T > operator-( const Vector< T > & vector ) const;

        /**
         * @brief Subtract corresponding elements of vector from calling vector
         * 
         * @param vector Vector to subtract by
         * 
         * @pre {vector}.size() must equal this->size()
         *      Type T must have operator-=() defined
         * @post Calling vector's elements are updated as difference
         * @return Reference to calling vector
         * @throws Error< SIZE_MISMATCH_ERROR > if vector sizes differ
         */
        Vector< T > & operator-=( const Vector< T > & vector );

        /**
         * @brief Negate all elements of vector
         * 
         * @pre Type T must have unary operator-() defined
         * @return Vector with elements as negations of calling vector's elements
         */
        Vector< T > operator-() const;

        /**
         * @brief Dot product
         * 
         * @param vector Vector to use as 2nd operand of dot product
         * 
         * @pre {vector}.size() must equal this->size()
         *      Type T must have binary operator*() defined
         * @return Sum of products of corresponding vector elements
         * @throws Error< SIZE_MISMATCH_ERROR > if vector sizes differ
         */
        T operator*( const Vector< T > & vector ) const;

        /**
         * @brief Scalar multiplication
         * 
         * @param scalar Scalar to multiply vector by
         * 
         * @pre None
         * @return Result vector from multiplying vector elements by scalar
         */
        Vector< T > operator*( const T & scalar ) const;

        /**
         * @brief Scalar multiplication
         * 
         * @param scalar Scalar to multiply vector by
         * 
         * @pre Type T must have operator*=() defined
         * @post Calling vector's elements are scaled
         * @return Reference to calling vector
         */
        Vector< T > & operator*=( const T & scalar );

        /**
         * @brief Scalar division
         * 
         * @param scalar Scalar to divide vector by
         * 
         * @pre None
         * @return Result vector from dividing vector elements by scalar
         */
        Vector< T > operator/( const T & scalar ) const;

        /**
         * @brief Scalar division
         * 
         * @param scalar Scalar to divide vector by
         * 
         * @pre Type T must have operator/=() defined
         * @post Calling vector's elements are scaled
         * @return Reference to calling vector
         */
        Vector< T > & operator/=( const T & scalar );

        /**
         * @brief Get number of elements in vector
         * 
         * @pre None
         * @return Size of vector 
         */
        int size() const;

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
        friend void swap( Vector< U > & vector1, Vector< U > & vector2 );
};

/**
 * @brief Stream insertion operator
 * 
 * @param os Stream to insert vector representation into
 * @param vector Vector to insert into stream
 * 
 * @pre Type T must have operator<<() defined
 * @post Stream has representation of vector elements inserted
 * @return Reference to stream 
 */
template< typename T >
std::ostream & operator<<( std::ostream & os, const Vector< T > & vector );

/**
 * @brief Stream extraction operator
 * 
 * @param is Stream to extract vector data from
 * @param vector Vector to fill with data from stream
 * 
 * @pre Vector must have data constructed and numElements set
 * @return Reference to stream
 * @throws Error< PARSE_ERROR > if invalid input is read
 */
template< typename T >
std::istream & operator>>( std::istream & is, Vector< T > & vector );

#include "vector.hpp"

#endif