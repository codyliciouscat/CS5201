/**
 * @file solver.hpp
 * @author Caleb Berg (cabr29@mst.edu)
 * @author Cody Moore (cmmyv8@mst.edu)
 * @brief Implementation header file for Solver class
 * @date 2019-05-12
 */

#ifndef SOLVER_HPP
#define SOLVER_HPP

#include "solver.h"
#include <iostream>
#include <chrono>
using namespace std;

template< typename T >
Vector< T > Solver< T >::gaussianElimination( const AbstractMatrix< T > & a, const Vector< T > & b ) const
{
    // Error check
    if( a.getColumns() != b.size() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrix and vector are not compatible sizes" );
    }
    int size = b.size();
    Vector< T > row_max(size);
    Vector< double > solution(size);
    Vector< int > perm(size);
    DenseMatrix< T > a_copy( size, size );
    Vector< T > b_copy( b );

    // copy a to a_copy
    for( int y = 0; y < size; ++y )
    {
        for( int x = 0; x < size; ++x )
        {
            a_copy.set( y, x, a.get( y, x ));
        }
    }

    // fill row_max and perm vectors
    for( int i = 0; i < size; ++i )
    {
        perm[ i ] = i;
        row_max[ i ] = getMaxValueInRow(a_copy, i);
    }

    for( int k = 0; k < size - 1; ++k )
    {
        double pivot_element = 0, // current pivot element
               max_pivot = 0;     // max pivot element
        int max_pos = 0;          // position of max pivot element

        // compute pivot elements and find max pivot element
        for( int i = k; i < size; ++i )
        {
            pivot_element = fabs( a_copy.get( perm[ i ], k ) / row_max[ i ] );
            if( pivot_element > max_pivot )
            {
                max_pivot = pivot_element;
                max_pos = i;
            }
        }

        // swap perm[ max_pos ] and perm[ k ]
        int tmp = perm[ max_pos ];
        perm[ max_pos ] = perm[ k ];
        perm[ k ] = tmp;

        // forward elimination
        for( int i = k + 1; i < size; ++i )
        {
            double x_mult = a_copy.get( perm[ i ], k ) / a_copy.get( perm[ k ], k );
            for( int j = k; j < size; ++j )
            {
                a_copy.set( perm[i], j, a_copy.get( perm[i], j ) - ( x_mult * a_copy.get( perm[ k ], j ) ) );
            }
            b_copy[ perm[ i ] ] = b_copy[ perm[ i ] ] - ( x_mult * b_copy[ perm[ k ] ] );
        }
    }

    // back substitution
    // calculates last solution (x) value
    solution[ size - 1 ] = b_copy[ perm[ size - 1 ] ] / a_copy.get( perm[ size - 1 ], size - 1 );
    // calculates solution (x) value for i = 0 : size - 2
    for( int i = size - 2; i >= 0; --i )
    {
        double sum = b_copy[ perm[ i ] ];
        for( int j = i + 1; j < size; ++j )
        {
            sum = sum - ( a_copy.get( perm[ i ], j ) ) * solution[ j ];
        }
        solution[ i ] = sum / a_copy.get( perm[ i ], i );
    }

    return solution;

}   /* Solver< T >::gaussian() */

template< typename T >
Vector< T > Solver< T >::FiniteDifferenceGauss( const AbstractMatrix< T > & a, const Vector< T > & b ) const
{
    // Error check
    if( a.getColumns() != b.size() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrix and vector are not compatible sizes" );
    }
    int size = b.size();
    Vector< double > solution(size);
    DenseMatrix< T > a_copy( size, size );
    Vector< T > b_copy( b );
    int meshTightness = static_cast<int>(sqrt(size) + 1);

    // copy a to a_copy
    for( int y = 0; y < size; ++y )
    {
        for( int x = 0; x < size; ++x )
        {
            a_copy.set( y, x, a.get( y, x ));
        }
    }

    for( int k = 0; k < size - 1; ++k )
    {

        // compute pivot elements and find max pivot element
        int i_end = k + meshTightness;
        if( i_end > size )
          i_end = size;

        // forward elimination
        for( int i = k + 1; i < i_end; ++i )
        {   
            double x_mult = a_copy.get( i , k ) / a_copy.get( k , k );
            // calculates end of nonzero elements
            int j_end = i + meshTightness;
            if( j_end > size )
              j_end = size;
            for( int j = k; j < j_end; ++j )
            {
                a_copy.set( i, j, a_copy.get( i, j ) - ( x_mult * a_copy.get( k , j ) ) );
            }
            b_copy[ i ] = b_copy[ i ] - ( x_mult * b_copy[ k ] );
        }
    }

    // back substitution
    // calculates last solution (x) value
    solution[ size - 1 ] = b_copy[ size - 1  ] / a_copy.get( size - 1 , size - 1 );
    // calculates solution (x) value for i = 0 : size - 2
    for( int i = size - 2; i >= 0; --i )
    {
        double sum = b_copy[ i  ];
        int j_end = i + meshTightness;
        if( j_end > size )
          j_end = size;
        for( int j = i + 1; j < j_end; ++j )
        {
            sum = sum - ( a_copy.get( i , j ) ) * solution[ j ];
        }
        solution[ i ] = sum / a_copy.get( i , i );
    }

    return solution;
}   /* Solver< T >::FiniteDifferenceGaussian() */

template< typename T >
Vector< T > Solver< T >::forwardSubstitution( const AbstractMatrix< T > & a, const Vector< T > & b ) const
{
    // Error check
    if( a.getColumns() != b.size() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrix and vector are not compatible sizes" );
    }
    Vector< T > x( b.size() );

    for( int i = 0; i < x.size(); ++i )
    {
        x[ i ] = b[ i ];
        for( int xi = i - 1; xi >= 0 ; --xi )
        {
            x[ i ] -= a.get( i, xi ) * x[ xi ];
        }
        x[ i ] /= a.get( i, i );
    }

    return( x );

}   /* Solver< T >::forwardSubstitution() */

template< typename T >
Vector< T > Solver< T >::backwardSubstitution( const AbstractMatrix< T > & a, const Vector< T > & b ) const
{
    // Error check
    if( a.getColumns() != b.size() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrix and vector are not compatible sizes" );
    }
    Vector< T > x( b.size() );

    for( int i = x.size() - 1; i >= 0; --i )
    {
        x[ i ] = b[ i ];
        for( int xi = i + 1; xi < x.size(); ++xi )
        {
            x[ i ] -= a.get( i, xi ) * x[ xi ];
        }
        x[ i ] /= a.get( i, i );
    }

    return( x );

}   /* Solver< T >::backwardSubstitution() */

template< typename T >
Vector< T > Solver< T >::thomas( const AbstractMatrix< T > & a, const Vector< T > & b ) const
{
    // Error check
    if( a.getColumns() != b.size() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrix and vector are not compatible sizes" );
    }
    int n = b.size();
    Vector< T > x( n ); /* Solution */
    Vector< T > aDiagonal( n ); /* Lower diagonal   */
    Vector< T > bDiagonal( n ); /* Middle diagonal  */
    Vector< T > cDiagonal( n ); /* Upper diagonal   */
    Vector< T > d( b );

    /* Fill diagonal vectors */
    for( int i = 0; i < n; ++i )
    {
        if( i > 0 )
        {
            aDiagonal[ i ] = a.get( i, i - 1 );
        }
        bDiagonal[ i ] = a.get( i, i );
        if( i < ( n - 1 ) )
        {
            cDiagonal[ i ] = a.get( i, i + 1 );
        }
    }
    
    for( int k = 1; k < n; ++k )
    {
        T m = ( aDiagonal[ k ] / bDiagonal[ k - 1 ] );
        bDiagonal[ k ] = ( bDiagonal[ k ] - ( m * cDiagonal[ k - 1 ] ) );
        d[ k ] = ( d[ k ] - ( m * d[ k - 1 ] ) );
    }

    x[ n - 1 ] = ( d[ n - 1 ] / bDiagonal[ n - 1 ] );
    for( int k = n - 2; k >= 0; --k )
    {
        x[ k ] = ( ( d[ k ] - ( cDiagonal[ k ] * x[ k + 1 ] ) ) / bDiagonal[ k ] );
    }

    return( x );

}   /* Solver< T >::thomas() */

template< typename T >
Vector< T > Solver< T >::choleskyDecomposition( const AbstractMatrix< T > & a, const Vector< T > & b ) const
{
    // Error check
    if( a.getColumns() != b.size() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrix and vector are not compatible sizes" );
    }
    int n = b.size();
    SymmetricMatrix< T > l( n );    /* Stores L and L transpose automatically */

    for( int i = 0; i < n; i++ )
    {
        for( int j = 0; j <= i; j++ )
        {
            if( i == j )
            {
                /* Diagonal formula */
                T sum = 0;
                for( int k = 0; k < j; ++k )
                {
                    sum += ( l.get( j, k ) * l.get( j, k ) );
                }
                l.set( j, j, sqrt( a.get( j, j ) - sum ) );
            }
            else
            {
                /* Non-diagonal formula */
                T sum = 0;
                for( int k = 0; k < j; ++k )
                {
                    sum += ( l.get( i, k ) * l.get( j, k ) );
                }
                l.set( i, j, ( ( a.get( i, j ) - sum ) / l.get( j, j ) ) );
            }
        }
    }

    Vector< T > y = this->forwardSubstitution( l, b );
    Vector< T > x = this->backwardSubstitution( l, y );

    return( x );

}   /* Solver< T >::choleskyDecomposition() */

template< typename T >
Vector< T > Solver< T >::FiniteDifferenceCholesky( const AbstractMatrix< T > & a, const Vector< T > & b ) const
{
    // Error check
    if( a.getColumns() != b.size() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrix and vector are not compatible sizes" );
    }
    int n = b.size();
    SymmetricMatrix< T > l( n );    /* Stores L and L transpose automatically */
    int meshTightness = static_cast<int>(sqrt(n) + 1);

    for( int i = 0; i < n; i++ )
    {
        int begin = 0;
        if( i >= meshTightness )
          begin = i - meshTightness + 1;
        for( int j = begin; j <= i; j++ )
        {
            if( i == j )
            {
                /* Diagonal formula */
                T sum = 0;
                for( int k = begin; k < j; ++k )
                {
                    sum += ( l.get( j, k ) * l.get( j, k ) );
                }
                l.set( j, j, sqrt( a.get( i, j ) - sum ) );
            }
            else
            {
                /* Non-diagonal formula */
                T sum = 0;
                for( int k = begin; k < j; ++k )
                {
                    sum += ( l.get( i, k ) * l.get( j, k ) );
                }
                l.set( i, j, ( ( a.get( i, j ) - sum ) / l.get( j, j ) ) );
            }
        }
    }

    // forward substitution
    Vector< T > y( n );
    y[ 0 ] = b[ 0 ] / l.get(0, 0);
    for( int i = 1; i < n; i++ )
    {
        double sum = b[ i ];
        int begin = 0;
        if( i >= meshTightness )
          begin = i - meshTightness + 1;
        // calculates solultion (x) value for j = begin : i
        for( int j = begin; j < i; ++j )
        {
            sum = sum - ( l.get( i, j ) * y[ j ] );
        }
        y[ i ] = sum / l.get( i, i );
    }

    // backward substitution
    Vector< T > x( n );
    x[ n - 1 ] = y[ n - 1  ] / l.get( n - 1, n - 1 );
    // calculates solution (x) value for i = 0 : size - 2
    for( int i = n - 2; i >= 0; --i )
    {
        double sum = y[ i  ];
        int j_end = i + meshTightness;
        if( j_end > n )
          j_end = n;
        for( int j = i + 1; j < j_end; ++j )
        {
            sum = sum - ( l.get( i, j ) ) * x[ j ];
        }
        x[ i ] = sum / l.get( i, i );
    }

    return( x );

}   /* Solver< T >::FiniteDifferenceCholesky() */

template< typename T >
T Solver< T >::getMaxValueInRow( const AbstractMatrix< T > & a, const int row ) const
{
    T maxValue = 0;
    for( int j = 0; j < a.getColumns(); ++j )
    {
        if( abs( a.get( row, j ) ) > maxValue )
        {
            maxValue = abs( a.get( row, j ) );
        }
    }
    return( maxValue );

}   /* Solver< T >::getMaxValueInRow() */

template< typename T >
int Solver< T >::getNextPivotRow( const AbstractMatrix< T > & a, const int iteration ) const
{
    int nextRow = iteration;
    for( int i = iteration + 1; i < a.getRows(); ++i )
    {
        if( abs( static_cast< double >( a.get( i, iteration ) ) / this->getMaxValueInRow( a, i ) )
            > abs( static_cast< double >( a.get( nextRow, iteration ) ) / this->getMaxValueInRow( a, nextRow ) ) )
        {
            nextRow = i;
        }
    }
    return( nextRow );

}   /* Solver< T >::getNextPivotRow() */

template< typename T >
Vector< T > Solver< T >::operator()( const AbstractMatrix< T > & a, const Vector< T > & b ) const
{
    // Error check
    if( a.getColumns() != b.size() )
    {
        throw Error< SIZE_MISMATCH_ERROR >( "Matrix and vector are not compatible sizes" );
    }

    Vector< T > x( b.size() );
    if( a.isDiagonal() )
    {
        for( int i = 0; i < x.size(); ++i )
        {
            x[ i ] = ( b[ i ] / a.get( i, i ) );
        }
    }
    else if( a.isTridiagonal() )
    {
        x = this->thomas( a, b );
    }
    else if( a.isSymmetric() )
    {
        x = this->choleskyDecomposition( a, b );
    }
    else if( a.isLowerTriangular() )
    {
        x = this->forwardSubstitution( a, b );
    }
    else if( a.isUpperTriangular() )
    {
        x = this->backwardSubstitution( a, b );
    }
    else
    {
        x = this->gaussianElimination( a, b );
    }
    return( x );

}   /* Solver< T >::operator()() */

#include "solver.hpp"

#endif  /* SOLVER_H */