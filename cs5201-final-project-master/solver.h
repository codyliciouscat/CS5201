/**
 * @file solver.h
 * @author Caleb Berg (cabr29@mst.edu)
 * @author Cody Moore (cmmyv8@mst.edu)
 * @brief Header file for Solver class
 * @date 2019-03-17
 */

#ifndef SOLVER_H
#define SOLVER_H

#include <cmath>
#include <utility>
#include "abstractmatrix.h"
#include "densematrix.h"
#include "lowertriangularmatrix.h"
#include "uppertriangularmatrix.h"
#include "vector.h"

/**
 * @brief Solver for Ax=b with gaussian elimination
 */
template< typename T >
class Solver
{
    public:
        /**
         * @brief Gaussian elimination solver
         * 
         * @param a Matrix that represents linear system
         * @param b Vector that represents linear system
         * 
         * @pre a.getColumns() must equal b.size()
         * @return Solution to system
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        Vector< T > gaussianElimination( const AbstractMatrix< T > & a, const Vector< T > & b ) const;

        /**
         * @brief Optimized Gaussian elimination solver for Finite Difference Method
         * 
         * @param a 
         * @param b 
         * 
         * @pre a.getColumns() must equal b.size()
         * @return Solution to system
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        Vector< T > FiniteDifferenceGauss( const AbstractMatrix< T > & a, const Vector< T > & b) const;

    public:
        /**
         * @brief Forward substitution
         * 
         * @param a Matrix that represents linear system
         * @param b Vector that represents linear system
         * 
         * @pre a.getColumns() must equal b.size()
         * @return Solution to system
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        Vector< T > forwardSubstitution( const AbstractMatrix< T > & a, const Vector< T > & b ) const;

        /**
         * @brief Backward substitution
         * 
         * @param a Matrix that represents linear system
         * @param b Vector that represents linear system
         * 
         * @pre a.getColumns() must equal b.size()
         * @return Solution to system
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        Vector< T > backwardSubstitution( const AbstractMatrix< T > & a, const Vector< T > & b ) const;

        /**
         * @brief Thomas algorithm
         * 
         * @param a Matrix that represents linear system
         * @param b Vector that represents linear system
         * 
         * @pre a.getColumns() must equal b.size()
         * @return Solution to system
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        Vector< T > thomas( const AbstractMatrix< T > & a, const Vector< T > & b ) const;

        /**
         * @brief Cholesky decomposition
         * 
         * @param a Matrix that represents linear system
         * @param b Vector that represents linear system
         * 
         * @pre a.getColumns() must equal b.size()
         * @return Solution to system
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        Vector< T > choleskyDecomposition( const AbstractMatrix< T > & a, const Vector< T > & b ) const;

        /**
         * @brief Optimized Cholesky decomposition for Finite Difference Method
         * 
         * @param a 
         * @param b 
         * 
         * @pre a.getColumns() must equal b.size()
         * @return Solution to system
         * @throws Error< SiZE_MISMATCH_ERROR > if precondition is not met
         */
        Vector< T > FiniteDifferenceCholesky( const AbstractMatrix< T > & a, const Vector< T > & b) const;

        /**
         * @brief Get the maximum absolute value in a row of a matrix
         * 
         * @param a Matrix to get value from
         * @param row Index of row in matrix
         * 
         * @pre row must be in range of [0, a.getRows())
         * @return Maximum value in specified row
         * @throws Error< INDEX_ERROR > if precondition is not met
         */
        T getMaxValueInRow( const AbstractMatrix< T > & a, const int row ) const;

        /**
         * @brief Get the row that will become the next pivot row
         * 
         * @param a Matrix to use
         * @param iteration Which iteration of the algorithm
         * 
         * @pre iteration must be in range [0, a.getRows())
         * @return Index of row
         */
        int getNextPivotRow( const AbstractMatrix< T > & a, const int iteration ) const;

    public:
        /**
         * @brief Pseudo function operator
         * 
         * @param a Matrix that represents linear system
         * @param b Vector that represents linear system
         * 
         * @pre a.getColumns() must equal b.size()
         * @return Solution to system
         * @throws Error< SIZE_MISMATCH_ERROR > if precondition is not met
         */
        Vector< T > operator()( const AbstractMatrix< T > & a, const Vector< T > & b ) const;

};  /* class Solver */

#include "solver.hpp"

#endif  /* SOLVER_H */