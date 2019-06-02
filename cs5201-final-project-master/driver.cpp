/**
 * @file driver.cpp
 * @author Caleb Berg (cabr29@mst.edu)
 * @author Cody Moore (cmmyv8@mst.edu)
 * @brief Program main file
 * @date 2019-05-12
 */

/**
 * @brief Get M_PI from cmath
 * 
 */
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <chrono>
#include <iomanip>
#include "abstractmatrix.h"
#include "densematrix.h"
#include "lowertriangularmatrix.h"
#include "uppertriangularmatrix.h"
#include "diagonalmatrix.h"
#include "symmetricmatrix.h"
#include "tridiagonalmatrix.h"
#include "vector.h"
#include "error.h"
#include "solver.h"

using namespace std;

/**
 * @brief Suppress unused parameter warning
 * 
 */
#define UNUSED( var ) (void)( var )

/**
 * @brief Maps a value from one range to another
 * 
 * @param value 
 * @param fromLow 
 * @param fromHigh 
 * @param toLow 
 * @param toHigh
 * 
 * @pre fromLow <= fromHigh; toLow <= toHigh 
 * @return Returns value proportional to new range
 */
double map( const double value, const double fromLow, const double fromHigh, const double toLow, const double toHigh )
{
    return( toLow + ( ( toHigh - toLow ) / ( fromHigh - fromLow ) ) * ( value - fromLow ) );

}   /* map() */

/**
 * @brief Forcing function for Poisson equation
 * 
 * @param x 
 * @param y 
 * @pre None
 * @return Returns forcing function evaluated at x and y
 */
double f( const double x, const double y )
{
    UNUSED( x );
    UNUSED( y );
    return( 0.0 );

}   /* f() */

/**
 * @brief Exact solution to problem
 * 
 * @param x 
 * @param y 
 * @pre None
 * @return Returns exact solution for x and y
 */
double u( const double x, const double y )
{
    return( ( 1.0 / sinh( M_PI ) ) * ( ( sin( x ) * sinh( M_PI - y ) ) + ( sin( y ) * sinh( M_PI - x ) ) ) );

}   /* u() */

/**
 * @brief Solves Dirichlet problem
 * 
 * @tparam forcingFunc( double, double ) 
 * @tparam bottomBoundFunc( double ) 
 * @tparam topBoundFunc( double ) 
 * @tparam leftBoundFunc( double ) 
 * @tparam rightBoundFunc( double ) 
 * @param meshTightness 
 * @param bottomBound 
 * @param topBound 
 * @param leftBound 
 * @param rightBound 
 * 
 * @pre bottomBound <= topBound; leftBound <= rightBound; meshTightness > 1
 * 
 * @post Generates and solves system of equations for Fininte Difference Method
 * 
 * @return Solution to system of equations for Finite Difference Method
 * 
 * @throw Throws Error< PARAMETER_ERROR > when meshTightness is less than 1
 */
template< double forcingFunc( double, double ), double bottomBoundFunc( double ), double topBoundFunc( double ), double leftBoundFunc( double ), double rightBoundFunc( double ) >
Vector< double > finiteDifferenceMethod( const int meshTightness /* N */, const double bottomBound, const double topBound, const double leftBound, const double rightBound )
{
    if( meshTightness < 1 )
    {
        throw Error< PARAMETER_ERROR >( "meshTightness must be greater than one" );
    }
    
    int N = meshTightness;
    double h = ( ( rightBound - leftBound ) / N );
    DenseMatrix< double > A( ( N - 1 ) * ( N - 1 ), ( N - 1 ) * ( N - 1 ) );
    Vector< double > b( ( N - 1 ) * ( N - 1 ) );

    /* Iterate over interior points of mesh */
    int directions[ 4 ][ 2 ] =
    {
        { -1, 0 },
        { 1, 0 },
        { 0, -1 },
        { 0, 1 }
    };

    // generates start time
    chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();
    // generates coefficient matrix
    for( int k = 1; k <= ( N - 1 ); ++k )
    {
        for( int j = 1; j <= ( N - 1 ); ++j )
        {
            int i = ( ( ( N - 1 ) * ( k - 1 ) ) + ( j - 1 ) );
            A.set( i, i, 1.0 );
            for( int d = 0; d < 4; ++d )
            {
                int j_new = j + directions[ d ][ 0 ];
                int k_new = k + directions[ d ][ 1 ];
                if( j_new == 0 )
                {
                    b[ i ] += ( 0.25 * leftBoundFunc( map( k_new, 0, N, bottomBound, topBound ) ) );
                }
                else if( j_new == N )
                {
                    b[ i ] += ( 0.25 * rightBoundFunc( map( k_new, 0, N, bottomBound, topBound ) ) );
                }
                else if( k_new == 0 )
                {
                    b[ i ] += ( 0.25 * bottomBoundFunc( map( j_new, 0, N, leftBound, rightBound ) ) );
                }
                else if( k_new == N )
                {
                    b[ i ] += ( 0.25 * topBoundFunc( map( j_new, 0, N, leftBound, rightBound ) ) );
                }
                else
                {
                    int i_new = ( ( ( N - 1 ) * ( k_new - 1 ) ) + ( j_new - 1 ) );
                    A.set( i, i_new, -0.25 );
                }
            }
            b[ i ] -= ( ( -( h * h ) / 4.0 ) * forcingFunc( map( j, 0, N, leftBound, rightBound ), map( k, 0, N, bottomBound, topBound ) ) );
        }
    }
    // generates stop time
    chrono::high_resolution_clock::time_point stopTime = chrono::high_resolution_clock::now();
    // calculates coefficient matrix generation time
    chrono::duration< double > timeSpan = chrono::duration_cast< chrono::duration< double > >( stopTime - startTime );
    double seconds = timeSpan.count();

    cout << endl << "N = " << N << endl;
    cout << "Matrix Generation Time: " << seconds << "s" << endl;

    Solver< double > solver;

    // calculate start time
    startTime = chrono::high_resolution_clock::now();
    // find x using Cholesky
    Vector< double > x = solver.FiniteDifferenceCholesky( A, b );
    // calculate stop time
    stopTime = chrono::high_resolution_clock::now();
    // calculate execution time for Cholesky
    timeSpan = chrono::duration_cast< chrono::duration< double > >( stopTime - startTime );
    seconds = timeSpan.count();

    cout << "Cholesky Execution Time: " << seconds << "s" << endl;
    cout << "CHOLESKY: (A * x) == b: " << ( A * x == b ? "True" : "False" ) << endl;

    // calculate start time
    startTime = chrono::high_resolution_clock::now();
    // find x using Gaussian
    x = solver.FiniteDifferenceGauss( A, b );
    // calculate stop time
    stopTime = chrono::high_resolution_clock::now();
    // calculate execution time for Cholesky
    timeSpan = chrono::duration_cast< chrono::duration< double > >( stopTime - startTime );
    seconds = timeSpan.count();

    cout << "Gaussian Execution Time: " << seconds << "s" << endl;
    cout << "GAUSSIAN: (A * x) == b: " << ( A * x == b ? "True" : "False" ) << endl;
 
    return( x );

}   /* finiteDifferenceMethod() */

/**
 * @brief Return zero
 * 
 * @param n Unused
 * @return Zero
 */
double zero( const double n )
{
    UNUSED( n );
    return( 0 );

}   /* zero() */

/**
 * @brief Main program execution
 * 
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main( int argc, char ** argv )
{
    try
    {
        // Check command line usage
        if( argc != 2 )
        {
            throw( Error< COMMAND_LINE_ARGUMENTS_ERROR >( string( "Invalid usage\nUsage: " ) + string( argv[ 0 ] ) + string( " <N>" ) ) );
        }

        int N = stoi( argv[ 1 ] );

        Vector< double > solution = finiteDifferenceMethod< f, sin, zero, sin, zero >( N, 0, M_PI, 0, M_PI );

        // calculate error
        double error = 0.0;
        double sum = 0.0;
        for( int i = 0; i < solution.size(); ++i )
        {
            double x_coord = map( ( i % ( N - 1 ) ), 0, N, 0, M_PI );
            double y_coord = map( ( i / ( N - 1 ) ), 0, N, 0, M_PI );
            double expected = u( x_coord, y_coord );
            double estimate = solution[ i ];
            double tmp = estimate - expected;
            sum += tmp * tmp;
        }
        error = sqrt( ( M_PI / N ) * sum );

        cout << "Error: " << error << endl << endl;
    }
    catch( Error< COMMAND_LINE_ARGUMENTS_ERROR > error )
    {
        cerr << "ERROR( COMMAND_LINE_ARGUMENTS_ERROR ): " << error.message() << endl;
    }
    catch( Error< INDEX_ERROR > error )
    {
        cerr << "ERROR( INDEX_ERROR ): " << error.message() << endl;
    }
    catch( Error< SIZE_MISMATCH_ERROR > error )
    {
        cerr << "ERROR( SIZE_MISMATCH_ERROR ): " << error.message() << endl;
    }
    catch( Error< PARAMETER_ERROR > error )
    {
        cerr << "ERROR( PARAMETER_ERROR ): " << error.message() << endl;
    }
    catch( Error< MATRIX_SIZE_ERROR > error )
    {
        cerr << "ERROR( MATRIX_SIZE_ERROR ): " << error.message() << endl;
    }
    catch( Error< PARSE_ERROR > error )
    {
        cerr << "ERROR( PARSE_ERROR ): " << error.message() << endl;
    }
    catch( ... )
    {
        cerr << "No clue what happened :/" << endl;
    }

    return( 0 );

}   /* main() */
