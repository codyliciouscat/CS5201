/**
 * @file error.h
 * @author Caleb Berg (cabr29@mst.edu)
 * @brief Simple error class
 * @date 2019-02-26
 */

#ifndef ERROR_H
#define ERROR_H

#include <string>

using namespace std;

/**
 * @brief Error enumerations for templated error class
 * 
 */
enum ErrorEnum
{
    INDEX_ERROR,
    PARSE_ERROR,
    MATRIX_SIZE_ERROR,
    COMMAND_LINE_ARGUMENTS_ERROR,
    PARAMETER_ERROR,
    SIZE_MISMATCH_ERROR
};

/**
 * @brief Basic error class
 */
template< ErrorEnum E >
class Error
{
    private:
        string msg;

    public:
        /**
         * @brief Default constructor
         * 
         * @pre None
         * @post Error is constructed with no message
         */
        Error() {}

        /**
         * @brief Parameterized constructor
         * 
         * @param msg Message associated with error
         * 
         * @pre None
         * @post Error is constructed with passed message
         */
        Error( const string & msg ) : msg( msg ) {}

        /**
         * @brief 
         * 
         * @pre None
         * @post Returns the message describing the error as a string
         * @return Message of error
         */
        const string message() const { return( this->msg ); }
};

#endif