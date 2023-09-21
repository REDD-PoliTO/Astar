#ifndef ___IMPORTEXPORT_HPP
#define ___IMPORTEXPORT_HPP
#include "ctime"
#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdarg.h>
#include <string.h>
#include <vector>
#include "auxMacros.hpp"


using namespace std;

namespace IOApp {

class Output;
class Input;

class Output
{
public:

    /// \brief Exit codes for methods throughout Application
    enum ExitCodes
    {
        Success = 1, ///< Success flag.
        Abort = -1, ///< Abort flag.
        MpiError = -2, ///< MPI error flag.
        PartitionError = -3, ///< Partitioning error flag.
        FileError = -4, ///< File error flag.
        GenericError = -5, ///< Generic error flag.
        UnimplementedMethod = -6 ///< Flag for an unimplemented feature.
    };

};

class Input
{
public:
    //---------------------------- For string delimiter ------------------------
    static vector<string> Split(string s, string delimiter);


};

}
#endif // IMPORTEXPORT_HPP
