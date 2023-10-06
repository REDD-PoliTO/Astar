#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include "Astar.hpp"
#include "ImportExport.hpp"
#include "Points.hpp"
#include "auxMacros.hpp"
#include <unistd.h>

#if USE_MPI ///Parallel native libraries
//#include <mpi.h>
#endif

#if USE_OPENMP
//#include <ompi_config.h>
#include <omp.h>
#endif

using namespace std;
using namespace IOApp;


int main(int argc, char** argv)
{
    int numThreads=1;
    Output::ExitCodes code;
    Points pointCollection;
    double vm, rss;


    int numProperties=4;
    string inputFolder = INPUTFOLDER;
    string testName=argv[1];
    string fileName=inputFolder+FILE_SEPARATOR+testName+FILE_SEPARATOR+"test.txt";
    unsigned int numInletOut=atoi(argv[2]);

#if USE_OPENMP
#pragma omp parallel
    {
        printf("Thread %i from %i\n", omp_get_thread_num(), omp_get_num_threads());
    }

#endif

    pointCollection.SetNumberProperties(numProperties);
    ///Grain value set as 1
    pointCollection.ImportFromFile(fileName);
    /// Define all properties of the collections
    /// dimensions, porosity, inlet-outlet...
    Properties(pointCollection,numInletOut,testName);

    string nameFold=OUTPUTFOLDER;
    int status;

    status = mkdir(nameFold.c_str(), 0777);// For Windows environment (in some cases) remove the code number

    Astar astar(&pointCollection);
    astar.ComputeNeighborhood();
    /// Running the AStar Loop
    astar.AStar3D();
    cout << endl << "End of main"<< endl;
    return 0;
}
