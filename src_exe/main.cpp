#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <fstream>
#include <iostream>
#include "Astar.hpp"
//#include "CustomTimer.hpp"
#include "ImportExport.hpp"
#include "Points.hpp"
#include "auxMacros.hpp"
#include <unistd.h>

void process_mem_usage(double& vm_usage, double& resident_set);

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

    //#if USE_OPENMP
//    status = mkdir(nameFold.c_str(), 0777);
    //#else
        status = mkdir(nameFold.c_str());
    //#endif

    Astar astar(&pointCollection);

    astar.ComputeNeighborhood();

    astar.AStar3D();
    cout << endl << "End of main"<< endl;
    return 0;
}

void process_mem_usage(double& vm_usage, double& resident_set)
{
    vm_usage     = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = 2/1024;//sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
}
