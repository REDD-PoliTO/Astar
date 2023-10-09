#ifndef ASTAR_HPP
#define ASTAR_HPP
#include <algorithm>
#include "auxMacros.hpp"
#include "Points.hpp"
#include <string>
#if USE_OPENMP
#include <omp.h>
#endif


typedef vector<vector<SinglePoint>> PathVec; ///Vector of the paths, ordered by couple index

using namespace IOApp;
class Astar
{
private:
    unsigned int _connectingDistance;


protected:
    PathVec _optimalpaths;

    vector<unsigned int> _Np; ///Number of Paths
    unsigned int _numDirections; /// Numbero of directions
    vector<double> _tortuosity; ///Tortuositi Values



public:
    Points _grid; ///Points of the pore medium
    Astar(Points* pointCollection);
    ~Astar();

    Output::ExitCodes ComputePathsDimension(); /// number of paths per direction
    double ComputePathLenght(unsigned int& pathIdx, bool serial=false); /// Compute the lenght of a path from the inlet point to the outlet point
    Output::ExitCodes ComputeNeighborhood(); ///Compute neighbours at beginning of computation
    Output::ExitCodes ComputeTortuosity(bool asmean=true); ///If TRUE tortuosity is computed by Mean Value, else as minimum value
    Output::ExitCodes AStar3D(); 
    Output::ExitCodes GetNeighsCD(SinglePoint& A, vector<SinglePoint*> &neighs); ///Ask for the neighborhood of a point
    Output::ExitCodes MinOpenFScore(SinglePoint& current, float &cost, int& id, vector<unsigned int> &_open_close, vector<float>& hcosts, vector<float> &fcosts); ///select the point of minimal fcost
    Output::ExitCodes ReversePath(vector<SinglePoint>& path, SinglePoint& current, SinglePoint &start, vector<unsigned int> &parents); ///Build the path starting from the target
    inline void SetConnectingDistance(const unsigned int& cd){_connectingDistance=cd; return;}
/// Print functions
    Output::ExitCodes PrintDeadPath(int numpo, int onfile=-1); 
    Output::ExitCodes PrintPath(vector<SinglePoint> &path, bool complete=false, int onfile=-1);

};

Output::ExitCodes Properties(Points& pointCollection, unsigned int &numInlets, string &testName);

#endif // ASTAR_HPP
