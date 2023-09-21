#ifndef ASTAR_HPP
#define ASTAR_HPP
#include <algorithm>
#include "auxMacros.hpp"
#include "Points.hpp"
#include <string>
#if USE_OPENMP
#include <omp.h>
#endif


typedef vector<vector<SinglePoint>> PathVec;

using namespace IOApp;
class Astar
{
private:
    unsigned int _connectingDistance;


protected:
    PathVec _optimalpaths;

    vector<unsigned int> _Np;
    unsigned int _numDirections;
    vector<double> _tortuosity;



public:
    Points _grid;
    Astar(Points* pointCollection);
    ~Astar();

    Output::ExitCodes ComputePathsDimension();
    double ComputePathLenght(unsigned int& pathIdx, bool serial=false);
    Output::ExitCodes ComputeNeighborhood();///Compute neighbours at beginning of computation
    Output::ExitCodes ComputeTortuosity(bool asmean=true); ///If TRUE tortuosity is computed by Mean Value, else as minimum value
    Output::ExitCodes AStar3D();
    Output::ExitCodes GetNeighsCD(SinglePoint& A, vector<SinglePoint*> &neighs);
    Output::ExitCodes MinOpenFScore(SinglePoint& current, float &cost, int& id, vector<unsigned int> &_open_close, vector<float>& hcosts, vector<float> &fcosts);
    Output::ExitCodes ReversePath(vector<SinglePoint>& path, SinglePoint& current, SinglePoint &start, vector<unsigned int> &parents);
    inline void SetConnectingDistance(const unsigned int& cd){_connectingDistance=cd; return;}
    Output::ExitCodes PrintDeadPath(int numpo, int onfile=-1);
    Output::ExitCodes PrintPath(vector<SinglePoint> &path, bool complete=false, int onfile=-1);

};

Output::ExitCodes Properties(Points& pointCollection, unsigned int &numInlets, string &testName);

#endif // ASTAR_HPP
