#ifndef POINTS_HPP
#define POINTS_HPP

#include "ImportExport.hpp"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <bits/stdc++.h>
#include <math.h>

using namespace IOApp;

struct SinglePoint{
    int x;
    int y;
    int z;
    bool medium; /// 0 to pores, 1 to medium
    unsigned int globidx; ///Global index of the point
    int pourous_index; /// Index of the point in the pore structure -1 values for medium points
    vector<SinglePoint*> _neighbourhood; ///Vector of pointers to the neighbors
    bool _deadSuspect = false; ///True when the point is suspected to be isolated
    bool operator == (const SinglePoint& A) const {
        return ((A.x == x)&&(A.y==y)&&(A.z==z));
    }
};

/// Compute Distance Between two points
double ComputeDistance(SinglePoint &A, SinglePoint &B);
bool find(SinglePoint &point, vector<SinglePoint *> &list);
void equal(SinglePoint &A, SinglePoint& B);
bool SamePoint(SinglePoint A, SinglePoint *B);
bool SamePoint(SinglePoint A, SinglePoint B);
/// Print function
void PrintPoint(SinglePoint A, bool matlabFormat=false);


class Points
{
protected:
    int _numberProperties;
    vector<SinglePoint> _gridFull; ///Points collection, a vector long as the number of points, each element of the vector has the properties x , y , z , value
    vector<SinglePoint> _grid; ///Points collection, a vector long as the number of pore points, each element of the vector has the properties x , y , z , value
    unsigned int _numberOfPoints;
    unsigned int _numberOfGrains; ///Number of grains in the collection
    unsigned int _numberOfPores; ///Number of Pores
    double _porosity;
    vector <vector<unsigned int>> inlet_indices;
    vector <vector<unsigned int>> target_indices;
    vector <unsigned int> path_couples; ///couples of indices inlet/target 
    vector<unsigned int> _startSizes; 
    vector<unsigned int> _targetSizes;
    vector<unsigned int> _globalIndexing;///It keeps in memory the global index from porous points
    unsigned int _nx, _ny,_nz;

public:
    Points();
    ~Points();

    friend class Astar;///Astar can access to members of Points class
    ///Getters
    ///Get number of grains in collection
    inline unsigned int GetNumberOfGrains(){return _numberOfGrains;}
    ///Get number of points in collection
    inline unsigned int GetNumberOfPoints(){return _numberOfPoints;}
    /// Get the index point by coordinates
    unsigned int GetPointIndex( unsigned int &x, unsigned int &y, unsigned int &z);
    /// Load File for import Points
    ///Get number of points in collection
    inline double GetPorosity(){return _porosity;}
    ///Passing Points Collection
    inline const vector<SinglePoint>& GetGrid(){return _gridFull;}




    ///Number of medium grains
    Output::ExitCodes ComputeNumberGrains();
    ///Porosity of the examined sample
    Output::ExitCodes ComputePorosity();
    ///Get global idx from porous one
    int FromPorousToGlob(int &porous);
      ///Get porous idx from global one
    int FromGlobToPor(unsigned int &globIdx);
    ///Create the Neighborhood of a Point
    Output::ExitCodes GetNeighs(SinglePoint& A, vector<SinglePoint*> &neighs);


    /// 4 properties for 3D simulation
    /// 3 properties for 2D (not debugged yet)
    void SetNumberProperties(int &numProp){_numberProperties=numProp;}

    Output::ExitCodes ImportFromFile(const string &filename);
    /// Load Inlet Outlet points from file, establishing them inlet=true if thei are inlet, false if they are outlet
    Output::ExitCodes ConfigureInletOutlet(const string &filename, const unsigned int &dir, bool inlet, unsigned int numData=100);
    /// Counting the couples avoiding direction
    Output::ExitCodes ConfigureInletOutletVector();
    /// Print the point at index idPoint
    Output::ExitCodes PrintPoint(unsigned int &idPoint);

};

#endif // POINTS_HPP
