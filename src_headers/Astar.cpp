#include "Astar.hpp"
#include <float.h>
using namespace IOApp;

Astar::Astar(Points *pointCollection)
{
    this->_grid=*pointCollection;

    _numDirections=_grid._numberProperties-1;
    _Np.resize(_numDirections,0);
    _tortuosity.resize(_numDirections,0.0);
    _connectingDistance=1;
    _optimalpaths.resize(_grid.path_couples.size()/2);
    ComputePathsDimension();

}

//*************************************************************************************
Astar::~Astar()
{
    _Np.clear();

 }

//*************************************************************************************
Output::ExitCodes Astar::ComputePathsDimension()
{
    for(unsigned int i =0;i<_numDirections;i++)
        _Np.at(i)=_grid._startSizes.at(i)*_grid._targetSizes.at(i);

    return Output::Success;
}
//*************************************************************************************

double Astar::ComputePathLenght(unsigned int &pathIdx, bool serial)
{
    double lenght=0.0;
    unsigned int size= _optimalpaths.at(pathIdx).size();

    for(unsigned int i =1; i<size;i++)
    {
        lenght+=ComputeDistance(_optimalpaths.at(pathIdx).at(i-1),_optimalpaths.at(pathIdx).at(i));
    }

    return lenght;
}

//*************************************************************************************

Output::ExitCodes Astar::ComputeTortuosity(bool asmean)
{
    double minTort=FLT_MAX;
    double globminTort;
    vector<int> dims(3,0);
    dims.at(0)=_grid._nx;
    dims.at(1)=_grid._ny;
    dims.at(2)=_grid._nz;
    unsigned int offSet=0, numel=0;
    for(unsigned int i =0;i<_numDirections;i++)
    {
        numel=0;
        if(asmean)
        {
            minTort=0.0;
        }
        else
        {
            minTort=FLT_MAX;
        }
#if USE_OPENMP
#pragma omp parallel for num_threads(THREAD_NUM)
#endif
        for(unsigned int j=0; j< _Np.at(i);j++)
        {

            unsigned int idx=j+offSet;
            double ln=ComputePathLenght(idx);

            if(ln>1)
            {
                if(asmean)
                {
                    minTort=minTort+ln;
                    numel++;
                }
                else if(minTort>ln)
                {
                    minTort=ln;
                }
            }
        }
        if(asmean)
        {
            minTort=minTort/numel;///Mean Value
        }
        _tortuosity.at(i)=minTort/(dims.at(i)-1);
        offSet+=_Np.at(i);
    }
    return Output::Success;
}

//*************************************************************************************

Output::ExitCodes Astar::ComputeNeighborhood()
{
    unsigned int numPo=_grid._numberOfPores;
    int dim= (_grid._numberProperties-1)*(_grid._numberProperties-1)*(_grid._numberProperties-1)*_connectingDistance*_connectingDistance;


#if USE_OPENMP
#pragma omp parallel for num_threads(THREAD_NUM)
#endif
    for(unsigned int i=0; i< numPo;i++)
    {
        SinglePoint* current=&_grid._grid.at(i);
        if(current->medium==false)
        {
            current->_neighbourhood.reserve(dim);
            GetNeighsCD(*current, current->_neighbourhood);
        }


    }
    ///destroy grid
    _grid._gridFull.clear();
    vector<SinglePoint>().swap(_grid._gridFull);

#if USE_OPENMP
    printf("Thread %i stopped computing neighboors\n", omp_get_thread_num());
    cout<<flush;
#endif
    return Output::Success;
}

//*************************************************************************************
Output::ExitCodes Astar::AStar3D()
{

    unsigned int np=0,numpaths=0,deadpaths=0;
    unsigned int numCouple=_grid.path_couples.size();
    bool reverse=false;


    string outputFold=OUTPUTFOLDER;


#if USE_OPENMP
#pragma omp parallel for num_threads(THREAD_NUM)
#endif
    for(unsigned int couple=0; couple<numCouple;couple+=2)
    {

        unsigned int currentNumberPoints=0,numPoro=_grid._numberOfPores, numberPointsReversing;

#if USE_OPENMP
        if((couple%10==0&&numCouple<1000)||(couple%100==0 && numCouple>=1000))
        {
            float missing= ((numCouple/THREAD_NUM)/2 -couple);
            printf("Thread %i in loop couple %d\n", omp_get_thread_num(),couple);
        }
        cout<<flush;
#endif


        string direct="Couple"+to_string(couple/2);
#if VERBOSE
        cout << _grid.path_couples.at(couple+1) << "\t "<< _grid.path_couples.at(couple+1)<< endl;
        cout<<flush;
#endif
        SinglePoint target=_grid._grid.at(_grid.path_couples.at(couple+1));

        SinglePoint inlet=_grid._grid.at(_grid.path_couples.at(couple));


reverse:
        vector<unsigned int> open_close;
        open_close.resize(numPoro,0); /// 1 assigned to open - 2 assigned to close
        ///costs vector
        vector<float> hcosts;
        vector<float> fcosts;
        vector<float> gcosts;
        vector<unsigned int> parentVector;
        hcosts.resize(numPoro,0.0);
        gcosts.resize(numPoro,0.0);
        fcosts.resize(numPoro,FLT_MAX);
        parentVector.resize(numPoro,0);
        currentNumberPoints=0;

        SinglePoint current=inlet;
        open_close.at(current.pourous_index)=1;
        double dist=ComputeDistance(inlet,target);
        fcosts.at(current.pourous_index)=dist;
        hcosts.at(current.pourous_index)=dist;

        float cost;
        currentNumberPoints=0;
        bool cond=true;
        if(inlet._deadSuspect==true || target._deadSuspect==true)
        {
            cond=false;
            deadpaths++;
            np++;
            numpaths++;
        }
rewhile:
        while(cond)
        {

            cost=FLT_MAX;
            int id;
            MinOpenFScore(current,cost,id,open_close,hcosts,fcosts);
            if(cost>=1.0e+12 )
            {

                _optimalpaths.at(np).reserve(1);
                PrintDeadPath(currentNumberPoints, (int)(couple/2));
                if(currentNumberPoints<=(int)(numPoro/2))
                {
                    inlet._deadSuspect=true;
                    _grid._grid.at(inlet.pourous_index)._deadSuspect=true;

                }
                else
                {
                    numberPointsReversing=currentNumberPoints;
                    SinglePoint cmd;//=(&target);
                    equal(cmd,target);
                    equal(target,inlet);
                    equal(inlet,cmd);
                    open_close.at(target.pourous_index)=1;
                    open_close.at(inlet.pourous_index)=1;

                    goto reverse;
                }

                deadpaths++;
                np++;
                numpaths++;


                cond=false;
                goto rewhile;
            }
            currentNumberPoints++;
            if(current==target)
            {
                _optimalpaths.at(np).reserve(_grid._nz);
                ReversePath(_optimalpaths.at(couple/2),current, inlet, parentVector);
#if VERBOSE
                PrintPath(_optimalpaths.at(couple/2), false, -1);

#endif

                cout<< flush;

                np++;
                numpaths++;
                cond=false;
                goto rewhile;
            }


            open_close.at(id)=2;
            fcosts.at(id)=FLT_MAX;
#if VERBOSE
#if USE_OPENMP
            if(currentNumberPoints%1000==0)
            {
             printf("Thread %i in loop explored %d points \n", omp_get_thread_num(),currentNumberPoints);
            }
#endif
#endif
            unsigned int nsize=current._neighbourhood.size();
            for(unsigned int i=0;i<nsize;i++)
            {
                SinglePoint& point=*current._neighbourhood.at(i);

                double newMovementCost=gcosts.at(current.pourous_index)+ComputeDistance(current,point);

                if(point.medium==1 || (gcosts.at(point.pourous_index)<=newMovementCost && open_close.at(point.pourous_index)!=0 ))//(open_close.at(point.pourous_index)==2 && gcosts.at(point.pourous_index)<=newMovementCost))
                {
                    continue;
                }

               if(newMovementCost< gcosts.at(point.pourous_index) || open_close.at(point.pourous_index)!=2)
                {


                    gcosts.at(point.pourous_index)=newMovementCost;
                    hcosts.at(point.pourous_index)=ComputeDistance(point,target);
                    fcosts.at(point.pourous_index)=hcosts.at(point.pourous_index)+gcosts.at(point.pourous_index);
                    parentVector.at(point.pourous_index)=current.pourous_index;
                    if(open_close.at(point.pourous_index)!=1)
                        open_close.at(point.pourous_index)=1;

                }

            }///Neighs loop

        }///While end
        cout << flush;
        double lenght;
        unsigned int idx=couple/2;
        lenght = ComputePathLenght(idx,true);
        unsigned int size= _optimalpaths.at(idx).size();



    }///Paths couple

    ComputeTortuosity(true);


#if USE_OPENMP
#pragma omp parallel for num_threads(THREAD_NUM)
#endif
    for(unsigned int couple=0; couple<numCouple;couple+=2)///Print Paths On File
    {
#if VERBOSE
        PrintPath(_optimalpaths.at(couple/2), false, (int)(couple/2));
#endif
    }


    return Output::Success;
}

//*************************************************************************************

Output::ExitCodes Astar::GetNeighsCD(SinglePoint &A, vector<SinglePoint *> &neighs)
{
    neighs.reserve(27);
    int start= _connectingDistance* (-1);
    int end= (int)(_connectingDistance);
    for(int x=start; x<=end; x++)
    {
        for(int y=start;y<=end;y++)
        {
            for(int z=start;z<=end;z++)
            {
                if(x==0&& y==0 && z==0)
                    continue;

                int checkx=A.x+x;
                int checky=A.y+y;
                int checkz=A.z+z;
                if(checkx>=0 && checkx < _grid._nx && checky>=0 && checky < _grid._ny && checkz>=0 && checkz < _grid._nz)
                {
                    unsigned int vicino=checkx+checky*_grid._nx+checkz*_grid._nx*_grid._ny;///_grid.GetPointIndex((unsigned int&)(checkx),(unsigned int&)checky,(unsigned int&)checkz);
                    int id=_grid._gridFull.at(vicino).pourous_index;
                    if(id!=-1)
                        neighs.push_back(&_grid._grid.at(id));


                }
            }
        }
    }
    return Output::Success;
}
//*************************************************************************************

Output::ExitCodes Astar::MinOpenFScore(SinglePoint &current, float &cost, int& id, vector<unsigned int>& _open_close, vector<float>& hcosts, vector<float> &fcosts)
{



    auto ite= min_element(fcosts.begin(),fcosts.end());
    int index= distance(fcosts.begin(), ite);

    if(ite!= fcosts.end())
    {
        cost=fcosts.at(index);
        equal(current,_grid._grid.at(index));
        id= index;

    }
    else
        return Output::Abort;


    return Output::Success;
}


//*************************************************************************************


Output::ExitCodes Astar::ReversePath(vector<SinglePoint> &path, SinglePoint &current, SinglePoint &start, vector<unsigned int> &parents)
{
    SinglePoint point= current;
    path.push_back(point);
    while(!SamePoint(point,start))
    {
        int id=point.pourous_index;
        point=_grid._grid.at(parents.at(id));
        path.push_back(point);
    }
    return Output::Success;
}

//*************************************************************************************

Output::ExitCodes Astar::PrintDeadPath(int numpo,int onfile)
{
    SinglePoint& target=_grid._grid.at(_grid.path_couples.at(onfile+1));
    SinglePoint& inlet=_grid._grid.at(_grid.path_couples.at(onfile));
    if(onfile>=0)
    {
#if VERBOSE
        string outputFold=OUTPUTFOLDER;
        ofstream outFile;
        outFile.open(outputFold+FILE_SEPARATOR+"Couple_"+to_string(onfile)+".mat",ios_base::out);


        outFile << " Couple " <<onfile<<" from "<< inlet.x << " , " <<  inlet.y << " , "<< inlet.z << " to " << target.x << " , "<< target.y << " , "<< target.z<< endl;
        outFile<< " The path is dead, Astar explored "<< numpo<< " points.";

        outFile.close();
#endif
    }else
    {
#if VERBOSE
        cout << " Couple " <<onfile<<" from "<< inlet.x << " , " <<  inlet.y << " , "<< inlet.z << " to " << target.x << " , "<< target.y << " , "<< target.z<< endl;
        cout<< " The path is dead, Astar explored "<< numpo<< " points.";
#endif
    }

    return Output::Success;
}


//*************************************************************************************

Output::ExitCodes Astar::PrintPath(vector<SinglePoint> &path,  bool complete, int onfile)
{
    unsigned int size=path.size();

    if(onfile>=0)
    {
        string outputFold=OUTPUTFOLDER;
        ofstream outFile;
        outFile.open(outputFold+FILE_SEPARATOR+"Couple_"+to_string(onfile)+".mat",ios_base::out);
#if ONFILE
        int a =1;
        for(unsigned int i =0 ; i< size; i++)
        {
            outFile<<path.at(i).y +a << " "<< path.at(i).x +a <<" "<< path.at(i).z+a<<endl;
        }
#else
        SinglePoint& target=_grid._grid.at(_grid.path_couples.at(onfile+1));
        SinglePoint& inlet=_grid._grid.at(_grid.path_couples.at(onfile));
        outFile << " Couple " <<onfile<<" from "<< inlet.x << " , " <<  inlet.y << " , "<< inlet.z << " to " << target.x << " , "<< target.y << " , "<< target.z<< endl << "Size: "<<size<< endl;

#endif
        outFile.close();
    }



    return Output::Success;
}
//*************************************************************************************


Output::ExitCodes Properties(Points &pointCollection,unsigned int &numInlets,string& testName)
{
    pointCollection.ComputeNumberGrains();
    pointCollection.ComputePorosity();

    string inputFolder = INPUTFOLDER;
    unsigned int numdata=numInlets;


    string prefix=inputFolder+FILE_SEPARATOR+testName+FILE_SEPARATOR;
    string inletFile=prefix+"test_inlet_x.txt";
    pointCollection.ConfigureInletOutlet(inletFile,1,true,numdata);
    inletFile=prefix+"test_outlet_x.txt";
    pointCollection.ConfigureInletOutlet(inletFile,1,false,numdata);
    inletFile=prefix+"test_inlet_y.txt";
    pointCollection.ConfigureInletOutlet(inletFile,2,true,numdata);
    inletFile=prefix+"test_outlet_y.txt";
    pointCollection.ConfigureInletOutlet(inletFile,2,false,numdata);
    inletFile=prefix+"test_inlet_z.txt";
    pointCollection.ConfigureInletOutlet(inletFile,3,true,numdata);
    inletFile=prefix+"test_outlet_z.txt";
    pointCollection.ConfigureInletOutlet(inletFile,3,false,numdata);
    pointCollection.ConfigureInletOutletVector();

    return Output::Success;
}
