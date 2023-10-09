#include "Points.hpp"

using namespace IOApp;
Points::Points()
{
    _numberOfGrains=0;
    _numberProperties=4;
    _nx=0;
    _ny=0;
    _nz=0;
    _targetSizes.resize(_numberProperties-1,0);
    target_indices.resize(_numberProperties-1);
    inlet_indices.resize(_numberProperties-1);
    _startSizes.resize(_numberProperties-1,0);
}
//*************************************************************************************
Points::~Points()
{

    _gridFull.clear();
    _grid.clear();


}


//*************************************************************************************
double ComputeDistance(SinglePoint &A, SinglePoint &B)
{
    double val=0.0;
    val=sqrt((B.x-A.x)*(B.x-A.x)+(B.y-A.y)*(B.y-A.y)+(B.z-A.z)*(B.z-A.z));
    return val;
}

//*************************************************************************************
Output::ExitCodes Points::ComputeNumberGrains()
{
    unsigned int idx=0;
    for(unsigned int i=0; i< _numberOfPoints;i++)
    {

        if(_gridFull.at(i).medium==false)
        {
            _gridFull.at(i).pourous_index=idx;
            _gridFull.at(i).globidx=i;
            idx++;
        }
        else
        {
            _gridFull.at(i).pourous_index=-1;
            _gridFull.at(i).globidx=i;
        }
        _numberOfGrains+=_gridFull.at(i).medium;
    }
    _numberOfPores=_gridFull.size()-_numberOfGrains;
    _globalIndexing.resize(_numberOfPores);
    _grid.reserve(_numberOfPores);
    idx=0;
    for(unsigned int i=0; i< _numberOfPoints;i++)
    {
        if(_gridFull.at(i).medium==false)
        {
            SinglePoint cmd=_gridFull.at(i);
            _grid.push_back(cmd);
            _globalIndexing.at(idx)=i;
            idx++;
        }
    }
    return Output::Success;
}

//*************************************************************************************
Output::ExitCodes Points::ComputePorosity()
{
    _porosity=(1.0- (double)((double)_numberOfGrains/(double)_numberOfPoints))*100;
    return Output::Success;
}
//*************************************************************************************

int Points::FromPorousToGlob(int &porous)
{
    int id=-1;
    id=_globalIndexing.at(porous);
    return id;
}

//*************************************************************************************

int Points::FromGlobToPor(unsigned int &globIdx)
{
    auto it=find(_globalIndexing.begin(),_globalIndexing.end(),globIdx);
    if(it!=_globalIndexing.end())
    {
        return it-_globalIndexing.begin();
    }
    else
        return -1;
}


//*************************************************************************************

Output::ExitCodes Points::GetNeighs(SinglePoint &A, vector<SinglePoint *> &neighs)
{
      for(int x= -1; x<=1; x++)
    {
        for(int y=-1;y<=1;y++)
        {
            for(int z=-1;z<=1;z++)
            {
                if(x==0&& y==0 && z==0)
                    continue;

                int checkx=A.x+x;
                int checky=A.y+y;
                int checkz=A.z+z;
                if(checkx>=0 && checkx < (int)(_nx) && checky>=0 && checky < (int)(_ny) && checkz>=0 && checkz < (int)(_nz))
                {
                    unsigned int vicino=checkx+checky*_nx+checkz*_nx*_ny;
                    int id=FromGlobToPor(vicino);
                    if(id!=-1)
                        neighs.push_back(&_grid.at(id));
                }
            }
        }
    }
    return Output::Success;
}

//*************************************************************************************
unsigned int Points::GetPointIndex(unsigned int &x, unsigned int &y, unsigned int &z)
{
    unsigned int index=_numberOfPoints;
    for(unsigned int i=0; i< _numberOfPoints;i++)
    {
        if(_gridFull.at(i).x==(int)(x) && _gridFull.at(i).y ==(int)(y) && _gridFull.at(i).z==(int)(z))
        {
            index=i;
            break;
        }

    }
    if(index==_numberOfPoints)
    {
        cout <<" The Point is not in the collection "<<x<<" "<<y<<" "<<z<<" "<<endl;
    }
    return index;
}

//*************************************************************************************

Output::ExitCodes Points::ImportFromFile( const string &filename)
{
    ifstream inFile;
    inFile.open(filename.c_str());
    string line;
    if(!inFile.is_open())
    {
        cout<< " Cannot open file"<< endl;
        return Output::GenericError;
    }
    else
    {
        int numberLines=0;
        int memX=0, memY=0, memZ=0;
        _nx=1,_ny=1,_nz=1;
        while(getline(inFile,line))
        {
            numberLines++;
            vector<string> parser;
            parser=Input::Split(line,",");
            if(memX<atoi(parser.at(0).c_str()))
            {
                memX++;
                _nx++;
            }
            if(memY<atoi(parser.at(1).c_str()))
            {
                memY++;
                _ny++;
            }
            if(memZ<atoi(parser.at(2).c_str()))
            {
                memZ++;
                _nz++;
            }
        }
        _numberOfPoints=numberLines;

        _gridFull.reserve(_numberOfPoints);
        _gridFull.resize(_numberOfPoints);

        inFile.close();
        inFile.open(filename.c_str());
        ///Number of Header Lines has to be removed
        numberLines=0;

        ///Storing the points collection
        while(getline(inFile,line))
        {
            if(numberLines>=0)
            {
                vector<string> parser;

                parser=Input::Split(line,",");

                int X=atoi(parser.at(0).c_str()),Y=atoi(parser.at(1).c_str()),Z=atoi(parser.at(2).c_str()),index;

                index=X+Y*_nx+Z*_nx*_ny;
                _gridFull.at(index).x=X;
                _gridFull.at(index).y=Y;

                _gridFull.at(index).z=Z;
                _gridFull.at(index).medium=atoi(parser.at(3).c_str());

            }
            numberLines++;

        }

        inFile.close();
    }

    return Output::Success;
}

//**************************************************************************************************
Output::ExitCodes Points::ConfigureInletOutlet(const string &filename, const unsigned int& dir, bool inlet, unsigned int numData)
{

    ifstream inFile;
    inFile.open(filename.c_str());
    string line;
    if(!inFile.is_open())
    {
        return Output::GenericError;
    }
    else
    {
        vector<unsigned int> val(3,0);
        int numberLines=0;



        ///Storing the points collection
        while(getline(inFile,line) && numberLines<(int)(numData))
        {

            vector<string> parser;
            parser=Input::Split(line,"\t");
            for(int i=0;i<_numberProperties-1;i++)
            {
                val.at(i)=atoi(parser.at(i).c_str());
                val.at(i)--;

            }
            //- index remind
            //- In the inlet outlet file x replaces y
            unsigned int index = val[1]+val[0]*_nx+val[2]*_nx*_ny;
            if(inlet)
            {
                if(dir==1)
                {
                    _startSizes.at(0)=_startSizes.at(0)+1;
                    inlet_indices.at(0).push_back(index);
                }else if(dir==2)
                {
                    _startSizes.at(1)=_startSizes.at(1)+1;
                    inlet_indices.at(1).push_back(index);
                }else if(dir==3)
                {
                    _startSizes.at(2)=_startSizes.at(2)+1;
                    inlet_indices.at(2).push_back(index);
                }
            }
            else
            {
                if(dir==1)
                {
                    _targetSizes.at(0)=_targetSizes.at(0)+1;
                    target_indices.at(0).push_back(index);
                }else if(dir==2)
                {
                    _targetSizes.at(1)=_targetSizes.at(1)+1;
                    target_indices.at(1).push_back(index);
                }else if(dir==3)
                {
                    _targetSizes.at(2)=_targetSizes.at(2)+1;
                    target_indices.at(2).push_back(index);
                }

            }
            numberLines++;

        }
        inFile.close();
    }

    return Output::Success;
}

//**************************************************************************************************
Output::ExitCodes Points::ConfigureInletOutletVector()
{
    unsigned int couples=0;
    for(int i=0; i<_numberProperties-1;i++)
    {
        couples= couples+ _startSizes.at(i)*_targetSizes.at(i);
    }
    path_couples.reserve(couples*2);
    for(int i=0; i<_numberProperties-1;i++)
    {
        for(unsigned int j=0; j<_startSizes.at(i);j++)
        {
            for(unsigned int k=0; k<_targetSizes.at(i);k++)
            {
                int id=_gridFull.at(inlet_indices.at(i).at(j)).pourous_index;
                path_couples.push_back(id);
                id=_gridFull.at(target_indices.at(i).at(k)).pourous_index;

                path_couples.push_back(id);
            }
        }

    }

    return Output::Success;
}


//**************************************************************************************************
Output::ExitCodes Points::PrintPoint(unsigned int &idPoint)
{
    cout << _gridFull.at(idPoint).x<< " , "<<  _gridFull.at(idPoint).y << " , " << _gridFull.at(idPoint).z << " , " << _gridFull.at(idPoint).medium <<
             " , " << _gridFull.at(idPoint).pourous_index <<";";
    cout << endl;
    return Output::Success;
}


//**************************************************************************************************

bool find(SinglePoint &point, vector<SinglePoint*> &list)
{
    bool res=false;
    for(int it=list.size()-1;it>=0;it--)
    {
        if(point==*list.at(it))
        {
            res=true;
            break;
        }
    }
    return res;
}
//**************************************************************************************************

void equal(SinglePoint &A, SinglePoint &B)
{
    A.x=B.x;
    A.y=B.y;
    A.z=B.z;
    A.medium=B.medium;
    A.globidx=B.globidx;
    A.pourous_index=B.pourous_index;
    A._neighbourhood=B._neighbourhood;
    A._deadSuspect=B._deadSuspect;

    return;
}

//**************************************************************************************************

bool SamePoint(SinglePoint A, SinglePoint* B)
{

    if(A.x==B->x && A.y==B->y && A.z==B->z)
    {
        return true;
    }
    else return false;

}
bool SamePoint(SinglePoint A, SinglePoint B)
{
    if(A.x==B.x && A.y==B.y && A.z==B.z)
    {
        return true;
    }
    else return false;
}
//**************************************************************************************************

void PrintPoint(SinglePoint A, bool matlabFormat)
{
    int a=0;
    if(matlabFormat)
    {
        a=1;
        cout << A.y +a << " "<< A.x +a <<" "<< A.z+a /*<< " "<<A.medium*/<< endl;
    }
    else
        cout << A.x +a << " "<< A.y +a <<" "<< A.z+a /*<< " "<<A.medium*/<< endl;

    return;
}
//**************************************************************************************************


