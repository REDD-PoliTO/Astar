# Astar

This library is related to the results presented in the submitted article 'A parallel programming application of the A* algorithm in Digital Rock Physics', Computers & Geosciences https://www.sciencedirect.com/journal/computers-and-geosciences

GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007

```rb
git clone https://github.com/REDD-PoliTO/Astar.git
```

Compile and Run the first test. 

```rb
cd Astar
pwd
cwd=$(pwd)
sed -i "s@FOLDER_PATH_TO_REPLACE@$cwd@g" src_headers/auxMacros.hpp
```
src_headers/auxMacros.hpp: INPUTFOLDER has been modified with the current folder.

```rb
mkdir build
cd build
cmake ../src_exe/ .
```
Change the flag In CMake Environment if necessary
```rb
make
./astarSolver test161 5
```

replace ``` 5 ``` with the number of test inlet/outlet per direction (maximum 50). The given example will run 5 inlet points and 5 outlet points per direction. All the available threads will be used.

The test case ```test161 ``` is related to the Sample161 presented in the article. 
