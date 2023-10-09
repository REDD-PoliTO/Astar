# Astar

This code has been submitted to 
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

replace ``` 5 ``` with the number of test inlet/outlet per direction (maximum 50)
The example will run 5 inlet points and 5 outlet points per direction. All the available threads will be used.
