# Astar

```rb
git clone https://github.com/REDD-PoliTO/Astar.git
```

Compile and Run the first test. 

```rb
cd Astar
pwd
```

copy the current folder path and replace it in
../src_headers/auxMacros.hpp: INPUTFOLDER

```rb
mkdir build
cd build
cmake ../src_exe/ .
```
Change the flag In CMake Environment if necessary
```rb
make
./astarSolver test161 #n
```

replace ``` #n ``` with the number of test inlet/outlet per direction (maximum 50)
