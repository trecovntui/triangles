LINUX:

Dependencies (Debian):
sudo apt install libsfml-dev
sudo apt install libcgal-dev

Compile:
mkdir build
cd build
cmake ../
make

Run (view distance, scale, path to obj):
./ren -100 1 path_to.obj


WINDOWS (Using MSVC):

Dependencies:
1. Download SFML, CGAL and Boost (required by CGAL) source.
2. Set include, lib (.lib) and bin (.dll) paths correctly in CMakeLists.txt.

Compile:
mkdir build
cd build
cmake ..\
msbuild ren.sln

Run (view distance, scale, path to obj):
ren.exe -100 1 path_to.obj