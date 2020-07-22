dependencies (Debian):
sudo apt install libsfml-dev
sudo apt install libcgal-dev

compile:
mkdir build
cd build
cmake ../
make

run example (view distance, scale, path to obj):
./ren -100 1 path_to.obj
