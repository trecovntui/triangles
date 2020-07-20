dependencies (Debian):
sudo apt install libsfml-dev
sudo apt install libcgal-dev

compile:
g++ ren.cpp -lCGAL -lgmp -lsfml-graphics -lsfml-window -lsfml-system
[or]
make
