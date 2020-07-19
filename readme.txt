dependencies (Debian):
sudo apt install libsfml-dev
sudo apt install libcgal-dev

compile:
g++ ren.cpp -lCGAL -lgmp -lsfml-graphics -lsfml-window -lsfml-system

for eigen (ren_e.cpp):
g++ -I ~/repos/eigen ren_e.cpp -lsfml-graphics -lsfml-window -lsfml-system

gnuplot:
set yrange [*:*] reverse
plot 'test_plot.txt' matrix with image reverse
