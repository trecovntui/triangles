MAINSRC = ren.cpp
MAINLIB = -lsfml-graphics -lsfml-window -lsfml-system

SRC = $(MAINSRC)

OBJ = $(SRC:.c=.o)

CC = g++

LIBS = -lCGAL -lgmp $(MAINLIB)

all: ren

.c.o:
	$(CC) -c $<

ren: $(OBJ)
	g++ -o $@ $^ $(LIBS)

ren_e:
	g++ -o ren_e -I ~/repos/eigen ren_e.cpp $(MAINLIB)

clean:
	rm -f ren *.o
	rm -f ren_e *.o
