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

clean:
	rm -f ren *.o
