CC	 = clang++
#CFLAGS	 = -I
#LDFLAGS = -l

all: sfunk 

sfunk: main.o sfunk.o rdf.o sim.o
	$(CC) -o $@ $^ #$(LDFLAGS)

main.o: main.cpp
	$(CC) -c -g main.cpp

sfunk.o: SFunk.cpp SFunk.h sfunk_const.h 
	$(CC) -c -g SFunk.cpp

rdf.o: rdf.cpp rdf.h
	$(CC) -c -g rdf.cpp

sim.0: sim.cpp sim.h
	$(CC) -c -g sim.cpp


clean:
	rm *.o
	


