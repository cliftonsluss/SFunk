CC	 = clang++
CFLAGS	 = -std=c++11
#CFLAGS	 = -I
#LDFLAGS = -l

all: SFunk_liquid

SFunk_liquid: main.o read_config.o rdf_reader.o funks.o
	#sfunk.o rdf.o sim.o
	$(CC) -o $@ $^

main.o: main.cpp
	$(CC) $(CFLAGS) -c -g main.cpp

rdf_reader.o: rdf_reader.cpp rdf_reader.h
	$(CC) $(CFLAGS) -c -g rdf_reader.cpp

funks.o: funks.cpp funks.h
	$(CC) $(CFLAGS) -c -g funks.cpp

read_config.o: ../solid/read_config.cpp ../solid/read_config.h
	$(CC) $(CFLAGS) -c -g ../solid/read_config.cpp

#sfunk.o: SFunk.cpp SFunk.h sfunk_const.h
#	$(CC) -c -g SFunk.cpp

#rdf.o: rdf.cpp rdf.h
#	$(CC) -c -g rdf.cpp

#sim.0: sim.cpp sim.h
#	$(CC) -c -g sim.cpp


clean:
	rm *.o SFunk_liquid
