#CC = g++
#CFLAGS = -std-c++11
CC	 = clang++
CFLAGS	 = -std=c++11 -I/Users/clifton/Dropbox/Grad_school/entropy/lammps_EF/source/include/
#LDFLAGS = -l

SFunk_solid: traj.o traj_reader.o read_config.o neighbor_list_generator.o PBCPointCloud_generator.o PBC.o
	$(CC) -o $@ $^

traj.o: traj.cpp test.h
	$(CC) $(CFLAGS) -c -g traj.cpp

traj_reader.o: traj_reader.cpp traj_reader.h
	$(CC) $(CFLAGS) -c -g traj_reader.cpp

neighbor_list_generator.o: neighbor_list_generator.cpp neighbor_list_generator.h
	$(CC) $(CFLAGS) -c -g neighbor_list_generator.cpp

PBC.o: PBC.cpp PBC.h
	$(CC) $(CFLAGS) -c -g PBC.cpp

PBCPointCloud_generator.o: PBCPointCloud_generator.cpp PBCPointCloud_generator.h
	$(CC) $(CFLAGS) -c -g PBCPointCloud_generator.cpp

read_config.o: read_config.cpp read_config.h
	$(CC) $(CFLAGS) -c -g read_config.cpp

clean:
	rm *.o SFunk_solid
