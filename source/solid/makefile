CC	 = clang++
CFLAGS	 = -std=c++11 -I/Users/clifton/Dropbox/Grad_school/entropy/lammps_EF/source/include/
#LDFLAGS = -l

SFunk_solid: traj.o traj_reader.o read_config.o
	$(CC) -o $@ $^

traj.o: traj.cpp test.h
	$(CC) $(CFLAGS) -c -g traj.cpp

traj_reader.o: traj_reader.cpp traj_reader.h 
	$(CC) $(CFLAGS) -c -g traj_reader.cpp

read_config.o: read_config.cpp read_config.h
	$(CC) $(CFLAGS) -c -g read_config.cpp

clean:
	rm *.o
