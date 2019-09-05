CC	 = clang++
CFLAGS	 = -std=c++11 -I/Users/clifton/Dropbox/Grad_school/entropy/lammps_EF/source/include/
#LDFLAGS = -l

traj: traj.o traj_reader.o
	$(CC) -o $@ $^

traj.o: traj.cpp test.h
	$(CC) $(CFLAGS) -c -g traj.cpp

traj_reader.o: traj_reader.cpp test.h
	$(CC) $(CFLAGS) -c -g traj_reader.cpp

clean:
	rm *.o
