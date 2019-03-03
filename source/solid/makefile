CC	 = clang++
CFLAGS	 = -std=c++11 -I/Users/clifton/Dropbox/Grad_school/entropy/lammps_EF/source/include/
#LDFLAGS = -l

all: test 

#test: pointcloud.o 
#	$(CC) -o $@ $^ #$(LDFLAGS)

main.o: main.cpp
	$(CC) -c -g main.cpp

test: test.o 
	$(CC) -o $@ $^ #$(LDFLAGS)

test.o: test.cpp test.h traj.h 
	$(CC) $(CFLAGS) -c -g test.cpp

traj: traj.o traj_reader.o
	$(CC) -o $@ $^

traj.o: traj.cpp test.h 
	$(CC) $(CFLAGS) -c -g traj.cpp 

traj_reader.o: traj_reader.cpp test.h
	$(CC) $(CFLAGS) -c -g traj_reader.cpp


gauss_fit.o: gauss_fit.cpp
	$(CC) $(CFLAGS) -c -g gauss_fit.cpp

a.out: gauss_fit.cpp
	$(CC) $(CFLAGS) -g gauss_fit.cpp

clean:
	rm *.o
	


