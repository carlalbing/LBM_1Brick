MPI_CC=CC
MPI_FLAGS=-fast -acc -Minfo=acc -Mnoopenmp
#MPI_FLAGS=-O3 -hnoomp -hacc -hlist=m
#MPI_FLAGS=-O2


SOURCES= wmBrick3D.cpp OpenChannel3D.cpp WallMountedBrick.cpp
OBJECTS=wmBrick3D.o OpenChannel3D.o WallMountedBrick.o

TARGET=WMBrick3D

%.o: %.cpp
	$(MPI_CC) $(MPI_FLAGS) -c $^

$(TARGET): $(OBJECTS)
	$(MPI_CC) $(MPI_FLAGS) -o $@ $^

clean:
	rm *.o $(TARGET) *~
