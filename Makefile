MPI_CC=CC
ifeq ($(PE_ENV),PGI)
	MPI_FLAGS=-O3 -fast -acc -Minfo=acc -Mnoopenmp -Mcuda=maxregcount:72
	ifeq ($(CRAYPAT_COMPILER_OPTIONS),1)
		MPI_FLAGS+= -DCRAYPAT
	endif
else
	MPI_FLAGS=-O3 -hnoomp -hacc -hlist=m 
endif



SOURCES= wmBrick3D.cpp OpenChannel3D.cpp WallMountedBrick.cpp
OBJECTS=wmBrick3D.o OpenChannel3D.o WallMountedBrick.o workArounds.o
LIBS=

ifeq ($(USE_NVTX),1)
	LIBS+=-lnvToolsExt
	MPI_FLAGS+= -DUSE_NVTX
endif

TARGET=WMBrick3D

%.o: %.cpp
	$(MPI_CC) $(MPI_FLAGS) -c $^

$(TARGET): $(OBJECTS)
	$(MPI_CC) $(MPI_FLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(TARGET) *~
