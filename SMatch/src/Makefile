CXX =		g++ #/usr/local/mpich2/bin/mpicxx
CXXFLAGS =	-fopenmp -Ofast -Wall -I/usr/local/boost/include -I/usr/local/mpich2/include -IC:/ASN/libs/nlopt/include -IC:/ASN/libs/boost/include
OBJS =		SMatch.o Mol.o Parser.o Printer.o Optimization.o Coord.o Engine.o
LIBS =		-LC:/ASN/libs/nlopt/lib -L/usr/local/boost/lib -lnlopt_cxx -LC:/ASN/libs/zlib -lz -static -lgomp # -lboost_mpi -lboost_serialization
TARGET =	SMatch
$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)
all:	$(TARGET)

clean:
	del $(OBJS) $(TARGET)
