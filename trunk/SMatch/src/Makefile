CXX =		/usr/local/mpich2/bin/mpicxx
CXXFLAGS =	-O3 -Wall -I/usr/local/boost/include -I/usr/local/mpich2/include -IC:/ASN/libs/nlopt/include

OBJS =		SMatch.o Mol.o Parser.o Printer.o Optimization.o Coord.o Engine.o

LIBS =		-LC:/ASN/libs/nlopt/lib -L/usr/local/boost/lib -lnlopt_cxx -lboost_mpi -lboost_serialization -lz #-static

TARGET =	SMatch

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
