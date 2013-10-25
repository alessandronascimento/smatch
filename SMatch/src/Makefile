CXX = 	g++ #/opt/intel/composer_xe_2013_sp1.0.080/bin/intel64/icpc
CXXFLAGS =	-O2 -g3 -Wall -fmessage-length=0

OBJS =		SMatch.o Mol.o Gaussian.o Coord.o Optimization.o Printer.o

LIBS =		/usr/local/intel/nlopt/lib/libnlopt_cxx.a -lz

TARGET =	SMatch

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
