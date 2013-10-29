CXX =		/opt/intel/composer_xe_2013_sp1.0.080/bin/intel64/icpc -fast#g++ 
CXXFLAGS =	-O2 -g3 -Wall -fmessage-length=0 -IC:/ASN/libs/nlopt/include

OBJS =		SMatch.o Mol.o Parser.o Printer.o Optimization.o Coord.o

LIBS =		-LC:/ASN/libs/nlopt/lib -L/usr/local/intel/nlopt/lib -lnlopt_cxx -lz

TARGET =	SMatch

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
