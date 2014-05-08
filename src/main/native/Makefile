#Makefile for cbc4j

$(LIBCBC4J): cbc4j.o cbc4j_wrap.o
	$(TARGET)-g++ -shared -DNDEBUG cbc4j.o cbc4j_wrap.o -L$(CBCDIR)/lib -lCbcSolver -lCbc -lCgl -lOsiClp -lOsiCbc -lOsi -lClp -lCoinUtils -lm $(LDFLAGS) -o $(LIBCBC4J)

cbc4j.o: cbc4j.cpp cbc4j.hpp
	$(TARGET)-g++ -fPIC -c cbc4j.cpp -I$(CBCDIR)/include/coin -I$(JAVADIR)/include/ -I$(JAVADIR)/include/linux -I.

cbc4j_wrap.o: cbc4j_wrap.cxx
	$(TARGET)-g++ -fPIC -c cbc4j_wrap.cxx -I$(CBCDIR)/include/coin -I$(JAVADIR)/include/ -I$(JAVADIR)/include/linux -I.

cbc4j_wrap.cxx: cbc4j.i cbc4j.hpp
	swig -outdir "$(SWIG_OUTDIR)" -package "$(SWIG_PACKAGE)" -c++ -java cbc4j.i

all: $(LIBCBC4J)
