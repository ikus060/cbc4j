#Makefile for javaCbc

$(LIBJAVACBC): javacbc.o javacbc_wrap.o
	$(TARGET)-g++ -shared -DNDEBUG javacbc.o javacbc_wrap.o -L$(CBCDIR)/lib -lCbcSolver -lCbc -lCgl -lOsiClp -lOsiCbc -lOsi -lClp -lCoinUtils -lm $(LDFLAGS) -o $(LIBJAVACBC)

javacbc.o: javacbc.cpp javacbc.hpp
	$(TARGET)-g++ -fPIC -c javacbc.cpp -I$(CBCDIR)/include/coin -I$(JAVADIR)/include/ -I$(JAVADIR)/include/linux -I.

javacbc_wrap.o: javacbc_wrap.cxx
	$(TARGET)-g++ -fPIC -c javacbc_wrap.cxx -I$(CBCDIR)/include/coin -I$(JAVADIR)/include/ -I$(JAVADIR)/include/linux -I.

javacbc_wrap.cxx: javacbc.i javacbc.hpp
	swig -outdir "$(SWIG_OUTDIR)" -package "$(SWIG_PACKAGE)" -c++ -java javacbc.i

all: $(LIBJAVACBC)
