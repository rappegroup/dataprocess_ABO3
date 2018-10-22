CXX=g++
CXXFLAGS= -std=c++11
DEPS= space.h atom.h interface.h polarconfig.h autospeed.h
LIBPATH =
vpath %.cpp src
vpath %.h include
vpath %.o obj
OBJDIR=obj
SRCDIR=src
INCDIR=include
CXXFLAGS+=-I$(INCDIR)
FFT_INC= -I/usr/local/include
FFT_PATH=-L/usr/local/lib
FFT_LIB=-lfftw3
ana.x: atom.o main.o space.o interface.o polarconfig.o autospeed.o
	mkdir -p obj bin
	$(CXX) -o ana.x $(LIBPATH) atom.o main.o space.o interface.o polarconfig.o autospeed.o
	mv *.o obj
%.o: %.c $(DEPS)
	$(CXX)  $(CXXFLAGS) $(FFT_INC) -c -o $@ $^
clean:
	rm -rf *.o obj bin
