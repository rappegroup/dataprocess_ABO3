CXX=icpc
CXXFLAGS= -std=c++11
DEPS= space.h atom.h interface.h polarconfig.h autospeed.h
LIBPATH =
vpath %.cpp src
vpath %.h include
vpath %.o obj
OBJDIR=./obj
SRCDIR=./src
INCDIR=./include
CXXFLAGS+=-I$(INCDIR)
FFT_INC = -I/opt/cray/pe/fftw/3.3.6.3/sandybridge/include
FFT_PATH = -L/opt/cray/pe/fftw/3.3.6.3/sandybridge/lib
FFT_LIB = -lfftw3
CXXFLAGS+= $(FFT_INC)
LIBPATH+= $(FFT_PATH)
ana.x: atom.o main.o space.o interface.o polarconfig.o autospeed.o
	mkdir -p obj bin
	$(CXX) -o ana.x $(LIBPATH) $(FFT_LIB) atom.o main.o space.o interface.o polarconfig.o autospeed.o
	mv *.o obj
%.o: %.c $(DEPS)
	$(CXX)  $(CXXFLAGS) -c -o $@ $^
clean:
	rm -rf *.o obj bin
