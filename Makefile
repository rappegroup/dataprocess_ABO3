CXX=icc
CXXFLAGS= -std=c++11
DEPS= space.h atom.h
LIBPATH =-L/opt/gcc/7.3.0/snos/lib64
CXXFLAGS += -I/opt/gcc/7.3.0/snos/include
CXXFLAGS +=-I./include
VPATH=./src
ana.x: atom.o main.o space.o
	mkdir -p obj bin
	$(CXX) -o ana.x $(LIBPATH) atom.o main.o space.o
	mv *.o obj
	mv ana.x bin
%.o: %.c $(DEPS)
	$(CXX)  $(CXXFLAGS) -c -o $@ $^
clean:
	rm -rf *.o obj bin
