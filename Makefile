CXX=g++
CXXFLAGS= -std=c++11
DEPS= space.h atom.h interface.h polarconfig.h autospeed.h
LIBPATH =
CXXFLAGS +=
CXXFLAGS +=-I./include
VPATH=./src
ana.x: atom.o main.o space.o interface.o polarconfig.o autospeed.o
	mkdir -p obj bin
	$(CXX) -o ana.x $(LIBPATH) atom.o main.o space.o interface.o polarconfig.o autospeed.o
	mv *.o obj
	mv ana.x bin
%.o: %.c $(DEPS)
	$(CXX)  $(CXXFLAGS) -c -o $@ $^
clean:
	rm -rf *.o obj bin
