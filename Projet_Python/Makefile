PROGS= main

all: $(PROGS)

CXXFLAGS: -Wall


main: main.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

clean:
	-rm *.o $(PROGS)

# les dependances

main.o:maillage.hpp
