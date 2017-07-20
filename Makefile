CXX=/opt/apps/compilers/gcc/6.2.0/bin/g++
CXXFLAGS=-Wall -O2 -std=c++11 -DNDEBUG
OBJ=misc.o main.o

%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@

GreedyRelated: $(OBJ)
		$(CXX) -static $^ -o $@
