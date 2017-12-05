CXX=g++
CXXFLAGS=-std=c++11 -O2
CPPFLAGS=
INCLUDE=/home/James/codes/cpp-opt/include
LIBS=-lnlopt

all: howarth-dorodnitsyn-coupled

howarth-dorodnitsyn-coupled: howarth-dorodnitsyn-coupled.cpp 
	$(CXX) -o $@ $(CXXFLAGS) -I$(INCLUDE) $(CPPFLAGS) $< $(LIBS)
