# Makefile for compiling sim.o
CXX := g++

sim.o: sim.cpp sim.h
	$(CXX) $< -o $@

clean:
	rm -f *.o *.txt

.PHONY: clean