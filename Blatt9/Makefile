CFLAGS = -Wall -Wextra
INCPATH = -I ~/usr/include/eigen3

all:  one plot1

one:	Aufgabe1.cpp | Data
		g++ $(CFLAGS) $(INCPATH) -O3 -ffast-math -march=native -std=c++11 Aufgabe1.cpp -o Aufgabe1
		./Aufgabe1

plot1: one | Plots
	python Aufgabe1.py

Plots:
	mkdir -p Plots

Data:
	mkdir -p Data

clean:
	
	rm -f Aufgabe1
	rm -rf Data/
	rm -rf Plots/

.PHONY: clean