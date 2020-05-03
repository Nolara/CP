CFLAGS = -Wall -Wextra
INCPATH = -I ~/usr/include/eigen3
all:	one plot1 two plot2 three plot3

one:	Aufgabe1.cpp | Data
		g++ $(CFLAGS) $(INCPATH) -O3 -ffast-math -march=native -std=c++11 Aufgabe1.cpp -o Aufgabe1
		./Aufgabe1

two:	Aufgabe2.cpp Profiler.cpp Profiler.h | Data
		g++ $(CFLAGS) $(INCPATH) -O3 -ffast-math -march=native -std=c++11 Aufgabe2.cpp profiler.cpp -o Aufgabe2
		./Aufgabe2

three:	Aufgabe3.cpp | Data
				g++ $(CFLAGS) $(INCPATH) -O3 -ffast-math -march=native -std=c++11 Aufgabe3.cpp -o Aufgabe3
				./Aufgabe3

plot1: one | Plots
	python Aufgabe1.py

plot2: two | Plots
	python Aufgabe2.py

plot3: three | Plots
	python Aufgabe3.py

Plots:
	mkdir -p Plots

Data:
	mkdir -p Data

clean:
	rm -f Aufgabe1
	rm -f Aufgabe2
	rm -f Aufgabe3
	rm -rf Data/
	rm -rf Plots/

.PHONY: clean
