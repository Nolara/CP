CFLAGS = -Wall -Wextra
all:	three plot3

three:	Aufgabe3.cpp
		g++ $(CFLAGS) -O3 -ffast-math -march=native Aufgabe3.cpp -o Aufgabe3
		./Aufgabe3

two:	Aufgabe2.cpp
		g++ $(CFLAGS) -O3 -ffast-math -march=native Aufgabe2.cpp -o Aufgabe2
		./Aufgabe2

plot3: three
	python Aufgabe3.py

plot2: two
	python Aufgabe2.py

clean:
	rm -f Aufgabe3
	rm -f Data/*
	rm -f Plots/*
