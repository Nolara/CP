CFLAGS = -Wall -Wextra
all:	three plot3

three:	Aufgabe3.cpp

		g++ $(CFLAGS) -O3 -ffast-math -march=native Aufgabe3.cpp -o Aufgabe3
		./Aufgabe3
plot3: Aufgabe3.py
	python Aufgabe3.py


clean:
	rm -f Aufgabe3
	rm -f Data/three.txt
