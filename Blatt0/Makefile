CFLAGS = -Wall -Wextra
all:	three plot3 two plot2

three:	Aufgabe3.cpp | Data
		g++ $(CFLAGS) -O3 -ffast-math -march=native Aufgabe3.cpp -o Aufgabe3
		./Aufgabe3

two:	Aufgabe2.cpp | Data
		g++ $(CFLAGS) -O3 -ffast-math -march=native Aufgabe2.cpp -o Aufgabe2
		./Aufgabe2

plot3: three | Plots
	python Aufgabe3.py

plot2: two | Plots 
	python Aufgabe2.py

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
