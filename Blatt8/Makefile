CFLAGS = -Wall -Wextra -g
INCPATH = -I ~/usr/include/eigen3
all:	one plot1 Animations

one:	Aufgabe1.cpp | Data
		g++ $(CFLAGS) $(INCPATH) -O3 -ffast-math -march=native -std=c++11 Aufgabe1.cpp -o Aufgabe1
		./Aufgabe1

plot1: one | Plots
	python Aufgabe1.py


Animations: one | Anim
		python Animation_Programs/Animate_0_01_aqu.py
		python Animation_Programs/Animate_0_01.py
		python Animation_Programs/Animate_1_aqu.py
		python Animation_Programs/Animate_1.py
		python Animation_Programs/Animate_10_aqu.py
		python Animation_Programs/Animate_10.py
		python Animation_Programs/Animate_isokin_aqu.py
		python Animation_Programs/Animate_isokin.py

Plots:
	mkdir -p Plots

Data:
	mkdir -p Data

Anim:
	mkdir -p Animations

clean:
	rm -f Aufgabe1
	rm -rf Data/
	rm -rf Plots/
	rm -rf Animations/

.PHONY: clean
