#Make

hist.pdf: data5 distr.py
	python3 distr.py

data5: cuenta.x
	./distr.x

distr.x: distr.c
	gcc -fopenmp distr.c -o distr.x
