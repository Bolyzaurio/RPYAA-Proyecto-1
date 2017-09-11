main:
	gcc nasa_shuttle.c -lgsl -lm -Wall -pedantic -o Main
	./Main
clean :
	rm Main 
