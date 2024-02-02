CC = g++
CFLAGS = -Wall -O3 --std=c++11

Exercise9 : main.o random.o Genetic.o City.o
	$(CC) random.o main.o Genetic.o City.o -o Exercise9.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
City.o: City.cpp City.h
	$(CC) -c City.cpp -o City.o $(CFLAGS)
Genetic.o: Genetic.cpp Genetic.h
	$(CC) -c Genetic.cpp -o Genetic.o $(CFLAGS)	
clean :
	rm *.o Exercise9.exe
