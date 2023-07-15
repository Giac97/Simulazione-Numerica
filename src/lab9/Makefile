CC = g++
CFLAGS = -Wall -O3 --std=c++11

Exercise9 : main.o random.o
	$(CC) random.o main.o -o Exercise8.exe
Exercise_8.o : Exercise_8.cpp
	$(CC) -c Exercise_8.cpp -o Exercise_8.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Exercise8.exe
