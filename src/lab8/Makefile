CC = g++
CFLAGS = -Wall -O3 --std=c++11

Exercise8.exe : Exercise_8.o random.o
	$(CC) random.o Exercise_8.o -o Exercise8.exe
Exercise_8.o : Exercise_8.cpp
	$(CC) -c Exercise_8.cpp -o Exercise_8.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Exercise8.exe