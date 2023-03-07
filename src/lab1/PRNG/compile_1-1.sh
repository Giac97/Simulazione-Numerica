g++ -c exercise1_1.cpp -o exercise1_1.o
g++ -c random.cpp -o random.o
g++ exercise1_1.o random.o -o ex1_1 -larmadillo
