g++ -c exercise1_3.cpp -o exercise1_3.o
g++ -c random.cpp -o random.o
g++ exercise1_3.o random.o -o ex1_3 -larmadillo
