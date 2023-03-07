g++ -c exercise1_2.cpp -o exercise1_2.o
g++ -c random.cpp -o random.o
g++ exercise1_2.o random.o -o ex1_2 -larmadillo
