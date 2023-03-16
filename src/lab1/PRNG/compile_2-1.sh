g++ -c exercise2_1.cpp -o exercise2_1.o
g++ -c random.cpp -o random.o
g++ exercise2_1.o random.o -o ex2_1 -larmadillo
