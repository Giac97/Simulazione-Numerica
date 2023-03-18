g++ -c exercise2_2.cpp -o exercise2_2.o
g++ -c random.cpp -o random.o
g++ exercise2_2.o random.o -o ex2_2 -larmadillo
rm *.o
