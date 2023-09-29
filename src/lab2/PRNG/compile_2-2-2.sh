g++ -c exercise2_2_2.cpp -o exercise2_2_2.o
g++ -c random.cpp -o random.o
g++ -c utils.cpp -o utils.o
g++ exercise2_2_2.o utils.o random.o -o ex2_2_2 -larmadillo
rm *.o
