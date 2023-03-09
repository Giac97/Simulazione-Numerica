g++ -c exercise2_1_1.cpp -o exercise2_1_1.o
g++ -c random.cpp -o random.o
g++ -c utils.cpp -o utils.o
g++ exercise2_1_1.o utils.o random.o -o ex2_1_1 -larmadillo
rm *.o
