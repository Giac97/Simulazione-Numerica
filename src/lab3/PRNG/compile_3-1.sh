g++ -c exercise3_1.cpp -o exercise3_1.o
g++ -c random.cpp -o random.o
g++ -c utils.cpp -o utils.o
g++ exercise3_1.o utils.o random.o -o ex3_1 -larmadillo
rm *.o
