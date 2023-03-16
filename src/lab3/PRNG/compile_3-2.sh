g++ -c exercise3_2.cpp -o exercise3_2.o
g++ -c random.cpp -o random.o
g++ -c utils.cpp -o utils.o
g++ exercise3_2.o utils.o random.o -o ex3_2 -larmadillo
rm *.o
