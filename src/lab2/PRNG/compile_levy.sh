g++ -c exercise2_2_2_extra.cpp -o exercise2_2_2_e.o
g++ -c random.cpp -o random.o
g++ -c utils.cpp -o utils.o
g++ exercise2_2_2_e.o utils.o random.o -o levy -larmadillo
rm *.o
