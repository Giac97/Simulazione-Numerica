g++ -c exercise_5_1_2p.cpp -o exercise5.o
g++ -c random.cpp -o random.o
g++ exercise5.o random.o -o ex5_2p -larmadillo
rm *.o
