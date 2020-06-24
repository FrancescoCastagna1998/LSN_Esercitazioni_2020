CC = g++
CFLAGS = -Wall -O3 --std=c++11

Es1.2.exe : Es1.2.o random.o Vettore.o funzioni.o
	$(CC) random.o Es1.2.o Vettore.o funzioni.o -o Es1.2.exe
Es1.2.o : Es1.2.cpp
	$(CC) -c Es1.2.cpp -o Es1.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Vettore.o : Vettore.cxx Vettore.h
	$(CC) -c Vettore.cxx -o Vettore.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h Vettore.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
clean :
	rm *.o Es1.2.exe seed.out

clean_result:
	rm st1.txt st2.txt st10.txt st100.txt
	rm exp1.txt exp2.txt exp10.txt exp100.txt
	rm cl1.txt cl2.txt cl10.txt cl100.txt
