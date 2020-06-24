CC = g++
CFLAGS = -Wall -O3 --std=c++11

Es1.1.exe : Es1.1.o random.o Vettore.o funzioni.o
	$(CC) random.o Es1.1.o Vettore.o funzioni.o -o Es1.1.exe
Es1.1.o : Es1.1.cpp
	$(CC) -c Es1.1.cpp -o Es1.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Vettore.o : Vettore.cxx Vettore.h
	$(CC) -c Vettore.cxx -o Vettore.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h Vettore.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
clean :
	rm *.o Es1.1.exe seed.out

clean_result:
	rm ave.txt
	rm variance.txt
	rm chi_2.txt
