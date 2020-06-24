CC = g++
CFLAGS = -Wall -O3 --std=c++11

Random_Walk_3D.exe : Random_Walk_3D.o random.o Vettore.o funzioni.o Posizione.o
	$(CC) random.o Random_Walk_3D.o Vettore.o funzioni.o Posizione.o -o Random_Walk_3D.exe
Random_Walk_3D.o : Random_Walk_3D.cpp
	$(CC) -c Random_Walk_3D.cpp -o Random_Walk_3D.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Vettore.o : Vettore.cxx Vettore.h
	$(CC) -c Vettore.cxx -o Vettore.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h Vettore.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
Posizione.o : Posizione.cpp Posizione.h 
	$(CC) -c Posizione.cpp -o Posizione.o $(CFLAGS)
clean :
	rm *.o Random_Walk_3D.exe seed.out

clean_result:
	rm RW3Ddiscreto.txt
	rm RW3Dcontinuo.txt
