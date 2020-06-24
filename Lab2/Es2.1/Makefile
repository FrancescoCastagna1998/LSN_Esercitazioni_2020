CC = g++
CFLAGS = -Wall -O3 --std=c++11

Integral_MC.exe : Integral_MC.o random.o Vettore.o funzioni.o
	$(CC) random.o Integral_MC.o Vettore.o funzioni.o -o Integral_MC.exe
Integral_MC.o : Integral_MC.cpp FunzioneBase.h Coseno.h CosenoIS.h
	$(CC) -c Integral_MC.cpp -o Integral_MC.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Vettore.o : Vettore.cxx Vettore.h
	$(CC) -c Vettore.cxx -o Vettore.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h Vettore.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
clean :
	rm *.o Integral_MC.exe seed.out

clean_result:
	rm int.txt int_IS.txt
