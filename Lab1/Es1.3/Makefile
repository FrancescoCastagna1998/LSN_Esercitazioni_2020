CC = g++
CFLAGS = -Wall -O3 --std=c++11

Buffon.exe : Buffon.o random.o
	$(CC) random.o Buffon.o -o Buffon.exe
Buffon.o : Buffon.cpp Buffon.h
	$(CC) -c Buffon.cpp -o Buffon.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Buffon.exe seed.out
clean_result:
	rm pigreco.txt

