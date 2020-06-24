CC = g++
CFLAGS = -Wall -O3 --std=c++11

Hydrogen_Atom.exe : Hydrogen_Atom.o random.o Posizione.o FunOndaHI2P.o FunOndaHIGS.o
	$(CC) random.o Hydrogen_Atom.o Posizione.o FunOndaHI2P.o FunOndaHIGS.o -o Hydrogen_Atom.exe
Hydrogen_Atom.o : Hydrogen_Atom.cpp Hydrogen_Atom.h Posizione.h 
	$(CC) -c Hydrogen_Atom.cpp -o Hydrogen_Atom.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Posizione.o : Posizione.cpp Posizione.h 
	$(CC) -c Posizione.cpp -o Posizione.o $(CFLAGS)
FunOndaHI2P.o : FunOndaHI2P.cpp FunOndaHI2P.h  FunzioneBase.h
	$(CC) -c FunOndaHI2P.cpp -o FunOndaHI2P.o $(CFLAGS)
FunOndaHIGS.o : FunOndaHIGS.cpp FunOndaHIGS.h  FunzioneBase.h
	$(CC) -c FunOndaHIGS.cpp -o FunOndaHIGS.o $(CFLAGS)

clean :
	rm *.o Hydrogen_Atom.exe seed.out

input_GS_Unif:
	cp input.GS_Unif input.dat

input_2P_Unif:
	cp input.2P_Unif input.dat

input_GS_Gauss:
	cp input.GS_Gauss input.dat

input_2P_Gauss:
	cp input.2P_Gauss input.dat

clean_GS_Unif:
	rm output.r_equi_GS_Unif
	rm output.rAve_GS_Unif
	rm output.XYZsampled_GS_Unif

clean_2P_Unif:
	rm output.r_equi_2P_Unif
	rm output.rAve_2P_Unif
	rm output.XYZsampled_2P_Unif

clean_GS_Gauss:
	rm output.r_equi_GS_Gauss
	rm output.rAve_GS_Gauss
	rm output.XYZsampled_GS_Gauss

clean_2P_Gauss:
	rm output.r_equi_2P_Gauss
	rm output.rAve_2P_Gauss
	rm output.XYZsampled_2P_Gauss
