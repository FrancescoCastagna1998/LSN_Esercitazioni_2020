CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_ISING_1D.exe : Monte_Carlo_ISING_1D.o random.o
	$(CC) random.o Monte_Carlo_ISING_1D.o -o Monte_Carlo_ISING_1D.exe
Monte_Carlo_ISING_1D.o : Monte_Carlo_ISING_1D.cpp
	$(CC) -c Monte_Carlo_ISING_1D.cpp -o Monte_Carlo_ISING_1D.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_ISING_1D.exe seed.out

clean_Gibbs :
	rm *.o Monte_Carlo_ISING_1D.exe
	rm output.ene_0_equi
	rm output.mag_0_equi
	rm output.ene.0
	rm output.ene.0_final
	rm output.heat.0
	rm output.heat.0_final
	rm output.mag.0
	rm output.mag.0_final
	rm output.chi.0
	rm output.chi.0_final

clean_Metropolis :
	rm *.o Monte_Carlo_ISING_1D.exe
	rm output.ene_1_equi
	rm output.mag_1_equi
	rm output.ene.1
	rm output.ene.1_final
	rm output.heat.1
	rm output.heat.1_final
	rm output.mag.1
	rm output.mag.1_final
	rm output.chi.1
	rm output.chi.1_final
input_Gibbs :
	cp input.Gibbs input.dat
input_Metropolis :
	cp input.Metropolis input.dat

