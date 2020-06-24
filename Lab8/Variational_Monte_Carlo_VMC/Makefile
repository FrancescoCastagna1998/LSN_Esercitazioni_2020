CC = g++
CFLAGS = -Wall -O3 --std=c++11

Variational_Monte_Carlo.exe : Variational_Monte_Carlo.o random.o TrialWaveFunction.o 
	$(CC) random.o Variational_Monte_Carlo.o TrialWaveFunction.o -o Variational_Monte_Carlo.exe
Variational_Monte_Carlo.o : Variational_Monte_Carlo.cpp 
	$(CC) -c Variational_Monte_Carlo.cpp -o Variational_Monte_Carlo.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
TrialWaveFunction.o : TrialWaveFunction.cpp TrialWaveFunction.h  FunzioneBase.h
	$(CC) -c TrialWaveFunction.cpp -o TrialWaveFunction.o $(CFLAGS)
clean :
	rm *.o Variational_Monte_Carlo.exe seed.out

clean_result:
	rm output.expection_values_hamiltonian_equi
	rm output.expection_values_hamiltonian
	rm ProbabilityGS.dat





