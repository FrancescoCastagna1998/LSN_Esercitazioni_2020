CC = g++
CFLAGS = -Wall -O3 --std=c++11

Travelling_Salesman_Problem_TSP.exe : Travelling_Salesman_Problem_TSP.o random.o city.o Individual.o Population.o funzioni.o
	mpicxx random.o city.o Individual.o Population.o funzioni.o Travelling_Salesman_Problem_TSP.o -o Travelling_Salesman_Problem_TSP.exe
Travelling_Salesman_Problem_TSP.o : Travelling_Salesman_Problem_TSP.cpp Travelling_Salesman_Problem_TSP.h 
	mpicxx -c Travelling_Salesman_Problem_TSP.cpp -o Travelling_Salesman_Problem_TSP.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Individual.o : Individual.cpp Individual.h city.h 
	$(CC) -c Individual.cpp -o Individual.o $(CFLAGS)
city.o : city.cpp city.h 
	$(CC) -c city.cpp -o city.o $(CFLAGS)
Population.o : Population.cpp Population.h Individual.h 
	$(CC) -c Population.cpp -o Population.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h  
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
clean :
	rm *.o Travelling_Salesman_Problem_TSP.exe seed.out

input_square:
	cp input.square input.dat
	
clean_result:	
	rm Starting_Population_square_0.dat
	rm Starting_Population_square_1.dat
	rm Starting_Population_square_2.dat
	rm Starting_Population_square_3.dat
	rm Best_Individual_Evolution_square_0.dat
	rm Best_Individual_Evolution_square_1.dat
	rm Best_Individual_Evolution_square_2.dat
	rm Best_Individual_Evolution_square_3.dat
	rm Best_Individual_square_0.dat
	rm Best_Individual_square_1.dat
	rm Best_Individual_square_2.dat
	rm Best_Individual_square_3.dat
	rm Average_Individuals_Evolution_square_0.dat
	rm Average_Individuals_Evolution_square_1.dat
	rm Average_Individuals_Evolution_square_2.dat
	rm Average_Individuals_Evolution_square_3.dat
