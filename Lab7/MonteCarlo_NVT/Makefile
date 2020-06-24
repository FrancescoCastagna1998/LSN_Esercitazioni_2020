CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_NVT.exe : Monte_Carlo_NVT.o random.o
	$(CC) random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT.exe
Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_NVT.exe

input_solid :
	cp input.solid input.dat

input_liquid :
	cp input.liquid input.dat

input_gas :
	cp input.gas input.dat

import_fcc :
	cp config.fcc config.0

import_solid :
	cp input.solid input.dat
	cp config.solid config.0

import_liquid :
	cp input.liquid input.dat
	cp config.liquid config.0

import_gas :
	cp input.gas input.dat
	cp config.gas config.0

save_solid :
	cp config.final config.solid

save_liquid :
	cp config.final config.liquid

save_gas :
	cp config.final config.gas

clean_solid : 
	rm output.epot_equi_solid
	rm output.epot_solid
	rm output.pres_solid
	output.gofr_solid_NVT
	output.gave_solid_NVT

clean_liquid : 
	rm output.epot_equi_liquid
	rm output.epot_liquid
	rm output.pres_liquid
	output.gofr_liquid_NVT

clean_gas : 
	rm output.epot_equi_gas
	rm output.epot_gas
	rm output.pres_gas
	output.gofr_gas_NVT
	output.gave_solid_NVT
