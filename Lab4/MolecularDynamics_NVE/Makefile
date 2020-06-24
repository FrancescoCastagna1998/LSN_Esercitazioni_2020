CC = g++
CFLAGS = -Wall -O3 --std=c++11

MolDyn_NVE.exe : MolDyn_NVE.o Vettore.o funzioni.o
	$(CC) MolDyn_NVE.o Vettore.o funzioni.o -o MolDyn_NVE.exe
MolDyn_NVE.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
Vettore.o : Vettore.cxx Vettore.h
	$(CC) -c Vettore.cxx -o Vettore.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h Vettore.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
clean :
	rm *.o MolDyn_NVE.exe 

import_fcc :
	cp config.fcc config.0
import_liquid :
	cp input.liquid input.dat
import_solid :
	cp input.solid input.dat
import_Argon_solid :
	cp input_Argon.solid input.dat
import_Argon_liquid :
	cp input_Argon.liquid input.dat
import_Argon_gas :
	cp input_Argon.gas input.dat
clean_solid :
	rm output_temp_istant_solid.dat
	rm output_temp_istant_equi_solid.dat
	rm ave_epot_solid.out
	rm ave_ekin_solid.out
	rm ave_temp_solid.out
	rm ave_etot_solid.out
clean_liquid :
	rm output_temp_istant_liquid.dat
	rm output_temp_istant_equi_liquid.dat
	rm ave_epot_liquid.out
	rm ave_ekin_liquid.out
	rm ave_temp_liquid.out
	rm ave_etot_liquid.out
clean_Argon_solid :
	rm output_temp_istant_Argon_solid.dat
	rm output_temp_istant_equi_Argon_solid.dat
	rm ave_epot_Argon_solid.out
	rm ave_ekin_Argon_solid.out
	rm ave_temp_Argon_solid.out
	rm ave_etot_Argon_solid.out
save_Argon_solid :
	cp config.0 config.Argon_solid
	cp old.0 old.Argon_solid
config_Argon_solid :
	cp config.Argon_solid config.0
	cp old.Argon_solid old.0 
clean_Argon_liquid :
	rm output_temp_istant_Argon_liquid.dat
	rm output_temp_istant_equi_Argon_liquid.dat
	rm ave_epot_Argon_liquid.out
	rm ave_ekin_Argon_liquid.out
	rm ave_temp_Argon_liquid.out
	rm ave_etot_Argon_liquid.out
save_Argon_liquid :
	cp config.0 config.Argon_liquid
	cp old.0 old.Argon_liquid
config_Argon_liquid :
	cp config.Argon_liquid config.0
	cp old.Argon_liquid old.0 
clean_Argon_gas :
	rm output_temp_istant_Argon_gas.dat
	rm output_temp_istant_equi_Argon_gas.dat
	rm ave_epot_Argon_gas.out
	rm ave_ekin_Argon_gas.out
	rm ave_temp_Argon_gas.out
	rm ave_etot_Argon_gas.out
save_Argon_gas :
	cp config.0 config.Argon_gas
	cp old.0 old.Argon_gas
config_Argon_gas :
	cp config.Argon_gas config.0
	cp old.Argon_gas old.0 



