CC = g++
CFLAGS = -Wall -O3 --std=c++11

European_Option.exe : European_Option.o random.o
	$(CC) random.o European_Option.o -o European_Option.exe
European_Option.o : European_Option.cpp European_Option.h
	$(CC) -c European_Option.cpp -o European_Option.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o European_Option.exe
clean_result:
	rm Call_sampling_final_asset_price.txt Put_sampling_final_asset_price.txt
	rm Call_sampling_dis.txt Put_sampling_dis.txt
