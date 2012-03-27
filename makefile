BIN = SkinSeg
OBJ = main.o skinFilter.o

CINCLUDE = `pkg-config --cflags opencv`
CLIB = `pkg-config --libs opencv`

all: $(BIN)

$(BIN): $(OBJ)
	gcc -g -o $(BIN) $(CLIB) $(OBJ)

main.o: main.c
	gcc -g $(CINCLUDE) -c main.c
skinFilter.o: skinFilter.c
	gcc -g $(CINCLUDE) -c skinFilter.c
clean: 
	rm $(OBJ) $(BIN)
