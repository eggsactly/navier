EXE:=navier
CC:=gcc
CFLAGS:=-O3 -Wall -pg -DAPP_NAME=\"$(EXE)\" -IlibBitmapWag/include/

.c.o:
	$(CC) -c $(CFLAGS) $<

.PHONY: all
all:    $(EXE)

.PHONY: clean
clean:
	rm -f *.o navier

$(EXE): alloc.o boundary.o init.o main.o output.o simulation.o \
libBitmapWag/lib/libBitmapWag.a
	$(CC) $(CFLAGS) -o $@ $^ -lm

boundary.o       : datadef.h
init.o           : datadef.h libBitmapWag/include/libBitmapWag.h
main.o           : alloc.h boundary.h datadef.h init.h simulation.h \
libBitmapWag/include/libBitmapWag.h

simulation.o     : datadef.h init.h

libBitmapWag/lib/libBitmapWag.a: libBitmapWag/libBitmapWag.c 
libBitmapWag/lib/libBitmapWag.a: libBitmapWag/libBitmapWag.h
	$(MAKE) -C libBitmapWag/

