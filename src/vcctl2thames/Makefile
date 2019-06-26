MAKE= make
CC= gcc
CFLAGS= -c -Wall -g -O2 -I/usr/local/include -I/usr/include -I.
LD_LIBS= -L/usr/local/lib -L/usr/lib
RM= /bin/rm -f
CP= /bin/cp -f

all:	vcctl2thames

%.o:	%.c
	$(CC) $(CFLAGS) $<

vcctl2thames: vcctl2thames.o vcctl.h
	$(CC) -o $@ vcctl2thames.o

clean:
	${RM} -r -R vcctl2thames
	${RM} -r -R vcctl2thames.o

cleanobj:
	${RM} -r -R vcctl2thames.o
