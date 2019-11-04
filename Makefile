CC = gcc
CFLAGS= -Wall -Ofast -std=gnu99
LDFLAGS = -L/usr/lib
LIBS = -lusb-1.0 -lairspy -lpthread -lfftw3f -lcurl -lm

OBJS = airspy_wsprd.o wsprd.o wsprsim_utils.o wsprd_utils.o tab.o fano.o nhash.o wsprnet.o

%.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

airspy_wsprd: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

airspy_wsprd.o: filter.h airspy_wsprd.c

clean:
	rm -f *.o airspy_wsprd wspr_wisdom.dat hashtable.txt
