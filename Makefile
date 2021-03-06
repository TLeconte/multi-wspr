CC = gcc
CFLAGS= -Wall -std=gnu99 -Ofast -march=native
LDFLAGS = -g -L/usr/lib
LIBS = -lusb-1.0 -lairspy -lpthread -lfftw3f -lcurl -lm 
OBJS = multi_wspr.o wsprd.o wsprsim_utils.o wsprd_utils.o tab.o fano.o nhash.o wsprnet.o airspy.o jelinek.o 

%.o: %.c
	${CC} ${CFLAGS} -c $< -o $@

multi_wspr: $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

multi_wspr.o: filter.h freqsets.h multi_wspr.c

clean:
	rm -f *.o multi_wspr
