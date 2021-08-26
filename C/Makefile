CFLAGS=-O3 -Wall -Wextra -Wpedantic -std=c99 \
	-Wundef -Wshadow -Wcast-align -Wpointer-arith -Wmissing-prototypes \
	-fstrict-aliasing -fno-common -pipe -Wno-unused-function

SOURCES=randombytes.c poly.c common.c zq.c
HEADERS=randombytes.h poly.h common.h zq.h

.PHONY: all clean
all: test_01schoolbook test_02karatsuba test_03toom test_04ntt test_05fft \
     test_06radix3fft test_07incomplete test_08goods

test_%: %.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $< $(SOURCES) -D_16BIT_COEFFICIENTS -lm

clean:
	$(RM) test_*
