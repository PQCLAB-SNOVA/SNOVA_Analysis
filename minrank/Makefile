ifeq ($(ARCH),)
ARCH = native
endif

CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -Wredundant-decls -Wshadow -Wvla -Wpointer-arith -O3 -march=$(ARCH) -mtune=$(ARCH)

OBJS += ./build/fips202.o
HDRS =
LIBS =


all: build beullens_16

clean:
	rm -rf build __pycache__

build:
	mkdir -p build

./build/%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) $(P) -o $@ -c $<

./build/beullens_q: $(OBJS) ./build/beullens_q.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

beullens_q: build ./build/beullens_q
	@./build/beullens_q -1 1

./build/beullens_16: $(OBJS) ./build/beullens_16.o $(LIBS)

beullens_16: build ./build/beullens_16
	@./build/beullens_16 -1 1
