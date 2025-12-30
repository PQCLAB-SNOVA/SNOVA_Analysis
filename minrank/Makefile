ifeq ($(ARCH),)
ARCH = native
endif

CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -Wredundant-decls -Wshadow -Wvla -Wpointer-arith -O3 -march=$(ARCH) -mtune=$(ARCH)

OBJS += ./build/fips202.o
HDRS = 
LIBS =


all: build analyze

clean:
	rm -rf build __pycache__

build:
	mkdir -p build

./build/%.o: %.c $(HDRS)
	$(CC) $(CFLAGS) $(P) -o $@ -c $<

./build/analyze: $(OBJS) ./build/beullens_attack.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

analyze: build ./build/analyze
	@./build/analyze 0 100000

./build/ranktest: $(OBJS) ./build/ranktest.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

ranktest: build ./build/ranktest
	@./build/ranktest
