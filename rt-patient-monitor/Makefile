# Real-time patient-monitor simulator.
#   make            - native build (Linux or self-hosted QNX)
#   make verify     - off-target verification suite -> results/verification.csv
#   make qnx        - cross build for QNX Neutrino (source qnxsdp-env.sh first); QNX_ARCH=aarch64le|x86_64
CC      ?= gcc
CFLAGS  ?= -std=gnu11 -O2 -g -Wall -Wextra -Wno-misleading-indentation
LDLIBS  := -lm -lpthread
CORE    := src/siggen.c src/ecg_proc.c src/spo2_proc.c src/alarm.c
RT      := src/monitor.c src/rt_util.c src/chan.c
HDRS    := $(wildcard src/*.h)

all: build/monitor build/loadgen build/test_detect
build/monitor: $(RT) $(CORE) $(HDRS) | build
	$(CC) $(CFLAGS) -o $@ $(RT) $(CORE) $(LDLIBS)
build/loadgen: src/loadgen.c | build
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)
build/test_detect: tests/test_detect.c $(CORE) $(HDRS) | build
	$(CC) $(CFLAGS) -o $@ tests/test_detect.c $(CORE) $(LDLIBS)
build:
	mkdir -p build results
verify: build/test_detect
	./build/test_detect results/verification.csv

QNX_ARCH ?= aarch64le
qnx:
	$(MAKE) CC="qcc -Vgcc_nto$(QNX_ARCH)" LDLIBS="-lm" all
clean:
	rm -rf build
.PHONY: all verify qnx clean
