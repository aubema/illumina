FORT = gfortran
FLAGS = -Wunused-parameter -mcmodel=medium -O3

version := $(lastword $(shell head -1 illum/__init__.py))

SRCS := $(wildcard illum/kernel/*.f)
BINS := $(SRCS:illum/kernel/%.f=bin/%)
LIBS := $(wildcard illum/kernel/libs/*.f)

.PHONY: all clean echo

all: ${BINS}

bin/illumina: illum/kernel/illumina.f
	sed -i "s/__version__/${version}/" $^
	${FORT} ${FLAGS} $^ ${LIBS} -o $@
	sed -i "s/${version}/__version__/" $^

bin/%: illum/kernel/%.f
	${FORT} ${FLAGS} $^ -o $@

clean:
	rm ${BINS}
