FORT = gfortran
FLAGS = -Wunused-parameter -mcmodel=medium -O3

version := $(lastword $(shell head -1 illum/__init__.py))

SRCS := $(wildcard illum/kernel/*.f*)
BINS := $(SRCS:illum/kernel/%.f=bin/%)
BINS := $(BINS:illum/kernel/%.f90=bin/%)
LIBS := $(wildcard illum/kernel/libs/*.f*)

.PHONY: all clean echo

all: ${BINS}

echo:
	@echo ${BINS}

bin/illum-health: illum/kernel/illum-health.f90
	@sed -i "s/__version__/${version}/" $^
	${FORT} ${FLAGS} $^ ${LIBS} -o $@
	@sed -i "s/${version}/__version__/" $^

bin/%: illum/kernel/%.f
	${FORT} ${FLAGS} $^ -o $@

bin/%: illum/kernel/%.f90
	${FORT} ${FLAGS} $^ -o $@

clean:
	rm ${BINS}
