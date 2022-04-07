SHELL := /bin/bash -e
LIBS = Quaternion

include ./projects/include.mk

all: sdfGen main

sdfGen:
	cd projects/sdfGen; make

main:
	cd projects/main; make

clean:
	@-for d in $(LIBS); do (echo -e "cd ./lib/$$d; rm *.o";cd ./lib/$$d; rm *.o; cd ../..); done
	cd projects/sdfGen; make clean
	cd projects/main; make clean
