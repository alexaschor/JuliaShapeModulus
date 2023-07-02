SHELL := /bin/bash -e
LIBS = Quaternion

include ./projects/include.mk

all: sdfGen main genPoly

sdfGen:
	cd projects/sdfGen; make

main:
	cd projects/main; make

genPoly:
	cd projects/genPoly; make

clean:
	@-for d in $(LIBS); do (echo -e "cd ./lib/$$d; rm *.o";cd ./lib/$$d; rm *.o; cd ../..); done
	cd projects/main; make clean
	cd projects/genPoly; make clean
	cd projects/sdfGen; make clean
