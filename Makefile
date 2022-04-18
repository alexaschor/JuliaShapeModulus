SHELL := /bin/bash -e
LIBS = Quaternion

include ./projects/include.mk

all: sdfGen main rootMorph genPoly

sdfGen:
	cd projects/sdfGen; make

main:
	cd projects/main; make

rootMorph:
	cd projects/rootMorph; make

genPoly:
	cd projects/genPoly; make

clean:
	@-for d in $(LIBS); do (echo -e "cd ./lib/$$d; rm *.o";cd ./lib/$$d; rm *.o; cd ../..); done
	cd projects/sdfGen; make clean
	cd projects/main; make clean
	cd projects/rootMorph; make clean
	cd projects/genPoly; make clean
