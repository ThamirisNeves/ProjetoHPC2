#!/bin/bash

#
#	Script para mudar a temperatura no código variando de 0.5 até 8
#	de 0.5 em 0.5 e depois compilar o código com as flags -OX 
#

for i in $(seq 0.1 0.1 8)
  do 
	j=3
	k=5
	sed "s/XX/$i/g" isingteste.c > isingparalelo.c
        mpicc -O$j -ffast-math isingparalelo.c -o isingparalelo-T$i.x -lm
	mpirun -np $k ./isingparalelo-T$i.x  
	cp tempo.dat time-core$k-T$i.dat
	cp media.dat media-core$k-T$i.dat
	
  done

