#!/bin/bash

#GNU- no flags
for i in $(seq 0.5 0.5 8)
  do    
        x=$(date "+%s")
        mpirun -np $i ../../../siesta < benzene.fdf > saida-$i.out
        rm *.DM *.XV *.CG *.ZM
        y=$(date "+%s")
	echo $i  $(($y-$x)) >> time.txt
  done


