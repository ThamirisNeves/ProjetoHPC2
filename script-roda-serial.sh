#!/bin/bash

#
#	Script para rodar os executáveis gerados no outro script
#	que são os compilados da temp 0-8 com as flags -OX
#

for j in $(seq 0 8)
  do
        for i in $(seq 0 3)
          do
                   isingO$i-T$j.x > saidaO$i-T$j.dat
          done
  done

