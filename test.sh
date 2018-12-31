#!/bin/bash
mpicxx -Wall -c ./dns.cpp
mpicxx -Wall -c ./matrix.cpp
mpicxx -Wall ./main.cpp ./dns.o ./matrix.o -o DNS_mul
rm -rf ./dns.o ./matrix.o
mpirun -np $1 ./DNS_mul 
