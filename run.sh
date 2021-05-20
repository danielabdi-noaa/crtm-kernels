#!/bin/sh
echo "Running version 1"
echo "======================="
PGI_COMPARE=abs=13 ./kernel_v1
echo "Running version 2"
echo "======================="
PGI_COMPARE=abs=13 ./kernel_v2
echo "Running on cpu"
echo "======================="
./kernel_cpu
