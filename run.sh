#!/bin/sh

export CUDA_VISIBLE_DEVICES="0"

echo
echo "Running version 1 of GPU kernels"
echo "================================="
PGI_COMPARE=abs=13 ./kernel_v1
echo
echo "Running version 2 of GPU kernel"
echo "==============================="
PGI_COMPARE=abs=13 ./kernel_v2
echo
echo "Running on cpu"
echo "======================="
./kernel_cpu
