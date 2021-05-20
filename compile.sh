#!/bin/sh
pgf95 -O3 -acc -Minfo=accel -ta=tesla:cc60,cuda10.1 kernel_v1.f90 -o kernel_v1
pgf95 -O3 -acc -Minfo=accel -ta=tesla:cc60,cuda10.1 kernel_v2.f90 -o kernel_v2
pgf95 -O3 kernel_v2.f90 -o kernel_cpu
