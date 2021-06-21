#!/bin/sh
pgf95 -O3 -mp -acc -Minfo=accel -ta=tesla:cc60,cuda10.1 kernel_v1.F90 -o kernel_v1
pgf95 -O3 -mp -acc -Minfo=accel -ta=tesla:cc60,cuda10.1 kernel_v2.F90 -o kernel_v2
pgf95 -O3 -mp kernel_v2.F90 -o kernel_cpu
