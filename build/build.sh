#!/bin/bash
export SDKROOT=$(xcrun --show-sdk-path)
g++ -I ../third-party/eigen-3.4.0 -I ../third-party/spectra-1.1.0/include -I ../third-party/MshIO-0.0.1/include -l mshio -L ../third-party/MshIO-0.0.1/build -O3 ../src/eigenfem/main.cpp -o main.o
exit 0