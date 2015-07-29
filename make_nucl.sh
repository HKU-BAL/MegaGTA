#!/bin/bash
c++ -std=c++11 -c nuclkmer_tester.cpp -o nucl.o
g++ nucl.o city.o -o nucl
