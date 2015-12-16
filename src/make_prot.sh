#!/bin/bash
c++ -std=c++11 -c protkmer_tester.cpp -o prot.o
g++ prot.o city.o -o prot
