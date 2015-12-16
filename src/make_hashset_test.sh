#!/bin/bash
g++ -c hash_set_tester.cpp -o hash_set_tester.o -fopenmp
g++ -c city.cpp -o city.o
g++ -o hash_test hash_set_tester.o city.o -fopenmp

