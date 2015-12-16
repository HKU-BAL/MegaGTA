#!/bin/bash
g++ -c kmer_tester.cpp -o tester.o
g++ -c city.cpp -o city.o
g++ -o test tester.o city.o
