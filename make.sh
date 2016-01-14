#!/bin/bash
cd src
make && rm -f ../bin/megahit_gt.py ../bin/megahit_gt && cp megahit_gt.py megahit_gt ../bin
cd ../share/RDPTools/
make

