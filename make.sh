#!/bin/bash
cd src
make && rm -f ../bin/kingassembler.py ../bin/megahit_gt && cp kingassembler.py megahit_gt ../bin
cd ../share/RDPTools/
make

