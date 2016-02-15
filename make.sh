#!/bin/bash
cd src
make && rm -f ../bin/megagta.py ../bin/megagta && cp megagta.py megagta ../bin
cd ../share/RDPTools/
make

