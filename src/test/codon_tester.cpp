#include "codon.h"
#include <iostream>
#include <stdint.h>

using namespace std;

int main(int argc, char **argv) {
    uint8_t a = 4, b = 4, c = 4;
    char emission = Codon::codonTable[a][b][c];

    cout << "emission	= " << emission << '\n';

    return 0;
}