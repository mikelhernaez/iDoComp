//
//  iDoComp.h
//  iDoComp_v1
//
//  Created by Mikel Hernaez on 8/7/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#ifndef iDoComp_v1_iDoComp_h

#define iDoComp_v1_iDoComp_h

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "idc_load_chr.h"

#include "fasta_compressor.h"


uint32_t start_fasta_compression(chromosome chrRoot, char* osPath, unsigned int numChr, unsigned int bpPerLine);
void start_fasta_decompression(char* osPath, FILE *inputFile, FILE * headerFile);
BASEPAIR char2BP(char c);
BASEPAIR BP2char(BASEPAIR c);
uint8_t** generate_char_transitions(chromosome chr);

#endif
