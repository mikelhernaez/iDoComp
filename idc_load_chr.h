//
//  idc_load_chr.h
//  iDoComp_v1
//
//  Created by Mikel Hernaez on 8/7/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU Affero General Public License version 3,
//as published by the Free Software Foundation.
//
//This program is distributed in the hope that it will be useful, but
//WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
//Affero General Public License for more details.
//
//You should have received a copy of the GNU Affero General Public
//License along with this program. If not, see
//<http://www.gnu.org/licenses/>.

#ifndef iDoComp_v1_idc_load_chr_h

#define iDoComp_v1_idc_load_chr_h

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


// Define memory allocation parameters
#define MAX_BP_CHR 300000000
#define MAX_HEADER 500

// Define the Substitution parameters
// first filter
#define MAX_LEN_SUB 400
#define MIN_LEN_GAP 10

//second filter
#define MIN_SNP 25
#define lenAlwaysSubs 20
#define MIN_GAP 0

typedef struct instruction_t{
    uint32_t pos;
    uint32_t length;
    char refChar;
    char targetChar;
}*instruction;

typedef struct substitution_t{
    uint32_t pos;
    char refChar;
    char targetChar;
}*substitution;

typedef struct insertion_t{
    uint32_t pos;
    char refChar;
    char targetChar;
}*insertion;

typedef struct chromosome_t{
    instruction instructions;
    substitution substitutions;
    insertion insertions;
    struct chromosome_t *next;
    
    uint32_t numInst;
    uint32_t numSubs;
    uint32_t numInse;
    
}*chromosome;

typedef struct llnode_t{
    char value;
    struct llnode_t * next;
} *llnode;

int array_compare(char*S, char* T, uint32_t* T_len);

int store_reference_in_memory(FILE* refFile, char** S);
int store_target_in_memory(FILE* targetFile, char** T, FILE* headerFile, uint32_t *bpPerLine);
chromosome generate_chromosome_sequence(FILE *inputFile, FILE *headerFile, unsigned int* chrCtr, unsigned int *numTotalInt);

chromosome create_chromosome(instruction inst, substitution subs, insertion inse, uint32_t numInst, uint32_t numSubs, uint32_t numInse);
instruction generate_instructions(char *S, char *T, uint32_t m, uint32_t n, int *SA, uint32_t *numInstructions);
uint32_t create_substitution(substitution substitutions, uint32_t currentSubPos, char prevRefChar, char prevTargetChar, uint32_t firstSub);
insertion generate_insertions(instruction* oldIns, uint32_t* numInstructions, uint32_t *numInsertions);
chromosome generate_chromosome_mapping(uint8_t chrCtr, char *S, char *T, uint32_t tagetLength, uint32_t refLength, int * SA);


#endif
