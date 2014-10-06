//
//  fasta_decompressor.c
//  iDoComp_v1
//
//  Created by Mikel Hernaez on 8/18/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//
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

#include "iDoComp.h"


chromosome idc_decompress_chr_ints(fasta_compressor fc, unsigned int* numChr, unsigned int* bpPerLine, uint32_t **targetLength){
    
    unsigned int i = 0, prevInt = 0, ctrChr = 0;
    
    chromosome chrRoot, chr;
    
    chrRoot = create_chromosome(NULL, NULL, NULL, 0, 0, 0);
    
    
    // Start compressing the Ints
    
    // number of Chromosomes
    *numChr = decompress_Ints(fc->Ints);
    
    *targetLength = (uint32_t *) calloc(*numChr, sizeof(uint32_t));
    
    // number of base pairs per line
    *bpPerLine = decompress_Ints(fc->Ints);
    
    // [numInst numSubs numInse] x numChr
    chr = chrRoot;
    while (ctrChr++ < *numChr) {
        
        chr->numInst = decompress_Ints(fc->Ints);
        chr->instructions = (instruction) calloc(chr->numInst, sizeof(struct instruction_t));
        
        chr->numInse = decompress_Ints(fc->Ints);
        chr->insertions = (insertion) calloc(chr->numInse, sizeof(struct insertion_t));
        
        chr->numSubs = decompress_Ints(fc->Ints);
        chr->substitutions = (substitution) calloc(chr->numSubs, sizeof(struct substitution_t));
        
        chr->next = create_chromosome(NULL, NULL, NULL, 0, 0, 0);
        chr = chr->next;
    }
    chr = chrRoot;
    ctrChr = 0;
    //printf("num chars to decompress: %d\n",numChars);
    // Store the rest of the numbers related to the edit distance between each of the pair of strings
    while (ctrChr < *numChr) {
        
        //Decompress instructions
        prevInt = 0;
        for (i = 0; i < chr->numInst; i++) {
            
            // decompress the position using the info of the sign
            chr->instructions[i].pos =  (decompress_Signs(fc->Signs))? prevInt - decompress_Ints(fc->Ints): prevInt + decompress_Ints(fc->Ints);
            
            // decompress the length of the instruction
            chr->instructions[i].length = decompress_Ints(fc->Ints);
            
            prevInt = chr->instructions[i].pos + chr->instructions[i].length;
            
            targetLength[0][ctrChr] += chr->instructions[i].length, (targetLength[0][ctrChr]++);
        }
        
        // Decompress insertions
        prevInt = 0;
        for (i = 0; i < chr->numInse; i++) {
            
            chr->insertions[i].pos = prevInt + decompress_Ints(fc->Ints);
            
            prevInt = chr->insertions[i].pos;
            
            targetLength[0][ctrChr]++;
        }
        
        // Decompress substitutions
        prevInt = 0;
        for (i = 0; i < chr->numSubs; i++) {
            
            chr->substitutions[i].pos = prevInt + decompress_Ints(fc->Ints);
            
            prevInt = chr->substitutions[i].pos;
        }
        
        chr = chr->next;
        
        ctrChr++;
    }
    
    return chrRoot;
    
}

//**************************************************************//
//                                                              //
//                      CHROMOSOMES                             //
//                                                              //
//**************************************************************//

uint32_t reconstruct_instructions(chromosome chr, llnode T, char *Ref, fasta_compressor fc){
    
    uint32_t i, j, t_ctr = 0;
    char targetChar;
    
    // Do the instructions, except the last one (in the last one we do not reconstruct the chars)
    for (i = 0; i < chr->numInst-1; i++) {
        
        // Reconstruct the "matched" letters
        for (j = chr->instructions[i].pos; j < chr->instructions[i].pos + chr->instructions[i].length; j++){
            T[t_ctr].value = Ref[j-1];
            T[t_ctr].next = &T[t_ctr+1], t_ctr++;
        }
        
        // The next character of the reference needs to be modified by the new char
        targetChar = decompress_Chars(fc->Chars, Ref[j-1]);
        
        T[t_ctr].value = targetChar;
        T[t_ctr].next = &T[t_ctr+1], t_ctr++;
        
    }
    // Do the last instruction, no chars
    for (j = chr->instructions[i].pos; j < chr->instructions[i].pos + chr->instructions[i].length; j++){
        T[t_ctr].value = Ref[j-1];
        T[t_ctr].next = &T[t_ctr+1], t_ctr++;
    }
    
    return t_ctr;
}


char* reconstruct_insertions(chromosome chr, llnode T, uint32_t *currentLength, fasta_compressor fc){
    
    uint32_t i, j, targetLength, index, currentPos;
    char targetChar;
    llnode Tidx;
    char *target;
    
    Tidx = T;
    
    // Create an array of length equal to the target length
    targetLength = *currentLength + chr->numInse;
    *currentLength = targetLength;
    
    
    target = calloc(targetLength, sizeof(char));
    
    // Traverse the current T and do the insertions while copying the target to the new array
    currentPos = 0;
    for (i = 0; i < chr->numInse; i++) {
        
        j = chr->insertions[i].pos;
        
        for (index = currentPos; index<j-1; index++){
            target[index] = Tidx->value, Tidx = Tidx->next;
        }
        currentPos = index;
        
        // We add the insertion
        targetChar = decompress_Chars(fc->Chars, Tidx->value);
        
        // Do the insertion
        target[currentPos++] = targetChar;
    }
    
    // traverse from the position of the last insertion to the end of the target
    for (index = currentPos; index<targetLength; index++){
        target[index] = Tidx->value, Tidx = Tidx->next;
    }
    
    return target;
    
    
}

void reconstruct_substitutions(chromosome chr, char* T, fasta_compressor fc){
    
    uint32_t i, j;
    char targetChar;
    
    // Do the substitutions, except the last one (in the last one we do not reconstruct the chars)
    
    for (i = 0; i < chr->numSubs; i++) {
        
        j = chr->substitutions[i].pos;
        
        // The next character of the reference needs to be modified by the new char
        targetChar = decompress_Chars(fc->Chars, T[j-1]);
        
        // Do the substitution
        T[j-1] = targetChar;
    }
    
    
}

char* reconstruct_target(char *Ref, chromosome chr, fasta_compressor fc, uint32_t *targetSize){
    
    uint32_t currentLength;
    llnode T_ll;
    char *target;
    
    T_ll = (llnode)calloc(*targetSize, sizeof(struct llnode_t));
    
    // printf("Reconstruct instructions \n");
    currentLength = reconstruct_instructions(chr, T_ll, Ref, fc);
    
    
    // printf("Reconstruct insertions \n");
    target = reconstruct_insertions(chr, T_ll, &currentLength, fc), free(T_ll);
    
    // printf("Reconstruct substitutions \n");
    reconstruct_substitutions(chr, target, fc);
    
    *targetSize = currentLength;
    
    return target;
}

uint8_t write_chromosomes_to_file(char *T, FILE *headerFile, FILE * targetFile, uint32_t bpPerLine, uint32_t targetLength){
    
    uint32_t numBps;
    
    char header[MAX_HEADER];
    
    // Read the header from the headerFile
    fgets(header, MAX_HEADER, headerFile);
    
    if (header[0] != '\n'){
        fputs(header, targetFile);
    }
    
    if (bpPerLine == 0) bpPerLine = targetLength;
    
    for (numBps = 0; numBps < targetLength; numBps++){
        fputc(T[numBps], targetFile);
        if ( ((numBps + 1)%bpPerLine) == 0) fputc('\n', targetFile);
    }
    
    return 1;
}


void start_fasta_decompression(char* osPath, FILE *inputFile, FILE * headerFile){
    
    unsigned int ctrChr = 0, numChr = 0, bpPerLine = 0;
    
    char refPath[256], targetPath[256];
    
    char *S = NULL, *T = NULL;
    
    unsigned int *targetSize = NULL;
    
    fasta_compressor fc;
    
    chromosome chr;
    
    FILE *refFile, *targetFile;
    
    // Initialize the decompressor
    
    fc = initialize_fasta_compressor(osPath, DECOMPRESSION, NULL);
    
    chr = idc_decompress_chr_ints(fc, &numChr, &bpPerLine, &targetSize);
    
    // Go through the Chromosomes and reconstruct the target
    while (ctrChr < numChr){
        
        // printf("Decompressing chromosome %02d\n",chrCtr++);
        // Read reference and target from file
        fscanf(inputFile, "%s %s", refPath, targetPath);
        while (refPath[0] == '#') { // This is the input to the mapping function
            fscanf(inputFile, "%s %s", refPath, targetPath);}
        
        // Open the Ref file
        refFile = fopen ( refPath , "r" );
        if (refFile==NULL) {fputs ("Reference File error\n",stderr); exit (1);}
        // Store Ref sequence in memory
        store_reference_in_memory(refFile, &S);
        fclose (refFile);
        
        // Open the Target file
        targetFile = fopen ( targetPath , "w" );
        if (targetFile==NULL) {fputs ("Target File error\n",stderr); exit (1);}
        
        
        // Reconstruct target chromosome and read chars at the same time
        T = reconstruct_target(S, chr, fc, &(targetSize[ctrChr]));
        
        // Reconstruct the Chromosome
        write_chromosomes_to_file(T, headerFile, targetFile, bpPerLine, (targetSize[ctrChr]));
        
        //free and close file
        fclose(targetFile);
        free(S);
        free(T);
        
        chr = chr->next;
        
        ctrChr++;
    }
}
