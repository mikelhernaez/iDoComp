//
//  idc_load_chr.c
//  iDoComp_v1
//
//  Created by Mikel Hernaez on 8/7/14.
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


#include "idc_load_chr.h"

//**************************************************************//
//                                                              //
//                  STORE REFERENCE IN MEMORY                   //
//                                                              //
//**************************************************************//
int store_reference_in_memory(FILE* refFile, char** S){
    uint32_t ch, letterCount;
    char *H;
    
    *S = (char *) malloc(MAX_BP_CHR*sizeof(char));
    
    // ******* Read and Store Reference****** //
    letterCount = 0;
    
    ch = getc(refFile);
    if (ch == '>' || ch == '@'){
        H = (char*)malloc(MAX_HEADER * sizeof(char));
        fgets(H, MAX_HEADER, refFile);
        free(H);
    }
    else{
        (*S)[letterCount++] = ch;
    }
    
    while (EOF != (ch=getc(refFile)))
        if (ch !='\n')
            (*S)[letterCount++] = ch;
    
    (*S)[letterCount] = '#';
    (*S)[letterCount+1] = '\0';
    
    (*S) = (char *) realloc(*S, letterCount+2);
    
    return letterCount+1;
    
}

//**************************************************************//
//                                                              //
//                  STORE TARGET IN MEMORY                      //
//                                                              //
//**************************************************************//
int store_target_in_memory(FILE* targetFile, char** T, FILE* headerFile, uint32_t *bpPerLine){
    
    uint32_t ch, letterCount;
    char *H;
    
    *T = (char *) malloc(MAX_BP_CHR*sizeof(char));
    
    // ******* Read and Store Target ****** //
    letterCount = 0;
    
    // Check to see if there is a header in the file
    ch = getc(targetFile);
    if (ch == '>' || ch == '@'){
        H = (char*)malloc(MAX_HEADER * sizeof(char));
        *H = ch;
        fgets(H+1, MAX_HEADER, targetFile);
        fputs(H, headerFile);
        free(H);
    }
    else{
        fputs("\n", headerFile);
        (*T)[letterCount++] = ch;
    }
    
    
    // Read and Store the rest
    while (EOF != (ch=getc(targetFile)))
        if (ch !='\n')
            (*T)[letterCount++] = ch;
        else{
            *bpPerLine = letterCount;
            break;
        }
    
    while (EOF != (ch=getc(targetFile)))
        if (ch !='\n')
            (*T)[letterCount++] = ch;
    
    // Terminate the target with '$
    (*T)[letterCount] = '$';
    (*T)[letterCount+1] = '\0';
    
    (*T) = (char *) realloc(*T, letterCount+2);
    
    return letterCount+2;
    
    
}



chromosome generate_chromosome_sequence(FILE *inputFile, FILE *headerFile, unsigned int *chrCtr, unsigned int *bpPerLine){
    
    char refPath[256];
    char targetPath[256];
    char saPath[256];
    
    FILE *targetFile, *refFile, *saFile;
    
    unsigned int refLength, targetLength;
    
    uint64_t totalDNAsize = 0;
    
    char *S = NULL, *T = NULL;
    
    int *SA;
    
    chromosome chrRoot, chr;
    
    chr = create_chromosome(NULL, NULL, NULL, 0, 0, 0);
    
    chrRoot = chr;
    
    *bpPerLine = 0, *chrCtr = 0;
    
    while (fscanf(inputFile, "%s %s %s", refPath, targetPath, saPath) != EOF){
        
        if (refPath[0] != '#') { // This is the input to the mapping function
            
            printf("//********* CHROMOSOME %02u ************//\n", *chrCtr);
            
            // Open the Ref file
            refFile = fopen ( refPath , "r" );
            if (refFile==NULL) {fputs ("Reference File error\n",stderr); return 0;}
            // Store Ref sequence in memory
            refLength = store_reference_in_memory(refFile, &S);
            fclose (refFile);
            
            // Open the Target file
            targetFile = fopen ( targetPath , "r" );
            if (targetFile==NULL) {fputs ("Target File error\n",stderr); return 0;}
            // Store Target sequence in memory
            targetLength = store_target_in_memory(targetFile, &T, headerFile, bpPerLine);
            fclose (targetFile);
            
            // Open the SA file
            saFile = fopen ( saPath , "r" );
            if (saFile==NULL) {fputs ("SA File error\n",stderr); return 0;}
            // Store SA in memory and return a pointer to the root
            SA = (int *)malloc((size_t)refLength * sizeof(int));
            fread(SA, refLength, sizeof(int), saFile);
            fclose(saFile);
            
            // Generate the mapping between the Target and the Reference
            chr->next = generate_chromosome_mapping(*chrCtr, S, T, targetLength, refLength, SA);
            chr = chr->next;
                
            // Free the rest of the memory that we dont need
            free(S); free(T); free(SA);
            
            // Update the total number of integers to compress
            //*numTotalInt += 2*chr->numInst + chr->numSubs + chr->numInse;
            
            totalDNAsize += targetLength;
            
            (*chrCtr)++;
            
            //end = clock();
            //printf("Time spent computing chromosome %02u: %.2f seconds\n", chrCtr, (double)(end - begin)/ CLOCKS_PER_SEC);
        }
    }
    
    //(*numTotalInt) += 2 + 3 * (*chrCtr); // 1 for numChr, 1 for bp/line, and each chromosome numIns, numInse, numSub
    
    return chrRoot->next;
}
