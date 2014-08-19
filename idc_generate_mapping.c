//
//  idc_generate_mapping.c
//  iDoComp_v1
//
//  Created by Mikel Hernaez on 8/7/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "idc_load_chr.h"

//**************************************************************//
//                                                              //
//                  CREATE CHROMOSOME                           //
//                                                              //
//**************************************************************//
chromosome create_chromosome(instruction inst, substitution subs, insertion inse, uint32_t numInst, uint32_t numSubs, uint32_t numInse){
    
    chromosome chr;
    
    chr = (chromosome ) malloc(sizeof(struct chromosome_t));
    
    chr->insertions = inse;
    chr->instructions = inst;
    chr->substitutions = subs;
    
    chr->numSubs = numSubs;
    chr->numInse = numInse;
    chr->numInst = numInst;
    
    chr->next = NULL;
    
    
    return chr;
}

//**************************************************************//
//                                                              //
//                  ARRAY COMPARE FUNCTION                      //
//                                                              //
//**************************************************************//

int array_compare(char*S, char* T, uint32_t* T_len){
    
    int ctr = 0;
    
    while ( (S[ctr]^T[ctr]) == 0){ctr++;}
    
    *T_len = ctr;
    
    if (T[ctr] > S[ctr]) {
        return 1;
    }
    
    return 0;
    
}

//**************************************************************//
//                                                              //
//                  GENERATE INSTRUCTIONS                       //
//                                                              //
//**************************************************************//
instruction generate_instructions(char *S, char *T, uint32_t m, uint32_t n, int *SA, uint32_t *numInstructions){
    
    uint32_t ctr, i, T_len, prevCTR, insCtr, l, r, mid, T_len_1;
    
    int32_t T_ini = 0, LAST_LOW;
    
    uint8_t letterComp, stopCause;
    
    char ch, ch1;
    
    //printf("Size of the target string: %d Mbp\n", m / 1000000 + 1 );
    
    // Allocate enought memory for the instructions
    instruction instructions;
    instructions = (instruction) malloc(m*sizeof(struct instruction_t));
    
    ctr = 0;
    prevCTR = 0;
    i = 0;
    letterComp = 0;
    stopCause = 0;
    
    insCtr = 0;
    
    while (ctr < m ) {
        
        l = 0; r = n-1; mid = 0, T_len = 0, T_len_1 = 0, LAST_LOW = 0;
        
        // We are choosing the first instruction
        while (l < r) {
            mid = (l+r)/2;
            if (array_compare( (S+SA[mid]),(T+ctr), &T_len) == 1) {
                l = mid + 1;
                LAST_LOW = 1;
            }
            else{
                r = mid;
                LAST_LOW = 0;
            }
            
        }
        
        // We need to make sure that we are taking the LONGEST substring
        if(l){
            if(LAST_LOW){
                array_compare( (S+SA[l]),(T+ctr), &T_len_1 );
                T_ini = (T_len>T_len_1)? SA[l-1]: SA[l];
            }
            else{
                array_compare( (S+SA[l-1]),(T+ctr), &T_len_1);
                T_ini = (T_len>T_len_1)? SA[l]: SA[l-1];
            }
        }
        
        ctr += (T_len>T_len_1)? T_len: T_len_1;
        T_len = (T_len>T_len_1)? T_len: T_len_1;
        
        
        // Store the Instruction
        //     printf("%c %d\n", T[ctr], T_len);
        
        if (T[ctr] == '$' && T_len == 0) {
            instructions[insCtr].pos = instructions[insCtr-1].pos + instructions[insCtr-1].length + 1;
            instructions[insCtr].length = 0;
            instructions[insCtr].targetChar = 'N';
            instructions[insCtr].refChar = 'N';
            insCtr++;
            break;
        }
        else if (T[ctr] == '$') {// is the last instruction
            instructions[insCtr].pos = T_ini+1;
            instructions[insCtr].length = T_len;
            instructions[insCtr].targetChar = 'N';
            instructions[insCtr].refChar = 'N';
            insCtr++;
            break;
        }
        
        else if (T_len == 0){ // the target letter doesn't exist in the reference
            
            // Check if we are in the case where the target has N's, but the reference has none
            instructions[insCtr].pos = insCtr ? instructions[insCtr-1].pos + instructions[insCtr-1].length + 1 : 1;
            instructions[insCtr].length = 0;
            instructions[insCtr].targetChar = T[ctr];
            instructions[insCtr].refChar = S[instructions[insCtr].pos - 1];
            insCtr++;
        }
        
        else{
            
            ch = T[ctr];
            ch1 = S[T_ini + T_len];
            
            instructions[insCtr].pos = T_ini+1;
            instructions[insCtr].length = T_len;
            instructions[insCtr].targetChar = ch;
            instructions[insCtr].refChar = ch1;
            
            insCtr++;
        }
        
        
        // start all over again but with the updated T counter (ctr)
        ctr++; // we need to skip the mismatched letter as we explicitily saved it
        prevCTR = ctr;
    }
    
    *numInstructions = insCtr;
    
    // Realloc the Instructions array to free unused memory
    instructions = (instruction) realloc(instructions, insCtr*sizeof(struct instruction_t));
    
    return instructions;
}


//**************************************************************//
//                                                              //
//                  GENERATE SUBSTITUTIONS                      //
//                                                              //
//**************************************************************//
uint32_t create_substitution(substitution substitutions, uint32_t currentSubPos, char prevRefChar, char prevTargetChar, uint32_t firstSub){
    
    static uint32_t subCtr = 0;
    
    if (firstSub) {
        subCtr = 0;
    }
    
    else{
        substitutions[subCtr].pos = currentSubPos + 1;
        substitutions[subCtr].refChar = prevRefChar;
        substitutions[subCtr++].targetChar = prevTargetChar;
    }
    return subCtr;
}


substitution generate_substitutions(instruction* oldIns, uint32_t* numInstructions, uint32_t *numSubs, char* S, uint32_t S_length){
    
    uint32_t numNewInstructions, prevPos, prevLen, currNewInsLength, i, currentSubPos, currentPos, currentLen, posNewInstruction, ctr = 0;
    
    char prevRefChar, prevTargetChar;
    uint8_t doInstruction;
    uint32_t tempNumSubs, futurePos;
    
    // Allocate the Array for the substitutions
    substitution substitutions;
    substitutions = (substitution) malloc((MAX_LEN_SUB + 1)*(*numInstructions)*sizeof(struct substitution_t));
    
    // Allocate the Array for the new instructions
    instruction newInstructions;
    newInstructions = (instruction) malloc((*numInstructions)*sizeof(struct instruction_t));
    
    
    create_substitution(NULL, 0, 0, 0,1);
    numNewInstructions = 0;
    
    
    // Read first old instruction
    prevPos = (*oldIns)[0].pos;
    prevLen = (*oldIns)[0].length;
    
    //Start of the first new instruction
    posNewInstruction = prevPos;
    currNewInsLength = 0;
    
    //Start with a new substitution
    currentSubPos = prevLen;
    
    for (i = 1; i < *numInstructions; i ++) {
        
        doInstruction = 1;
        tempNumSubs = 0;
        
        // Read the previous chars
        prevRefChar = (*oldIns)[i-1].refChar;
        prevTargetChar = (*oldIns)[i-1].targetChar;
        
        //move to the current old instruction
        currentPos = (*oldIns)[i].pos;
        currentLen = (*oldIns)[i].length;
        
        //future pos and len
        futurePos = (i == *numInstructions -1) ? currentPos + currentLen + MIN_LEN_GAP + 10 : (*oldIns)[i+1].pos;
        
        
        // Check if there is a substitution at the end of the prev instruction, and update the length of the new instruction
        
        if ( (prevPos + prevLen + 1) == currentPos) { //only a gap of one between instructions
            
            // we concatenate instructions
            currNewInsLength = currNewInsLength + prevLen + 1;
            
            // Create the new substitution
            if (prevTargetChar != prevRefChar)
                *numSubs = create_substitution(substitutions, currentSubPos, prevRefChar, prevTargetChar,0);
            
            doInstruction = 0;
        }
        
        // Type 2 substitution
        else if ( (abs(currentPos - (prevPos + prevLen)) > MIN_LEN_GAP) && (currentLen < MAX_LEN_SUB) && (prevPos + prevLen + currentLen < S_length))
            // else if ( (abs(currentPos - (prevPos + prevLen)) > MIN_LEN_GAP) && abs(futurePos - (currentPos + currentLen)) > MIN_LEN_GAP && (currentLen < MAX_LEN_SUB) && (prevPos + prevLen + currentLen < S_length))
        {
            // Count the number of tempSubs that we would add
            tempNumSubs = 0;
            for (ctr = 0; ctr<currentLen; ctr++)
                if (S[currentPos - 1 + ctr] != S[prevPos + prevLen + ctr]) tempNumSubs++;
            
            // tempNumSubs <  floor(percentageAlwaysSubs*(double)currentLen)
            // (tempNumSubs< 1 + floor(percentageSecondSubs*(double)currentLen))
            if ( (tempNumSubs < MIN_SNP) || (currentLen < lenAlwaysSubs && (abs(currentPos - (prevPos + prevLen)) > MIN_GAP) )){
                
                // we concatenate instructions
                currNewInsLength = currNewInsLength + prevLen + 1;
                
                // Create a substitution with the char of the previous instructions
                if (prevTargetChar != prevRefChar)
                    *numSubs = create_substitution(substitutions, currentSubPos, prevRefChar, prevTargetChar,0);
                
                // Create as many substitutions as length of current instruction (except in those positions where target and ref coincide)
                for (ctr = 0; ctr<currentLen; ctr++){
                    if (S[currentPos - 1 + ctr] != S[prevPos + prevLen + ctr]){
                        *numSubs = create_substitution(substitutions, currentSubPos + 1 + ctr, S[prevPos + prevLen + ctr], S[currentPos - 1 + ctr],0);
                    }
                }
                
                // Now currentPos needs to be modified, because it should now point to prevPos + prevLen + 1
                currentPos = prevPos + prevLen + 1;
                // Ref char also needs to be modified
                (*oldIns)[i].refChar = S[prevPos + prevLen + ctr];
                doInstruction = 0;
            }
        }
        
        if (doInstruction == 1){
        	// There is no substitution
            
            // Current instruction cannot be added to the new-instruction concatenation, so we finish the ongoing new instruction and start a new one
            if (prevTargetChar != prevRefChar){
                // Store the new instruction
                newInstructions[numNewInstructions].pos = posNewInstruction;
                newInstructions[numNewInstructions].length = currNewInsLength + prevLen;
                newInstructions[numNewInstructions].targetChar = prevTargetChar;
                newInstructions[numNewInstructions].refChar = prevRefChar;
                
                numNewInstructions++;
                
                // Start a new instruction
                posNewInstruction = currentPos;
                currNewInsLength = 0;
            }
            
            else{
                // We extend the ongoing instruction using the next instructions until a mismatch is found
                while (S[prevPos + prevLen] == S[currentPos - 1]){
                    //keep extending the ongoing instruction and reducing the next one
                    prevLen++;
                    currentPos++;
                    currentLen--;
                    ctr++;
                    if (currentLen == 0) break;
                }
                prevRefChar = S[prevPos + prevLen];
                if (currentLen != 0){ // We still have some letters remainding in the new instruction
                    prevTargetChar = S[currentPos - 1];
                    
                    prevLen++;
                    currentPos++;
                    currentLen--;
                    
                    // Store the new instruction
                    newInstructions[numNewInstructions].pos = posNewInstruction;
                    newInstructions[numNewInstructions].length = currNewInsLength + prevLen;
                    newInstructions[numNewInstructions].targetChar = prevTargetChar;
                    newInstructions[numNewInstructions].refChar = prevRefChar;
                    
                    numNewInstructions++;
                    
                    // Start a new instruction
                    posNewInstruction = currentPos;
                    currNewInsLength = 0;
                    
                }
                else
                {// we have merge completely the two instructions
                    prevLen++;
                    currentPos = prevPos;
                    currentLen = prevLen;
                    (*oldIns)[i].refChar = S[currentPos + currentLen - 1];
                }
                
            }
            
        }
        
        currentSubPos = currentSubPos + (*oldIns)[i].length + 1;
        prevPos = currentPos;
        prevLen = currentLen;
    }
    
    // Finish the entry of the previous new instruction
    
    prevRefChar = 'N';
    prevTargetChar = 'N';
    
    newInstructions[numNewInstructions].pos = posNewInstruction;
    newInstructions[numNewInstructions].length = currNewInsLength + prevLen;
    newInstructions[numNewInstructions].targetChar = prevTargetChar;
    newInstructions[numNewInstructions].refChar = prevRefChar;
    numNewInstructions++;
    
    /*if (prevDelete == 1){
     (*numSubs)++;
     newSub = create_substitution(currentSubPos + 1, prevTargetChar, prevRefChar);
     oldSub->next = newSub;
     oldSub = newSub;
     
     }*/
    
    *numInstructions = numNewInstructions;
    
    // Free the old Instructions and free the unsued array of the new allocations
    free(*oldIns);
    
    //newInstructions = (instruction *) realloc(newInstructions, numNewInstructions*sizeof(instruction));
    //substitutions = (substitution*)realloc(substitutions, subCtr*sizeof(substitution));
    
    // Replace the old instructions with the new ones
    *oldIns = newInstructions;
    
    
    return substitutions;
}

//**************************************************************//
//                                                              //
//          FUNCTION FOR GENERATING THE INSERTIONS              //
//                                                              //
//**************************************************************//

insertion generate_insertions(instruction* oldIns, uint32_t* numInstructions, uint32_t *numInsertions){
    
    uint32_t numNewInstructions, prevPos, prevLen, currNewInsLength, i, currentInsertionPos, currentPos, currentLen, posNewInstruction, inseCtr;
    
    char prevRefChar, prevTargetChar;
    
    // Allocate the Array for the insertions
    insertion insertions;
    insertions = (insertion ) malloc((*numInstructions)*sizeof(struct insertion_t));
    
    // Allocate the Array for the new instructions
    instruction newInstructions;
    newInstructions = (instruction ) malloc((*numInstructions)*sizeof(struct instruction_t));
    
    numNewInstructions = 0;
    inseCtr = 0;
    
    // Read first old instruction
    prevPos = (*oldIns)[0].pos;
    prevLen = (*oldIns)[0].length;
    
    //Start of the first new instruction
    posNewInstruction = prevPos;
    currNewInsLength = 0;
    
    //Start with a new insertion
    currentInsertionPos = prevLen;
    
    for (i = 1; i < *numInstructions; i ++) {
        
        // get the previous chars
        prevRefChar = (*oldIns)[i-1].refChar;
        prevTargetChar = (*oldIns)[i-1].targetChar;
        
        //move to the next old instruction
        currentPos = (*oldIns)[i].pos;
        currentLen = (*oldIns)[i].length;
        
        // Check if there is a insertion at the end of the prev instruction, and update the length of the new instruction
        if ( (prevPos + prevLen) == currentPos) { //no gap between instructions
            
            // we concatenate instructions
            currNewInsLength = currNewInsLength + prevLen;
            
            // Create the new insertion
            insertions[inseCtr].pos = currentInsertionPos + 1;
            insertions[inseCtr].targetChar = prevTargetChar;
            insertions[inseCtr].refChar = prevRefChar;
            inseCtr++;
        }
        else{// There is no substitution (no length-one gap between instructions)
            
            // Current instruction cannot be added to the new-instruction concatenation, so we finish the ongoing new instruction and start a new one
            
            // Store the new instruction
            newInstructions[numNewInstructions].pos = posNewInstruction;
            newInstructions[numNewInstructions].length = currNewInsLength + prevLen;
            newInstructions[numNewInstructions].targetChar = prevTargetChar;
            newInstructions[numNewInstructions].refChar = prevRefChar;
            numNewInstructions++;
            
            // Start a new instruction
            posNewInstruction = currentPos;
            currNewInsLength = 0;
        }
        
        currentInsertionPos = currentInsertionPos + currentLen + 1;
        prevPos = currentPos;
        prevLen = currentLen;
    }
    
    
    
    // Finish the entry of the previous new instruction
    prevRefChar = (*oldIns)[i-1].refChar;
    prevTargetChar = (*oldIns)[i-1].targetChar;
    
    
    newInstructions[numNewInstructions].pos = posNewInstruction;
    newInstructions[numNewInstructions].length = currNewInsLength + prevLen;
    newInstructions[numNewInstructions].targetChar = prevTargetChar;
    newInstructions[numNewInstructions].refChar = prevRefChar;
    numNewInstructions++;
    
    /*if (prevDelete == 1){
     (*numSubs)++;
     newSub = create_substitution(currentSubPos + 1, prevTargetChar, prevRefChar);
     oldSub->next = newSub;
     oldSub = newSub;
     
     }*/
    
    *numInstructions = numNewInstructions;
    *numInsertions = inseCtr;
    
    // Free the old Instructions and free unused memory of insertions
    free(*oldIns);
    newInstructions = (instruction) realloc(newInstructions, numNewInstructions*sizeof(struct instruction_t));
    insertions = (insertion) realloc(insertions, inseCtr*sizeof(struct insertion_t));
    
    // Replace the old instructions with the new ones
    *oldIns = newInstructions;
    
    
    return insertions;
}


//**************************************************************//
//                                                              //
//          FUNCTION FOR GENERATING CHROMOSOME MAPPING          //
//                                                              //
//**************************************************************//

chromosome generate_chromosome_mapping(uint8_t chrCtr, char *S, char *T, uint32_t tagetLength, uint32_t refLength, int * SA){
    
    chromosome chr;
    
    chr = create_chromosome(NULL, NULL, NULL, 0, 0, 0);
    
    // Do mapping between ref and target
    chr->instructions = generate_instructions(S, T, tagetLength, refLength, SA, &(chr->numInst));
    
    // Do Substitutions
    chr->substitutions = generate_substitutions(&(chr->instructions), &(chr->numInst), &(chr->numSubs), S,refLength);
    // printf("Substitutions created\n");
    
    chr->insertions = generate_insertions(&(chr->instructions), &(chr->numInst), &(chr->numInse));
    // printf("Insertions created\n");
    
    return chr;
    
}




