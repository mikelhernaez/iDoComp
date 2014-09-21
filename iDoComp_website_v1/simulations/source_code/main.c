//
//  main.c
//  iDoComp_v1
//
//  Created by Mikel Hernaez on 8/7/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "iDoComp.h"

BASEPAIR char2BP(char c){
    switch (c) {
        case 'A':
            return bp_A;
        case 'C':
            return bp_C;
        case 'G':
            return bp_G;
        case 'T':
            return bp_T;
        case 'N':
            return bp_N;
        case 'U':
            return bp_U;
        case 'R':
            return bp_R;
        case 'Y':
            return bp_Y;
        case 'K':
            return bp_K;
        case 'M':
            return bp_M;
        case 'S':
            return bp_S;
        case 'W':
            return bp_W;
        case 'B':
            return bp_B;
        case 'D':
            return bp_D;
        case 'H':
            return bp_H;
        case 'V':
            return bp_V;
        case 'X':
            return bp_X;
        case 'a':
            return bp_a;
        case 'c':
            return bp_c;
        case 'g':
            return bp_g;
        case 't':
            return bp_t;
        case 'n':
            return bp_n;
        case 'u':
            return bp_u;
        case 'r':
            return bp_r;
        case 'y':
            return bp_y;
        case 'k':
            return bp_k;
        case 'm':
            return bp_m;
        case 's':
            return bp_s;
        case 'w':
            return bp_w;
        case 'b':
            return bp_b;
        case 'd':
            return bp_d;
        case 'h':
            return bp_h;
        case 'v':
            return bp_v;
        case 'x':
            return bp_x;
        case '#':
            return bp_N;
        default:
            printf("%c is not a valid FASTA base pair in a file, using N instead.\n", c);
            return bp_N;
    }
};

BASEPAIR BP2char(BASEPAIR c){
    switch (c) {
        case bp_A:
            return 'A';
        case bp_C:
            return 'C';
        case bp_G:
            return 'G';
        case bp_T:
            return 'T';
        case bp_N:
            return 'N';
        case bp_U:
            return 'U';
        case bp_R:
            return 'R';
        case bp_Y:
            return 'Y';
        case bp_K:
            return 'K';
        case bp_M:
            return 'M';
        case bp_S:
            return 'S';
        case bp_W:
            return 'W';
        case bp_B:
            return 'B';
        case bp_D:
            return 'D';
        case bp_H:
            return 'H';
        case bp_V:
            return 'V';
        case bp_X:
            return 'X';
        case bp_a:
            return 'a';
        case bp_c:
            return 'c';
        case bp_g:
            return 'g';
        case bp_t:
            return 't';
        case bp_n:
            return 'n';
        case bp_u:
            return 'u';
        case bp_r:
            return 'r';
        case bp_y:
            return 'y';
        case bp_k:
            return 'k';
        case bp_m:
            return 'm';
        case bp_s:
            return 's';
        case bp_w:
            return 'w';
        case bp_b:
            return 'b';
        case bp_d:
            return 'd';
        case bp_h:
            return 'h';
        case bp_v:
            return 'v';
        case bp_x:
            return 'x';
        default:
            printf("Something wrong compressing chars...\n");
            return 0;
    }
};

int compression_main(int argc, const char * argv[])
{
    
    FILE *inputFile, *headerFile;
    
    char outputPath[256];
    
    clock_t begin, end;
    
    unsigned int chrCtr = 0, bpPerLine = 0, compressedSize = 0;
    
    chromosome chr;
    
    if (argc != 4) {
        printf("ERROR: The input argument is 1) c for compression, 2) a .txt file with the reference, the target and the SA, and 3) output (i.e., compressed) file path. \n");
        return 9;
    }
    
    begin = clock();
    
    //open the input file
    inputFile = fopen ( argv[2] , "r" );
    if (inputFile==NULL) {fputs ("Input File error\n",stderr); exit (1);}
    
    strcpy(outputPath, argv[3]);
    strcat(outputPath, "_headers.ido");
    
    //open the Header file
    headerFile = fopen ( outputPath , "w" );
    if (headerFile==NULL) {fputs ("Header File error\n",stderr); exit (1);}
    
    //open the dataFile
    strcpy(outputPath, argv[3]);
    
    chr = generate_chromosome_sequence(inputFile, headerFile, &chrCtr, &bpPerLine);
    
    compressedSize = start_fasta_compression(chr, outputPath, chrCtr, bpPerLine);
    
    end = clock();
    
    
    printf("Compression CPU Time: %.2f seconds\n", (double)(end - begin)/ CLOCKS_PER_SEC);
    
    printf("Compressed Size: %d Bytes\n", compressedSize);
    
    return 0;
}

int decompression_main(int argc, const char *argv[]){
    
    
    FILE *inputFile, *headerFile;
    
    char outputPath[256];
    
    clock_t begin, end;
    
    if (argc != 4) {
        printf("ERROR: The input argument is 1) d for decompression, 2) a .txt file with the reference, and the reconstructed file, and 3) compressed file path.\n");
        return 9;
    }
    
    begin = clock();
    
    //open the input file
    inputFile = fopen ( argv[2] , "r" );
    if (inputFile==NULL) {fputs ("Input File error\n",stderr); exit (1);}
    
    strcpy(outputPath, argv[3]);
    strcat(outputPath, "_headers.ido");
    
    //open the Header file
    headerFile = fopen ( outputPath , "r" );
    if (headerFile==NULL) {fputs ("Header File error\n",stderr); exit (1);}
    
    //open the dataFile
    strcpy(outputPath, argv[3]);
    
    start_fasta_decompression(outputPath, inputFile, headerFile);
    
    end = clock();
    
    
    printf("Decompression CPU Time: %.2f seconds\n", (double)(end - begin)/ CLOCKS_PER_SEC);
    
    return 0;
}


int main(int argc, const char * argv[])
{

    //Check that the input argument are ok
    
    switch ((int)argv[1][0]) {
        case 'c':
            compression_main(argc, argv);
            return 1;
            
        case 'd':
            decompression_main(argc, argv);
            return 1;
            
            
        default:
            printf("ERROR: use 'c' for compression or 'd' for decompression as first argument.\n");
            return 0;
    }
    
    
    return 0;
}


