//
//  fasta_stream.c
//  iDoComp_v1
//
//  Created by Mikel Hernaez on 8/7/14.
//  Copyright (c) 2014 Mikel Hernaez. All rights reserved.
//

#include "fasta_compressor.h"



arithStream initialize_arithStream_Ints(char* osPath, uint8_t decompressor_flag){
    
    arithStream as;
    FILE *fos;
    
    uint32_t osPathLength = (uint32_t)strlen(osPath);
    
    strcat(osPath, "_ints.ido");
    fos = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    as = (arithStream) calloc(1, sizeof(struct arithStream_t));
    as->stats = initialize_stream_stats_Ints();
    as->a = initialize_arithmetic_encoder(m_INTS);
    as->os = initialize_osStream(1, fos, NULL, decompressor_flag);
    as->a->t = (decompressor_flag)? read_uint32_from_stream(as->a->m, as->os):0;
    
    *(osPath + osPathLength) = 0;
    
    return as;
    
}

arithStream initialize_arithStream_Signs(char* osPath, uint8_t decompressor_flag){
    
    arithStream as;
    FILE *fos;
    
    uint32_t osPathLength = (uint32_t)strlen(osPath);
    
    strcat(osPath, "_signs.ido");
    fos = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    as = (arithStream) calloc(1, sizeof(struct arithStream_t));
    as->stats = initialize_stream_stats_Ints();
    as->a = initialize_arithmetic_encoder(m_SIGNS);
    as->os = initialize_osStream(1, fos, NULL, decompressor_flag);
    as->a->t = (decompressor_flag)? read_uint32_from_stream(as->a->m, as->os):0;
    
    *(osPath + osPathLength) = 0;
    
    return as;
    
}

arithStream initialize_arithStream_Chars(char* osPath, uint8_t decompressor_flag, uint8_t **BP_trans){
    
    arithStream as;
    FILE *fos;
    
    uint32_t osPathLength = (uint32_t)strlen(osPath);
    
    strcat(osPath, "_char.ido");
    fos = (decompressor_flag)? fopen(osPath, "r"):fopen(osPath, "w");
    
    as = (arithStream) calloc(1, sizeof(struct arithStream_t));
    as->stats = initialize_stream_stats_Chars(BP_trans);
    as->a = initialize_arithmetic_encoder(m_CHARS);
    as->os = initialize_osStream(1, fos, NULL, decompressor_flag);
    as->a->t = (decompressor_flag)? read_uint32_from_stream(as->a->m, as->os):0;
    
    *(osPath + osPathLength) = 0;
    
    return as;
    
}

fasta_compressor initialize_fasta_compressor(char osPath[], uint8_t streamDirection, uint8_t **BP_trans){
    
    fasta_compressor s;
    
    s = calloc(1, sizeof(struct fasta_compressor_t));
    
    s->Ints = initialize_arithStream_Ints(osPath, streamDirection);
    s->Chars = initialize_arithStream_Chars(osPath, streamDirection, BP_trans);
    s->Signs = initialize_arithStream_Signs(osPath, streamDirection);
    
    return s;
}
