
#ifndef RAWASM_H
#define RAWASM_H

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include "sdict.h"
#include "asg.h"
#include "miniasm.h"
#include "rsig.h"

//////////////////////////
//  Rawasm Structures   //
//////////////////////////

typedef struct {
	uint8_t commit;
    int key;
    ri_sig_t * sig;
} ra_ugl_ugr_t; // Rawasm Unitig List Read Signal

typedef struct {
	ra_ugl_ugr_t * read;
    int last_commit; // Index of the last committed signal write on file
    uint32_t n;
} ra_ugl_ug_t; // Rawasm Unitig List Read 

typedef struct {
	ra_ugl_ug_t * ug;
    uint32_t n;
} ra_ugl_t; // Rawasm Unitig List

// Redefine some local structures to asm.c
typedef struct {
	uint32_t utg:31, ori:1, start, len;
} utg_intv_t;


///////////////////////////////
//  Rawasm Main Extension    //
///////////////////////////////

void ra_ug_print(const ma_ug_t *ug, const sdict_t *d, const ma_sub_t *sub, FILE *fp);
int ra_ug_seq(ma_ug_t *g, const sdict_t *d, const ma_sub_t *sub, char **src_fname, const int src_fnum, const char * dst_dir);

//////////////////////////////////
//  Rawasm Unitig List Methods	//
//////////////////////////////////

ra_ugl_t * ra_ugl_create(const size_t ug_count);
void ra_ugl_readlist_create(ra_ugl_t * ugl, const uint32_t ug_id, const uint32_t signum);
void ra_ugl_readlist_insert(ra_ugl_t * ugl, const uint32_t ug_id, const uint64_t key);
void ra_ugl_readlist_add(ra_ugl_t * ugl, const uint32_t ug_id, const int key, ri_sig_t * sig);
void ra_ug_readlist_commit(ra_ugl_t * ugl, const uint32_t ug_id, ri_sig_file_t* sig_fw);
void ra_ugl_destroy(ra_ugl_t **);
void ra_ugl_fprint(ra_ugl_t * ugl);

//////////////////
//  Rawasm Misc //
//////////////////

char ** ra_get_fnum(const char * fn_reads, int * filenum);
void ra_sig_trim_copy(ri_sig_t* trimmed_sig, ri_sig_t* tmp_sig, const uint32_t sig_len);
void ra_ugl_check(ra_ugl_t * ugl);

#endif