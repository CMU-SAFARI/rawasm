#include "rawasm.h"

///////////////////////////////
//  Rawasm Main Extension    //
///////////////////////////////

// The function assumes to have a list of filenames. It can work with only 1 file if specified so.
int ra_ug_seq(ma_ug_t *g, const sdict_t *d, const ma_sub_t *sub, char **src_fname, const int src_fnum, const char * dst_dir)
{
	// TODO LIST:
	//	--	Make the patch (self-contained version (if possible) in order to make it run on the SAFARI server)
	//	--	Complete copy of the signal structure (including metadata) (atm we only copy the float signal)

	// NOTES:
    //  -- Reverse complement is an open problem atm
	//	-- We can process multiple files or a monolithic file. However using a mobolithic file is much faster. (overhead for close/open I guess)
	//	-- Tested on 16 input files for Fast5, Pod5, Slow5.

	////////////////////////////////////////////////////////////////////////////////////////

	// Raw signals variables
	ri_sig_file_t* ug_sig_fr;		// Unitigs signal file reader pointer
	ri_sig_file_t** ug_sig_fw;		// Unitigs signal file writer pointer
	size_t read_count;				// Number of reads from file reader
	size_t ug_count;
	size_t ug_radix;
	size_t ug_format;
	size_t ug_namelen;
	char * utg_fname;
	char utg_nbuff[64];
	ri_sig_t* tmp_sig;				// Temporary signal
	ri_sig_t* trimmed_sig;			// Temporary signal (trimmed)
	ra_ugl_t* ugl;

	// Dynamic allocation for local variables
	utg_intv_t * tmp = (utg_intv_t*)calloc(d->n_seq, sizeof(utg_intv_t));		// Unitig Vector
	ug_sig_fw = (ri_sig_file_t**)calloc(g->u.n, sizeof(ri_sig_file_t *));		// File Writer Vector
	tmp_sig = (ri_sig_t *) calloc(1, sizeof(ri_sig_t));							// Temporary Signal Structure
	trimmed_sig = (ri_sig_t *) calloc(1, sizeof(ri_sig_t));						// Trimmed Temporary Signal Structure

	// Building the Files (FAST5, POD5, SLOW5)	
	ug_count = g->u.n;
	ug_radix = 8;	// Utg num digits (fixed to 8, max 10^8 unitigs)
	ug_format = (strstr(src_fname[0], ".pod5")) ? strlen(".pod5") : strlen(".fast5");
	ug_namelen = strlen(dst_dir) + strlen("/utg") + ug_radix + ug_format;
	utg_fname = (char *) calloc(ug_namelen + 1, sizeof(char));

	memset(utg_fname, '0', ug_namelen);
	memcpy(utg_fname, dst_dir, strlen(dst_dir));
	memcpy(utg_fname + strlen(dst_dir), "/utg", strlen("/utg"));

	if (strstr(src_fname[0], ".fast5")) 	memcpy(utg_fname + strlen(dst_dir) + strlen("/utg") + ug_radix, ".fast5\0", strlen(".fast5\0"));
	else if (strstr(src_fname[0], ".pod5")) memcpy(utg_fname + strlen(dst_dir) + strlen("/utg") + ug_radix, ".pod5\0", strlen(".pod5\0"));
	else if (strstr(src_fname[0], ".slow5")) memcpy(utg_fname + strlen(dst_dir) + strlen("/utg") + ug_radix, ".slow5\0", strlen(".slow5\0"));
	else if (strstr(src_fname[0], ".blow5")) memcpy(utg_fname + strlen(dst_dir) + strlen("/utg") + ug_radix, ".blow5\0", strlen(".blow5\0"));

	fprintf(stderr, "[M::%s] Building %lu unitigs, please wait\n", __func__, ug_count);
	//fprintf(stderr, "[M::%s] Number of sequences in dictionary: %u\n", __func__, d->n_seq);

	for (int i = 0; i < ug_count; i++){
		sprintf(utg_nbuff, "%08u", i);
		memcpy(utg_fname + strlen(dst_dir) + strlen("/utg"), utg_nbuff, ug_radix);
		//fprintf(stderr, "%s\n", utg_fname);
		ug_sig_fw[i] = open_sig(utg_fname, SIG_WRITE_OP);
	}

	free(utg_fname);

	// Create the Unitig List Structure
	ugl = ra_ugl_create(ug_count);
	for (int i = 0; i < ug_count; i++){

		ma_utg_t *tmp_utg = &g->u.a[i];							// Select the unitig i
		uint32_t l = 0;
		//fprintf(stderr, "[M::%s] Unitig %u\n", __func__, i);
		ra_ugl_readlist_create(ugl, i, tmp_utg->n);

		for (int j = 0; j < tmp_utg->n; ++j) {	

			utg_intv_t *t = &tmp[tmp_utg->a[j]>>33];					// a[j] is a 64 bit value, so maybe it is the ID for the read
			assert(t->len == 0);
			t->utg = i, t->ori = tmp_utg->a[j]>>32&1;					// Assign the utg ID, strand, start and length. Then increase the global len of the utg.
			t->start = l, t->len = (uint32_t)tmp_utg->a[j];
			l += t->len;
			//fprintf(stderr, "\t\tRead %08x - index=%u : utg=%u : ori=%u : start=%u : len=%u\n", j, tmp_utg->a[j]>>33, t->utg, t->ori, t->start, t->len);
			ra_ugl_readlist_insert(ugl, i, t->start);
		}
	}
	
	// Write the Signal Files
	for(int k = 0; k < src_fnum; k++){
			
		ug_sig_fr = open_sig(src_fname[k], SIG_READ_OP);

		if(!ug_sig_fr) return 0;
		// get the total number of signals in the file
		read_count = ri_sig_count(ug_sig_fr);
		// For all the reads in this file
		for(int w = 0; w < read_count; w++){					

			utg_intv_t * ug_ptr;
			int32_t id;
			uint32_t ug_id;
			int ug_key;
			uint32_t ug_rlen;

			// Read the next signal from rawsignal file
			ri_read_sig(ug_sig_fr, tmp_sig);
			
			// Get the read ID 
			id = sd_get(d, tmp_sig->name);
			// If that read is valid go on (i.e. it is mapped on at least a unitig)
			if (id < 0 || tmp[id].len == 0) { free(tmp_sig->sig); free(tmp_sig->name); continue;	}	
			// Find the unitig that contains the read

			ug_ptr = &tmp[id];
			ug_id = ug_ptr->utg;
			ug_key = ug_ptr->start;
			ug_rlen = ug_ptr->len;	// This could be misleading, but it is the read length in the unitig
			// fprintf(stderr, "Signal found for utg %u in file %u, number %u :: ", ug_id, k, w);

			// Build the signal ( Reverse complement is not supported atm )
			if (sub) {									
				assert(sub[id].e - sub[id].s <= tmp_sig->l_sig);
				memmove(tmp_sig->sig, tmp_sig->sig + sub[id].s, sub[id].e - sub[id].s);
				tmp_sig->l_sig = sub[id].e - sub[id].s;
			}	

			// Here generate the trimmed signal
			ra_sig_trim_copy(trimmed_sig, tmp_sig, ug_rlen);
			
			// Write the signal into the correct position on destination file
			ra_ugl_readlist_add(ugl, ug_id, ug_key, trimmed_sig);
			ra_ug_readlist_commit(ugl, ug_id, ug_sig_fw[ug_id]);

			// Free the signal-related pointers after use (Signals sizes are supposed to change from one signal to another)
			free(tmp_sig->sig);		free(tmp_sig->name);
			free(trimmed_sig->sig);	free(trimmed_sig->name);

		}

		ri_sig_close(ug_sig_fr);

	}

	ra_ugl_check(ugl);
	// Free section
	for (int i = 0; i < ug_count; i++) ri_sig_close(ug_sig_fw[i]);
	free(ug_sig_fw);
	free(tmp);
	free(tmp_sig);
	free(trimmed_sig);
	ra_ugl_destroy(&ugl);


	return 0;
}

void ra_ug_print(const ma_ug_t *ug, const sdict_t *d, const ma_sub_t *sub, FILE *fp)
{
	uint32_t i, j, l;
	char name[128];

	for (i = 0; i < ug->u.n; ++i) { // the Segment lines in GFA

		ma_utg_t *p = &ug->u.a[i];
		sprintf(name, "utg%.6d%c", i + 1, "lc"[p->circ]);
		fprintf(fp, "S\t%s\tUGID-%08d\tLN:i:%d\n", name, i, p->len);

		if (p->circ) {  // make circularising links (both forward and reverse directions) for circular unitigs
			fprintf(fp, "L\t%s\t+\t%s\t+\t0M\n", name, name);
			fprintf(fp, "L\t%s\t-\t%s\t-\t0M\n", name, name);
		}
		for (j = l = 0; j < p->n; l += (uint32_t)p->a[j++]) {
			uint32_t x = p->a[j]>>33;
			if (sub) fprintf(fp, "a\t%s\t%d\t%s:%d-%d\t%c\t%d\n", name, l, d->seq[x].name, sub[x].s + 1, sub[x].e, "+-"[p->a[j]>>32&1], (uint32_t)p->a[j]);
			else fprintf(fp, "a\t%s\t%d\t%s\t%c\t%d\n", name, l, d->seq[x].name, "+-"[p->a[j]>>32&1], (uint32_t)p->a[j]);
		}

	}

	for (i = 0; i < ug->g->n_arc; ++i) { // the Link lines in GFA
		uint32_t u = ug->g->arc[i].ul>>32, v = ug->g->arc[i].v;
		fprintf(fp, "L\tutg%.6d%c\t%c\tutg%.6d%c\t%c\t%dM\tSD:i:%d\n", (u>>1)+1, "lc"[ug->u.a[u>>1].circ], "+-"[u&1],
				(v>>1)+1, "lc"[ug->u.a[v>>1].circ], "+-"[v&1], ug->g->arc[i].ol, asg_arc_len(ug->g->arc[i]));
	}
	for (i = 0; i < ug->u.n; ++i) { // summary of unitigs
		uint32_t cnt[2];
		ma_utg_t *u = &ug->u.a[i];
		if (u->start == UINT32_MAX) {
			fprintf(fp, "x\tutg%.6dc\t%d\t%d\n", i + 1, u->len, u->n);
		} else {
			for (j = 0; j < 2; ++j) cnt[j] = asg_arc_n(ug->g, i<<1|j);
			if (sub)
				fprintf(fp, "x\tutg%.6dl\t%d\t%d\t%d\t%d\t%s:%d-%d\t%c\t%s:%d-%d\t%c\n", i + 1, u->len, u->n, cnt[1], cnt[0],
						d->seq[u->start>>1].name, sub[u->start>>1].s + 1, sub[u->start>>1].e, "+-"[u->start&1],
						d->seq[u->end>>1].name, sub[u->end>>1].s + 1, sub[u->end>>1].e, "+-"[u->end&1]);
			else
				fprintf(fp, "x\tutg%.6dl\t%d\t%d\t%d\t%d\t%s\t%c\t%s\t%c\n", i + 1, u->len, u->n, cnt[1], cnt[0],
						d->seq[u->start>>1].name, "+-"[u->start&1], d->seq[u->end>>1].name, "+-"[u->end&1]);
		}
	}

}

//////////////////////////////////
// Rawasm Unitig List Methods	//
//////////////////////////////////

ra_ugl_t * ra_ugl_create(const size_t ug_count){

	ra_ugl_t * ugl = (ra_ugl_t *)calloc(1, sizeof(ra_ugl_t));
	ugl->ug = (ra_ugl_ug_t *)calloc(ug_count, sizeof(ra_ugl_ug_t));
	ugl->n = ug_count;
	return ugl;
}

void ra_ugl_readlist_create(ra_ugl_t * ugl, const uint32_t ug_id, const uint32_t signum){

	ugl->ug[ug_id].read = (ra_ugl_ugr_t *)calloc(signum, sizeof(ra_ugl_ugr_t));
	ugl->ug[ug_id].n = signum;
	ugl->ug[ug_id].last_commit = -1;

	for(int i = 0; i < ugl->ug[ug_id].n; i++){
		ugl->ug[ug_id].read[i].commit = 0;
		ugl->ug[ug_id].read[i].key = -1;
	}
}

void ra_ugl_readlist_insert(ra_ugl_t * ugl, const uint32_t ug_id, const uint64_t key){

	// Inserting an empty signal slot ordered by key value
	for(int i = 0; i < ugl->ug[ug_id].n; i++){
		// If we reach an empty value, then init that value to key.
		if(ugl->ug[ug_id].read[i].key == -1){
			ugl->ug[ug_id].read[i].key = key;
			return;
		}

		if(ugl->ug[ug_id].read[i].key > key){
			for(int j = ugl->ug[ug_id].n-1; j > i; j--) ugl->ug[ug_id].read[j].key = ugl->ug[ug_id].read[j-1].key;
			ugl->ug[ug_id].read[i].key = key;
			return;
		}
	}

	
}

void ra_ugl_readlist_add(ra_ugl_t * ugl, const uint32_t ug_id, const int key, ri_sig_t * sig){

	ra_ugl_ug_t * ug = &ugl->ug[ug_id];
	uint32_t low, high, val;

	// Let's use a Binary Search
	low = 0;
	high = ug->n - 1;

	//fprintf(stderr, "Key %d, size is %u, let's test\n", key, high + 1);

	while(low <= high){
		val = low + (high - low) / 2;
		//fprintf(stderr, "Val:%u, high:%u, low:%u -- %d\n", val, high, low, ug->read[val].key);
		if(ug->read[val].key == key){

			// Copy here the signal
			ug->read[val].sig = (ri_sig_t *) calloc(1, sizeof(ri_sig_t));
			ug->read[val].sig->l_sig = sig->l_sig;
			ug->read[val].sig->offset = sig->offset;
			ug->read[val].sig->rid = sig->rid;
			ug->read[val].sig->name = (char *)calloc(37, sizeof(char));
			ug->read[val].sig->sig = (float *)calloc(sig->l_sig, sizeof(float));
			memcpy(ug->read[val].sig->name, sig->name, 37);
			memcpy(ug->read[val].sig->sig, sig->sig, sig->l_sig);
			return;

		} else if(ug->read[val].key < key) {
			low = val + 1;
		} else {
			high = val - 1;
		}
	}

}

void ra_ug_readlist_commit(ra_ugl_t * ugl, const uint32_t ug_id, ri_sig_file_t* sig_fw){

	ra_ugl_ug_t * ug = &ugl->ug[ug_id];

	int index_commit = ug->last_commit + 1;

	// Check if the first signal to commit is available
	if(ug->read[index_commit].sig == NULL){
		//fprintf(stderr, "\tCannot Commit<%08u>::%d :: ", ug_id, index_commit);
		return;
	}

	while(ug->read[index_commit].sig != NULL){

		//fprintf(stderr, "\tCommit Unitig<%08u>::%d :: ", ug_id, index_commit);
		// Commit the signal on file
		ri_write_sig(sig_fw, ug->read[index_commit].sig, 1, 0);
		ug->read[index_commit].commit = 1;
		ug->last_commit++;

		// Once the signal is written on file, release its memory.
		free(ug->read[index_commit].sig->sig);
		free(ug->read[index_commit].sig->name);
		free(ug->read[index_commit].sig);

		// Increase the counter
		index_commit++;

		if(index_commit == ug->n) break;

	}
	
}

void ra_ugl_destroy(ra_ugl_t ** ugl){

	// Signals are freed in the commit function
	ra_ugl_t * tmp = *ugl;
	for(int i = 0; i < tmp->n; i++) free(tmp->ug[i].read);
	free(tmp->ug);
	free(*ugl);
}

// Do not use this function unless commit is not used
void ra_ugl_fprint(ra_ugl_t * ugl){
	for(int i = 0; i < ugl->n; i++){
		fprintf(stderr, "Unitig<%08d>:\n", i);
		for(int j = 0; j < ugl->ug[i].n; j++){
			fprintf(stderr, "\t\tRead<%08d>: [%u, %d]\n\t\t\t{", j, ugl->ug[i].read[j].commit, ugl->ug[i].read[j].key);
			int m = (ugl->ug[i].read[j].sig->l_sig < 16) ? ugl->ug[i].read[j].sig->l_sig : 16;
			for(int k = 0; k < m; k++){
				fprintf(stderr,"%02.4f ", ugl->ug[i].read[j].sig->sig[k]);
			}
			if(m == 16) fprintf(stderr, "... }\n");
			else fprintf(stderr, "}\n");
		}
	}
}

//////////////////
//  Rawasm Misc //
//////////////////

// Given the filename, check if it is a single file or a directory (containing multiple files of the same type)
char ** ra_get_fnum(const char * fn_reads, int * filenum){

    fprintf(stderr, "[M::%s] Counting HD5 files ... \n", __func__);
	ri_sig_file_t** fr;		// Unitigs signal file reader
    char ** fn_list;

    // First, check if the path is a Single File or a Directory
    if(strstr(fn_reads, ".fast5") || strstr(fn_reads, ".pod5") || strstr(fn_reads, ".slow5") || strstr(fn_reads, ".blow5")){
        // One single file, copy the path into the fn_list
        *filenum = 1;
        fn_list = (char **) calloc(*filenum, sizeof(char *));
        fn_list[0] = (char *) calloc(strlen(fn_reads), sizeof(char));
        strcpy(fn_list[0], fn_reads);
    } else {
        // Open the directory with signals (it MUST only contains rawsignal files)
        struct dirent* in_file;
        DIR * fd = opendir(fn_reads);

        if(!fd){
            fprintf(stderr, "\x1B[31m[M::%s] CRITICAL ERROR - DIRECTORY NOT FOUND, INVALID OR CORRUPTED\x1B[0m\n", __func__);
            exit(-1);
        }
        
        // Count all the signals files
        *filenum = 0;
        while ((in_file = readdir(fd))) *filenum = (!strcmp(in_file->d_name,".") || !strcmp(in_file->d_name,"..")) ? *filenum : *filenum + 1;
        rewinddir(fd);

        // Save all files in a char list
        fn_list = (char **) calloc(*filenum, sizeof(char *));
        int local_cnt = 0;
        while ((in_file = readdir(fd))){
            if(strcmp(in_file->d_name,".") && strcmp(in_file->d_name,"..")){
                fn_list[local_cnt] = (char *) calloc(strlen(fn_reads) + 1 + strlen(in_file->d_name), sizeof(char));
                strcpy(fn_list[local_cnt], fn_reads);
                strcat(fn_list[local_cnt], "/");
                strcat(fn_list[local_cnt], in_file->d_name);
                local_cnt++;
                
            }
        }
    }
	
	fprintf(stderr, "[M::%s] Ready to process %d files \n", __func__, *filenum);

    return fn_list;

}

void ra_sig_trim_copy(ri_sig_t* trimmed_sig, ri_sig_t* tmp_sig, const uint32_t sig_len){

	trimmed_sig->name = (char *)calloc(37,sizeof(char));
	memcpy(trimmed_sig->name, tmp_sig->name, 37);
	trimmed_sig->sig = (float *)calloc(sig_len, sizeof(float));
	trimmed_sig->l_sig = sig_len;
	
	for (int f = 0; f < trimmed_sig->l_sig; ++f)
		trimmed_sig->sig[f] = tmp_sig->sig[f];
}

void ra_ugl_check(ra_ugl_t * ugl){

	fprintf(stderr, "[M::%s] Unitig files commit notification\n", __func__);

	for(int i = 0; i < ugl->n; i++){
		fprintf(stderr, "\t\tUnitig<%08d> Committed %d/%d\n", i, ugl->ug[i].last_commit + 1, ugl->ug[i].n);
	}

}