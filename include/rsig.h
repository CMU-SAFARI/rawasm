#ifndef RSIG_H
#define RSIG_H

#include "rutils.h"
#ifndef NHDF5RH
#include "hdf5_tools.hpp"
#endif
#ifndef NPOD5RH
#include "pod5_format/c_api.h"
#endif
#ifndef NSLOW5RH
#include <slow5/slow5.h>
#endif

#include "uuid.h"


#ifdef __cplusplus
extern "C" {
#endif

// Allow to open files in both write and read mode
// Implement the write file function

#define SIG_READ_OP		0
#define SIG_WRITE_OP 	1

#ifndef NPOD5RH

// POD5-related data types
struct pod5_file_t{
	Pod5FileReader_t* reader; 
	Pod5FileWriter_t* writer;
};


#endif

typedef struct ri_sig_s{
	uint32_t rid, l_sig; //read id and length of the signal values
	char *name; //name of the read

	uint64_t offset; // offset in ri_idx_t::S

	float dig, ran; //digitalisation, range
	float* sig; //signal values of a read
} ri_sig_t;

typedef struct { size_t n, m; ri_sig_t **a; } rhsig_v;

typedef struct ri_sig_file_s {
	// gzFile fp;
	// kseq_t *ks;
	// ri_sig_t s;
	int num_read; //Number of reads
	int cur_read; //Number of processed reads by RawHash (shows the id of the next read to process)
	
	//HDF5-related
	char** raw_path; //List of paths to raw values
	char** ch_path; //List of paths to channel
	#ifndef NHDF5RH
	hdf5_tools::File* fp; //FAST5 file pointer
	#endif
	
	//POD5-related
	unsigned long int pod5_batch_count;
	unsigned long int pod5_batch;
	unsigned long int pod5_row_count;
	unsigned long int pod5_row;
	#ifndef NPOD5RH
	Pod5ReadRecordBatch_t* batch;
	pod5_file_t * pp;	//POD5 file pointer
	#endif

	#ifndef NSLOW5RH
	slow5_file_t* sp; //SLOW5 file pointer
	#endif
} ri_sig_file_t;

/**
 * Closes the signal file.
 *
 * @param fp	a struct that includes the file pointer to a signal file
 * 				Allocated memories in the struct are destroyed here.
 */
void ri_sig_close(ri_sig_file_t *fp);

/**
 * Opens the signal file (e.g., a FAST5 file)
 *
 * @param fn	path to the signal file
 * @param mode	file operation mode (0 for Read or 1 for Write)
 * 				POD5
 * 					-	opening in write mode means creating a new file. Opening in write mode an existing pod5 file is not supported atm.
 * 						only standard write opening options are supported. 
 * 
 * @return		a struct that includes the file pointer to the opened signal file
 * 				Returned struct (and its variables) is allocated in this function.
 */
ri_sig_file_t *open_sig(const char *fn, uint8_t const mode);

/**
 * Opens all the signal files (e.g., FAST5 files)
 *
 * @param n		number of files
 * @param fn	list of paths to the files
 * @param mode	files operation mode (0 for Read or 1 for Write)
 * 
 * @return		List of structs that include the file pointers to each opened signal file
 * 				Returned structs (and their variables) are allocated in this function.
 */
ri_sig_file_t **open_sigs(int n, const char **fn, const uint8_t mode);

/**
 * Converts the sequence into its expected event values
 *
 * @param str		sequence to convert to its expected event values
 * @param len		length of the $str
 * @param pore_vals	expected event values for each possible k-mer
 * @param pore_kmer	k-mer size of each event
 * @param strand	directin of strand. 1 is forward, 0 is reverse direction.
 * @param s_len		length of $s_values
 * @param s_values	expected event values of each k-mer in $str
 */
void ri_seq_to_sig(const char *str, int len, const float* pore_vals, const int pore_kmer, const int strand, uint32_t* s_len, float* s_values);

/**
 * Write a new signal values in the file
 *
 * @param fp					file pointer to the signal file (i.e., either FAST5 or SLOW5)
 * @param s						attribute of the read and the signal values.
 * @param calibration_scale		used to turn the signal floating point representation to int16
 * @param calibration_offset	used to turn the signal floating point representation to int16
 * 
 */

void ri_write_sig(ri_sig_file_t* fp, ri_sig_t* s, const float calibration_scale, const float calibration_offset);

/**
 * Recursively find all files that ends with "fast5" under input directory const char *A
 *
 * @param A			path to a directory where the fast5 files are searched
 * @param fnames	list of fast5 or pod5 files
 * 					fnames->a = List of file names
 * 					fnames->n = Number of fast5 files
 */

void ri_read_sig(ri_sig_file_t* fp, ri_sig_t* s);

/**
 * Recursively find all files that ends with "fast5" under input directory const char *A
 *
 * @param A			path to a directory where the fast5 files are searched
 * @param fnames	list of fast5 or pod5 files
 * 					fnames->a = List of file names
 * 					fnames->n = Number of fast5 files
 */

void find_sfiles(const char *A, ri_char_v *fnames);

/**
 * Get the total number of signal in the file
 *
 * @param 
 * @param 
 * 
 */

// Counting the number of signals (doesn't affect the current record)
int ri_sig_count(ri_sig_file_t* hd5_fp);

#ifdef __cplusplus
}
#endif
#endif //RSIG_H

