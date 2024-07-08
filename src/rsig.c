#include "rsig.h"
#include "rh_kvec.h"
#include <math.h>
#ifndef NHDF5RH
#include <cstring>
#include <cassert>
#endif
#include <errno.h>
#include <assert.h>
#include <sys/stat.h>
#include <dirent.h>

void ri_seq_to_sig(const char *str, int len, const float* pore_vals, const int k, const int strand, uint32_t* s_len, float* s_values){

	int i, j, l, pos, n = 0;
	// uint64_t shift1 = 2 * (k - 1);
	uint64_t mask = (1ULL<<2*k) - 1, kmer = 0;
	double mean = 0, std_dev = 0, sum = 0, sum2 = 0, curval = 0;

	for (i = l = j = n = 0; i < len; ++i) {
		if(strand) pos = len - i -1;
		else pos = i;
		int c = seq_nt4_table[(uint8_t)str[pos]];
		if (c < 4) { // not an ambiguous base
			if(!strand) kmer = (kmer << 2 | c) & mask;    // forward k-mer
			// else kmer = (kmer >> 2) | (3ULL^c) << shift1; // reverse k-mer
			//TODO: this is currently based on the ordering in the original ordering in the sigmap implementation. Change later to above
			else kmer = ((kmer << 2) | (3ULL^c)) & mask; // reverse k-mer
		}else
      		kmer = (kmer << 2) & mask; //TODO: This is not the best approach. We are basically inserting 00 (A?) to kmer whenever c >= 4. Mask it instead

		if(i+1 < k) continue;

		curval = pore_vals[kmer];
		s_values[j++] = curval;
		sum += curval;
		sum2 += curval*curval;
	}

	mean = sum/j;
	std_dev = sqrt(sum2/j - (mean)*(mean));

	for(i = 0; i < j; ++i)
		s_values[i] = (s_values[i]-mean)/std_dev;

	*s_len = j;
}

#ifndef NHDF5RH
static inline ri_sig_file_t *ri_sig_open_fast5(const char *fn, const uint8_t mode)
{
	ri_sig_file_t *fp;	
	fp = (ri_sig_file_t*)calloc(1, sizeof(ri_sig_file_t));

	#ifndef NPOD5RH
	fp->pp = NULL;
	#endif

	#ifndef NSLOW5RH
	fp->sp = NULL;
	#endif

	hdf5_tools::File* fast5_file = new hdf5_tools::File();

	if(mode == SIG_WRITE_OP){
		// Create file template
		std::string channel_id = "rawasm_channel";
		fast5_file->create(fn, 1); // Truncate if file is existing

		// Set Root group attributes
		std::string file_type, file_version;
		file_type = "multi-read";
		file_version = "2.2";
		fast5_file->write_attribute("/file_type", file_type);
		fast5_file->write_attribute("/file_version", file_type);

		if (!fast5_file->is_open()) return 0;

		fp->fp = fast5_file;

		fp->num_read = 0;
		fp->cur_read = 0;
	} else if(mode == SIG_READ_OP) {
		fast5_file->open(std::string(fn));

		if (!fast5_file->is_open()) return 0;

		fp->fp = fast5_file;

		bool is_single = false;
		std::vector<std::string> fast5_file_groups = fast5_file->list_group("/");
		fp->num_read = fast5_file_groups.size();
		fp->ch_path = (char**)calloc(fp->num_read, sizeof(char*));
		fp->raw_path = (char**)calloc(fp->num_read, sizeof(char*));

		for (std::string &group : fast5_file_groups) {
			if (group == "Raw") {
				is_single = true;
				break;
			}
		}

		std::string raw_path;
		std::string ch_path;
		int i = 0;

		if (is_single) {
			ch_path = "/UniqueGlobalKey/channel_id";
			for (std::string &read : fast5_file->list_group("/Raw/Reads")) {
				raw_path = "/Raw/Reads/" + read;
				if(i == fp->num_read){
					fprintf(stderr, "ERROR: More reads than previously predicted (%d). Stopped reading the reads here.\n", fp->num_read);
					break;
				}
				fp->ch_path[i] = strdup(ch_path.c_str());
				fp->raw_path[i++] = strdup(raw_path.c_str());
			}
		} else {
			for (std::string &read : fast5_file_groups) {
				raw_path = "/" + read + "/Raw";
				ch_path = "/" + read + "/channel_id";
				fp->ch_path[i] = strdup(ch_path.c_str());
				fp->raw_path[i++] = strdup(raw_path.c_str());
			}
		}

		fp->num_read = i;
		fp->cur_read = 0;
	}
	
	return fp;
}
#endif

#ifndef NPOD5RH
static inline ri_sig_file_t *ri_sig_open_pod5(const char *fn, const uint8_t mode){

	ri_sig_file_t *fp;

	// Allocate signal file
	fp = (ri_sig_file_t*) calloc(1, sizeof(ri_sig_file_t));
	fp->pp = (pod5_file_t *) calloc(1, sizeof(pod5_file_t));

	#ifndef NHDF5RH
	fp->fp = NULL;
	#endif

	#ifndef NSLOW5RH
	fp->sp = NULL;
	#endif


	// Initialize pod5 library
	pod5_init();

	if(!mode){
		// Open file in READ MODE, ready for walking
		Pod5FileReader_t* pod5_file = pod5_open_file(fn);

		if (!pod5_file) {
			fprintf(stderr, "ERROR: Failed to open file (%s) %s\n", fn, pod5_get_error_string());
			free(fp->pp);
			free(fp);
			return 0;
		}

		long unsigned int batch_count = 0;
		if (pod5_get_read_batch_count(&batch_count, pod5_file) != POD5_OK) {
			fprintf(stderr, "Failed to query batch count: %s\n", pod5_get_error_string());
			pod5_close_and_free_reader(pod5_file);
			free(fp->pp);
			free(fp);
			return 0;
		}

		fp->pp->reader = pod5_file;

		fp->num_read = batch_count; //num_read is the number of batches for pod5 (even if it should be the number of rows for batch)
		fp->cur_read = 0; //cur_read is cur_batch for pod5

		Pod5ReadRecordBatch_t* batch = NULL;
		if(pod5_get_read_batch(&batch, pod5_file, fp->cur_read) != POD5_OK){
			fprintf(stderr, "Failed to get batch: %s\n", pod5_get_error_string());
			pod5_close_and_free_reader(pod5_file);
			free(fp->pp);
			free(fp);
			return 0;
		}

		long unsigned int batch_row_count = 0;
		if(pod5_get_read_batch_row_count(&batch_row_count, batch) != POD5_OK) {
			fprintf(stderr, "Failed to get batch row count\n");
			pod5_close_and_free_reader(pod5_file);
			pod5_free_read_batch(batch);
			free(fp->pp);
			free(fp);
			return 0;
		}

		fp->pod5_row_count = batch_row_count;
		fp->pod5_row = 0;

		fp->pod5_batch_count = batch_count;
		fp->pod5_batch = 0;
		
		fp->batch = batch;
		fp->num_read = batch_count;

	} else {

		// Create a file in WRITE MODE
		Pod5FileReader_t* pod5_reader;
		Pod5FileWriter_t* pod5_writer;
		Pod5WriterOptions_t pod5_writer_options = {0, DEFAULT_SIGNAL_COMPRESSION, 0, 0}; // Default options
		int16_t pd5_pore_index;
		int16_t pod5_run_info_index;

		if(pod5_reader = pod5_open_file(fn)){
			fprintf(stderr, "Error, File %s already existing, please remove it\n", fn);
			pod5_close_and_free_reader(pod5_reader);
			free(fp->pp);
			free(fp);
			return 0;
		}

		if(!(pod5_writer = pod5_create_file(fn, "writer", &pod5_writer_options))){
			fprintf(stderr, "Failed to create file %s\n", pod5_get_error_string());
			free(fp->pp);
			free(fp);
			return 0;
		}

		// Create a Pore Dictionary
		pod5_add_pore(&pd5_pore_index, pod5_writer, "rsig_default_pore");
		
		// Create a Run Info Dictionary
		pod5_add_run_info(
			&pod5_run_info_index,
			pod5_writer, 
			"ffffffffffffffffffffffffffffffffffffffff",
			0,
			4095,
			-4096,
			0,
			nullptr,
			nullptr,
			"",
			"",
			"",
			"",
			"",
			0,
			"test_sample",
			4000,
			"",
			"",
			"",
			"",
			"",
			"",
			0,
			nullptr,
			nullptr);


		// Init the new file
		fp->pp->writer = pod5_writer;
		fp->num_read = 0;
		fp->cur_read = 0; 

		fp->pod5_row_count = 0;
		fp->pod5_row = 0;
		fp->pod5_batch_count = 0;
		fp->pod5_batch = 0;
		
		fp->batch = nullptr;

	}


	return fp;
}
#endif

#ifndef NSLOW5RH
static inline ri_sig_file_t *ri_sig_open_slow5(const char *fn, const uint8_t mode){
	
	ri_sig_file_t *fp;

	//open the SLOW5 file
	slow5_file_t *sp = (mode == SIG_WRITE_OP) ? slow5_open(fn, "w") : slow5_open(fn, "r");

	if(sp==NULL){
       fprintf(stderr, "ERROR: Failed to open file %s\n", fn);
	   return 0;
    }

	fp = (ri_sig_file_t*)calloc(1, sizeof(ri_sig_file_t));

	#ifndef NHDF5RH
	fp->fp = NULL;
	#endif

	#ifndef NPOD5RH
	fp->pp = NULL;
	#endif
	
	fp->sp = sp;
	fp->num_read = (mode == SIG_WRITE_OP) ? 0 : ri_sig_count(fp);
	fp->cur_read = 0;

	return fp;
}
#endif

ri_sig_file_t *ri_sig_open(const char *fn, const uint8_t mode){
	if (strstr(fn, ".fast5")) {
		#ifndef NHDF5RH
		return ri_sig_open_fast5(fn, mode);
		#endif
	} else if (strstr(fn, ".pod5") || strstr(fn, ".pod")) {
		#ifndef NPOD5RH
		return ri_sig_open_pod5(fn, mode);
		#endif
	} else if (strstr(fn, ".slow5") || strstr(fn, ".blow5")) {
		#ifndef NSLOW5RH
		return ri_sig_open_slow5(fn, mode);
		#endif
	}

	return 0;
}

void ri_sig_close(ri_sig_file_t *fp)
{
	if(!fp) return;
	// gzclose(fp->fp);

	#ifndef NHDF5RH
	if(fp->fp){
		fp->fp->close();
		// ch_path and raw_path are initialized only in Read Mode (WR mode is not supported atm)
		if(fp->ch_path && fp->raw_path)
			for(int i = 0; i < fp->num_read; ++i){
				if(fp->ch_path[i])free(fp->ch_path[i]);
				if(fp->raw_path[i])free(fp->raw_path[i]);
			}
	}
	#endif
	#ifndef NPOD5RH
	if(fp->pp){
		if(fp->pp->reader) pod5_close_and_free_reader(fp->pp->reader);
		if(fp->pp->writer) pod5_close_and_free_writer(fp->pp->writer);
	}
	#endif

	#ifndef NSLOW5RH
	if(fp->sp){
		slow5_close(fp->sp);
	}
	#endif

	free(fp);
}

ri_sig_file_t *open_sig(const char *fn, const uint8_t mode) //TODO: make this a part of the pipeline. Sequntially reading from many FAST5 files creates an overhead
{
	ri_sig_file_t *fp;
	fp = (ri_sig_file_t*)calloc(1,sizeof(ri_sig_file_t));
	if ((fp = ri_sig_open(fn, mode)) == 0) {
		fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn, strerror(errno));
		ri_sig_close(fp);
		return 0;
	}
	return fp;
}

ri_sig_file_t **open_sigs(int n, const char **fn, const uint8_t mode) //TODO: make this a part of the pipeline. Sequntially reading from many FAST5 files creates an overhead
{
	ri_sig_file_t **fp;
	int i, j;
	fp = (ri_sig_file_t**)calloc(n, sizeof(ri_sig_file_t*));
	for (i = 0; i < n; ++i) {
		if ((fp[i] = ri_sig_open(fn[i], mode)) == 0) {
			fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn[i], strerror(errno));
			for (j = 0; j < i; ++j) ri_sig_close(fp[j]);
			free(fp);
			return 0;
		}
	}
	return fp;
}

//Check if the input const char* A is a directory
//Generated by GitHub Copilot
int is_dir(const char *A)
{
	struct stat st;
	if (stat(A, &st) == -1) return 0;
	return S_ISDIR(st.st_mode);
}

//Recursively find all files that ends with "fast5", "pod5", or "s/blow5" under input directory const char *A
//Generated by GitHub Copilot
void find_sfiles(const char *A, ri_char_v *fnames)
{
	if (!is_dir(A)) {
		if (strstr(A, ".fast5") || strstr(A, ".pod5") || strstr(A, ".pod") || strstr(A, ".slow5") || strstr(A, ".blow5")) {
			char** cur_fname;
			rh_kv_pushp(char*, 0, *fnames, &cur_fname);
			(*cur_fname) = strdup(A);
		}
		return;
	}

	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir(A)) != NULL) {
		while ((ent = readdir(dir)) != NULL) {
			char *tmp = (char*)malloc(strlen(A) + strlen(ent->d_name) + 2);
			sprintf(tmp, "%s/%s", A, ent->d_name);
			if (is_dir(tmp)) {
				if (strcmp(ent->d_name, ".") && strcmp(ent->d_name, ".."))
					find_sfiles(tmp, fnames);
			} else {
				if (strstr(ent->d_name, ".fast5") || strstr(ent->d_name, ".pod5") || strstr(ent->d_name, ".pod") || strstr(ent->d_name, ".slow5") || strstr(ent->d_name, ".blow5")) {
					char** cur_fname;
					rh_kv_pushp(char*, 0, *fnames, &cur_fname);
					(*cur_fname) = strdup(tmp);
				}
			}
			free(tmp);
		}
		closedir(dir);
	}
}

#ifndef NHDF5RH
static inline void ri_read_sig_fast5(ri_sig_file_t* fp, ri_sig_t* s){

	if(fp->cur_read >= fp->num_read) return;
	
	s->name = 0;
	for (auto a : fp->fp->get_attr_map(fp->raw_path[fp->cur_read])) {
		if (a.first == "read_id") {
			s->name = strdup(a.second.c_str());
		}
	}
	assert(s->name);

	float dig = 0, ran = 0, offset = 0;

	// float digitisation = 0, range = 0, offset = 0;
	for (auto a : fp->fp->get_attr_map(fp->ch_path[fp->cur_read])) {
		if (a.first == "channel_number") {
		// channel_idx = atoi(a.second.c_str()) - 1;
		} else if (a.first == "digitisation") {
			dig = atof(a.second.c_str());
		} else if (a.first == "range") {
			ran = atof(a.second.c_str());
		} else if (a.first == "offset") {	
			offset = atof(a.second.c_str());
		}
	}
	

	std::string sig_path = std::string(fp->raw_path[fp->cur_read]) + "/Signal";
	std::vector<float> sig;
	fp->fp->read(sig_path, sig);

	// convert to pA
	uint32_t l_sig = 0;
	float scale = ran/dig;
	float pa = 0;
	for (size_t i = 0; i < sig.size(); i++) {
		pa = (sig[i]+offset)*scale;
		if (pa > 30.0f && pa < 200.0f) {
			sig[l_sig++] = pa;
		}
	}

	s->sig = (float*)calloc(l_sig, sizeof(float));
	s->l_sig = l_sig;
	std::copy(sig.begin(), sig.begin() + l_sig, s->sig);
	s->dig = dig;
	s->ran = ran;
	fp->cur_read++;
}

static inline void ri_write_sig_fast5(ri_sig_file_t* fp, ri_sig_t* s, const float calibration_scale, const float calibration_offset){

	// Define reads attributes 
	std::vector<short int> sig;
	std::string read_name;
	std::string read_id = "00000000-0000-0000-0000-000000000000";
	std::string pore_id = "";
	int read_number;
	std::string run_id = "0000000000000000000000000000000000000000";
	int start_time = 0;
	int duration = 0;
	int start_mux = 0;
	int median_before = 0;
	int end_reason = 0; // UNKNOWN
	// Define channel attributes
	int channel_num = 0;
	float digitisation = 0;
	float range = 1;	// Scale should be range/digitalization 
	float offset = 0;
	float sampling_rate = 0;
	hdf5_tools::File* fast5_file = fp->fp;

	// Convert the signal to int16
	sig.resize(s->l_sig);
	for(int i = 0; i < s->l_sig; i++){
		sig[i] = s->sig[i]/calibration_scale - calibration_offset; 
	}
	
	// Read Attribute
	read_name = s->name;
	fast5_file->write_attribute("/read_" + read_name + "/run_id", run_id);
	fast5_file->write_attribute("/read_" + read_name + "/pore_type", pore_id);
	// Raw signal read signal
	fast5_file->write_dataset("/read_" + read_name + "/Raw/Signal", sig);
	fast5_file->write_attribute("/read_" + read_name + "/Raw/read_number", read_number);
	fast5_file->write_attribute("/read_" + read_name + "/Raw/start_time", start_time);
	fast5_file->write_attribute("/read_" + read_name + "/Raw/duration", duration);
	fast5_file->write_attribute("/read_" + read_name + "/Raw/start_mux", start_mux);
	fast5_file->write_attribute("/read_" + read_name + "/Raw/median_before", median_before);
	fast5_file->write_attribute("/read_" + read_name + "/Raw/end_reason", end_reason);
	fast5_file->write_attribute("/read_" + read_name + "/Raw/read_id", s->rid);
	// Adding channel ID attributes
	fast5_file->write_attribute("/read_" + read_name + "/channel_id/channel_num", channel_num);
	fast5_file->write_attribute("/read_" + read_name + "/channel_id/digitisation", digitisation);
	fast5_file->write_attribute("/read_" + read_name + "/channel_id/range", calibration_scale);
	fast5_file->write_attribute("/read_" + read_name + "/channel_id/offset", calibration_offset);
	fast5_file->write_attribute("/read_" + read_name + "/channel_id/sampling_rate", sampling_rate);
	// Ignoring Context Tags ATM
	// Ignoring Tracking ID ATM

	// Now update file pointer state
	// Here we do not update raw and ch path, since we do not expect to use the file in read and write mode at the same time.
	fp->num_read++;
	sig.clear();

}
#endif

#ifndef NPOD5RH
static inline void ri_read_sig_pod5(ri_sig_file_t* fp, ri_sig_t* s){

	if(fp->cur_read >= fp->num_read) return;
	if(fp->pod5_row >= fp->pod5_row_count) return;

	pod5_init();

    // Open the file ready for walking:
    Pod5FileReader_t* pod5_file = fp->pp->reader;

    if (!pod5_file) {
        fprintf(stderr, "Failed to open file %s\n", pod5_get_error_string());
		return;
    }

	uint16_t read_table_version = 0;
	ReadBatchRowInfo_t read_data;
	if (pod5_get_read_batch_row_info_data(fp->batch, fp->pod5_row, READ_BATCH_ROW_INFO_VERSION, &read_data, &read_table_version) != POD5_OK) {
		fprintf(stderr, "Failed to get read %lu\n", fp->pod5_row);
		return;
	}

	//Retrieve global information for the run
	RunInfoDictData_t* run_info_data;
	if (pod5_get_run_info(fp->batch, read_data.run_info, &run_info_data) != POD5_OK) {
		fprintf(stderr, "Failed to get Run Info %lu %s\n", fp->pod5_row, pod5_get_error_string());
		pod5_release_run_info(run_info_data);
		return;
	}

	// Retrieves the digitisation and range for the channel
	CalibrationExtraData_t calibration_extra_info;
	pod5_get_calibration_extra_info(fp->batch, fp->pod5_row, &calibration_extra_info);
	s->ran = calibration_extra_info.range;
	s->dig = (float)calibration_extra_info.digitisation;

	//This is same as read_data.num_samples
	// unsigned long int sample_count = 0;
	// pod5_get_read_complete_sample_count(fp->pp, fp->batch, fp->pod5_row, &sample_count);

	int16_t *sig = (int16_t*)malloc(read_data.num_samples * sizeof(int16_t));
	float *sigF = (float*)malloc(read_data.num_samples * sizeof(float));
	
	if (pod5_get_read_complete_signal(pod5_file, fp->batch, fp->pod5_row, read_data.num_samples, sig) != POD5_OK) {
		fprintf(stderr, "Failed to get read %lu signal: %s\n", fp->pod5_row, pod5_get_error_string());
		pod5_release_run_info(run_info_data);
		if(sig) free(sig);
		return;
	}

	char read_id_tmp[37];
	pod5_format_read_id(read_data.read_id, read_id_tmp);

	uint32_t l_sig = 0;
	float pa = 0.0f;
	for(uint64_t i = 0; i < read_data.num_samples; i++){
		pa = (sig[i]+read_data.calibration_offset)*read_data.calibration_scale;
		if (pa > 30 && pa < 200) {
			sigF[l_sig++] = pa;
		}
	}

	free(sig);
	s->sig = (float*)calloc(l_sig, sizeof(float));
	s->l_sig = l_sig;
	memcpy(s->sig, sigF, l_sig*sizeof(float));
	s->name = strdup(read_id_tmp);
	free(sigF);

	// fprintf(stderr, "%s %lu\n", read_id_tmp, s->l_sig);
	// for(int i = 0; i < s->l_sig; i++){
	// 	fprintf(stderr, "%f ", s->sig[i]);
	// }
	// fprintf(stderr, "\n");

	pod5_release_run_info(run_info_data);

	fp->pod5_row++;
	if(fp->pod5_row >= fp->pod5_row_count){

		if (pod5_free_read_batch(fp->batch) != POD5_OK) {
            fprintf(stderr, "Failed to release batch\n");
			// pod5_close_and_free_reader(fp->pp);
			return;
        }

		fp->cur_read++;
		if(fp->cur_read < fp->num_read){
			Pod5ReadRecordBatch_t* batch = NULL;
			if(pod5_get_read_batch(&batch, fp->pp->reader, fp->cur_read) != POD5_OK){
				fprintf(stderr, "Failed to get batch: %s\n", pod5_get_error_string());
				return;
			}

			long unsigned int batch_row_count = 0;
			if(pod5_get_read_batch_row_count(&batch_row_count, batch) != POD5_OK) {
				fprintf(stderr, "Failed to get batch row count\n");
				return;
			}

			fp->pod5_row_count = batch_row_count;
			fp->pod5_row = 0;

			fp->batch = batch;
		}
	}

}

static inline void ri_write_sig_pod5(ri_sig_file_t* fp, ri_sig_t* s, const float calibration_scale, const float calibration_offset){

	// We should technically provide all this info from the outside
	// Pointers for the row info
	ReadBatchRowInfoArray_t * pod5_row_data;
	read_id_t * pod5_read_id;
    uint32_t * pod5_read_number;
    uint64_t * pod5_start_sample;
    float * pod5_median_before;
    uint16_t * pod5_channel;
    uint8_t * pod5_well;
    int16_t * pod5_pore_type;
    float * pod5_calibration_offset;
    float * pod5_calibration_scale;
    pod5_end_reason_t * pod5_end_reason;
    uint8_t * pod5_end_reason_forced;
    int16_t * pod5_run_info_id;
    uint64_t * pod5_num_minknow_events;
    float * pod5_tracked_scaling_scale;
    float * pod5_tracked_scaling_shift;
    float * pod5_predicted_scaling_scale;
    float * pod5_predicted_scaling_shift;
    uint32_t * pod5_num_reads_since_mux_change;
    float * pod5_time_since_mux_change;

	// Pointers to signal info
	uint32_t pod5_read_count;
	int16_t ** pod5_signal;
	uint32_t * pod5_signal_size;	

	// Initialize pod5 lib
	pod5_init();

    // Open the file 
    Pod5FileWriter_t* pod5_file = fp->pp->writer;

    if (!pod5_file) {
        fprintf(stderr, "Failed to open file %s\n", pod5_get_error_string());
		return;
    }

	// Add the signal to the POD5 file
	pod5_read_count = 1; // Right now we allow only for one signal write
	pod5_row_data = (ReadBatchRowInfoArray_t *) calloc(pod5_read_count, sizeof(ReadBatchRowInfoArray_t));

	// Allocate the Batch Info Array
	pod5_read_id = 						(read_id_t *)			calloc (pod5_read_count,sizeof(read_id_t));			
	pod5_read_number = 					(uint32_t *)			calloc (pod5_read_count,sizeof(uint32_t));			
	pod5_start_sample = 				(uint64_t *)			calloc (pod5_read_count,sizeof(uint64_t));		
	pod5_median_before = 				(float *)				calloc (pod5_read_count,sizeof(float));				
	pod5_channel =	 					(uint16_t *)			calloc (pod5_read_count,sizeof(uint16_t));		
	pod5_well =		 					(uint8_t *)				calloc (pod5_read_count,sizeof(uint8_t));			
	pod5_pore_type = 					(int16_t *)				calloc (pod5_read_count,sizeof(int16_t));		
	pod5_calibration_offset = 			(float *)				calloc (pod5_read_count,sizeof(float));			
	pod5_calibration_scale = 			(float *)				calloc (pod5_read_count,sizeof(float));		
	pod5_end_reason = 					(pod5_end_reason_t *)	calloc (pod5_read_count,sizeof(pod5_end_reason_t));			
	pod5_end_reason_forced = 			(uint8_t *)				calloc (pod5_read_count,sizeof(uint8_t));				
	pod5_run_info_id = 					(int16_t *)				calloc (pod5_read_count,sizeof(int16_t));			
	pod5_num_minknow_events = 			(uint64_t *)			calloc (pod5_read_count,sizeof(uint64_t));			
	pod5_tracked_scaling_scale = 		(float *)				calloc (pod5_read_count,sizeof(float));				
	pod5_tracked_scaling_shift = 		(float *)				calloc (pod5_read_count,sizeof(float));					
	pod5_predicted_scaling_scale = 		(float *)				calloc (pod5_read_count,sizeof(float));					
	pod5_predicted_scaling_shift = 		(float *)				calloc (pod5_read_count,sizeof(float));				
	pod5_num_reads_since_mux_change =	(uint32_t *)			calloc (pod5_read_count,sizeof(uint32_t));				
	pod5_time_since_mux_change = 		(float *)				calloc (pod5_read_count,sizeof(float));				
	
	for(int i = 0; i < pod5_read_count; i++){

		// Assign to constant pointer
		pod5_row_data[i].read_id = &pod5_read_id[i];
		pod5_row_data[i].read_number = &pod5_read_number[i]; 
		pod5_row_data[i].start_sample = &pod5_start_sample[i]; 
		pod5_row_data[i].median_before = &pod5_median_before[i];
		pod5_row_data[i].channel = &pod5_channel[i];
		pod5_row_data[i].well = &pod5_well[i];
		pod5_row_data[i].pore_type = &pod5_pore_type[i];
		pod5_row_data[i].calibration_offset = &pod5_calibration_offset[i];
		pod5_row_data[i].calibration_scale = &pod5_calibration_scale[i];
		pod5_row_data[i].end_reason = &pod5_end_reason[i];
		pod5_row_data[i].end_reason_forced = &pod5_end_reason_forced[i];
		pod5_row_data[i].run_info_id = &pod5_run_info_id[i];
		pod5_row_data[i].num_minknow_events = &pod5_num_minknow_events[i];
		pod5_row_data[i].tracked_scaling_scale = &pod5_tracked_scaling_scale[i];
		pod5_row_data[i].tracked_scaling_shift = &pod5_tracked_scaling_shift[i];
		pod5_row_data[i].predicted_scaling_scale = &pod5_predicted_scaling_scale[i];
		pod5_row_data[i].predicted_scaling_shift = &pod5_predicted_scaling_shift[i];
		pod5_row_data[i].num_reads_since_mux_change = &pod5_num_reads_since_mux_change[i];
		pod5_row_data[i].time_since_mux_change = &pod5_time_since_mux_change[i];
	
		// Innitialize the array
		pod5_read_number[i] = fp->num_read++;
		pod5_start_sample[i] = 0;
		pod5_median_before[i] = 0;
		pod5_channel[i] = 0;
		pod5_well[i] = 0;
		pod5_pore_type[i] = 0;
		pod5_calibration_offset[i] = calibration_offset;
		pod5_calibration_scale[i] = calibration_scale;
		pod5_end_reason[i] = POD5_END_REASON_UNKNOWN;
		pod5_end_reason_forced[i] = 0;
		pod5_run_info_id[i] = 0;
		pod5_num_minknow_events[i] = 0;
		pod5_tracked_scaling_scale[i] = 0;
		pod5_tracked_scaling_shift[i] = 0;
		pod5_predicted_scaling_scale[i] = 0;
		pod5_predicted_scaling_shift[i] = 0;
		pod5_num_reads_since_mux_change[i] = 0;
		pod5_time_since_mux_change[i] = 0;

		// We need to generate the byte-vector UUID from the string
		uuid_parse(s[i].name, pod5_read_id[i]);
	}

	// Initialize the signal
	pod5_signal_size = (uint32_t *) calloc(pod5_read_count, sizeof(uint32_t));
	pod5_signal = (int16_t **) calloc(pod5_read_count, sizeof(int16_t *));

	// Convert signal from float to int16
	for(int i = 0; i < pod5_read_count; i++) {
		pod5_signal_size[i] = s->l_sig;
		pod5_signal[i] = (int16_t *) calloc(pod5_signal_size[i], sizeof(int16_t));
		for(int j = 0; j < pod5_signal_size[i]; j++){
			pod5_signal[i][j] = s->sig[j]/pod5_calibration_scale[i] - pod5_calibration_offset[i];
		}
	}
	
	// Finally add the signal
	pod5_add_reads_data(pod5_file, pod5_read_count, READ_BATCH_ROW_INFO_VERSION_3, pod5_row_data, (const int16_t **)pod5_signal, pod5_signal_size);

	// Free the allocated arrays
	free(pod5_read_id);
	free(pod5_read_number);	
	free(pod5_start_sample);
	free(pod5_median_before);		
	free(pod5_channel);
	free(pod5_well);		
	free(pod5_pore_type);	
	free(pod5_calibration_offset);	
	free(pod5_calibration_scale);	
	free(pod5_end_reason);
	free(pod5_end_reason_forced);		
	free(pod5_run_info_id);		
	free(pod5_num_minknow_events);	
	free(pod5_tracked_scaling_scale);	
	free(pod5_tracked_scaling_shift);				
	free(pod5_predicted_scaling_scale);				
	free(pod5_predicted_scaling_shift);			
	free(pod5_num_reads_since_mux_change);			
	free(pod5_time_since_mux_change);		
	free(pod5_row_data);

	// Update file info
	fp->pod5_row_count += pod5_read_count; // Update the number of rows
	fp->num_read += pod5_read_count;

}
#endif

#ifndef NSLOW5RH
static inline void ri_read_sig_slow5(ri_sig_file_t* fp, ri_sig_t* s){
	
	if(fp->cur_read >= fp->num_read) return;

	slow5_file_t *sp = fp->sp;
	slow5_rec_t *rec = NULL;
	int ret = 0;
	if ((ret = slow5_get_next(&rec, sp)) < 0) {
		fp->cur_read = fp->num_read;
		return;
	}

	s->name = strdup(rec->read_id);
	float *sigF = (float*)malloc(rec->len_raw_signal * sizeof(float));
	uint32_t l_sig = 0;
	float pa = 0.0f;
	float scale = rec->range/rec->digitisation;

	s->ran = rec->range;
	s->dig = rec->digitisation;
	
	for(int i = 0; i < rec->len_raw_signal; ++i){
		pa = (rec->raw_signal[i]+rec->offset)*scale;
		if (pa > 30.0f && pa < 200.0f) {
			sigF[l_sig++] = pa;
		}
	}

	fp->cur_read++;

	if(ret != SLOW5_ERR_EOF){  //check if proper end of file has been reached
		fp->num_read++;
    }else{
		fp->cur_read = fp->num_read;
	}

	s->sig = (float*)calloc(l_sig, sizeof(float));
	s->l_sig = l_sig;
	memcpy(s->sig, sigF, l_sig*sizeof(float));
	free(sigF);
	slow5_rec_free(rec);
}

static inline void ri_write_sig_slow5(ri_sig_file_t* fp, ri_sig_t* s, const float calibration_scale, const float calibration_offset){
	
	slow5_file_t *sp = fp->sp;
	
	// SLOW5 files needs 3 info: Index (idx), Header (hdx) and Record (rec)
	// What is the Index about?
	// Is the header single for a whole slow5 or we have one for each read?
	// Headers can be associated with at least one auxiliary data. But what are Aux?

	// Records seems to contain all significant information
	slow5_rec_t *read;

	read = slow5_rec_init();
	read->read_group = 0;
	read->read_id_len = 36; // Since it doesn't include the null character
	read->read_id = (char *) calloc(read->read_id_len, sizeof(char));
	memcpy(read->read_id, s->name, read->read_id_len);
	read->digitisation = 0; // remember that digit and range ar functions of calibration_scale
	read->offset = calibration_offset;
	read->range = 0;
	read->sampling_rate = 0;
	read->len_raw_signal = s->l_sig;

	//what about auxiliary fields? (ignored atm)
	read->raw_signal = (short int *) calloc(read->len_raw_signal, sizeof(short int));
	// Convert signal from float to int16
	for(int i = 0; i < read->len_raw_signal; i++) {
		read->raw_signal[i] = s->sig[i]/calibration_scale - calibration_offset;
	}

	read->aux_map = NULL;
	slow5_write(read, sp);
	fp->num_read++;

	// slow5_rec_free also frees read_id and raw_signal
	slow5_rec_free(read);

}
#endif

void ri_read_sig(ri_sig_file_t* fp, ri_sig_t* s){

	assert(fp->cur_read < fp->num_read);

	#ifndef NHDF5RH
	if(fp->fp) ri_read_sig_fast5(fp, s);
	#endif
	#ifndef NPOD5RH
	if(fp->pp) ri_read_sig_pod5(fp, s);
	#endif
	#ifndef NSLOW5RH
	if(fp->sp) ri_read_sig_slow5(fp, s);
	#endif
}

void ri_write_sig(ri_sig_file_t* fp, ri_sig_t* s, const float calibration_scale, const float calibration_offset){
	#ifndef NHDF5RH
	if(fp->fp) ri_write_sig_fast5(fp, s, calibration_scale, calibration_offset);
	#endif
	#ifndef NPOD5RH
	if(fp->pp) ri_write_sig_pod5(fp, s, calibration_scale, calibration_offset);
	#endif
	#ifndef NSLOW5RH
	if(fp->sp) ri_write_sig_slow5(fp, s, calibration_scale, calibration_offset);
	#endif
}

int ri_sig_count(ri_sig_file_t* hd5_fp){

	int sig_cnt = 0;

	#ifndef NHDF5RH
	if(hd5_fp->fp){

		return hd5_fp->num_read;
	}
	#endif
	#ifndef NPOD5RH
	if(hd5_fp->pp){
		Pod5ReadRecordBatch_t* batch = NULL;
		long unsigned int batch_row_count = hd5_fp->pod5_row_count;
		sig_cnt = batch_row_count;

		for(int i = 1; i < hd5_fp->num_read; i++){
			pod5_get_read_batch(&batch, hd5_fp->pp->reader, i);
			pod5_get_read_batch_row_count(&batch_row_count, batch);
			pod5_free_read_batch(batch);
			sig_cnt += batch_row_count;
		}
		return sig_cnt;
	}
	#endif
	#ifndef NSLOW5RH
	if(hd5_fp->sp){
		slow5_rec_t *rec = NULL;
		slow5_file_t* tmp_sp;
		int ret = 0;

		tmp_sp = slow5_open(hd5_fp->sp->meta.pathname, "r");
		if(tmp_sp){
			while ((ret = slow5_get_next(&rec, tmp_sp)) >= 0) sig_cnt++;
			slow5_close(tmp_sp);
		} else return -1;

		hd5_fp->cur_read = 0;
		return sig_cnt;
	}
	#endif

	return -1;

}
