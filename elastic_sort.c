/*
  Copyright (C) 2013- The University of Notre Dame
  This software is distributed under the GNU General Public License.
  See the file COPYING for details.
*/

/* 
  Distributed sort using Work Queue.
*/

#include "debug.h"
#include <work_queue.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

#define LINE_SIZE 2048

//Defaults
#define PARTITION_DEFAULT 20
#define BW_DEFAULT 100 //BW in Mbps

#define PARTITION_COEFF_A_DEFAULT 175
#define MERGE_COEFF_A_DEFAULT 10
#define MERGE_COEFF_B_DEFAULT 380
#define PER_RECORD_SORT_TIME_DEFAULT 0.000003

static unsigned long long total_records = 0;

double partition_overhead_coefficient_a = PARTITION_COEFF_A_DEFAULT;
double merge_overhead_coefficient_a = MERGE_COEFF_A_DEFAULT;
double merge_overhead_coefficient_b = MERGE_COEFF_B_DEFAULT;
double per_record_execution_time = PER_RECORD_SORT_TIME_DEFAULT;
	
static int created_partitions = 0;	

unsigned long long get_total_lines(char *infile) {
	FILE *input_file = fopen(infile, "r");
	unsigned long long line_count = 0; 
	int ch;
	
	while ((ch=fgetc(input_file)) != EOF) {
		if (ch=='\n') 
	        ++line_count;
	}
	fclose(input_file);

	return line_count;
}

//Returns the end byte offset for a given line number in file
off_t get_file_line_end_offset(FILE *fp, off_t start_offset, unsigned long long line_number) {
	unsigned long long line_count = 0; 
	int ch;
	long end_offset = -1;
	
	if (fp == NULL)
		return -1;	
	
	fseek(fp, start_offset, SEEK_SET);	
	
	while ((ch=fgetc(fp)) != EOF && line_count < line_number) {
		if (ch == '\n') 
	        ++line_count;
	}

	if(line_count == line_number) {
		end_offset = ftell(fp);	
	}
	
	return (end_offset-2); //subtract two to rewind back to newline at end of line	
}

int submit_task(struct work_queue *q, const char *command, const char *executable, const char *infile, off_t infile_offset_start, off_t infile_offset_end, const char *outfile) {
	struct work_queue_task *t;
	int taskid;

	char *infile_dup = strdup(infile); //basename() modifies its arguments. So we need to pass a duplicate.
	char *executable_dup = strdup(executable);

	t = work_queue_task_create(command);
	if (!work_queue_task_specify_file_piece(t, infile, basename(infile_dup), infile_offset_start, infile_offset_end, WORK_QUEUE_INPUT, WORK_QUEUE_NOCACHE)) {
		printf("task_specify_file_piece() failed for %s: start offset %ld, end offset %ld.\n", infile, infile_offset_start, infile_offset_end);
		return 0;	
	}
	if (!work_queue_task_specify_file(t, executable, basename(executable_dup), WORK_QUEUE_INPUT, WORK_QUEUE_CACHE)) {
		printf("task_specify_file() failed for %s: check if arguments are null or remote name is an absolute path.\n", executable);
		return 0;	
	}
	if (!work_queue_task_specify_file(t, outfile, outfile, WORK_QUEUE_OUTPUT, WORK_QUEUE_NOCACHE)) {
		printf("task_specify_file() failed for %s: check if arguments are null or remote name is an absolute path.\n", outfile);
		return 0;	
	}

	taskid = work_queue_submit(q, t);
	printf("submitted task (id# %d): %s\n", taskid, t->command_line);

	free(infile_dup);
	free(executable_dup);

	return taskid;
}

/* Partition the input file according to the number of partitions specified and
 * create tasks that sort each of these partitions.
 */
off_t partition_tasks(struct work_queue *q, const char *executable, const char *executable_args, const char *infile, int infile_offset_start, const char *outfile_prefix, int partitions, int records_to_sort) {
	char outfile[256], remote_infile[256], command[256];
	
	off_t file_offset_start = infile_offset_start;
	off_t file_offset_end;
	unsigned long long task_end_line = 0;
	unsigned long long lines_to_submit;
	FILE *infile_fs;

	unsigned long long lines_per_task = (unsigned long long)ceil((double)records_to_sort/partitions); 

	char *infile_dup = strdup(infile);
	strcpy(remote_infile, basename(infile_dup));
	free(infile_dup);

	infile_fs = fopen(infile, "r");
	if (infile_fs == NULL) {
		printf ("Opening %s file failed: %s!\n", infile, strerror(errno)); 
		return 0;	
	}	
	
	while(task_end_line < records_to_sort) {
		//we partition input into pieces by tracking the file offset of the lines in it.
		lines_to_submit = (records_to_sort - task_end_line) < lines_per_task ? (records_to_sort - task_end_line) : lines_per_task;	
		task_end_line += lines_to_submit;
		file_offset_end = get_file_line_end_offset(infile_fs, file_offset_start, lines_to_submit);		
		if (file_offset_end < 0) {
			printf ("End file offset for line %llu is:%ld\n", task_end_line, file_offset_end);
			return 0;	
		}
		
		//create and submit tasks for sorting the pieces.
		sprintf(outfile, "%s.%d", outfile_prefix, created_partitions);
		if (executable_args){	
			sprintf(command, "./%s %s %s > %s", executable, executable_args, remote_infile, outfile);
		} else {
			sprintf(command, "./%s %s > %s", executable, remote_infile, outfile);
		}

		if(!submit_task(q, command, executable, infile, file_offset_start, file_offset_end, outfile))
			return 0;

		created_partitions++;
		file_offset_start = file_offset_end + 1;
	}
	
	fclose(infile_fs);	
	return file_offset_start;
}

int get_file_line_value(FILE *fp) {
	char *line = (char*) malloc(sizeof(char) * LINE_SIZE);	
	int line_value;	
	
	if (!fgets(line, LINE_SIZE, fp)) {
			if(feof(fp)) {
				return -1;	
			}	
	}
	line_value = atoi(line);
	free(line);	
	return line_value;
}

//compute min of array and also return the position of min.
int find_min(int *vals, int length, int *min_pos) {
	int i;
	int min = INT_MAX;
	for (i = 0; i < length; i++) {
		if(vals[i] >= 0 && vals[i] <= min) {
			min = vals[i];
			*min_pos = i;	
		}
	}
	return min;
}

// Do k-way merge of the sorted outputs returned by tasks. 
int merge_sorted_outputs(const char *outfile, const char *partition_file_prefix, int partitions) {
	FILE *outfile_fp;
	
	char partition_file[256];
	FILE **partition_file_fps;
	int *partition_file_line_vals;	
	
	int min_pos, min_value;
	int merged_partitions = 0;
	int i;	

	outfile_fp = fopen(outfile, "w");
	if(!outfile_fp) {
		fprintf(stderr, "Opening file %s failed: %s!\n", outfile, strerror(errno));
		return -1;	
	}	
	
	partition_file_line_vals = malloc(sizeof(int) * partitions);
	partition_file_fps = malloc(sizeof(FILE *) * partitions);

	for(i = 0; i < partitions; i++) {
		sprintf(partition_file, "%s.%d", partition_file_prefix, i);	
		partition_file_fps[i] = fopen(partition_file, "r");
		if(!partition_file_fps[i]) {
			fprintf(stderr, "Opening file %s failed: %s!\n", partition_file, strerror(errno));
			goto cleanup;	
			return -1;	
		}
	}

	//read the first lines of each output file into the array
	for(i = 0; i < partitions; i++) {
		partition_file_line_vals[i] = get_file_line_value(partition_file_fps[i]);
	}

	//compute the minimum of array and load a new value from the contributing
	//file into the array index of the minimum.
	while (merged_partitions < partitions) {
		min_value = find_min(partition_file_line_vals, partitions, &min_pos);
		fprintf(outfile_fp, "%d\n", min_value); //write current min value to output file
		partition_file_line_vals[min_pos] = get_file_line_value(partition_file_fps[min_pos]);
		if (partition_file_line_vals[min_pos] < 0) {	
			merged_partitions++;	
		}
	}

  cleanup:
	for(i = 0; i < partitions; i++) {
		fclose(partition_file_fps[i]);
		sprintf(partition_file, "%s.%d", partition_file_prefix, i);	
		unlink(partition_file);	
	}
	free(partition_file_line_vals);	
	free(partition_file_fps);	
	fclose(outfile_fp);	
	
	return 1;
}

double get_partition_coefficient(char *input_file) {
	
	//Sample the time to partition to empirically compute the partition coefficient 'a' in the model.
	struct timeval current;
	FILE *fp = fopen(input_file, "r");	
	long long unsigned int sample_partition_start_time, sample_partition_end_time, sample_partition_time; 
	
	gettimeofday(&current, 0);
	sample_partition_start_time = ((long long unsigned int) current.tv_sec) * 1000000 + current.tv_usec;
	get_file_line_end_offset(fp, 0, total_records/10); 
	gettimeofday(&current, 0);
	sample_partition_end_time = ((long long unsigned int) current.tv_sec) * 1000000 + current.tv_usec;
	sample_partition_time = sample_partition_end_time - sample_partition_start_time;

	fclose(fp);
	return 10*(sample_partition_time/1000000.0); 
}

double* sort_estimate_runtime(char *input_file, char *executable, int bandwidth, int resources, int tasks) {
	//Model: T(n,k,r) = [T_part + T_merge] + [(t*n)/k * ceil(k/r)] + [(d_n + (d_r * r))/BW_Bps]
	
	double partition_overhead;
	double merge_overhead;
	double parallel_execution_time;
	double transfer_overhead;	
	double total_execution_time;
	double *estimated_times;

	int BW_Bps = bandwidth * 1000000/8; //assume 100Mbps = 100000000/8 Bytes per sec to start
	
	double total_records_in_billion;	
	long long record_bytes = 0;
	long long sw_bytes = 0;
	
	struct stat stat_buf;

	if(!stat(input_file, &stat_buf)){
		record_bytes = stat_buf.st_size;
	}
	
	if(!stat(executable, &stat_buf)){
		sw_bytes = stat_buf.st_size;
	}
		
	if(total_records == 0) {
		total_records = get_total_lines(input_file);
	} 

	total_records_in_billion = total_records/1000000000.0;
	
	//we transfer the records twice - for input and output.
	transfer_overhead = ((double)((2*record_bytes) + (sw_bytes * resources))) / BW_Bps;

	parallel_execution_time = (total_records * per_record_execution_time) / tasks;	
	parallel_execution_time *= ceil((double)tasks/(double)resources);
	
	/* Model of partition is based on the partitioning done in partition_tasks():
	 * Its asymptotic runtime is O(n) where n is number of records in billions and m is number of partitions.
	 * Its actual runtime is modeled as: (a*n). The values of a and b are found by sampling.
	 */	
	if(!partition_overhead_coefficient_a) 
		partition_overhead_coefficient_a = get_partition_coefficient(input_file);

	partition_overhead = (partition_overhead_coefficient_a * total_records_in_billion); 
	
	/* Model of merge is based on the running time of merge_sorted_outputs():
	 * Its asymptotic runtime is O(n*m) where n is number of records in billions and m is number of partitions.
	 * Its actual runtime is modeled as: (a*n*m + b*n). The values of a and b are found sampling.
	 */
	merge_overhead = (merge_overhead_coefficient_a * total_records_in_billion * tasks) + (merge_overhead_coefficient_b * total_records_in_billion);	
	
	total_execution_time = partition_overhead + merge_overhead + parallel_execution_time + transfer_overhead;

	estimated_times = (double *)malloc(sizeof(double) * 5);
	if (estimated_times == NULL) {
		printf ("Allocating memory for estimated_times failed!\n");
		return NULL;
	}
	estimated_times[0] = total_execution_time;
	estimated_times[1] = partition_overhead;
	estimated_times[2] = merge_overhead;
	estimated_times[3] = parallel_execution_time;
	estimated_times[4] = transfer_overhead;
	return estimated_times;
}

int print_optimal_runtimes(char *input_file, char *executable, int bandwidth, int resources, double *optimal_times) {
	double *estimated_times;
	double optimal_execution_time = -1;
	double execution_time = -1;
	int optimal_partitions;
	int i;
	for (i = 1; i <= 5*resources; i++) { 	
		estimated_times = sort_estimate_runtime(input_file, executable, bandwidth, resources, i);
		execution_time = estimated_times[0];
		if (optimal_execution_time < 0 || execution_time < optimal_execution_time) {
			optimal_execution_time = optimal_times[0] = execution_time;
			optimal_partitions = i;
			optimal_times[1] = estimated_times[1];
			optimal_times[2] = estimated_times[2];
			optimal_times[3] = estimated_times[3];
			optimal_times[4] = estimated_times[4];
		}	
		free(estimated_times);	
	}
	return optimal_partitions;
}

static void show_help(const char *cmd) {
    fprintf(stdout, "Use: %s [options] <sort program> <file 1>\n", cmd);
	fprintf(stdout, "where options are:\n");
	fprintf(stdout, " %-30s Specify a project name for the Work Queue master. (default = none)\n", "-N <string>");
	fprintf(stdout, " %-30s Specify the number of partitions to create of the input data. (default = 20)\n", "-k <int>");
	fprintf(stdout, " %-30s Automatically determine the optimal partition size. (default = 20)\n", "-A <int>");
	fprintf(stdout, " %-30s Empirically estimate the model coefficients by sampling the execution environment. (default = off)\n", "-S <int>");
	fprintf(stdout, " %-30s Specify the arguments for the sort program.\n", "-p <string>");
	fprintf(stdout, " %-30s Estimate and print the optimal number of partitions for different resource sizes and exit.\n", "-M");
	fprintf(stdout, " %-30s Specify the number of records in the input file.(default=auto).\n", "-L <int>");
	fprintf(stdout, " %-30s Specify the keepalive interval for WQ.(default=300).\n", "-I <int>");
	fprintf(stdout, " %-30s Specify the keepalive timeout for WQ.(default=30).\n", "-T <int>");
	fprintf(stdout, " %-30s Estimate and print the runtime for specified partition and exit.\n", "-R <int>");
	fprintf(stdout, " %-30s Set the estimated bandwidth to workers for estimating optimal paritions. (default=%d)\n", "-B <int>", BW_DEFAULT);
	fprintf(stdout, " %-30s Show this help screen\n", "-h,--help");
}

int main(int argc, char *argv[])
{
	struct work_queue *q;
	int port = 9000;
	int c;

	int bandwidth = BW_DEFAULT; 
	
	char *sort_arguments = NULL;
	const char *proj_name = NULL;
	int auto_partition = 0;	
	int sample_env = 0;	
	int print_runtime_estimates = 0;
	int estimate_partition= 0;
	struct timeval current;
	long long unsigned int execn_start_time, execn_time;
	int keepalive_interval = 300;
	int keepalive_timeout = 30;

	int partitions = PARTITION_DEFAULT;

	gettimeofday(&current, 0);
	execn_start_time = ((long long unsigned int) current.tv_sec) * 1000000 + current.tv_usec;

	debug_flags_set("all");
	if(argc < 3) {
		show_help(argv[0]);
		return 0;
	}
		
	while((c = getopt(argc, argv, "N:k:ASp:MR:L:I:T:B:h")) != (char) -1) {
		switch (c) {
		case 'N':
			proj_name = strdup(optarg);
			break;
		case 'k':
			partitions = atoi(optarg);
			break;
		case 'A':
			auto_partition = 1;
			break;
		case 'S':
			sample_env = 1;
			break;
		case 'p':
			sort_arguments = strdup(optarg);
			break;
		case 'M':
			print_runtime_estimates = 1;
			break;
		case 'R':
			estimate_partition = atoi(optarg);
			break;
		case 'L':
			total_records = atoll(optarg);
			break;
		case 'I':
			keepalive_interval = atoi(optarg);
			break;
		case 'T':
			keepalive_timeout = atoi(optarg);
			break;
		case 'B':
			bandwidth = atoi(optarg);
			break;
		case 'h':
			show_help(argv[0]);
			return 0;
		default:
			show_help(argv[0]);
			return -1;
		}
	}

	char sort_executable[256]; 
	char infile[256], outfile_prefix[256]; 
	struct work_queue_task *t;	
	off_t last_partition_offset_end = 0;
	char *infile_temp = NULL;
	int i;
	int records;	
	int optimal_partitions;
	int optimal_resources; 
	int current_optimal_partitions;
	double current_optimal_time = DBL_MAX;
	double optimal_times[5];

	sprintf(sort_executable, "%s", argv[optind]);
	sprintf(infile, "%s", argv[optind+1]);

	infile_temp = strdup(infile);		
	strcpy(outfile_prefix, basename(infile_temp));
	sprintf(outfile_prefix, "%s.sorted", outfile_prefix);
	free(infile_temp);

	if(estimate_partition) {
		double *estimated_runtimes = (double *)malloc(sizeof(double) * 5); 
		for (i = 1; i <= 2*estimate_partition; i++) {
			estimated_runtimes = sort_estimate_runtime(infile, sort_executable, bandwidth, i, estimate_partition); 
			if(estimated_runtimes[0] < current_optimal_time) {
				current_optimal_time = estimated_runtimes[0];
				optimal_times[0] = estimated_runtimes[0];
				optimal_times[1] = estimated_runtimes[1];
				optimal_times[2] = estimated_runtimes[2];
				optimal_times[3] = estimated_runtimes[3];
				optimal_times[4] = estimated_runtimes[4];
				optimal_resources = i;
			}
		}	
		printf("For partition %d: %d %f %f %f %f %f\n", estimate_partition, optimal_resources, optimal_times[0], optimal_times[1], optimal_times[2], optimal_times[3], optimal_times[4]);	
		free(estimated_runtimes);	
		return 1;	
	}

	if(print_runtime_estimates) {
		printf("Resources \t Partitions \t Runtime \t Part time \t Merge time\n");
		for (i = 1; i <= 100; i++) {
			optimal_partitions = print_optimal_runtimes(infile, sort_executable, bandwidth, i, optimal_times); 
			printf("%d \t \t %d \t %f \t %f \t %f \t %f \t %f\n", i, optimal_partitions, optimal_times[0], optimal_times[1], optimal_times[2], optimal_times[3], optimal_times[4]);	
		}
		return 1;	
	}

	if(auto_partition) {
		printf("Determining optimal partition size for %s\n", infile);
		for (i = 1; i <= 100; i++) {
			current_optimal_partitions = print_optimal_runtimes(infile, sort_executable, bandwidth, i, optimal_times); 
			if (optimal_times[0] < current_optimal_time) {
				current_optimal_time = optimal_times[0];	
				optimal_partitions = current_optimal_partitions;
				optimal_resources = i;	
			}
		}
		printf("Optimal partition size is %d that runs the workload in %f\n", optimal_partitions, current_optimal_time);	
		printf("--> Please allocate %d resources for running this workload in a cost-efficient manner.\n", optimal_resources);	
		partitions = optimal_partitions;	
	}

	q = work_queue_create(port);
	if(!q) {
		printf("couldn't listen on port %d: %s\n", port, strerror(errno));
		return 1;
	}

	printf("listening on port %d...\n", work_queue_port(q));
	
	if(proj_name){
		work_queue_specify_master_mode(q, WORK_QUEUE_MASTER_MODE_CATALOG);	
		work_queue_specify_name(q, proj_name);
	}
	work_queue_specify_keepalive_interval(q, keepalive_interval);
	work_queue_specify_keepalive_timeout(q, keepalive_timeout);
	//work_queue_activate_fast_abort(q, 2.5);

	free((void *)proj_name);

	printf("%s will be run to sort contents of %s\n", sort_executable, infile);

	if(total_records == 0) {
		total_records = get_total_lines(infile);
	} 
	records = total_records;

	if(sample_env) {
		//Reset the coefficients to 0 so we use the empirically determined values from sampling the execution environment.
		partition_overhead_coefficient_a = 0;
		merge_overhead_coefficient_a = 0;
		merge_overhead_coefficient_b = 0;
		per_record_execution_time = 0;
	}

	long long unsigned int part_start_time, part_end_time, part_time;
	gettimeofday(&current, 0);
	part_start_time = ((long long unsigned int) current.tv_sec) * 1000000 + current.tv_usec;

	last_partition_offset_end = partition_tasks(q, sort_executable, sort_arguments, infile, 0, outfile_prefix, partitions, records);
    	
	gettimeofday(&current, 0);
	part_end_time = ((long long unsigned int) current.tv_sec) * 1000000 + current.tv_usec;
	part_time = part_end_time - part_start_time;
	printf("Partition time is %llu\n", part_time);
	
	free(sort_arguments);

	FILE *task_file = fopen("wq_sort.tasktimes", "w");
   	if (!task_file) {
        	printf("Opening of wq_sort.tasktimes file failed!\n");
        	return 1;
    	}
	
	printf("Waiting for tasks to complete...\n");
	long long unsigned int parallel_start_time, parallel_end_time, parallel_time;
	gettimeofday(&current, 0);
	parallel_start_time = ((long long unsigned int) current.tv_sec) * 1000000 + current.tv_usec;
	while(!work_queue_empty(q)) {
		t = work_queue_wait(q, 5);
		if(t) {
			printf("Task (taskid# %d) complete: %s (return code %d)\n", t->taskid, t->command_line, t->return_status);
			fprintf(task_file, "%d: %llu\n", t->taskid, (long long unsigned) t->cmd_execution_time);	
			work_queue_task_delete(t);
		}
	}
	fclose(task_file);
	
	gettimeofday(&current, 0);
	parallel_end_time = ((long long unsigned int) current.tv_sec) * 1000000 + current.tv_usec;
	parallel_time = parallel_end_time - parallel_start_time;
	printf("Parallel execution time is %llu\n", parallel_time);
	
	long long unsigned int merge_start_time, merge_end_time, merge_time;
	gettimeofday(&current, 0);
	merge_start_time = ((long long unsigned int) current.tv_sec) * 1000000 + current.tv_usec;

	merge_sorted_outputs(outfile_prefix, outfile_prefix, partitions);	
	
	gettimeofday(&current, 0);
	merge_end_time = ((long long unsigned int) current.tv_sec) * 1000000 + current.tv_usec;
	merge_time = merge_end_time - merge_start_time;
	printf("Merge time is %llu\n", merge_time);
	
	printf("Sorting complete. Output is at: %s!\n", outfile_prefix);

	execn_time = merge_end_time - execn_start_time;
	printf("Execn time is %llu\n", execn_time);

	FILE *time_file = fopen("wq_sort.times", "w");
	if (time_file) {
		fprintf(time_file, "Partition time: %llu\n", part_time);
		fprintf(time_file, "Parallel time: %llu\n", parallel_time);
		fprintf(time_file, "Merge time: %llu\n", merge_time);
		fprintf(time_file, "Execution time: %llu\n", execn_time);
	}
	fclose(time_file);

	work_queue_delete(q);
	return 0;
}
