/**********************************************************************
  Copyright(c) 2011-2013 Intel Corporation All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.
    * Neither the name of Intel Corporation nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************/
#define __USE_GNU 
#include <sched.h>
#include <pthread.h>
#include <getopt.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "raid.h"
#include "cauchy.h"
#include "galois.h"
#include "coded.h"
#include "vandermonde.h"

//#include "types.h"
#include <time.h>
#define TEST_SOURCES 15
#define TEST_LEN     16*1024
#define MMAX 255
#define KMAX 255
int codenum=0;
int alldata=0;
int decodedata=0;
long long int encode_time_arrs[15500];
long long int decode_time_arrs[15500];
int encode_test(int k, int p, int w, int len, int packetsize,int e);
int encode_test_split(int k, int p, int w, int len, int packetsize,int e,int n);
int usage(void)
{
	fprintf(stderr,
		"Usage: ec_simple_example [options]\n"
		"  -h        Help\n"
		"  -k <val>  Number of source fragments\n"
		"  -p <val>  Number of parity fragments\n"
		"  -w <val>  \n"
		"  -len <val>  Length of fragments\n"
		"  -packetsize \n");
	exit(0);
}

int en_n=1000;
int de_n=1000;
int *b;



int main(int argc, char *argv[])
{
	int i,j;
	int c;
	int k=10, p=2, w=8;
	int packetsize=1024;
	int e=1;
	int d=p;
	int n=1;
	while ((c = getopt(argc, argv, "k:p:s:e:w:d:n:")) != -1) {
		switch (c) {
		case 'k':
			k = atoi(optarg);
			break;
		case 'p':
			p = atoi(optarg);
			break;
		case 's':
			packetsize=atoi(optarg);
			break;
		case 'e':
			e=atoi(optarg);
			break;
		case 'w':
			w=atoi(optarg);
			break;
		case 'd':
			d=atoi(optarg);
			break;
		case 'n':
			n=atoi(optarg);
			break;
		break;
		}
	}
	
	int len=10*1024*1024;

	encode_test(k,p,w,len,packetsize,e);
	//encode_test_split(k,p,w,len,packetsize,e,n);
	long long int encode_time_all=0;
	long long int decode_time_all=0;
	for(i=0;i < en_n;i++)
	{
		encode_time_all += encode_time_arrs[i];
	}

	for(i=0;i < de_n;i++)
	{
		decode_time_all+=decode_time_arrs[i];
	}
	encode_time_all=encode_time_all/en_n;
	decode_time_all=decode_time_all/de_n;
	double encode_v=((((double)alldata))/1024/1024/1024)/((double)encode_time_all/1000/1000/1000);
	double decode_v=((((double)decodedata))/1024/1024/1024)/((double)decode_time_all/1000/1000/1000);

	printf("%lf %lf\n",encode_v,decode_v);
	return 0;
}

int encode_test(int k, int p, int w, int len, int packetsize,int e)
{
	struct timespec time1 = {0, 0};
    struct timespec time2 = {0, 0};
	long encode_time=0;

	int i, j, m,a;
	int alignlen;
	// Fragment buffer pointers
	u8 *frag_ptrs[MMAX];
	u8 *copy_ptrs[MMAX];
	u8 *tmp[MMAX];
	unsigned int* encode_matrix;
	int* encode_bitmatrix;
	int* slp;
	int** schedule;
	m=k+p;

	int datasize=10*1024*1024;
	alignlen = datasize/k;
	if (packetsize != 0) {
		if (alignlen%(w*packetsize) != 0) { 
			while (alignlen%(w*packetsize) != 0) 
				alignlen++;
		}
	}
	else {
		if (alignlen%w != 0) {
			while (alignlen%w != 0) 
				alignlen++;
		}
	}
    alldata = alignlen*(k);
	decodedata = alignlen*(k);
	// Allocate the src & parity buffers
	for (i = 0; i < m; i++) {
		void *buf;
        if (posix_memalign(&buf, 32, alignlen)) {
            printf("alloc error: Fail");
            return 1;
        }
        frag_ptrs[i] = buf;
		//frag_ptrs_copy[i] = frag_ptrs[i];
	}

	for (i = 0; i < m; i++)
		for (j = 0; j < alignlen; j++)
			frag_ptrs[i][j] = rand();

	encode_matrix=cauchy_good_general_coding_matrix(k,p,w);
	//encode_matrix=gen_best_matrix_all(k,p,w);
	encode_bitmatrix=matrix_to_bitmatrix(k, p, w, encode_matrix);
	
	schedule=bitmatrix_to_schedule(k,p,w,encode_bitmatrix);

	for(i=0;i<en_n;i++){
		clock_gettime(CLOCK_REALTIME, &time1);
		encode_best(k, p, w, schedule, frag_ptrs, &frag_ptrs[k],alignlen, packetsize);
		clock_gettime(CLOCK_REALTIME, &time2);
		encode_time=(time2.tv_sec-time1.tv_sec)*1000000000+(time2.tv_nsec-time1.tv_nsec);
		encode_time_arrs[i]=encode_time;
	}

	u8 **ptrs;
	int *decode_bitmatrix;
	int **de_operations;
	int *eraseds;
	int *survives;
	int *recover;
	int *de_recover;

	eraseds = (int *)malloc(sizeof(int)*(m));
	survives = (int *)malloc(sizeof(int)*(m));
	recover = (int *)malloc(sizeof(int)*(m));
	de_recover = (int *)malloc(sizeof(int)*(m));
	int erro_num=e;

	for(i=0;i < erro_num;i++)
	{
		eraseds[i]=i;
	}
	eraseds[erro_num] = -1;
	for(i=0;i<m-erro_num;i++)
	{
		survives[i]=i+erro_num;
		de_recover[i] = i+ erro_num;
	}

	decode_bitmatrix = generate_decoding_bitmatrix(k, p, w, encode_bitmatrix, eraseds,de_recover);
	ptrs = set_up_ptrs_for_opertaions(k, erro_num, eraseds, frag_ptrs, &frag_ptrs[k],de_recover);
	schedule=bitmatrix_to_schedule(k,erro_num,w,decode_bitmatrix);

	for(i=0;i<de_n;i++)
	{
		clock_gettime(CLOCK_REALTIME, &time1);
		encode_best(k, erro_num, w, schedule, ptrs, &ptrs[k],alignlen, packetsize);
		clock_gettime(CLOCK_REALTIME, &time2);
		encode_time=(time2.tv_sec-time1.tv_sec)*1000000000+(time2.tv_nsec-time1.tv_nsec);
		decode_time_arrs[i]=encode_time;
	}
	
	return 1;
}


int encode_test_split(int k, int p, int w, int len, int packetsize,int e,int n)
{	
	struct timespec time1 = {0, 0};
    struct timespec time2 = {0, 0};
	long encode_time=0;

	int i, j, m,a;
	int alignlen;
	// Fragment buffer pointers
	u8 *frag_ptrs[MMAX];
	u8 *copy_ptrs[MMAX];
	u8 *tmp[MMAX];
	unsigned int* encode_matrix;
	int* encode_bitmatrix;
	int* slp;
	int** schedule;
 	int ***operations;
	m=k+p;

	int datasize=10*1024*1024;
	//int datasize=10321920;
	alignlen = datasize/k;
	//align len with w*packetsize
	if (packetsize != 0) {
		if (alignlen%(w*packetsize) != 0) { 
			while (alignlen%(w*packetsize) != 0) 
				alignlen++;
		}
	}
	else {
		if (alignlen%w != 0) {
			while (alignlen%w != 0) 
				alignlen++;
		}
	}
    alldata = alignlen*(k);
	decodedata = alignlen*(k);
	// Allocate the src & parity buffers
	for (i = 0; i < m; i++) {
		void *buf;
        if (posix_memalign(&buf, 32, alignlen)) {
            printf("alloc error: Fail");
            return 1;
        }
        frag_ptrs[i] = buf;
	}

	for (i = 0; i < m; i++)
		for (j = 0; j < alignlen; j++)
			frag_ptrs[i][j] = rand();

	//encode_matrix=gen_best_matrix_all(k,p,w);
	encode_matrix=cauchy_good_general_coding_matrix(k,p,w);

	encode_bitmatrix=matrix_to_bitmatrix(k, p, w, encode_matrix);
	operations=bitmatrix_to_schedule_split(k,p,w,n,encode_bitmatrix);

	void **data_indexs;
	void **parity_indexs;
	long *buf;
	posix_memalign(&buf, 32, p*w*2*k*w*sizeof(long));
	data_indexs=buf;
	posix_memalign(&buf, 32, 2*2*k*w*sizeof(long));
	parity_indexs=buf;
	void **data_buffs;
	void **parity_buffs;
	data_buffs=malloc(sizeof(void *) * (2*k*w));
	parity_buffs=malloc(sizeof(void *));
	for(i=0;i<en_n;i++){
		clock_gettime(CLOCK_REALTIME, &time1);

		encode_best_spilt(k, p, w, operations, frag_ptrs, &frag_ptrs[k],n,alignlen, packetsize);

		clock_gettime(CLOCK_REALTIME, &time2);
		encode_time=(time2.tv_sec-time1.tv_sec)*1000000000+(time2.tv_nsec-time1.tv_nsec);
		encode_time_arrs[i]=encode_time;
	}

	u8 **ptrs;
	int *decode_bitmatrix;
	int **de_operations;
	int *eraseds;
	int *survives;
	int *recover;
	int *de_recover;

	eraseds = (int *)malloc(sizeof(int)*(m));
	survives = (int *)malloc(sizeof(int)*(m));
	recover = (int *)malloc(sizeof(int)*(m));
	de_recover = (int *)malloc(sizeof(int)*(m));
	int erro_num=e;

	for(i=0;i < erro_num;i++)
	{
		eraseds[i]=i;
	}
	eraseds[erro_num] = -1;
	for(i=0;i<m-erro_num;i++)
	{
		survives[i]=i+erro_num;
		de_recover[i] = i+ erro_num;
	}
	decode_bitmatrix = generate_decoding_bitmatrix(k, erro_num, w, encode_bitmatrix, eraseds,de_recover);
	ptrs = set_up_ptrs_for_opertaions(k, erro_num, eraseds, frag_ptrs, &frag_ptrs[k],de_recover);
	operations=bitmatrix_to_schedule_split(k,erro_num,w,n,decode_bitmatrix);
	
	for(i=0;i<de_n;i++)
	{
		clock_gettime(CLOCK_REALTIME, &time1);
	
		encode_best_spilt(k, erro_num, w, operations, ptrs, &ptrs[k],n,alignlen, packetsize);

		clock_gettime(CLOCK_REALTIME, &time2);
		encode_time=(time2.tv_sec-time1.tv_sec)*1000000000+(time2.tv_nsec-time1.tv_nsec);
		decode_time_arrs[i]=encode_time;
	}
	return 1;
}


