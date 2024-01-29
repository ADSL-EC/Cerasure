/**********************************************************************
  Copyright(c) 2011-2018 Intel Corporation All rights reserved.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include "erasure_code.h"	// use <isa-l.h> instead when linking against installe
#include <time.h>
#include "test.h"
#define MMAX 255
#define KMAX 255

typedef unsigned char u8;
int encode_time_arrs[15500];

int usage(void)
{
	fprintf(stderr,
		"Usage: ec_simple_example [options]\n"
		"  -h        Help\n"
		"  -k <val>  Number of source fragments\n"
		"  -p <val>  Number of parity fragments\n"
		"  -l <val>  Length of fragments\n"
		"  -e <val>  Simulate erasure on frag index val. Zero based. Can be repeated.\n"
		"  -r <seed> Pick random (k, p) with seed\n");
	exit(0);
}

static int gf_gen_decode_matrix_simple(u8 * encode_matrix,
				       u8 * decode_matrix,
				       u8 * invert_matrix,
				       u8 * temp_matrix,
				       u8 * decode_index,
				       u8 * frag_err_list, int nerrs, int k, int m);
//--------------------------测试编解码性能----------------------------
int main(int argc, char *argv[])
{
	cpu_set_t mask;  //CPU核的集合
    CPU_ZERO(&mask);    //置空
    CPU_SET(30,&mask);

	struct timespec time1 = {0, 0};
    struct timespec time2 = {0, 0};
	long encode_time=0;

	void *buf;
	int i, j, m, c, e, ret;
	int k = 10, p = 4, len = 8 * 1024;	// Default params
	int nerrs = 0;
	// Fragment buffer pointers
	u8 *frag_ptrs[MMAX];
	u8 *recover_srcs[KMAX];
	u8 *recover_outp[KMAX];
	u8 frag_err_list[MMAX];

	// Coefficient matrices
	u8 *encode_matrix, *decode_matrix;
	u8 *invert_matrix, *temp_matrix;
	u8 *g_tbls;
	u8 decode_index[MMAX];
	// if (argc == 1)


	while ((c = getopt(argc, argv, "k:p:l:e:r:h:")) != -1) {
		switch (c) {
		case 'k':
			k = atoi(optarg);
			break;
		case 'p':
			p = atoi(optarg);
			break;
			case 'e':
			e = atoi(optarg);
			break;
		}
	}
	m = k + p;
	int datasize=10*1024*1024;
	len=datasize/k;
	for (i = 0; i < e; i++)
		frag_err_list[nerrs++] = i;
	// Check for valid parameters
	if (m > MMAX || k > KMAX || m < 0 || p < 1 || k < 1) {
		printf(" Input test parameter error m=%d, k=%d, p=%d, erasures=%d\n",
		       m, k, p, nerrs);
		usage();
	}
	if (nerrs > p) {
		printf(" Number of erasures chosen exceeds power of code erasures=%d p=%d\n",
		       nerrs, p);
		usage();
	}
	for (i = 0; i < nerrs; i++) {
		if (frag_err_list[i] >= m) {
			printf(" fragment %d not in range\n", frag_err_list[i]);
			usage();
		}
	}

	//printf("ec_simple_example:\n");

	// Allocate coding matrices
	encode_matrix = malloc(m * k);
	decode_matrix = malloc(m * k);
	invert_matrix = malloc(m * k);
	temp_matrix = malloc(m * k);
	g_tbls = malloc(k * p * 32);

	if (encode_matrix == NULL || decode_matrix == NULL
	    || invert_matrix == NULL || temp_matrix == NULL || g_tbls == NULL) {
		printf("Test failure! Error with malloc\n");
		return -1;
	}
	// Allocate the src & parity buffers
	for (i = 0; i < m; i++) {
		if (posix_memalign(&buf, 32, len)) {
			printf("alloc error: Fail\n");
			return -1;
		}
		frag_ptrs[i] = buf;
	}
	// for (i = 0; i < m; i++) {
	// 	if (NULL == (frag_ptrs[i] = malloc(len))) {
	// 		printf("alloc error: Fail\n");
	// 		return -1;
	// 	}
	// }

	// Allocate buffers for recovered data
	for (i = 0; i < p; i++) {
		if (NULL == (recover_outp[i] = malloc(len))) {
			printf("alloc error: Fail\n");
			return -1;
		}
	}

	// Fill sources with random data
	for (i = 0; i < k; i++)
		for (j = 0; j < len; j++)
			frag_ptrs[i][j] = rand();
	

	//printf(" encode (m,k,p)=(%d,%d,%d) len=%d\n", m, k, p, len);

	// Pick an encode matrix. A Cauchy matrix is a good choice as even
	// large k are always invertable keeping the recovery rule simple.
	gf_gen_cauchy1_matrix(encode_matrix, m, k);
    //u8 jerasure_matrix[24] =  {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,1,1,1,1,1,2,142,4};
    // printf("encode_matrix_begin\n");
    // for(i = 0;i < m*k; i++)
    // {
    //     //encode_matrix[i] = jerasure_matrix[i];
    //     printf("%d ",encode_matrix[i]);
    //     if((i+1) % k == 0)
    //     {
    //         printf("\n");
    //     }
    // }

	// Initialize g_tbls from encode matrix
	ec_init_tables(k, p, &encode_matrix[k * k], g_tbls);

	int n=1000;
	int de_n=1;
	printf("%d %d\n",k,p);
	clock_gettime(CLOCK_REALTIME, &time1);
	//ec_encode_data(len, k, p, g_tbls, frag_ptrs, &frag_ptrs[k]);
	int a=0;
	for(i=0;i<n;i++){
		//clock_gettime(CLOCK_REALTIME, &time1);
		ec_encode_data(len, k, p, g_tbls, frag_ptrs, &frag_ptrs[k]);
		//clock_gettime(CLOCK_REALTIME, &time2);
		//encode_time=(time2.tv_sec-time1.tv_sec)*1000000000+(time2.tv_nsec-time1.tv_nsec);
		//encode_time_arrs[i]=encode_time;
		//printf("encode_time is %ld   ", encode_time);
		//printf("encode throughput(MB/s): %0.10f\n", (((double)len*m)/1024/1024)/((double)encode_time/1000/1000/1000));
	}
	clock_gettime(CLOCK_REALTIME, &time2);
	encode_time=(time2.tv_sec-time1.tv_sec)*1000000000+(time2.tv_nsec-time1.tv_nsec);
	encode_time/=n;
	//printf("encode_time is %ld \n",encode_time);
	//printf("encode throughput(MB/s): %d\n", (int)((((double)len*m)/1024/1024)/((double)encode_time/1000/1000/1000)));
	printf("%lf\n", ((((double)len*k)/1024/1024/1024)/((double)encode_time/1000/1000/1000)));
	if (nerrs <= 0)
		return 0;

	//printf(" recover %d fragments\n", nerrs);

	// Find a decode matrix to regenerate all erasures from remaining frags
	ret = gf_gen_decode_matrix_simple(encode_matrix, decode_matrix,
					  invert_matrix, temp_matrix, decode_index,
					  frag_err_list, nerrs, k, m);
	if (ret != 0) {
		printf("Fail on generate decode matrix\n");
		return -1;
	}
	// Pack recovery array pointers as list of valid fragments
	for (i = 0; i < k; i++)
		recover_srcs[i] = frag_ptrs[decode_index[i]];

	// Recover data
	ec_init_tables(k, nerrs, decode_matrix, g_tbls);
	clock_gettime(CLOCK_REALTIME, &time1);
	for(i=0;i<de_n;i++)
	{
		ec_encode_data(len, k, nerrs, g_tbls, recover_srcs, recover_outp);
	}
	clock_gettime(CLOCK_REALTIME, &time2);
	encode_time=(time2.tv_sec-time1.tv_sec)*1000000000+(time2.tv_nsec-time1.tv_nsec);
	encode_time/=de_n;
	printf("%lf\n", ((((double)len*(k))/1024/1024/1024.0)/((double)encode_time/1000/1000/1000)));
	// Check that recovered buffers are the same as original
	for(i=0;i<p;i++)
	{
		for(j=0;j<len;j++)
		{
			if(frag_ptrs[k+i][j]!=recover_outp[i][j])
			{
				printf("erro");
			}
		}
	}
	return 0;
}

// /*
//  * Generate decode matrix from encode matrix and erasure list
//  *
//  */


//---------------------------------------------测试LOAD XOR STORE GF 时间--------------------------------------------------
// int main(int argc, char *argv[])
// {
// 	cpu_set_t mask;  //CPU核的集合
//     CPU_ZERO(&mask);    //置空
//     CPU_SET(30,&mask);

// 	struct timespec time1 = {0, 0};
//     struct timespec time2 = {0, 0};
// 	long encode_time=0;

// 	void *buf;
// 	int i, j, m, c, e, ret;
// 	int k = 100, p = 4, len = 8 * 1024;	// Default params
// 	int nerrs = 0;
// 	// Fragment buffer pointers
// 	u8 *frag_ptrs[MMAX];
// 	u8 *recover_srcs[KMAX];
// 	u8 *recover_outp[KMAX];
// 	u8 frag_err_list[MMAX];

// 	// Coefficient matrices
// 	u8 *encode_matrix, *decode_matrix;
// 	u8 *invert_matrix, *temp_matrix;
// 	u8 *g_tbls;
// 	u8 decode_index[MMAX];
// 	// if (argc == 1)



// 	m = k+p;
// 	int datasize=10*1024*1024;
// 	len=datasize/k;
// 	e=0;

// 	// Check for valid parameters
	
// 	//printf("ec_simple_example:\n");

// 	// Allocate coding matrices
// 	encode_matrix = malloc(m * k);
// 	decode_matrix = malloc(m * k);
// 	invert_matrix = malloc(m * k);
// 	temp_matrix = malloc(m * k);
// 	g_tbls = malloc(k * p * 32);
// 	// Allocate the src & parity buffers

// 	for (i = 0; i < m; i++) {
// 		if (posix_memalign(&buf, 32, len)) {
// 			printf("alloc error: Fail\n");
// 			return -1;
// 		}
// 		frag_ptrs[i] = buf;
// 	}
// 	// for (i = 0; i < m; i++) {
// 	// 	if (NULL == (frag_ptrs[i] = malloc(len))) {
// 	// 		printf("alloc error: Fail\n");
// 	// 		return -1;
// 	// 	}
// 	// }
// 	// Fill sources with random data
// 	for (i = 0; i < m; i++)
// 		for (j = 0; j < len; j++)
// 			frag_ptrs[i][j] = rand();
	

// 	//printf(" encode (m,k,p)=(%d,%d,%d) len=%d\n", m, k, p, len);

// 	// Pick an encode matrix. A Cauchy matrix is a good choice as even
// 	// large k are always invertable keeping the recovery rule simple.
// 	gf_gen_cauchy1_matrix(encode_matrix, m, k);
// 	// Initialize g_tbls from encode matrix
// 	ec_init_tables(k, p, &encode_matrix[k * k], g_tbls);

// 	int n=1000;
// 	clock_gettime(CLOCK_REALTIME, &time1);
// 	//ec_encode_data(len, k, p, g_tbls, frag_ptrs, &frag_ptrs[k]);
// 	for(i=0;i<n;i++){
// 		//clock_gettime(CLOCK_REALTIME, &time1);
// 		ec_encode_data(len, k, p, g_tbls, frag_ptrs, &frag_ptrs[k]);
// 		//clock_gettime(CLOCK_REALTIME, &time2);
// 		//encode_time=(time2.tv_sec-time1.tv_sec)*1000000000+(time2.tv_nsec-time1.tv_nsec);
// 		//encode_time_arrs[i]=encode_time;
// 		//printf("encode_time is %ld   ", encode_time);
// 		//printf("encode throughput(MB/s): %0.10f\n", (((double)len*m)/1024/1024)/((double)encode_time/1000/1000/1000));
// 	}
// 	clock_gettime(CLOCK_REALTIME, &time2);
// 	encode_time=(time2.tv_sec-time1.tv_sec)*1000000000+(time2.tv_nsec-time1.tv_nsec);
// 	encode_time/=n;
// 	//printf("encode_time is %ld \n",encode_time);
// 	printf("encode throughput(MB/s): %d\n", (int)((((double)len*m)/1024/1024)/((double)encode_time/1000/1000/1000)));
// 	//printf("%ld\n", encode_time);
// 	return 0;
// }

/*
 * Generate decode matrix from encode matrix and erasure list
 *
 */

static int gf_gen_decode_matrix_simple(u8 * encode_matrix,
				       u8 * decode_matrix,
				       u8 * invert_matrix,
				       u8 * temp_matrix,
				       u8 * decode_index, u8 * frag_err_list, int nerrs, int k,
				       int m)
{
	int i, j, p, r;
	int nsrcerrs = 0;
	u8 s, *b = temp_matrix;
	u8 frag_in_err[MMAX];

	memset(frag_in_err, 0, sizeof(frag_in_err));

	// Order the fragments in erasure for easier sorting
	for (i = 0; i < nerrs; i++) {
		if (frag_err_list[i] < k)
			nsrcerrs++;
		frag_in_err[frag_err_list[i]] = 1;
	}

	// Construct b (matrix that encoded remaining frags) by removing erased rows
	for (i = 0, r = 0; i < k; i++, r++) {
		while (frag_in_err[r])
			r++;
		for (j = 0; j < k; j++)
			b[k * i + j] = encode_matrix[k * r + j];
		decode_index[i] = r;
	}

	// Invert matrix to get recovery matrix
	if (gf_invert_matrix(b, invert_matrix, k) < 0)
		return -1;

	// Get decode matrix with only wanted recovery rows
	for (i = 0; i < nerrs; i++) {
		if (frag_err_list[i] < k)	// A src err
			for (j = 0; j < k; j++)
				decode_matrix[k * i + j] =
				    invert_matrix[k * frag_err_list[i] + j];
	}

	// For non-src (parity) erasures need to multiply encode matrix * invert
	for (p = 0; p < nerrs; p++) {
		if (frag_err_list[p] >= k) {	// A parity err
			for (i = 0; i < k; i++) {
				s = 0;
				for (j = 0; j < k; j++)
					s ^= gf_mul(invert_matrix[j * k + i],
						    encode_matrix[k * frag_err_list[p] + j]);
				decode_matrix[k * p + i] = s;
			}
		}
	}
	return 0;
}
