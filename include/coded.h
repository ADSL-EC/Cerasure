//
// Created by nty on 2022/7/8.
//
#ifndef ISA_L_XOR_MATCH_CODE_H
#define ISA_L_XOR_MATCH_CODE_H

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "raid.h"
#include "cauchy.h"
#include "galois.h"
//#include "types.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

typedef unsigned char u8;

#include"coded.h"
extern int mceil(int n,int m);

//-----------------decode---------------------------------------------------------------------------------------------------/
extern int invert_bitmatrix(int *mat, int *inv, int rows);

extern u8 **set_up_ptrs_for_opertaions(int k, int p, int *eraseds, char **data_ptrs, char **coding_ptrs,int *recover);

extern int set_up_ids_for_operations(int k, int p, int *eraseds, int *row_ids, int *ind_to_row,int *recover);

extern int *eraseds_to_erased(int k, int p, int *eraseds);

extern int count_bitmatrix_one_number(int k,int p,int w,int *bitmatrix,int begin,int end);

extern int*** bitmatrix_to_schedule_split(int k,int p, int w,int n,int *bitmarix);
extern int max_var(int *slp);

extern int rows_of_schedule(int **schedule);


extern void encode_base_deforestation(int k, int p, int w, int **schedule, u8 **src, u8 **dest,int len, int packetsize);

extern void encode_base_once_deforestation(int k, int p, int w, int **schedule,u8 **copy_ptr,int packetsize);

extern void encode_best(int k, int p, int w, int **schedule, u8 **src, u8 **dest,int len, int packetsize);
extern void  encode_best_spilt(int k, int p, int w, int ***operations, u8 **src, u8 **dest,int n,int len, int packetsize);

typedef unsigned char u8;

typedef struct {
	int k;
	int p;
	int w;
	int **schedule;
	u8 **src;
	u8 **dest;
	int len; 
	int packetsize;
	int value;
	int begin;
	int n;
}ptInfo;



#endif


