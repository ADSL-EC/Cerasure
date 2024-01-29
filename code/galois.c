//
// Created by nty on 2022/6/2.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include "galois.h"

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))
#define MAX_GF_INSTANCES 64
gf_t *gfp_array[MAX_GF_INSTANCES] = { 0 };
int  gfp_is_composite[MAX_GF_INSTANCES] = { 0 };

gf_t *galois_get_field_ptr(int w)
{
    if (gfp_array[w] != NULL) {
        return gfp_array[w];
    }

    return NULL;
}

gf_t* galois_init_field(int w,
                        int mult_type,
                        int region_type,
                        int divide_type,
                        uint64_t prim_poly,
                        int arg1,
                        int arg2)
{
    int scratch_size;
    void *scratch_memory;
    gf_t *gfp;

    if (w <= 0 || w > 32) {
        fprintf(stderr, "ERROR -- cannot init default Galois field for w=%d\n", w);
        assert(0);
    }

    gfp = (gf_t *) malloc(sizeof(gf_t));
    if (!gfp) {
        fprintf(stderr, "ERROR -- cannot allocate memory for Galois field w=%d\n", w);
        assert(0);
    }

    scratch_size = gf_scratch_size(w, mult_type, region_type, divide_type, arg1, arg2);
    if (!scratch_size) {
        fprintf(stderr, "ERROR -- cannot get scratch size for base field w=%d\n", w);
        assert(0);
    }

    scratch_memory = malloc(scratch_size);
    if (!scratch_memory) {
        fprintf(stderr, "ERROR -- cannot get scratch memory for base field w=%d\n", w);
        assert(0);
    }

    if(!gf_init_hard(gfp,
                     w,
                     mult_type,
                     region_type,
                     divide_type,
                     prim_poly,
                     arg1,
                     arg2,
                     NULL,
                     scratch_memory))
    {
        fprintf(stderr, "ERROR -- cannot init default Galois field for w=%d\n", w);
        assert(0);
    }

    gfp_is_composite[w] = 0;
    return gfp;
}

gf_t* galois_init_composite_field(int w,
                                  int region_type,
                                  int divide_type,
                                  int degree,
                                  gf_t* base_gf)
{
    int scratch_size;
    void *scratch_memory;
    gf_t *gfp;

    if (w <= 0 || w > 32) {
        fprintf(stderr, "ERROR -- cannot init composite field for w=%d\n", w);
        assert(0);
    }

    gfp = (gf_t *) malloc(sizeof(gf_t));
    if (!gfp) {
        fprintf(stderr, "ERROR -- cannot allocate memory for Galois field w=%d\n", w);
        assert(0);
    }

    scratch_size = gf_scratch_size(w, GF_MULT_COMPOSITE, region_type, divide_type, degree, 0);
    if (!scratch_size) {
        fprintf(stderr, "ERROR -- cannot get scratch size for composite field w=%d\n", w);
        assert(0);
    }

    scratch_memory = malloc(scratch_size);
    if (!scratch_memory) {
        fprintf(stderr, "ERROR -- cannot get scratch memory for composite field w=%d\n", w);
        assert(0);
    }

    if(!gf_init_hard(gfp,
                     w,
                     GF_MULT_COMPOSITE,
                     region_type,
                     divide_type,
                     0,
                     degree,
                     0,
                     base_gf,
                     scratch_memory))
    {
        fprintf(stderr, "ERROR -- cannot init default composite field for w=%d\n", w);
        assert(0);
    }
    gfp_is_composite[w] = 1;
    return gfp;
}

int galois_init_default_field(int w)
{
    if (gfp_array[w] == NULL) {
        gfp_array[w] = (gf_t*)malloc(sizeof(gf_t));
        if(gfp_array[w] == NULL)
            return ENOMEM;
        if (!gf_init_easy(gfp_array[w], w))
            return EINVAL;
    }
    return 0;
}

int galois_uninit_field(int w)
{
    int ret = 0;
    if (gfp_array[w] != NULL) {
        int recursive = 1;
        ret = gf_free(gfp_array[w], recursive);
        free(gfp_array[w]);
        gfp_array[w] = NULL;
    }
    return ret;
}

 void galois_init(int w)
{
    if (w <= 0 || w > 32) {
        fprintf(stderr, "ERROR -- cannot init default Galois field for w=%d\n", w);
        assert(0);
    }

    switch (galois_init_default_field(w)) {
        case ENOMEM:
            fprintf(stderr, "ERROR -- cannot allocate memory for Galois field w=%d\n", w);
            assert(0);
            break;
        case EINVAL:
            fprintf(stderr, "ERROR -- cannot init default Galois field for w=%d\n", w);
            assert(0);
            break;
    }
}


static int is_valid_gf(gf_t *gf, int w)
{
    // TODO: I assume we may eventually
    // want to do w=64 and 128, so w
    // will be needed to perform this check
    (void)w;

    if (gf == NULL) {
        return 0;
    }
    if (gf->multiply.w32 == NULL) {
        return 0;
    }
    if (gf->multiply_region.w32 == NULL) {
        return 0;
    }
    if (gf->divide.w32 == NULL) {
        return 0;
    }
    if (gf->inverse.w32 == NULL) {
        return 0;
    }
    if (gf->extract_word.w32 == NULL) {
        return 0;
    }

    return 1;
}

void galois_change_technique(gf_t *gf, int w)
{
    if (w <= 0 || w > 32) {
        fprintf(stderr, "ERROR -- cannot support Galois field for w=%d\n", w);
        assert(0);
    }

    if (!is_valid_gf(gf, w)) {
        fprintf(stderr, "ERROR -- overriding with invalid Galois field for w=%d\n", w);
        assert(0);
    }

    if (gfp_array[w] != NULL) {
        gf_free(gfp_array[w], gfp_is_composite[w]);
    }

    gfp_array[w] = gf;
}

int galois_single_multiply(int x, int y, int w)
{
    if (x == 0 || y == 0) return 0;

    if (gfp_array[w] == NULL) {
        galois_init(w);
    }

    if (w <= 32) {
        return gfp_array[w]->multiply.w32(gfp_array[w], x, y);
    } else {
        fprintf(stderr, "ERROR -- Galois field not implemented for w=%d\n", w);
        return 0;
    }
}

int galois_single_divide(int x, int y, int w)
{
    if (x == 0) return 0;
    if (y == 0) return -1;

    if (gfp_array[w] == NULL) {
        galois_init(w);
    }

    if (w <= 32) {
        return gfp_array[w]->divide.w32(gfp_array[w], x, y);
    } else {
        fprintf(stderr, "ERROR -- Galois field not implemented for w=%d\n", w);
        return 0;
    }
}

void galois_w08_region_multiply(char *region,      /* Region to multiply */
                                int multby,       /* Number to multiply by */
                                int nbytes,        /* Number of bytes in region */
                                char *r2,          /* If r2 != NULL, products go here */
                                int add)
{
    if (gfp_array[8] == NULL) {
        galois_init(8);
    }
    gfp_array[8]->multiply_region.w32(gfp_array[8], region, r2, multby, nbytes, add);
}

void galois_w16_region_multiply(char *region,      /* Region to multiply */
                                int multby,       /* Number to multiply by */
                                int nbytes,        /* Number of bytes in region */
                                char *r2,          /* If r2 != NULL, products go here */
                                int add)
{
    if (gfp_array[16] == NULL) {
        galois_init(16);
    }
    gfp_array[16]->multiply_region.w32(gfp_array[16], region, r2, multby, nbytes, add);
}


void galois_w32_region_multiply(char *region,      /* Region to multiply */
                                int multby,       /* Number to multiply by */
                                int nbytes,        /* Number of bytes in region */
                                char *r2,          /* If r2 != NULL, products go here */
                                int add)
{
    if (gfp_array[32] == NULL) {
        galois_init(32);
    }
    gfp_array[32]->multiply_region.w32(gfp_array[32], region, r2, multby, nbytes, add);
}

void galois_w8_region_xor(void *src, void *dest, int nbytes)
{
    if (gfp_array[8] == NULL) {
        galois_init(8);
    }
    gfp_array[8]->multiply_region.w32(gfp_array[32], src, dest, 1, nbytes, 1);
}

void galois_w16_region_xor(void *src, void *dest, int nbytes)
{
    if (gfp_array[16] == NULL) {
        galois_init(16);
    }
    gfp_array[16]->multiply_region.w32(gfp_array[16], src, dest, 1, nbytes, 1);
}

void galois_w32_region_xor(void *src, void *dest, int nbytes)
{
    if (gfp_array[32] == NULL) {
        galois_init(32);
    }
    gfp_array[32]->multiply_region.w32(gfp_array[32], src, dest, 1, nbytes, 1);
}

void galois_region_xor(char *src, char *dest, int nbytes)
{
    if (nbytes >= 16) {
        galois_w32_region_xor(src, dest, nbytes);
    } else {
        int i = 0;
        for (i = 0; i < nbytes; i++) {
            *dest ^= *src;
            dest++;
            src++;
        }
    }
}

int galois_inverse(int y, int w)
{
    if (y == 0) return -1;
    return galois_single_divide(1, y, w);
}

void print_bitmatrix(int *m, int rows, int cols, int w)
{
    int i, j;

    for (i = 0; i < rows; i++) {
        if (i != 0 && i%w == 0) printf("\n");
        for (j = 0; j < cols; j++) {
            if (j != 0 && j%w == 0) printf(" ");
            printf("%d,", m[i*cols+j]);
        }
        printf("\n");
    }
}

int *matrix_to_bitmatrix(int k, int m, int w, int *matrix)
{
    int *bitmatrix;
    int rowelts, rowindex, colindex, elt, i, j, l, x;

    if (matrix == NULL) { return NULL; }

    bitmatrix = talloc(int, k*m*w*w);
    if (!bitmatrix) return NULL;

    rowelts = k * w;
    rowindex = 0;

    for (i = 0; i < m; i++) {
        colindex = rowindex;
        for (j = 0; j < k; j++) {
            elt = matrix[i*k+j];
            for (x = 0; x < w; x++) {
                for (l = 0; l < w; l++) {
                    bitmatrix[colindex+x+l*rowelts] = ((elt & (1 << l)) ? 1 : 0);
                }
                elt = galois_single_multiply(elt, 2, w);
            }
            colindex += w;
        }
        rowindex += rowelts * w;
    }
    return bitmatrix;
}

