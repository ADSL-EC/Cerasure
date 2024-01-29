#pragma once

#ifdef __cplusplus
extern "C" {
#endif

extern int *cauchy_original_coding_matrix(int k, int m, int w);
extern int *cauchy_xy_coding_matrix(int k, int m, int w, int *x, int *y);
extern void cauchy_improve_coding_matrix(int k, int m, int w, int *matrix);
extern int *cauchy_good_general_coding_matrix(int k, int m, int w);
extern int cauchy_n_ones(int n, int w);

#ifdef __cplusplus
}
#endif