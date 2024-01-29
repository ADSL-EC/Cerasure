#include "vandermonde.h"
#include "galois.h"
#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

#define MAX_CHECK 63		/* Size is limited by using uint64_t to represent subsets */
#define M_MAX 0x3F
#define K_MAX 0x3F
#define ROWS M_MAX
#define COLS K_MAX

int **bitmatrix_to_schedule(int k, int p, int w, int *bitmatrix)
{
	int **operations;
	int oper_index, oper_count;
	int i, j;

	operations=malloc(sizeof(int *) * (p*w+1));
	if (!operations) 
		return NULL;
  	oper_index=0;
	oper_count=0;

	for(i=0;i<p*w;i++){
		for(j=0;j<k*w;j++){
			if(bitmatrix[i*k*w+j]){
				oper_count++;
			}
		}
		operations[i]=malloc(sizeof(int)*(2*oper_count+3));
		operations[i][0]=oper_count;
		oper_index++;
		operations[i][oper_index]=k+i/w;
		oper_index++;
		operations[i][oper_index]=i%w;
		oper_index++;
		for(j=0;j<k*w;j++){
			if(bitmatrix[i*k*w+j]){
				operations[i][oper_index]=j/w;
				oper_index++;
				operations[i][oper_index]=j%w;
				oper_index++;
			}
		}

		oper_index=0;
		oper_count=0;
	}
	operations[i]=malloc(sizeof(int));
	operations[p*w][0]=-1;

	return operations;	
}

int gf_invert_matrix_w(unsigned int *in_mat, unsigned int *out_mat, const int n,int w)
{
	int i, j, k;
	unsigned int temp;
	// Set out_mat[] to the identity matrix
	for (i = 0; i < n * n; i++)	// memset(out_mat, 0, n*n)
		out_mat[i] = 0;

	for (i = 0; i < n; i++)
		out_mat[i * n + i] = 1;
	// Inverse
	for (i = 0; i < n; i++) {
		// Check for 0 in pivot element
		if (in_mat[i * n + i] == 0) {
			// Find a row with non-zero in current column and swap
			for (j = i + 1; j < n; j++)
				if (in_mat[j * n + i])
					break;
			if (j == n)	// Couldn't find means it's singular
				return -1;
			for (k = 0; k < n; k++) {	// Swap rows i,j
				temp = in_mat[i * n + k];
				in_mat[i * n + k] = in_mat[j * n + k];
				in_mat[j * n + k] = temp;

				temp = out_mat[i * n + k];
				out_mat[i * n + k] = out_mat[j * n + k];
				out_mat[j * n + k] = temp;
			}
		}
		temp = galois_inverse(in_mat[i * n + i],w);	// 1/pivot
        
		for (j = 0; j < n; j++) {	// Scale row i by 1/pivot
			in_mat[i * n + j] = galois_single_multiply(in_mat[i * n + j], temp,w);
			out_mat[i * n + j] = galois_single_multiply(out_mat[i * n + j], temp,w);
		}

		for (j = 0; j < n; j++) {
			if (j == i)
				continue;

			temp = in_mat[j * n + i];
			for (k = 0; k < n; k++) {
				out_mat[j * n + k] ^= galois_single_multiply(temp, out_mat[i * n + k],w);
				in_mat[j * n + k] ^= galois_single_multiply(temp, in_mat[i * n + k],w);
			}
		}
	}
	return 0;
}

static inline int min(int a, int b)
{
	if (a <= b)
		return a;
	else
		return b;
}

void gen_sub_matrix(unsigned int *out_matrix, unsigned int *in_matrix, int rows,
		    int cols, unsigned int row_indicator, unsigned int col_indicator)
{
	int i, j, r, s,dim;
    int p=rows-cols;
    dim=0;

    for(i = 0; i < cols; i++)
    {
        if (!(col_indicator & ((uint64_t) 1 << i)))
        {
            for (j = 0 ; j < cols; j++) {
                out_matrix[dim * cols + j] = in_matrix[cols * i + j];
            }
            dim++;
        }
    }

    for(i=0;i<p;i++)
    {
        if (row_indicator & ((uint64_t) 1 << i))
        {
            for (j = 0 ; j < cols; j++) {
			    out_matrix[dim * cols + j] = in_matrix[cols * (i+cols) + j];
            }
            dim++;
        }
    }

}

/* Gosper's Hack */
uint64_t next_subset(uint64_t * subset, uint64_t element_count, uint64_t subsize)
{
	uint64_t tmp1 = *subset & -*subset;
	uint64_t tmp2 = *subset + tmp1;
	*subset = (((*subset ^ tmp2) >> 2) / tmp1) | tmp2;
	if (*subset & (((uint64_t) 1 << element_count))) {
		/* Overflow on last subset */
		*subset = ((uint64_t) 1 << subsize) - 1;
		return 1;
	}

	return 0;
}



void gen_sub_matrix_w(unsigned int *out_matrix, int dim, unsigned int *in_matrix, int rows,
		    int cols, uint64_t row_indicator, uint64_t col_indicator)
{
	int i, j, r, s;
	for (i = 0, r = 0; i < rows; i++) {
		if (!(row_indicator & ((uint64_t) 1 << i)))
			continue;
		for (j = 0, s = 0; j < cols; j++) {
			if (!(col_indicator & ((uint64_t) 1 << j)))
				continue;
			out_matrix[dim * r + s] = in_matrix[cols * i + j];
			s++;
		}
		r++;
	}
}

int are_submatrices_singular_w(unsigned int *vmatrix, long long int rows, long long int cols,int w)
{
	unsigned int matrix[COLS * COLS];
	unsigned int invert_matrix[COLS * COLS];
	uint64_t row_indicator, col_indicator, subset_init, subsize;

	/* Check all square subsize x subsize submatrices of the rows x cols
	 * vmatrix for singularity*/
    subsize=min(rows, cols);
    subset_init = (1 << subsize) - 1;
    int number=0;
    int i=0;
    int j=0;
    col_indicator = subset_init;
    do {
        row_indicator = subset_init;
        do {
            gen_sub_matrix_w(matrix, subsize, vmatrix, rows,
                        cols, row_indicator, col_indicator);

            if (gf_invert_matrix_w(matrix, invert_matrix, subsize,w)==-1)
                return 1;
        } while (next_subset(&row_indicator, rows, subsize) == 0);
    } while (next_subset(&col_indicator, cols, subsize) == 0);

	return 0;
}

//检查矩阵是否满足MDS性质
int check_coding_matrix_w(int k, int m, int w,unsigned int* encode_matrix)
{

    unsigned int *matrix;
    int i,j;
    matrix=malloc(k * (k+m) * sizeof(unsigned int));
	memset(matrix, 0, k * (k+m) );
	for (i = 0; i < k; i++)
    {
        matrix[k * i + i] = 1;
    }
	for (i = k; i < (m+k); i++) {
		for (j = 0; j < k; j++) {
			matrix[k * i + j] = encode_matrix[(i-k)*k+j];
		}
	}

    if(are_submatrices_singular_w(matrix,k+m,k,w))
    {
        return 1;
    }
    return 0;
}

unsigned int *vandermonde_sub_coding_matrix(int k, int m, int w){
    int i,j,tmp,elem;
    unsigned int *matrix;
    matrix=malloc(sizeof(unsigned int)*k*m);
    elem=1;
    for(i=0;i<m;i++)
    {
        tmp=1;
        for(j=0;j<k;j++)
        {
            matrix[i*k+j]=tmp;
            tmp=galois_single_multiply(tmp,elem,w);
        }
        elem=galois_single_multiply(elem,2,w);
    }
    return matrix;
}

int count_one_numbers(int k, int m, int w, unsigned int* bitmatrix)
{
    int count=0;
    int i,j;
    for(i=0;i<m*w;i++)
    {
        for(j=0;j<k*w;j++)
        {
            if(bitmatrix[i*k*w+j])
            {
                count++;
            }
        }
    }
    return count;
}


void vandermonde_R(int n,int w,int a)
{
    int i;
    int tmp=1;
    for(i=0;i<n-1;i++)
    {
        printf("%d,",tmp);
        tmp=galois_single_multiply(tmp,a,w);
    }
    printf("%d",tmp);
}

void vandermonde_all(w)
{
    int i;
    int n = 1;
    for(i=0;i<w;i++)
    {
        n=n*2;
    } 
    printf("{");
    for(i=1;i<n-1;i++)
    {
        printf("{");
        vandermonde_R(n,w,i);
        printf("},");
    }
    printf("{");
    vandermonde_R(n,w,i);
    printf("}");
    printf("}");
}

void count_one_R(int n,int w,int a)
{
    int i;
    int tmp=a;
    int count=0;
    int* encode_bitmatrix;
    for(i=0;i<n-1;i++)
    {
        encode_bitmatrix=matrix_to_bitmatrix(1, 1, w, &R_w4[a-1][i]);
        count+=count_one_numbers(1,1,w,encode_bitmatrix);
        printf("%d,",count);
    }
    encode_bitmatrix=matrix_to_bitmatrix(1, 1, w, &R_w4[a-1][i]);
    count+=count_one_numbers(1,1,w,encode_bitmatrix);
    printf("%d",count);
}

void count_one_R_all(int w)
{
    int i;
    int n = 1;
    for(i=0;i<w;i++)
    {
        n=n*2;
    } 
    printf("{");
    for(i=1;i<n-1;i++)
    {
        printf("{");
        count_one_R(n,w,i);
        printf("},");
    }
    printf("{");
    count_one_R(n,w,i);
    printf("}");
    printf("}");
}
int comp(const R_one* a,const R_one* b)//用来做比较的函数。
{
    return (*a).num-(*b).num;
}

unsigned int *vandermonde_sub_coding_matrix_greedy(int k, int m, int w,int count_cauchy,int* flag){
    int loop=0;
    int i,j,tmp,elem,index;
    int n = 1;
    for(i=0;i<w;i++)
    {
        n=n*2;
    } 
    R_one order[n];
    unsigned int *matrix;
    int* encode_bitmatrix;
    int parents[n];
    int (*p)[n];
    int (*p_ones)[n];
    index_order indexs[m];
    int count=0;
    matrix=malloc(sizeof(unsigned int)*k*m);
    for(i=0;i<n;i++)
    {
        parents[i]=0;
    }
    if(w==8)
    {
        p=R_w8;
        p_ones=R_ones_w8;
    }
    else if(w==4)
    {
        p=R_w4;
        p_ones=R_ones_w4;
    }
    for(i=0;i<n-1;i++)
    {
        order[i].index = i;
        order[i].num = p_ones[i][k-1];
    }
    qsort(order,n-1,sizeof(R_one),comp);
    for(i=0;i<m;i++)
    {
        indexs[i].index=order[i].index;
        indexs[i].order_index=i;
        parents[indexs[i].index]=1;
    }
    for(i=0;i<m;i++)
    {
        for(j=0;j<k;j++)
        {
            matrix[i*k+j]=p[indexs[i].index][j];
        }
    }

    encode_bitmatrix=matrix_to_bitmatrix(k, m, w, matrix);
    count=count_one_numbers(k,m,w,encode_bitmatrix);
    int tmp_begin_index=-1;
    int tmp_end_index=-1;
    int tmp_num;
    int tmp_min;
    while(count < count_cauchy && check_coding_matrix_w(k,m,w,matrix))
    {
        tmp_min=1<<30;
        for(i=0;i<m;i++)
        {
            j= 1 + indexs[i].order_index;
            while( j < n-1 && parents[order[j].index]==1)
            {
                j++;
            }
            if(j < n-1)
            {
                tmp_num=p_ones[order[j].index][k-1]-p_ones[indexs[i].index][k-1];
                if(tmp_min >= tmp_num)
                {
                    tmp_min=tmp_num;
                    tmp_begin_index=i;
                    tmp_end_index=j;
                }
            }
        }
        parents[indexs[tmp_begin_index].index]=0;
        parents[order[tmp_end_index].index]=1;
        indexs[tmp_begin_index].order_index = tmp_end_index;
        indexs[tmp_begin_index].index=order[tmp_end_index].index;
        for(i=0;i<m;i++)
        {
            for(j=0;j<k;j++)
            {
                matrix[i*k+j]=p[indexs[i].index][j];
            }
        }
        encode_bitmatrix=matrix_to_bitmatrix(k, m, w, matrix);
        count=count_one_numbers(k,m,w,encode_bitmatrix);
        loop++;
        // if(loop>nums)
        // {
        //     break;
        // }
    }
    if(count >= count_cauchy )
    {
        *flag = 0;
        return matrix;
    }
    *flag=1;
    return matrix;
}

unsigned int* v_serach(int k, int m, int w)
{
    unsigned int* encode_matrix_sub;
    unsigned int* encode_matrix_cauchy;
    int flag=0;
	int* encode_bitmatrix_cauchy;
    int count_sub=0;
    int count_cauchy=0;
    encode_matrix_cauchy= cauchy_good_general_coding_matrix(k,m,w);
    encode_bitmatrix_cauchy=matrix_to_bitmatrix(k, m, w, encode_matrix_cauchy);
    count_cauchy=count_one_numbers(k,m,w,encode_bitmatrix_cauchy);
    encode_matrix_sub=vandermonde_sub_coding_matrix_greedy(k, m, w,count_cauchy,&flag);
    //printf("%d\n",flag);
    if(flag==1)
    {
        return encode_matrix_sub;
    }
    return encode_matrix_cauchy;
}


int comp_desc(const matrix_num_index* a,const matrix_num_index* b)//用来做比较的函数。
{
    return (*b).num-(*a).num;
}

int comp_asce(const matrix_num_index* a,const matrix_num_index* b)//用来做比较的函数。
{
    return (*a).num-(*b).num;
}

