#include"coded.h"
int mceil(int n,int m)
{
	if(n%m==0)
	{
		return n/m;
	}
	return n/m+1;
}

//-----------------decode---------------------------------------------------------------------------------------------------/
int invert_bitmatrix(int *mat, int *inv, int rows)
{
  int cols, i, j, k;
  int tmp;
 
  cols = rows;

  k = 0;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      inv[k] = (i == j) ? 1 : 0;
      k++;
    }
  }

  /* First -- convert into upper triangular */

  for (i = 0; i < cols; i++) {

    /* Swap rows if we have a zero i,i element.  If we can't swap, then the 
       matrix was not invertible */

    if ((mat[i*cols+i]) == 0) { 
      for (j = i+1; j < rows && (mat[j*cols+i]) == 0; j++) ;
      if (j == rows) return -1;
      for (k = 0; k < cols; k++) {
        tmp = mat[i*cols+k]; mat[i*cols+k] = mat[j*cols+k]; mat[j*cols+k] = tmp;
        tmp = inv[i*cols+k]; inv[i*cols+k] = inv[j*cols+k]; inv[j*cols+k] = tmp;
      }
    }
 
    /* Now for each j>i, add A_ji*Ai to Aj */
    for (j = i+1; j != rows; j++) {
      if (mat[j*cols+i] != 0) {
        for (k = 0; k < cols; k++) {
          mat[j*cols+k] ^= mat[i*cols+k]; 
          inv[j*cols+k] ^= inv[i*cols+k];
        }
      }
    }
  }

  /* Now the matrix is upper triangular.  Start at the top and multiply down */

  for (i = rows-1; i >= 0; i--) {
    for (j = 0; j < i; j++) {
      if (mat[j*cols+i]) {
        for (k = 0; k < cols; k++) {
          mat[j*cols+k] ^= mat[i*cols+k]; 
          inv[j*cols+k] ^= inv[i*cols+k];
        }
      }
    }
  } 
  return 0;
}

u8 **set_up_ptrs_for_opertaions(int k, int p, int *eraseds, char **data_ptrs, char **coding_ptrs,int *recover)
{
  int ddf, cdf;
  int *erased;
  int *visted;
  char **ptrs;
  int i, j, x;

  ddf = 0;
  cdf = 0;
  for (i = 0; eraseds[i] != -1; i++) {
    if (eraseds[i] < k) ddf++; else cdf++;
  }
  
  erased = eraseds_to_erased(k, p, eraseds);
  if (erased == NULL) return NULL;

  /* Set up ptrs.  It will be as follows:

       - If data drive i has not failed, then ptrs[i] = data_ptrs[i].
       - If data drive i has failed, then ptrs[i] = coding_ptrs[j], where j is the 
            lowest unused non-failed coding drive.
       - Elements k to k+ddf-1 are data_ptrs[] of the failed data drives.
       - Elements k+ddf to k+ddf+cdf-1 are coding_ptrs[] of the failed data drives.

       The array row_ids contains the ids of ptrs.
       The array ind_to_row_ids contains the row_id of drive i.
  
       However, we're going to set row_ids and ind_to_row in a different procedure.
   */      
	ptrs = talloc(u8 *, k+p);
	if (!ptrs) {
		free(erased);
		return NULL;
	}
	/* 重新实现指针的赋值*/
	visted = talloc(int, k+p);
	for(i=0;i<k+p;i++)
	{
		visted[i]=0;
	}
	j = 0;

	for(i=0;i<k;i++)
	{
		if(recover[i]<k)
		{
			ptrs[recover[i]] = data_ptrs[recover[i]];
			visted[recover[i]]=recover[i];
		}
		else
		{
			while(visted[j])j++;
			ptrs[j]=coding_ptrs[recover[i]-k];
			visted[j]=j;
			j++;
		}
	}
	x=k;
	for (i = 0; i < k+p; i++) {
		if (erased[i]==1) 
		{
			if(erased[i]<k)
			{
				ptrs[x] = data_ptrs[i];
			}
			else
			{
				ptrs[x] = coding_ptrs[i-k];
			}
			x++;
		}
	}
  
  free(erased);
  return ptrs;
}

int set_up_ids_for_operations(int k, int p, int *eraseds, int *row_ids, int *ind_to_row,int *recover)
{
	int ddf, cdf;
	int *erased;
	int i, j, x;

	ddf = 0;
	cdf = 0;
	for (i = 0; eraseds[i] != -1; i++) {
		if (eraseds[i] < k) ddf++; else cdf++;
	}
	
	erased = eraseds_to_erased(k, p, eraseds);
	if (erased == NULL) return -1;
	for(i=0;i<k+p;i++)
	{
		row_ids[i]=0;
	}
	/* See set_up_ptrs_for_opertaions for how these are set */

	j = 0;
	for(i=0;i<k;i++)
	{
		if(recover[i] < k)
		{
			row_ids[recover[i]]=recover[i];
			ind_to_row[recover[i]] = recover[i];
		}
		else
		{
			while(row_ids[j]) j++;
			row_ids[j]=recover[i];
			ind_to_row[j] = recover[i];
			j++;
		}
	}
	x=k;
	for (i = 0; i < k+p; i++) {
		if (erased[i]==1) 
		{
			row_ids[x] = i;
			ind_to_row[x] = i;
			x++;
		}
	}
	
	free(erased);
	return 0;
}

int *eraseds_to_erased(int k, int p, int *eraseds)
{
  int td;
  int t_non_erased;
  int *erased;
  int i;

  td = k+p;
  erased = talloc(int, td);
  if (erased == NULL) return NULL;
  t_non_erased = td;

  for (i = 0; i < td; i++) erased[i] = 0;

  for (i = 0; eraseds[i] != -1; i++) {
    if (erased[eraseds[i]] == 0) {
      erased[eraseds[i]] = 1;
      t_non_erased--;
      if (t_non_erased < k) {
        free(erased);
        return NULL;
      }
    }
  }
  return erased;
}

int count_bitmatrix_one_number(int k,int p,int w,int *bitmatrix,int begin,int end)
{
	int one_num=0;
	int i=0;
	for(i=0;i<k * w * (end-begin) * w;i++)
	{
		if(bitmatrix[i+begin*w*w*k]==1)
		{
			one_num++;
		}
	}
	return one_num;
}

int*** bitmatrix_to_schedule_split(int k,int p, int w,int n,int *bitmarix)
{
	int ***schedule;
	int i;
	schedule=(int **)malloc(sizeof(int**)*(n+1));
	for(i=0;i<n;i++)
	{
		int offset=i*(k*w/n);
		schedule[i]=(int *)malloc(sizeof(int *) * (p*w+1));
		int j,a,count;
		j=0,a=0;
		for(j=0;j<p*w;j++)
		{
			schedule[i][j]=(int )malloc(sizeof(int *) * (2*k*w+3));
			count=0;
			schedule[i][j][1] = j/w+k;
			schedule[i][j][2] = j%w;
			for(a=0;a<k*w/n;a++)
			{
				if(bitmarix[j*k*w+a+offset]==1)
				{
					schedule[i][j][count*2+3]= (a+offset)/w;
					schedule[i][j][count*2+4]= (a+offset)%w;
					count++;
				}
			}
			schedule[i][j][0]=count;
		}
		schedule[i][j]=malloc(sizeof(int *) * (2*k*w+3));
		schedule[i][j][0]=-1;
	}
	return schedule;
}

int max_var(int *schedule)
{
    int i=0;
    int value=0;
    while(schedule[i]!=-1)
    {
        if(schedule[i] > value)
        {
            value=schedule[i];
        }
        i++;
    }
    return value;
}

int rows_of_schedule(int **schedule)
{
	int index=0;
    while(schedule[index][0]!=-1)
    {
        index++;
    }
    return index;
}


void encode_base_deforestation(int k, int p, int w, int **schedule, u8 **src, u8 **dest,int len, int packetsize)
{
	int i,j;
    u8 *copy_ptr[p*w+1];
	for(i=0;i<k;i++)
	{
		copy_ptr[i] = src[i];
	}
	for(i=0;i<p;i++)
	{
		copy_ptr[i+k] = dest[i];
	}

	for(i=0;i<len;i += packetsize*w)
	{
		encode_base_once_deforestation(k, p, w, schedule, copy_ptr ,packetsize);
		for(j=0;j<k+p;j++)
		{
			copy_ptr[j] += (packetsize*w);
		}
	}
}

void encode_base_once_deforestation(int k, int p, int w, int **schedule,u8 **copy_ptr,int packetsize)
{
    void **data_buffs;
	void **parity_buffs;
	int j;
	data_buffs=malloc(sizeof(void *) * (2*k*w));
	parity_buffs=malloc(sizeof(void *));
	int index=0;

	while(schedule[index][0] >= 0){
        for(j=0;j<schedule[index][0];j++)
        {
            data_buffs[j] = copy_ptr[schedule[index][2*j+3]] + schedule[index][2*j+4] * packetsize ;
        }
        parity_buffs[0]=copy_ptr[schedule[index][1]] + schedule[index][2]*packetsize;
        xor_gen(schedule[index][0],packetsize, data_buffs, parity_buffs);
		index++;
	}
}

void encode_best(int k, int p, int w, int **schedule, u8 **src, u8 **dest,int len, int packetsize)
{
	int i,j;
    u8 *copy_ptr[k+p];

    long *add_offset;
    long *buf;
	if (posix_memalign(&buf, 32, 16*sizeof(long))) {
		return 1;
	}
	add_offset = buf;
    long add_value= w * packetsize;
	for(i=0;i<16;i++)
	{
		add_offset[i]= add_value;
	}
    
	for(i=0;i<k;i++)
	{
		copy_ptr[i] = src[i];
	}
	for(i=0;i<p;i++)
	{
		copy_ptr[i+k] = dest[i];
	}

	void **indexs;
	void **data_buffs[p*w];
	void **parity_buffs[p*w];
	for(i=0;i<p*w;i++)
	{
		posix_memalign(&buf, 32, k*w*sizeof(long));
		data_buffs[i]=buf;
		parity_buffs[i]=malloc(sizeof(void *));
	}
	posix_memalign(&buf, 32, p*w*k*w*sizeof(long));
	indexs=buf;
	int ptr_index=0;
	int index=0;
    while(schedule[index][0] >= 0){
        data_buffs[index]=&indexs[ptr_index];
        for(j=0;j<schedule[index][0];j++)
        {
            indexs[ptr_index++] = copy_ptr[schedule[index][2*j+3]] + schedule[index][2*j+4] * packetsize;
        }
        parity_buffs[index]=&indexs[ptr_index];
        indexs[ptr_index++] = copy_ptr[schedule[index][1]] + schedule[index][2]*packetsize;
		index++;
	}
	int loop=mceil(ptr_index,8)*64;
	for(i=0;i<len;i+=add_value)
	{
		for(j=0;j < index;j++)
		{
			xor_gen(schedule[j][0], packetsize, data_buffs[j] , parity_buffs[j]);
		}
		buff_finally(loop,add_offset,indexs);
	}
}

void  encode_best_spilt(int k, int p, int w, int ***operations, u8 **src, u8 **dest,int n,int len, int packetsize)
{
	int i,j,a;
    u8 *copy_ptr[k+p];

    long *add_offset;
    long *buf;
	if (posix_memalign(&buf, 32, 16*sizeof(long))) {
		return 1;
	}
	add_offset = buf;
    long add_value= w * packetsize;
	for(i=0;i<16;i++)
	{
		add_offset[i]= add_value;
	}
    
	for(i=0;i<k;i++)
	{
		copy_ptr[i] = src[i];
	}
	for(i=0;i<p;i++)
	{
		copy_ptr[i+k] = dest[i];
	}


	void **indexs;
	void **data_buffs[n][p*w];
	void **parity_buffs[n][p*w];
	int oper_group_indexs[n];

	for(i=0;i<n;i++)
	{
		for(j=0;j<p*w;j++)
		{
			data_buffs[i][j]=malloc(k*w*sizeof(long)/n);
			parity_buffs[i][j]=malloc(sizeof(long) * 2);
		}
	}
	posix_memalign(&buf, 32,p*w*k*w*sizeof(long));
	indexs=buf;
	int ptr_index=0;
	int index;
	for(i=0;i<n;i++)
	{
		index=0;
		while(operations[i][index][0] >= 0){
			data_buffs[i][index]=&indexs[ptr_index];
		
			for(j=0;j<operations[i][index][0];j++)
			{
				indexs[ptr_index++] = copy_ptr[operations[i][index][2*j+3]] + operations[i][index][2*j+4] * packetsize;
			}
			if(i==0)
			{
				parity_buffs[i][index]=&indexs[ptr_index];
				indexs[ptr_index++] = copy_ptr[operations[i][index][1]] + operations[i][index][2]*packetsize;
			}else
			{
				parity_buffs[i][index]=parity_buffs[0][index];
			}
			index++;
		}	
	}	
	int loop=mceil(ptr_index,8)*64;
	for(i=0;i<len;i+=add_value)
	{
		for(a=0;a < p*w;a++)
		{
			split_genf(operations[0][a][0], packetsize, data_buffs[0][a], parity_buffs[0][a]);
		}

		for(j=1;j<n;j++)
		{			
			for(a=0;a < p*w;a++)
			{
				split_gen(operations[j][a][0], packetsize, data_buffs[j][a], parity_buffs[j][a]);
			}
		}
		buff_finally(loop,add_offset,indexs);
	}

}


int *generate_decoding_bitmatrix(int k, int p, int w, int *bitmatrix, int *eraseds,int *recover)
{
  int i, j, x, drive, y, index, z;
  int *decoding_matrix, *inverse, *real_decoding_matrix;
  int *ptr;
  int *row_ids;
  int *ind_to_row;
  int ddf, cdf;
  int **schedule;
  int *b1, *b2;
 
 /* First, figure out the number of data drives that have failed, and the
    number of coding drives that have failed: ddf and cdf */

  ddf = 0;
  cdf = 0;
  for (i = 0; eraseds[i] != -1; i++) {
    if (eraseds[i] < k) ddf++; else cdf++;
  }
  
  row_ids = talloc(int, k+p);
  if (!row_ids) return NULL;
  ind_to_row = talloc(int, k+p);
  if (!ind_to_row) {
    free(row_ids);
    return NULL;
  }

  if (set_up_ids_for_operations(k, p, eraseds, row_ids, ind_to_row,recover) < 0) {
    free(row_ids);
    free(ind_to_row);
    return NULL;
  }
  /* Now, we're going to create one decoding matrix which is going to 
     decode everything with one call.  The hope is that the scheduler
     will do a good job.    This matrix has w*e rows, where e is the
     number of eraseds (ddf+cdf) */
  real_decoding_matrix = talloc(int, k*w*(cdf+ddf)*w);
  if (!real_decoding_matrix) {
    free(row_ids);
    free(ind_to_row);
    return NULL;
  }

  /* First, if any data drives have failed, then initialize the first
     ddf*w rows of the decoding matrix from the standard decoding
     matrix inversion */

  if (ddf > 0) {
    
    decoding_matrix = talloc(int, k*k*w*w);
    if (!decoding_matrix) {
      free(row_ids);
      free(ind_to_row);
      return NULL;
    }
    ptr = decoding_matrix;
    for (i = 0; i < k; i++) {
      if (row_ids[i] == i) {
        bzero(ptr, k*w*w*sizeof(int));
        for (x = 0; x < w; x++) {
          ptr[x+i*w+x*k*w] = 1;
        } 
      } else {
        memcpy(ptr, bitmatrix+k*w*w*(row_ids[i]-k), k*w*w*sizeof(int));
      }
      ptr += (k*w*w);
    }
    inverse = talloc(int, k*k*w*w);
    if (!inverse) {
      free(row_ids);
      free(ind_to_row);
      free(decoding_matrix);
      return NULL;
    }

    invert_bitmatrix(decoding_matrix, inverse, k*w);

    free(decoding_matrix);
    ptr = real_decoding_matrix;
    for (i = 0; i < ddf; i++) {
      memcpy(ptr, inverse+k*w*w*row_ids[k+i], sizeof(int)*k*w*w);
      ptr += (k*w*w);
    }
    free(inverse);
  } 

  /* Next, here comes the hard part.  For each coding node that needs
     to be decoded, you start by putting its rows of the distribution
     matrix into the decoding matrix.  If there were no failed data
     nodes, then you're done.  However, if there have been failed
     data nodes, then you need to modify the columns that correspond
     to the data nodes.  You do that by first zeroing them.  Then
     whereever there is a one in the distribution matrix, you XOR
     in the corresponding row from the failed data node's entry in
     the decoding matrix.  The whole process kind of makes my head
     spin, but it works.
   */

  for (x = 0; x < cdf; x++) {
    drive = row_ids[x+ddf+k]-k;
    ptr = real_decoding_matrix + k*w*w*(ddf+x);
    memcpy(ptr, bitmatrix+drive*k*w*w, sizeof(int)*k*w*w);

    for (i = 0; i < k; i++) {
      if (row_ids[i] != i) {
        for (j = 0; j < w; j++) {
          bzero(ptr+j*k*w+i*w, sizeof(int)*w);
        }
      }  
    }

    /* There's the yucky part */

    index = drive*k*w*w;
    for (i = 0; i < k; i++) {
      if (row_ids[i] != i) {
        b1 = real_decoding_matrix+(ind_to_row[i]-k)*k*w*w;
        for (j = 0; j < w; j++) {
          b2 = ptr + j*k*w;
          for (y = 0; y < w; y++) {
            if (bitmatrix[index+j*k*w+i*w+y]) {
              for (z = 0; z < k*w; z++) {
                b2[z] = b2[z] ^ b1[z+y*k*w];
              }
            }
          }
        }
      }  
    }
  }

  free(row_ids);
  free(ind_to_row);
  return real_decoding_matrix;
}