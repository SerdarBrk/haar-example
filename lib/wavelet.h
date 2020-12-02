#ifndef H_WAVELET_H
#define H_WAVELET_H
#define WT_FWD 0
#define WT_REV 1

void haar_transform_vector(double *vector,int n,int dir);
void haar_transform_matrix(double **matrix,int m,int n,int dir);

#endif