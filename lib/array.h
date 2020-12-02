#ifndef H_ARRAY_H
#define H_ARRAY_H
#include <stdio.h>
#include <stdlib.h>


#define make_vector(v,n) ((v)=malloc((n)*sizeof *(v)))
#define free_vector(v) do{ free(v); v=NULL;}while(0)

int **make_imatrix(size_t m,size_t n){
    int **matrix;
    size_t i;
    make_vector(matrix,m+1);
    for (i=0;i<m;i++)
        make_vector(matrix[i],n);
    return matrix;
}

double **make_dmatrix(size_t m,size_t n){
    double **matrix;
    size_t i;
    make_vector(matrix,m+1);
    for (i=0;i<m;i++)
        make_vector(matrix[i],n);
    return matrix;
}


char **make_cmatrix(size_t m,size_t n){
    char **matrix;
    size_t i;
    make_vector(matrix,m+1);
    for (i=0;i<m;i++)
        make_vector(matrix[i],n);
    return matrix;
}



void free_imatrix(int **matrix){
    if(matrix != NULL){
        size_t i;
        for(i=0; matrix[i] != NULL; i++)
            free_vector(matrix[i]);
        free(matrix);
        matrix=NULL;
    }
}

void free_dmatrix(double **matrix){
    if(matrix != NULL){
        size_t i;
        for(i=0; matrix[i] != NULL; i++)
            free_vector(matrix[i]);
        free(matrix);
        matrix=NULL;
    }
}

void free_cmatrix(char **matrix){
    if(matrix != NULL){
        size_t i;
        for(i=0; matrix[i] != NULL; i++)
            free_vector(matrix[i]);
        free(matrix);
        matrix=NULL;
    }
}

#endif
