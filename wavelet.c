#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "lib/wavelet.h"
#include "lib/array.h"

#define SQRT1_2   0.70710678118654752440


static void haar_transform_vector_forward(double *vector, int n){
    int i;
    int d;
    double h;
    double *qvector;

    h = sqrt (n);
    
    for(i=0; i<n; i++){
        vector[i] /= h;
    }

    qvector=malloc (n*sizeof(qvector));

    for (i=0;i<n;i++ )
    {
        qvector[i] = 0.0;
    }

    d = 1;
    while (d*2<=n)
    {
        d*=2;
    }
  
    while ( 1 < d )
    {
        d /= 2;
        for ( i = 0; i < d; i++ )
        {
            qvector[i]=(vector[2*i] + vector[2*i+1])*SQRT1_2;
            qvector[i+d]=( vector[2*i] - vector[2*i+1])*SQRT1_2;
        }
        for (i = 0; i < d * 2; i++)
        {
            vector[i] = qvector[i];
        }
    }

    free (qvector);
    
}

static void haar_transform_vector_reverse(double *vector, int n){

    int i;
    int d;
    double h;
    double *qvector;

    h = sqrt (n);
    
    for(i=0; i<n; i++){
        vector[i] *= h;
    }

    qvector =malloc(n*sizeof(qvector));

    for ( i = 0; i < n; i++ )
    {
        qvector[i] = 0.0;
    }

    d = 1;
    while (d*2 <= n)
    {
        for ( i = 0; i < d; i++ )
        {
            qvector[2*i]=(vector[i] + vector[i+d])*SQRT1_2;
            qvector[2*i+1]=(vector[i] - vector[i+d])*SQRT1_2;
        }
        for ( i = 0; i < d * 2; i++ )
        {
            vector[i] = qvector[i];
        }
        d*=2;
    }

    free (qvector);
    


}

void printTransformedVector(double *vector,int n){
    int i;
    printf("----------TRANSFORMED VECTOR----------\n");
    for(i=0;i<n;i++){
        printf("\t%f",vector[i]);
    }
    printf("\n");
}

void printReconstructedVector(double *vector,int n){
    int i;
    printf("----------RECONSTRUCTED VECTOR----------\n");
    for(i=0;i<n;i++){
        printf("\t%f",vector[i]);
    }
    printf("\n");
}

static void haar_transform_matrix_forward(double **matrix, int m,int n){

    int i,j;
    double *temp;
    temp=malloc(m*sizeof(temp));
    for(i = 0; i<m; i++){
        haar_transform_vector(matrix[i],n,WT_FWD);
    }

    for(j=0; j<n;j++){

        for(i=0; i<m;i++){
            temp[i]=matrix[i][j];
        }
        haar_transform_vector(temp,m,WT_FWD);
        for(i=0; i<m;i++){
            matrix[i][j]=temp[i];
        }
    }

    printf("----------TRANSFORMED MATRIX----------\n");
    for(i=0; i<m; i++){
        for(j=0;j<n;j++){
            printf("\t%f",matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

static void haar_transform_matrix_reverse(double **matrix, int m,int n){
    int i,j;
    double *temp;
    temp=malloc(m*sizeof(temp));
    for(i = 0; i<m; i++){
        haar_transform_vector(matrix[i],n,WT_REV);
    }

    for(j=0; j<n;j++){

        for(i=0; i<m;i++){
            temp[i]=matrix[i][j];
        }
        haar_transform_vector(temp,m,WT_REV);
        for(i=0; i<m;i++){
            matrix[i][j]=temp[i];
        }
    }

    printf("----------RECONSTRUCTED MATRIX----------\n");
    for(i=0; i<m; i++){
        for(j=0;j<n;j++){
            printf("\t%f",matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void haar_transform_vector(double *vector,int n,int dir){


    if( dir == WT_FWD ){
        haar_transform_vector_forward(vector,n);
    }else if( dir == WT_REV ){
        haar_transform_vector_reverse(vector,n);
    }else {
        
        exit(EXIT_FAILURE);
    }

}

void haar_transform_matrix(double **matrix,int m,int n,int dir){
    if( dir == WT_FWD ){
        haar_transform_matrix_forward(matrix,m,n);
    }else if( dir == WT_REV ){
        haar_transform_matrix_reverse(matrix,m,n);
    }else {
        exit(EXIT_FAILURE);
    }

}

int main(void){
    int i,j,m=4,n=8;
    double *vector,**matrix;
    make_vector(vector,n);
    matrix=make_dmatrix(m,n);

    printf("-----------------ORIGINAL VECTOR---------\n");
    for(i=0;i<n;i++){
        vector[i]=sqrt((i+1.0)/2000);
    }

    for(i=0;i<n;i++){
        printf("\t%f",vector[i]);
    }
    printf("\n");
    haar_transform_vector(vector,n,0);
    printTransformedVector(vector,n);
    haar_transform_vector(vector,n,1);
    printReconstructedVector(vector,n);
    
    

    for(i=0; i<m; i++){
        for(j=0;j<n;j++){
            matrix[i][j] = 1.0/(1+i+j);
        }
    }
    
    free(vector);
    printf("-----------------ORIGINAL MATRIX---------\n");
    for(i=0; i<m; i++){
        for(j=0;j<n;j++){
            printf("\t%f",matrix[i][j]);
        }
        printf("\n");
    }
    haar_transform_matrix(matrix,m,n,0);
    haar_transform_matrix(matrix,m,n,1);
    
    free(matrix);
}
