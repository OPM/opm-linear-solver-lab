//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "mat.hpp"
#include "vec.hpp"
#include <stdio.h>

void mat_show(const double *A, int m, int n, const char *name)
{
    printf("%s=[\n",name);
    for(int i=0;i<m;i++)
    { 
        for(int j=0;j<n;j++) printf("%+.2f ",A[n*i+j]);
        printf("\n");
    }
    printf("]\n\n");
}

void mat_identity(double *A, int m, int n)
{
    vec_fill(A,m*n,0);
    for(int i=0;i<n;i++) A[(n+1)*i]=1.0;
}

void mat_inv(const double *A, double *B,int n)
{
//    [a b c]
//    [d e f]
//    [g h i]
//    det = a*(ei-fh) - d*(bi-ch) + g*(bf-ce)

    double det= A[0]*( A[4]*A[8] - A[5]*A[7] ) - A[3]*( A[1]*A[8] - A[2]*A[7] ) + A[6]*( A[1]*A[5] - A[2]*A[4] ); 

    B[0] = ( A[4]*A[8] - A[5]*A[7] )/det;  
    B[1] =-( A[1]*A[8] - A[2]*A[7] )/det;
    B[2] = ( A[1]*A[5] - A[2]*A[4] )/det;
    B[3] =-( A[3]*A[8] - A[5]*A[6] )/det;
    B[4] = ( A[0]*A[8] - A[2]*A[6] )/det;
    B[5] =-( A[0]*A[5] - A[2]*A[3] )/det;
    B[6] = ( A[3]*A[7] - A[4]*A[6] )/det;
    B[7] =-( A[0]*A[7] - A[1]*A[6] )/det;
    B[8] = ( A[0]*A[4] - A[1]*A[3] )/det; 
}

void mat_matmul(const double *A, const double *B, double *C, int n)
{
    double bj[3];
    double b_ij;

    for(int j=0;j<n;j++)
    {
        for(int k=0;k<n;k++) bj[k]=B[k*n+j];
        for(int i=0;i<n;i++)
        {
            b_ij=0.0;
            for(int k=0;k<n;k++)
            {
                b_ij+=A[i*n+k]*bj[k];
            }
            C[i*n+j]=b_ij;
        }
    }
}

void mat_matfms(const double *A, const double *B, double *C, int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            for(int k=0;k<n;k++)
            {
                C[i*n+j]-=A[i*n+k]*B[k*n+j];
            }
        }
    }
}

void mat_vecmul(const double *A, double *b, int n)
{
    double bj[3];
    double bi;

    for(int j=0;j<n;j++) bj[j]=b[j];
    for(int i=0;i<n;i++)
    {
        bi=0.0;
        for(int j=0;j<n;j++)
        {
            bi+=A[i*n+j]*bj[j];
        }
        b[i]=bi;
    }
}

void mat_vecfms(const double *A, const double *b, double *c, int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            c[i]-=A[i*n+j]*b[j];
        }
    }
}


void mat_fma(const double *A, const double *x, double *y, int n)
{
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            y[i]+=A[i*n+j]*x[j];
        }
    }
}

void mat_dia(double *A, double c, int n)
{
    vec_fill(A,n*n,0);
    for(int i=0;i<n;i++) A[(n+1)*i]=c;
}
