//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "bilu.hpp"
#include<stdio.h>
#include "mat.hpp"
#include "vec.hpp"


void bilu_apply(bilu_prec *P, double *x)
{
    bsr_matrix *L= &(P->L);
    bsr_matrix *U= &(P->U);

    int b=L->bsize;
    int bb=b*b;
//    csr_matrix *R= &(P->R);

    // Lower triangular solve assuming ones on diagonal
    for(int i=0;i<L->ncols;i++)
    {
        for(int k=L->rowptr[i];k<L->rowptr[i+1];k++)
        {
            int j=L->colidx[k];
            //x[j]-=L->data[k]*x[i];
            mat_vecfms(L->data+k*bb,x+b*i,x+b*j,b);
        }
    }

    // Upper triangular solve assuming nonzeros stored in original order
    double scale[9];
    for(int i=U->ncols;i>0;i--)
    {
        for(int k=U->rowptr[i]-1;k>U->rowptr[i-1];k--)
        {
            int j=U->colidx[k];
            //x[i-1]-=U->data[k]*x[j];
            mat_vecfms(U->data+k*bb,x+b*j,x+b*(i-1),b);
        }
        int d=U->rowptr[i-1];
        //double scale=1.0/U->data[d];
        //x[i-1]*=scale;
        mat_inv(U->data+d*bb,scale,b);
        mat_vecmul(scale,x+b*(i-1),b);
    }

/*
    // Upper triangular solve assuming nonzeros stored in reverse order
    double *y=x+R->nrows-1;
    for(int i=0;i<R->ncols;i++)
    {
        int d=R->rowptr[i];
        double scale=1.0/R->data[d];
        y[-i]*=scale;
        for(int k=R->rowptr[i]+1;k<R->rowptr[i+1];k++)
        {
            int j=R->colidx[k];
            y[-j]-=R->data[k]*y[-i];
        }
    }
*/
}

void bilu_create(bsr_matrix *M, bilu_prec *P)
{

    int b=M->bsize;
    int bb=b*b;

    //csr_matrix *U= &(P->U);
    bsr_matrix *L= &(P->L);
    bsr_matrix *U= &(P->U);
//    csr_matrix *R= &(P->R);

    bsr_matrix l;           // could temporarily use P->U, because P->L is strictly lower triangular
    bsr_lower(M,&l,true);
    bsr_transpose(&l,L);    // should reuse memory allocation

//    csr_matrix U;           // can we temporarily use P->R?
    bsr_upper(M,U,false);

    double scale[9];
    for(int i=0;i<M->nrows;i++)
    {
        int d=U->rowptr[i];
        //double scale =1.0/U->data[d];
        mat_inv(U->data+d*bb,scale,b);
        for(int k=L->rowptr[i];k<L->rowptr[i+1];k++)
        {
            //scale column i of L
            int j=L->colidx[k];
            //L->data[k]*=scale;
            mat_matmul(scale,L->data+k*bb,L->data+k*bb,b);

            //update diagonal of U
            d=U->rowptr[j];
            //U->data[d]-=U->data[i+k+1]*L->data[k];
            mat_matfms(U->data+(i+k+1)*bb,L->data+k*bb,U->data+d*bb,b);

            //NOT IMPLEMENTED!
            //update off-diagonal entries of U and L
            for(int m=L->rowptr[j];m<L->rowptr[j+1];m++)
            {
                if(L->colidx[m]==j)
                {
                    printf("ILU OFF_DIAGONALS NOT IMPLEMENTED!\n");
                    printf("(%d,%d)",m,j);
                    getchar();
                }
            }

        }

        // let a,b,c be off-diagonal indices of column i
        // assume a < b < c
        // coordinates (a,i), (b,i), (c,i)
        // structural symmetry implies (i,a), (i,b), (i,c)
        // diagonals (a,i)(i,a), (b,i)(i,b), (c,i)(i,c) always exist
        // do off-diagonals (a,i)(i,b), (a,i)(i,c), (b,i)(i,c) exist?
        // only if off-diagonals (a,b), (a,c), (b,c) exist
        // goto column b check for row a
        // goto column c, check for row a and b

//        break;
    }

    //assume no neigbors are connected to each other
    //works for two-point flux approximations on cartesian grids
/*
    csr_matrix T;
    csr_transpose(U,&T);  //transpose in-place? necessary?
    csr_reverse(&T,R);     //reverse in-place? necessary? bake into factorization?
*/
}

void bildu_apply(bildu_prec *P, double *x)
{
    bsr_matrix *L= &(P->L);
    bsr_matrix *D= &(P->D);
    bsr_matrix *U= &(P->U);

    int b=L->bsize;
    int bb=b*b;

    // Lower triangular solve assuming ones on diagonal
    for(int i=0;i<L->ncols;i++)
    {
        for(int k=L->rowptr[i];k<L->rowptr[i+1];k++)
        {
            int j=L->colidx[k];
            //x[j]-=L->data[k]*x[i];
            mat_vecfms(L->data+k*bb,x+b*i,x+b*j,b);
        }
    }
    
    //mat_show(x,L->nrows,b,"xx");
     
   //  return;
    // Muliply by (inverse) diagonal matrix
    
    
    for(int i=0;i<D->ncols;i++) mat_vecmul(D->data+bb*i,x+b*i,b);
    

    // Upper triangular solve assuming nonzeros stored in original order
    for(int i=U->ncols;i>0;i--)
    {
        for(int k=U->rowptr[i]-1;k>U->rowptr[i-1]-1;k--)
        {
            int j=U->colidx[k];
            //x[i-1]-=U->data[k]*x[j];
            mat_vecfms(U->data+k*bb,x+b*j,x+b*(i-1),b);
        }
    }
    
}


void bildu_create(bsr_matrix *M, bildu_prec *P)
{
    int b=M->bsize;
    int bb=b*b;

    bsr_matrix *L= &(P->L);
    bsr_matrix *D= &(P->D);
    bsr_matrix *U= &(P->U);

    bsr_matrix l;           // could temporarily use P->U, because P->L is strictly lower triangular
    bsr_lower(M,&l,true);
    bsr_transpose(&l,L);    // should reuse memory allocation

    bsr_diagonal(M,D);
    bsr_upper(M,U,true);

    double scale[9];
    for(int i=0;i<M->nrows;i++)
    {
        mat_inv(D->data+i*bb,scale,b);
        vec_copy(scale, D->data+bb*i, bb); //store inverse instead to simplify application
        for(int k=L->rowptr[i];k<L->rowptr[i+1];k++)
        {
            //scale column i of L
            mat_matmul(L->data+k*bb,scale,L->data+k*bb,b);

            //TODO: update diagonal of U
            int j=L->colidx[k];
            mat_matfms(U->data+k*bb,L->data+k*bb,D->data+j*bb,b);

            //TODO: scale row i of U
            mat_matmul(scale,U->data+k*bb,U->data+k*bb,b);

            //NOT IMPLEMENTED!
            //update off-diagonal entries of U and L
            for(int m=L->rowptr[j];m<L->rowptr[j+1];m++)
            {
                if(L->colidx[m]==j)
                {
                    printf("ILU OFF_DIAGONALS NOT IMPLEMENTED!\n");
                    printf("(%d,%d)",m,j);
                    getchar();
                }
            }

        }
    }

}
