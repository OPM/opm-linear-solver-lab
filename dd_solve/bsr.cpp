//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "bsr.hpp"
#include "vec.hpp"
#include "mat.hpp"
#include <stdlib.h>
#include <stdio.h>

void bsr_dirichlet(bsr_matrix *M, int i)
{
    int b = M->bsize;
    int bb=b*b;

    for(int j=M->rowptr[i];j<M->rowptr[i+1];j++)
    {
        double c = (i==M->colidx[j] ? 1.0 : 0.0);
        mat_dia(M->data+bb*j,c,b);
    }
}

void bsr_to_graph(bsr_matrix *M, bsr_graph *G)
{
    int nrows=M->nrows;
    int ncols=M->ncols;
    int nnz=M->nnz;

    bsr_graph_create(G, nrows, ncols, nnz);
    for(int i=0;i<nrows+1;i++) G->rowptr[i]=M->rowptr[i];
    for(int i=0;i<nnz;i++) G->colidx[i]=M->colidx[i];
    for(int i=0;i<nnz;i++) G->weight[i]=i;
}

void bsr_reorder(bsr_matrix *M, bsr_graph *G, bsr_matrix *T)
{
    T->nrows=G->nrows;
    T->ncols=G->ncols;
    T->nnz=G->nnz;
	T->bsize=M->bsize;
	int bsize = T->bsize;
    T->rowptr=G->rowptr;
    T->colidx=(int*)malloc(G->nnz*sizeof(int));
    T->data=(double*)malloc(bsize*bsize*G->nnz*sizeof(double));
    for(int i=0;i<G->nnz;i++)
    {
        int k=G->weight[i];
		for (int j=0; j<bsize*bsize;j++)
		{
			T->data[i*bsize*bsize+j]=M->data[k*bsize*bsize+j];
		}
        T->colidx[i]=G->colidx[i];
    }
}


void  bsr_dot(bsr_matrix *M, const double *x, double *y)
{
    int b=M->bsize;
    int bb=b*b;

    int k=0;
    int j;
    for(int i=0;i<M->nrows;i++)
    {
        //y[i]=0;
        vec_fill(y+i*b,b,0.0);
        while(k<M->rowptr[i+1])
        {
            j=M->colidx[k];
            mat_fma(M->data+k*bb,x+j*b,y+i*b,b);
            k++;
        }
    }
}

void bsr_create(bsr_matrix *M, int nrows, int ncols, int nnz, int bsize)
{
    M->bsize=bsize;
    M->nrows=nrows;
    M->ncols=ncols;
    M->nnz=nnz;

    M->colidx=(int*)malloc(nnz*sizeof(int));
    M->rowptr=(int*)malloc((nrows+1)*sizeof(int));
    M->data=(double*)malloc(bsize*bsize*nnz*sizeof(double));
}




void bsr_transpose(bsr_matrix *M, bsr_matrix *T)
{
    int b   = M->bsize;
    int bb  = b*b;
    int nrows=M->ncols;
    int ncols=M->nrows;
    int nnz=M->nnz;
    bsr_create(T,nrows,ncols,nnz,b);

    for(int i=0;i<nrows+1;i++) T->rowptr[i]=0;
    for(int k=0;k<nnz;k++)
    {
        int j=M->colidx[k]+1;
        //printf("%d\n",j);
        T->rowptr[j]+=1;
    }
    for(int i=0;i<nrows;i++) T->rowptr[i+1]+=T->rowptr[i];
    //assert T->rowptr[nrows]==nnz;

    int k=0;
    for(int i=0;i<M->nrows;i++)
    {
        while(k<M->rowptr[i+1])
        {
            int j = M->colidx[k];
            int m = T->rowptr[j]++;
            T->colidx[m]=i;
            //T->data[m]=M->data[k];
            vec_copy(M->data+bb*k, T->data+bb*m, bb);
            k++;
        }
    }

    for(int i=nrows;i>0;i--) T->rowptr[i]=T->rowptr[i-1];
    T->rowptr[0]=0;

    // sorting may be required
}

void  bsr_diagonal(bsr_matrix *M, bsr_matrix *D)
{
    int n   = M->nrows;
    int b   = M->bsize;
    int bb  = b*b;
    bsr_create(D,n,n,n,b);

    int k=0;
    D->rowptr[0]=0;
    for(int i=0;i<n;i++)
    {
        while(k<M->rowptr[i+1])
        {
            int j=M->colidx[k];
            if(j==i)
            {
                D->colidx[i]=j;
                vec_copy(M->data+bb*k, D->data+bb*i, bb);
            }
            k++;
        }
        D->rowptr[i+1]=i+1;
    }
}



void  bsr_lower(bsr_matrix *M, bsr_matrix *L, bool strict)
{
    // assume structural symmetry
    int b   = M->bsize;
    int bb  = b*b;
    int n   = M->nrows;
    int nnz = (M->nnz - n)/2 + (strict ? 0 : n);
    bsr_create(L,n,n,nnz,b);

    int m=0;
    int k=0;
    L->rowptr[0]=0;
    for(int i=0;i<n;i++)
    {
        while(k<M->rowptr[i+1])
        {
            int j=M->colidx[k];
            if(j< (strict ? i : i+1))
            {
                L->colidx[m]=j;
                vec_copy(M->data+bb*k, L->data+bb*m, bb);
                m++;
            }
            k++;
        }
        L->rowptr[i+1]=m;
    }
}


void  bsr_upper(bsr_matrix *M, bsr_matrix *U, bool strict)
{
    // assume structural symmetry
    int b   = M->bsize;
    int bb  = b*b;
    int n   = M->nrows;
    int nnz = (M->nnz - n)/2 + (strict ? 0 : n);
    bsr_create(U,n,n,nnz,b);

    int m=0;
    int k=0;
    U->rowptr[0]=0;
    for(int i=0;i<n;i++)
    {
        while(k<M->rowptr[i+1])
        {
            int j=M->colidx[k];
            if(j> (strict ? i : i-1))
            {
                U->colidx[m]=j;
                vec_copy(M->data+bb*k, U->data+bb*m, bb);
                m++;
            }
            k++;
        }
        U->rowptr[i+1]=m;
    }
}

void  bsr_cartesian(bsr_matrix *M, int nx, int ny, int nz, int bsize)
{
    int n=nx*ny*nz;
    int nnz= nx*ny*nz +2*(nx-1)*ny*nz + 2*nx*(ny-1)*nz + 2*nx*ny*(nz-1);

    M->bsize=bsize;
    M->nrows=n;
    M->ncols=n;
    M->nnz=nnz;

    int b=bsize;
    int bb=b*b;

    M->colidx=(int*)malloc(nnz*sizeof(int));
    M->rowptr=(int*)malloc((n+1)*sizeof(int));
    M->data=(double*)malloc(bb*nnz*sizeof(double));

    vec_fill(M->data,bb*nnz,0.0);

    int m=0;
    M->rowptr[0]=0;
    for(int i=0;i<nz;i++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int k=0;k<nx;k++)
            {
                if(i>0) M->colidx[m++] = (i-1)*nx*ny +j*nx + k;
                if(j>0) M->colidx[m++] = i*nx*ny +(j-1)*nx + k;
                if(k>0) M->colidx[m++] = i*nx*ny + j*nx + k-1;

                M->colidx[m++] = i*nx*ny + j*nx + k;

                if((nx-k)>1) M->colidx[m++] = i*nx*ny + j*nx + k+1;
                if((ny-j)>1) M->colidx[m++] = i*nx*ny +(j+1)*nx + k;
                if((nz-i)>1) M->colidx[m++] = (i+1)*nx*ny +j*nx + k;

                M->rowptr[i*nx*ny+j*nx+k+1]=m;
            }
        }
    }
}

void bsr_flux_balance(bsr_matrix *M)
{
    int b=M->bsize;
    int bb=b*b;
    

    int k=0;
    int d=-1;
    double f;
    for(int i=0;i<M->nrows;i++)
    {
        f=0.0;
        while(k<M->rowptr[i+1])
        {
            if(M->colidx[k]==i) d=k++;   
            else
            {
                //M->data[k++]=-1.0;
                mat_dia(M->data+bb*k++,-1.0,b);
                f+=1.0;
            }
        }
        //M->data[d]=f;
        mat_dia(M->data+bb*d,f,b);
    }
}


void  bsr_laplacian(bsr_matrix *M, int nx, int ny, int nz, int bsize)
{
    int n=nx*ny*nz;
    int nnz= nx*ny*nz +2*(nx-1)*ny*nz + 2*nx*(ny-1)*nz + 2*nx*ny*(nz-1);

    M->bsize=bsize;
    M->nrows=n;
    M->ncols=n;
    M->nnz=nnz;

    int b=bsize;
    int bb=b*b;

    M->colidx=(int*)malloc(nnz*sizeof(int));
    M->rowptr=(int*)malloc((n+1)*sizeof(int));
    M->data=(double*)malloc(bb*nnz*sizeof(double));

    int m=0;
    M->rowptr[0]=0;
    for(int i=0;i<nz;i++)
    {
        for(int j=0;j<ny;j++)
        {
            for(int k=0;k<nx;k++)
            {
                if(i>0)
                {
                    M->colidx[m] = (i-1)*nx*ny +j*nx + k;
                    //M->data[m]=-1.0;
                    mat_dia(M->data+bb*m,-1.0,b);
                    m++;
                }
                if(j>0)
                {
                    M->colidx[m] = i*nx*ny +(j-1)*nx + k;
                    //M->data[m]=-1.0;
                    mat_dia(M->data+bb*m,-1.0,b);
                    m++;
                }
                if(k>0)
                {
                    M->colidx[m] = i*nx*ny + j*nx + k-1;
                    //M->data[m]=-1.0;
                    mat_dia(M->data+bb*m,-1.0,b);
                    m++;
                }
                {
                    M->colidx[m] = i*nx*ny + j*nx + k;
                    //M->data[m]=6.0;
                    mat_dia(M->data+bb*m,6.0,b);
                    m++;
                }
                if((nx-k)>1)
                {
                    M->colidx[m] = i*nx*ny + j*nx + k+1;
                    //M->data[m]=-1.0;
                    mat_dia(M->data+bb*m,-1.0,b);
                    m++;
                }
                if((ny-j)>1)
                {
                    M->colidx[m] = i*nx*ny +(j+1)*nx + k;
                    //M->data[m]=-1.0;
                    mat_dia(M->data+bb*m,-1.0,b);
                    m++;
                }
                if((nz-i)>1)
                {
                    M->colidx[m] = (i+1)*nx*ny +j*nx + k;
                    //M->data[m]=-1.0;
                    mat_dia(M->data+bb*m,-1.0,b);
                    m++;
                }
                M->rowptr[i*nx*ny+j*nx+k+1]=m;
            }
        }
    }
}

void bsr_identity(bsr_matrix* M, int n, int bsize)
{
    bsr_create(M,n,n,n,bsize);

    for(int i=0;i<n;i++)
    {
        M->colidx[i]=i;
        M->rowptr[i]=i;
        mat_identity(M->data+i*bsize*bsize, bsize, bsize);
    }
    M->rowptr[n]=n;
}

void  bsr_show(bsr_matrix *M, const char *name)
{
    int b=M->bsize;
    int bb=b*b;

    printf("%s=[\n",name);
    int k=0;
    int j=0;
    for(int n=0;n<M->nrows;n++)
    {
        for(int ii=0;ii<b;ii++)
        {
            j=0;
            if(ii>0) k=M->rowptr[n];
            while(k<M->rowptr[n+1])
            {
                while(j<M->colidx[k])
                {
                    for(int jj=0;jj<b;jj++) printf("   -  ");
                    j++;
                }
                for(int jj=0;jj<b;jj++) printf("%+5.2f ",M->data[bb*k+b*ii+jj]);
                k++;
                j++;
            }
            while(j<M->ncols)
            {
                for(int jj=0;jj<b;jj++) printf("   -  ");
                j++;
            }
            printf("\n");
        }
    }
    printf("];\n\n");
}
