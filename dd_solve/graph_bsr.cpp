//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "bsr.hpp"
#include <stdlib.h>
#include <stdio.h>
#include "vec.hpp"



void bsr_graph_create(bsr_graph *M, int nrows, int ncols, int nnz)
{
    M->nrows=nrows;
    M->ncols=ncols;
    M->nnz=nnz;


    M->colidx=(int*)malloc(nnz*sizeof(int));
    M->rowptr=(int*)malloc((nrows+1)*sizeof(int));
    M->weight=(int*)malloc(nnz*sizeof(int));
}



void bsr_graph_part(
    int nx, int ny, int nz,
    int nblkx, int nblky, int nblkz,
    int* part_idxs
)
{
    int bx = nx/nblkx;
    int by = ny/nblky;
    int ibx, jby, kbz;
    int gidx, pidx;

    for(int i=0; i<nx; i++)
    {
        ibx = i/nblkx;
        for(int j=0; j<ny; j++)
        {
            jby = j/nblky;
            for(int k=0; k<nz; k++)
            {
                kbz = k/nblkz;
                gidx = i + nx*(j + ny*k);
//                gidx = i*nx*ny + j*nx + k;
                pidx = ibx + bx*(jby + by*kbz);
                part_idxs[gidx] = pidx;
            }
        }
    }
}


void bsr_graph_tags(bsr_graph *M, int *tag, int n)
{
    int nnz=n;
    int nrows=n;
    int ncols=0;
    for(int i=0;i<nnz;i++) ncols = tag[i]>ncols? tag[i] : ncols;
    ncols+=1;

    bsr_graph_create(M,nrows,ncols,nnz);
    for(int i=0;i<nnz;i++)
    {
        M->rowptr[i]=i;
        M->colidx[i]=tag[i];
        M->weight[i]=i;
    }
    M->rowptr[nnz]=nnz;
}

void  bsr_graph_show(bsr_graph *M, const char *name)
{
    printf("%s=[\n",name);
    int k=0;
    int j=0;
    for(int n=0;n<M->nrows;n++)
    {
        j=0;
        while(k<M->rowptr[n+1])
        {
            while(j<M->colidx[k])
            {
                //printf("%.1f ",0.0);
                printf("  -  ");
                j++;
            }
            printf(" %3d ",M->weight[k]);
            k++;
            j++;
        }
        while(j<M->ncols)
        {
            //printf("%.1f ",0.0);
            printf("  -  ");
            j++;
        }
        printf("\n");
    }
    printf("];\n\n");
}

void bsr_graph_transpose(bsr_graph *M, bsr_graph *T)
{
    int nrows=M->ncols;
    int ncols=M->nrows;
    int nnz=M->nnz;
    bsr_graph_create(T,nrows,ncols,nnz);

    for(int i=0;i<nrows+1;i++) T->rowptr[i]=0;
    for(int k=0;k<nnz;k++)
    {
        int j=M->colidx[k]+1;
        T->rowptr[j]++;
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
            T->weight[m]=M->weight[k];
            k++;
        }
    }

    for(int i=nrows;i>0;i--) T->rowptr[i]=T->rowptr[i-1];
    T->rowptr[0]=0;

    // sorting may be required
}
void bsr_graph_decompose(bsr_graph *M, bsr_graph *P, bsr_graph *D)
{
    int nparts   = P->nrows;
    int *partptr = P->rowptr;

    int nrows   = M->nrows;
    int ncols   = M->ncols;
    int nnz     = M->nnz;
    int *colidx = M->colidx;
    int *rowptr = M->rowptr;
    int *weight = M->weight;

    bsr_graph_create(D,nrows,ncols,nnz);
    for(int i=0;i<nnz;i++) D->weight[i]=i;

    int n=0;
    D->rowptr[0]=0;
    for(int i=0;i<nparts;i++)
    {
        for(int j=partptr[i]; j<partptr[i+1]; j++)
        {
            for(int k=rowptr[j]; k<rowptr[j+1]; k++)
            {
                if( colidx[k]>=partptr[i] && colidx[k]<partptr[i+1])
                {
                    D->colidx[n] = colidx[k];
                    D->weight[n] = weight[k];
                    n++;
                }
            }
            D->rowptr[j+1]=n;
        }
    }
    D->post_dd_edges=n;
    //D->nnz=n;
}

void bsr_graph_reorder(bsr_graph *M, int **pmap,bsr_graph *P)
{
    int nnz=M->nnz;
    int nrows=M->nrows;
    int ncols=M->ncols;

    bsr_graph_create(P,nrows,ncols,nnz);

    int *rowptr=P->rowptr;
    int *colidx=P->colidx;
    int *weight=P->weight;


    int *imap = (int*)malloc(nrows*sizeof(int));
    for(int i=0;i<nrows;i++) imap[(*pmap)[i]]=i;

    // permute row pointer offsets
    rowptr[0]=0;
    for(int i=0;i<nrows;i++)
    {
            int k=(*pmap)[i];
            rowptr[i+1]=M->rowptr[k+1]-M->rowptr[k];
    }
    for(int i=0;i<nrows;i++) rowptr[i+1]+=rowptr[i];

    int n=0;
    for(int i=0;i<nrows;i++)
    {
        // permute rows, columns, and nonzeros
        int m=(*pmap)[i];
        for(int j=M->rowptr[m];j<M->rowptr[m+1];j++)
        {
            int k=imap[M->colidx[j]];
            colidx[n]=k;
            weight[n]=M->weight[j];
            n++;
        }

        // sort row entries by colum index if necessary
        bool sorted=false;
        while(!sorted)
        {
            sorted=true;
            for(int j=rowptr[i];j<rowptr[i+1]-1;j++)
            {
                if (colidx[j]>colidx[j+1])
                {
                    int tmp;
                    tmp=colidx[j];
                    colidx[j]=colidx[j+1];
                    colidx[j+1]=tmp;

                    tmp=weight[j];
                    weight[j]=weight[j+1];
                    weight[j+1]=tmp;

                    sorted=false;
                }
            }
        }
    }
    free(imap);
};

void  bsr_graph_rowcount(bsr_graph *M, int **count)
{
    int nrows= M->nrows;
    *count=(int*)malloc(nrows*sizeof(int));
    for(int i=0;i<nrows;i++) (*count)[i]=M->rowptr[i+1]-M->rowptr[i];
}


