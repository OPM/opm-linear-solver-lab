//**************************************************************************
//* Copyright (c) 2025, Advanced Micro Devices, Inc. All rights reserved.
//**************************************************************************

#include "vec.hpp"
#include <stdio.h>

void vec_show(const double *x, int n, const char *name)
{
    printf("%s=[\n",name);
    for(int i=0;i<n;i++) printf("%+.2f ",x[i]);
    printf("\n]\n\n");

}

void vec_ishow(const int *x, int n, const char *name)
{
    printf("%s=[\n",name);
    for(int i=0;i<n;i++) printf("%2d ",x[i]);
    printf("\n]\n\n");

}

void vec_fill(double *v, int n, double a)
{
    for(int i=0;i<n;i++) v[i]=a;
}

void vec_ifill(int *v, int n, int a)
{
    for(int i=0;i<n;i++) v[i]=a;
}

void vec_copy(const double *x, double *y, int n)
{
    for(int i=0;i<n;i++) y[i]=x[i];
}

double vec_inner(const double *x, const double *y, int n)
{
    double ans=0.0;
    for(int i=0;i<n;i++) ans+=x[i]*y[i];
    return ans;
}
