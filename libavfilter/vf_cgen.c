/*
 * Copyright (C) 2018 Christian Bacher
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with FFmpeg; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

/**
 * @file
 * Video generation based on complex functions
 */

#include <float.h>
#include <stdio.h>
#include <string.h>
#include "libavutil/avassert.h"
#include "libavutil/avstring.h"
#include "libavutil/opt.h"
#include "libavutil/pixdesc.h"
#include "internal.h"
#include <complex.h>

#define PARAMSIZE 10

typedef struct CGenContext {
    const AVClass *class;
    const char *f1;
    double p1[PARAMSIZE];
    const char *f2;
    double p2[PARAMSIZE];
    const char *f3;
    double p3[PARAMSIZE];
    const char *nf;
    double np[PARAMSIZE];
    const char *wf;
    double wp[PARAMSIZE];
    int cp3; // copy params from 0=none, 1=p1, 2=p2
    int cp2; // copy params from 0=none, 1=p1, 3=p3
    uint8_t (*funcs[3])(int,int,int,int,int,double*);
    uint8_t (*funcs2[3])(int,int,int,int,int,struct CGenContext*);
    double (*nfunc)(int,double*);
    double (*wfunc)(int,double*);
    int hsub, vsub;             ///< chroma subsampling
    int planes;                 ///< number of planes
    int is_rgb;
    int bps;
    int offset;
} CGenContext;

typedef struct Func2 {
    const char *name;
    uint8_t (*f)(int x, int y, int w, int h, int n, CGenContext *ctx);
} Func2;

typedef struct Func {
    const char *name;
    uint8_t (*f)(int x, int y, int w, int h, int n, double *p);
} Func;

typedef struct NFunc {
    const char *name;
    double (*f)(int n, double *p);
} NFunc;

static double complex rotate(double complex z, double angle) {
     double phi = angle*(M_PI/180);
     return (cos(phi) + sin(phi)*I) * z;
}

static double complex poly0(double complex z) {
    return z * z;
}


static double complex poly1(double complex z) {
    return z * z * z - (8+0*I);
}

static double complex poly2(double complex z) {
    return ((z*z-(1+0*I)) * (z-(2-I)) * (z-(2-I))) / (z*z + (2+2*I));
}

static double complex poly3(double complex z) {
    return (z*z-(2+0*I)) / (z*z*z + (1+0*I));
}

static double complex poly4(double complex z) {
    return (z - 2 + 0*I) * (z + 0+2*I) * (z + 2+0*I) / (z + 0-0*I) * (z + 2-2*I) * (z - 2-2*I);
}

static double complex poly5(double complex z) {
    double complex z1 = rotate(z,-30);
    return 1/((z1 * z1) * (z1 * z1) * (z1 * z1) *z1 - z1);
}

static double complex id(double complex z) {
    return z;
}

static double complex expinv(double complex z) {
    return cexp(1/z);
}

static double complex expinv1(double complex z) {
    return 0.85 + 0.9 * z * cexp( I * -7.3/(1/pow(cabs(z),2)));
}

static double complex fexpinv2(double complex z, double *p, int n) {
    //return 0.85 + 0.9 * z * cexp( I * -7.3/(1/pow(cabs(z),2)));
    //return 0.999 * z * cexp( I * -7.3/(1/pow(cabs(z),2)));
    return p[0] * z * cexp( I * p[1]/pow(cabs(z),p[3]));
}

static double complex fexpinv1(double complex z, double *p, int n) {
    return p[0] + (p[1]?p[1]:1) * z * cexp( I * (p[2]?p[2]:1)/(1/pow(cabs(z),2)));
}

static double complex fpoly4(double complex z, double *p, int n) {
    return (z - 2 + 0*I) * (z + 0+2*I) * (z + 2+0*I) / (z + 0+0.002*n*I) * (z + 2-2*I) * (z - 2-2*I);
}

static double complex sininv(double complex z) {
    return csin(1/z);
}

static double complex sinz(double complex z) {
    return csin(z);
}

static double complex tanz(double complex z) {
    return ctan(z);
}
static void get_hsv2rgb_params(double phi, double S, double V, uint8_t *h, double *p, double *q, double *t) {
    double H = phi<0?((180/M_PI)*phi)+360:(180/M_PI)*phi;
    *h = floor(H/60.0);
    double f = H/60.0 - *h;
    *p = V*(1-S);
    *q = V*(1-S*f);
    *t = V*(1-S*(1-f)); 
}

static uint8_t hsv_r(double phi, double S, double V) {
    double p,q,t;
    uint8_t h;
    get_hsv2rgb_params(phi, S, V, &h, &p, &q, &t);
    uint8_t retval = 0;
    switch(h) {
        case 0:
        case 5: retval = V * 255; break;
        case 1: retval = q * 255; break;
        case 2:
        case 3: retval = p * 255; break;
        case 4: retval = t * 255; break;        
    }
    return retval;
}

static uint8_t hsv_g(double phi, double S, double V) {
    double p,q,t;
    uint8_t h;
    get_hsv2rgb_params(phi, S, V, &h, &p, &q, &t);
    uint8_t retval = 0;
    switch(h) {
        case 1:
        case 2: retval = V * 255; break;
        case 3: retval = q * 255; break;
        case 4:
        case 5: retval = p * 255; break;
        case 0: retval = t * 255; break;        
    }
    return retval;
}

static uint8_t hsv_b(double phi, double S, double V) {
    double p,q,t;
    uint8_t h;
    get_hsv2rgb_params(phi, S, V, &h, &p, &q, &t);
    uint8_t retval = 0;
    switch(h) {
        case 3:
        case 4: retval = V * 255; break;
        case 5: retval = q * 255; break;
        case 0:
        case 1: retval = p * 255; break;
        case 2: retval = t * 255; break;        
    }
    return retval;
}

/* ============================= NFUNCS *******************************************/

static double poly(int n, double *p) {
    double x = (double)n * p[0];
    return p[1] + p[2]*x + p[3]*x+x + p[4]*x*x*x + p[5]*x*x*x*x;
}

static double nsin(int n, double *p) {
    return (sin(p[0]*n - M_PI/2)+1)*p[1]; 
}

static double idn(int n, double *p) {
    return n;
}

static double constant(int n, double *p) {
    return p[0]?p[0]:1;
}

static double lin(int n, double *p) {
    return p[0] + p[1]*n;
}

static double inv(int n, double *p) {
    double x = (double)n/p[0];
    return p[1]*(1/((double)x+1));
}

static double pchain(int n, double *p) {
    int len = p[0];
    int section = n/len + 1;
    if((section+1) >= PARAMSIZE) return 0;
    double x = (double)(n%len)/(double)len;
    return p[section] + (p[section+1]-p[section])*x;
}

static double gauss(int n, double *p) {
    double m = p[1]==0?1:p[1];
    double fac = p[2]==0?0.5:p[2];
    double x = (double)n/p[0];
    return m * exp(-fac * pow(x-p[3],2)); 
}

static NFunc nfuncs[] = {
    {"idn",idn},
    {"constant",constant},
    {"lin",lin},
    {"poly",poly},
    {"pchain",pchain},
    {"sin",nsin},
    {"inv",inv},
    {"gauss",gauss},
    {NULL,NULL}
};

static double (*getNFunc(const char *name))(int,double*) {
    int k=0;
    while(nfuncs[k].name) {
        if(!strcmp(name, nfuncs[k].name)) {
            return nfuncs[k].f;
        }
        k++;
    }
    return NULL;
}

uint8_t (*hsv[3])(double phi, double S, double V) = { hsv_g, hsv_b, hsv_r };


double complex **last;

/* ============================= FUNCS2 =========================================== */
static uint8_t ffunc(double complex (*func)(double complex z, double *p, int n), int plane, int x, int y, int w, int h, int n, CGenContext *ctx) {
    double f = ctx->wfunc(n-ctx->offset,ctx->wp)/(double)w;
    double complex z = func((x*f+ctx->p1[4]+ctx->p1[2]*n) + ((y*f+ctx->p1[5]+ctx->p1[3]*n) * I), &ctx->p1[6],n-ctx->offset);
    return ctx->nfunc(n,ctx->np) * ctx->p1[0] * cabs(z) * hsv[plane](carg(z), ctx->p1[1], 1);
}


static uint8_t ffuncg(double complex (*func)(double complex z, double *p, int n), int x, int y, int w, int h, int n, CGenContext *ctx) {
    double f = ctx->wfunc(n-ctx->offset,ctx->wp)/(double)w;
    double complex z = func((x*f+ctx->p1[4]+ctx->p1[2]*n) + ((y*f+ctx->p1[5]+ctx->p1[3]*n) * I), &ctx->p1[6],n-ctx->offset);
    return ctx->nfunc(n,ctx->np) * ctx->p1[0] * cabs(z) * hsv_g(carg(z), ctx->p1[1], 1);
}

static uint8_t ffuncb(double complex (*func)(double complex z, double *p, int n), int x, int y, int w, int h, int n, CGenContext *ctx) {
    double f = ctx->wfunc(n-ctx->offset,ctx->wp)/(double)w;
    double complex z = func((x*f+ctx->p2[4]+ctx->p2[2]*n) + ((y*f+ctx->p2[5]+ctx->p2[3]*n) * I), &ctx->p2[6],n-ctx->offset);
    return ctx->nfunc(n,ctx->np) * ctx->p2[0] * cabs(z) * hsv_b(carg(z), ctx->p2[1], 1);
}

static uint8_t ffuncr(double complex (*func)(double complex z, double *p, int n), int x, int y, int w, int h, int n, CGenContext *ctx) {
    double f = ctx->wfunc(n-ctx->offset,ctx->wp)/(double)w;
    double complex z = func((x*f+ctx->p3[4]+ctx->p3[2]*n) + ((y*f+ctx->p3[5]+ctx->p3[3]*n) * I), &ctx->p3[6],n-ctx->offset);
    return ctx->nfunc(n,ctx->np) * ctx->p3[0] * cabs(z) * hsv_r(carg(z), ctx->p3[1], 1);
}

static uint8_t fexpinvg(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return ffuncg(fexpinv1,x,y,w,h,n,ctx);
}

static uint8_t fexpinvb(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return ffuncb(fexpinv1,x,y,w,h,n,ctx);
}

static uint8_t fexpinvr(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return ffuncr(fexpinv1,x,y,w,h,n,ctx);
}

static uint8_t nfuncg(double complex (*func)(double complex z), int x, int y, int w, int h, int n, CGenContext *ctx) {
    double f = ctx->wfunc(n-ctx->offset,ctx->wp)/(double)w;
    double complex z = func((x*f+ctx->p1[4]+ctx->p1[2]*n) + ((y*f+ctx->p1[5]+ctx->p1[3]*n) * I));
    return ctx->nfunc(n,ctx->np) * ctx->p1[0] * cabs(z) * hsv_g(carg(z), ctx->p1[1], 1);
}

static uint8_t nfuncb(double complex (*func)(double complex z), int x, int y, int w, int h, int n, CGenContext *ctx) {
    double f = ctx->wfunc(n-ctx->offset,ctx->wp)/(double)w;
    double complex z = func((x*f+ctx->p2[4]+ctx->p2[2]*n) + ((y*f+ctx->p2[5]+ctx->p2[3]*n) * I));
    return ctx->nfunc(n,ctx->np) * ctx->p2[0] * cabs(z) * hsv_b(carg(z), ctx->p2[1], 1);
}

static uint8_t nfuncr(double complex (*func)(double complex z), int x, int y, int w, int h, int n, CGenContext *ctx) {
    double f = ctx->wfunc(n-ctx->offset,ctx->wp)/(double)w;
    double complex z = func((x*f+ctx->p3[4]+ctx->p3[2]*n) + ((y*f+ctx->p3[5]+ctx->p3[3]*n) * I));
    return ctx->nfunc(n,ctx->np) * ctx->p3[0] * cabs(z) * hsv_r(carg(z), ctx->p3[1], 1);
}

static uint8_t nsininvg(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return nfuncg(sininv,x,y,w,h,n,ctx);
}

static uint8_t nsininvb(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return nfuncb(sininv,x,y,w,h,n,ctx);
}

static uint8_t nsininvr(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return nfuncr(sininv,x,y,w,h,n,ctx);
}

static uint8_t nexpinvg(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return nfuncg(expinv1,x,y,w,h,n,ctx);
}

static uint8_t nexpinvb(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return nfuncb(expinv1,x,y,w,h,n,ctx);
}

static uint8_t nexpinvr(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return nfuncr(expinv1,x,y,w,h,n,ctx);
}


static uint8_t ntang(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return nfuncg(tanz,x,y,w,h,n,ctx);
}

static uint8_t ntanb(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return nfuncb(tanz,x,y,w,h,n,ctx);
}

static uint8_t ntanr(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return nfuncr(tanz,x,y,w,h,n,ctx);
}

static uint8_t npoly5g(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return nfuncg(poly5,x,y,w,h,n,ctx);
}

static uint8_t npoly5b(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return nfuncb(poly5,x,y,w,h,n,ctx);
}

static uint8_t npoly5r(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return nfuncr(poly5,x,y,w,h,n,ctx);
}

static uint8_t fpoly4g(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return ffuncg(fpoly4,x,y,w,h,n,ctx);
}

static uint8_t fpoly4b(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return ffuncb(fpoly4,x,y,w,h,n,ctx);
}

static uint8_t fpoly4r(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return ffuncr(fpoly4,x,y,w,h,n,ctx);
}

static uint8_t nzero(int x, int y, int w, int h, int n, CGenContext *ctx) {
    return 0;
}

static uint8_t sqasing(int x, int y, int w, int h, int n, CGenContext *ctx) {
    double f = ctx->wfunc(n-ctx->offset,ctx->wp)/(double)h;
    double x1=f*x;
    double y1=f*y;
    return asin(sqrt(x1*x1+y1*y1)/(h/2))*n*ctx->p1[0];
}

static uint8_t sqasinb(int x, int y, int w, int h, int n, CGenContext *ctx) {
    double f = ctx->wfunc(n-ctx->offset,ctx->wp)/(double)h;
    double x1=f*x;
    double y1=f*y;
    return asin(sqrt(x1*x1+y1*y1)/(h/2))*n*ctx->p2[0];
}

static uint8_t sqasinr(int x, int y, int w, int h, int n, CGenContext *ctx) {
    double f = ctx->wfunc(n-ctx->offset,ctx->wp)/(double)h;
    double x1=f*x;
    double y1=f*y;
    return asin(sqrt(x1*x1+y1*y1)/(h/2))*n*ctx->p3[0];
}


static uint8_t map1g(int x, int y, int w, int h, int n, CGenContext *ctx) {
    int nn = n - ctx->offset;
    double f = ctx->p1[0];
    if(nn == 0) {
        double complex z = x*f + y*f * I;
        last[x+960][y+540] = z;
        return hsv_g(carg(z),ctx->p1[1], cabs(z)/(h/2));
    } else {
        double complex z = expinv1(last[x+960][y+540]);
        last[x+960][y+540] = z;
        return hsv_g(carg(z),ctx->p1[1], cabs(z)/(h/2));
    }
}

static uint8_t map1b(int x, int y, int w, int h, int n, CGenContext *ctx) {
    int nn = n - ctx->offset;
    double f = ctx->p2[0];
    if(nn == 0) {
        double complex z = x*f + y*f * I;
        last[x+960][y+540] = z;
        return hsv_b(carg(z),ctx->p2[1], cabs(z)/(h/2));
    } else {
        double complex z = expinv1(last[x+960][y+540]);
        last[x+960][y+540] = z;
        return hsv_b(carg(z),ctx->p2[1], cabs(z)/(h/2));
    }
}

static uint8_t map1r(int x, int y, int w, int h, int n, CGenContext *ctx) {
    int nn = n - ctx->offset;
    double f = ctx->p3[0];
    if(nn == 0) {
        double complex z = x*f + y*f * I;
        last[x+960][y+540] = z;
        return hsv_r(carg(z),ctx->p3[1], cabs(z)/(h/2));
    } else {
        double complex z = expinv1(last[x+960][y+540]);
        last[x+960][y+540] = z;
        return hsv_r(carg(z),ctx->p3[1], cabs(z)/(h/2));
    }
}


static Func2 func2s[] = {
    { "nzero", nzero },
    { "fpoly4r", fpoly4r },
    { "fpoly4g", fpoly4g },
    { "fpoly4b", fpoly4b },
    { "npoly5r", npoly5r },
    { "npoly5g", npoly5g },
    { "npoly5b", npoly5b },
    { "ntanb", ntanb },
    { "ntanr", ntanr },
    { "ntang", ntang },
    { "nsininvb", nsininvb },
    { "nsininvr", nsininvr },
    { "nsininvg", nsininvg },
    { "nexpinvb", nexpinvb },
    { "nexpinvr", nexpinvr },
    { "nexpinvg", nexpinvg },
    { "sqasinb", sqasinb },
    { "sqasinr", sqasinr },
    { "sqasing", sqasing },
    { "fexpinvb", fexpinvb },
    { "fexpinvr", fexpinvr },
    { "fexpinvg", fexpinvg },
    { "map1b", map1b },
    { "map1r", map1r },
    { "map1g", map1g },
    { NULL, NULL }
};

static uint8_t (*getFunc2(const char *name))(int,int,int,int,int,CGenContext*) {
    int k=0;
    while(func2s[k].name) {
        if(!strcmp(name, func2s[k].name)) {
            return func2s[k].f;
        }
        k++;
    }
    return NULL;
}

/* ============================= FUNCS ===========================================*/

static uint8_t zero(int x, int y, int w, int h, int n, double *p) {
    return 0;
}
static uint8_t sq(int x, int y, int w, int h, int n, double *p) {
    return (sqrt(x*x+y*y)/(w/2))*n*p[0];
}

static uint8_t poly1re(int x, int y, int w, int h, int n, double *p) {
    double f = (double)n/((double)w*(p[0]==0?1:p[0]));
    double complex z = x*f + (y*f) * I;
    return 128 + creal(poly1(z));
}
static uint8_t poly1img(int x, int y, int w, int h, int n, double *p) {
    double f = (double)n/((double)w*(p[0]==0?1:p[0]));
    double complex z = x*f + (y*f) * I;
    return 128 + cimag(poly1(z));
}
static uint8_t poly1abs(int x, int y, int w, int h, int n, double *p) {
    double f = (double)n/((double)w*(p[0]==0?1:p[0]));
    double complex z = x*f + (y*f) * I;
    return cabs(poly1(z));
}

static uint8_t idre(int x, int y, int w, int h, int n, double *p) {
    double f = (double)n/((double)w*(p[0]==0?1:p[0]));
    double complex z = x*f + (y*f) * I;
    return 128 + creal(id(z));
}
static uint8_t idimg(int x, int y, int w, int h, int n, double *p) {
    double f = (double)n/((double)w*(p[0]==0?1:p[0]));
    double complex z = x*f + (y*f) * I;
    return 128 + cimag(id(z));
}

static uint8_t idabs(int x, int y, int w, int h, int n, double *p) {
    double f = (double)n/((double)w*(p[0]==0?1:p[0]));
    double complex z = x*f + (y*f) * I;
    return cabs(id(z));
}

static uint8_t funcg(double complex (*func)(double complex z),int x, int y, int w, int h, int n, double *p) {
    double f = (p[0]+p[5]*n)/(double)w;
    double complex z = func(x*f + (y*f) * I);
    return ((double)n + 1 + p[3] * pow(n,p[4])) * p[2] * cabs(z) * hsv_g(carg(z), p[1], 1);
}

static uint8_t funcb(double complex (*func)(double complex z),int x, int y, int w, int h, int n, double *p) {
    double f = (p[0]+p[5]*n)/(double)w;
    double complex z = func(x*f + (y*f) * I);
    return ((double)n + 1 + p[3] * pow(n,p[4])) * p[2] * cabs(z) * hsv_b(carg(z), p[1], 1);
}

static uint8_t funcr(double complex (*func)(double complex z),int x, int y, int w, int h, int n, double *p) {
    double f = (p[0]+p[5]*n)/(double)w;
    double complex z = func(x*f + (y*f) * I);
    return ((double)n + 1 + p[3] * pow(n,p[4])) * p[2] * cabs(z) * hsv_r(carg(z), p[1], 1);
}

static uint8_t tang(int x, int y, int w, int h, int n, double *p) {
    return funcg(tanz,x,y,w,h,n,p);
}

static uint8_t tanb(int x, int y, int w, int h, int n, double *p) {
    return funcb(tanz,x,y,w,h,n,p);
}

static uint8_t tanr(int x, int y, int w, int h, int n, double *p) {
    return funcr(tanz,x,y,w,h,n,p);
}

static uint8_t sing(int x, int y, int w, int h, int n, double *p) {
    return funcg(sinz,x,y,w,h,n,p);
}

static uint8_t sinb(int x, int y, int w, int h, int n, double *p) {
    return funcb(sinz,x,y,w,h,n,p);
}

static uint8_t sinr(int x, int y, int w, int h, int n, double *p) {
    return funcr(sinz,x,y,w,h,n,p);
}

static uint8_t sininvg(int x, int y, int w, int h, int n, double *p) {
    return funcg(sininv,x,y,w,h,n,p);
}

static uint8_t sininvb(int x, int y, int w, int h, int n, double *p) {
    return funcb(sininv,x,y,w,h,n,p);
}

static uint8_t sininvr(int x, int y, int w, int h, int n, double *p) {
    return funcr(sininv,x,y,w,h,n,p);
}

static uint8_t expinvg(int x, int y, int w, int h, int n, double *p) {
    return funcg(expinv,x,y,w,h,n,p);
}

static uint8_t expinvb(int x, int y, int w, int h, int n, double *p) {
    return funcb(expinv,x,y,w,h,n,p);
}

static uint8_t expinvr(int x, int y, int w, int h, int n, double *p) {
    return funcr(expinv,x,y,w,h,n,p);
}

static uint8_t poly5g(int x, int y, int w, int h, int n, double *p) {
    return funcg(poly5,x,y,w,h,n,p);
}

static uint8_t poly5b(int x, int y, int w, int h, int n, double *p) {
    return funcb(poly5,x,y,w,h,n,p);
}

static uint8_t poly5r(int x, int y, int w, int h, int n, double *p) {
    return funcr(poly5,x,y,w,h,n,p);
}

static uint8_t poly4g(int x, int y, int w, int h, int n, double *p) {
    return funcg(poly4,x,y,w,h,n,p);
}

static uint8_t poly4b(int x, int y, int w, int h, int n, double *p) {
    return funcb(poly4,x,y,w,h,n,p);
}

static uint8_t poly4r(int x, int y, int w, int h, int n, double *p) {
    return funcr(poly4,x,y,w,h,n,p);
}

static uint8_t poly3g(int x, int y, int w, int h, int n, double *p) {
    return funcg(poly3,x,y,w,h,n,p);
}

static uint8_t poly3b(int x, int y, int w, int h, int n, double *p) {
    return funcb(poly3,x,y,w,h,n,p);
}

static uint8_t poly3r(int x, int y, int w, int h, int n, double *p) {
    return funcr(poly3,x,y,w,h,n,p);
}

static uint8_t poly2g(int x, int y, int w, int h, int n, double *p) {
    return funcg(poly2,x,y,w,h,n,p);
}

static uint8_t poly2b(int x, int y, int w, int h, int n, double *p) {
    return funcb(poly2,x,y,w,h,n,p);
}

static uint8_t poly2r(int x, int y, int w, int h, int n, double *p) {
    return funcr(poly2,x,y,w,h,n,p);
}

static uint8_t poly1g(int x, int y, int w, int h, int n, double *p) {
    return funcg(poly1,x,y,w,h,n,p);
}

static uint8_t poly1b(int x, int y, int w, int h, int n, double *p) {
    return funcb(poly1,x,y,w,h,n,p);
}

static uint8_t poly1r(int x, int y, int w, int h, int n, double *p) {
    return funcr(poly1,x,y,w,h,n,p);
}

static uint8_t idrgbg(int x, int y, int w, int h, int n, double *p) {
    double complex z = x + y * I;
    return hsv_g(carg(z),1, cabs(z)/(h/2));
}

static uint8_t idrgbb(int x, int y, int w, int h, int n, double *p) {
    double complex z = x + y * I;
    return hsv_b(carg(z),1, cabs(z)/(h/2));
}

static uint8_t idrgbr(int x, int y, int w, int h, int n, double *p) {
    double complex z = x + y * I;
    return hsv_r(carg(z),1, cabs(z)/(h/2));
}

static uint8_t sqg(int x, int y, int w, int h, int n, double *p) {
    return funcg(poly0,x,y,w,h,n,p);
}

static uint8_t sqb(int x, int y, int w, int h, int n, double *p) {
    return funcb(poly0,x,y,w,h,n,p);
}

static uint8_t sqr(int x, int y, int w, int h, int n, double *p) {
    return funcr(poly0,x,y,w,h,n,p);
}


static Func funcs[] = {
    { "zero", zero },
    { "sq", sq },
    { "idre", idre },
    { "idimg", idimg },
    { "idabs", idabs },
    { "poly1re", poly1re },
    { "poly1img", poly1img },
    { "poly1abs", poly1abs },
    { "poly5r", poly5r },
    { "poly5g", poly5g },
    { "poly5b", poly5b },
    { "poly4r", poly4r },
    { "poly4g", poly4g },
    { "poly4b", poly4b },
    { "poly3r", poly3r },
    { "poly3g", poly3g },
    { "poly3b", poly3b },
    { "poly2r", poly2r },
    { "poly2g", poly2g },
    { "poly2b", poly2b },
    { "poly1r", poly1r },
    { "poly1g", poly1g },
    { "poly1b", poly1b },
    { "expinvr", expinvr },
    { "expinvg", expinvg },
    { "expinvb", expinvb },
    { "sininvr", sininvr },
    { "sininvg", sininvg },
    { "sininvb", sininvb },
    { "tanr", tanr },
    { "tang", tang },
    { "tanb", tanb },
    { "sinr", sinr },
    { "sing", sing },
    { "sinb", sinb },
    { "idrgbr", idrgbr },
    { "idrgbg", idrgbg },
    { "idrgbb", idrgbb },
    { "sqr", sqr },
    { "sqg", sqg },
    { "sqb", sqb },
    {NULL, NULL},
};

static uint8_t (*getFunc(const char *name))(int,int,int,int,int,double*) {
    int k=0;
    while(funcs[k].name != NULL) {
        if(!strcmp(name, funcs[k].name)) {
            return funcs[k].f;
        }
        k++;
    }
    return NULL;
}
/* ================================== FILTER ========================================== */


#define OFFSET(x) offsetof(CGenContext, x)
#define FLAGS AV_OPT_FLAG_VIDEO_PARAM|AV_OPT_FLAG_FILTERING_PARAM

static const AVOption cgen_options[] = {
    { "f1",  "f1",   OFFSET(f1),  AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "p11","p11",   OFFSET(p1[0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p12","p12",   OFFSET(p1[1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p13","p13",   OFFSET(p1[2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p14","p14",   OFFSET(p1[3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p15","p15",   OFFSET(p1[4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p16","p16",   OFFSET(p1[5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p17","p17",   OFFSET(p1[6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p18","p18",   OFFSET(p1[7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p19","p19",   OFFSET(p1[8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p1a","p1a",   OFFSET(p1[9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "f2",  "f2",   OFFSET(f2),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "cp2","cp2",   OFFSET(cp2),   AV_OPT_TYPE_INT,    {.i64=0},    0,        3,        FLAGS },
    { "p21","p21",   OFFSET(p2[0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p22","p22",   OFFSET(p2[1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p23","p23",   OFFSET(p2[2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p24","p24",   OFFSET(p2[3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p25","p25",   OFFSET(p2[4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p26","p26",   OFFSET(p2[5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p27","p27",   OFFSET(p2[6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p28","p28",   OFFSET(p2[7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p29","p29",   OFFSET(p2[8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p2a","p2a",   OFFSET(p2[9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "f3",  "f3",   OFFSET(f3),  AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "cp3","cp3",   OFFSET(cp3),   AV_OPT_TYPE_INT,    {.i64=0},    0,        2,        FLAGS },
    { "p31","p31",   OFFSET(p3[0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p32","p32",   OFFSET(p3[1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p33","p33",   OFFSET(p3[2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p34","p34",   OFFSET(p3[3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p35","p35",   OFFSET(p3[4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p36","p36",   OFFSET(p3[5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p37","p37",   OFFSET(p3[6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p38","p38",   OFFSET(p3[7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p39","p39",   OFFSET(p3[8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p3a","p3a",   OFFSET(p3[9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "nf",  "nf",   OFFSET(nf),  AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf1","nf1",   OFFSET(np[0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "nf2","nf2",   OFFSET(np[1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "nf3","nf3",   OFFSET(np[2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "nf4","nf4",   OFFSET(np[3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "nf5","nf5",   OFFSET(np[4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "nf6","nf6",   OFFSET(np[5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "nf7","nf7",   OFFSET(np[6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "nf8","nf8",   OFFSET(np[7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "nf9","nf9",   OFFSET(np[8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "nfa","nfa",   OFFSET(np[9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "wf",  "wf",   OFFSET(wf),  AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "wf1","wf1",   OFFSET(wp[0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "wf2","wf2",   OFFSET(wp[1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "wf3","wf3",   OFFSET(wp[2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "wf4","wf4",   OFFSET(wp[3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "wf5","wf5",   OFFSET(wp[4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "wf6","wf6",   OFFSET(wp[5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "wf7","wf7",   OFFSET(wp[6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "wf8","wf8",   OFFSET(wp[7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "wf9","wf9",   OFFSET(wp[8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "wfa","wfa",   OFFSET(wp[9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "rgb","rgb",OFFSET(is_rgb), AV_OPT_TYPE_INT,    {.i64=0},    0,        1,        FLAGS },
    { "ofs","ofs",OFFSET(offset), AV_OPT_TYPE_INT,    {.i64=0},    0,     INT_MAX,        FLAGS },
    {NULL},
};

AVFILTER_DEFINE_CLASS(cgen);

static av_cold int cgen_init(AVFilterContext *ctx)
{
    int k;
    CGenContext *cgen = ctx->priv;
    av_log(ctx, AV_LOG_INFO, "rgb=%d ofs=%d \n", cgen->is_rgb, cgen->offset);
    if(cgen->wf) {
        cgen->wfunc = getNFunc(cgen->wf);
        if(!cgen->wfunc) {
            av_log(ctx, AV_LOG_WARNING, "function for wf not found %s\n", cgen->wf);
            cgen->wfunc = constant;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for wf is %s %f %f %f %f %f %f\n", cgen->wf, cgen->wp[0], cgen->wp[1], cgen->wp[2], cgen->wp[3], cgen->wp[4], cgen->wp[5]);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for nf\n");
        cgen->nfunc = constant;
    }

    if(cgen->nf) {
        cgen->nfunc = getNFunc(cgen->nf);
        if(!cgen->nfunc) {
            av_log(ctx, AV_LOG_WARNING, "function for nf not found %s\n", cgen->nf);
            cgen->nfunc = idn;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for nf is %s %f %f %f %f %f %f\n", cgen->nf, cgen->np[0], cgen->np[1], cgen->np[2], cgen->np[3], cgen->np[4], cgen->np[5]);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for nf\n");
        cgen->nfunc = idn;
    }

    if(cgen->f1) {
        av_log(ctx, AV_LOG_DEBUG, "F1\n");
        cgen->funcs[0] = getFunc(cgen->f1);
        if(!cgen->funcs[0]) {
            av_log(ctx, AV_LOG_DEBUG, "GET FUNC2 F1\n");
            cgen->funcs2[0] = getFunc2(cgen->f1);
            if(!cgen->funcs2[0]) {
                av_log(ctx, AV_LOG_WARNING, "function for f1 not found %s\n", cgen->f1);
                cgen->funcs[0] = zero;
            } else {
                av_log(ctx, AV_LOG_DEBUG, "FUNC2 found\n");
            }
        }

        if(cgen->funcs2[0] || cgen->funcs[0]) {
            av_log(ctx, AV_LOG_INFO, "function for f1 is %s %f %f %f %f %f %f\n", cgen->f1, cgen->p1[0], cgen->p1[1], cgen->p1[2], cgen->p1[3], cgen->p1[4], cgen->p1[5]);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for f1\n");
        cgen->funcs[0] = zero;
    }

    if(cgen->f2) {
        cgen->funcs[1] = getFunc(cgen->f2);
        if(!cgen->funcs[1]) {
            cgen->funcs2[1] = getFunc2(cgen->f2);
            if(!cgen->funcs2[1]) {
                av_log(ctx, AV_LOG_WARNING, "function not found %s\n", cgen->f2);
                cgen->funcs[1] = zero;
            }
        } 
        
        if(cgen->funcs2[1] || cgen->funcs[1])  {
            if(cgen->cp2 == 1) {
                memcpy(cgen->p2,cgen->p1,sizeof(double)*PARAMSIZE);
            } else if(cgen->cp2 == 3) {
                memcpy(cgen->p2,cgen->p3,sizeof(double)*PARAMSIZE);
            }
            av_log(ctx, AV_LOG_INFO, "function for f2 is %s %f %f %f %f %f %f\n", cgen->f2, cgen->p2[0], cgen->p2[1], cgen->p2[2], cgen->p2[3], cgen->p2[4], cgen->p2[5]);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for f2\n");
        cgen->funcs[1] = zero;
    }

    if(cgen->f3) {
        cgen->funcs[2] = getFunc(cgen->f3);
        if(!cgen->funcs[2]) {
            cgen->funcs2[2] = getFunc2(cgen->f3);
            if(!cgen->funcs2[2]) {
                av_log(ctx, AV_LOG_WARNING, "function not found %s\n", cgen->f3);
                cgen->funcs[2] = zero;
            }
        } 
        
        if(cgen->funcs2[2] || cgen->funcs[2]) {
            if(cgen->cp3 == 1) {
                memcpy(cgen->p3,cgen->p1,sizeof(double)*PARAMSIZE);
            } else if(cgen->cp2 == 2) {
                memcpy(cgen->p3,cgen->p2,sizeof(double)*PARAMSIZE);
            }
            av_log(ctx, AV_LOG_INFO, "function for f3 is %s %f %f %f %f %f %f\n", cgen->f3, cgen->p3[0], cgen->p3[1], cgen->p3[2], cgen->p3[3], cgen->p3[4], cgen->p3[5]);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for f3\n");
        cgen->funcs[2] = zero;
    }
    last = (double complex **)malloc(1920 * sizeof(double complex *));
    for(k=0;k<1920;k++) {
        last[k] = (double complex *)malloc(1080 * sizeof(double complex));
    }

    return 0;
}

static int cgen_query_formats(AVFilterContext *ctx)
{
    CGenContext *cgen = ctx->priv;
    static const enum AVPixelFormat yuv_pix_fmts[] = {
        AV_PIX_FMT_YUV444P,  AV_PIX_FMT_YUV422P,  AV_PIX_FMT_YUV420P,
        AV_PIX_FMT_YUV411P,  AV_PIX_FMT_YUV410P,  AV_PIX_FMT_YUV440P,
        AV_PIX_FMT_YUVA444P, AV_PIX_FMT_YUVA422P, AV_PIX_FMT_YUVA420P,
        AV_PIX_FMT_GRAY8,
        AV_PIX_FMT_YUV444P9,  AV_PIX_FMT_YUV422P9,  AV_PIX_FMT_YUV420P9,
        AV_PIX_FMT_YUVA444P9, AV_PIX_FMT_YUVA422P9, AV_PIX_FMT_YUVA420P9,
        AV_PIX_FMT_YUV444P10,  AV_PIX_FMT_YUV422P10,  AV_PIX_FMT_YUV420P10,
        AV_PIX_FMT_YUV440P10,
        AV_PIX_FMT_YUVA444P10, AV_PIX_FMT_YUVA422P10, AV_PIX_FMT_YUVA420P10,
        AV_PIX_FMT_GRAY9, AV_PIX_FMT_GRAY10,
        AV_PIX_FMT_YUV444P12,  AV_PIX_FMT_YUV422P12,  AV_PIX_FMT_YUV420P12,
        AV_PIX_FMT_GRAY12, AV_PIX_FMT_GRAY14,
        AV_PIX_FMT_YUV444P14,  AV_PIX_FMT_YUV422P14,  AV_PIX_FMT_YUV420P14,
        AV_PIX_FMT_YUV444P16,  AV_PIX_FMT_YUV422P16,  AV_PIX_FMT_YUV420P16,
        AV_PIX_FMT_YUVA444P16, AV_PIX_FMT_YUVA422P16, AV_PIX_FMT_YUVA420P16,
        AV_PIX_FMT_GRAY16,
        AV_PIX_FMT_NONE
    };
    static const enum AVPixelFormat rgb_pix_fmts[] = {
        AV_PIX_FMT_GBRP, AV_PIX_FMT_GBRAP,
        AV_PIX_FMT_GBRP9,
        AV_PIX_FMT_GBRP10, AV_PIX_FMT_GBRAP10,
        AV_PIX_FMT_GBRP12, AV_PIX_FMT_GBRAP12,
        AV_PIX_FMT_GBRP14,
        AV_PIX_FMT_GBRP16, AV_PIX_FMT_GBRAP16,
        AV_PIX_FMT_NONE
    };
    AVFilterFormats *fmts_list;

    if (cgen->is_rgb) {
        fmts_list = ff_make_format_list(rgb_pix_fmts);
    } else
        fmts_list = ff_make_format_list(yuv_pix_fmts);
    if (!fmts_list)
        return AVERROR(ENOMEM);
    return ff_set_common_formats(ctx, fmts_list);
}

static int cgen_config_props(AVFilterLink *inlink)
{
    CGenContext *cgen = inlink->dst->priv;
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(inlink->format);

    av_assert0(desc);

    cgen->hsub = desc->log2_chroma_w;
    cgen->vsub = desc->log2_chroma_h;
    cgen->bps = desc->comp[0].depth;
    cgen->planes = desc->nb_components;
    return 0;
}

typedef struct ThreadData {
    int height;
    int width;
    int plane;
    int linesize;
    int n;
    AVFrame *in;
} ThreadData;

static int slice_cgen_filter(AVFilterContext *ctx, void *arg, int jobnr, int nb_jobs)
{
    CGenContext *cgen = ctx->priv;
    ThreadData *td = arg;
    const int height = td->height;
    const int width = td->width;
    const int plane = td->plane;
    const int linesize = td->linesize;
    const int slice_start = (height *  jobnr) / nb_jobs;
    const int slice_end = (height * (jobnr+1)) / nb_jobs;
    int x, y;
    uint8_t *ptr;
    
    if(plane == 0) {
        for (y = slice_start; y < slice_end; y++) {
            ptr = td->in->data[plane] + linesize * y;

            for (x = 0; x < width; x++) {
                ptr[x] = cgen->funcs2[0]?cgen->funcs2[0](x-width/2,y-height/2,width,height,td->n + cgen->offset, cgen):cgen->funcs[0](x-width/2,y-height/2,width,height,td->n + cgen->offset,cgen->p1);
            }
        }
    }
    if(plane == 1) {
        for (y = slice_start; y < slice_end; y++) {
            ptr = td->in->data[plane] + linesize * y;

            for (x = 0; x < width; x++) {
                ptr[x] = cgen->funcs2[1]?cgen->funcs2[1](x-width/2,y-height/2,width,height,td->n + cgen->offset, cgen):cgen->funcs[1](x-width/2,y-height/2,width,height,td->n + cgen->offset,cgen->p2);
            }
        }
    }
    if(plane == 2) {
        for (y = slice_start; y < slice_end; y++) {
            ptr = td->in->data[plane] + linesize * y;

            for (x = 0; x < width; x++) {
                ptr[x] = cgen->funcs2[2]?cgen->funcs2[2](x-width/2,y-height/2,width,height,td->n + cgen->offset, cgen):cgen->funcs[2](x-width/2,y-height/2,width,height,td->n + cgen->offset,cgen->p3);
            }
        }
    }
    return 0;
}

static int cgen_filter_frame(AVFilterLink *inlink, AVFrame *in)
{
    int plane;
    AVFilterContext *ctx = inlink->dst;
    const int nb_threads = ff_filter_get_nb_threads(ctx);
    CGenContext *cgen = ctx->priv;
    AVFilterLink *outlink = inlink->dst->outputs[0];
    av_frame_make_writable(in);
    double n = inlink->frame_count_out;
    //double t = in->pts == AV_NOPTS_VALUE ? NAN : in->pts * av_q2d(inlink->time_base);
    

    for (plane = 0; plane < cgen->planes && in->data[plane]; plane++) {
        const int width = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(inlink->w, cgen->hsub) : inlink->w;
        const int height = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(inlink->h, cgen->vsub) : inlink->h;
        const int linesize = in->linesize[plane];
        ThreadData td;
        td.width = width;
        td.height = height;
        td.plane = plane;
        td.linesize = linesize;
        td.in = in;
        td.n = n;
        ctx->internal->execute(ctx, slice_cgen_filter, &td, NULL, FFMIN(height, nb_threads));
    }

    return ff_filter_frame(outlink, in);
}

static av_cold void cgen_uninit(AVFilterContext *ctx)
{
    
    av_log(ctx, AV_LOG_INFO, "uninit\n");
    int k;
    for(k=0;k<1920;k++) {
        free(last[k]);
    }
    free(last);
}

static const AVFilterPad cgen_inputs[] = {
    {
        .name         = "default",
        .type         = AVMEDIA_TYPE_VIDEO,
        .config_props = cgen_config_props,
        .filter_frame = cgen_filter_frame,
    },
    { NULL }
};

static const AVFilterPad cgen_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO,
    },
    { NULL }
};

AVFilter ff_vf_cgen = {
    .name          = "cgen",
    .description   = NULL_IF_CONFIG_SMALL("Apply generic equation to each pixel."),
    .priv_size     = sizeof(CGenContext),
    .init          = cgen_init,
    .uninit        = cgen_uninit,
    .query_formats = cgen_query_formats,
    .inputs        = cgen_inputs,
    .outputs       = cgen_outputs,
    .priv_class    = &cgen_class,
    .flags         = AVFILTER_FLAG_SUPPORT_TIMELINE_GENERIC | AVFILTER_FLAG_SLICE_THREADS,
};
