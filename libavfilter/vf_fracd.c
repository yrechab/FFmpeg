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
#include <quadmath.h>
#include <stdio.h>
#include <string.h>
#include "libavutil/avassert.h"
#include "libavutil/avstring.h"
#include "libavutil/opt.h"
#include "libavutil/pixdesc.h"
#include "libavutil/genutil.h"
#include "libavutil/plot.h"
#include "internal.h"
#include <complex.h>
//#define PRECLD
#define PRECD
//#define PREC128
#ifdef PRECLD
#define MYFLT long double
#define MYCOMPLEX long double complex
#define CA(z_,r_,i_) (z_) = (r_) + I * (i_)  
#define MYSIN sinl
#define MYCSIN csinl
#define MYSINH sinhl
#define MYCSINH csinhl
#define MYCOS cosl
#define MYCCOS ccosl
#define MYCOSH coshl
#define MYCCOSH ccoshl
#define MYTAN tanl
#define MYCTAN ctanl
#define MYATAN atanl
#define MYCATAN catanl
#define MYEXP expl
#define MYCEXP cexpl
#define MYSQRT sqrtl
#define MYCSQRT csqrtl
#define MYLOG logl
#define MYCLOG clogl
#define MYPOW powl
#define MYCPOW cpowl
#define MYREAL creall
#define MYIMAG cimagl
#define MYFABS fabsl
#define MYCABS cabsl
#define MYSTRTOF(s_) strtold((s_),0)
#endif

#ifdef PRECD
#define MYFLT double
#define MYCOMPLEX double complex
#define CA(z_,r_,i_) (z_) = (r_) + I * (i_)  
#define MYSIN sin
#define MYCSIN csin
#define MYSINH sinh
#define MYCSINH csinh
#define MYCOS cos
#define MYCCOS ccos
#define MYCOSH cosh
#define MYCCOSH ccosh
#define MYTAN tan
#define MYCTAN ctan
#define MYATAN atan
#define MYCATAN catan
#define MYEXP exp
#define MYCEXP cexp
#define MYSQRT sqrt
#define MYCSQRT csqrt
#define MYLOG log
#define MYCLOG clog
#define MYPOW pow
#define MYCPOW cpow
#define MYREAL creal
#define MYIMAG cimag
#define MYFABS fabs
#define MYCABS cabs
#define MYSTRTOF(s_) atof((s_))
#endif

#ifdef PREC128
#define MYFLT __float128
#define MYCOMPLEX __complex128
#define CA(z_, r_, i_) {__real__(z_) = (r_); __imag__(z_) = (i_);}
#define MYSIN sinq
#define MYCSIN csinq
#define MYSINH sinhq
#define MYCSINH csinhq
#define MYCOS cosq
#define MYCCOS ccosq
#define MYCOSH coshq
#define MYCCOSH ccoshq
#define MYTAN tanq
#define MYCTAN ctanq
#define MYATAN atanq
#define MYCATAN catanq
#define MYEXP expq
#define MYCEXP cexpq
#define MYSQRT sqrtq
#define MYCSQRT csqrtq
#define MYLOG logq
#define MYCLOG clogq
#define MYPOW powq
#define MYCPOW cpowq
#define MYREAL crealq
#define MYIMAG cimagq
#define MYFABS fabsq
#define MYCABS cabsq
#define MYSTRTOF(s_) strtoflt128((s_),0)
#endif


typedef struct FracDFuncParams {
    MYFLT x;
    MYFLT y;
    double w;
    double h;
    MYFLT fx;
    MYFLT fy;
    double rot;
    double sat;
    double fac;
    double *p;
    int ifcmode;
    MYCOMPLEX (*ifunc)(MYCOMPLEX,MYFLT*);
    MYFLT ip[40];
    double (*cfunc[3])(double,double*);
    double cp[3][10];
    int cmod;
    int is_rgb;
} FracDFuncParams;

typedef struct FracDContext {
    const AVClass *class;
    const char *f;
    void (*ffunc)(FracDFuncParams *, AVFrame*,int,int,int);
    double p[40];
    const char *ifc;
    MYCOMPLEX (*ifunc)(MYCOMPLEX,MYFLT*);
    MYFLT ip[40];
    int ifcmode;
    const char *nf[10];
    MYFLT np[10][100];
    int nmod[10]; // modulo n
    MYFLT (*nfunc[10])(int,MYFLT*);
    const char *cf[3];
    double (*cfunc[3])(double,double*);
    double cp[3][10];
    int cmod;
    MYFLT x;
    MYFLT y;
    MYFLT fx;
    MYFLT fy;
    char *str_x;
    char *str_y;
    char *str_fx;
    char *str_fy;
    double w;
    double h;
    double rot;
    double sat;
    double fac;
    int hsub, vsub;             ///< chroma subsampling
    int planes;                 ///< number of planes
    int is_rgb;
    int bps;
    int offset;
    int dbg;
    char *rf[20];

} FracDContext;


#define OFFSET(x) offsetof(FracDContext, x)
#define FLAGS AV_OPT_FLAG_VIDEO_PARAM|AV_OPT_FLAG_FILTERING_PARAM

static const AVOption fracd_options[] = {
    { "f",  "f",   OFFSET(f),  AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    //{ "x","x",       OFFSET(x),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    //{ "y","y",       OFFSET(y),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    //{ "fx","fx",   OFFSET(fx),   AV_OPT_TYPE_DOUBLE, {.dbl =  1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    //{ "fy","fy",   OFFSET(fy),   AV_OPT_TYPE_DOUBLE, {.dbl =  1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "x","x",       OFFSET(str_x),     AV_OPT_TYPE_STRING, {.str = "0"}, CHAR_MIN, CHAR_MAX,  FLAGS },
    { "y","y",       OFFSET(str_y),     AV_OPT_TYPE_STRING, {.str = "0"}, CHAR_MIN, CHAR_MAX,  FLAGS },
    { "fx","fx",   OFFSET(str_fx),   AV_OPT_TYPE_STRING, {.str =  "1"}, CHAR_MIN, CHAR_MAX,  FLAGS },
    { "fy","fy",   OFFSET(str_fy),   AV_OPT_TYPE_STRING, {.str =  "1"}, CHAR_MIN, CHAR_MAX,  FLAGS },
    { "rot","rot",   OFFSET(rot),   AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "sat","sat",   OFFSET(sat),   AV_OPT_TYPE_DOUBLE, {.dbl = 1.0}, 0,  1,  FLAGS },
    { "fac","fac",   OFFSET(fac),   AV_OPT_TYPE_DOUBLE, {.dbl = 1.0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "ifc",  "ifc",   OFFSET(ifc),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "ifcmode","ifcmode",OFFSET(ifcmode), AV_OPT_TYPE_INT,    {.i64=0},    0,        10,        FLAGS },
    { "nf0",  "nf0",     OFFSET(nf[0]),     AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf1",  "nf1",     OFFSET(nf[1]),     AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf2",  "nf2",     OFFSET(nf[2]),     AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf3",  "nf3",     OFFSET(nf[3]),     AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf4",  "nf4",     OFFSET(nf[4]),     AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf5",  "nf5",     OFFSET(nf[5]),     AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf6",  "nf6",     OFFSET(nf[6]),     AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf7",  "nf7",     OFFSET(nf[7]),     AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf8",  "nf8",     OFFSET(nf[8]),     AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf9",  "nf9",     OFFSET(nf[9]),     AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nm0","nm0",   OFFSET(nmod[0]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },
    { "nm1","nm1",   OFFSET(nmod[1]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },
    { "nm2","nm2",   OFFSET(nmod[2]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },
    { "nm3","nm3",   OFFSET(nmod[3]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },
    { "nm4","nm4",   OFFSET(nmod[4]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },
    { "nm5","nm5",   OFFSET(nmod[5]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },
    { "nm6","nm6",   OFFSET(nmod[6]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },
    { "nm7","nm7",   OFFSET(nmod[7]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },
    { "nm8","nm8",   OFFSET(nmod[8]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },
    { "nm9","nm9",   OFFSET(nmod[9]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },

    { "cf0",  "cf0",   OFFSET(cf[0]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "cf1",  "cf1",   OFFSET(cf[1]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "cf2",  "cf2",   OFFSET(cf[2]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "cmod","cmod",OFFSET(cmod), AV_OPT_TYPE_INT,    {.i64=0},    0,        2,        FLAGS },
 
    { "rgb","rgb",OFFSET(is_rgb), AV_OPT_TYPE_INT,    {.i64=1},    0,        1,        FLAGS },
    { "dbg","dbg",OFFSET(dbg), AV_OPT_TYPE_INT,    {.i64=0},    0,        1,        FLAGS },
    { "ofs","ofs",OFFSET(offset), AV_OPT_TYPE_INT,    {.i64=0},    0,     INT_MAX,        FLAGS },
    { "r0", "r0",  OFFSET(rf[0]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "r1", "r1",  OFFSET(rf[1]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "r2", "r2",  OFFSET(rf[2]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "r3", "r3",  OFFSET(rf[3]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "r4", "r4",  OFFSET(rf[4]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "r5", "r5",  OFFSET(rf[5]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "r6", "r6",  OFFSET(rf[6]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "r7", "r7",  OFFSET(rf[7]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "r8", "r8",  OFFSET(rf[8]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "r9", "r9",  OFFSET(rf[9]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },

    {NULL},
};

AVFILTER_DEFINE_CLASS(fracd);


static MYCOMPLEX rotate(MYCOMPLEX z, MYFLT phi) {
    return phi!=0.0?(MYCOS(phi) + MYSIN(phi)*I) * z:z;
}
/* ============================= NFUNCS *******************************************/

/*** spline ***/
static MYFLT p3(int n, MYFLT a,MYFLT b,MYFLT c,MYFLT d) {
    return a + b*n + c*n*n + d*n*n*n;
}

static MYFLT pcubic(int n, MYFLT *p) {
    int len = p[0];
    int end = p[1];
    int k,index;
    for(k=1;k<len;k++) {
        if(n<p[k*5+2]) {
            index = (k-1)*5 +2;
            return p3(n-p[index],p[index+1],p[index+2],p[index+3],p[index+4]);
        } 
    }
    if(n<end) {
        index = (len-1)*5+2;
        return p3(n-p[index],p[index+1],p[index+2],p[index+3],p[index+4]);
    } else {
        index = (len-1)*5+2;
        return p3(end-p[index],p[index+1],p[index+2],p[index+3],p[index+4]);
    }
}

static MYFLT pspline(int n, MYFLT *p) {
    int len;
    MYFLT *out,ret;

    len = p[0];
    out = calloc(len*5+2, sizeof(MYFLT));
    av_genutil_spline(p,out);
    ret =  pcubic(n,out);
    free(out);
    return ret;
}

static MYFLT npoly(int n, MYFLT *p) {
    MYFLT x;
    x = (MYFLT)n * p[0];
    return p[1] + p[2]*x + p[3]*x+x + p[4]*x*x*x + p[5]*x*x*x*x;
}

static MYFLT nsin(int n, MYFLT *p) {
    return (MYSIN(p[0]*n - M_PI/2 + p[3]))*p[1]+p[2]; 
}

static MYFLT rsin(int n, MYFLT *p) {
    return (MYSIN(p[0]*n + p[3]))*(p[1]==0.0?1:p[1])+p[2]; 
}

static MYFLT rcos(int n, MYFLT *p) {
    return (MYCOS(p[0]*n + p[3]))*(p[1]==0.0?1:p[1])+p[2]; 
}

static MYFLT linsin(int n, MYFLT *p) {
    return p[0]*n+(MYSIN(p[1]*n + p[3]))*p[2]; 
}

static MYFLT idn(int n, MYFLT *p) {
    return n;
}

static MYFLT nexp(int n, MYFLT *p) {
    return p[2] + p[1]*(MYEXP(p[0]*n-p[3]));
}

static MYFLT nexpinv(int n, MYFLT *p) {
    return p[2] + p[1]/(MYEXP(p[0]*n-p[3]));
}

static MYFLT nexp2(int n, MYFLT *p) {
    return p[2] + p[1]*(MYEXP(p[0]*n*n-p[3]));
}

static MYFLT constant(int n, MYFLT *p) {
    return p[0];
}

static MYFLT sqn(int n, MYFLT *p) {
    MYFLT x;
    x = (MYFLT)n*p[0];
    return p[1] * (MYSQRT(x+0.25)-0.5);
}

static MYFLT ln(int n, MYFLT *p) {
    MYFLT x;
    x = (MYFLT)n*p[0];
    return p[1] * MYLOG(x+1);
}

static MYFLT lin(int n, MYFLT *p) {
    return p[0] + p[1]*n;
}

static MYFLT inv(int n, MYFLT *p) {
    MYFLT x;
    x = (MYFLT)n*p[0];
    return p[1]*(1/(p[2]*((MYFLT)x+1)));
}

static MYFLT pval(int n, MYFLT *p) {
    int k;
    if(n<p[0]) return 0;
    if(p[NMAXPARAMS-2]&&(n>p[NMAXPARAMS-2])) return p[NMAXPARAMS-1];
    for(k=0;k<NMAXPARAMS-2;k+=2) {
        if(k>0 && p[k]==0) return p[k-1]; 
        if(n>p[k] && n<=p[k+2])
            return p[k+1];
    }
    return 0;
}

static MYFLT pchain(int n, MYFLT *p) {
    int len,section;
    MYFLT x;
    len = p[0];
    section = n/len + 1;
    if(section > NMAXPARAMS) return 0;
    x = (MYFLT)(n%len)/(MYFLT)len;
    return p[section] + (p[section+1]-p[section])*x;
}

static MYFLT gauss(int n, MYFLT *p) {
    MYFLT m,fac,x;
    m = p[1]==0?1:p[1];
    fac = p[2]==0?0.5:p[2];
    x = (MYFLT)n*p[0];
    return m * MYEXP(-fac * MYPOW(x-p[3],2)); 
}

typedef struct NFunc {
    const char *name;
    MYFLT (*f)(int n, MYFLT *p);
} NFunc;

static NFunc nfuncs[] = {
    {"idn",idn},
    {"constant",constant},
    {"lin",lin},
    {"linsin",linsin},
    {"poly",npoly},
    {"pchain",pchain},
    {"pcubic",pcubic},
    {"pspline",pspline},
    {"pval",pval},
    {"sin",nsin},
    {"rsin",rsin},
    {"rcos",rcos},
    {"inv",inv},
    {"gauss",gauss},
    {"sq",sqn},
    {"ln",ln},
    {"exp",nexp},
    {"expinv",nexpinv},
    {"exp2",nexp2},
    {NULL,NULL}
};

static MYFLT (*get_nfunc(const char *name))(int,MYFLT*) {
    int k=0;
    while(nfuncs[k].name) {
        if(!strcmp(name, nfuncs[k].name)) {
            return nfuncs[k].f;
        }
        k++;
    }
    return NULL;
}

static void parse_nfunc(const char *nf, MYFLT *p, MYFLT (**f)(int,MYFLT*)) {
    char *saveptr,*token,*str;
    int j;
    str = strdup(nf);
    for (j = 0; ; j++, str = NULL) {
        token = av_strtok(str, " |", &saveptr);
        if (token == NULL)
            break;
        if(j==0) {
            *f = get_nfunc(token);
        } else {
            p[j-1] = MYSTRTOF(token);
        }
    }
    free(str);
}



// ********************* inner functions **************************/
typedef struct IFunc {
    const char *name;
    MYCOMPLEX (*f)(MYCOMPLEX,MYFLT*);
} IFunc;

static MYCOMPLEX izero(MYCOMPLEX z, MYFLT *p) {
    return 0;
}

static MYCOMPLEX hopalong(MYCOMPLEX z, MYFLT *p) {
    MYFLT a,b,c,x0,y0,x,y;
    MYCOMPLEX ret;
    a = p[0];
    b = p[1];
    c = p[2];
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    x =  y0 - SIGN(x0)*MYSQRT(MYFABS(b*x0-c));
    y = a - x0;
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX hopalongsin(MYCOMPLEX z, MYFLT *p) {
    MYFLT a,b,c,d,x0,y0,xx,x,y;
    MYCOMPLEX ret;
    a = p[0];
    b = p[1];
    c = p[2];
    d = p[3];
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    xx = MYFABS(b*x0-c);
    x =  y0 - SIGN(x0)*MYSQRT(xx + MYSIN(d*xx));
    y = a - x0;
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX hopalongp(MYCOMPLEX z, MYFLT *p) {
    MYFLT a,b,c,d,x0,y0,x,y;
    MYCOMPLEX ret;
    a = p[0];
    b = p[1];
    c = p[2];
    d = p[3];
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    x =  y0 - SIGN(x0)*MYPOW(MYFABS(b*x0-c),d);
    y = a - x0;
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX hopalonglog(MYCOMPLEX z, MYFLT *p) {
    MYFLT a,b,c,x0,y0,x,y;
    MYCOMPLEX ret;
    a = p[0];
    b = p[1];
    c = p[2];
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    x =  y0 - SIGN(x0)*MYLOG(MYFABS(b*x0-c));
    y = a - x0;
    CA(ret,x,y);
    return ret;
}


static MYCOMPLEX hopalong_s(MYCOMPLEX z, MYFLT *p) {
    MYFLT a,b,c,x0,y0,x,y;
    MYCOMPLEX ret;
    a = p[0];
    b = p[1];
    c = p[2];
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    x =  y0 - MYSQRT(MYFABS(b*x0-c));
    y = a - x0;
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX hopatest(MYCOMPLEX z, MYFLT *p) {
    MYFLT a,b,c,x0,y0,x,y;
    MYCOMPLEX ret;
    a = p[0];
    b = p[1];
    c = p[2];
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    if(x0<0) x0 = -x0;
    x =  y0 - SIGN(x0)*MYSQRT(MYFABS(b*x0-c));
    y = a - x0;
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX bedhead(MYCOMPLEX z, MYFLT *p) {
    MYFLT a,b,x0,y0,x,y;
    MYCOMPLEX ret;
    a = p[0];
    b = p[1];
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    x = MYSIN(x0*y0/b)*y0+MYCOS(a*x0-y0);
    y = x0+MYSIN(y0)/b;
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX _clifford(MYCOMPLEX z, MYFLT a, MYFLT b, MYFLT c, MYFLT d) {
    MYFLT x0,y0,x,y;
    MYCOMPLEX ret;
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    x = MYSIN(a*y0)+c*MYCOS(a*x0);
    y = MYSIN(b*x0)+d*MYCOS(b*y0);
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX clifford1(MYCOMPLEX z, MYFLT *p) {
    return _clifford(z,p[0],p[1],p[2],p[3]);
}

static MYCOMPLEX clifford2(MYCOMPLEX z, MYFLT *p) {
    return _clifford(z,p[2],p[3],p[0],p[1]);
}

static MYCOMPLEX clifford3(MYCOMPLEX z, MYFLT *p) {
    return _clifford(z,p[2],p[1],p[0],p[3]);
}

static MYCOMPLEX ginger(MYCOMPLEX z, MYFLT *p) {
    MYFLT x,y;
    MYCOMPLEX ret;
    x = p[0] - MYIMAG(z) + MYFABS(MYREAL(z));
    y = p[1] * MYREAL(z);
    CA(ret,x,y);
    return ret;
}


static MYFLT miraf(MYFLT x, MYFLT a) {
    return a*x-(1-a)*(2*x*x/(1+x*x));
}

static MYCOMPLEX mira(MYCOMPLEX z, MYFLT *p) {
    MYFLT a,b,x,y;
    MYCOMPLEX ret;
    a = p[0];
    b = p[1];
    x = b * MYIMAG(z) + miraf(MYREAL(z),a);
    y = -MYREAL(z) + miraf(x,a);
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX kaneko(MYCOMPLEX z, MYFLT *p) {
    MYFLT a,b,x,y;
    MYCOMPLEX ret;
    a = p[0];
    b = p[1];
    x = a * MYREAL(z) + (1-a)*(1-b*MYIMAG(z)*MYIMAG(z));
    y = MYREAL(z);
    CA(ret,x,y);
    return ret;
}


static MYCOMPLEX mandel(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX c;
    CA(c,p[0],p[1]);
    return z*z + c;
}

static MYCOMPLEX ship(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX c;
    MYCOMPLEX _z;
    MYFLT x,y;
    x = MYFABS(MYREAL(z));
    y = MYFABS(MYIMAG(z));
    CA(c,p[0],p[1]);
    CA(_z,x,y);
    return _z*_z + c;
}

static MYCOMPLEX exp4(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX c;
    CA(c,p[0],p[1]);
    return MYCEXP(z/(c*c*c*c));
}

static MYCOMPLEX expc(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX c;
    CA(c,p[0],p[1]);
    return MYCEXP((z*z+p[2]*z)/MYCSQRT(c*c*c));
}

static MYCOMPLEX sin_h(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX c;
    CA(c,p[0],p[1]);
    return MYCSINH(z/c);
}

static MYCOMPLEX exp3(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX c;
    CA(c,p[0],p[1]);
    return MYCEXP((z*z*z)/(c*c*c));
}

static MYCOMPLEX exp_2(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX c;
    CA(c,p[0],p[1]);
    return MYCEXP((z*z-p[2]*z)/(c*c*c));
}

static MYCOMPLEX cosinv(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX c;
    CA(c,p[0],p[1]);
    return MYCCOS(z/c);
}

static MYCOMPLEX taninv(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX c;
    CA(c,p[0],p[1]);
    return MYCTAN(z/c);
}

static MYCOMPLEX id(MYCOMPLEX z, MYFLT *p) {
    return z;
}

static MYCOMPLEX rat(MYCOMPLEX z, MYFLT *p) {
    int zg,k,ng;
    MYCOMPLEX zl,nen;
    zg = p[0];
    CA(zl,p[1],0);
    for(k = 0; k<(zg*2); k+=2) {
	MYCOMPLEX zk;
	CA(zk,p[k+2],p[k+3]);
        zl *= (z-zk);
    }
    ng = p[20];
    if(ng == 0) return zl;

    CA(nen,p[1],0);
    for(k = 20; k<(ng*2+20); k+=2) {
	MYCOMPLEX zk;
	CA(zk,p[k+2],p[k+3]);
        nen *= (z-zk);
    }
    return zl/nen;
}

static MYCOMPLEX fexp(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX ex;
    CA(ex,0,(p[2]?p[2]:1)/(1/MYPOW(MYCABS(z),p[3])));
    return p[0] + (p[1]?p[1]:1) * z * MYEXP(ex);
}

static MYCOMPLEX expr(MYCOMPLEX z, MYFLT *p) {
    return MYEXP(rat(z,p));
}

static MYCOMPLEX sinr(MYCOMPLEX z,  MYFLT *p) {
    return MYSIN(rat(z,p));
}

static MYCOMPLEX tanr(MYCOMPLEX z, MYFLT *p) {
    return MYTAN(rat(z,p));
}

static MYCOMPLEX hypno(MYCOMPLEX z, MYFLT *p) {
    MYFLT r;
    MYCOMPLEX ret;
    r = MYCABS(z);
    CA(ret,p[0]*r*MYCOS(p[1]*r),p[2]*r*MYSIN(p[3]*r));
    return ret;
}

static MYCOMPLEX hypno1(MYCOMPLEX z, MYFLT *p) {
    MYFLT r,x0,y0,phi;
    MYCOMPLEX ret;
    r = MYCABS(z);
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    if(y0<0) x0 = -x0;
    phi = MYATAN(x0/y0);
    CA(ret,p[0]*phi*r*MYCOS(p[1]*r),p[2]*r*phi*MYSIN(p[3]*r));
    return ret;
}

static MYCOMPLEX v3(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX ret;
    MYFLT x0,y0,r,x,y;
    r = MYCABS(z);
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    x = x0*MYSIN(r*r)-y0*MYCOS(r*r);
    y = x0*MYCOS(r*r)-y0*MYSIN(r*r);
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX v1(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX ret;
    MYFLT x0,y0;
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    CA(ret,MYSIN(x0),MYSIN(y0));
    return ret;
}

static MYCOMPLEX v6(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX ret;
    MYFLT x0,y0,r,phi,x,y;
    r = MYCABS(z);
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    if(y0<0) x0 = -x0;
    phi = MYATAN(x0/y0);
    x = MYSIN(phi+r);
    y = MYCOS(phi-r);
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX v8(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX ret;
    MYFLT x0,y0,r,phi,x,y;
    r = MYCABS(z);
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    if(y0<0) x0 = -x0;
    phi = MYATAN(x0/y0);
    x = (phi/p[0])*MYSIN(p[1]*r);
    y = (phi/p[2])*MYCOS(p[3]*r);
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX v9(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX ret;
    MYFLT x0,y0,r,phi,x,y;
    r = MYCABS(z);
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    if(y0<0) x0 = -x0;
    phi = MYATAN(x0/y0);
    x = (MYCOS(phi)+MYSIN(r))/r;
    y = (MYSIN(phi)-MYCOS(r))/r;
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX v10(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX ret;
    MYFLT x0,y0,r,phi,x,y;
    r = MYCABS(z);
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    if(y0<0) x0 = -x0;
    phi = MYATAN(x0/y0);
    x = MYSIN(phi)/r;
    y = MYCOS(phi)*r;
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX v11(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX ret;
    MYFLT x0,y0,r,phi,x,y;
    r = MYCABS(z);
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    if(y0<0) x0 = -x0;
    phi = MYATAN(x0/y0);
    x = MYSIN(phi)*MYCOS(r);
    y = MYCOS(phi)*MYSIN(r);
    CA(ret,x,y);
    return ret;
}
#define P3(x) (x)*(x)*(x)
#define P2(x) (x)*(x)
static MYCOMPLEX v12(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX ret;
    MYFLT p0,p1,x0,y0,r,phi,x,y;
    r = MYCABS(z);
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    if(y0<0) x0 = -x0;
    phi = MYATAN(x0/y0);
    p0 = MYSIN(phi+r);
    p1 = MYCOS(phi-r);
    x = P3(p0)+P3(p1);
    y = P3(p0)-P3(p1);
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX v612(MYCOMPLEX z, MYFLT *p) {
    MYCOMPLEX z0,ret;
    MYFLT p0,p1,x,y;
    z0 = v6(z,p);
    p0 = MYREAL(z0);
    p1 = MYIMAG(z0);
    x = p[0]*MYPOW(p0,p[1]) + p[2]*MYPOW(p1,p[3]);
    y = p[4]*MYPOW(p0,p[5]) + p[6]*MYPOW(p1,p[7]);
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX v19(MYCOMPLEX z, MYFLT *p) {
    MYFLT x0,y0,x,y,r,phi,ep;
    MYCOMPLEX ret;
    r = MYCABS(z);
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    if(y0<0) x0 = -x0;
    phi = MYATAN(x0/y0);
    ep = MYPOW(r,MYSIN(phi));
    x = ep*MYCOS(phi);
    y = ep*MYSIN(phi);
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX v20(MYCOMPLEX z, MYFLT *p) {
    MYFLT x,y,x0,y0;
    MYCOMPLEX ret;
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    x = MYCOS((M_PI)*x0)*MYCOSH(y0);
    y = MYSIN((M_PI)*x0)*MYSINH(y0);
    CA(ret,x,y);
    return ret;
}

static MYCOMPLEX v48(MYCOMPLEX z, MYFLT *p) {
    MYFLT x0,y0,r;
    MYCOMPLEX ret;
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    r = MYSQRT(1/(P2((P2(x0)-P2(y0)))));
    CA(ret,r*x0,r*y0);
    return ret;
}

static MYCOMPLEX vrnd(MYCOMPLEX z, MYFLT *p) {
    MYFLT x0,y0,x,y;
    MYCOMPLEX ret;
    double rnd;
    x0 = MYREAL(z);
    y0 = MYIMAG(z);
    rnd = (double)rand()/RAND_MAX;
    x = x0 + rnd/20;
    y = y0 + rnd/20;
    CA(ret,x,y);
    return ret;
}



static IFunc ifuncs[] = {
    { "zero", izero },
    { "hopalong", hopalong },
    { "hopalongp", hopalongp },
    { "hopalongsin", hopalongsin },
    { "hopalonglog", hopalonglog },
    { "hopatest", hopatest },
    { "hopalong_s", hopalong_s },
    { "mira", mira },
    { "kaneko", kaneko },
    { "ginger", ginger },
    { "bedhead", bedhead },
    { "cliff1", clifford1 },
    { "cliff2", clifford2 },
    { "cliff3", clifford3 },
    { "mandel", mandel },
    { "ship", ship },
    { "hypno", hypno },
    { "hypno1", hypno1 },
    { "exp4", exp4 },
    { "exp3", exp3 },
    { "exp2", exp_2 },
    { "expc", expc },
    { "sinh", sin_h },
    { "cosinv", cosinv },
    { "tan", taninv },
    { "rat", rat },
    { "ikeda", fexp },
    { "fexpr", expr },
    { "fsin", sinr },
    { "ftan", tanr},
    { "fv3", v3 },
    { "fv1", v1 },
    { "fv6", v6 },
    { "fv8", v8 },
    { "fv9", v9 },
    { "fv10", v10 },
    { "fv11", v11 },
    { "fv12", v12 },
    { "fv19", v19 },
    { "fv20", v20 },
    { "fv48", v48 },
    { "fv612", v612 },
    { "fvrnd", vrnd },
    { "id", id },
    { NULL, NULL }
};

static MYCOMPLEX (*get_ifunc(const char *name))(MYCOMPLEX, MYFLT*) {
    int k=0;
    while(ifuncs[k].name) {
        if(!strcmp(name, ifuncs[k].name)) {
            return ifuncs[k].f;
        }
        k++;
    }
    return NULL;
}

static void parse_ifunc(const char *ff, MYFLT *p, MYCOMPLEX (**f)(MYCOMPLEX, MYFLT*)) {
    char *saveptr,*token;
    char *str = strdup(ff);
    int j;
    for (j = 0; ; j++, str = NULL) {
        token = av_strtok(str, " |", &saveptr);
        if (token == NULL)
            break;
        if(j==0) {
            *f = get_ifunc(token);
        } else {
            p[j-1] = MYSTRTOF(token);
        }
    }
    free(str);
}

// ************************** hsv util **************************************/
static void get_hsv2rgb_params(double phi, double S, double V, uint8_t *h, double *p, double *q, double *t) {
    double H,f;
    H = phi<0?((180/M_PI)*phi)+360:(180/M_PI)*phi;
    *h = floor(H/60.0);
    f = H/60.0 - *h;
    *p = V*(1-S);
    *q = V*(1-S*f);
    *t = V*(1-S*(1-f)); 
}

static uint8_t hsv_r(double phi, double S, double V) {
    double p,q,t;
    uint8_t h;
    uint8_t retval = 0;

    get_hsv2rgb_params(phi, S, V, &h, &p, &q, &t);
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
    uint8_t retval = 0;

    get_hsv2rgb_params(phi, S, V, &h, &p, &q, &t);
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
    uint8_t retval = 0;

    get_hsv2rgb_params(phi, S, V, &h, &p, &q, &t);
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

static uint8_t (*hsv[3])(double phi, double S, double V) = { hsv_g, hsv_b, hsv_r };

static double fhsv_r(double phi, double S, double V) {
    double p,q,t;
    uint8_t h;
    double retval = 0;

    get_hsv2rgb_params(phi, S, V, &h, &p, &q, &t);
    switch(h) {
        case 0:
        case 5: retval = V; break;
        case 1: retval = q; break;
        case 2:
        case 3: retval = p; break;
        case 4: retval = t; break;        
    }
    return retval;
}

static double fhsv_g(double phi, double S, double V) {
    double p,q,t;
    uint8_t h;
    double retval = 0;

    get_hsv2rgb_params(phi, S, V, &h, &p, &q, &t);
    switch(h) {
        case 1:
        case 2: retval = V; break;
        case 3: retval = q; break;
        case 4:
        case 5: retval = p; break;
        case 0: retval = t; break;        
    }
    return retval;
}

static double fhsv_b(double phi, double S, double V) {
    double p,q,t;
    uint8_t h;
    double retval = 0;

    get_hsv2rgb_params(phi, S, V, &h, &p, &q, &t);
    switch(h) {
        case 3:
        case 4: retval = V; break;
        case 5: retval = q; break;
        case 0:
        case 1: retval = p; break;
        case 2: retval = t; break;        
    }
    return retval;
}

static double (*fhsv[3])(double phi, double S, double V) = { fhsv_g, fhsv_b, fhsv_r };

// ************************** functions *************************************/


typedef struct Func {
    const char *name;
    void (*f)(FracDFuncParams *params, AVFrame *in, int x, int y, int n);
} Func;

static void zero(FracDFuncParams *params, AVFrame *in, int x, int y, int n) {
}

static void sq(FracDFuncParams *params, AVFrame *in, int x, int y, int n) {
    int plane,x0,y0;
    uint8_t *ptr;
    x0 = x-params->w/2;
    y0 = y-params->h/2;
    for(plane=0;plane<3;plane++) {
        ptr = in->data[plane] + in->linesize[plane] * y;
        ptr[x] = (sqrt(x0*x0+y0*y0)/(params->w/2))*n*params->p[plane];
    }
}

static void sqs(FracDFuncParams *params, AVFrame *in, int x, int y, int n) {
    int plane,x0,y0;
    uint8_t *ptr;

    x0 = (x-params->w/2)*params->fx + params->x;
    y0 = (y-params->h/2)*params->fy + params->y;
    for(plane=0;plane<3;plane++) {
        ptr = in->data[plane] + in->linesize[plane] * y;
        ptr[x] = asin(sqrt(x0*x0+y0*y0)/(params->w/2))*n*params->p[plane];
    }
}

static void c(FracDFuncParams *params, AVFrame *in, int x, int y, int n) {
    double xi0,yi0;
    MYCOMPLEX z,_z;
    MYFLT x0,y0;
    uint8_t *ptr;
    int plane;

    xi0 = ((double)x-params->w/2)/(double)params->w;
    yi0 = ((double)y-params->h/2)/(double)params->w;
    CA(_z,xi0,yi0);
    _z = av_genutil_rotate(_z,params->rot);
    x0 = (MYREAL(_z) * params->fx + params->x); 
    y0 = (MYIMAG(_z) * params->fy + params->y);
    CA(_z,x0,y0);
    z = params->ifunc(_z,params->ip);
    for(plane=0;plane<3;plane++) {
        double color = params->fac * MYCABS(z) * hsv[plane](cargl(z),params->sat,1);
        ptr = in->data[plane] + in->linesize[plane] * y;
        ptr[x] = floor(color);
    }
}

static void cf(FracDFuncParams *params, AVFrame *in, int x, int y, int n) {
    double xi0,yi0;
    MYCOMPLEX z,_z;
    MYFLT x0,y0;
    uint8_t *ptr;
    int plane;

    xi0 = ((double)x-params->w/2)/(double)params->w;
    yi0 = ((double)y-params->h/2)/(double)params->w;
    CA(_z,xi0,yi0);
    _z = av_genutil_rotate(_z,params->rot);
    x0 = (MYREAL(_z) * params->fx + params->x); 
    y0 = (MYIMAG(_z) * params->fy + params->y);
    CA(_z,x0,y0);
    z = params->ifunc(_z,params->ip);

    for(plane=0;plane<3;plane++) {
        ptr = in->data[plane] + in->linesize[plane] * y;
        ptr[x] = (sin(params->fac * MYCABS(z) * fhsv[plane](cargl(z),params->sat,1))+1)*127;
    }
}

static void cs(FracDFuncParams *params, AVFrame *in, int x, int y, int n) {
    double xi0,yi0;
    MYCOMPLEX z,_z;
    MYFLT x0,y0;
    uint8_t *ptr;
    int plane;

    xi0 = ((double)x-params->w/2)/(double)params->w;
    yi0 = ((double)y-params->h/2)/(double)params->w;
    CA(_z,xi0,yi0);
    _z = av_genutil_rotate(_z,params->rot);
    x0 = (MYREAL(_z) * params->fx + params->x); 
    y0 = (MYIMAG(_z) * params->fy + params->y);
    CA(_z,x0,y0);
    z = params->ifunc(_z,params->ip);

    for(plane=0;plane<3;plane++) {
        ptr = in->data[plane] + in->linesize[plane] * y;
        ptr[x] = (sin(params->fac * MYCABS(z) * hsv[plane](cargl(z),params->sat,1))+1)*127;
    }
}


static void j(FracDFuncParams *params, AVFrame *in, int x, int y, int n) {
    int plane,len,k;
    double colors[3];
    MYCOMPLEX z;
    uint8_t *ptr;

    len = floor(params->p[0]);
    CA(z,(x-params->w/2)*params->fx+params->x, (y-params->h/2)*params->fy+params->y);
    for(k=0;k<len;k++) {
        z = params->ifunc(z,params->ip);
        if(params->ifcmode == 0 && MYCABS(z) > params->p[1]) {
            break;
        }
    }
    av_genutil_get_color(params->cfunc, params->cp,k,params->cmod, params->is_rgb, colors);
    for(plane=0;plane<3;plane++) {
        ptr = in->data[plane] + in->linesize[plane] * y;
        ptr[x] = (k>=len)?0:floor(colors[plane]*255);
    }
}

static void m(FracDFuncParams *params, AVFrame *in, int x, int y, int n) {
    int plane,len,k;
    MYFLT xr0,yr0;
    MYCOMPLEX zr,z; 
    double max;
    double colors[3];
    uint8_t *ptr;
    MYFLT ip[2],modulus,mu;

    len = params->p[0];
    xr0 = x-params->w/2;
    yr0 = y-params->h/2;
    CA(zr,xr0,yr0);
    zr = rotate(zr,params->rot);
    xr0 = MYREAL(zr);
    yr0 = MYIMAG(zr);
    ip[0] = xr0*params->fx+params->x; 
    ip[1] = yr0*params->fy+params->y;
    z = 0;
    max=0;
    for(k=0;k<len;k++) {
        z = params->ifunc(z,ip);
        if(params->ifcmode == 0 && MYCABS(z) > params->p[1]) {
            max = k;
            break;
        }
        if(params->ifcmode == 2) {
            if(MYFABS(MYREAL(z)) > max) max = (MYFABS(MYREAL(z)));
        }
    }
    modulus = MYCABS(z);
    mu = max - (MYLOG (MYLOG (modulus)))/ MYLOG (2.0);
    max = (float)mu;
    av_genutil_get_color(params->cfunc, params->cp,max,params->cmod, params->is_rgb, colors);
    for(plane=0;plane<3;plane++) {
        ptr = in->data[plane] + in->linesize[plane] * y;
        ptr[x] = (k>=len)?0:floor(colors[plane]*255);
    }
}

static void hopa(FracDFuncParams *params, AVFrame *in, int x, int y, int n) {
    int plane,len,k;
    MYCOMPLEX z;
    double max;
    double colors[3];
    uint8_t *ptr;
    MYFLT ip[10];

    len = params->p[0];

    ip[0] = params->ip[0] + params->p[2] * (y+params->y)*params->fy;
    ip[1] = params->ip[1] + params->p[3] * (x+params->x)*params->fx;
    ip[2] = params->ip[2];
    z = 0;
    max = 0;
    for(k=0;k<len;k++) {
        z = params->ifunc(z,ip);
        if(params->ifcmode == 0 && MYCABS(z) > params->p[1]) {
            max = k;
            break;
        }
        if(params->ifcmode == 1 && (MYFABS(MYREAL(z)) > params->p[1])) {
            max = k;
            break;
        }
        if(params->ifcmode == 2) {
            if(MYFABS(MYREAL(z)) > max) max = (MYFABS(MYREAL(z)));
        }
    }
    av_genutil_get_color(params->cfunc, params->cp,max,params->cmod, params->is_rgb, colors);
    for(plane=0;plane<3;plane++) {
        ptr = in->data[plane] + in->linesize[plane] * y;
        ptr[x] = (k>=len&&params->ifcmode!=2)?0:floor(colors[plane]*255);
    }
}

static void hopa1(FracDFuncParams *params, AVFrame *in, int x, int y, int n) {
    int plane,k,len;
    double max;
    MYFLT y00,x0,y0;
    MYCOMPLEX zr,z; 
    double colors[3];
    uint8_t *ptr;
    MYFLT ip[10];

    x0 = x-params->w/2;
    y0 = y-params->h/2;
    CA(zr,x0,y0);
    zr = av_genutil_rotate(zr,params->rot);
    x0 = MYREAL(zr);
    y0 = MYIMAG(zr);
    len = params->p[0];
    y00 = y0*params->fy+params->y;
    ip[0] = params->ip[0] + params->p[3] * (y00==0?params->fy:y00);
    ip[1] = params->ip[1] + params->p[2] * (x0*params->fx+params->x);
    ip[2] = params->ip[2];
    ip[3] = params->ip[3];
    CA(z,params->p[4], params->p[5]);
    for(k=0;k<len;k++) {
        z = params->ifunc(z,ip);
        if(params->ifcmode == 0 && MYCABS(z) > params->p[1]) {
            max = k;
            break;
        }
        if(params->ifcmode == 1 && (MYFABS(MYREAL(z)) > params->p[1])) {
            max = k;
            break;
        }
        if(params->ifcmode == 2) {
            if(fabsl(MYREAL(z)) > max) max = (MYFABS(MYREAL(z)));
        } else if(params->ifcmode ==3) {
            if(MYCABS(z) > max) max = MYCABS(z);
        }
    }
    av_genutil_get_color(params->cfunc, params->cp,max,params->cmod, params->is_rgb, colors);
    for(plane=0;plane<3;plane++) {
        ptr = in->data[plane] + in->linesize[plane] * y;
        ptr[x] = (k>=len&&params->ifcmode!=2)?0:floor(colors[plane]*255);
    }
}

static Func funcs[] = {
    { "zero", zero },
    { "sq", sq },
    { "sqs", sqs },
    { "hopa", hopa },
    { "hopa1", hopa1 },
    { "j", j },
    { "m", m },
    { "c", c },
    { "cs", cs },
    { "cf", cf },
    { NULL, NULL }
};

static void (*get_func(const char *name))(FracDFuncParams*,AVFrame*,int,int,int) {
    int k=0;
    while(funcs[k].name) {
        if(!strcmp(name, funcs[k].name)) {
            return funcs[k].f;
        }
        k++;
    }
    return NULL;
}

static void parse_ffunc(const char *ff, double *p, void (**f)(FracDFuncParams*,AVFrame*,int,int,int)) {
    char *saveptr,*token;
    char *str = strdup(ff);
    int j;
    for (j = 0; ; j++, str = NULL) {
        token = av_strtok(str, " |", &saveptr);
        if (token == NULL)
            break;
        if(j==0) {
            *f = get_func(token);
        } else {
            p[j-1] = atof(token);
        }
    }
    free(str);
}

/* ================================== FILTER ========================================== */
static av_cold void logParameters(AVFilterContext *ctx,double* p, int len) {
    int k;
    for(k=0;k<len;k++) {
        av_log(ctx, AV_LOG_INFO," %f",p[k]); 
    }
    av_log(ctx, AV_LOG_INFO,"\n"); 
}

static av_cold void logParametersL(AVFilterContext *ctx,MYFLT* p, int len) {
    int k;
    for(k=0;k<len;k++) {
        av_log(ctx, AV_LOG_INFO," %Lf",p[k]); 
    }
    av_log(ctx, AV_LOG_INFO,"\n"); 
}

static av_cold int fracd_init(AVFilterContext *ctx)
{
    int k;
    FracDContext *fracd = ctx->priv;
    av_log(ctx, AV_LOG_INFO, "x=%s y=%s fx=%s fy=%s\n",fracd->str_x, fracd->str_y, fracd->str_fx, fracd->str_fy);
    fracd->fx = MYSTRTOF(fracd->str_fx);
    fracd->fy = MYSTRTOF(fracd->str_fy);
    fracd->x = MYSTRTOF(fracd->str_x);
    fracd->y = MYSTRTOF(fracd->str_y);
    av_log(ctx, AV_LOG_INFO, "rgb=%d ofs=%d x=%Lf y=%Lf fx=%Lf fy=%Lf\n",fracd->is_rgb, fracd->offset, fracd->x, fracd->y, fracd->fx, fracd->fy);

    for(k=0;k<10;k++) {
        if(fracd->nf[k]) {
            parse_nfunc(fracd->nf[k],fracd->np[k],&fracd->nfunc[k]); 
            if(!fracd->nfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for nf%d not found %s\n", k, fracd->nf[k]);
            } else {
                av_log(ctx, AV_LOG_INFO, "function for nf%d is [%s]", k, fracd->nf[k]);
                logParametersL(ctx,fracd->np[k],10);
            }
        }
    }
    if(fracd->f) {
        parse_ffunc(fracd->f,fracd->p,&fracd->ffunc); 
        if(!fracd->ffunc) {
            av_log(ctx, AV_LOG_WARNING, "function for ff not found %s\n", fracd->f);
        } else {
            av_log(ctx, AV_LOG_INFO, "function for ff is [%s]", fracd->f);
            logParameters(ctx,fracd->p,40);
        }
    } 
    if(!fracd->ffunc) {
        fracd->ffunc = zero;
    }

    if(fracd->ifc) {
        parse_ifunc(fracd->ifc,fracd->ip,&fracd->ifunc); 
        if(!fracd->ifunc) {
            av_log(ctx, AV_LOG_WARNING, "function for if not found %s\n", fracd->ifc);
        } else {
            av_log(ctx, AV_LOG_INFO, "function for if is [%s]", fracd->ifc);
            logParametersL(ctx,fracd->ip,40);
        }
    } 
    if(!fracd->ifunc) {
        fracd->ifunc = izero;
    }

    for(k=0;k<3;k++) {
        if(fracd->cf[k]) {
            av_genutil_parse_cfunc(fracd->cf[k],fracd->cp[k],&fracd->cfunc[k]); 
            if(!fracd->cfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for cf%d not found %s\n", k, fracd->cf[k]);
            } else {
                av_log(ctx, AV_LOG_INFO, "function for cf%d is [%s]", k, fracd->cf[k]);
                logParameters(ctx,fracd->cp[k],10);
            }
        } 
     }

    return 0;
}

static int fracd_query_formats(AVFilterContext *ctx)
{
    FracDContext *fracd = ctx->priv;
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

    if (fracd->is_rgb) {
        fmts_list = ff_make_format_list(rgb_pix_fmts);
    } else
        fmts_list = ff_make_format_list(yuv_pix_fmts);
    if (!fmts_list)
        return AVERROR(ENOMEM);
    return ff_set_common_formats(ctx, fmts_list);
}

static int fracd_config_props(AVFilterLink *inlink)
{
    FracDContext *fracd = inlink->dst->priv;
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(inlink->format);

    av_assert0(desc);
    fracd->w = inlink->w;
    fracd->h = inlink->h;
    fracd->hsub = desc->log2_chroma_w;
    fracd->vsub = desc->log2_chroma_h;
    fracd->bps = desc->comp[0].depth;
    fracd->planes = desc->nb_components;
    return 0;
}

typedef struct ThreadData {
    int n;
    AVFrame *in;
    FracDFuncParams *params;
} ThreadData;

static int slice_fracd_filter(AVFilterContext *ctx, void *arg, int jobnr, int nb_jobs)
{
    FracDContext *fracd = ctx->priv;
    ThreadData *td = arg;
    const int slice_start = (td->in->height *  jobnr) / nb_jobs;
    const int slice_end = (td->in->height * (jobnr+1)) / nb_jobs;
    int x, y;
    
    for (y = slice_start; y < slice_end; y++) {
        for (x = 0; x < td->in->width; x++) {
            fracd->ffunc(td->params,td->in,x,y, td->n + fracd->offset);
        }
    }
    return 0;
}

static void make_params(FracDContext *fracd, FracDFuncParams *params, int frame_number) {
    int k,j,i;
    const char *format = "%s %d %s %s %d";
    MYFLT np[10][100];

    for(k=0;k<40;k++) params->p[k] = fracd->p[k];
    for(k=0;k<3;k++) params->cfunc[k] = fracd->cfunc[k];
    for(k=0;k<3;k++) {
        for(j=0;j<10;j++) params->cp[k][j] = fracd->cp[k][j];
    }
    params->cmod = fracd->cmod;
    params->ifcmode = fracd->ifcmode;
    for(j=0;j<40;j++) params->ip[j] = fracd->ip[j];
    params->ifunc = fracd->ifunc;
    params->w = fracd->w;
    params->h = fracd->h;
    params->is_rgb = fracd->is_rgb;
    params->fx = fracd->fx;
    params->fy = fracd->fy;
    params->x = fracd->x;
    params->y = fracd->y;
    params->rot = fracd->rot;
    params->sat = fracd->sat;
    params->fac = fracd->fac;

    for(k=0;k<10;k++) {
        for(j=0;j<100;j++) {
            np[k][j] = fracd->np[k][j];
        }
    }

    i=0;
    for(j=0;j<10;j++) {
        if(fracd->rf[j]) {
            char input[40];
            char target[4];
            char src[4];
            char mode[4];
            int index,m;
            MYFLT value;
            strcpy(input, fracd->rf[j]);
            sscanf(input, format, target, &index, mode, src, &m);
            switch(src[0]) {
                case 'n': {
                    value = fracd->nfunc[m](fracd->nmod[m]?frame_number%fracd->nmod[m]:frame_number,np[m]);
                    break;
                }                             
            }
            switch(target[0]) {
                case 'p': {
                    if(mode[0] == 'o') params->p[index] = value;
                    if(mode[0] == 'a') params->p[index] = fracd->p[index] + value;
                    if(mode[0] == 's') params->p[index] = fracd->p[index] - value;
                    break;
                }
                case 'n': {
                    int f = index/10;
                    int i = index % 10;
                    if(mode[0] == 'o') np[f][i] = value;
                    if(mode[0] == 'a') np[f][i] = fracd->np[f][i] + value;
                    if(mode[0] == 's') np[f][i] = fracd->np[f][i] - value;
                    break;
                }
                case 'c': {
                    int f = index/10;
                    int i = index % 10;
                    if(mode[0] == 'o') params->cp[f][i] = value;
                    if(mode[0] == 'a') params->cp[f][i] = fracd->cp[f][i] + value;
                    if(mode[0] == 's') params->cp[f][i] = fracd->cp[f][i] - value;
                    break;
                }
                case 'i': {
                    if(mode[0] == 'o') params->ip[index] = value;
                    if(mode[0] == 'a') params->ip[index] = fracd->ip[index] + value;
                    if(mode[0] == 's') params->ip[index] = fracd->ip[index] - value;
                    break;
                }

                case 'x': {
                    if(mode[0] == 'o') params->x = value;
                    if(mode[0] == 'a') params->x = fracd->x + value;
                    if(mode[0] == 's') params->x = fracd->x - value;
                    break;
                }
                case 'y': {
                    if(mode[0] == 'o') params->y = value;
                    if(mode[0] == 'a') params->y = fracd->y + value;
                    if(mode[0] == 's') params->y = fracd->y - value;
                    break;
                }
                case 'w': {
                    if(mode[0] == 'o') params->fx = value;
                    if(mode[0] == 'a') params->fx = fracd->fx + value;
                    if(mode[0] == 's') params->fx = fracd->fx - value;
                    break;
                }
                case 'h': {
                    if(mode[0] == 'o') params->fy = value;
                    if(mode[0] == 'a') params->fy = fracd->fy + value;
                    if(mode[0] == 's') params->fy = fracd->fy - value;
                    break;
                }
                case 'r': {
                    if(mode[0] == 'o') params->rot = value;
                    if(mode[0] == 'a') params->rot = fracd->rot + value;
                    if(mode[0] == 's') params->rot = fracd->rot - value;
                    break;
                }
                case 'f': {
                    if(mode[0] == 'o') params->fac = value;
                    if(mode[0] == 'a') params->fac = fracd->fac + value;
                    if(mode[0] == 's') params->fac = fracd->fac - value;
                    break;
                }
                case 's': {
                    if(mode[0] == 'o') params->sat = value;
                    if(mode[0] == 'a') params->sat = fracd->sat + value;
                    if(mode[0] == 's') params->sat = fracd->sat - value;
                    break;
                }
            }
        }
        i++;
    }
 
}

static void fdebug(FracDFuncParams *params, int frame_number, AVFrame *in) {
    int k;
    for(k=0;k<4;k++) {
        av_genutil_draw_number_c(k*100+20, params->h-16, params->p[k], 0, in, 0);
        av_genutil_draw_number_c(k*100+20, params->h-16, params->p[k], 0, in, 1);
        av_genutil_draw_number_c(k*100+20, params->h-16, params->p[k], 0, in, 2);
    }
    for(k=4;k<8;k++) {
        av_genutil_draw_number_c(k*100+20, params->h-16, params->ip[k-4], 0, in, 0);
        av_genutil_draw_number_c(k*100+20, params->h-16, params->ip[k-4], 0, in, 1);
        av_genutil_draw_number_c(k*100+20, params->h-16, params->ip[k-4], 0, in, 2);
    }
    av_genutil_draw_number_c(k*100+20, params->h-16, params->fac, 0, in, 0);
    av_genutil_draw_number_c(k*100+20, params->h-16, params->fac, 0, in, 1);
    av_genutil_draw_number_c(k*100+20, params->h-16, params->fac, 0, in, 2);
    k++;
    av_genutil_draw_number_c(k*100+20, params->h-16, params->sat, 0, in, 0);
    av_genutil_draw_number_c(k*100+20, params->h-16, params->sat, 0, in, 1);
    av_genutil_draw_number_c(k*100+20, params->h-16, params->sat, 0, in, 2);
    
    av_genutil_draw_int_c(params->w-50, params->h-16, frame_number, 0, in, 0);
    av_genutil_draw_int_c(params->w-50, params->h-16, frame_number, 0, in, 1);
    av_genutil_draw_int_c(params->w-50, params->h-16, frame_number, 0, in, 2);
    av_plot_form(1,params->w/2,params->h/2,1,params->w,params->h,in->linesize[0],in->data[0]);

    av_genutil_draw_number_c(20, params->h-32, params->fx, 0, in, 0);
    av_genutil_draw_number_c(20, params->h-32, params->fx, 0, in, 1);
    av_genutil_draw_number_c(20, params->h-32, params->fx, 0, in, 2);
    av_genutil_draw_number_c(120, params->h-32, params->fy, 0, in, 0);
    av_genutil_draw_number_c(120, params->h-32, params->fy, 0, in, 1);
    av_genutil_draw_number_c(120, params->h-32, params->fy, 0, in, 2);
    av_genutil_draw_int_c(220, params->h-32, params->x, 0, in, 0);
    av_genutil_draw_int_c(220, params->h-32, params->x, 0, in, 1);
    av_genutil_draw_int_c(220, params->h-32, params->x, 0, in, 2);
    av_genutil_draw_int_c(320, params->h-32, params->y, 0, in, 0);
    av_genutil_draw_int_c(320, params->h-32, params->y, 0, in, 1);
    av_genutil_draw_int_c(320, params->h-32, params->y, 0, in, 2);
    av_genutil_draw_number_c(420, params->h-32, params->rot, 0, in, 0);
    av_genutil_draw_number_c(420, params->h-32, params->rot, 0, in, 1);
    av_genutil_draw_number_c(420, params->h-32, params->rot, 0, in, 2);
}


static int fracd_filter_frame(AVFilterLink *inlink, AVFrame *in)
{
    int nb_threads;
    double *p;
    double n;
    FracDFuncParams params;
    ThreadData td;
    AVFilterContext *ctx = inlink->dst;
    FracDContext *fracd = ctx->priv;
    AVFilterLink *outlink = inlink->dst->outputs[0];

    srand(0);
    nb_threads = ff_filter_get_nb_threads(ctx);
    av_frame_make_writable(in);
    n = inlink->frame_count_out;
    p = calloc(40,sizeof(double));
    params.p = p;
    make_params(fracd,&params,n+fracd->offset);
    
    td.in = in;
    td.n = n;
    td.params = &params;
    ctx->internal->execute(ctx, slice_fracd_filter, &td, NULL, FFMIN(in->height, nb_threads));
    if(fracd->dbg) fdebug(&params,n,in);
    free(p);
    return ff_filter_frame(outlink, in);
}

static av_cold void fracd_uninit(AVFilterContext *ctx)
{
    av_log(ctx, AV_LOG_INFO, "uninit\n");
}

static const AVFilterPad fracd_inputs[] = {
    {
        .name         = "default",
        .type         = AVMEDIA_TYPE_VIDEO,
        .config_props = fracd_config_props,
        .filter_frame = fracd_filter_frame,
    },
    { NULL }
};

static const AVFilterPad fracd_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO,
    },
    { NULL }
};

AVFilter ff_vf_fracd = {
    .name          = "fracd",
    .description   = NULL_IF_CONFIG_SMALL("Apply generic equation to each pixel."),
    .priv_size     = sizeof(FracDContext),
    .init          = fracd_init,
    .uninit        = fracd_uninit,
    .query_formats = fracd_query_formats,
    .inputs        = fracd_inputs,
    .outputs       = fracd_outputs,
    .priv_class    = &fracd_class,
    .flags         = AVFILTER_FLAG_SUPPORT_TIMELINE_GENERIC | AVFILTER_FLAG_SLICE_THREADS,
};
