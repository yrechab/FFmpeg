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


typedef struct CGen2Context {
    const AVClass *class;
    const char *f;
    double p[40];
    const char *n;
    double np[20];
    const char *w;
    double wp[20];
    uint8_t (*ffunc)(int,int,int,int,int,int,struct CGen2Context*);
    double (*nfunc)(int,double*);
    double (*wfunc)(int,double*);
    double x;
    double xn;
    double y;
    double yn;
    double sat;
    double fac;
    double div;
    double rot;
    double rotn;
    const char *mf;
    double complex (*cfunc)(double complex, double*, double);
    int hsub, vsub;             ///< chroma subsampling
    int planes;                 ///< number of planes
    int is_rgb;
    int bps;
    int offset;
} CGen2Context;


#define OFFSET(x) offsetof(CGen2Context, x)
#define FLAGS AV_OPT_FLAG_VIDEO_PARAM|AV_OPT_FLAG_FILTERING_PARAM

static const AVOption cgen2_options[] = {
    { "f",  "f",   OFFSET(f),  AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "p00","p00",   OFFSET(p[0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p01","p01",   OFFSET(p[1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p02","p02",   OFFSET(p[2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p03","p03",   OFFSET(p[3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p04","p04",   OFFSET(p[4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p05","p05",   OFFSET(p[5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p06","p06",   OFFSET(p[6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p07","p07",   OFFSET(p[7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p08","p08",   OFFSET(p[8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p09","p09",   OFFSET(p[9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p10","p10",   OFFSET(p[10]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p11","p11",   OFFSET(p[11]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p12","p12",   OFFSET(p[12]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p13","p13",   OFFSET(p[13]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p14","p14",   OFFSET(p[14]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p15","p15",   OFFSET(p[15]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p16","p16",   OFFSET(p[16]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p17","p17",   OFFSET(p[17]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p18","p18",   OFFSET(p[18]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p19","p19",   OFFSET(p[19]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p20","p20",   OFFSET(p[20]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p21","p21",   OFFSET(p[21]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p22","p22",   OFFSET(p[22]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p23","p23",   OFFSET(p[23]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p24","p24",   OFFSET(p[24]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p25","p25",   OFFSET(p[25]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p26","p26",   OFFSET(p[26]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p27","p27",   OFFSET(p[27]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p28","p28",   OFFSET(p[28]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p29","p29",   OFFSET(p[29]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p30","p30",   OFFSET(p[30]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p31","p31",   OFFSET(p[31]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p32","p32",   OFFSET(p[32]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p33","p33",   OFFSET(p[33]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p34","p34",   OFFSET(p[34]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p35","p35",   OFFSET(p[35]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p36","p36",   OFFSET(p[36]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p37","p37",   OFFSET(p[37]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p38","p38",   OFFSET(p[38]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p39","p39",   OFFSET(p[39]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "x","x",       OFFSET(x),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "y","y",       OFFSET(y),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "xn","xn",     OFFSET(xn),    AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "yn","yn",     OFFSET(yn),    AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "sat","sat",   OFFSET(sat),   AV_OPT_TYPE_DOUBLE, {.dbl =  1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "fac","fac",   OFFSET(fac),   AV_OPT_TYPE_DOUBLE, {.dbl =  1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "div","div",   OFFSET(div),   AV_OPT_TYPE_DOUBLE, {.dbl =  1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "rot","rot",   OFFSET(rot),   AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "rotn","rotn", OFFSET(rotn),  AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "mf",  "mf",   OFFSET(mf),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "n",  "n",     OFFSET(n),     AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "n00","n00",   OFFSET(np[0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n01","n01",   OFFSET(np[1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n02","n02",   OFFSET(np[2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n03","n03",   OFFSET(np[3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n04","n04",   OFFSET(np[4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n05","n05",   OFFSET(np[5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n06","n06",   OFFSET(np[6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n07","n07",   OFFSET(np[7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n08","n08",   OFFSET(np[8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n09","n09",   OFFSET(np[9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n10","n10",   OFFSET(np[10]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n11","n11",   OFFSET(np[11]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n12","n12",   OFFSET(np[12]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n13","n13",   OFFSET(np[13]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n14","n14",   OFFSET(np[14]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n15","n15",   OFFSET(np[15]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n16","n16",   OFFSET(np[16]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n17","n17",   OFFSET(np[17]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n18","n18",   OFFSET(np[18]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n19","n19",   OFFSET(np[19]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w",  "w",     OFFSET(w),      AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "w00","w00",   OFFSET(wp[0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w01","w01",   OFFSET(wp[1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w02","w02",   OFFSET(wp[2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w03","w03",   OFFSET(wp[3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w04","w04",   OFFSET(wp[4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w05","w05",   OFFSET(wp[5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w06","w06",   OFFSET(wp[6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w07","w07",   OFFSET(wp[7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w08","w08",   OFFSET(wp[8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w09","w09",   OFFSET(wp[9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w10","w10",   OFFSET(wp[10]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w11","w11",   OFFSET(wp[11]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w12","w12",   OFFSET(wp[12]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w13","w13",   OFFSET(wp[13]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w14","w14",   OFFSET(wp[14]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w15","w15",   OFFSET(wp[15]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w16","w16",   OFFSET(wp[16]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w17","w17",   OFFSET(wp[17]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w18","w18",   OFFSET(wp[18]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "w19","w19",   OFFSET(wp[19]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "rgb","rgb",OFFSET(is_rgb), AV_OPT_TYPE_INT,    {.i64=0},    0,        1,        FLAGS },
    { "ofs","ofs",OFFSET(offset), AV_OPT_TYPE_INT,    {.i64=0},    0,     INT_MAX,        FLAGS },
    {NULL},
};

AVFILTER_DEFINE_CLASS(cgen2);


typedef struct Func {
    const char *name;
    uint8_t (*f)(int plane, int x, int y, int w, int h, int n, CGen2Context *ctx);
} Func;

typedef struct NFunc {
    const char *name;
    double (*f)(int n, double *p);
} NFunc;

typedef struct CFunc {
    const char *name;
    double complex (*f)(double complex z, double *p, double t);
} CFunc;

#define CFUNC(x) double complex (*x)(double complex z, double *p, double t)

static double complex rotate(double complex z, double angle, double t) {
    if(angle == 0.0 && t == 0.0) return z;
    double phi = angle*(M_PI/180)+t;
    return (cos(phi) + sin(phi)*I) * z;
}

static double complex rat(double complex z, double *p, double t) {
    int zg = p[0];
    double complex zl = p[1];
    int k;
    for(k = 0; k<(zg*2); k+=2) {
        zl *= (z-(p[k+2] + p[k+3]*I));
    }
    int ng = p[20];
    if(ng == 0) return zl;

    double complex nen = p[21];
    for(k = 20; k<(ng*2+20); k+=2) {
        nen *= (z-(p[k+2] + p[k+3]*I));
    }
    return zl/nen;
}


static double complex id(double complex z, double *p, double t) {
    return z;
}

static double complex fexp(double complex z, double *p, double t) {
    return p[0] + (p[1]?p[1]:1) * z * cexp( I * (p[2]?p[2]:1)/(1/pow(cabs(z),p[3])));
}

static double complex expr(double complex z, double *p, double t) {
    return cexp(rat(z,p,t));
}

static double complex sinr(double complex z,  double *p, double t) {
    return csin(rat(z,p,t));
}

static double complex tanr(double complex z, double *p, double t) {
    return ctan(rat(z,p,t));
}


static CFunc cfuncs[] = {
    {"rat",rat},
    {"exp",fexp},
    {"expr",expr},
    {"sinr",sinr},
    {"tanr",tanr},
    {NULL,NULL}
};

static double complex (*getCFunc(const char *name))(double complex, double*, double) {
    int k=0;
    while(cfuncs[k].name) {
        if(!strcmp(name, cfuncs[k].name)) {
            return cfuncs[k].f;
        }
        k++;
    }
    return NULL;
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

static double npoly(int n, double *p) {
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
    if((section+1) >= 20) return 0;
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
    {"poly",npoly},
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

static uint8_t (*hsv[3])(double phi, double S, double V) = { hsv_g, hsv_b, hsv_r };


static double complex **last;

/* ============================= FUNCS =========================================== */
static uint8_t ffunc(double complex (*func)(double complex z, double *p, double t), int plane, int x, int y, int w, int h, int n, CGen2Context *ctx) {
    double f = ctx->wfunc(n-ctx->offset,ctx->wp)/(double)w;
    double complex _z = rotate((x*f+ctx->x+ctx->xn*n) + ((y*f+ctx->y+ctx->yn*n) * I), ctx->rot, ctx->rotn * n);
    //double complex z = func((x*f+ctx->x+ctx->xn*n) + ((y*f+ctx->y+ctx->yn*n) * I), ctx->p,(n-ctx->offset)*ctx->div);
    double complex z = func(_z, ctx->p,(n-ctx->offset)*ctx->div);
    return ctx->nfunc(n,ctx->np) * ctx->fac * cabs(z) * hsv[plane](carg(z), ctx->sat, 1);
}

static uint8_t ffunc2(double complex (*func)(double complex z, double *p, double t), int plane, int x, int y, int w, int h, int n, CGen2Context *ctx) {
    double f = ctx->wfunc(n-ctx->offset,ctx->wp)/(double)w;
    double complex _z = rotate((x*f+ctx->x+ctx->xn*n) + ((y*f+ctx->y+ctx->yn*n) * I), ctx->rot, ctx->rotn * n);
    //double complex z = func((x*f+ctx->x+ctx->xn*n) + ((y*f+ctx->y+ctx->yn*n) * I), ctx->p,(n-ctx->offset)*ctx->div);
    double complex z = func(_z, ctx->p,(n-ctx->offset)*ctx->div);
    return (sin(ctx->nfunc(n,ctx->np) * ctx->fac * cabs(z) * hsv[plane](carg(z), ctx->sat, 1))+1)*127;
}

static uint8_t fexp1(int plane,int x, int y, int w, int h, int n, CGen2Context *ctx) {
    return ffunc(fexp,plane,x,y,w,h,n,ctx);
}

static uint8_t frat(int plane,int x, int y, int w, int h, int n, CGen2Context *ctx) {
    return ffunc(rat,plane,x,y,w,h,n,ctx);
}

static uint8_t frats(int plane,int x, int y, int w, int h, int n, CGen2Context *ctx) {
    return ffunc2(rat,plane,x,y,w,h,n,ctx);
}

static uint8_t fid(int plane,int x, int y, int w, int h, int n, CGen2Context *ctx) {
    return ffunc(id,plane,x,y,w,h,n,ctx);
}

static uint8_t sqasin(int plane,int x, int y, int w, int h, int n, CGen2Context *ctx) {
    double f = ctx->wfunc(n-ctx->offset,ctx->wp)/(double)h;
    double x1=f*x;
    double y1=f*y;
    return asin(sqrt(x1*x1+y1*y1)/(h/2))*n*ctx->p[plane];
}

static uint8_t fsin(int plane,int x, int y, int w, int h, int n, CGen2Context *ctx) {
    return ffunc(sinr,plane,x,y,w,h,n,ctx);
}

static uint8_t ftan(int plane,int x, int y, int w, int h, int n, CGen2Context *ctx) {
    return ffunc(tanr,plane,x,y,w,h,n,ctx);
}

static uint8_t fexpr(int plane,int x, int y, int w, int h, int n, CGen2Context *ctx) {
    return ffunc(expr,plane,x,y,w,h,n,ctx);
}


static uint8_t map(int plane,int x, int y, int w, int h, int n, CGen2Context *ctx) {
    int nn = n - ctx->offset;
    double f = ctx->p[4];
    if(nn == 0) {
        double complex z = x*f + y*f * I;
        last[x+960][y+540] = z;
        return hsv[plane](carg(z),ctx->p[plane], cabs(z)/(h/2));
    } else {
        double complex z = ctx->cfunc(last[x+960][y+540], ctx->p, n*ctx->div);
        last[x+960][y+540] = z;
        return hsv[plane](carg(z),ctx->p[plane], cabs(z)/(h/2));
    }
}

static uint8_t zero(int plane,int x, int y, int w, int h, int n, CGen2Context *ctx) {
    return 0;
}

static Func funcs[] = {
    { "zero", zero },
    { "sqasin", sqasin },
    { "frat", frat },
    { "frats", frats },
    { "fexp", fexp1 },
    { "fexpr", fexpr },
    { "fsin", fsin },
    { "ftan", ftan },
    { "fid", fid },
    { NULL, NULL }
};

static uint8_t (*getFunc(const char *name))(int,int,int,int,int,int,CGen2Context*) {
    int k=0;
    while(funcs[k].name) {
        if(!strcmp(name, funcs[k].name)) {
            return funcs[k].f;
        }
        k++;
    }
    return NULL;
}

/* ================================== FILTER ========================================== */
static av_cold void logParameters(AVFilterContext *ctx,double* p, int len) {
    int k;
    for(k=0;k<len;k++) {
        av_log(ctx, AV_LOG_INFO," %f",p[k]); 
    }
    av_log(ctx, AV_LOG_INFO,"\n"); 
}

static av_cold int cgen2_init(AVFilterContext *ctx)
{
    int k;
    CGen2Context *cgen2 = ctx->priv;
    av_log(ctx, AV_LOG_INFO, "rgb=%d ofs=%d x=%f y=%f xn=%f yn=%f, sat=%f fac=%f div=%f\n", cgen2->is_rgb, cgen2->offset, cgen2->x, cgen2->y, cgen2->xn, cgen2->yn, cgen2->sat, cgen2->fac, cgen2->div);
    if(cgen2->w) {
        cgen2->wfunc = getNFunc(cgen2->w);
        if(!cgen2->wfunc) {
            av_log(ctx, AV_LOG_WARNING, "function for wf not found %s\n", cgen2->w);
            cgen2->wfunc = constant;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for wf is %s", cgen2->w);
            logParameters(ctx,cgen2->wp,20);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for w\n");
        cgen2->wfunc = constant;
    }

    if(cgen2->n) {
        cgen2->nfunc = getNFunc(cgen2->n);
        if(!cgen2->nfunc) {
            av_log(ctx, AV_LOG_WARNING, "function for nf not found %s\n", cgen2->n);
            cgen2->nfunc = idn;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for nf is %s", cgen2->n);
            logParameters(ctx,cgen2->np,20);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for n\n");
        cgen2->nfunc = idn;
    }

    if(cgen2->f) {
        cgen2->ffunc = getFunc(cgen2->f);
        if(!cgen2->ffunc) {
            av_log(ctx, AV_LOG_WARNING, "function for f not found %s\n", cgen2->f);
            cgen2->ffunc = zero;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for f is %s\n", cgen2->f);
            logParameters(ctx,cgen2->p,40);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for f\n");
        cgen2->ffunc = zero;
    }

    if(cgen2->mf) {
        CFUNC(cf) = getCFunc(cgen2->mf);
        if(cf) {
            cgen2->cfunc = cf;
            cgen2->ffunc = map;
        } else {
            av_log(ctx, AV_LOG_WARNING, "function for mf not found %s\n", cgen2->mf);
            cgen2->ffunc = zero;
        }
    }

    last = (double complex **)malloc(1920 * sizeof(double complex *));
    for(k=0;k<1920;k++) {
        last[k] = (double complex *)malloc(1080 * sizeof(double complex));
    }

    return 0;
}

static int cgen2_query_formats(AVFilterContext *ctx)
{
    CGen2Context *cgen2 = ctx->priv;
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

    if (cgen2->is_rgb) {
        fmts_list = ff_make_format_list(rgb_pix_fmts);
    } else
        fmts_list = ff_make_format_list(yuv_pix_fmts);
    if (!fmts_list)
        return AVERROR(ENOMEM);
    return ff_set_common_formats(ctx, fmts_list);
}

static int cgen2_config_props(AVFilterLink *inlink)
{
    CGen2Context *cgen2 = inlink->dst->priv;
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(inlink->format);

    av_assert0(desc);

    cgen2->hsub = desc->log2_chroma_w;
    cgen2->vsub = desc->log2_chroma_h;
    cgen2->bps = desc->comp[0].depth;
    cgen2->planes = desc->nb_components;
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

static int slice_cgen2_filter(AVFilterContext *ctx, void *arg, int jobnr, int nb_jobs)
{
    CGen2Context *cgen2 = ctx->priv;
    ThreadData *td = arg;
    const int height = td->height;
    const int width = td->width;
    const int plane = td->plane;
    const int linesize = td->linesize;
    const int slice_start = (height *  jobnr) / nb_jobs;
    const int slice_end = (height * (jobnr+1)) / nb_jobs;
    int x, y;
    uint8_t *ptr;
    
    for (y = slice_start; y < slice_end; y++) {
        ptr = td->in->data[plane] + linesize * y;

        for (x = 0; x < width; x++) {
            ptr[x] = cgen2->ffunc(plane, x-width/2,y-height/2, width, height, td->n + cgen2->offset, cgen2);
        }
    }
    return 0;
}

static int cgen2_filter_frame(AVFilterLink *inlink, AVFrame *in)
{
    int plane;
    AVFilterContext *ctx = inlink->dst;
    const int nb_threads = ff_filter_get_nb_threads(ctx);
    CGen2Context *cgen2 = ctx->priv;
    AVFilterLink *outlink = inlink->dst->outputs[0];
    av_frame_make_writable(in);
    double n = inlink->frame_count_out;
    //double t = in->pts == AV_NOPTS_VALUE ? NAN : in->pts * av_q2d(inlink->time_base);
    

    for (plane = 0; plane < cgen2->planes && in->data[plane]; plane++) {
        const int width = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(inlink->w, cgen2->hsub) : inlink->w;
        const int height = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(inlink->h, cgen2->vsub) : inlink->h;
        const int linesize = in->linesize[plane];
        ThreadData td;
        td.width = width;
        td.height = height;
        td.plane = plane;
        td.linesize = linesize;
        td.in = in;
        td.n = n;
        ctx->internal->execute(ctx, slice_cgen2_filter, &td, NULL, FFMIN(height, nb_threads));
    }

    return ff_filter_frame(outlink, in);
}

static av_cold void cgen2_uninit(AVFilterContext *ctx)
{
    
    av_log(ctx, AV_LOG_INFO, "uninit\n");
    int k;
    for(k=0;k<1920;k++) {
        free(last[k]);
    }
    free(last);
}

static const AVFilterPad cgen2_inputs[] = {
    {
        .name         = "default",
        .type         = AVMEDIA_TYPE_VIDEO,
        .config_props = cgen2_config_props,
        .filter_frame = cgen2_filter_frame,
    },
    { NULL }
};

static const AVFilterPad cgen2_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO,
    },
    { NULL }
};

AVFilter ff_vf_cgen2 = {
    .name          = "cgen2",
    .description   = NULL_IF_CONFIG_SMALL("Apply generic equation to each pixel."),
    .priv_size     = sizeof(CGen2Context),
    .init          = cgen2_init,
    .uninit        = cgen2_uninit,
    .query_formats = cgen2_query_formats,
    .inputs        = cgen2_inputs,
    .outputs       = cgen2_outputs,
    .priv_class    = &cgen2_class,
    .flags         = AVFILTER_FLAG_SUPPORT_TIMELINE_GENERIC | AVFILTER_FLAG_SLICE_THREADS,
};
