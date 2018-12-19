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
 * Video generation based on a function plotter
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
#include <math.h>


typedef struct PlotContext {
    const AVClass *class;
    const char *f;
    double p[40];
    const char *n;
    double np[20];
    const char *w;
    double wp[20];
    void (*ffunc)(struct PlotContext*,int, int, int, int, uint8_t*);
    double (*nfunc)(int,double*);
    double (*wfunc)(int,double*);
    double x;
    double xn;
    double y;
    double yn;
    double sat;
    double speed;
    double length;
    int form;
    int count;
    int layers;
    int dim;
    double rot;
    double rotn;
    uint8_t *last;
    const char *mf;
    double complex (*cfunc)(double complex, double*, double);
    int hsub, vsub;             ///< chroma subsampling
    int planes;                 ///< number of planes
    int is_rgb;
    int bps;
    int offset;
} PlotContext;


#define OFFSET(x) offsetof(PlotContext, x)
#define FLAGS AV_OPT_FLAG_VIDEO_PARAM|AV_OPT_FLAG_FILTERING_PARAM

static const AVOption plot_options[] = {
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
    { "p38","p38",   OFFSET(p[38]), AV_OPT_TYPE_DOUBLE, {.dbl =  1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p39","p39",   OFFSET(p[39]), AV_OPT_TYPE_DOUBLE, {.dbl =  1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "x","x",       OFFSET(x),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "y","y",       OFFSET(y),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "xn","xn",     OFFSET(xn),    AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "yn","yn",     OFFSET(yn),    AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "sat","sat",   OFFSET(sat),   AV_OPT_TYPE_DOUBLE, {.dbl =  1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "form","form", OFFSET(form),  AV_OPT_TYPE_INT,    {.i64 =  1}, 0,  20,  FLAGS },
    { "layers","l",  OFFSET(layers),AV_OPT_TYPE_INT,    {.i64 =  4}, 1,  INT_MAX,  FLAGS },
    { "count","c",   OFFSET(count), AV_OPT_TYPE_INT,    {.i64 =  1}, 1,  INT_MAX,  FLAGS },
    { "dim","d",     OFFSET(dim),   AV_OPT_TYPE_INT,    {.i64 =  4}, 1,  255,  FLAGS },
    { "speed","s",   OFFSET(speed), AV_OPT_TYPE_DOUBLE, {.dbl =  0.01}, 0,  DBL_MAX,  FLAGS },
    { "length","len",OFFSET(length), AV_OPT_TYPE_DOUBLE, {.dbl =  2*M_PI}, 0,  DBL_MAX,  FLAGS },
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

AVFILTER_DEFINE_CLASS(plot);

typedef struct FFunc {
    const char *name;
    void (*f)(PlotContext*, int, int, int, int, uint8_t*);
} FFunc;

typedef struct NFunc {
    const char *name;
    double (*f)(int n, double *p);
} NFunc;

#define PFUNC(x) double complex (*x)(double t, double *p)

static double complex rotate(double complex z, double phi) {
    return (cos(phi) + sin(phi)*I) * z;
}


static double complex psin(double t,  double *p) {
    return p[0] * cos(t) + I * p[1] * sin(t);
}

static double complex kardio(double t,  double *p) {
    return p[0] * (cos(t)*cos(t)) + p[1] * cos(t) + I * (p[0] * cos(t) * sin(t) + p[1] * sin(t));
}

static double complex lemniskateG(double t,  double *p) {
    double x = cos(t);
    double y = cos(t)*sin(t);
    return p[0] * x + I * p[1] * y;
}

static double complex lemniskateB(double t,  double *p) {
    double x = sqrt(cos(2*t))*cos(t);
    double y = sqrt(cos(2*t))*sin(t);
    return p[0] * x + I * p[1] * y;
}

static double complex epizycloid(double t, double *p) {
    double a = p[0];
    double b = p[1];
    double x = (a+b)*cos(t) - a * cos((1+b/a)*t);
    double y = (a+b)*sin(t) - a * sin((1+b/a)*t);
    return x + I * y;
}

static double complex lissajous(double t, double *p) {
    if(p[3]==0 && p[4] == 0) p[4] = 1.0;
    double a = p[1];
    double b = p[2];
    double c = p[3];
    double d = p[4];
    double x = (c * t + d) * sin(t);
    double y = sin(a * t + b);
    return p[0] * x + I * p[0] * y;
}

static double complex lissajousG(double t, double *p) {
    double a = p[1];
    double b = p[2];
    double c = p[3];
    double d = p[4];
    double e = p[5];
    double x = sin(t);
    double y = sin(a*t+b) + c * sin(d*t+e);
    return p[0] * x + I * p[0] * y;
}

static double complex lissajousQ(double t, double *p) {
    double a = p[1];
    double b = p[2];
    double c = p[3];
    double d = p[4];
    double x = cos(t) + a * sin(b*t) * sin(b*t);
    double y = sin(t) + c * sin(d*t) * sin(d*t);
    return p[0] * x + I * p[0] * y;
}


static double poly3(double a,double b,double x) {
    return a*x*x*x + b*x;
}

static double complex legendre(double t, double *p) {
    double x = p[0] * poly3(p[1],p[2],cos(t));
    double y = p[0] * poly3(p[1],p[2],sin(t));
    return x + I * y;
}

static double complex hypocycloid(double t, double *p) {
    double a = p[1];
    double b = p[2];
    double x = p[0]*((1-a)*cos(a*t)+a*b*cos((1-a)*t));
    double y = p[0]*((1-a)*sin(a*t)-a*b*sin((1-a)*t));
    return x + I * y;
}

static double complex cb(double t, double *p) {
    double a = p[2];
    double b = p[3];
    double x = p[0]*sin(a*t)*cos(t);
    double y = p[1]*sin(b*t)*sin(t);
    return x + I * y;
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

static double sqn(int n, double *p) {
    double x = (double)n/p[0];
    return p[1] * (sqrt(x+0.25)-0.5);
}

static double ln(int n, double *p) {
    double x = (double)n/p[0];
    return p[1] * log(x+1);
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
    {"sq",sqn},
    {"ln",ln},
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


/* ================================= PLOTTER ================================================ */
static void plotPix(int x, int y, double a, int w, int h,int linesize, uint8_t *buf) {
    if(x>=0 && x<w && y>=0 && y<h) {
        buf[x + linesize * y] = floor(a*255);
    }
}

static uint8_t numbers[12][8] = {
   {
       0b00000000,
       0b01111110,
       0b01100110,
       0b01100110,
       0b01100110,
       0b01100110,
       0b01100110,
       0b01111110,
   },
   {
       0b00000000,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
   },
   {
       0b00000000,
       0b01111110,
       0b00000110,
       0b00000110,
       0b01111110,
       0b01100000,
       0b01100000,
       0b01111110,
   },
   {
       0b00000000,
       0b01111110,
       0b00000110,
       0b00000110,
       0b01111110,
       0b00000110,
       0b00000110,
       0b01111110,
   },
   {
       0b00000000,
       0b01100110,
       0b01100110,
       0b01100110,
       0b01111110,
       0b00000110,
       0b00000110,
       0b00000110,
   },
   {
       0b00000000,
       0b01111110,
       0b01100000,
       0b01100000,
       0b01111110,
       0b00000110,
       0b00000110,
       0b01111110,
   },
   {
       0b00000000,
       0b01111110,
       0b01100000,
       0b01100000,
       0b01111110,
       0b01100110,
       0b01100110,
       0b01111110,
   },
   {
       0b00000000,
       0b01111110,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
   },
   {
       0b00000000,
       0b01111110,
       0b01100110,
       0b01100110,
       0b01111110,
       0b01100110,
       0b01100110,
       0b01111110,
   },
   {
       0b00000000,
       0b01111110,
       0b01100110,
       0b01100110,
       0b01111110,
       0b00000110,
       0b00000110,
       0b01111110,
   },
   {
       0b00000000,
       0b00000000,
       0b00000000,
       0b00000000,
       0b00000000,
       0b00000000,
       0b00111000,
       0b00111000,
   },
   {
       0b00000000,
       0b00000000,
       0b00000000,
       0b00000000,
       0b01111110,
       0b00000000,
       0b00000000,
       0b00000000,
   },

};

static void plotCh(char ch, int x, int y, int w, int h, int linesize, uint8_t *buf) {
    uint8_t *bmp;
    if(ch == '.') {
        bmp = numbers[10];
    } else if(ch == '-') {
        bmp = numbers[11];
    } else if(ch>=48 && ch<58) {
        bmp = numbers[ch-48];
    } else {
        return;
    }
    int k,j;
    for(k=0;k<8;k++) {
        uint8_t b = bmp[k];
        for(j=0;j<8;j++) {
            if( (1 << 7-j)&b) {
                plotPix(x+j,y+k,1,w,h,linesize,buf);
            }
        }
    }
}

static void drawNumber(int x, int y, double z, int w, int h, int linesize, uint8_t *buf) {
    char s[50];
    sprintf(s,"%lf",z);
    char *ptr = s;
    int k = 0;
    while(*ptr) {
        plotCh(*ptr,x+k,y,w,h,linesize,buf);
        ptr++;
        k+=9;
    }
}




static void plotForm(int form, int x, int y, double a, int w, int h,int linesize, uint8_t *buf) {
    switch(form) {
        case 0: plotPix(x,y,a,w,h,linesize,buf);
                break;
        case 1: plotPix(x,y,a,w,h,linesize,buf);
                plotPix(x+1,y,a,w,h,linesize,buf);
                plotPix(x,y+1,a,w,h,linesize,buf);
                plotPix(x+1,y+1,a,w,h,linesize,buf);
                break;
        case 2: plotPix(x,y,a,w,h,linesize,buf);
                plotPix(x+1,y,a,w,h,linesize,buf);
                break;
        case 3: plotPix(x,y,a,w,h,linesize,buf);
                plotPix(x,y+1,a,w,h,linesize,buf);
                break;
        case 4: plotPix(x,y,a,w,h,linesize,buf);
                plotPix(x-1,y,a,w,h,linesize,buf);
                plotPix(x+1,y,a,w,h,linesize,buf);
                plotPix(x,y-1,a,w,h,linesize,buf);
                plotPix(x,y+1,a,w,h,linesize,buf);
                break;
        case 5: plotPix(x,y,a,w,h,linesize,buf);
                plotPix(x-1,y,a,w,h,linesize,buf);
                plotPix(x+1,y,a,w,h,linesize,buf);
                plotPix(x,y-1,a,w,h,linesize,buf);
                plotPix(x,y+1,a,w,h,linesize,buf);
                plotPix(x+1,y+1,a,w,h,linesize,buf);
                plotPix(x-1,y-1,a,w,h,linesize,buf);
                break;
        case 6: plotPix(x,y,a,w,h,linesize,buf);
                plotPix(x+1,y,a,w,h,linesize,buf);
                plotPix(x-1,y,a,w,h,linesize,buf);
                plotPix(x+2,y,a,w,h,linesize,buf);
                plotPix(x-2,y,a,w,h,linesize,buf);
                break;
        case 7: plotPix(x,y,a,w,h,linesize,buf);
                plotPix(x+1,y,a,w,h,linesize,buf);
                plotPix(x-1,y,a,w,h,linesize,buf);
                plotPix(x+2,y,a,w,h,linesize,buf);
                plotPix(x-2,y,a,w,h,linesize,buf);
                plotPix(x,y-1,a,w,h,linesize,buf);
                plotPix(x+1,y-1,a,w,h,linesize,buf);
                plotPix(x-1,y-1,a,w,h,linesize,buf);
                plotPix(x,y-2,a,w,h,linesize,buf);
                plotPix(x,y+1,a,w,h,linesize,buf);
                plotPix(x+1,y+1,a,w,h,linesize,buf);
                plotPix(x-1,y+1,a,w,h,linesize,buf);
                plotPix(x,y+2,a,w,h,linesize,buf);
                break;
        case 8: plotPix(x,y,a,w,h,linesize,buf);
                plotPix(x+4,y,a,w,h,linesize,buf);
                plotPix(x-4,y,a,w,h,linesize,buf);
                plotPix(x+8,y,a,w,h,linesize,buf);
                plotPix(x-8,y,a,w,h,linesize,buf);
                plotPix(x+12,y,a,w,h,linesize,buf);
                plotPix(x-12,y,a,w,h,linesize,buf);
                break;
        default:plotPix(x,y,a,w,h,linesize,buf);
                break;
    }
}

static void drive(double t, PFUNC(f), double *p, int w, int h, int linesize, int ofsX, int ofsY, double a, double phi, int form, uint8_t *buf) {
    double complex _z = f(t,p);
    double complex z = rotate(_z,phi);
    int x = floor(p[38]*creal(z) + w/2) + ofsX;
    int y = floor(p[39]*cimag(z) + h/2) + ofsY;
    plotForm(form,x,y,a,w,h,linesize,buf);
} 

static void circs(PlotContext *plot, int n, int w, int h, int linesize, uint8_t *buf) {
    double count = plot->count;
    double speed = plot->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int ofsX = floor(cos(alpha) * plot->p[0]);
        int ofsY = floor(sin(alpha) * plot->p[1]);
        double start = 2*M_PI*k/count;
        int j;
        int layers = plot->layers;
        for(j=0;j<layers;j++) {
            double s = start + plot->length*j/layers;
            double t = s + n*speed;
            drive(t,psin,plot->p,w,h,linesize,ofsX,ofsY,1,0,plot->form,buf);
        }
    }
} 

static void kardioids(PlotContext *plot, int n, int w, int h, int linesize, uint8_t *buf) {
    double count = plot->count;
    double speed = plot->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int ofsX = 0;
        int ofsY = 0;
        int j;
        double start = 0;
        int layers = plot->layers;
        for(j=0;j<layers;j++) {
            double s = start + plot->length*j/layers;
            double t = s + n*speed;
            drive(t,kardio,plot->p,w,h,linesize,ofsX,ofsY,1,alpha,plot->form,buf);
        }
    }
}

static void lemG(PlotContext *plot, int n, int w, int h, int linesize, uint8_t *buf) {
    double count = plot->count;
    double speed = plot->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int ofsX = 0;
        int ofsY = 0;
        int j;
        double start = 0;
        int layers = plot->layers;
        for(j=0;j<layers;j++) {
            double s = start + plot->length*j/layers;
            double t = s + n*speed;
            drive(t,lemniskateG,plot->p,w,h,linesize,ofsX,ofsY,1,alpha,plot->form,buf);
        }
    }
}

static void lemB(PlotContext *plot, int n, int w, int h, int linesize, uint8_t *buf) {
    double count = plot->count;
    double speed = plot->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int ofsX = 0;
        int ofsY = 0;
        int j;
        double start = 0;
        int layers = plot->layers;
        for(j=0;j<layers;j++) {
            double s = start + plot->length*j/layers;
            double t = s + n*speed;
            drive(t,lemniskateB,plot->p,w,h,linesize,ofsX,ofsY,1,alpha,plot->form,buf);
        }
    }
}

static void epi(PlotContext *plot, int n, int w, int h, int linesize, uint8_t *buf) {
    double count = plot->count;
    double speed = plot->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int ofsX = 0;
        int ofsY = 0;
        int j;
        double start = 0;
        int layers = plot->layers;
        for(j=0;j<layers;j++) {
            double s = start + plot->length*j/layers;
            double t = s + n*speed;
            drive(t,epizycloid,plot->p,w,h,linesize,ofsX,ofsY,1,alpha,plot->form,buf);
        }
    }
}

static void lissGP(PlotContext *plot, int n, int w, int h, int linesize, uint8_t *buf) {
    double t = 0;
    double tmp = plot->p[1];
    double tmp5 = plot->p[4];
    double x = plot->nfunc(n,plot->np);
    plot->p[1]+=x;
    plot->p[4]+=x;
    while(t<plot->length) {
        int ofsX = plot->x;
        int ofsY = plot->y;
        drive(t,lissajousG,plot->p,w,h,linesize,ofsX,ofsY,1,0,plot->form,buf);
        t+=plot->speed;
    }
    drawNumber(w-200, h-15, plot->p[1], w, h, linesize, buf);
    plot->p[1] = tmp;
    plot->p[4] = tmp5;
}


static void lissQP(PlotContext *plot, int n, int w, int h, int linesize, uint8_t *buf) {
    double t = 0;
    double tmp = plot->p[1];
    double x = plot->nfunc(n,plot->np);
    plot->p[1]+=x;
    while(t<plot->length) {
        int ofsX = plot->x;
        int ofsY = plot->y;
        drive(t,lissajousQ,plot->p,w,h,linesize,ofsX,ofsY,1,0,plot->form,buf);
        t+=plot->speed;
    }
    drawNumber(w-200, h-15, plot->p[1], w, h, linesize, buf);
    plot->p[1] = tmp;
}

static void lissP(PlotContext *plot, int n, int w, int h, int linesize, uint8_t *buf) {
    double t = 0;
    double tmp = plot->p[1];
    double x = plot->nfunc(n,plot->np);
    plot->p[1]+=x;
    while(t<plot->length) {
        int ofsX = 0;
        int ofsY = 0;
        double alpha = 0;
        drive(t,lissajous,plot->p,w,h,linesize,ofsX,ofsY,1,alpha,plot->form,buf);
        t+=plot->speed;
    }
    drawNumber(w-200, h-15, plot->p[1], w, h, linesize, buf);
    plot->p[1] = tmp;
}

static void liss(PlotContext *plot, int n, int w, int h, int linesize, uint8_t *buf) {
    double count = plot->count;
    double speed = plot->speed;
    int k;
    double tmp = plot->p[1];
    double x = plot->nfunc(n,plot->np);
    plot->p[1]-=x;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int ofsX = 0;
        int ofsY = 0;
        int j;
        double start = 0;
        int layers = plot->layers;
        for(j=0;j<layers;j++) {
            double s = start + plot->length*j/layers;
            double t = s + n*speed;
            drive(t,lissajous,plot->p,w,h,linesize,ofsX,ofsY,1,alpha,plot->form,buf);
        }
    }
    plot->p[1] = tmp;
}

static void leg(PlotContext *plot, int n, int w, int h, int linesize, uint8_t *buf) {
    double count = plot->count;
    double speed = plot->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int ofsX = 0;
        int ofsY = 0;
        int j;
        double start = 0;
        int layers = plot->layers;
        for(j=0;j<layers;j++) {
            double s = start + plot->length*j/layers;
            double t = s + n*speed;
            drive(t,legendre,plot->p,w,h,linesize,ofsX,ofsY,1,alpha,plot->form,buf);
        }
    }
}

static void hypo(PlotContext *plot, int n, int w, int h, int linesize, uint8_t *buf) {
    double count = plot->count;
    double speed = plot->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int ofsX = 0;
        int ofsY = 0;
        int j;
        double start = 0;
        int layers = plot->layers;
        for(j=0;j<layers;j++) {
            double s = start + plot->length*j/layers;
            double t = s + n*speed;
            drive(t,hypocycloid,plot->p,w,h,linesize,ofsX,ofsY,1,alpha,plot->form,buf);
        }
    }
}


static void cbP(PlotContext *plot, int n, int w, int h, int linesize, uint8_t *buf) {
    double t = plot->p[39];
    double x = plot->nfunc(n,plot->np);
    double tmp = plot->p[2];
    plot->p[2]+=x;
    while(t<plot->length) {
        int ofsX = plot->x;
        int ofsY = plot->y;
        double alpha = 0;
        drive(t,cb,plot->p,w,h,linesize,ofsX,ofsY,1,alpha,plot->form,buf);
        t+=plot->speed;
    }
    drawNumber(w-200, h-15, plot->p[2], w, h, linesize, buf);
    plot->p[2] = tmp;
}

static void zero(PlotContext *plot, int n, int w, int h, int linesize, uint8_t *buf) {
}

static FFunc ffuncs[] = {
    {"circs",circs},
    {"kardioids",kardioids},
    {"lemG",lemG},
    {"lemB",lemB},
    {"epi",epi},
    {"liss",liss},
    {"lissP",lissP},
    {"lissQP",lissQP},
    {"lissGP",lissGP},
    {"cbP",cbP},
    {"leg",leg},
    {"hypo",hypo},
    {NULL,NULL}
};

static void (*getFunc(const char *name))(PlotContext *, int, int, int, int, uint8_t*) {
    int k=0;
    while(ffuncs[k].name) {
        if(!strcmp(name, ffuncs[k].name)) {
            return ffuncs[k].f;
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

static av_cold int plot_init(AVFilterContext *ctx)
{
    PlotContext *plot = ctx->priv;
    av_log(ctx, AV_LOG_INFO, "rgb=%d ofs=%d x=%f y=%f xn=%f yn=%f, layers=%d speed=%f length=%f count=%d dim=%d\n", 
            plot->is_rgb, plot->offset, plot->x, plot->y, plot->xn, plot->yn, plot->layers, plot->speed, plot->length, plot->count, plot->dim);
    if(plot->w) {
        plot->wfunc = getNFunc(plot->w);
        if(!plot->wfunc) {
            av_log(ctx, AV_LOG_WARNING, "function for wf not found %s\n", plot->w);
            plot->wfunc = constant;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for wf is %s", plot->w);
            logParameters(ctx,plot->wp,20);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for w\n");
        plot->wfunc = constant;
    }

    if(plot->n) {
        plot->nfunc = getNFunc(plot->n);
        if(!plot->nfunc) {
            av_log(ctx, AV_LOG_WARNING, "function for nf not found %s\n", plot->n);
            plot->nfunc = idn;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for nf is %s", plot->n);
            logParameters(ctx,plot->np,20);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for n\n");
        plot->nfunc = idn;
    }

    if(plot->f) {
        plot->ffunc = getFunc(plot->f);
        if(!plot->ffunc) {
            av_log(ctx, AV_LOG_WARNING, "function for f not found %s\n", plot->f);
            plot->ffunc = zero;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for f is %s\n", plot->f);
            logParameters(ctx,plot->p,40);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for f\n");
        plot->ffunc = zero;
    }

    plot->last = calloc(1920*1080,sizeof(uint8_t));

    return 0;
}

static int plot_query_formats(AVFilterContext *ctx)
{
    PlotContext *plot = ctx->priv;
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

    if (plot->is_rgb) {
        fmts_list = ff_make_format_list(rgb_pix_fmts);
    } else
        fmts_list = ff_make_format_list(yuv_pix_fmts);
    if (!fmts_list)
        return AVERROR(ENOMEM);
    return ff_set_common_formats(ctx, fmts_list);
}

static int plot_config_props(AVFilterLink *inlink)
{
    PlotContext *plot = inlink->dst->priv;
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(inlink->format);

    av_assert0(desc);

    plot->hsub = desc->log2_chroma_w;
    plot->vsub = desc->log2_chroma_h;
    plot->bps = desc->comp[0].depth;
    plot->planes = desc->nb_components;
    return 0;
}

typedef struct ThreadData {
    int height;
    int width;
    int linesize;
    int n;
    AVFrame *in;
} ThreadData;



static void dim(PlotContext *plot, AVFrame *in, int width, int height, int linesize) {
    uint8_t *src;
    uint8_t *dst;
    int x,y;
    for (y = 0; y < height; y++) {
        src = plot->last + linesize * y;
        dst = in->data[0] + linesize * y;
        for (x = 0; x < width; x++) {
            dst[x] = src[x]<plot->dim?0:src[x]-plot->dim;
        }
    }
}

static int plot_filter_frame(AVFilterLink *inlink, AVFrame *in)
{
    AVFilterContext *ctx = inlink->dst;
    PlotContext *plot = ctx->priv;
    AVFilterLink *outlink = inlink->dst->outputs[0];
    av_frame_make_writable(in);
    double n = inlink->frame_count_out;
    int plane;
    for (plane = 0; plane < plot->planes && in->data[plane]; plane++) {
        const int width = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(inlink->w, plot->hsub) : inlink->w;
        const int height = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(inlink->h, plot->vsub) : inlink->h;
        const int linesize = in->linesize[plane];
        int x, y;
        uint8_t *ptr;
        uint8_t *last;
        if(plane == 0) {
            dim(plot,in,width,height,linesize);
            plot->ffunc(plot, n, width, height, linesize, in->data[plane]);
            for (y = 0; y < height; y++) {
                ptr = in->data[plane] + linesize * y;
                last = plot->last + linesize * y;
                for (x = 0; x < width; x++) {
                    last[x] = ptr[x];
                }
            }
        } else {
            if(in->data[plane]) {
                for (y = 0; y < height; y++) {
                    ptr = in->data[plane] + linesize * y;
                    for (x = 0; x < width; x++) {
                        ptr[x] = 128;
                    }
                }
            }
        }
    }
    return ff_filter_frame(outlink, in);
}

static av_cold void plot_uninit(AVFilterContext *ctx)
{
    av_log(ctx, AV_LOG_INFO, "uninit\n");
    PlotContext *plot = ctx->priv;
    free(plot->last);
    
}

static const AVFilterPad plot_inputs[] = {
    {
        .name         = "default",
        .type         = AVMEDIA_TYPE_VIDEO,
        .config_props = plot_config_props,
        .filter_frame = plot_filter_frame,
    },
    { NULL }
};

static const AVFilterPad plot_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO,
    },
    { NULL }
};

AVFilter ff_vf_plot = {
    .name          = "plot",
    .description   = NULL_IF_CONFIG_SMALL("Apply generic equation to each pixel."),
    .priv_size     = sizeof(PlotContext),
    .init          = plot_init,
    .uninit        = plot_uninit,
    .query_formats = plot_query_formats,
    .inputs        = plot_inputs,
    .outputs       = plot_outputs,
    .priv_class    = &plot_class,
    .flags         = AVFILTER_FLAG_SUPPORT_TIMELINE_GENERIC | AVFILTER_FLAG_SLICE_THREADS,
};
