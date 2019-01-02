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
#include "libavutil/genutil.h"
#include "libavutil/plot.h"
#include "internal.h"
#include "drawutils.h"
#include <complex.h>
#include <math.h>

enum { Y, U, V, A };
typedef struct PlotContext {
    const AVClass *class;
    const char *f;
    double p[40];
    double _p[40];
    double (*nfunc[4])(int,double*);
    double np[4][10];
    double _np[4][10];
    const char *n[4];
    int nmod[4]; // modulo n
    void (*ffunc)(struct PlotContext*,int, AVFrame *in);
    int w,h;
    double fx;
    double fy;
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
    int dbg;
    double rot;
    double rotn;
    uint8_t *last;

    const char *mf;
    double complex (*cfunc)(double complex, double*, double);
    int offset;
    int is_packed_rgb;
    uint8_t rgba_map[4];
    char *rf[20];
 
    AVFilterContext *ctx;
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

    { "fx","fx",       OFFSET(fx),     AV_OPT_TYPE_DOUBLE, {.dbl = 100}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "fy","fy",       OFFSET(fy),     AV_OPT_TYPE_DOUBLE, {.dbl = 100}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "x","x",       OFFSET(x),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "y","y",       OFFSET(y),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "xn","xn",     OFFSET(xn),    AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "yn","yn",     OFFSET(yn),    AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "sat","sat",   OFFSET(sat),   AV_OPT_TYPE_DOUBLE, {.dbl =  1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "form","form", OFFSET(form),  AV_OPT_TYPE_INT,    {.i64 =  1}, 0,  20,  FLAGS },
    { "layers","l",  OFFSET(layers),AV_OPT_TYPE_INT,    {.i64 =  4}, 1,  INT_MAX,  FLAGS },
    { "count","c",   OFFSET(count), AV_OPT_TYPE_INT,    {.i64 =  1}, 1,  INT_MAX,  FLAGS },
    { "dim","d",     OFFSET(dim),   AV_OPT_TYPE_INT,    {.i64 =  4}, 0,  255,  FLAGS },
    { "speed","s",   OFFSET(speed), AV_OPT_TYPE_DOUBLE, {.dbl =  0.01}, 0,  DBL_MAX,  FLAGS },
    { "length","len",OFFSET(length), AV_OPT_TYPE_DOUBLE, {.dbl =  2*M_PI}, 0,  DBL_MAX,  FLAGS },
    { "rot","rot",   OFFSET(rot),   AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "rotn","rotn", OFFSET(rotn),  AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "mf",  "mf",   OFFSET(mf),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "ofs","ofs",OFFSET(offset), AV_OPT_TYPE_INT,    {.i64=0},    0,     INT_MAX,        FLAGS },
    { "dbg","dbg",OFFSET(dbg), AV_OPT_TYPE_INT,    {.i64=0},    0,     1,        FLAGS },
        { "n0",  "n0",   OFFSET(n[0]),     AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nm0","nm0",   OFFSET(nmod[0]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },
    { "n00","n00",   OFFSET(np[0][0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n01","n01",   OFFSET(np[0][1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n02","n02",   OFFSET(np[0][2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n03","n03",   OFFSET(np[0][3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n04","n04",   OFFSET(np[0][4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n05","n05",   OFFSET(np[0][5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n06","n06",   OFFSET(np[0][6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n07","n07",   OFFSET(np[0][7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n08","n08",   OFFSET(np[0][8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n09","n09",   OFFSET(np[0][9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n1",  "n1",   OFFSET(n[1]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nm1","nm1",   OFFSET(nmod[1]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },
    { "n10","n10",   OFFSET(np[1][0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n11","n11",   OFFSET(np[1][1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n12","n12",   OFFSET(np[1][2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n13","n13",   OFFSET(np[1][3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n14","n14",   OFFSET(np[1][4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n15","n15",   OFFSET(np[1][5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n16","n16",   OFFSET(np[1][6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n17","n17",   OFFSET(np[1][7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n18","n18",   OFFSET(np[1][8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n19","n19",   OFFSET(np[1][9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n2",  "n2",  OFFSET(n[2]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nm2","nm2",   OFFSET(nmod[2]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },
    { "n20","n20",   OFFSET(np[2][0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n21","n21",   OFFSET(np[2][1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n22","n22",   OFFSET(np[2][2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n23","n23",   OFFSET(np[2][3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n24","n24",   OFFSET(np[2][4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n25","n25",   OFFSET(np[2][5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n26","n26",   OFFSET(np[2][6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n27","n27",   OFFSET(np[2][7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n28","n28",   OFFSET(np[2][8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n29","n29",   OFFSET(np[2][9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n3",  "n3",   OFFSET(n[3]),      AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nm3","nm3",   OFFSET(nmod[3]),  AV_OPT_TYPE_INT,    {.i64=0},    0,        INT_MAX,        FLAGS },
    { "n30","n30",   OFFSET(np[3][0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n31","n31",   OFFSET(np[3][1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n32","n32",   OFFSET(np[3][2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n33","n33",   OFFSET(np[3][3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n34","n34",   OFFSET(np[3][4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n35","n35",   OFFSET(np[3][5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n36","n36",   OFFSET(np[3][6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n37","n37",   OFFSET(np[3][7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n38","n38",   OFFSET(np[3][8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n39","n39",   OFFSET(np[3][9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },


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

AVFILTER_DEFINE_CLASS(plot);

typedef struct FFunc {
    const char *name;
    void (*f)(PlotContext*, int,AVFrame *in);
} FFunc;


static void draw_number(int x, int y, double z, AVFrame *in, int packed) {
    uint8_t color[] = {255,255,255,255};
    if(packed) {
        av_plot_number_p(x,y,z,color,in->width,in->height,in->linesize[0],in->data[0]);
    } else {
        av_plot_number(x,y,z,in->width,in->height,in->linesize[0],in->data[0]);
    }
}

static void set_alpha(PlotContext *plot, AVFrame *in, double t, PFUNC(f), double a,double alpha,int ofsX, int ofsY) {
    double complex _z = f(t,plot->_p);
    double complex z = av_genutil_rotate(_z,alpha);
    int x = floor(plot->fx*creal(z) + plot->w/2) + plot->x + ofsX;
    int y = floor(plot->fy*cimag(z) + plot->h/2) + plot->y + ofsY;
    if(plot->is_packed_rgb) {
        av_plot_form_p(plot->form,x,y,a,plot->w,plot->h,in->linesize[0],in->data[0]+plot->rgba_map[A]);
    } else {
        av_plot_form(plot->form,x,y,a,plot->w,plot->h,in->linesize[3],in->data[3]);
    }
}


static void black_yuv(PlotContext *plot,AVFrame *in) {
    int x,y;
    const int width = in->width;
    const int height = in->height;
    uint8_t *ptr;
    
    if(plot->dim==0) {
        if(plot->is_packed_rgb) {
            for (y = 0; y < height; y++) {
                ptr = in->data[0] + in->linesize[0] * y + plot->rgba_map[A];
                for (x = 0; x < width; x+=4) {
                    ptr[x] = 0;
                }
            }
        } else {
            if(in->data[3]) {
                for (y = 0; y < height; y++) {
                    ptr = in->data[3] + in->linesize[3] * y;
                    for (x = 0; x < width; x++) {
                        ptr[x] = 0;
                    }
                }
            }
        }
     }
}

static void circs(PlotContext *plot, int n, AVFrame *in) {
    black_yuv(plot,in);
    double count = plot->count;
    double speed = plot->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int ofsX = floor(cos(alpha) * plot->_p[0]);
        int ofsY = floor(sin(alpha) * plot->_p[1]);
        double start = 2*M_PI*k/count;
        int j;
        int layers = plot->layers;
        for(j=0;j<layers;j++) {
            double s = start + plot->length*j/layers;
            double t = s + n*speed;
            set_alpha(plot,in,t,av_genutil_circ,1,0,ofsX,ofsY);
        }
    }
} 

static void curve(PlotContext *plot, PFUNC(f),int n, AVFrame *in) {
    black_yuv(plot,in);
    double count = plot->count;
    double speed = plot->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int j;
        double start = 0;
        int layers = plot->layers;
        for(j=0;j<layers;j++) {
            double s = start + plot->length*j/layers;
            double t = s + n*speed;
            set_alpha(plot,in,t,f,1,alpha,0,0);
        }
    }
}

static void curve2(PlotContext *plot, PFUNC(f),int n, AVFrame *in) {
    black_yuv(plot,in);
    double count = plot->count;
    double speed = plot->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int j;
        double start = 0;
        int layers = plot->layers;
        for(j=0;j<layers;j++) {
            double s = start + plot->length*j/layers;
            double t = s + n*speed;
            set_alpha(plot,in,t,f,1,alpha,0,0);
            set_alpha(plot,in,-t,f,1,alpha,0,0);
        }
    }
}


static void kardioids(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_kardio,n,in);
}

static void lemG(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_lemniskateG,n,in);
}

static void lemB(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_lemniskateB,n,in);
}

static void epi(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_epicycloid,n,in);
}

static void lissGP(PlotContext *plot, int n, AVFrame *in) {
    black_yuv(plot,in);
    double t = 0;
    while(t<plot->length) {
        int ofsX = plot->x;
        int ofsY = plot->y;
        set_alpha(plot,in,t,av_genutil_lissajousG,1,0,ofsX,ofsY);
        t+=plot->speed;
    }
    if(plot->dbg) draw_number(plot->w-200, plot->h-15, plot->_p[1], in, plot->is_packed_rgb);
}


static void lissQP(PlotContext *plot, int n, AVFrame *in) {
    black_yuv(plot,in);
    double t = 0;
    while(t<plot->length) {
        int ofsX = plot->x;
        int ofsY = plot->y;
        set_alpha(plot,in,t,av_genutil_lissajousQ,1,0,ofsX,ofsY);
        t+=plot->speed;
    }
    if(plot->dbg) draw_number(plot->w-200, plot->h-15, plot->_p[1], in, plot->is_packed_rgb);
}

static void lissP(PlotContext *plot, int n, AVFrame *in) {
    black_yuv(plot,in);
    double t = 0;
    while(t<plot->length) {
        int ofsX = 0;
        int ofsY = 0;
        set_alpha(plot,in,t,av_genutil_lissajous,1,0,ofsX,ofsY);
        t+=plot->speed;
    }
    if(plot->dbg) draw_number(plot->w-200, plot->h-15, plot->_p[1], in, plot->is_packed_rgb);
}

static void liss(PlotContext *plot, int n, AVFrame *in) {
    black_yuv(plot,in);
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
            set_alpha(plot,in,t,av_genutil_lissajous,1,alpha,ofsX,ofsY);
        }
    }
}

static void leg(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_legendre,n,in);
}

static void hypo(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_hypocycloid,n,in);
}

static void rhodo(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_rhodonea,n,in);
    if(plot->dbg) {
        draw_number(plot->w-100, plot->h-15, plot->_p[1], in, plot->is_packed_rgb);
        draw_number(plot->w-200, plot->h-15, plot->_p[0], in, plot->is_packed_rgb);
    }
}

static void nodal(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_nodal,n,in);
    if(plot->dbg) draw_number(plot->w-100, plot->h-15, plot->_p[0], in, plot->is_packed_rgb);
}

static void talbot(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_talbot,n,in);
    if(plot->dbg) draw_number(plot->w-100, plot->h-15, plot->_p[0], in, plot->is_packed_rgb);
}

static void folium(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_folium,n,in);
    if(plot->dbg) draw_number(plot->w-100, plot->h-15, plot->_p[0], in, plot->is_packed_rgb);
}

static void gielis(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_gielis,n,in);
    int k;
    if(plot->dbg) {
        for(k=1;k<5;k++) {
            draw_number(plot->w-k*100, plot->h-15, plot->_p[k-1], in, plot->is_packed_rgb);
        }
    }
}

static void super_spiral(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_super_spiral,n,in);
    int k;
    if(plot->dbg) {
        for(k=1;k<6;k++) {
            draw_number(plot->w-k*100, plot->h-15, plot->_p[k-1], in, plot->is_packed_rgb);
        }
    }
}

static void super_rose(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_super_rose,n,in);
    int k;
    if(plot->dbg) {
        for(k=1;k<6;k++) {
            draw_number(plot->w-k*100, plot->h-15, plot->_p[k-1], in, plot->is_packed_rgb);
        }
    }
}

static void cbP(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_cb,n,in);
}

static void epi_spiral(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_epi_spiral,n,in);
    int k;
    if(plot->dbg) {
        for(k=1;k<2;k++) {
            draw_number(plot->w-k*100, plot->h-15, plot->_p[k-1], in, plot->is_packed_rgb);
        }
    }
}

static void spiral(PlotContext *plot, int n, AVFrame *in) {
    curve2(plot,av_genutil_spiral,n,in);
    int k;
    if(plot->dbg) {
        for(k=1;k<4;k++) {
            draw_number(plot->w-k*100, plot->h-15, plot->_p[k-1], in, plot->is_packed_rgb);
        }
    }
}

static void atom_spiral(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_atom_spiral,n,in);
    int k;
    if(plot->dbg) {
        for(k=1;k<2;k++) {
            draw_number(plot->w-k*100, plot->h-15, plot->_p[k-1], in, plot->is_packed_rgb);
        }
    }
}

static void cotes_spiral(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_cotes_spiral,n,in);
    int k;
    if(plot->dbg) {
        for(k=1;k<2;k++) {
            draw_number(plot->w-k*100, plot->h-15, plot->_p[k-1], in, plot->is_packed_rgb);
        }
    }
}

static void sin_spiral(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_sin_spiral,n,in);
    int k;
    if(plot->dbg) {
        for(k=1;k<2;k++) {
            draw_number(plot->w-k*100, plot->h-15, plot->_p[k-1], in, plot->is_packed_rgb);
        }
    }
}

static void maclaurin(PlotContext *plot, int n, AVFrame *in) {
    curve(plot,av_genutil_maclaurin,n,in);
    int k;
    if(plot->dbg) {
        for(k=1;k<2;k++) {
            draw_number(plot->w-k*100, plot->h-15, plot->_p[k-1], in, plot->is_packed_rgb);
        }
    }
}

static void zero(PlotContext *plot, int n, AVFrame *in) {
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
    {"nodal",nodal},
    {"talbot",talbot},
    {"folium",folium},
    {"rhodo",rhodo},
    {"gielis",gielis},
    {"superspiral",super_spiral},
    {"superrose",super_rose},
    {"epispiral",epi_spiral},
    {"spiral",spiral},
    {"atomspiral",atom_spiral},
    {"cotesspiral",cotes_spiral},
    {"sinspiral",sin_spiral},
    {"maclaurin",maclaurin},
    {NULL,NULL}
};

static void (*get_func(const char *name))(PlotContext *, int, AVFrame *in) {
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

static void init_last(PlotContext *plot) {
    if(plot->dim <= 0) return;
    int x, y;
    uint8_t *last;
    int height = plot->h;
    int width = plot->w;
    for (y = 0; y < height; y++) {
        last = plot->last + width * y;
        for (x = 0; x < width; x++) {
             last[x] = 0;
        }
    }
}

static av_cold int plot_init(AVFilterContext *ctx)
{
    PlotContext *plot = ctx->priv;
    av_log(ctx, AV_LOG_INFO, "form=%d ofs=%d x=%f y=%f xn=%f yn=%f, layers=%d speed=%f length=%f count=%d dim=%d\n", 
            plot->form, plot->offset, plot->x, plot->y, plot->xn, plot->yn, plot->layers, plot->speed, plot->length, plot->count, plot->dim);
    int k;
    for(k=0;k<4;k++) {
        if(plot->n[k]) {
            plot->nfunc[k] = av_genutil_get_nfunc(plot->n[k]);
            if(!plot->nfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for n%d not found %s\n", k, plot->n[k]);
            } else {
                av_log(ctx, AV_LOG_INFO, "function for n%d is %s", k, plot->n[k]);
                logParameters(ctx,plot->np[k],10);
            }
        } else {
            av_log(ctx, AV_LOG_WARNING, "no function given for n%d\n",k);
        }
    }

    if(plot->f) {
        plot->ffunc = get_func(plot->f);
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
    init_last(plot);
    plot->ctx = ctx;
    return 0;
}

static int plot_query_formats(AVFilterContext *ctx)
{
    static const enum AVPixelFormat main_fmts[] = {
        AV_PIX_FMT_YUVA444P, AV_PIX_FMT_YUVA422P, AV_PIX_FMT_YUVA420P,
        AV_PIX_FMT_GBRAP,
        AV_PIX_FMT_RGBA, AV_PIX_FMT_BGRA, AV_PIX_FMT_ARGB, AV_PIX_FMT_ABGR,
        AV_PIX_FMT_NONE
    };

    AVFilterFormats *fmts_list;
    fmts_list = ff_make_format_list(main_fmts);
    if (!fmts_list)
        return AVERROR(ENOMEM);
    return ff_set_common_formats(ctx, fmts_list);
}

static int plot_config_props(AVFilterLink *inlink)
{
    PlotContext *plot = inlink->dst->priv;
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(inlink->format);

    av_assert0(desc);
    plot->w = inlink->w;
    plot->h = inlink->h;
    plot->is_packed_rgb =
        ff_fill_rgba_map(plot->rgba_map, inlink->format) >= 0 &&
        inlink->format != AV_PIX_FMT_GBRAP;
    av_log(inlink->dst, AV_LOG_INFO, "packed rgb %d\n", plot->is_packed_rgb);
    av_log(inlink->dst, AV_LOG_INFO, "A index %d\n", plot->rgba_map[A]);
 
    return 0;
}

static void dim(PlotContext *plot, AVFrame *out) {
    if(plot->dim <= 0) return;
    int plane = 3;
    uint8_t *src;
    uint8_t *dst;
    int x,y;
    int height = plot->h;
    int width = plot->w;
    if(plot->is_packed_rgb) {
        for (y = 0; y < height; y++) {
            src = plot->last + width * y;
            dst = out->data[0] + out->linesize[0] * y + plot->rgba_map[A];
            int z=0;
            for (x = 0; x < width; x++) {
                    dst[z] = src[x]<plot->dim?0:src[x]-plot->dim;
                    z+=4;
            }
        }
    } else {
        for (y = 0; y < height; y++) {
            src = plot->last + out->linesize[plane] * y;
            dst = out->data[plane] + out->linesize[plane] * y;
            for (x = 0; x < width; x++) {
                    dst[x] = src[x]<plot->dim?0:src[x]-plot->dim;
            }
        }
    }
}

static void copy0(PlotContext *plot, AVFrame *out) {
    if(plot->dim <= 0) return;
    int x, y;
    uint8_t *ptr;
    uint8_t *last;
    int height = plot->h;
    int width = plot->w;
    if(plot->is_packed_rgb) {
        for (y = 0; y < height; y++) {
            ptr = out->data[0] + out->linesize[0] * y + plot->rgba_map[A];
            last = plot->last + width * y;
            int z = 0;
            for (x = 0; x < width; x++) {
                 last[x] = ptr[z];
                 z+=4;
            }
        }
    } else {
        for (y = 0; y < height; y++) {
            ptr = out->data[3] + out->linesize[3] * y;
            last = plot->last + width * y;
            for (x = 0; x < width; x++) {
                 last[x] = ptr[x];
            }
        }
    }
}

static int modify(PlotContext *s, int frameNumber) {
    int j,i;
    const char *format = "%s %d %s %s %d";
    i=0;
    for(j=0;j<10;j++) {
        if(s->rf[j]) {
            char input[40];
            char target[4];
            char src[4];
            char mode[4];
            int index,m;
            strcpy(input, s->rf[j]);
            sscanf(input, format, target, &index, mode, src, &m);
            double value;
            switch(src[0]) {
                case 'n': {
                    value = s->nfunc[m](s->nmod[m]?frameNumber%s->nmod[m]:frameNumber,s->_np[m]);
                    break;
                }                             
            }
            switch(target[0]) {
                case 'p': {
                    if(mode[0] == 'o') s->_p[index] = value;
                    if(mode[0] == 'a') s->_p[index] = s->p[index] + value;
                    if(mode[0] == 's') s->_p[index] = s->p[index] - value;
                    break;
                }
                case 'n': {
                    int f = index/10;
                    int i = index % 10;
                    if(mode[0] == 'o') s->_np[f][i] = value;
                    if(mode[0] == 'a') s->_np[f][i] = s->np[f][i] + value;
                    if(mode[0] == 's') s->_np[f][i] = s->np[f][i] - value;
                }
            }
        }
        i++;
    }
    return 0;
}

static void copy_params(PlotContext *genvideo) {
    int k,j;
    for(k=0;k<40;k++) {
        genvideo->_p[k] = genvideo->p[k];
    }
    for(k=0;k<4;k++) {
        for(j=0;j<10;j++) {
            genvideo->_np[k][j] = genvideo->np[k][j];
        }
    }
}

static int plot_filter_frame(AVFilterLink *inlink, AVFrame *in)
{
    AVFilterContext *ctx = inlink->dst;
    //av_log(ctx, AV_LOG_INFO, "linesize %d\n", in->linesize[0]);
    PlotContext *plot = ctx->priv;
    AVFilterLink *outlink = inlink->dst->outputs[0];
    av_frame_make_writable(in);
    double n = inlink->frame_count_out;
    dim(plot,in);
    copy_params(plot);
    modify(plot,inlink->frame_count_out);
    plot->ffunc(plot, n, in);
    copy0(plot,in);
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
