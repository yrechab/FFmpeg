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

typedef struct PlotContext {
    const AVClass *class;
    const char *f;
    const char *ff;
    double p[40];
    const char *nf[20];

    double (*nfunc[20])(int,double*);
    double np[20][NMAXPARAMS];
    const char *n[4];
    int nmod[20]; // modulo n
    void (*ffunc)(GenutilFuncParams*,int,AVFrame *in);
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
    { "ff",  "ff",   OFFSET(ff),  AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
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
    { "nf0",  "nf0",   OFFSET(nf[0]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf1",  "nf1",   OFFSET(nf[1]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf2",  "nf2",   OFFSET(nf[2]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf3",  "nf3",   OFFSET(nf[3]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf4",  "nf4",   OFFSET(nf[4]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf5",  "nf5",   OFFSET(nf[5]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf6",  "nf6",   OFFSET(nf[6]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf7",  "nf7",   OFFSET(nf[7]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf8",  "nf8",   OFFSET(nf[8]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nf9",  "nf9",   OFFSET(nf[9]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },

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
    
    for(k=0;k<10;k++) {
        if(plot->nf[k]) {
            av_genutil_parse_nfunc(plot->nf[k],plot->np[k],&plot->nfunc[k]); 
            if(!plot->nfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for nf%d not found %s\n", k, plot->n[k]);
            } else {
                av_log(ctx, AV_LOG_INFO, "function for nf%d is [%s]", k, plot->nf[k]);
                logParameters(ctx,plot->np[k],10);
            }
        }
        if(k<4 && plot->n[k] && !plot->nfunc[k]) {
            plot->nfunc[k] = av_genutil_get_nfunc(plot->n[k]);
            if(!plot->nfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for nf%d not found %s\n", k, plot->n[k]);
            } else {
                av_log(ctx, AV_LOG_INFO, "function for n%d is %s", k, plot->n[k]);
                logParameters(ctx,plot->np[k],10);
            }
        }
    }

    if(plot->ff) {
        av_genutil_parse_ffunc(plot->ff,plot->p,&plot->ffunc); 
        if(!plot->ffunc) {
            av_log(ctx, AV_LOG_WARNING, "function for ff not found %s\n", plot->ff);
        } else {
            av_log(ctx, AV_LOG_INFO, "function for ff is [%s]", plot->ff);
            logParameters(ctx,plot->p,40);
        }
    } 
    if(plot->f && !plot->ffunc) {
        plot->ffunc = av_genutil_get_ffunc(plot->f);
        if(!plot->ffunc) {
            av_log(ctx, AV_LOG_WARNING, "function for f not found %s\n", plot->f);
            plot->ffunc = gzero;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for f is %s\n", plot->f);
            logParameters(ctx,plot->p,40);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for f\n");
    }
    if(!plot->ffunc) {
        plot->ffunc = gzero;
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

static void make_params(PlotContext *plot, GenutilFuncParams *params, int frame_number) {
    int k,j;
    for(k=0;k<40;k++) params->p[k] = plot->p[k];
    params->length = plot->length;
    params->w = plot->w;
    params->h = plot->h;
    params->fx = plot->fx;
    params->fy = plot->fy;
    params->x = plot->x;
    params->y = plot->y;
    params->rot = plot->rot;
    params->form = plot->form;
    params->speed = plot->speed;
    params->count = plot->count;
    params->layers = plot->layers;

    for(k=0;k<4;k++) params->rgba_map[k]=plot->rgba_map[k];
    params->is_packed_rgb = plot->is_packed_rgb;

    double np[10][NMAXPARAMS];
    for(k=0;k<10;k++) {
        for(j=0;j<NMAXPARAMS;j++) {
            np[k][j] = plot->np[k][j];
        }
    }

    const char *format = "%s %d %s %s %d";
    int i=0;
    for(j=0;j<10;j++) {
        if(plot->rf[j]) {
            char input[40];
            char target[4];
            char src[4];
            char mode[4];
            int index,m;
            strcpy(input, plot->rf[j]);
            sscanf(input, format, target, &index, mode, src, &m);
            double value;
            switch(src[0]) {
                case 'n': {
                    value = plot->nfunc[m](plot->nmod[m]?frame_number%plot->nmod[m]:frame_number,np[m]);
                    break;
                }                             
            }
            switch(target[0]) {
                case 'p': {
                    if(mode[0] == 'o') params->p[index] = value;
                    if(mode[0] == 'a') params->p[index] = plot->p[index] + value;
                    if(mode[0] == 's') params->p[index] = plot->p[index] - value;
                    break;
                }
                case 'n': {
                    int f = index/10;
                    int i = index % 10;
                    if(mode[0] == 'o') np[f][i] = value;
                    if(mode[0] == 'a') np[f][i] = plot->np[f][i] + value;
                    if(mode[0] == 's') np[f][i] = plot->np[f][i] - value;
                    break;
                }
                case 'x': {
                    if(mode[0] == 'o') params->x = floor(value);
                    if(mode[0] == 'a') params->x = floor(plot->x + value);
                    if(mode[0] == 's') params->x = floor(plot->x - value);
                    break;
                }
                case 'y': {
                    if(mode[0] == 'o') params->y = floor(value);
                    if(mode[0] == 'a') params->y = floor(plot->y + value);
                    if(mode[0] == 's') params->y = floor(plot->y - value);
                    break;
                }
                case 'w': {
                    if(mode[0] == 'o') params->fx = floor(value);
                    if(mode[0] == 'a') params->fx = floor(plot->fx + value);
                    if(mode[0] == 's') params->fx = floor(plot->fx - value);
                    break;
                }
                case 'h': {
                    if(mode[0] == 'o') params->fy = floor(value);
                    if(mode[0] == 'a') params->fy = floor(plot->fy + value);
                    if(mode[0] == 's') params->fy = floor(plot->fy - value);
                    break;
                }
                case 'r': {
                    if(mode[0] == 'o') params->rot = value;
                    if(mode[0] == 'a') params->rot = plot->rot + value;
                    if(mode[0] == 's') params->rot = plot->rot - value;
                    break;
                }
                case 'f': {
                    if(mode[0] == 'o') if(value<9) params->form = floor(value);
                    if(mode[0] == 'a') if(value<9) params->form = floor(plot->form + value);
                    if(mode[0] == 's') if(value<9) params->form = floor(plot->form - value);
                    break;
                }
                case 'l': {
                    if(mode[0] == 'o') params->length = value;
                    if(mode[0] == 'a') params->length = plot->length + value;
                    if(mode[0] == 's') params->length = plot->length - value;
                    break;
                }
            }
        }
        i++;
    }
 
}


static int plot_filter_frame(AVFilterLink *inlink, AVFrame *in)
{
    AVFilterContext *ctx = inlink->dst;
    srand(0);
    //av_log(ctx, AV_LOG_INFO, "linesize %d\n", in->linesize[0]);
    PlotContext *plot = ctx->priv;
    AVFilterLink *outlink = inlink->dst->outputs[0];
    av_frame_make_writable(in);
    black_yuv(plot,in);
    double n = inlink->frame_count_out;
    dim(plot,in);
    GenutilFuncParams params;
    double *p = calloc(40,sizeof(double));
    params.p = p;
    make_params(plot,&params,n);
    if(plot->dbg) debug(&params,n,in,plot->is_packed_rgb?-1:0);
    plot->ffunc(&params, n, in);
    copy0(plot,in);
    free(p);
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
