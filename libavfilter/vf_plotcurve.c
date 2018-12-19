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
 * Video generation based on a function plotcurveter
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
#include <complex.h>
#include <math.h>


typedef struct PlotCurveContext {
    const AVClass *class;
    const char *f;
    void (*ffunc)(struct PlotCurveContext*,int,AVFrame *in);
    double p[40];
    double (*nfunc[4])(int,double*);
    int nx[4];
    double np[4][10];
    const char *n[4];
    double (*cfunc[3])(double,double*);
    double cp[3][10];
    const char *c[3];
    int cmod; // YUV = 0, RGB = 1, HSV = 2
    double fx;
    double fy;
    double x;
    double xn;
    double y;
    double yn;
    double delta;
    double start;
    double length;
    int form;
    int dim;
    double rot;
    double rotn;
    int hsub, vsub;             ///< chroma subsampling
    int w, h;            
    int planes;                 ///< number of planes
    int is_rgb;
    int bps;
    int dbg;
    int offset;
    AVFilterContext *ctx;
} PlotCurveContext;


#define OFFSET(x) offsetof(PlotCurveContext, x)
#define FLAGS AV_OPT_FLAG_VIDEO_PARAM|AV_OPT_FLAG_FILTERING_PARAM

static const AVOption plotcurve_options[] = {
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
    { "delta","d",   OFFSET(delta),   AV_OPT_TYPE_DOUBLE, {.dbl =  1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "form","form", OFFSET(form),  AV_OPT_TYPE_INT,    {.i64 =  1}, 0,  20,  FLAGS },
    { "dim","d",     OFFSET(dim),   AV_OPT_TYPE_INT,    {.i64 =  4}, 1,  255,  FLAGS },
    { "length","len",OFFSET(length), AV_OPT_TYPE_DOUBLE, {.dbl =  2*M_PI}, 0,  DBL_MAX,  FLAGS },
    { "start","s",   OFFSET(start), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, 0,  DBL_MAX,  FLAGS },
    { "rot","rot",   OFFSET(rot),   AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "rotn","rotn", OFFSET(rotn),  AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "n0",  "n0",   OFFSET(n[0]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "nx0","nx0",   OFFSET(nx[0]),   AV_OPT_TYPE_INT,    {.i64 =  0}, 0,  39,  FLAGS },
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
    { "nx1","nx1",   OFFSET(nx[1]),   AV_OPT_TYPE_INT,    {.i64 =  1}, 0,  39,  FLAGS },
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
    { "nx2","nx2",   OFFSET(nx[2]),   AV_OPT_TYPE_INT,    {.i64 =  2}, 0,  39,  FLAGS },
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
    { "nx3","nx3",   OFFSET(nx[3]),   AV_OPT_TYPE_INT,    {.i64 =  3}, 0,  39,  FLAGS },
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
    { "c0",  "c0",   OFFSET(c[0]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "c00","c00",   OFFSET(cp[0][0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c01","c01",   OFFSET(cp[0][1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c02","c02",   OFFSET(cp[0][2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c03","c03",   OFFSET(cp[0][3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c04","c04",   OFFSET(cp[0][4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c05","c05",   OFFSET(cp[0][5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c06","c06",   OFFSET(cp[0][6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c07","c07",   OFFSET(cp[0][7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c08","c08",   OFFSET(cp[0][8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c09","c09",   OFFSET(cp[0][9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c1",  "c1",   OFFSET(c[1]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "c10","c10",   OFFSET(cp[1][0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c11","c11",   OFFSET(cp[1][1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c12","c12",   OFFSET(cp[1][2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c13","c13",   OFFSET(cp[1][3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c14","c14",   OFFSET(cp[1][4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c15","c15",   OFFSET(cp[1][5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c16","c16",   OFFSET(cp[1][6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c17","c17",   OFFSET(cp[1][7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c18","c18",   OFFSET(cp[1][8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c19","c19",   OFFSET(cp[1][9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c2",  "c2",   OFFSET(c[2]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "c20","c20",   OFFSET(cp[2][0]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c21","c21",   OFFSET(cp[2][1]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c22","c22",   OFFSET(cp[2][2]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c23","c23",   OFFSET(cp[2][3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c24","c24",   OFFSET(cp[2][4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c25","c25",   OFFSET(cp[2][5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c26","c26",   OFFSET(cp[2][6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c27","c27",   OFFSET(cp[2][7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c28","c28",   OFFSET(cp[2][8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c29","c29",   OFFSET(cp[2][9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "cmod","cmod",OFFSET(cmod), AV_OPT_TYPE_INT,    {.i64=0},    0,        2,        FLAGS },
    { "rgb","rgb",OFFSET(is_rgb), AV_OPT_TYPE_INT,    {.i64=0},    0,        1,        FLAGS },
    { "dbg","dbg",OFFSET(dbg), AV_OPT_TYPE_INT,    {.i64=0},    0,        1,        FLAGS },
    { "ofs","ofs",OFFSET(offset), AV_OPT_TYPE_INT,    {.i64=0},    0,     INT_MAX,        FLAGS },
    {NULL},
};

AVFILTER_DEFINE_CLASS(plotcurve);

typedef struct FFunc {
    const char *name;
    void (*f)(PlotCurveContext*, int,AVFrame *in);
} FFunc;

static void drawNumber(int x, int y, double z, AVFrame *in) {
    av_plot_number(x,y,z,in->width,in->height,in->linesize[0],in->data[0]);
}

static void drawPoint(PlotCurveContext *plotcurve, AVFrame *in, int plane, double t, PFUNC(f), double a) {
    double complex _z = f(t,plotcurve->p);
    double complex z = av_genutil_rotate(_z,plotcurve->rot);
    int _x = floor(plotcurve->fx*creal(z) + plotcurve->w/2) + plotcurve->x;
    int _y = floor(plotcurve->fy*cimag(z) + plotcurve->h/2) + plotcurve->y;
    int x,y,w,h;
    if((plane==1 || plane==2)&&(plotcurve->is_rgb==0)) {
        w = AV_CEIL_RSHIFT(plotcurve->w, plotcurve->hsub);
        x = AV_CEIL_RSHIFT(_x, plotcurve->hsub);
        h = AV_CEIL_RSHIFT(plotcurve->h, plotcurve->vsub);
        y = AV_CEIL_RSHIFT(_y, plotcurve->vsub);
    } else {
        x=_x;y=_y;w=plotcurve->w;h=plotcurve->h;
    }
    av_plot_form(plotcurve->form,x,y,a,w,h,in->linesize[plane],in->data[plane]);
}

static void blackYUV(PlotCurveContext *plotcurve, AVFrame *in) {
    int plane,x,y;
    const int width = AV_CEIL_RSHIFT(in->width, plotcurve->hsub);
    const int height = AV_CEIL_RSHIFT(in->height, plotcurve->vsub);
    uint8_t *ptr;
    for(plane=1;plane<=2;plane++) {
        if(in->data[plane]) {
            for (y = 0; y < height; y++) {
                ptr = in->data[plane] + in->linesize[plane] * y;
                for (x = 0; x < width; x++) {
                    ptr[x] = 128;
                }
            }
        }
    }
}


static void curve(PlotCurveContext *plotcurve, PFUNC(f),int n, AVFrame *in) {
    blackYUV(plotcurve,in);
    double t = plotcurve->start;
    int k;
    for(k=0;k<4;k++) {
        if(plotcurve->nfunc[k]) {
            plotcurve->p[plotcurve->nx[k]] = plotcurve->nfunc[k](n,plotcurve->np[k]);
            if(plotcurve->dbg) drawNumber(plotcurve->w-k*100-100, plotcurve->h-35, plotcurve->p[plotcurve->nx[k]], in);
        }
    }
    while(t<plotcurve->start+plotcurve->length) {
        double colors[3];
        av_genutil_get_color(plotcurve->cfunc, plotcurve->cp,t,plotcurve->cmod, plotcurve->is_rgb, colors);
        drawPoint(plotcurve,in,0,t,f,colors[0]);
        drawPoint(plotcurve,in,1,t,f,colors[1]);
        drawPoint(plotcurve,in,2,t,f,colors[2]);
        t+=plotcurve->delta;
    }
}

static void lissG(PlotCurveContext *plotcurve, int n, AVFrame *in) {
    curve(plotcurve,av_genutil_lissajousG,n,in);
}

static void lissQ(PlotCurveContext *plotcurve, int n, AVFrame *in) {
    curve(plotcurve,av_genutil_lissajousQ,n,in);
}

static void liss(PlotCurveContext *plotcurve, int n, AVFrame *in) {
    curve(plotcurve,av_genutil_lissajous,n,in);
}

static void leg(PlotCurveContext *plotcurve, int n, AVFrame *in) {
    curve(plotcurve,av_genutil_legendre,n,in);
}

static void hypo(PlotCurveContext *plotcurve, int n, AVFrame *in) {
    curve(plotcurve,av_genutil_hypocycloid,n,in);
}

static void cbP(PlotCurveContext *plotcurve, int n, AVFrame *in) {
    curve(plotcurve,av_genutil_cb,n,in);
}

static void sinoid(PlotCurveContext *plotcurve, int n, AVFrame *in) {
    curve(plotcurve,av_genutil_sinusoid,n,in);
}

static void capri(PlotCurveContext *plotcurve, int n, AVFrame *in) {
    curve(plotcurve,av_genutil_capricornoid,n,in);
}

static void scara(PlotCurveContext *plotcurve, int n, AVFrame *in) {
    curve(plotcurve,av_genutil_scarabaeus,n,in);
}

static void tro(PlotCurveContext *plotcurve, int n, AVFrame *in) {
    curve(plotcurve,av_genutil_trochoid,n,in);
}

static void circ(PlotCurveContext *plotcurve, int n, AVFrame *in) {
    curve(plotcurve,av_genutil_circoid,n,in);
}

static void zero(PlotCurveContext *plotcurve, int n, AVFrame *in) {
}

static FFunc ffuncs[] = {
    {"liss",liss},
    {"lissQ",lissQ},
    {"lissG",lissG},
    {"cbP",cbP},
    {"leg",leg},
    {"hypo",hypo},
    {"sinoid",sinoid},
    {"capri",capri},
    {"circ",circ},
    {"scara",scara},
    {"tro",tro},
    {NULL,NULL}
};

static void (*getFunc(const char *name))(PlotCurveContext*, int, AVFrame*) {
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

static av_cold int plotcurve_init(AVFilterContext *ctx)
{
    PlotCurveContext *plotcurve = ctx->priv;
    plotcurve->ctx = ctx;
    av_log(ctx, AV_LOG_INFO, "rgb=%d ofs=%d x=%f y=%f xn=%f yn=%f, delta=%f start=%f length=%f fx=%f fy=%f dbg=%d\n", 
            plotcurve->is_rgb, plotcurve->offset, plotcurve->x, plotcurve->y, plotcurve->xn, plotcurve->yn, 
            plotcurve->delta, plotcurve->start, plotcurve->length, plotcurve->fx,plotcurve->fy,plotcurve->dbg);
    int k;
    for(k=0;k<4;k++) {
        if(plotcurve->n[k]) {
            plotcurve->nfunc[k] = av_genutil_get_nfunc(plotcurve->n[k]);
            if(!plotcurve->nfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for nf%d not found %s\n", k, plotcurve->n[k]);
            } else {
                av_log(ctx, AV_LOG_INFO, "function for nf%d is %s", k, plotcurve->n[k]);
                logParameters(ctx,plotcurve->np[k],10);
            }
        } else {
            av_log(ctx, AV_LOG_WARNING, "no function given for n\n");
        }
    }

    for(k=0;k<3;k++) {
        if(plotcurve->c[k]) {
            plotcurve->cfunc[k] = av_genutil_get_cfunc(plotcurve->c[k]);
            if(!plotcurve->cfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for cf%d not found %s\n", k, plotcurve->c[k]);
                plotcurve->cfunc[k] = czero;
            } else {
                av_log(ctx, AV_LOG_INFO, "function for cf%d is %s", k, plotcurve->c[k]);
                logParameters(ctx,plotcurve->cp[k],10);
            }
        } else {
            av_log(ctx, AV_LOG_WARNING, "no function given for n\n");
            plotcurve->cfunc[k] = czero;
        }
    }


    if(plotcurve->f) {
        plotcurve->ffunc = getFunc(plotcurve->f);
        if(!plotcurve->ffunc) {
            av_log(ctx, AV_LOG_WARNING, "function for f not found %s\n", plotcurve->f);
            plotcurve->ffunc = zero;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for f is %s\n", plotcurve->f);
            logParameters(ctx,plotcurve->p,40);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for f\n");
        plotcurve->ffunc = zero;
    }

    return 0;
}

static int plotcurve_query_formats(AVFilterContext *ctx)
{
    PlotCurveContext *plotcurve = ctx->priv;
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

    if (plotcurve->is_rgb) {
        fmts_list = ff_make_format_list(rgb_pix_fmts);
    } else
        fmts_list = ff_make_format_list(yuv_pix_fmts);
    if (!fmts_list)
        return AVERROR(ENOMEM);
    return ff_set_common_formats(ctx, fmts_list);
}

static int plotcurve_config_props(AVFilterLink *inlink)
{
    PlotCurveContext *plotcurve = inlink->dst->priv;
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(inlink->format);

    av_assert0(desc);

    plotcurve->hsub = desc->log2_chroma_w;
    plotcurve->vsub = desc->log2_chroma_h;
    plotcurve->w = inlink->w;
    plotcurve->h = inlink->h;
    plotcurve->bps = desc->comp[0].depth;
    plotcurve->planes = desc->nb_components;
    av_log(plotcurve->ctx, AV_LOG_INFO, "%d %d %d %d\n", plotcurve->hsub, plotcurve->vsub, plotcurve->w, plotcurve->h);
    return 0;
}

typedef struct ThreadData {
    int height;
    int width;
    int linesize;
    int n;
    AVFrame *in;
} ThreadData;



static int plotcurve_filter_frame(AVFilterLink *inlink, AVFrame *in)
{
    AVFilterContext *ctx = inlink->dst;
    PlotCurveContext *plotcurve = ctx->priv;
    AVFilterLink *outlink = inlink->dst->outputs[0];
    av_frame_make_writable(in);
    double n = inlink->frame_count_out;
    plotcurve->ffunc(plotcurve,n+plotcurve->offset,in);
    return ff_filter_frame(outlink, in);
}

static av_cold void plotcurve_uninit(AVFilterContext *ctx)
{
}

static const AVFilterPad plotcurve_inputs[] = {
    {
        .name         = "default",
        .type         = AVMEDIA_TYPE_VIDEO,
        .config_props = plotcurve_config_props,
        .filter_frame = plotcurve_filter_frame,
    },
    { NULL }
};

static const AVFilterPad plotcurve_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO,
    },
    { NULL }
};

AVFilter ff_vf_plotcurve = {
    .name          = "plotcurve",
    .description   = NULL_IF_CONFIG_SMALL("Apply generic equation to each pixel."),
    .priv_size     = sizeof(PlotCurveContext),
    .init          = plotcurve_init,
    .uninit        = plotcurve_uninit,
    .query_formats = plotcurve_query_formats,
    .inputs        = plotcurve_inputs,
    .outputs       = plotcurve_outputs,
    .priv_class    = &plotcurve_class,
    .flags         = AVFILTER_FLAG_SUPPORT_TIMELINE_GENERIC | AVFILTER_FLAG_SLICE_THREADS,
};
