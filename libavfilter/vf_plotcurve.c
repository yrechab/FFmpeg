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
 * Video generation based on drawing 2d curves
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
    const char *ff;
    //void (*ffunc)(struct PlotCurveContext*,int,AVFrame *in);
    void (*ffunc)(GenutilFuncParams*,int,AVFrame *in);
    double p[40];
    char *s[40];
    const char *nf[20];
    double (*nfunc[20])(int,double*);
    double np[20][NMAXPARAMS];
    const char *n[4];
    int nmod[20]; // modulo n
    const char *cf[3];
    double (*cfunc[3])(double,double*);
    double cp[3][CMAXPARAMS];
    const char *c[3];
    int cmod; // YUV = 0, RGB = 1, HSV = 2
    uint8_t colors[10][4];
    double fx;
    double fy;
    double x;
    double y;
    double delta;
    double start;
    double length;
    int mode;
    int count;
    int form;
    int dim;
    double rot;
    int hsub, vsub;             ///< chroma subsampling
    int w, h;            
    int planes;                 ///< number of planes
    int is_rgb;
    int bps;
    int dbg;
    int offset;
    char *rf[20];

    AVFilterContext *ctx;
} PlotCurveContext;


#define OFFSET(x) offsetof(PlotCurveContext, x)
#define FLAGS AV_OPT_FLAG_VIDEO_PARAM|AV_OPT_FLAG_FILTERING_PARAM

static const AVOption plotcurve_options[] = {
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
    { "color0", "color0",  OFFSET(colors[0]),   AV_OPT_TYPE_COLOR, {.str="white"}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "color1", "color1",  OFFSET(colors[1]),   AV_OPT_TYPE_COLOR, {.str="white"}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "color2", "color2",  OFFSET(colors[2]),   AV_OPT_TYPE_COLOR, {.str="white"}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "color3", "color3",  OFFSET(colors[3]),   AV_OPT_TYPE_COLOR, {.str="white"}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "color4", "color4",  OFFSET(colors[4]),   AV_OPT_TYPE_COLOR, {.str="white"}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "color5", "color5",  OFFSET(colors[5]),   AV_OPT_TYPE_COLOR, {.str="white"}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "color6", "color6",  OFFSET(colors[6]),   AV_OPT_TYPE_COLOR, {.str="white"}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "color7", "color7",  OFFSET(colors[7]),   AV_OPT_TYPE_COLOR, {.str="white"}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "color8", "color8",  OFFSET(colors[8]),   AV_OPT_TYPE_COLOR, {.str="white"}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "color9", "color9",  OFFSET(colors[9]),   AV_OPT_TYPE_COLOR, {.str="white"}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "s0", "s0",  OFFSET(s[0]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "s1", "s1",  OFFSET(s[1]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "s2", "s2",  OFFSET(s[2]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "s3", "s3",  OFFSET(s[3]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "s4", "s4",  OFFSET(s[4]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "s5", "s5",  OFFSET(s[5]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "s6", "s6",  OFFSET(s[6]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "s7", "s7",  OFFSET(s[7]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "s8", "s8",  OFFSET(s[8]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "s9", "s9",  OFFSET(s[9]),   AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "fx","fx",       OFFSET(fx),     AV_OPT_TYPE_DOUBLE, {.dbl = 100}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "fy","fy",       OFFSET(fy),     AV_OPT_TYPE_DOUBLE, {.dbl = 100}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "x","x",       OFFSET(x),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "y","y",       OFFSET(y),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "delta","d",   OFFSET(delta),   AV_OPT_TYPE_DOUBLE, {.dbl =  1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "form","form", OFFSET(form),  AV_OPT_TYPE_INT,    {.i64 =  1}, 0,  20,  FLAGS },
    { "count","count", OFFSET(count),  AV_OPT_TYPE_INT,    {.i64 =  1}, 0,  200,  FLAGS },
    { "dim","d",     OFFSET(dim),   AV_OPT_TYPE_INT,    {.i64 =  4}, 1,  255,  FLAGS },
    { "length","len",OFFSET(length), AV_OPT_TYPE_DOUBLE, {.dbl =  2*M_PI}, 0,  DBL_MAX,  FLAGS },
    { "start","s",   OFFSET(start), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, 0,  DBL_MAX,  FLAGS },
    { "rot","rot",   OFFSET(rot),   AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
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
    { "n0",  "n0",   OFFSET(n[0]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
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
    { "cf0",  "cf0",   OFFSET(cf[0]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "cf1",  "cf1",   OFFSET(cf[1]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "cf2",  "cf2",   OFFSET(cf[2]),    AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
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
    { "mode","mode",OFFSET(mode), AV_OPT_TYPE_INT,    {.i64=0},    0,        100,        FLAGS },
    { "cmod","cmod",OFFSET(cmod), AV_OPT_TYPE_INT,    {.i64=0},    0,        2,        FLAGS },
    { "rgb","rgb",OFFSET(is_rgb), AV_OPT_TYPE_INT,    {.i64=0},    0,        1,        FLAGS },
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

AVFILTER_DEFINE_CLASS(plotcurve);

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
    av_log(ctx, AV_LOG_INFO, "rgb=%d ofs=%d x=%f y=%f delta=%f start=%f length=%f fx=%f fy=%f dbg=%d\n", 
            plotcurve->is_rgb, plotcurve->offset, plotcurve->x, plotcurve->y, 
            plotcurve->delta, plotcurve->start, plotcurve->length, plotcurve->fx,plotcurve->fy,plotcurve->dbg);
    int k;
    for(k=0;k<10;k++) {
        if(plotcurve->nf[k]) {
            av_genutil_parse_nfunc(plotcurve->nf[k],plotcurve->np[k],&plotcurve->nfunc[k]); 
            if(!plotcurve->nfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for nf%d not found %s\n", k, plotcurve->n[k]);
            } else {
                av_log(ctx, AV_LOG_INFO, "function for nf%d is [%s]", k, plotcurve->nf[k]);
                logParameters(ctx,plotcurve->np[k],10);
            }
        }
        if(k<4 && plotcurve->n[k] && !plotcurve->nfunc[k]) {
            plotcurve->nfunc[k] = av_genutil_get_nfunc(plotcurve->n[k]);
            if(!plotcurve->nfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for nf%d not found %s\n", k, plotcurve->n[k]);
            } else {
                av_log(ctx, AV_LOG_INFO, "function for n%d is %s", k, plotcurve->n[k]);
                logParameters(ctx,plotcurve->np[k],10);
            }
        }
    }

    for(k=0;k<3;k++) {
        if(plotcurve->cf[k]) {
            av_genutil_parse_cfunc(plotcurve->cf[k],plotcurve->cp[k],&plotcurve->cfunc[k]); 
            if(!plotcurve->cfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for cf%d not found %s\n", k, plotcurve->cf[k]);
            } else {
                av_log(ctx, AV_LOG_INFO, "function for cf%d is [%s]", k, plotcurve->cf[k]);
                logParameters(ctx,plotcurve->cp[k],10);
            }
        } 

        if(plotcurve->c[k] && !plotcurve->cfunc[k]) {
            plotcurve->cfunc[k] = av_genutil_get_cfunc(plotcurve->c[k]);
            if(!plotcurve->cfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for c%d not found %s\n", k, plotcurve->c[k]);
                plotcurve->cfunc[k] = czero;
            } else {
                av_log(ctx, AV_LOG_INFO, "function for c%d is %s", k, plotcurve->c[k]);
                logParameters(ctx,plotcurve->cp[k],10);
            }
        } 
    }
    
    for(k=0;k<10;k++) {
        av_log(ctx, AV_LOG_INFO, "color %d is %d %d %d\n", k, plotcurve->colors[k][0],plotcurve->colors[k][1],plotcurve->colors[k][2]);
    }

    if(plotcurve->ff) {
        av_genutil_parse_ffunc(plotcurve->ff,plotcurve->p,&plotcurve->ffunc); 
        if(!plotcurve->ffunc) {
            av_log(ctx, AV_LOG_WARNING, "function for ff not found %s\n", plotcurve->ff);
        } else {
            av_log(ctx, AV_LOG_INFO, "function for ff is [%s]", plotcurve->ff);
            logParameters(ctx,plotcurve->p,40);
        }
    } 
    if(plotcurve->f && !plotcurve->ffunc) {
        plotcurve->ffunc = av_genutil_get_ffunc(plotcurve->f);
        if(!plotcurve->ffunc) {
            av_log(ctx, AV_LOG_WARNING, "function for f not found %s\n", plotcurve->f);
            plotcurve->ffunc = gzero;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for f is %s\n", plotcurve->f);
            logParameters(ctx,plotcurve->p,40);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for f\n");
    }
    if(!plotcurve->ffunc) {
        plotcurve->ffunc = gzero;
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

static void convert_colors(PlotCurveContext *plotcurve,GenutilFuncParams *params) {
    int k;
    for(k=0;k<10;k++) {
        double r,g,b;
        r = ((double)plotcurve->colors[k][0])/255.0;
        g = ((double)plotcurve->colors[k][1])/255.0;
        b = ((double)plotcurve->colors[k][2])/255.0;
        if(plotcurve->is_rgb) {
            params->colors[k][0] = r;
            params->colors[k][1] = g;
            params->colors[k][2] = b;
        } else {
            params->colors[k][0] = r * 0.299000  + g * 0.587000  + b * 0.114000;
            params->colors[k][1] = r * -0.168736 + g * -0.331264 + b * 0.500000  + 0.5;
            params->colors[k][2] = r * 0.500000  + g * -.418688  + b * -0.081312 + 0.5;  
        }
    }

}

static void make_params(PlotCurveContext *plotcurve, GenutilFuncParams *params, int frame_number) {
    int k,j;
    for(k=0;k<40;k++) params->p[k] = plotcurve->p[k];
    for(k=0;k<10;k++) params->s[k] = plotcurve->s[k];
    for(k=0;k<3;k++) params->cfunc[k] = plotcurve->cfunc[k];
    for(k=0;k<3;k++) {
        for(j=0;j<CMAXPARAMS;j++) params->cp[k][j] = plotcurve->cp[k][j];
    }
    params->mode = plotcurve->mode;
    params->cmod = plotcurve->cmod;
    params->count = plotcurve->count;
    params->delta = plotcurve->delta;
    params->start = plotcurve->start;
    params->length = plotcurve->length;
    params->w = plotcurve->w;
    params->h = plotcurve->h;
    params->is_rgb = plotcurve->is_rgb;
    params->fx = plotcurve->fx;
    params->fy = plotcurve->fy;
    params->x = plotcurve->x;
    params->y = plotcurve->y;
    params->rot = plotcurve->rot;
    params->form = plotcurve->form;
    params->ctx = plotcurve->ctx;

    double np[10][NMAXPARAMS];
    for(k=0;k<10;k++) {
        for(j=0;j<NMAXPARAMS;j++) {
            np[k][j] = plotcurve->np[k][j];
        }
    }
    convert_colors(plotcurve,params);

    const char *format = "%s %d %s %s %d";
    int i=0;
    for(j=0;j<10;j++) {
        if(plotcurve->rf[j]) {
            char input[40];
            char target[4];
            char src[4];
            char mode[4];
            int index,m;
            strcpy(input, plotcurve->rf[j]);
            sscanf(input, format, target, &index, mode, src, &m);
            double value;
            switch(src[0]) {
                case 'n': {
                    value = plotcurve->nfunc[m](plotcurve->nmod[m]?frame_number%plotcurve->nmod[m]:frame_number,np[m]);
                    break;
                }                             
            }
            switch(target[0]) {
                case 'p': {
                    if(mode[0] == 'o') params->p[index] = value;
                    if(mode[0] == 'a') params->p[index] = plotcurve->p[index] + value;
                    if(mode[0] == 's') params->p[index] = plotcurve->p[index] - value;
                    break;
                }
                case 'n': {
                    int f = index/10;
                    int i = index % 10;
                    if(mode[0] == 'o') np[f][i] = value;
                    if(mode[0] == 'a') np[f][i] = plotcurve->np[f][i] + value;
                    if(mode[0] == 's') np[f][i] = plotcurve->np[f][i] - value;
                    break;
                }
                case 'c': {
                    int f = index/10;
                    int i = index % 10;
                    if(mode[0] == 'o') params->cp[f][i] = value;
                    if(mode[0] == 'a') params->cp[f][i] = plotcurve->cp[f][i] + value;
                    if(mode[0] == 's') params->cp[f][i] = plotcurve->cp[f][i] - value;
                    break;
                }
                case 'x': {
                    if(mode[0] == 'o') params->x = floor(value);
                    if(mode[0] == 'a') params->x = floor(plotcurve->x + value);
                    if(mode[0] == 's') params->x = floor(plotcurve->x - value);
                    break;
                }
                case 'y': {
                    if(mode[0] == 'o') params->y = floor(value);
                    if(mode[0] == 'a') params->y = floor(plotcurve->y + value);
                    if(mode[0] == 's') params->y = floor(plotcurve->y - value);
                    break;
                }
                case 'w': {
                    if(mode[0] == 'o') params->fx = value;
                    if(mode[0] == 'a') params->fx = plotcurve->fx + value;
                    if(mode[0] == 's') params->fx = plotcurve->fx - value;
                    break;
                }
                case 'h': {
                    if(mode[0] == 'o') params->fy = value;
                    if(mode[0] == 'a') params->fy = plotcurve->fy + value;
                    if(mode[0] == 's') params->fy = plotcurve->fy - value;
                    break;
                }
                case 'r': {
                    if(mode[0] == 'o') params->rot = value;
                    if(mode[0] == 'a') params->rot = plotcurve->rot + value;
                    if(mode[0] == 's') params->rot = plotcurve->rot - value;
                    break;
                }
                case 'f': {
                    if(mode[0] == 'o') if(value<9) params->form = floor(value);
                    if(mode[0] == 'a') if(value<9) params->form = floor(plotcurve->form + value);
                    if(mode[0] == 's') if(value<9) params->form = floor(plotcurve->form - value);
                    break;
                }
                case 'l': {
                    if(mode[0] == 'o') params->length = value;
                    if(mode[0] == 'a') params->length = plotcurve->length + value;
                    if(mode[0] == 's') params->length = plotcurve->length - value;
                    break;
                }
                case 's': {
                    if(mode[0] == 'o') params->start = value;
                    if(mode[0] == 'a') params->start = plotcurve->start + value;
                    if(mode[0] == 's') params->start = plotcurve->start - value;
                    break;
                }
            }
        }
        i++;
    }
 
}
static void black_yuv(PlotCurveContext *plotcurve, AVFrame *in) {
    int plane,x,y;
    //const int width = AV_CEIL_RSHIFT(in->width, plotcurve->hsub);
    //const int height = AV_CEIL_RSHIFT(in->height, plotcurve->vsub);
    const int width = in->width;
    const int height = in->height;
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

static int plotcurve_filter_frame(AVFilterLink *inlink, AVFrame *in)
{
    AVFilterContext *ctx = inlink->dst;
    srand(0);
    PlotCurveContext *plotcurve = ctx->priv;
    AVFilterLink *outlink = inlink->dst->outputs[0];
    av_frame_make_writable(in);
    black_yuv(plotcurve,in);
    double n = inlink->frame_count_out;
    GenutilFuncParams params;
    double *p = calloc(40,sizeof(double));
    params.p = p;
    make_params(plotcurve,&params,n);
    if(plotcurve->dbg) debug(&params,n,in,0);
    plotcurve->ffunc(&params,n+plotcurve->offset,in);
    free(p);
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
