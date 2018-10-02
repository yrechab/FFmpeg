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

typedef struct Func {
    const char *name;
    uint8_t (*f)(int x, int y, int w, int h, int n, double p1, double p2, double p3);
} Func;

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
    return 5 * (z - 2*I) * (z - 1+1*I) * (z - 1+1*I) / (z + 3-2*I) * (z + 3-2*I) * (z + 3-2*I);
}

static double complex id(double complex z) {
    return z;
}

static double complex expinv(double complex z) {
    return cexp(1/z);
}

/* ============================= FUNCS ===========================================*/

static uint8_t zero(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    return 0;
}

static uint8_t poly1re(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = (double)n/((double)w*(p1==0?1:p1));
    double complex z = x*f + (y*f) * I;
    return 128 + creal(poly1(z));
}
static uint8_t poly1img(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = (double)n/((double)w*(p1==0?1:p1));
    double complex z = x*f + (y*f) * I;
    return 128 + cimag(poly1(z));
}
static uint8_t poly1abs(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = (double)n/((double)w*(p1==0?1:p1));
    double complex z = x*f + (y*f) * I;
    return cabs(poly1(z));
}

static uint8_t idre(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = (double)n/((double)w*(p1==0?1:p1));
    double complex z = x*f + (y*f) * I;
    return 128 + creal(id(z));
}
static uint8_t idimg(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = (double)n/((double)w*(p1==0?1:p1));
    double complex z = x*f + (y*f) * I;
    return 128 + cimag(id(z));
}

static uint8_t idabs(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = (double)n/((double)w*(p1==0?1:p1));
    double complex z = x*f + (y*f) * I;
    return cabs(id(z));
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

static uint8_t polyg(double complex (*poly)(double complex z),int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly(x*f + (y*f) * I);
    return ((double)n+ 1) * p3 * cabs(z) * hsv_g(carg(z), p2, 1);
}

static uint8_t polyb(double complex (*poly)(double complex z),int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly(x*f + (y*f) * I);
    return ((double)n + 1 ) * p3 * cabs(z) * hsv_b(carg(z), p2, 1);
}

static uint8_t polyr(double complex (*poly)(double complex z),int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly(x*f + (y*f) * I);
    return ((double)n + 1) * p3 * cabs(z) * hsv_r(carg(z), p2, 1);
}

static uint8_t expinvg(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    return polyg(expinv,x,y,w,h,n,p1,p2,p3);
}

static uint8_t expinvb(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    return polyb(expinv,x,y,w,h,n,p1,p2,p3);
}

static uint8_t expinvr(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    return polyr(expinv,x,y,w,h,n,p1,p2,p3);
}

static uint8_t poly4g(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    return polyg(poly4,x,y,w,h,n,p1,p2,p3);
}

static uint8_t poly4b(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    return polyb(poly4,x,y,w,h,n,p1,p2,p3);
}

static uint8_t poly4r(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    return polyr(poly4,x,y,w,h,n,p1,p2,p3);
}

static uint8_t poly3g(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly3(x*f + (y*f) * I);
    return ((double)n+ 1 + 0.001*n*n) * p3 * cabs(z) * hsv_g(carg(z), p2, 1);
}

static uint8_t poly3b(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly3(x*f + (y*f) * I);
    return ((double)n + 1 + 0.001*n*n) * p3 * cabs(z) * hsv_b(carg(z), p2, 1);
}

static uint8_t poly3r(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly3(x*f + (y*f) * I);
    return ((double)n + 1 + 0.001*n*n) * p3 * cabs(z) * hsv_r(carg(z), p2, 1);
}

static uint8_t poly2g(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly2(x*f + (y*f) * I);
    return ((double)n+1) * p3 * cabs(z) * hsv_g(carg(z), p2, 1);
}

static uint8_t poly2b(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly2(x*f + (y*f) * I);
    return ((double)n+1) * p3 * cabs(z) * hsv_b(carg(z), p2, 1);
}

static uint8_t poly2r(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly2(x*f + (y*f) * I);
    return ((double)n+1) * p3 * cabs(z) * hsv_r(carg(z), p2, 1);
}

static uint8_t poly1g(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly1(x*f + (y*f) * I);
    return ((double)n+1) * p3 * cabs(z) * hsv_g(carg(z), p2, 1);
}

static uint8_t poly1b(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly1(x*f + (y*f) * I);
    return ((double)n+1) * p3 * cabs(z) * hsv_b(carg(z), p2, 1);
}

static uint8_t poly1r(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly1(x*f + (y*f) * I);
    return ((double)n+1) * p3 * cabs(z) * hsv_r(carg(z), p2, 1);
}

static uint8_t idrgbg(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double complex z = x + y * I;
    return hsv_g(carg(z),1, cabs(z)/(h/2));
}

static uint8_t idrgbb(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double complex z = x + y * I;
    return hsv_b(carg(z),1, cabs(z)/(h/2));
}

static uint8_t idrgbr(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double complex z = x + y * I;
    return hsv_r(carg(z),1, cabs(z)/(h/2));
}

static uint8_t sqg(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly0(x*f + (y*f) * I);
    return ((double)n+1) * p3 * cabs(z) * hsv_g(carg(z), p2, 1);
}

static uint8_t sqb(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly0(x*f + (y*f) * I);
    return ((double)n+1) * p3 * cabs(z) * hsv_b(carg(z), p2, 1);
}

static uint8_t sqr(int x, int y, int w, int h, int n, double p1, double p2, double p3) {
    double f = p1/(double)w;
    double complex z = poly0(x*f + (y*f) * I);
    return ((double)n+1) * p3 * cabs(z) * hsv_r(carg(z), p2, 1);
}


static Func funcs[] = {
    { "zero", zero },
    { "idre", idre },
    { "idimg", idimg },
    { "idabs", idabs },
    { "poly1re", poly1re },
    { "poly1img", poly1img },
    { "poly1abs", poly1abs },
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
    { "idrgbr", idrgbr },
    { "idrgbg", idrgbg },
    { "idrgbb", idrgbb },
    { "sqr", sqr },
    { "sqg", sqg },
    { "sqb", sqb },
    {NULL, NULL},
};

static uint8_t (*getFunc(const char *name))(int,int,int,int,int,double,double,double) {
    int k=0;
    while(funcs[k++].name) {
        if(!strcmp(name, funcs[k].name)) {
            return funcs[k].f;
        }
    }
    return NULL;
}
/* ================================== FILTER ========================================== */

typedef struct CGenContext {
    const AVClass *class;
    const char *f1;
    double p11;
    double p12;
    double p13;
    const char *f2;
    double p21;
    double p22;
    double p23;
    const char *f3;
    double p31;
    double p32;
    double p33;
    uint8_t (*funcs[3])(int,int,int,int,int,double,double,double);
    int hsub, vsub;             ///< chroma subsampling
    int planes;                 ///< number of planes
    int is_rgb;
    int bps;
} CGenContext;


#define OFFSET(x) offsetof(CGenContext, x)
#define FLAGS AV_OPT_FLAG_VIDEO_PARAM|AV_OPT_FLAG_FILTERING_PARAM

static const AVOption cgen_options[] = {
    { "f1",  "f1",   OFFSET(f1),  AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "p11","p11",   OFFSET(p11), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p12","p12",   OFFSET(p12), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p13","p13",   OFFSET(p13), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "f2",  "f2",   OFFSET(f2),  AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "p21","p21",   OFFSET(p21), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p22","p22",   OFFSET(p22), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p23","p23",   OFFSET(p23), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "f3",  "f3",   OFFSET(f3),  AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "p31","p31",   OFFSET(p31), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p32","p32",   OFFSET(p32), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "p33","p33",   OFFSET(p33), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "rgb","rgb",OFFSET(is_rgb), AV_OPT_TYPE_INT,    {.i64=0},    0,        1,        FLAGS },
    {NULL},
};

AVFILTER_DEFINE_CLASS(cgen);

static av_cold int cgen_init(AVFilterContext *ctx)
{
    CGenContext *cgen = ctx->priv;
    if(cgen->f1) {
        cgen->funcs[0] = getFunc(cgen->f1);
        if(!cgen->funcs[0]) {
            av_log(ctx, AV_LOG_WARNING, "function for f1 not found %s\n", cgen->f1);
            cgen->funcs[0] = zero;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for f1 is %s %f %f %f\n", cgen->f1, cgen->p11, cgen->p12, cgen->p13);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for f1\n");
        cgen->funcs[0] = zero;
    }

    if(cgen->f2) {
        cgen->funcs[1] = getFunc(cgen->f2);
        if(!cgen->funcs[1]) {
            av_log(ctx, AV_LOG_WARNING, "function not found %s\n", cgen->f2);
            cgen->funcs[1] = zero;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for f2 is %s %f %f %f\n", cgen->f2, cgen->p21, cgen->p22, cgen->p23);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for f2\n");
        cgen->funcs[1] = zero;
    }

    if(cgen->f3) {
        cgen->funcs[2] = getFunc(cgen->f3);
        if(!cgen->funcs[2]) {
            av_log(ctx, AV_LOG_WARNING, "function not found %s\n", cgen->f3);
            cgen->funcs[2] = zero;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for f3 is %s %f %f %f\n", cgen->f3, cgen->p31, cgen->p32, cgen->p33);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for f3\n");
        cgen->funcs[2] = zero;
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
                ptr[x] = cgen->funcs[0](x-width/2,y-height/2,width,height,td->n,cgen->p11,cgen->p12,cgen->p13);
            }
        }
    }
    if(plane == 1) {
        for (y = slice_start; y < slice_end; y++) {
            ptr = td->in->data[plane] + linesize * y;

            for (x = 0; x < width; x++) {
                ptr[x] = cgen->funcs[1](x-width/2,y-height/2,width,height,td->n,cgen->p21,cgen->p22,cgen->p23);
            }
        }
    }
    if(plane == 2) {
        for (y = slice_start; y < slice_end; y++) {
            ptr = td->in->data[plane] + linesize * y;

            for (x = 0; x < width; x++) {
                ptr[x] = cgen->funcs[2](x-width/2,y-height/2,width,height,td->n,cgen->p31,cgen->p32,cgen->p33);
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
