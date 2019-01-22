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
#include "libavutil/genutil.h"
#include "libavutil/plot.h"
#include "internal.h"
#include <complex.h>

typedef struct FracFuncParams {
    double x;
    double y;
    double w;
    double h;
    double fx;
    double fy;
    double rot;
    double *p;
    int ifcmode;
    long double complex (*ifunc)(long double complex,double*);
    double ip[10];
    double (*cfunc[3])(double,double*);
    double cp[3][10];
    int cmod;
    int is_rgb;
} FracFuncParams;

typedef struct FracContext {
    const AVClass *class;
    const char *f;
    void (*ffunc)(FracFuncParams *, AVFrame*,int,int,int);
    double p[40];
    const char *ifc;
    long double complex (*ifunc)(long double complex,double*);
    double ip[10];
    int ifcmode;
    const char *nf[10];
    double np[10][10];
    int nmod[10]; // modulo n
    double (*nfunc[10])(int,double*);
    const char *cf[3];
    double (*cfunc[3])(double,double*);
    double cp[3][10];
    int cmod;
    double x;
    double y;
    double w;
    double h;
    double fx;
    double fy;
    double rot;
    int hsub, vsub;             ///< chroma subsampling
    int planes;                 ///< number of planes
    int is_rgb;
    int bps;
    int offset;
    int dbg;
    char *rf[20];

} FracContext;


#define OFFSET(x) offsetof(FracContext, x)
#define FLAGS AV_OPT_FLAG_VIDEO_PARAM|AV_OPT_FLAG_FILTERING_PARAM

static const AVOption frac_options[] = {
    { "f",  "f",   OFFSET(f),  AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "x","x",       OFFSET(x),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "y","y",       OFFSET(y),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "fx","fx",   OFFSET(fx),   AV_OPT_TYPE_DOUBLE, {.dbl =  1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "fy","fy",   OFFSET(fy),   AV_OPT_TYPE_DOUBLE, {.dbl =  1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "rot","rot",   OFFSET(rot),   AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
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

AVFILTER_DEFINE_CLASS(frac);
typedef struct IFunc {
    const char *name;
    long double complex (*f)(long double complex,double*);
} IFunc;

static long double complex izero(long double complex z, double *p) {
    return 0;
}

static long double complex hopalong(long double complex z, double *p) {
    double a = p[0];
    double b = p[1];
    double c = p[2];
    long double x0 = creal(z);
    long double y0 = cimag(z);
    long double x =  y0 - SIGN(x0)*sqrt(fabsl(b*x0-c));
    long double y = a - x0;
    return x + I * y;
}

static long double complex mandel(long double complex z, double *p) {
    long double complex c = p[0]+ I * p[1];
    return z*z + c;
}

static long double complex exp4(long double complex z, double *p) {
    long double complex c = p[0]+ I * p[1];
    return cexpl(z/(c*c*c*c));
}

static long double complex exp3(long double complex z, double *p) {
    long double complex c = p[0]+ I * p[1];
    return cexpl((z*z*z)/(c*c*c));
}

static long double complex exp_2(long double complex z, double *p) {
    long double complex c = p[0]+ I * p[1];
    return cexpl((z*z-p[2]*z)/(c*c*c));
}

static long double complex cosinv(long double complex z, double *p) {
    long double complex c = p[0]+ I * p[1];
    return ccosl(z/c);
}

static IFunc ifuncs[] = {
    { "zero", izero },
    { "hopalong", hopalong },
    { "mandel", mandel },
    { "exp4", exp4 },
    { "exp3", exp3 },
    { "exp2", exp_2 },
    { "cosinv", cosinv },
    { NULL, NULL }
};

static long double complex (*get_ifunc(const char *name))(long double complex, double*) {
    int k=0;
    while(ifuncs[k].name) {
        if(!strcmp(name, ifuncs[k].name)) {
            return ifuncs[k].f;
        }
        k++;
    }
    return NULL;
}

static void parse_ifunc(const char *ff, double *p, long double complex (**f)(long double complex, double*)) {
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
            p[j-1] = atof(token);
        }
    }
    free(str);
}


typedef struct Func {
    const char *name;
    void (*f)(FracFuncParams *params, AVFrame *in, int x, int y, int n);
} Func;

static void zero(FracFuncParams *params, AVFrame *in, int x, int y, int n) {
}

static void sq(FracFuncParams *params, AVFrame *in, int x, int y, int n) {
    int plane;
    int x0 = x-params->w/2;
    int y0 = y-params->h/2;
    uint8_t *ptr;
    for(plane=0;plane<3;plane++) {
        ptr = in->data[plane] + in->linesize[plane] * y;
        ptr[x] = (sqrt(x0*x0+y0*y0)/(params->w/2))*n*params->p[plane];
    }
}

static void sqs(FracFuncParams *params, AVFrame *in, int x, int y, int n) {
    int plane;
    int x0 = x-params->w/2;
    int y0 = y-params->h/2;
    uint8_t *ptr;
    for(plane=0;plane<3;plane++) {
        ptr = in->data[plane] + in->linesize[plane] * y;
        ptr[x] = asin(sqrt(x0*x0+y0*y0)/(params->w/2))*n*params->p[plane];
    }
}

static void j(FracFuncParams *params, AVFrame *in, int x, int y, int n) {
    int plane;
    double len = params->p[0];
    int k;
    long double complex z = (x-params->w/2)*params->fx+params->x + I * (y-params->h/2)*params->fy+params->y;
    for(k=0;k<len;k++) {
        z = params->ifunc(z,params->ip);
        if(params->ifcmode == 0 && cabsl(z) > params->p[1]) {
            break;
        }
    }
    double colors[3];
    av_genutil_get_color(params->cfunc, params->cp,k,params->cmod, params->is_rgb, colors);
    uint8_t *ptr;
    for(plane=0;plane<3;plane++) {
        ptr = in->data[plane] + in->linesize[plane] * y;
        ptr[x] = (k>=len)?0:floor(colors[plane]*255);
    }
}

static void m(FracFuncParams *params, AVFrame *in, int x, int y, int n) {
    int plane;
    double len = params->p[0];
    int k;
    long double x0 = (x-params->w/2)*params->fx+params->x; 
    long double y0 = (y-params->h/2)*params->fy+params->y;
    double ip[] = { x0, y0 };
    long double complex z = 0;
    double max=0;
    for(k=0;k<len;k++) {
        z = params->ifunc(z,ip);
        if(params->ifcmode == 0 && cabsl(z) > params->p[1]) {
            max = k;
            break;
        }
        if(params->ifcmode == 2) {
            if(fabsl(creal(z)) > max) max = (fabsl(creal(z)));
        }
    }
    double colors[3];
    av_genutil_get_color(params->cfunc, params->cp,max,params->cmod, params->is_rgb, colors);
    uint8_t *ptr;
    for(plane=0;plane<3;plane++) {
        ptr = in->data[plane] + in->linesize[plane] * y;
        ptr[x] = (k>=len)?0:floor(colors[plane]*255);
    }
}

static void hopa(FracFuncParams *params, AVFrame *in, int x, int y, int n) {
    int plane;
    //int x0 = x-params->w/2;
    //int y0 = y-params->h/2;
    double len = params->p[0];
    double a = params->ip[0] + params->p[2] * (y+params->y)*params->fy;
    double b = params->ip[1] + params->p[3] * (x+params->x)*params->fx;
    double c = params->ip[2];
    double ip[] = { a, b, c, 0, 0 ,0, 0, 0, 0, 0 };
    int k;
    long double complex z = 0;
    double max=0;
    for(k=0;k<len;k++) {
        z = params->ifunc(z,ip);
        if(params->ifcmode == 0 && cabsl(z) > params->p[1]) {
            max = k;
            break;
        }
        if(params->ifcmode == 1 && (fabsl(creal(z)) > params->p[1])) {
            max = k;
            break;
        }
        if(params->ifcmode == 2) {
            if(fabsl(creal(z)) > max) max = (fabsl(creal(z)));
        }
    }
    double colors[3];
    av_genutil_get_color(params->cfunc, params->cp,max,params->cmod, params->is_rgb, colors);
    uint8_t *ptr;
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
    { "j", j },
    { "m", m },
    { NULL, NULL }
};

static void (*get_func(const char *name))(FracFuncParams*,AVFrame*,int,int,int) {
    int k=0;
    while(funcs[k].name) {
        if(!strcmp(name, funcs[k].name)) {
            return funcs[k].f;
        }
        k++;
    }
    return NULL;
}

static void parse_ffunc(const char *ff, double *p, void (**f)(FracFuncParams*,AVFrame*,int,int,int)) {
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

static av_cold int frac_init(AVFilterContext *ctx)
{
    int k;
    FracContext *frac = ctx->priv;
    av_log(ctx, AV_LOG_INFO, "rgb=%d ofs=%d x=%f y=%f fx=%f fy=%f\n", frac->is_rgb, frac->offset, frac->x, frac->y, frac->fx, frac->fy);

    for(k=0;k<10;k++) {
        if(frac->nf[k]) {
            av_genutil_parse_nfunc(frac->nf[k],frac->np[k],&frac->nfunc[k]); 
            if(!frac->nfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for nf%d not found %s\n", k, frac->nf[k]);
            } else {
                av_log(ctx, AV_LOG_INFO, "function for nf%d is [%s]", k, frac->nf[k]);
                logParameters(ctx,frac->np[k],10);
            }
        }
    }
    if(frac->f) {
        parse_ffunc(frac->f,frac->p,&frac->ffunc); 
        if(!frac->ffunc) {
            av_log(ctx, AV_LOG_WARNING, "function for ff not found %s\n", frac->f);
        } else {
            av_log(ctx, AV_LOG_INFO, "function for ff is [%s]", frac->f);
            logParameters(ctx,frac->p,40);
        }
    } 
    if(!frac->ffunc) {
        frac->ffunc = zero;
    }

    if(frac->ifc) {
        parse_ifunc(frac->ifc,frac->ip,&frac->ifunc); 
        if(!frac->ifunc) {
            av_log(ctx, AV_LOG_WARNING, "function for if not found %s\n", frac->ifc);
        } else {
            av_log(ctx, AV_LOG_INFO, "function for if is [%s]", frac->ifc);
            logParameters(ctx,frac->ip,10);
        }
    } 
    if(!frac->ifunc) {
        frac->ifunc = izero;
    }

    for(k=0;k<3;k++) {
        if(frac->cf[k]) {
            av_genutil_parse_cfunc(frac->cf[k],frac->cp[k],&frac->cfunc[k]); 
            if(!frac->cfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for cf%d not found %s\n", k, frac->cf[k]);
            } else {
                av_log(ctx, AV_LOG_INFO, "function for cf%d is [%s]", k, frac->cf[k]);
                logParameters(ctx,frac->cp[k],10);
            }
        } 
     }


    return 0;
}

static int frac_query_formats(AVFilterContext *ctx)
{
    FracContext *frac = ctx->priv;
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

    if (frac->is_rgb) {
        fmts_list = ff_make_format_list(rgb_pix_fmts);
    } else
        fmts_list = ff_make_format_list(yuv_pix_fmts);
    if (!fmts_list)
        return AVERROR(ENOMEM);
    return ff_set_common_formats(ctx, fmts_list);
}

static int frac_config_props(AVFilterLink *inlink)
{
    FracContext *frac = inlink->dst->priv;
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(inlink->format);

    av_assert0(desc);
    frac->w = inlink->w;
    frac->h = inlink->h;
    frac->hsub = desc->log2_chroma_w;
    frac->vsub = desc->log2_chroma_h;
    frac->bps = desc->comp[0].depth;
    frac->planes = desc->nb_components;
    return 0;
}

typedef struct ThreadData {
    int n;
    AVFrame *in;
    FracFuncParams *params;
} ThreadData;

static int slice_frac_filter(AVFilterContext *ctx, void *arg, int jobnr, int nb_jobs)
{
    FracContext *frac = ctx->priv;
    ThreadData *td = arg;
    const int slice_start = (td->in->height *  jobnr) / nb_jobs;
    const int slice_end = (td->in->height * (jobnr+1)) / nb_jobs;
    int x, y;
    
    for (y = slice_start; y < slice_end; y++) {
        for (x = 0; x < td->in->width; x++) {
            frac->ffunc(td->params,td->in,x,y, td->n + frac->offset);
        }
    }
    return 0;
}

static void make_params(FracContext *frac, FracFuncParams *params, int frame_number) {
    int k,j;
    for(k=0;k<40;k++) params->p[k] = frac->p[k];
    for(k=0;k<3;k++) params->cfunc[k] = frac->cfunc[k];
    for(k=0;k<3;k++) {
        for(j=0;j<10;j++) params->cp[k][j] = frac->cp[k][j];
    }
    params->cmod = frac->cmod;
    params->ifcmode = frac->ifcmode;
    for(j=0;j<10;j++) params->ip[j] = frac->ip[j];
    params->ifunc = frac->ifunc;
    params->w = frac->w;
    params->h = frac->h;
    params->is_rgb = frac->is_rgb;
    params->fx = frac->fx;
    params->fy = frac->fy;
    params->x = frac->x;
    params->y = frac->y;
    params->rot = frac->rot;

    double np[10][10];
    for(k=0;k<10;k++) {
        for(j=0;j<10;j++) {
            np[k][j] = frac->np[k][j];
        }
    }

    const char *format = "%s %d %s %s %d";
    int i=0;
    for(j=0;j<10;j++) {
        if(frac->rf[j]) {
            char input[40];
            char target[4];
            char src[4];
            char mode[4];
            int index,m;
            strcpy(input, frac->rf[j]);
            sscanf(input, format, target, &index, mode, src, &m);
            double value;
            switch(src[0]) {
                case 'n': {
                    value = frac->nfunc[m](frac->nmod[m]?frame_number%frac->nmod[m]:frame_number,np[m]);
                    break;
                }                             
            }
            switch(target[0]) {
                case 'p': {
                    if(mode[0] == 'o') params->p[index] = value;
                    if(mode[0] == 'a') params->p[index] = frac->p[index] + value;
                    if(mode[0] == 's') params->p[index] = frac->p[index] - value;
                    break;
                }
                case 'n': {
                    int f = index/10;
                    int i = index % 10;
                    if(mode[0] == 'o') np[f][i] = value;
                    if(mode[0] == 'a') np[f][i] = frac->np[f][i] + value;
                    if(mode[0] == 's') np[f][i] = frac->np[f][i] - value;
                    break;
                }
                case 'c': {
                    int f = index/10;
                    int i = index % 10;
                    if(mode[0] == 'o') params->cp[f][i] = value;
                    if(mode[0] == 'a') params->cp[f][i] = frac->cp[f][i] + value;
                    if(mode[0] == 's') params->cp[f][i] = frac->cp[f][i] - value;
                    break;
                }
                case 'i': {
                    if(mode[0] == 'o') params->ip[index] = value;
                    if(mode[0] == 'a') params->ip[index] = frac->ip[index] + value;
                    if(mode[0] == 's') params->ip[index] = frac->ip[index] - value;
                    break;
                }

                case 'x': {
                    if(mode[0] == 'o') params->x = floor(value);
                    if(mode[0] == 'a') params->x = floor(frac->x + value);
                    if(mode[0] == 's') params->x = floor(frac->x - value);
                    break;
                }
                case 'y': {
                    if(mode[0] == 'o') params->y = floor(value);
                    if(mode[0] == 'a') params->y = floor(frac->y + value);
                    if(mode[0] == 's') params->y = floor(frac->y - value);
                    break;
                }
                case 'w': {
                    if(mode[0] == 'o') params->fx = value;
                    if(mode[0] == 'a') params->fx = frac->fx + value;
                    if(mode[0] == 's') params->fx = frac->fx - value;
                    break;
                }
                case 'h': {
                    if(mode[0] == 'o') params->fy = value;
                    if(mode[0] == 'a') params->fy = frac->fy + value;
                    if(mode[0] == 's') params->fy = frac->fy - value;
                    break;
                }
                case 'r': {
                    if(mode[0] == 'o') params->rot = value;
                    if(mode[0] == 'a') params->rot = frac->rot + value;
                    if(mode[0] == 's') params->rot = frac->rot - value;
                    break;
                }
            }
        }
        i++;
    }
 
}

static void fdebug(FracFuncParams *params, int frame_number, AVFrame *in) {
    int k;
    for(k=0;k<4;k++) {
        av_genutil_draw_number(k*100+20, params->h-16, params->p[k], in, 0);
        av_genutil_draw_number(k*100+20, params->h-16, params->p[k], in, 1);
        av_genutil_draw_number(k*100+20, params->h-16, params->p[k], in, 2);
    }
    for(k=4;k<8;k++) {
        av_genutil_draw_number(k*100+20, params->h-16, params->ip[k-4], in, 0);
        av_genutil_draw_number(k*100+20, params->h-16, params->ip[k-4], in, 1);
        av_genutil_draw_number(k*100+20, params->h-16, params->ip[k-4], in, 2);
    }
    av_genutil_draw_int(params->w-50, params->h-16, frame_number, in, 0);
    av_genutil_draw_int(params->w-50, params->h-16, frame_number, in, 1);
    av_genutil_draw_int(params->w-50, params->h-16, frame_number, in, 2);
}


static int frac_filter_frame(AVFilterLink *inlink, AVFrame *in)
{
    AVFilterContext *ctx = inlink->dst;
    const int nb_threads = ff_filter_get_nb_threads(ctx);
    FracContext *frac = ctx->priv;
    AVFilterLink *outlink = inlink->dst->outputs[0];
    av_frame_make_writable(in);
    double n = inlink->frame_count_out;
    double *p = calloc(40,sizeof(double));
    FracFuncParams params;
    params.p = p;
    make_params(frac,&params,n+frac->offset);
    
    ThreadData td;
    td.in = in;
    td.n = n;
    td.params = &params;
    ctx->internal->execute(ctx, slice_frac_filter, &td, NULL, FFMIN(in->height, nb_threads));
    if(frac->dbg) fdebug(&params,n,in);
    free(p);
    return ff_filter_frame(outlink, in);
}

static av_cold void frac_uninit(AVFilterContext *ctx)
{
    av_log(ctx, AV_LOG_INFO, "uninit\n");
}

static const AVFilterPad frac_inputs[] = {
    {
        .name         = "default",
        .type         = AVMEDIA_TYPE_VIDEO,
        .config_props = frac_config_props,
        .filter_frame = frac_filter_frame,
    },
    { NULL }
};

static const AVFilterPad frac_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO,
    },
    { NULL }
};

AVFilter ff_vf_frac = {
    .name          = "frac",
    .description   = NULL_IF_CONFIG_SMALL("Apply generic equation to each pixel."),
    .priv_size     = sizeof(FracContext),
    .init          = frac_init,
    .uninit        = frac_uninit,
    .query_formats = frac_query_formats,
    .inputs        = frac_inputs,
    .outputs       = frac_outputs,
    .priv_class    = &frac_class,
    .flags         = AVFILTER_FLAG_SUPPORT_TIMELINE_GENERIC | AVFILTER_FLAG_SLICE_THREADS,
};
