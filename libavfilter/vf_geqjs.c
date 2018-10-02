/*
 * Copyright (C) 2006 Michael Niedermayer <michaelni@gmx.at>
 * Copyright (C) 2012 Clément Bœsch <u pkh me>
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
 * Generic equation change filter with the duktape JS engine instead of the av util evaluator
 * Originally written by Michael Niedermayer for the MPlayer project, and
 * ported by Clément Bœsch for FFmpeg.
 * JS Engine integration by Christian Bacher (in short expressions much slower but if it gets complex it may be usable)
 */
#include "libavutil/avassert.h"
#include "libavutil/avstring.h"
#include "libavutil/eval.h"
#include "libavutil/opt.h"
#include "libavutil/file.h"
#include "libavutil/pixdesc.h"
#include "libavutil/duktape.h"
#include "internal.h"

static const char *const var_names[] = {   "X",   "Y",   "W",   "H",   "N",   "SW",   "SH",   "T",        NULL };
enum                                   { VAR_X, VAR_Y, VAR_W, VAR_H, VAR_N, VAR_SW, VAR_SH, VAR_T, VAR_VARS_NB };

typedef struct GEQJSContext {
    const AVClass *class;
    int mode;
    duk_context *duk_ctxs[4]; // javascript context for each plane
    char *expr_str[4+3];        ///< expression strings for each plane
    AVFrame *picref;            ///< current input buffer
    double values[VAR_VARS_NB]; ///< expression values
    int hsub, vsub;             ///< chroma subsampling
    int planes;                 ///< number of planes
    int is_rgb;
    int bps;
    char *tmp;
} GEQJSContext;

enum { Y = 0, U, V, A, G, B, R };

#define OFFSET(x) offsetof(GEQJSContext, x)
#define FLAGS AV_OPT_FLAG_VIDEO_PARAM|AV_OPT_FLAG_FILTERING_PARAM

static const AVOption geqjs_options[] = {
    { "lum_expr",   "set luminance expression",   OFFSET(expr_str[Y]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "lum",        "set luminance expression",   OFFSET(expr_str[Y]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "cb_expr",    "set chroma blue expression", OFFSET(expr_str[U]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "cb",         "set chroma blue expression", OFFSET(expr_str[U]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "cr_expr",    "set chroma red expression",  OFFSET(expr_str[V]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "cr",         "set chroma red expression",  OFFSET(expr_str[V]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "alpha_expr", "set alpha expression",       OFFSET(expr_str[A]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "a",          "set alpha expression",       OFFSET(expr_str[A]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "red_expr",   "set red expression",         OFFSET(expr_str[R]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "r",          "set red expression",         OFFSET(expr_str[R]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "green_expr", "set green expression",       OFFSET(expr_str[G]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "g",          "set green expression",       OFFSET(expr_str[G]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "blue_expr",  "set blue expression",        OFFSET(expr_str[B]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "b",          "set blue expression",        OFFSET(expr_str[B]), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    {NULL},
};

AVFILTER_DEFINE_CLASS(geqjs);

static inline double getpix(void *priv, double x, double y, int plane)
{
    int xi, yi;
    GEQJSContext *geqjs = priv;
    AVFrame *picref = geqjs->picref;
    const uint8_t *src = picref->data[plane];
    int linesize = picref->linesize[plane];
    const int w = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(picref->width,  geqjs->hsub) : picref->width;
    const int h = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(picref->height, geqjs->vsub) : picref->height;

    if (!src)
        return 0;

    xi = x = av_clipf(x, 0, w - 2);
    yi = y = av_clipf(y, 0, h - 2);

    x -= xi;
    y -= yi;

    if (geqjs->bps > 8) {
        const uint16_t *src16 = (const uint16_t*)src;
        linesize /= 2;

        return (1-y)*((1-x)*src16[xi +  yi    * linesize] + x*src16[xi + 1 +  yi    * linesize])
              +   y *((1-x)*src16[xi + (yi+1) * linesize] + x*src16[xi + 1 + (yi+1) * linesize]);
    } else {
        return (1-y)*((1-x)*src[xi +  yi    * linesize] + x*src[xi + 1 +  yi    * linesize])
              +   y *((1-x)*src[xi + (yi+1) * linesize] + x*src[xi + 1 + (yi+1) * linesize]);
    }
}

//TODO: cubic interpolate
//TODO: keep the last few frames
static double lum(void *priv, double x, double y) { return getpix(priv, x, y, 0); }
static double  cb(void *priv, double x, double y) { return getpix(priv, x, y, 1); }
static double  cr(void *priv, double x, double y) { return getpix(priv, x, y, 2); }
static double alpha(void *priv, double x, double y) { return getpix(priv, x, y, 3); }


static AVFilterContext *getAVFilterContext(duk_context *duk_ctx) {
    duk_push_global_stash(duk_ctx);
    duk_get_prop_string(duk_ctx, -1, "av_filter_context");
    AVFilterContext *ctx = (AVFilterContext *)duk_to_pointer(duk_ctx, -1);
    duk_pop(duk_ctx);
    duk_pop(duk_ctx);
    return ctx;
}

// functions for usage inside the javscript

static duk_ret_t print(duk_context *duk_ctx) {
    duk_push_string(duk_ctx, " ");
    duk_insert(duk_ctx, 0);
    duk_join(duk_ctx, duk_get_top(duk_ctx) - 1);
    AVFilterContext *ctx = getAVFilterContext(duk_ctx); 
    av_log(ctx, AV_LOG_INFO, "%s\n", duk_safe_to_string(duk_ctx, -1));
    return 0;
}

static duk_ret_t lum_or_g(duk_context *duk_ctx) {
    duk_get_top(duk_ctx);  /* #args not needed */
    AVFilterContext* ctx = getAVFilterContext(duk_ctx);
    double x = duk_to_number(duk_ctx, 0);
    double y = duk_to_number(duk_ctx, 1);
    double res = lum(ctx->priv, x, y);
    duk_push_number(duk_ctx, res);
    return 1;  /* one return value */
}

static duk_ret_t cb_or_b(duk_context *duk_ctx) {
    duk_get_top(duk_ctx);  /* #args not needed */
    AVFilterContext* ctx = getAVFilterContext(duk_ctx);
    double x = duk_to_number(duk_ctx, 0);
    double y = duk_to_number(duk_ctx, 1);
    double res = cb(ctx->priv, x, y);
    duk_push_number(duk_ctx, res);
    return 1;  /* one return value */
}

static duk_ret_t cr_or_r(duk_context *duk_ctx) {
    duk_get_top(duk_ctx);  /* #args not needed */
    AVFilterContext* ctx = getAVFilterContext(duk_ctx);
    double x = duk_to_number(duk_ctx, 0);
    double y = duk_to_number(duk_ctx, 1);
    double res = cr(ctx->priv, x, y);
    duk_push_number(duk_ctx, res);
    return 1;  /* one return value */
}

static duk_ret_t a(duk_context *duk_ctx) {
    duk_get_top(duk_ctx);  /* #args not needed */
    AVFilterContext* ctx = getAVFilterContext(duk_ctx);
    double x = duk_to_number(duk_ctx, 0);
    double y = duk_to_number(duk_ctx, 1);
    double res = alpha(ctx->priv, x, y);
    duk_push_number(duk_ctx, res);
    return 1;  /* one return value */
}


static av_cold void declareFunctions(duk_context *duk_ctx, GEQJSContext *geqjs, int plane) {
    static duk_ret_t (*p[])(duk_context *) = { lum_or_g, cb_or_b, cr_or_r, a };
    duk_push_c_function(duk_ctx, print, DUK_VARARGS);
    duk_put_global_string(duk_ctx, "print");
    duk_push_c_function(duk_ctx, lum_or_g, DUK_VARARGS);
    duk_put_global_string(duk_ctx, geqjs->is_rgb?"g":"lum");
    duk_push_c_function(duk_ctx, cb_or_b, DUK_VARARGS);
    duk_put_global_string(duk_ctx, geqjs->is_rgb?"b":"cb");
    duk_push_c_function(duk_ctx, cr_or_r, DUK_VARARGS);
    duk_put_global_string(duk_ctx, geqjs->is_rgb?"r":"cr");
    duk_push_c_function(duk_ctx, a, DUK_VARARGS);
    duk_put_global_string(duk_ctx, "alpha");
    duk_push_c_function(duk_ctx, p[plane], DUK_VARARGS);
    duk_put_global_string(duk_ctx, "p");
}

static av_cold int compileJS(AVFilterContext *ctx, duk_context *duk_ctx, const char* program_body)
{
    int success = 0;

    av_log(ctx, AV_LOG_DEBUG, "compiling %s\n", program_body);
    // Compile the JS into bytecode
    if (duk_pcompile_string(duk_ctx, 0, program_body) != 0)
    {
        // Error in program code
        av_log(ctx, AV_LOG_ERROR, "Compile failed: ");
        av_log(ctx, AV_LOG_ERROR, "%s\n", duk_safe_to_string(duk_ctx, -1));
        success = AVERROR(EINVAL);
    }
    else
    {
        // Actually evaluate it - this will push the compiled code into the global scope
        duk_pcall(duk_ctx, 0);
    }
    duk_pop(duk_ctx);
    return success;
}

static int load_jsfile(AVFilterContext *ctx, const char *file_name)
{
    GEQJSContext *geqjs = ctx->priv;
    int err;
    uint8_t *textbuf;
    uint8_t *tmp;
    size_t textbuf_size;

    geqjs->tmp = av_strdup("");
    if ((err = av_file_map(file_name, &textbuf, &textbuf_size, 0, ctx)) < 0) {
        av_log(ctx, AV_LOG_ERROR,
               "The js file '%s' could not be read or is empty\n",
               file_name);
        return err;
    }
    av_log(ctx, AV_LOG_TRACE, "loaded %s\n", textbuf);

    if (textbuf_size > SIZE_MAX - 1 || !(tmp = av_realloc(geqjs->tmp, textbuf_size + 1))) {
        av_file_unmap(textbuf, textbuf_size);
        return AVERROR(ENOMEM);
    }
    geqjs->tmp = tmp;
    memcpy(geqjs->tmp, textbuf, textbuf_size);
    geqjs->tmp[textbuf_size] = 0;
    av_log(ctx, AV_LOG_TRACE, "after memcpy %s\n", geqjs->tmp);
    av_file_unmap(textbuf, textbuf_size);

    return 0;
}

static av_cold int compileJSFile(AVFilterContext *ctx, duk_context *duk_ctx, const char* file_name)
{
    GEQJSContext *geqjs = ctx->priv;
    int result = load_jsfile(ctx, file_name);
    if(!result) {
        av_log(ctx, AV_LOG_DEBUG, "loaded %s\n", geqjs->tmp );
        return compileJS(ctx, duk_ctx, geqjs->tmp);
    }
    return result;
}

static int endsWith (char* base, const char* str) {
    int blen = strlen(base);
    int slen = strlen(str);
    return (blen >= slen) && (0 == strcmp(base + blen - slen, str));
}

static av_cold duk_context *get_duk_context(AVFilterContext *ctx, GEQJSContext *geqjs, int plane) {
    duk_context *duk_ctx = duk_create_heap_default();

    duk_push_global_stash(duk_ctx);
    duk_push_pointer(duk_ctx, (void *)ctx);
    duk_put_prop_string(duk_ctx, -2, "av_filter_context");
    duk_pop(duk_ctx);
    declareFunctions(duk_ctx, geqjs, plane);
    char *arg = geqjs->expr_str[plane < 3 && geqjs->is_rgb ? plane+4 : plane];
    if(endsWith(arg, ".js")) {
        av_log(ctx, AV_LOG_DEBUG, "compiling file %s for plane %d\n", arg ,plane);
        if(compileJSFile(ctx, duk_ctx, arg)==0) {
            return duk_ctx; 
        } else {
            return NULL;
        }
    } else {
        int length = strlen(arg);
        char function_body[length+37];
        snprintf(function_body, length+37, "function ev(X,Y,W,H,N,SW,SH,T) { %s }", arg);
        av_log(ctx, AV_LOG_DEBUG, "compiling %s for plane %d\n", function_body ,plane);
        if(compileJS(ctx, duk_ctx, function_body)==0) {
            return duk_ctx; 
        } else {
            return NULL;
        }
    }
}

static av_cold int geqjs_init(AVFilterContext *ctx)
{
    GEQJSContext *geqjs = ctx->priv;
    int plane;

    if (!geqjs->expr_str[Y] && !geqjs->expr_str[G] && !geqjs->expr_str[B] && !geqjs->expr_str[R]) {
        av_log(ctx, AV_LOG_ERROR, "A luminance or RGB expression is mandatory\n");
        return AVERROR(EINVAL);
    }
    geqjs->is_rgb = !geqjs->expr_str[Y];

    if ((geqjs->expr_str[Y] || geqjs->expr_str[U] || geqjs->expr_str[V]) && (geqjs->expr_str[G] || geqjs->expr_str[B] || geqjs->expr_str[R])) {
        av_log(ctx, AV_LOG_ERROR, "Either YCbCr or RGB but not both must be specified\n");
        return AVERROR(EINVAL);
    }

    if (!geqjs->expr_str[U] && !geqjs->expr_str[V]) {
        /* No chroma at all: fallback on luma */
        geqjs->expr_str[U] = av_strdup(geqjs->expr_str[Y]);
        geqjs->expr_str[V] = av_strdup(geqjs->expr_str[Y]);
    } else {
        /* One chroma unspecified, fallback on the other */
        if (!geqjs->expr_str[U]) geqjs->expr_str[U] = av_strdup(geqjs->expr_str[V]);
        if (!geqjs->expr_str[V]) geqjs->expr_str[V] = av_strdup(geqjs->expr_str[U]);
    }

    if (!geqjs->expr_str[A]) {
        if(geqjs->bps == 8) {
            geqjs->expr_str[A] = av_strdup("return 255;");
        } else {
            geqjs->expr_str[A] = av_strdup("return 65535;");
        }
    }
    if (!geqjs->expr_str[G])
        geqjs->expr_str[G] = av_strdup("return g(X,Y);");
    if (!geqjs->expr_str[B])
        geqjs->expr_str[B] = av_strdup("return b(X,Y);");
    if (!geqjs->expr_str[R])
        geqjs->expr_str[R] = av_strdup("return r(X,Y);");

    if (geqjs->is_rgb ?
            (!geqjs->expr_str[G] || !geqjs->expr_str[B] || !geqjs->expr_str[R])
                    :
            (!geqjs->expr_str[U] || !geqjs->expr_str[V] || !geqjs->expr_str[A])) {
        return AVERROR(ENOMEM);
    }

    for (plane = 0; plane < 4; plane++) {
        duk_context *duk_ctx = get_duk_context(ctx, geqjs, plane);
        if(duk_ctx) {
            geqjs->duk_ctxs[plane] = duk_ctx;
        } else {
            return AVERROR(EINVAL);
        }
    }

    return 0;
}

static int geqjs_query_formats(AVFilterContext *ctx)
{
    GEQJSContext *geqjs = ctx->priv;
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

    if (geqjs->is_rgb) {
        fmts_list = ff_make_format_list(rgb_pix_fmts);
    } else
        fmts_list = ff_make_format_list(yuv_pix_fmts);
    if (!fmts_list)
        return AVERROR(ENOMEM);
    return ff_set_common_formats(ctx, fmts_list);
}

static int geqjs_config_props(AVFilterLink *inlink)
{
    GEQJSContext *geqjs = inlink->dst->priv;
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(inlink->format);

    av_assert0(desc);

    geqjs->hsub = desc->log2_chroma_w;
    geqjs->vsub = desc->log2_chroma_h;
    geqjs->bps = desc->comp[0].depth;
    geqjs->planes = desc->nb_components;
    return 0;
}

typedef struct ThreadData {
    int height;
    int width;
    AVFrame *out;
} ThreadData;

static duk_int_t runJSFunction(AVFilterContext *ctx, duk_context *duk_ctx, const char* funcName, double *values, int length)
{
    duk_int_t returnVal;
    // Get a reference to the named JS function
    if (duk_get_global_string(duk_ctx, funcName))
    {
        // Function found, push the args
        av_log(ctx, AV_LOG_TRACE, "function found %s\n", funcName);
        int k = 0;
        for(;k<length;k++) {
            duk_push_number(duk_ctx, values[k]);
        }
        // Use pcall - this lets you catch and handle any errors
        if (duk_pcall(duk_ctx, length) != DUK_EXEC_SUCCESS)
        {
            // An error occurred - display a stack trace
            duk_get_prop_string(duk_ctx, -1, "stack");
            av_log(ctx, AV_LOG_ERROR, "Error in execution %s\n", duk_safe_to_string(duk_ctx, -1));
            returnVal = 0;
        }
        else
        {
            // function executed successfully - get result
            returnVal = duk_get_int(duk_ctx, -1);
        }
    }
    else
    {
        av_log(ctx, AV_LOG_ERROR, "The function %s was not found in input js\n", funcName );
        returnVal = 0;
    }
    duk_pop(duk_ctx); // pop result
    av_log(ctx, AV_LOG_TRACE, "function returned %d %p\n", returnVal, duk_ctx);
    return returnVal;
}


static int slice_geqjs_filter(AVFilterContext *ctx, void *arg, int jobnr, int nb_jobs)
{
    GEQJSContext *geqjs = ctx->priv;
    if(jobnr >= geqjs->planes) return 0;
    ThreadData *td = arg;
    const int plane = jobnr;
    const int width = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(td->width, geqjs->hsub) : td->width;
    const int height = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(td->height, geqjs->vsub) : td->height;
    const int slice_start = 0;
    const int slice_end = height;
    int x, y;
    uint8_t *ptr;
    uint16_t *ptr16;
    const int linesize = td->out->linesize[plane];
    uint8_t *plane_data = td->out->data[plane];
    uint16_t *plane_data_16 = (uint16_t*)td->out->data[plane];

    double values[VAR_VARS_NB];
    values[VAR_W] = width;
    values[VAR_H] = height;
    values[VAR_N] = geqjs->values[VAR_N];
    values[VAR_SW] = width/td->width;
    values[VAR_SH] = height/td->height;
    values[VAR_T] = geqjs->values[VAR_T];

    if (geqjs->bps == 8) {
        for (y = slice_start; y < slice_end; y++) {
            ptr = plane_data + linesize * y;
            values[VAR_Y] = y;

            for (x = 0; x < width; x++) {
                values[VAR_X] = x;
                av_log(ctx, AV_LOG_TRACE, "8bit running js with args %d %f %f\n", plane, values[VAR_X], values[VAR_Y]);
                ptr[x] = runJSFunction(ctx, geqjs->duk_ctxs[plane], "ev", values, VAR_VARS_NB);
            }
        }
    }
    else {
        for (y = slice_start; y < slice_end; y++) {
            ptr16 = plane_data_16 + (linesize/2) * y;
            values[VAR_Y] = y;
            for (x = 0; x < width; x++) {
                values[VAR_X] = x;
                av_log(ctx, AV_LOG_TRACE, "running js with args %d %f %f\n", plane, values[VAR_X], values[VAR_Y]);
                ptr16[x] = runJSFunction(ctx, geqjs->duk_ctxs[plane], "ev", values, VAR_VARS_NB);
            }
        }
    }

    return 0;
}

static int geqjs_filter_frame(AVFilterLink *inlink, AVFrame *in)
{
    AVFilterContext *ctx = inlink->dst;
    GEQJSContext *geqjs = ctx->priv;
    AVFilterLink *outlink = inlink->dst->outputs[0];
    AVFrame *out;

    geqjs->values[VAR_N] = inlink->frame_count_out;
    geqjs->values[VAR_T] = in->pts == AV_NOPTS_VALUE ? NAN : in->pts * av_q2d(inlink->time_base);

    geqjs->picref = in;
    av_log(ctx, AV_LOG_DEBUG, "filter fraem out: %d %d frameNr=%f\n", outlink->w, outlink->h, geqjs->values[VAR_N]);
    out = ff_get_video_buffer(outlink, outlink->w, outlink->h);
    if (!out) {
        av_frame_free(&in);
        return AVERROR(ENOMEM);
    }
    av_frame_copy_props(out, in);
    ThreadData td;
    td.out = out;
    td.width = inlink->w;
    td.height = inlink->h;
    // run one thread for each plane, we need a seperate heap for every thread
    // if we would slice then we may have a big overhead in compiling the js for every thread
    // see unsolved issue https://github.com/svaarala/duktape/issues/770
    // however maybe a TODO to checkout to run with more contexts and slicing
    ctx->internal->execute(ctx, slice_geqjs_filter, &td, NULL, 4);

    av_frame_free(&geqjs->picref);
    return ff_filter_frame(outlink, out);
}

static av_cold void geqjs_uninit(AVFilterContext *ctx)
{
    av_log(ctx, AV_LOG_DEBUG, "cleaning up\n");
    GEQJSContext *geqjs = ctx->priv;
    for(int k=0;k<4;k++) {
        duk_destroy_heap(geqjs->duk_ctxs[k]);
    }
}

static const AVFilterPad geqjs_inputs[] = {
    {
        .name         = "default",
        .type         = AVMEDIA_TYPE_VIDEO,
        .config_props = geqjs_config_props,
        .filter_frame = geqjs_filter_frame,
    },
    { NULL }
};

static const AVFilterPad geqjs_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO,
    },
    { NULL }
};

AVFilter ff_vf_geqjs = {
    .name          = "geqjs",
    .description   = NULL_IF_CONFIG_SMALL("Apply generic equation to each pixel via javascrpit."),
    .priv_size     = sizeof(GEQJSContext),
    .init          = geqjs_init,
    .uninit        = geqjs_uninit,
    .query_formats = geqjs_query_formats,
    .inputs        = geqjs_inputs,
    .outputs       = geqjs_outputs,
    .priv_class    = &geqjs_class,
    .flags         = AVFILTER_FLAG_SUPPORT_TIMELINE_GENERIC | AVFILTER_FLAG_SLICE_THREADS,
};
