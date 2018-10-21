/*
 * Copyright (c) 2018 Christian Bacher
 * 
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/**
 * @file
 * Javascript filter VERY SLOW - only usable for quick tests with small sizes
 * Based on the GEQ Filter
 */

#include "libavutil/avassert.h"
#include "libavutil/imgutils.h"
#include "libavutil/avstring.h"
#include "libavutil/opt.h"
#include "libavutil/duktape.h"
#include "libavutil/file.h"
#include "avfilter.h"
#include "formats.h"
#include "internal.h"
#include "video.h"

static const char *const var_names[] = { "P",   "W",   "H",   "L",   "N",   "T",   "SW",   "SH",        NULL };
enum                                   { VAR_P, VAR_W, VAR_H, VAR_L, VAR_N, VAR_T, VAR_SW, VAR_SH, VAR_VARS_NB };



typedef struct JSContext {
    const AVClass *class;
    char *jsfile;
    int is_rgb;
    duk_context *duk_ctxs[4]; // javascript context for each plane
    AVFrame *picref;            ///< current input buffer
    char *js;                   // the js code
    int hsub, vsub;             ///< chroma subsampling
    int planes;                 ///< number of planes
    int bps;

} JSContext;

#define OFFSET(x) offsetof(JSContext, x)
#define FLAGS AV_OPT_FLAG_FILTERING_PARAM|AV_OPT_FLAG_VIDEO_PARAM
static const AVOption js_options[] = {
    { "js",   "the js filter",   OFFSET(jsfile), AV_OPT_TYPE_STRING, {.str=NULL}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "rgb", "is rgb", OFFSET(is_rgb), AV_OPT_TYPE_INT, {.i64=1}, 0, 1, FLAGS },
    { NULL }
};

AVFILTER_DEFINE_CLASS(js);

static inline double getpix(void *priv, double x, double y, int plane)
{
    int xi, yi;
    JSContext *js = priv;
    AVFrame *picref = js->picref;
    const uint8_t *src = picref->data[plane];
    int linesize = picref->linesize[plane];
    const int w = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(picref->width,  js->hsub) : picref->width;
    const int h = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(picref->height, js->vsub) : picref->height;

    if (!src)
        return 0;

    xi = x = av_clipf(x, 0, w - 2);
    yi = y = av_clipf(y, 0, h - 2);

    x -= xi;
    y -= yi;

    if (js->bps > 8) {
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


static av_cold void declare_functions(duk_context *duk_ctx, JSContext *js, int plane) {
    static duk_ret_t (*p[])(duk_context *) = { lum_or_g, cb_or_b, cr_or_r, a };
    duk_push_c_function(duk_ctx, print, DUK_VARARGS);
    duk_put_global_string(duk_ctx, "print");
    duk_push_c_function(duk_ctx, lum_or_g, DUK_VARARGS);
    duk_put_global_string(duk_ctx, js->is_rgb?"g":"lum");
    duk_push_c_function(duk_ctx, cb_or_b, DUK_VARARGS);
    duk_put_global_string(duk_ctx, js->is_rgb?"b":"cb");
    duk_push_c_function(duk_ctx, cr_or_r, DUK_VARARGS);
    duk_put_global_string(duk_ctx, js->is_rgb?"r":"cr");
    duk_push_c_function(duk_ctx, a, DUK_VARARGS);
    duk_put_global_string(duk_ctx, "alpha");
    duk_push_c_function(duk_ctx, p[plane], DUK_VARARGS);
    duk_put_global_string(duk_ctx, "p");
}


static av_cold int compile_js(AVFilterContext *ctx, duk_context *duk_ctx, const char* program_body)
{
    int success = 0;

    av_log(ctx, AV_LOG_TRACE, "compiling %s\n", program_body);
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


static av_cold int load_js(AVFilterContext *ctx) {
    JSContext *js = ctx->priv;
    int err;
    uint8_t *textbuf;
    uint8_t *tmp;
    size_t textbuf_size;

    js->js = av_strdup("");
    if ((err = av_file_map(js->jsfile, &textbuf, &textbuf_size, 0, ctx)) < 0) {
        av_log(ctx, AV_LOG_ERROR,
               "The js file '%s' could not be read or is empty\n",
               js->jsfile);
        return err;
    }
    av_log(ctx, AV_LOG_TRACE, "loaded %s\n", textbuf);

    if (textbuf_size > SIZE_MAX - 1 || !(tmp = av_realloc(js->js, textbuf_size + 1))) {
        av_file_unmap(textbuf, textbuf_size);
        return AVERROR(ENOMEM);
    }
    js->js = tmp;
    memcpy(js->js, textbuf, textbuf_size);
    js->js[textbuf_size] = 0;
    av_file_unmap(textbuf, textbuf_size);

    av_log(ctx, AV_LOG_DEBUG, "loaded %s\n", js->js );
    return 0;
}

static av_cold duk_context *get_duk_context(AVFilterContext *ctx, JSContext *js, int plane) {
    duk_context *duk_ctx = duk_create_heap_default();
    duk_push_global_stash(duk_ctx);
    duk_push_pointer(duk_ctx, (void *)ctx);
    duk_put_prop_string(duk_ctx, -2, "av_filter_context");
    duk_pop(duk_ctx);
    declare_functions(duk_ctx, js, plane);
    return duk_ctx;
}

static av_cold int init(AVFilterContext *ctx)
{
    JSContext *js = ctx->priv;
    load_js(ctx);
    for(int k = 0; k<4; k++) {
        duk_context *duk_ctx = get_duk_context(ctx, js, k);
        compile_js(ctx, duk_ctx, js->js);
        js->duk_ctxs[k] = duk_ctx;
    }
    return 0;
}

static int js_query_formats(AVFilterContext *ctx)
{
    JSContext *js = ctx->priv;
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

    if (js->is_rgb) {
        fmts_list = ff_make_format_list(rgb_pix_fmts);
    } else
        fmts_list = ff_make_format_list(yuv_pix_fmts);
    if (!fmts_list)
        return AVERROR(ENOMEM);
    return ff_set_common_formats(ctx, fmts_list);
}

static int js_config_props(AVFilterLink *inlink)
{
    JSContext *js = inlink->dst->priv;
    const AVPixFmtDescriptor *desc = av_pix_fmt_desc_get(inlink->format);

    av_assert0(desc);

    js->hsub = desc->log2_chroma_w;
    js->vsub = desc->log2_chroma_h;
    js->bps = desc->comp[0].depth;
    js->planes = desc->nb_components;
    av_log(inlink->dst, AV_LOG_INFO, "function found %d %d\n", js->hsub, js->vsub);
    return 0;
}


static duk_int_t runJSFunction(AVFilterContext *ctx, duk_context *duk_ctx, const char* funcName, double *values, void *out, int out_len)
{
    duk_int_t returnVal;
    // Get a reference to the named JS function
    if (duk_get_global_string(duk_ctx, funcName))
    {
        // Function found, push the args
        av_log(ctx, AV_LOG_TRACE, "function found %s\n", funcName);
        int k = 0;
        for(;k<VAR_VARS_NB;k++) {
            duk_push_number(duk_ctx, values[k]);
        }
        duk_push_external_buffer(duk_ctx);
        duk_config_buffer(duk_ctx, -1, out, out_len);
        // Use pcall - this lets you catch and handle any errors
        if (duk_pcall(duk_ctx, VAR_VARS_NB + 1) != DUK_EXEC_SUCCESS)
        {
            // An error occurred - display a stack trace
            duk_get_prop_string(duk_ctx, -1, "stack");
            av_log(ctx, AV_LOG_ERROR, "Error in execution %s\n", duk_safe_to_string(duk_ctx, -1));
            returnVal = -1;
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
        returnVal = AVERROR(EINVAL);
    }
    duk_pop(duk_ctx); // pop result
    av_log(ctx, AV_LOG_DEBUG, "function returned %d %p\n", returnVal, duk_ctx);
    return returnVal;
}

typedef struct ThreadData {
    int height;
    int width;
    AVFrame *out;
    double n;
    double t;
} ThreadData;


static int slice_js_filter(AVFilterContext *ctx, void *arg, int jobnr, int nb_jobs)
{
    JSContext *js = ctx->priv;
    if(jobnr >= js->planes) return 0;
    ThreadData *td = arg;
    const int plane = jobnr;
    const int width = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(td->width, js->hsub) : td->width;
    const int height = (plane == 1 || plane == 2) ? AV_CEIL_RSHIFT(td->height, js->vsub) : td->height;
    const int linesize = td->out->linesize[plane];
    uint8_t *plane_data = td->out->data[plane];

    double values[VAR_VARS_NB];
    values[VAR_P] = jobnr;
    values[VAR_W] = width;
    values[VAR_H] = height;
    values[VAR_L] = linesize;
    values[VAR_N] = td->n;
    values[VAR_SW] = width/td->width;
    values[VAR_SH] = height/td->height;
    values[VAR_T] = td->t;
    runJSFunction(ctx, js->duk_ctxs[plane], "ev", values, plane_data, linesize * height );


    return 0;
}


static int filter_frame(AVFilterLink *inlink, AVFrame *in)
{
    AVFilterContext *ctx = inlink->dst;
    JSContext *js = ctx->priv;
    AVFilterLink *outlink = inlink->dst->outputs[0];
    AVFrame *out;
    const int nb_threads = ff_filter_get_nb_threads(ctx);
    //av_log(ctx, AV_LOG_INFO, "number of threads %d\n", nb_threads);

    av_log(ctx, AV_LOG_DEBUG, "filter fraem out: %d %d frameNr=%lld\n", outlink->w, outlink->h, inlink->frame_count_out);
    out = ff_get_video_buffer(outlink, outlink->w, outlink->h);
    if (!out) {
        av_frame_free(&in);
        return AVERROR(ENOMEM);
    }
    av_frame_copy_props(out, in);
    js->picref = in;

    ThreadData td;
    td.out = out;
    td.width = inlink->w;
    td.height = inlink->h;
    td.n = inlink->frame_count_out;
    td.t = in->pts == AV_NOPTS_VALUE ? NAN : in->pts * av_q2d(inlink->time_base);
    // run one thread for each plane, we need a seperate heap for every thread
    // if we would slice then we may have a big overhead in compiling the js for every thread
    // see unsolved issue https://github.com/svaarala/duktape/issues/770
    ctx->internal->execute(ctx, slice_js_filter, &td, NULL, 4);

    return ff_filter_frame(outlink, out);
    
}

static av_cold void uninit(AVFilterContext *ctx)
{
    int k;
    JSContext *js = ctx->priv;

    for (k = 0; k < 4; k++) {
         duk_destroy_heap(js->duk_ctxs[k]);
    }
}

static const AVFilterPad js_inputs[] = {
    {
        .name         = "default",
        .type         = AVMEDIA_TYPE_VIDEO,
        .config_props = js_config_props,
        .filter_frame = filter_frame,
    },
    { NULL }
};

static const AVFilterPad js_outputs[] = {
    {
        .name = "default",
        .type = AVMEDIA_TYPE_VIDEO,
    },
    { NULL }
};

AVFilter ff_vf_js = {
    .name          = "js",
    .description   = NULL_IF_CONFIG_SMALL("Detect and draw edge."),
    .priv_size     = sizeof(JSContext),
    .init          = init,
    .uninit        = uninit,
    .query_formats = js_query_formats,
    .inputs        = js_inputs,
    .outputs       = js_outputs,
    .priv_class    = &js_class,
    .flags         = AVFILTER_FLAG_SUPPORT_TIMELINE_GENERIC,
};
