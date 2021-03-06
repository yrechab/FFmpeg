/*
 * Copyright (c) 2019 Christian Bacher
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

#include "libavutil/imgutils.h"
#include "libavutil/intreadwrite.h"
#include "libavutil/opt.h"
#include "libavutil/pixdesc.h"
#include "framesync.h"
#include "avfilter.h"
#include "formats.h"
#include "internal.h"
#include "video.h"

typedef struct MartinContext {
    const AVClass *class;
    const AVPixFmtDescriptor *desc;
    FFFrameSync fs;
    int radius;
    float factor;
    float threshold;
    int planes;

    int llimit;
    int hlimit;
    int nb_inputs;
    int nb_frames;

    int depth;
    int nb_planes;
    int linesize[4];
    int height[4];

    AVFrame **frames;
    AVFrame **index_frames;
} MartinContext;

static int query_formats(AVFilterContext *ctx)
{
    static const enum AVPixelFormat pixel_fmts[] = {
        AV_PIX_FMT_GRAY8, AV_PIX_FMT_GRAY9,
        AV_PIX_FMT_GRAY10, AV_PIX_FMT_GRAY12, AV_PIX_FMT_GRAY14,
        AV_PIX_FMT_GRAY16,
        AV_PIX_FMT_YUV410P, AV_PIX_FMT_YUV411P,
        AV_PIX_FMT_YUV420P, AV_PIX_FMT_YUV422P,
        AV_PIX_FMT_YUV440P, AV_PIX_FMT_YUV444P,
        AV_PIX_FMT_YUVJ420P, AV_PIX_FMT_YUVJ422P,
        AV_PIX_FMT_YUVJ440P, AV_PIX_FMT_YUVJ444P,
        AV_PIX_FMT_YUVJ411P,
        AV_PIX_FMT_YUV420P9, AV_PIX_FMT_YUV422P9, AV_PIX_FMT_YUV444P9,
        AV_PIX_FMT_YUV420P10, AV_PIX_FMT_YUV422P10, AV_PIX_FMT_YUV444P10,
        AV_PIX_FMT_YUV440P10,
        AV_PIX_FMT_YUV444P12, AV_PIX_FMT_YUV422P12, AV_PIX_FMT_YUV420P12,
        AV_PIX_FMT_YUV440P12,
        AV_PIX_FMT_YUV444P14, AV_PIX_FMT_YUV422P14, AV_PIX_FMT_YUV420P14,
        AV_PIX_FMT_YUV420P16, AV_PIX_FMT_YUV422P16, AV_PIX_FMT_YUV444P16,
        AV_PIX_FMT_GBRP, AV_PIX_FMT_GBRP9, AV_PIX_FMT_GBRP10,
        AV_PIX_FMT_GBRP12, AV_PIX_FMT_GBRP14, AV_PIX_FMT_GBRP16,
        AV_PIX_FMT_NONE
    };
    AVFilterFormats *formats = ff_make_format_list(pixel_fmts);
    if (!formats)
        return AVERROR(ENOMEM);
    return ff_set_common_formats(ctx, formats);
}

typedef struct ThreadData {
    AVFrame **in, *out, **index;
    int xoff,yoff;
} ThreadData;

static int martin_frame(AVFilterContext *ctx, void *arg, int jobnr, int nb_jobs)
{
    MartinContext *s = ctx->priv;
    ThreadData *td = arg;
    AVFrame **in = td->in;
    AVFrame **idx = td->index;
    AVFrame *out = td->out;
    //const int nb_inputs = s->nb_inputs;
    //const int depth = s->depth;
    int x, y;

    const int slice_start = (s->height[0] * jobnr) / nb_jobs;
    const int slice_end = (s->height[0] * (jobnr+1)) / nb_jobs;
    uint8_t *dst0 = out->data[0] + slice_start * out->linesize[0];
    uint8_t *dst1 = out->data[1] + slice_start * out->linesize[1];
    uint8_t *dst2 = out->data[2] + slice_start * out->linesize[2];
    for (y = slice_start; y < slice_end; y++) {
      for (x = 0; x < s->linesize[0]; x++) {
	 int xi = x+td->xoff;
	 int yi = y+td->yoff;
	 uint8_t index=0;
	 if(xi>=0 && xi < idx[0]->width && yi>=0 && yi < idx[0]->height)
	 {
	   uint8_t i = idx[0]->data[0][yi * idx[0]->linesize[0] + xi];
	   float p = (float)i/256.0;
	   index = p*s->nb_frames;
	 }
         dst0[x] = in[index]->data[0][y * in[index]->linesize[0] + x];
         dst1[x] = in[index]->data[1][y * in[index]->linesize[1] + x];
         dst2[x] = in[index]->data[2][y * in[index]->linesize[2] + x];
      }
      dst0 += out->linesize[0];
      dst1 += out->linesize[1];
      dst2 += out->linesize[2];
    }
    return 0;
}
static int config_output(AVFilterLink *outlink)
{
    AVFilterContext *ctx = outlink->src;
    av_log(ctx, AV_LOG_DEBUG, "#################### config_output\n");
    MartinContext *s = ctx->priv;
    AVFilterLink *inlink = ctx->inputs[0];
    int ret;

    s->desc = av_pix_fmt_desc_get(outlink->format);
    if (!s->desc)
        return AVERROR_BUG;
    s->nb_planes = av_pix_fmt_count_planes(outlink->format);
    s->depth = s->desc->comp[0].depth;

    if ((ret = av_image_fill_linesizes(s->linesize, inlink->format, inlink->w)) < 0)
        return ret;

    s->height[1] = s->height[2] = AV_CEIL_RSHIFT(inlink->h, s->desc->log2_chroma_h);
    s->height[0] = s->height[3] = inlink->h;
    if ((ret = ff_framesync_init_dualinput(&s->fs, ctx)) < 0)
            return ret;
    ff_framesync_configure(&s->fs);
    return 0;
}

static av_cold void uninit(AVFilterContext *ctx)
{
    MartinContext *s = ctx->priv;
    int i;
    ff_framesync_uninit(&s->fs);
    if (s->frames) {
        for (i = 0; i < s->nb_frames; i++) {
           av_frame_free(&s->frames[i]);
           av_frame_free(&s->index_frames[i]);
	}
    }
    av_freep(&s->frames);
    av_freep(&s->index_frames);
}

static AVFrame *martin_filter_frame(AVFilterContext *ctx, AVFrame *in, AVFrame *index)
{
    av_log(ctx, AV_LOG_DEBUG, "#################### filterframe\n");
    AVFilterLink *outlink = ctx->outputs[0];
    MartinContext *s = ctx->priv;
    ThreadData td;
    AVFrame *out;
    if (s->nb_frames < s->nb_inputs) {
        s->frames[s->nb_frames] = av_frame_clone(in);
        s->index_frames[s->nb_frames] = av_frame_clone(index);
        s->nb_frames++;
    	av_log(ctx, AV_LOG_DEBUG, "#################### storing %d\n",s->nb_frames);
        //return NULL;
    } else {
    	av_log(ctx, AV_LOG_DEBUG, "#################### rotate 1\n");
        av_frame_free(&s->frames[0]);
        av_frame_free(&s->index_frames[0]);
        memmove(&s->frames[0], &s->frames[1], sizeof(*s->frames) * (s->nb_inputs - 1));
        memmove(&s->index_frames[0], &s->index_frames[1], sizeof(*s->index_frames) * (s->nb_inputs - 1));
        s->frames[s->nb_inputs - 1] = av_frame_clone(in);
        s->index_frames[s->nb_inputs - 1] = av_frame_clone(index);
    }
    av_log(ctx, AV_LOG_DEBUG, "#################### get output buffer %d\n",s->nb_frames);
    out = ff_get_video_buffer(outlink, outlink->w, outlink->h);
    if (!out)
        return NULL;
    out->pts = s->frames[0]->pts;
    //out->pts = in->pts;
    td.xoff = (index->width - in->width)/2; 
    td.yoff = (index->height - in->height)/2; 
    td.out = out;
    td.index = s->index_frames;
    td.in = s->frames;
    av_log(ctx, AV_LOG_DEBUG, "#################### calling internal execute\n");
    ctx->internal->execute(ctx, martin_frame, &td, NULL, FFMIN(s->height[1], ff_filter_get_nb_threads(ctx)));
    av_frame_free(&in);
    av_log(ctx, AV_LOG_DEBUG, "#################### returning out\n");
    return out;
}
static int filter_frame_for_dualinput(FFFrameSync *fs)
{   
    AVFilterContext *ctx = fs->parent;
    av_log(ctx, AV_LOG_DEBUG, "#################### dualinput\n");
    AVFrame *default_buf, *index_buf, *dst_buf;
    int ret;
    
    ret = ff_framesync_dualinput_get(fs, &default_buf, &index_buf);
    if (ret < 0)
        return ret;
    if (!index_buf)
        return ff_filter_frame(ctx->outputs[0], default_buf);
    dst_buf = martin_filter_frame(ctx, default_buf, index_buf);
    return ff_filter_frame(ctx->outputs[0], dst_buf);
}   

static av_cold int init(AVFilterContext *ctx)
{
    av_log(ctx, AV_LOG_INFO, "#################### init\n");
    MartinContext *s = ctx->priv;

    //s->nb_inputs = 256;

    s->frames = av_calloc(s->nb_inputs, sizeof(*s->frames));
    s->index_frames = av_calloc(s->nb_inputs, sizeof(*s->index_frames));
    if (!s->frames || !s->index_frames)
        return AVERROR(ENOMEM);
    s->fs.on_event = filter_frame_for_dualinput;
    return 0;
}


#define OFFSET(x) offsetof(MartinContext, x)
#define FLAGS AV_OPT_FLAG_VIDEO_PARAM | AV_OPT_FLAG_FILTERING_PARAM

static const AVOption martin_options[] = {
    { "frames", "set radius", OFFSET(nb_inputs), AV_OPT_TYPE_INT, {.i64=256}, 1, 256, .flags = FLAGS },
    { "factor", "set factor", OFFSET(factor), AV_OPT_TYPE_FLOAT, {.dbl=2}, 0, UINT16_MAX, .flags = FLAGS },
    { "threshold", "set threshold", OFFSET(threshold), AV_OPT_TYPE_FLOAT, {.dbl=10}, 0, UINT16_MAX, .flags = FLAGS },
    { "low", "set low limit for amplification", OFFSET(llimit), AV_OPT_TYPE_INT, {.i64=UINT16_MAX}, 0, UINT16_MAX, .flags = FLAGS },
    { "high", "set high limit for amplification", OFFSET(hlimit), AV_OPT_TYPE_INT, {.i64=UINT16_MAX}, 0, UINT16_MAX, .flags = FLAGS },
    { "planes", "set what planes to filter", OFFSET(planes), AV_OPT_TYPE_FLAGS, {.i64=7},    0, 15,  FLAGS },
    { NULL },
};

static int activate(AVFilterContext *ctx) {
    MartinContext *s = ctx->priv;
    av_log(ctx, AV_LOG_DEBUG, "#################### activate\n");
    return ff_framesync_activate(&s->fs);
}

static const AVFilterPad inputs[] = {
    {
        .name          = "default",
        .type          = AVMEDIA_TYPE_VIDEO,
    },
    {
        .name          = "index",
        .type          = AVMEDIA_TYPE_VIDEO,
    },
   { NULL }
};

static const AVFilterPad outputs[] = {
    {
        .name          = "default",
        .type          = AVMEDIA_TYPE_VIDEO,
        .config_props  = config_output,
    },
    { NULL }
};

FRAMESYNC_DEFINE_CLASS(martin, MartinContext, fs);

AVFilter ff_vf_martin = {
    .name          = "martin",
    .description   = NULL_IF_CONFIG_SMALL("Martin changes between successive video frames."),
    .priv_size     = sizeof(MartinContext),
    .priv_class    = &martin_class,
    .preinit       = martin_framesync_preinit,
    .query_formats = query_formats,
    .activate      = activate,
    .outputs       = outputs,
    .inputs        = inputs,
    .init          = init,
    .uninit        = uninit,
    .flags         = AVFILTER_FLAG_SUPPORT_TIMELINE_INTERNAL | AVFILTER_FLAG_SLICE_THREADS,
};
