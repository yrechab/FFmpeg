/*
 * Copyright (c) 2015 Paul B Mahol
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

#include <float.h>
#include <math.h>

#include "libavcodec/avfft.h"
#include "libavutil/audio_fifo.h"
#include "libavutil/avassert.h"
#include "libavutil/avstring.h"
#include "libavutil/channel_layout.h"
#include "libavutil/intreadwrite.h"
#include "libavutil/opt.h"
#include "libavutil/parseutils.h"
#include "audio.h"
#include "video.h"
#include "avfilter.h"
#include "internal.h"
#include "window_func.h"
#include "libavutil/plot.h"
#include "libavutil/genutil.h"

enum DisplayMode    { LINE, BAR, DOT, NB_MODES };
enum ChannelMode    { COMBINED, SEPARATE, NB_CMODES };
enum FrequencyScale { FS_LINEAR, FS_LOG, FS_RLOG, NB_FSCALES };
enum AmplitudeScale { AS_LINEAR, AS_SQRT, AS_CBRT, AS_LOG, NB_ASCALES };
#define PC 40
#define NC 10

typedef struct GenVideoContext {
    const AVClass *class;
    int w, h;
    int mode;
    int cmode;
    int fft_bits;
    int ascale, fscale;
    int avg;
    int win_func;
    FFTContext *fft;
    FFTComplex **fft_data;
    float **avg_data;
    float *window_func_lut;
    float overlap;
    float minamp;
    int hop_size;
    int nb_channels;
    int nb_freq;
    int win_size;
    float scale;
    char *colors;
    AVAudioFifo *fifo;
    int64_t pts;

    double p[PC];
    double _p[PC];
    double (*nfunc[4])(int,double*);
    double np[4][NC];
    double _np[4][NC];
    const char *n[4];
    int nmod[4]; // modulo n
    double (*cfunc[3])(double,double*);
    double cp[3][NC];
    double _cp[3][NC];
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
    double speed;
    int count;
    int layers;
    double rot;
    double rotn;
    int hsub, vsub;             ///< chroma subsampling
    int planes;                 ///< number of planes
    int is_rgb;
    int alpha;
    int dbg;
    int offset;
    double last[10][10];
    char *f;
    void (*ffunc)(struct GenVideoContext*,int,AVFrame *in);
    int ag;
    uint8_t *lastFrame;
    char *rf[20];
} GenVideoContext;

#define OFFSET(x) offsetof(GenVideoContext, x)
#define FLAGS AV_OPT_FLAG_FILTERING_PARAM|AV_OPT_FLAG_VIDEO_PARAM
#define RE(x, ch) s->fft_data[ch][x].re
#define IM(x, ch) s->fft_data[ch][x].im
#define M(a, b) (sqrt((a) * (a) + (b) * (b)))

typedef struct FFunc {
    const char *name;
    void (*f)(GenVideoContext*, int,AVFrame *in);
} FFunc;



static const AVOption genvideo_options[] = {
    { "size", "set video size", OFFSET(w), AV_OPT_TYPE_IMAGE_SIZE, {.str = "960x540"}, 0, 0, FLAGS },
    { "s",    "set video size", OFFSET(w), AV_OPT_TYPE_IMAGE_SIZE, {.str = "960x540"}, 0, 0, FLAGS },
    { "mode", "set display mode", OFFSET(mode), AV_OPT_TYPE_INT, {.i64=BAR}, 0, NB_MODES-1, FLAGS, "mode" },
        { "line", "show lines",  0, AV_OPT_TYPE_CONST, {.i64=LINE},   0, 0, FLAGS, "mode" },
        { "bar",  "show bars",   0, AV_OPT_TYPE_CONST, {.i64=BAR},    0, 0, FLAGS, "mode" },
        { "dot",  "show dots",   0, AV_OPT_TYPE_CONST, {.i64=DOT},    0, 0, FLAGS, "mode" },
    { "ascale", "set amplitude scale", OFFSET(ascale), AV_OPT_TYPE_INT, {.i64=AS_LOG}, 0, NB_ASCALES-1, FLAGS, "ascale" },
        { "lin",  "linear",      0, AV_OPT_TYPE_CONST, {.i64=AS_LINEAR}, 0, 0, FLAGS, "ascale" },
        { "sqrt", "square root", 0, AV_OPT_TYPE_CONST, {.i64=AS_SQRT},   0, 0, FLAGS, "ascale" },
        { "cbrt", "cubic root",  0, AV_OPT_TYPE_CONST, {.i64=AS_CBRT},   0, 0, FLAGS, "ascale" },
        { "log",  "logarithmic", 0, AV_OPT_TYPE_CONST, {.i64=AS_LOG},    0, 0, FLAGS, "ascale" },
    { "fscale", "set frequency scale", OFFSET(fscale), AV_OPT_TYPE_INT, {.i64=FS_LINEAR}, 0, NB_FSCALES-1, FLAGS, "fscale" },
        { "lin",  "linear",              0, AV_OPT_TYPE_CONST, {.i64=FS_LINEAR}, 0, 0, FLAGS, "fscale" },
        { "log",  "logarithmic",         0, AV_OPT_TYPE_CONST, {.i64=FS_LOG},    0, 0, FLAGS, "fscale" },
        { "rlog", "reverse logarithmic", 0, AV_OPT_TYPE_CONST, {.i64=FS_RLOG},   0, 0, FLAGS, "fscale" },
    { "win_size", "set window size", OFFSET(fft_bits), AV_OPT_TYPE_INT, {.i64=11}, 4, 16, FLAGS, "fft" },
        { "w16",    0, 0, AV_OPT_TYPE_CONST, {.i64=4},  0, 0, FLAGS, "fft" },
        { "w32",    0, 0, AV_OPT_TYPE_CONST, {.i64=5},  0, 0, FLAGS, "fft" },
        { "w64",    0, 0, AV_OPT_TYPE_CONST, {.i64=6},  0, 0, FLAGS, "fft" },
        { "w128",   0, 0, AV_OPT_TYPE_CONST, {.i64=7},  0, 0, FLAGS, "fft" },
        { "w256",   0, 0, AV_OPT_TYPE_CONST, {.i64=8},  0, 0, FLAGS, "fft" },
        { "w512",   0, 0, AV_OPT_TYPE_CONST, {.i64=9},  0, 0, FLAGS, "fft" },
        { "w1024",  0, 0, AV_OPT_TYPE_CONST, {.i64=10}, 0, 0, FLAGS, "fft" },
        { "w2048",  0, 0, AV_OPT_TYPE_CONST, {.i64=11}, 0, 0, FLAGS, "fft" },
        { "w4096",  0, 0, AV_OPT_TYPE_CONST, {.i64=12}, 0, 0, FLAGS, "fft" },
        { "w8192",  0, 0, AV_OPT_TYPE_CONST, {.i64=13}, 0, 0, FLAGS, "fft" },
        { "w16384", 0, 0, AV_OPT_TYPE_CONST, {.i64=14}, 0, 0, FLAGS, "fft" },
        { "w32768", 0, 0, AV_OPT_TYPE_CONST, {.i64=15}, 0, 0, FLAGS, "fft" },
        { "w65536", 0, 0, AV_OPT_TYPE_CONST, {.i64=16}, 0, 0, FLAGS, "fft" },
    { "win_func", "set window function", OFFSET(win_func), AV_OPT_TYPE_INT, {.i64=WFUNC_HANNING}, 0, NB_WFUNC-1, FLAGS, "win_func" },
        { "rect",     "Rectangular",      0, AV_OPT_TYPE_CONST, {.i64=WFUNC_RECT},     0, 0, FLAGS, "win_func" },
        { "bartlett", "Bartlett",         0, AV_OPT_TYPE_CONST, {.i64=WFUNC_BARTLETT}, 0, 0, FLAGS, "win_func" },
        { "hanning",  "Hanning",          0, AV_OPT_TYPE_CONST, {.i64=WFUNC_HANNING},  0, 0, FLAGS, "win_func" },
        { "hamming",  "Hamming",          0, AV_OPT_TYPE_CONST, {.i64=WFUNC_HAMMING},  0, 0, FLAGS, "win_func" },
        { "blackman", "Blackman",         0, AV_OPT_TYPE_CONST, {.i64=WFUNC_BLACKMAN}, 0, 0, FLAGS, "win_func" },
        { "welch",    "Welch",            0, AV_OPT_TYPE_CONST, {.i64=WFUNC_WELCH},    0, 0, FLAGS, "win_func" },
        { "flattop",  "Flat-top",         0, AV_OPT_TYPE_CONST, {.i64=WFUNC_FLATTOP},  0, 0, FLAGS, "win_func" },
        { "bharris",  "Blackman-Harris",  0, AV_OPT_TYPE_CONST, {.i64=WFUNC_BHARRIS},  0, 0, FLAGS, "win_func" },
        { "bnuttall", "Blackman-Nuttall", 0, AV_OPT_TYPE_CONST, {.i64=WFUNC_BNUTTALL}, 0, 0, FLAGS, "win_func" },
        { "bhann",    "Bartlett-Hann",    0, AV_OPT_TYPE_CONST, {.i64=WFUNC_BHANN},    0, 0, FLAGS, "win_func" },
        { "sine",     "Sine",             0, AV_OPT_TYPE_CONST, {.i64=WFUNC_SINE},     0, 0, FLAGS, "win_func" },
        { "nuttall",  "Nuttall",          0, AV_OPT_TYPE_CONST, {.i64=WFUNC_NUTTALL},  0, 0, FLAGS, "win_func" },
        { "lanczos",  "Lanczos",          0, AV_OPT_TYPE_CONST, {.i64=WFUNC_LANCZOS},  0, 0, FLAGS, "win_func" },
        { "gauss",    "Gauss",            0, AV_OPT_TYPE_CONST, {.i64=WFUNC_GAUSS},    0, 0, FLAGS, "win_func" },
        { "tukey",    "Tukey",            0, AV_OPT_TYPE_CONST, {.i64=WFUNC_TUKEY},    0, 0, FLAGS, "win_func" },
        { "dolph",    "Dolph-Chebyshev",  0, AV_OPT_TYPE_CONST, {.i64=WFUNC_DOLPH},    0, 0, FLAGS, "win_func" },
        { "cauchy",   "Cauchy",           0, AV_OPT_TYPE_CONST, {.i64=WFUNC_CAUCHY},   0, 0, FLAGS, "win_func" },
        { "parzen",   "Parzen",           0, AV_OPT_TYPE_CONST, {.i64=WFUNC_PARZEN},   0, 0, FLAGS, "win_func" },
        { "poisson",  "Poisson",          0, AV_OPT_TYPE_CONST, {.i64=WFUNC_POISSON},  0, 0, FLAGS, "win_func" },
    { "overlap",  "set window overlap", OFFSET(overlap), AV_OPT_TYPE_FLOAT, {.dbl=1.}, 0., 1., FLAGS },
    { "averaging", "set time averaging", OFFSET(avg), AV_OPT_TYPE_INT, {.i64=1}, 0, INT32_MAX, FLAGS },
    { "colors", "set channels colors", OFFSET(colors), AV_OPT_TYPE_STRING, {.str = "red|green|blue|yellow|orange|lime|pink|magenta|brown" }, 0, 0, FLAGS },
    { "cmode", "set channel mode", OFFSET(cmode), AV_OPT_TYPE_INT, {.i64=COMBINED}, 0, NB_CMODES-1, FLAGS, "cmode" },
        { "combined", "show all channels in same window",  0, AV_OPT_TYPE_CONST, {.i64=COMBINED}, 0, 0, FLAGS, "cmode" },
        { "separate", "show each channel in own window",   0, AV_OPT_TYPE_CONST, {.i64=SEPARATE}, 0, 0, FLAGS, "cmode" },
    { "minamp",  "set minimum amplitude", OFFSET(minamp), AV_OPT_TYPE_FLOAT, {.dbl=1e-6}, FLT_MIN, 1e-6, FLAGS },

    { "fx","fx",       OFFSET(fx),     AV_OPT_TYPE_DOUBLE, {.dbl = 100}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "fy","fy",       OFFSET(fy),     AV_OPT_TYPE_DOUBLE, {.dbl = 100}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "x","x",       OFFSET(x),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "y","y",       OFFSET(y),     AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "xn","xn",     OFFSET(xn),    AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "yn","yn",     OFFSET(yn),    AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "delta","d",   OFFSET(delta),   AV_OPT_TYPE_DOUBLE, {.dbl =  0.002}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "form","form", OFFSET(form),  AV_OPT_TYPE_INT,    {.i64 =  1}, 0,  20,  FLAGS },
    { "dim","d",     OFFSET(dim),   AV_OPT_TYPE_INT,    {.i64 =  0}, 0,  255,  FLAGS },
    { "ag","ag",     OFFSET(ag),   AV_OPT_TYPE_INT,    {.i64 =  4}, 1,  255,  FLAGS },
    { "length","len",OFFSET(length), AV_OPT_TYPE_DOUBLE, {.dbl =  2*M_PI}, 0,  DBL_MAX,  FLAGS },
    { "start","s",   OFFSET(start), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, 0,  DBL_MAX,  FLAGS },
    { "rot","rot",   OFFSET(rot),   AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "rotn","rotn", OFFSET(rotn),  AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "layers","l",  OFFSET(layers),AV_OPT_TYPE_INT,    {.i64 =  4}, 1,  INT_MAX,  FLAGS },
    { "count","c",   OFFSET(count), AV_OPT_TYPE_INT,    {.i64 =  1}, 1,  INT_MAX,  FLAGS },
    { "speed","s",   OFFSET(speed), AV_OPT_TYPE_DOUBLE, {.dbl =  0.01}, 0,  DBL_MAX,  FLAGS },
  
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


    { "c0",  "c0",   OFFSET(c[0]),    AV_OPT_TYPE_STRING, {.str="sin"}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "c00","c00",   OFFSET(cp[0][0]), AV_OPT_TYPE_DOUBLE, {.dbl =  5}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c01","c01",   OFFSET(cp[0][1]), AV_OPT_TYPE_DOUBLE, {.dbl =  3}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c02","c02",   OFFSET(cp[0][2]), AV_OPT_TYPE_DOUBLE, {.dbl =  2}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c03","c03",   OFFSET(cp[0][3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c04","c04",   OFFSET(cp[0][4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c05","c05",   OFFSET(cp[0][5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c06","c06",   OFFSET(cp[0][6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c07","c07",   OFFSET(cp[0][7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c08","c08",   OFFSET(cp[0][8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c09","c09",   OFFSET(cp[0][9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c1",  "c1",   OFFSET(c[1]),    AV_OPT_TYPE_STRING, {.str="sin"}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "c10","c10",   OFFSET(cp[1][0]), AV_OPT_TYPE_DOUBLE, {.dbl =  2}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c11","c11",   OFFSET(cp[1][1]), AV_OPT_TYPE_DOUBLE, {.dbl =0.1}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c12","c12",   OFFSET(cp[1][2]), AV_OPT_TYPE_DOUBLE, {.dbl =0.8}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c13","c13",   OFFSET(cp[1][3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c14","c14",   OFFSET(cp[1][4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c15","c15",   OFFSET(cp[1][5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c16","c16",   OFFSET(cp[1][6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c17","c17",   OFFSET(cp[1][7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c18","c18",   OFFSET(cp[1][8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c19","c19",   OFFSET(cp[1][9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c2",  "c2",   OFFSET(c[2]),    AV_OPT_TYPE_STRING, {.str="sin"}, CHAR_MIN, CHAR_MAX, FLAGS },
    { "c20","c20",   OFFSET(cp[2][0]), AV_OPT_TYPE_DOUBLE, {.dbl =  4}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c21","c21",   OFFSET(cp[2][1]), AV_OPT_TYPE_DOUBLE, {.dbl =0.5}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c22","c22",   OFFSET(cp[2][2]), AV_OPT_TYPE_DOUBLE, {.dbl =0.5}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c23","c23",   OFFSET(cp[2][3]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c24","c24",   OFFSET(cp[2][4]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c25","c25",   OFFSET(cp[2][5]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c26","c26",   OFFSET(cp[2][6]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c27","c27",   OFFSET(cp[2][7]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c28","c28",   OFFSET(cp[2][8]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },
    { "c29","c29",   OFFSET(cp[2][9]), AV_OPT_TYPE_DOUBLE, {.dbl =  0}, -DBL_MAX,  DBL_MAX,  FLAGS },

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


    { "cmod","cmod",OFFSET(cmod), AV_OPT_TYPE_INT,    {.i64=2},    0,        2,        FLAGS },
    { "rgb","rgb",OFFSET(is_rgb), AV_OPT_TYPE_INT,    {.i64=0},    0,        1,        FLAGS },
    { "alpha","alpha",OFFSET(alpha), AV_OPT_TYPE_INT,    {.i64=0},    0,        1,        FLAGS },
    { "dbg","dbg",OFFSET(dbg), AV_OPT_TYPE_INT,    {.i64=0},    0,        1,        FLAGS },
    { "ofs","ofs",OFFSET(offset), AV_OPT_TYPE_INT,    {.i64=0},    0,     INT_MAX,        FLAGS },

    { NULL }
};

AVFILTER_DEFINE_CLASS(genvideo);

static int query_formats(AVFilterContext *ctx)
{
    GenVideoContext *genvideo = ctx->priv;
    AVFilterFormats *formats = NULL;
    AVFilterChannelLayouts *layouts = NULL;
    AVFilterLink *inlink = ctx->inputs[0];
    AVFilterLink *outlink = ctx->outputs[0];
    static const enum AVSampleFormat sample_fmts[] = { AV_SAMPLE_FMT_FLTP, AV_SAMPLE_FMT_NONE };
    static const enum AVPixelFormat pix_fmts[] = { AV_PIX_FMT_YUV444P, AV_PIX_FMT_NONE };
    static const enum AVPixelFormat pix_fmts_a[] = { AV_PIX_FMT_YUVA444P, AV_PIX_FMT_NONE };
    int ret;

    /* set input audio formats */
    formats = ff_make_format_list(sample_fmts);
    if ((ret = ff_formats_ref(formats, &inlink->out_formats)) < 0)
        return ret;

    layouts = ff_all_channel_layouts();
    if ((ret = ff_channel_layouts_ref(layouts, &inlink->out_channel_layouts)) < 0)
        return ret;

    formats = ff_all_samplerates();
    if ((ret = ff_formats_ref(formats, &inlink->out_samplerates)) < 0)
        return ret;

    /* set output video format */
    formats = ff_make_format_list(genvideo->alpha?pix_fmts_a:pix_fmts);
    if ((ret = ff_formats_ref(formats, &outlink->in_formats)) < 0)
        return ret;

    return 0;
}


static void draw_number(int x, int y, double z, AVFrame *in) {
    av_plot_number(x,y,z,in->width,in->height,in->linesize[0],in->data[0]);
}

static void drawPoint(GenVideoContext *genvideo, AVFrame *in, int plane, double t, PFUNC(f), double a,double phi,int ofsX, int ofsY) {
    double complex _z = f(t,genvideo->_p);
    double complex z = av_genutil_rotate(_z,phi);
    int _x = floor(genvideo->fx*creal(z) + genvideo->w/2) + genvideo->x + ofsX;
    int _y = floor(genvideo->fy*cimag(z) + genvideo->h/2) + genvideo->y + ofsY;
    int x,y,w,h;
    if((plane==1 || plane==2)&&(genvideo->is_rgb==0)) {
        w = AV_CEIL_RSHIFT(genvideo->w, genvideo->hsub);
        x = AV_CEIL_RSHIFT(_x, genvideo->hsub);
        h = AV_CEIL_RSHIFT(genvideo->h, genvideo->vsub);
        y = AV_CEIL_RSHIFT(_y, genvideo->vsub);
    } else {
        x=_x;y=_y;w=genvideo->w;h=genvideo->h;
    }
    av_plot_form(genvideo->form,x,y,a,w,h,in->linesize[plane],in->data[plane]);
}

static void blackYUV(GenVideoContext *genvideo, AVFrame *in) {
    int plane,x,y;
    const int width = AV_CEIL_RSHIFT(in->width, genvideo->hsub);
    const int height = AV_CEIL_RSHIFT(in->height, genvideo->vsub);
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
    if(genvideo->alpha && genvideo->dim==0) {
        if(in->data[plane]) {
            for (y = 0; y < height; y++) {
                ptr = in->data[plane] + in->linesize[plane] * y;
                for (x = 0; x < width; x++) {
                    ptr[x] = 255;
                }
            }
        }
     }
}


static void curve(GenVideoContext *genvideo, PFUNC(f),int n, AVFrame *in) {
    blackYUV(genvideo,in);
    double t = genvideo->start;
    int k;
    if(genvideo->dbg) {
        for(k=0;k<4;k++) {
            draw_number(genvideo->w-k*100-100, genvideo->h-35, genvideo->_p[k], in);
        }
        draw_number(100, genvideo->h-35, n, in);
    }
    while(t<genvideo->start+genvideo->length) {
        double colors[3];
        av_genutil_get_color(genvideo->cfunc, genvideo->_cp,t,genvideo->cmod, genvideo->is_rgb, colors);
        drawPoint(genvideo,in,0,t,f,colors[0],0,0,0);
        drawPoint(genvideo,in,1,t,f,colors[1],0,0,0);
        drawPoint(genvideo,in,2,t,f,colors[2],0,0,0);
        t+=genvideo->delta;
    }
}




static void lissG(GenVideoContext *genvideo, int n, AVFrame *in) {
    curve(genvideo,av_genutil_lissajousG,n,in);
}

static void lissQ(GenVideoContext *genvideo, int n, AVFrame *in) {
    curve(genvideo,av_genutil_lissajousQ,n,in);
}

static void liss(GenVideoContext *genvideo, int n, AVFrame *in) {
    curve(genvideo,av_genutil_lissajous,n,in);
}

static void leg(GenVideoContext *genvideo, int n, AVFrame *in) {
    curve(genvideo,av_genutil_legendre,n,in);
}

static void hypo(GenVideoContext *genvideo, int n, AVFrame *in) {
    curve(genvideo,av_genutil_hypocycloid,n,in);
}

static void cbP(GenVideoContext *genvideo, int n, AVFrame *in) {
    curve(genvideo,av_genutil_cb,n,in);
}

static void sinoid(GenVideoContext *genvideo, int n, AVFrame *in) {
    curve(genvideo,av_genutil_sinusoid,n,in);
}

static void capri(GenVideoContext *genvideo, int n, AVFrame *in) {
    curve(genvideo,av_genutil_capricornoid,n,in);
}

static void scara(GenVideoContext *genvideo, int n, AVFrame *in) {
    curve(genvideo,av_genutil_scarabaeus,n,in);
}

static void tro(GenVideoContext *genvideo, int n, AVFrame *in) {
    curve(genvideo,av_genutil_trochoid,n,in);
}

static void circ(GenVideoContext *genvideo, int n, AVFrame *in) {
    curve(genvideo,av_genutil_circoid,n,in);
}

static void super_rose(GenVideoContext *genvideo, int n, AVFrame *in) {
    curve(genvideo,av_genutil_super_rose,n,in);
}

static void zero(GenVideoContext *genvideo, int n, AVFrame *in) {
}


static void circs(GenVideoContext *genvideo, int n, AVFrame *in) {
    double count = genvideo->count + (int)genvideo->_p[3];
    double speed = genvideo->speed;
    blackYUV(genvideo,in);
	    if(genvideo->dbg) {
        int k;
        for(k=0;k<4;k++) {
            draw_number(genvideo->w-k*100-100, genvideo->h-35, genvideo->_p[k], in);
        }
    }

    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int ofsX = floor(cos(alpha) * genvideo->_p[0]);
        int ofsY = floor(sin(alpha) * genvideo->_p[1]);
        double start = 2*M_PI*k/count;
        int j;
        //int layers = genvideo->layers;
        int layers = genvideo->layers + (int)genvideo->_p[2];
        for(j=0;j<layers;j++) {
            double s = start + genvideo->length*j/layers;
            double t = s + n*speed;
            drawPoint(genvideo,in,genvideo->alpha?3:0,t,av_genutil_circ,genvideo->alpha?0:1,0,ofsX,ofsY);
        }
    }
} 

static void plot(GenVideoContext *genvideo, PFUNC(f),int n, AVFrame *in) {
    blackYUV(genvideo,in);
    if(genvideo->dbg) {
        int k;
        for(k=0;k<2;k++) {
            draw_number(genvideo->w-k*100-100, genvideo->h-35, genvideo->_p[k], in);
        }
        draw_number(genvideo->w-2*100-100, genvideo->h-35, genvideo->_p[38], in);
        draw_number(genvideo->w-3*100-100, genvideo->h-35, genvideo->_p[39], in);
    }
    double count = genvideo->count + (int)genvideo->_p[39];
    double speed = genvideo->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int j;
        double start = 0;
        int layers = genvideo->layers + (int)genvideo->_p[38];
        for(j=0;j<layers;j++) {
            double s = start + genvideo->length*j/layers;
            double t = s + n*speed;
            drawPoint(genvideo,in,genvideo->alpha?3:0,t,av_genutil_circ,genvideo->alpha?0:1,alpha,0,0);
        }
    }
}
static void lemG(GenVideoContext *genvideo, int n, AVFrame *in) {
    plot(genvideo,av_genutil_lemniskateG,n,in);
}

static void lemB(GenVideoContext *genvideo, int n, AVFrame *in) {
    plot(genvideo,av_genutil_lemniskateB,n,in);
}

static void trop(GenVideoContext *genvideo, int n, AVFrame *in) {
    plot(genvideo,av_genutil_trochoid,n,in);
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
    {"circs",circs},
    {"scara",scara},
    {"tro",tro},
    {"trop",trop},
    {"lemG",lemG},
    {"lemB",lemB},
    {"superrose",super_rose},
    {NULL,NULL}
};

static void (*get_func(const char *name))(GenVideoContext*, int, AVFrame*) {
    int k=0;
    while(ffuncs[k].name) {
        if(!strcmp(name, ffuncs[k].name)) {
            return ffuncs[k].f;
        }
        k++;
    }
    return NULL;
}



static av_cold void log_parameters(AVFilterContext *ctx,double* p, int len) {
    int k;
    for(k=0;k<len;k++) {
        av_log(ctx, AV_LOG_INFO," %f",p[k]); 
    }
    av_log(ctx, AV_LOG_INFO,"\n"); 
}

static void initLast(GenVideoContext *genvideo) {
    if(genvideo->dim <= 0) return;
    int x, y;
    uint8_t *last;
    int height = genvideo->h;
    int width = genvideo->w;
    for (y = 0; y < height; y++) {
        last = genvideo->lastFrame + width * y;
        for (x = 0; x < width; x++) {
             last[x] = 255;
        }
    }
}

static av_cold int init(AVFilterContext *ctx)
{
    GenVideoContext *s = ctx->priv;
    int k;

    av_log(ctx, AV_LOG_INFO, "rgb=%d ofs=%d x=%f y=%f xn=%f yn=%f, delta=%f start=%f length=%f fx=%f fy=%f dbg=%d\n", 
            s->is_rgb, s->offset, s->x, s->y, s->xn, s->yn, 
            s->delta, s->start, s->length, s->fx,s->fy,s->dbg);

    for(k=0;k<4;k++) {
        if(s->n[k]) {
            s->nfunc[k] = av_genutil_get_nfunc(s->n[k]);
            if(!s->nfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for n%d not found %s\n", k, s->n[k]);
            } else {
                av_log(ctx, AV_LOG_INFO, "function for n%d is %s", k, s->n[k]);
                log_parameters(ctx,s->np[k],10);
            }
        } else {
            av_log(ctx, AV_LOG_WARNING, "no function given for n%d\n",k);
        }
    }
    

    for(k=0;k<3;k++) {
        if(s->c[k]) {
            s->cfunc[k] = av_genutil_get_cfunc(s->c[k]);
            if(!s->cfunc[k]) {
                av_log(ctx, AV_LOG_WARNING, "function for cf%d not found %s\n", k, s->c[k]);
                s->cfunc[k] = czero;
            } else {
                av_log(ctx, AV_LOG_INFO, "function for cf%d is %s", k, s->c[k]);
                log_parameters(ctx,s->cp[k],10);
            }
        } else {
            av_log(ctx, AV_LOG_WARNING, "no function given for c\n");
            s->cfunc[k] = czero;
        }
    }

    if(s->f) {
        s->ffunc = get_func(s->f);
        if(!s->ffunc) {
            av_log(ctx, AV_LOG_WARNING, "function for f not found %s\n", s->f);
            s->ffunc = zero;
        } else {
            av_log(ctx, AV_LOG_INFO, "function for f is %s\n", s->f);
            log_parameters(ctx,s->p,PC);
        }
    } else {
        av_log(ctx, AV_LOG_WARNING, "no function given for f\n");
        s->ffunc = zero;
    }

    s->lastFrame = av_calloc(1920*1080,sizeof(uint8_t));
    initLast(s);
    s->pts = AV_NOPTS_VALUE;
    /*
    for(k=0;k<10;k++) {
        if(s->rf[k]) {
            av_genutil_replace_char(s->rf[k],'|',' ');
        }
    }
    */
    return 0;
}

static int config_output(AVFilterLink *outlink)
{
    AVFilterContext *ctx = outlink->src;
    AVFilterLink *inlink = ctx->inputs[0];
    GenVideoContext *s = ctx->priv;
    float overlap;
    int i;

    s->nb_freq = 1 << (s->fft_bits - 1);
    s->win_size = s->nb_freq << 1;
    av_audio_fifo_free(s->fifo);
    av_fft_end(s->fft);
    s->fft = av_fft_init(s->fft_bits, 0);
    if (!s->fft) {
        av_log(ctx, AV_LOG_ERROR, "Unable to create FFT context. "
               "The window size might be too high.\n");
        return AVERROR(ENOMEM);
    }

    /* FFT buffers: x2 for each (display) channel buffer.
     * Note: we use free and malloc instead of a realloc-like function to
     * make sure the buffer is aligned in memory for the FFT functions. */
    for (i = 0; i < s->nb_channels; i++) {
        av_freep(&s->fft_data[i]);
        av_freep(&s->avg_data[i]);
    }
    av_freep(&s->fft_data);
    av_freep(&s->avg_data);
    s->nb_channels = inlink->channels;

    s->fft_data = av_calloc(s->nb_channels, sizeof(*s->fft_data));
    if (!s->fft_data)
        return AVERROR(ENOMEM);
    s->avg_data = av_calloc(s->nb_channels, sizeof(*s->avg_data));
    if (!s->fft_data)
        return AVERROR(ENOMEM);
    for (i = 0; i < s->nb_channels; i++) {
        s->fft_data[i] = av_calloc(s->win_size, sizeof(**s->fft_data));
        s->avg_data[i] = av_calloc(s->nb_freq, sizeof(**s->avg_data));
        if (!s->fft_data[i] || !s->avg_data[i])
            return AVERROR(ENOMEM);
    }

    /* pre-calc windowing function */
    s->window_func_lut = av_realloc_f(s->window_func_lut, s->win_size,
                                      sizeof(*s->window_func_lut));
    if (!s->window_func_lut)
        return AVERROR(ENOMEM);
    generate_window_func(s->window_func_lut, s->win_size, s->win_func, &overlap);
    if (s->overlap == 1.)
        s->overlap = overlap;
    s->hop_size = (1. - s->overlap) * s->win_size;
    if (s->hop_size < 1) {
        av_log(ctx, AV_LOG_ERROR, "overlap %f too big\n", s->overlap);
        return AVERROR(EINVAL);
    }

    for (s->scale = 0, i = 0; i < s->win_size; i++) {
        s->scale += s->window_func_lut[i] * s->window_func_lut[i];
    }

    outlink->frame_rate = av_make_q(inlink->sample_rate, s->win_size * (1.-s->overlap));
    av_log(ctx, AV_LOG_INFO, "overlap is %f\n ",   s->overlap );
    av_log(ctx, AV_LOG_INFO, "frame rate is %f\n ",   av_q2d(outlink->frame_rate) );

    outlink->sample_aspect_ratio = (AVRational){1,1};
    outlink->w = s->w;
    outlink->h = s->h;

    s->fifo = av_audio_fifo_alloc(inlink->format, inlink->channels, s->win_size);
    if (!s->fifo)
        return AVERROR(ENOMEM);



    s->hsub = 0;
    s->vsub = 0;
    s->planes = s->alpha?4:3;

    

    return 0;
}


static double freq_avg(double *buf, int len, double val) {
    int k;
    for(k=0;k<len;k++) {
        buf[k] = buf[k+1];
    }
    buf[len]=val;
    return av_genutil_avg(buf,len+1);
}

static void dim(GenVideoContext *genvideo, AVFrame *out) {
    if(genvideo->dim <= 0) return;
    int plane = genvideo->alpha?3:0;
    uint8_t *src;
    uint8_t *dst;
    int x,y;
    int height = genvideo->h;
    int width = genvideo->w;
    for (y = 0; y < height; y++) {
        src = genvideo->lastFrame + out->linesize[plane] * y;
        dst = out->data[plane] + out->linesize[plane] * y;
        for (x = 0; x < width; x++) {
            if(genvideo->alpha) {
                dst[x] = src[x]>(255-genvideo->dim)?255:src[x]+genvideo->dim;
            } else {
                dst[x] = src[x]<genvideo->dim?0:src[x]-genvideo->dim;
            }
        }
    }
}

static void copy0(GenVideoContext *genvideo, AVFrame *out) {
    if(genvideo->dim <= 0) return;
    int x, y;
    uint8_t *ptr;
    uint8_t *last;
    int plane = genvideo->alpha?3:0;
    int height = genvideo->h;
    int width = genvideo->w;
    for (y = 0; y < height; y++) {
        ptr = out->data[plane] + out->linesize[plane] * y;
        last = genvideo->lastFrame + out->linesize[plane] * y;
        for (x = 0; x < width; x++) {
             last[x] = ptr[x];
        }
    }
 
}

static int modify(GenVideoContext *s, int frameNumber) {
    int j,i,k;
    const char *format = "%s %d %s %s %d %lf %d %d %d";
    const int win_size = s->win_size;
    double *buf=calloc(win_size,sizeof(double));
    i=0;
    for(j=0;j<10;j++) {
        if(s->rf[j]) {
            char input[40];
            char target[4];
            char src[4];
            char mode[4];
            int index,m,a,b,ag;
            double amount;
            strcpy(input, s->rf[j]);
            sscanf(input, format, target, &index, mode, src, &m, &amount, &a, &b, &ag);
            double value;
            switch(src[0]) {
                case 'a': {
                    int len = b-a+1;
                    for(k=0;k<len;k++) {
                        buf[k] =  av_clipd(M(RE(k+a, 0), IM(k+a, 0)) / s->scale, 0, 1);
                    }
                    value = freq_avg(s->last[i],ag,amount *  av_genutil_avg(buf,len));
                    break;
                }
                case 'n': {
                    value = s->nfunc[m](s->nmod[m]?frameNumber%s->nmod[m]:frameNumber,s->_np[m]);
                    break;
                }                             
            }
            switch(target[0]) {
                case 'p': {
                    if(mode[0] == 'o') s->_p[index] = value;
                    if(mode[0] == 'a') s->_p[index] = s->_p[index] + value;
                    if(mode[0] == 's') s->_p[index] = s->_p[index] - value;
                    break;
                }
                case 'c': {
                    int f = index/10;
                    int i = index % 10;
                    if(mode[0] == 'o') s->_cp[f][i] = value;
                    if(mode[0] == 'a') s->_cp[f][i] = s->cp[f][i] + value;
                    if(mode[0] == 's') s->_cp[f][i] = s->cp[f][i] - value;
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
    free(buf);
    return 0;
}

static void copy_params(GenVideoContext *genvideo) {
    int k,j;
    for(k=0;k<PC;k++) {
        genvideo->_p[k] = genvideo->p[k];
    }
    for(k=0;k<3;k++) {
        for(j=0;j<NC;j++) {
            genvideo->_cp[k][j] = genvideo->cp[k][j];
        }
    }
    for(k=0;k<4;k++) {
        for(j=0;j<NC;j++) {
            genvideo->_np[k][j] = genvideo->np[k][j];
        }
    }
}
static int plot_freqs(AVFilterLink *inlink, AVFrame *in)
{
    AVFilterContext *ctx = inlink->dst;
    AVFilterLink *outlink = ctx->outputs[0];
    GenVideoContext *s = ctx->priv;
    const int win_size = s->win_size;
    AVFrame *out;
    int ch, n;

    out = ff_get_video_buffer(outlink, outlink->w, outlink->h);
    if (!out)
        return AVERROR(ENOMEM);

    for (n = 0; n < outlink->h; n++)
        memset(out->data[0] + out->linesize[0] * n, 0, outlink->w * 4);
    dim(s,out);
    /* fill FFT input with the number of samples available */
    for (ch = 0; ch < s->nb_channels; ch++) {
        const float *p = (float *)in->extended_data[ch];

        for (n = 0; n < in->nb_samples; n++) {
            s->fft_data[ch][n].re = p[n] * s->window_func_lut[n];
            s->fft_data[ch][n].im = 0;
        }
        for (; n < win_size; n++) {
            s->fft_data[ch][n].re = 0;
            s->fft_data[ch][n].im = 0;
        }
    }

    /* run FFT on each samples set */
    for (ch = 0; ch < s->nb_channels; ch++) {
        av_fft_permute(s->fft, s->fft_data[ch]);
        av_fft_calc(s->fft, s->fft_data[ch]);
    }
    copy_params(s);
    modify(s,inlink->frame_count_out);
    s->ffunc(s,inlink->frame_count_out,out);
    copy0(s,out);
    out->pts = in->pts;
    out->sample_aspect_ratio = (AVRational){1,1};
    return ff_filter_frame(outlink, out);
}

static int filter_frame(AVFilterLink *inlink, AVFrame *in)
{
    AVFilterContext *ctx = inlink->dst;
    GenVideoContext *s = ctx->priv;
    AVFrame *fin = NULL;
    int consumed = 0;
    int ret = 0;

    if (s->pts == AV_NOPTS_VALUE)
        s->pts = in->pts - av_audio_fifo_size(s->fifo);

    av_audio_fifo_write(s->fifo, (void **)in->extended_data, in->nb_samples);
    while (av_audio_fifo_size(s->fifo) >= s->win_size) {
        fin = ff_get_audio_buffer(inlink, s->win_size);
        if (!fin) {
            ret = AVERROR(ENOMEM);
            goto fail;
        }

        fin->pts = s->pts + consumed;
        consumed += s->hop_size;
        ret = av_audio_fifo_peek(s->fifo, (void **)fin->extended_data, s->win_size);
        if (ret < 0)
            goto fail;

        ret = plot_freqs(inlink, fin);
        av_frame_free(&fin);
        av_audio_fifo_drain(s->fifo, s->hop_size);
        if (ret < 0)
            goto fail;
    }

fail:
    s->pts = AV_NOPTS_VALUE;
    av_frame_free(&fin);
    av_frame_free(&in);
    return ret;
}

static av_cold void uninit(AVFilterContext *ctx)
{
    GenVideoContext *s = ctx->priv;
    int i;

    av_fft_end(s->fft);
    for (i = 0; i < s->nb_channels; i++) {
        if (s->fft_data)
            av_freep(&s->fft_data[i]);
        if (s->avg_data)
            av_freep(&s->avg_data[i]);
    }
    av_freep(&s->fft_data);
    av_freep(&s->avg_data);
    av_freep(&s->window_func_lut);
    av_audio_fifo_free(s->fifo);
    av_freep(&s->lastFrame);
}

static const AVFilterPad genvideo_inputs[] = {
    {
        .name         = "default",
        .type         = AVMEDIA_TYPE_AUDIO,
        .filter_frame = filter_frame,
    },
    { NULL }
};

static const AVFilterPad genvideo_outputs[] = {
    {
        .name          = "default",
        .type          = AVMEDIA_TYPE_VIDEO,
        .config_props  = config_output,
    },
    { NULL }
};

AVFilter ff_avf_genvideo = {
    .name          = "genvideo",
    .description   = NULL_IF_CONFIG_SMALL("Convert input audio to a frequencies video output."),
    .init          = init,
    .uninit        = uninit,
    .query_formats = query_formats,
    .priv_size     = sizeof(GenVideoContext),
    .inputs        = genvideo_inputs,
    .outputs       = genvideo_outputs,
    .priv_class    = &genvideo_class,
};
