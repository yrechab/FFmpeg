#ifndef AVUTIL_GENUTIL_H
#define AVUTIL_GENUTIL_H

#include <float.h>
#include <stdio.h>
#include <stdint.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include "libavutil/frame.h"

#define PFUNC(x) double complex (*x)(double t, double *p)

#define NMAXPARAMS 40
#define CMAXPARAMS 40

#define SIGN(x) (x==0?0:x/fabsl(x))

enum { Y, U, V, A };

typedef struct GenutilFuncParams {
    // func params
    double *p;
    //color
    double (*cfunc[3])(double,double*);
    double cp[3][10];
    int cmod;
    // curve params
    double delta;
    double start;
    double length;
    // frame
    int w,h;
    int alpha;
    int is_packed_rgb;
    int is_rgb;
    uint8_t rgba_map[4];

    // zoom/transpose
    double fx;
    double fy;
    double x;
    double y;
    double rot;
 
    // plot parms
    int dim;
    double speed;
    int count;
    int layers;
    int form;
} GenutilFuncParams;

double complex av_genutil_rotate(double complex z, double phi);
double complex av_genutil_circ(double t,  double *p);
double complex av_genutil_kardio(double t,  double *p);
double complex av_genutil_lemniskateG(double t,  double *p);
double complex av_genutil_lemniskateB(double t,  double *p);
double complex av_genutil_circoid(double t,  double *p);
double complex av_genutil_sinusoid(double t,  double *p);
double complex av_genutil_lissajous(double t, double *p);
double complex av_genutil_lissajousG(double t, double *p);
double complex av_genutil_lissajousQ(double t, double *p);
double complex av_genutil_trochoid(double t, double *p);
double complex av_genutil_legendre(double t, double *p);
double complex av_genutil_hypocycloid(double t, double *p);
double complex av_genutil_cb(double t, double *p);
double complex av_genutil_capricornoid(double t, double *p);
double complex av_genutil_scarabaeus(double t, double *p);
double complex av_genutil_epicycloid(double t, double *p);
double complex av_genutil_nodal(double t, double *p);
double complex av_genutil_talbot(double t, double *p);
double complex av_genutil_folium(double t, double *p);
double complex av_genutil_rhodonea(double t, double *p);
double complex av_genutil_super_spiral(double t, double *p);
double complex av_genutil_super_rose(double t, double *p);
double complex av_genutil_gielis(double t, double *p);
double complex av_genutil_epi_spiral(double t, double *p);
double complex av_genutil_spiral(double t, double *p);
double complex av_genutil_atom_spiral(double t, double *p);
double complex av_genutil_cotes_spiral(double t, double *p);
double complex av_genutil_sin_spiral(double t, double *p);
double complex av_genutil_maclaurin(double t, double *p);
double complex av_genutil_harmonic(double t, double *p);
double complex av_genutil_cornoid(double t, double *p);
double complex av_genutil_line(double t, double *p);
double complex av_genutil_poly(double t, double *p);
double complex av_genutil_addoid(double t, double *p);

double (*av_genutil_get_nfunc(const char *name))(int,double*);
void av_genutil_parse_nfunc(const char *nf, double *p, double (**f)(int,double*));
double (*av_genutil_get_cfunc(const char *name))(double,double*);
void av_genutil_parse_cfunc(const char *nf, double *p, double (**f)(double,double*));
void (*av_genutil_get_ffunc(const char *name))(GenutilFuncParams*,int,AVFrame*);
void av_genutil_parse_ffunc(const char *nf, double *p, void (**f)(GenutilFuncParams*,int,AVFrame*));

void av_genutil_get_color(double (*cfunc[3])(double,double*), double cp[3][10], double t, int cmod, int is_rgb, double *out);

double av_genutil_avg(double *p, int len);
void av_genutil_replace_char(char *str, char needle, char rep);
double czero(double t, double *p);
void gzero(GenutilFuncParams *params, int n, AVFrame *in);
void debug(GenutilFuncParams *params, int frame_number, AVFrame *in, int plane);
void av_genutil_draw_number(int x, int y, double z, AVFrame *in, int plane);
void av_genutil_draw_int(int x, int y, int n, AVFrame *in,int plane);
#endif
