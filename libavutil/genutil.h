#ifndef AVUTIL_GENUTIL_H
#define AVUTIL_GENUTIL_H

#include <float.h>
#include <stdio.h>
#include <stdint.h>
#include <complex.h>
#include <math.h>
#include <string.h>

#define PFUNC(x) double complex (*x)(double t, double *p)

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


double (*av_genutil_get_nfunc(const char *name))(int,double*);
double (*av_genutil_get_cfunc(const char *name))(double,double*);

void av_genutil_get_color(double (*cfunc[3])(double,double*), double cp[3][10], double t, int cmod, int is_rgb, double *out);

double av_genutil_avg(double *p, int len);

double czero(double t, double *p);

#endif
