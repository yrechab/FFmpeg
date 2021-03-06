#include "genutil.h"
#include "plot.h"
#include <string.h>
#include "libavutil/avstring.h"


typedef struct CFunc {
    const char *name;
    double (*f)(double t, double *p);
} CFunc;

#define P2C(r,t) r*cos(t)+I*r*sin(t)

double complex av_genutil_rotate(double complex z, double phi) {
    return phi!=0.0?(cos(phi) + sin(phi)*I) * z:z;
}

int av_genutil_spline(double *p, double *out) {
    /** Step 0 */
    int n, i, j, end;
    n = p[0];
    n--;
    double *x = calloc(n+1, sizeof(double));
    double *a = calloc(n+1, sizeof(double));
    double *h = calloc(n, sizeof(double));
    double *A = calloc(n, sizeof(double));
    double *l = calloc(n+1, sizeof(double));
    double *u = calloc(n+1, sizeof(double));
    double *z = calloc(n+1, sizeof(double));
    double *c = calloc(n+1, sizeof(double));
    double *b = calloc(n, sizeof(double));
    double *d = calloc(n, sizeof(double));
    for (i = 0; i < n + 1; ++i) {
        x[i] = p[i*2+1];
    }
    end = x[n];
    //scanf("%lf", &x[i]);
    for (i = 0; i < n + 1; ++i) {
        a[i] = p[i*2+2];
    }
    //scanf("%lf", &a[i]);

    /** Step 1 */
    for (i = 0; i <= n - 1; ++i) {
        h[i] = x[i + 1] - x[i];
        if(h[i] == 0) h[i] = 0.0001;
    }

    /** Step 2 */
    for (i = 1; i <= n - 1; ++i)
        A[i] = 3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1];

    /** Step 3 */
    l[0] = 1;
    u[0] = 0;
    z[0] = 0;

    /** Step 4 */
    for (i = 1; i <= n - 1; ++i) {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * u[i - 1];
        u[i] = h[i] / l[i];
        z[i] = (A[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    /** Step 5 */
    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    /** Step 6 */
    for (j = n - 1; j >= 0; --j) {
        c[j] = z[j] - u[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }
    
    /** Step 7 */
    //printf("%2s %8s %8s %8s %8s\n", "i", "ai", "bi", "ci", "di");
    out[0] = n;
    out[1] = end;
    for (i = 0; i < n; ++i) {
        out[i*5+2] = x[i];
        out[i*5+3] = a[i];
        out[i*5+4] = b[i];
        out[i*5+5] = c[i];
        out[i*5+6] = d[i];
    }
    //    printf("%2d %.16f %2.16f %2.16f %2.16f\n", i, a[i], b[i], c[i], d[i]);
    free(x);
    free(a);
    free(h);
    free(A);
    free(l);
    free(u);
    free(z);
    free(c);
    free(b);
    free(d);
    return 0;
}


static double complex av_genutil_line_p(double t, double x0, double y0, double x1, double y1,double *next) {
    int steep = fabs(y1 - y0) > fabs(x1 - x0);
    double tmp;
    if (steep) {
        tmp = x0;
        x0 = y0;
        y0 = tmp;
        tmp = x1;
        x1 = y1;
        y1 = tmp;
    }
    double dx = fabs(x1 - x0);
    double dy = y1 - y0;
    double gradient = dy / dx;
    double x,y;

    if(x0>x1) {
        if(x0-t<=x1) {
            y = y0 + fmod(t,dx) * gradient;
            x = x0 - fmod(t,dx);;
            *next = x0-x1;
        } else {
            y = y0 + t * gradient;
            x = x0 - t;
            *next = 0;
        }
     } else {
        if(x0+t>=x1) {
            y = y0 + fmod(t,dx) * gradient;
            x = x0 + fmod(t,dx);
            *next = x1-x0;
        } else {
            y = y0 + t * gradient;
            x = x0 + t;
            *next = 0;
        }
    }
    if(steep) {
        return y + I * x;
    } else {
        return x + I * y;
    }
}

static double line_length(double x0,double y0,double x1,double y1) {
    if(fabs(y1 - y0) > fabs(x1 - x0)) return fabs(y1-y0);
    return fabs(x1-x0);
}

static double poly_length(double *p) {
    int k;
    int len = floor(p[0]);
    double ret = 0;
    for(k=0;k<len*2;k+=2) {
        ret += line_length(p[k+1],p[k+2],p[k+3],p[k+4]);
    }
    return ret;
}

double complex av_genutil_poly(double t, double *p) {
    int len = floor(p[0]);
    int k;
    double next;
    double complex ret;
    double length = poly_length(p);
    double _t = fmod(t,length);
    //double _t = t;
    for(k=0;k<len*2;k+=2) {
        ret = av_genutil_line_p(_t,p[k+1],p[k+2],p[k+3],p[k+4],&next);
        if(next == 0) break;
        _t-=next;
    }
    return ret;
}

static double complex av_genutil_poly_s(double t, double *p) {
    int len = floor(p[0]);
    int k;
    double next;
    double complex ret = 0;
    double length = poly_length(p);
    if(t>length) return ret;
    for(k=0;k<len*2;k+=2) {
        ret = av_genutil_line_p(t,p[k+1],p[k+2],p[k+3],p[k+4],&next);
        if(next == 0) break;
        t-=next;
    }
    return ret;
}

double complex av_genutil_poly_u(double t, double *p) {
    int len = floor(p[0]);
    int k;
    double next;
    double complex ret;
    double length = poly_length(p);
    double _t = fmod(t,length);
    //double _t = t;
    for(k=0;k<len*4;k+=4) {
        ret = av_genutil_line_p(_t,p[k+1],p[k+2],p[k+3],p[k+4],&next);
        if(next == 0) break;
        _t-=next;
    }
    return ret;
}


double complex av_genutil_line(double t, double *p) {
    double dummy;
    double complex ret = av_genutil_line_p(t,p[0],p[1],p[2],p[3],&dummy);
    return ret;
}

static double p3d(double x, double a, double b, double c, double d) {
    return a + b*x + c*x*x + d*x*x*x;
}

double complex av_genutil_spline_2d(double t, double *p) {
    int len = p[0];
    double end = p[1];
    int k, index;
    double x,y,t0;
    double tm = fmod(t,end);
    for(k=1;k<len;k++) {
        if(tm<p[k*9+2]) {
            index = (k-1)*9+2;
            t0 = tm-p[index];
            x = p3d(t0,p[index+1],p[index+2],p[index+3],p[index+4]);
            y = p3d(t0,p[index+5],p[index+6],p[index+7],p[index+8]);
            return x + I * y;
        }    
    }
    index = (len-1)*9+2;
    t0 = tm-p[index];
    x = p3d(t0,p[index+1],p[index+2],p[index+3],p[index+4]);
    y = p3d(t0,p[index+5],p[index+6],p[index+7],p[index+8]);
    return x + I * y;
}

double complex av_genutil_circ(double t,  double *p) {
    return p[0] * cos(t) + I * p[1] * sin(t);
}

double complex av_genutil_kardio(double t,  double *p) {
    return p[0] * (cos(t)*cos(t)) + p[1] * cos(t) + I * (p[0] * cos(t) * sin(t) + p[1] * sin(t));
}

double complex av_genutil_lemniskateG(double t,  double *p) {
    double x = cos(t);
    double y = cos(t)*sin(t);
    return p[0] * x + I * p[1] * y;
}

double complex av_genutil_lemniskateB(double t,  double *p) {
    double x = sqrt(cos(2*t))*cos(t);
    double y = sqrt(cos(2*t))*sin(t);
    return p[0] * x + I * p[1] * y;
}


double complex av_genutil_sinusoid(double t,  double *p) {
    double x = pow((sin(p[0]*t)),p[1])*cos(t);
    double y = pow((sin(p[0]*t)),p[1])*sin(t);
    return x + I * y;
}

double complex av_genutil_circoid(double t,  double *p) {
    double x = p[0] * sin(p[1]*t);
    double y = p[2] * cos(p[3]*t);
    return x + I * y;
}

double complex av_genutil_epicycloid(double t, double *p) {
    double a = p[0];
    double b = p[1];
    double x = (a+b)*cos(t) - a * cos((1+b/a)*t);
    double y = (a+b)*sin(t) - a * sin((1+b/a)*t);
    return x + I * y;
}

double complex av_genutil_lissajous(double t, double *p) {
    if(p[3]==0 && p[4] == 0) p[4] = 1.0;
    double a = p[1];
    double b = p[2];
    double c = p[3];
    double d = p[4];
    double x = (c * t + d) * sin(t);
    double y = sin(a * t + b);
    return p[0] * x + I * p[0] * y;
}

double complex av_genutil_lissajousG(double t, double *p) {
    double a = p[1];
    double b = p[2];
    double c = p[3];
    double d = p[4];
    double e = p[5];
    double x = sin(t);
    double y = sin(a*t+b) + c * sin(d*t+e);
    return p[0] * x + I * p[0] * y;
}

double complex av_genutil_lissajousQ(double t, double *p) {
    double a = p[0];
    double b = p[1];
    double c = p[2];
    double d = p[3];
    double x = p[4] * cos(t) + a * sin(b*t) * sin(b*t) + p[6];
    double y = p[5] * sin(t) + c * sin(d*t) * sin(d*t) + p[7];
    return x + I * y;
}

double complex av_genutil_trochoid(double t, double *p) {
    double r1 = p[0];
    double r2 = p[1];
    double o1 = p[2];
    double o2 = p[3];
    double x = r1*cos(o1*t) + r2*cos(o2*t); 
    double y = r1*sin(o1*t) + r2*sin(o2*t); 
    return x + I * y;
}

double complex av_genutil_addoid(double t, double *p) {
    double r1 = p[0];
    double r2 = p[1];
    double r3 = p[2];
    double r4 = p[3];
    double o1 = p[4];
    double o2 = p[5];
    double o3 = p[6];
    double o4 = p[7];
    double x = r1*cos(o1*t) + r2*cos(o2*t); 
    double y = r3*sin(o3*t) + r4*sin(o4*t); 
    return x + I * y;
}


static double poly3(double a,double b,double x) {
    return a*x*x*x + b*x;
}

double complex av_genutil_legendre(double t, double *p) {
    double x = p[0] * poly3(p[1],p[2],cos(t));
    double y = p[0] * poly3(p[1],p[2],sin(t));
    return x + I * y;
}

double complex av_genutil_hypocycloid(double t, double *p) {
    double a = p[1];
    double b = p[2];
    double c = p[3];
    double x = p[0]*((1-c)*cos(a*t)+a*b*cos((1-a)*t));
    double y = p[0]*((1-c)*sin(a*t)-a*b*sin((1-a)*t));
    return x + I * y;
}

double complex av_genutil_rhodonea(double t, double *p) {
    double a = (1-1/p[1])/2;
    double b = 1/a - 1; 
    double x = p[0]*((1-a)*cos(a*t)+a*b*cos((1-a)*t));
    double y = p[0]*((1-a)*sin(a*t)-a*b*sin((1-a)*t));
    return x + I * y;
}

double complex av_genutil_cb(double t, double *p) {
    double a = p[2];
    double b = p[3];
    double x = p[0]*sin(a*t)*cos(t);
    double y = p[1]*sin(b*t)*sin(t);
    return x + I * y;
}

double complex av_genutil_capricornoid(double t, double *p) {
    double a = p[0];
    double x = (sin(t)/(cos(t)+a))*cos(t);
    double y = (sin(t)/(cos(t)+a))*sin(t);
    return x + I * y;
}

double complex av_genutil_scarabaeus(double t, double *p) {
    double a = p[0];
    double b = p[1];
    double r = a*cos(2*t) + b*cos(t);
    return r*cos(t) + I * r*sin(t);
}

double complex av_genutil_nodal(double t, double *p) {
    double a = p[0];
    double r = tan(a*t);
    return r*cos(t) + I * r*sin(t);
}

double complex av_genutil_talbot(double t, double *p) {
    double a = p[0];
    double x = (1+a*sin(t)*sin(t))*cos(t);
    double y = (1-a-a*cos(t)*cos(t))*sin(t);
    return x + I * y;
}

double complex av_genutil_folium(double t, double *p) {
    double a = p[0];
    double r = (sin(t)*sin(t)-a)*cos(t);
    return P2C(r,t);
}

static double super(double t, double a, double b, double c, double d) {
    return pow(pow(fabs(cos(d*t)),a)+pow(fabs(sin(d*t)),b),c);
}

double complex av_genutil_gielis(double t, double *p) {
    double r = super(t,p[0],p[1],p[2],p[3]);
    return P2C(r,t);
}

double complex av_genutil_super_spiral(double t, double *p) {
    double r = exp(p[4]*t) * super(t,p[0],p[1],p[2],p[3]);
    return P2C(r,t);
}

double complex av_genutil_super_rose(double t, double *p) {
    double r = sin(p[4]*t) * super(t,p[0],p[1],p[2],p[3]);
    return P2C(r,t);
}

double complex av_genutil_epi_spiral(double t, double *p) {
    double r = 1/cos(p[0]*t);
    return P2C(r,t);
}

double complex av_genutil_spiral(double t, double *p) {
    double phi = p[0]*t+p[2];
    double r;
    if(phi<0) {
        r = -pow(-phi,p[1]);
        return r*cos(t) - I * r * sin(t);
    }
    r = pow(phi, p[1]);
    return P2C(r,t);
}

double complex av_genutil_atom_spiral(double t, double *p) {
    double r = t/(t-p[0]);
    return P2C(r,t);
}

double complex av_genutil_cotes_spiral(double t, double *p) {
    double r = t/cosh(t*p[0]);
    return P2C(r,t);
}

double complex av_genutil_sin_spiral(double t, double *p) {
    double r = pow(sin(p[0]*t),p[0]);
    return P2C(r,t);
}

double complex av_genutil_maclaurin(double t, double *p) {
    double r = sin(p[0]*t+p[1])/sin(p[2]*t+p[3]);
    return P2C(r,t);
}

double complex av_genutil_harmonic(double t, double *p) {
    double a = p[0];
    double b = p[1];
    double r = (1+a*cos(b*t));
    return P2C(r,t);
}

double complex av_genutil_cornoid(double t, double *p) {
    double a = p[0];
    double b = p[1];
    double c = p[2];
    double x = cos(a*t)*cos(t);
    double y = (c+cos(b*t))*sin(t);
    return x + I * y;
}


/******************** UTIL ****************************/
void av_genutil_draw_number(int x, int y, double z, AVFrame *in, int plane) {
    av_genutil_draw_number_c(x,y,z,1,in,plane);
}
void av_genutil_draw_number_c(int x, int y, double z, double val, AVFrame *in, int plane) {
    if(plane == -1) { // packed rgb
        uint8_t color[] = {255,255,255,255};
        av_plot_number_p(x,y,z,color,in->width,in->height,in->linesize[0],in->data[0]);
    } else {
        av_plot_number_c(x,y,z,val, in->width,in->height,in->linesize[plane],in->data[plane]);
    }
}

void av_genutil_draw_int(int x, int y, int n, AVFrame *in,int plane) {
    av_genutil_draw_int_c(x,y,n,1,in,plane);
}
void av_genutil_draw_int_c(int x, int y, int n, double val, AVFrame *in,int plane) {
    if(plane == -1) { // packed rgb
        uint8_t color[] = {255,255,255,255};
        av_plot_int_p(x,y,n,color,in->width,in->height,in->linesize[0],in->data[0]);
    } else {
        av_plot_int_c(x,y,n,val,in->width,in->height,in->linesize[plane],in->data[plane]);
    }
}

static void drawPoint(GenutilFuncParams *params, AVFrame *in, int plane, double t, PFUNC(f), double a) {
    double complex _z = f(t,params->p);
    double complex z = av_genutil_rotate(_z,params->rot);
    int x = floor(params->fx*creal(z) + params->w/2) + params->x;
    int y = floor(params->fy*cimag(z) + params->h/2) + params->y;
    av_plot_form(params->form,x,y,a,params->w,params->h,in->linesize[plane],in->data[plane]);
}

static void plot_point(GenutilFuncParams *params, AVFrame *in, int plane, long double complex _z, double a) {
    double complex z = av_genutil_rotate(_z,params->rot);
    int x = floor(params->fx*creall(_z) + params->w/2) + params->x;
    int y = floor(params->fy*cimagl(_z) + params->h/2) + params->y;
    av_plot_form(params->form,x,y,a,params->w,params->h,in->linesize[plane],in->data[plane]);
}

void debug(GenutilFuncParams *params, int frame_number, AVFrame *in, int plane) {
    int k;
    for(k=0;k<9;k++) {
        av_genutil_draw_number(k*100+20, params->h-35, params->p[k], in, plane);
    }
    for(k=0;k<3;k++) {
        av_genutil_draw_number(k*100+20, params->h-16, params->colors[0][k], in, plane);
    }
    av_genutil_draw_int(params->w-50, params->h-35, frame_number, in, plane);
}
#define SWAP(x,y)   { \
       int j; \
       j=x; x=y; y=j; \
       }

static void c_line(GenutilFuncParams *params, double *r, double *color, AVFrame *in) {
    double dx;
    double dy;
    double gradient;
    int x0,x1,y0,y1;
    int steep;
    if(params->rot != 0) {
        double complex z0 = av_genutil_rotate(r[0]+I*r[1],params->rot);
        double complex z1 = av_genutil_rotate(r[2]+I*r[3],params->rot);
        r[0] = creal(z0);
        r[1] = cimag(z0);
        r[2] = creal(z1);
        r[3] = cimag(z1);
    }
    x0 = floor(params->fx*r[0] + params->w/2) + params->x;
    y0 = floor(params->fy*r[1] + params->h/2) + params->y;
    x1 = floor(params->fx*r[2] + params->w/2) + params->x;
    y1 = floor(params->fy*r[3] + params->h/2) + params->y;
    steep = abs(y1 - y0) > abs(x1 - x0);
    if (steep) {
        SWAP(x0,y0);
        SWAP(x1,y1);
    }
    if (x0 > x1) {
        SWAP(x0,x1);
        SWAP(y0,y1);
    }

    dx = x1 - x0;
    dy = y1 - y0;
    gradient = dy / dx;

    for(int k = x0;k<=x1;k++) {
        double y = y0 + (k-x0) * gradient;
        for(int plane = 0;plane <3;plane++) {
            if(steep) {
                av_plot_form(params->form,floor(y),k,color[plane],params->w,params->h,in->linesize[plane],in->data[plane]);
            } else {
                av_plot_form(params->form,k,floor(y),color[plane],params->w,params->h,in->linesize[plane],in->data[plane]);
            }
        }
    }
}

static void c_line_c(GenutilFuncParams *params, double *r, double *t, double delta, AVFrame *in) {
    double dx;
    double dy;
    double gradient;
    int x0,x1,y0,y1;
    int steep;
    
    if(params->rot != 0) {
        double complex z0 = av_genutil_rotate(r[0]+I*r[1],params->rot);
        double complex z1 = av_genutil_rotate(r[2]+I*r[3],params->rot);
        r[0] = creal(z0);
        r[1] = cimag(z0);
        r[2] = creal(z1);
        r[3] = cimag(z1);
    }
    x0 = floor(params->fx*r[0] + params->w/2) + params->x;
    y0 = floor(params->fy*r[1] + params->h/2) + params->y;
    x1 = floor(params->fx*r[2] + params->w/2) + params->x;
    y1 = floor(params->fy*r[3] + params->h/2) + params->y;
    steep = abs(y1 - y0) > abs(x1 - x0);
    if (steep) {
        SWAP(x0,y0);
        SWAP(x1,y1);
    }
    if (x0 > x1) {
        SWAP(x0,x1);
        SWAP(y0,y1);
    }

    dx = x1 - x0;
    dy = y1 - y0;
    //delta = delta*sqrt(dx*dx+dy*dy)/dx;
    gradient = dy / dx;
    double colors[3];
    for(int k = x0;k<=x1;k++) {
        av_genutil_get_color(params->cfunc, params->cp,*t,params->cmod, params->is_rgb, colors);
        double y = y0 + (k-x0) * gradient;
        for(int plane = 0;plane <3;plane++) {
            if(steep) {
                av_plot_form(params->form,floor(y),k,colors[plane],params->w,params->h,in->linesize[plane],in->data[plane]);
            } else {
                av_plot_form(params->form,k,floor(y),colors[plane],params->w,params->h,in->linesize[plane],in->data[plane]);
            }
        }
        *t = *t + delta;
    }
}

static void colored_curve(GenutilFuncParams *params, PFUNC(f), AVFrame *in) {
    double t = params->start;
    while(t<params->start+params->length) {
        double colors[3];
        double complex _z = f(t,params->p);
        double complex z = av_genutil_rotate(_z,params->rot);
        int x = floor(params->fx*creal(z) + params->w/2) + params->x;
        int y = floor(params->fy*cimag(z) + params->h/2) + params->y;
        av_genutil_get_color(params->cfunc, params->cp,t,params->cmod, params->is_rgb, colors);
        av_plot_form(params->form,x,y,colors[0],params->w,params->h,in->linesize[0],in->data[0]);
        av_plot_form(params->form,x,y,colors[1],params->w,params->h,in->linesize[1],in->data[1]);
        av_plot_form(params->form,x,y,colors[2],params->w,params->h,in->linesize[2],in->data[2]);
        //drawPoint(params,in,0,t,f,colors[0]);
        //drawPoint(params,in,1,t,f,colors[1]);
        //drawPoint(params,in,2,t,f,colors[2]);
        t+=params->delta;
    }
}

static void colored_curve2(GenutilFuncParams *params, PFUNC(f), AVFrame *in) {
    double t = params->start;
    while(t<params->start+params->length) {
        double colors[3];
        av_genutil_get_color(params->cfunc, params->cp,t,params->cmod, params->is_rgb, colors);
        drawPoint(params,in,0,t,f,colors[0]);
        drawPoint(params,in,1,t,f,colors[1]);
        drawPoint(params,in,2,t,f,colors[2]);
        drawPoint(params,in,0,-t,f,colors[0]);
        drawPoint(params,in,1,-t,f,colors[1]);
        drawPoint(params,in,2,-t,f,colors[2]);
        t+=params->delta;
    }
}
    
static void logParameters(AVFilterContext *ctx,double* p, int len, int debug) {
    int k;
    for(k=0;k<len;k++) {
        av_log(ctx, debug?AV_LOG_DEBUG:AV_LOG_INFO," %f",p[k]); 
    }
    av_log(ctx, debug?AV_LOG_DEBUG:AV_LOG_INFO,"\n"); 
}


static void colored_spline(GenutilFuncParams *params, AVFrame *in) {
    int len = params->p[0];
    double *tarr = calloc(len,sizeof(double));
    tarr[0] = 0;
    for(int i=1;i<len;i++) {
        double lx = params->p[1+i*2] - params->p[1+(i-1)*2];
        double ly = params->p[2+i*2] - params->p[2+(i-1)*2];
        if(lx == 0.0) tarr[i] = fabs(ly);
        else if(ly == 0.0) tarr[i] = fabs(lx);
        else tarr[i]=sqrt(lx*lx+ly*ly);
        tarr[i]+=tarr[i-1];
    }
    double *inputX = calloc(len*2+1,sizeof(double));
    double *inputY = calloc(len*2+1,sizeof(double));
    inputX[0]=len;
    inputY[0]=len;
    for(int i=0;i<len;i++) {
        inputX[1+i*2] = tarr[i];
        inputX[2+i*2] = params->p[1+i*2];
        inputY[1+i*2] = tarr[i];
        inputY[2+i*2] = params->p[2+i*2];
    }
    double *outputX = calloc(2+len*5,sizeof(double));
    double *outputY = calloc(2+len*5,sizeof(double));
    av_genutil_spline(inputX,outputX);
    av_genutil_spline(inputY,outputY);
    double *splineInput = calloc(2+inputX[0]*9,sizeof(double));
    splineInput[0]=outputX[0];
    splineInput[1]=outputX[1];
    for(int i=0;i<outputX[0];i++) {
        splineInput[2+i*9]=outputX[2+i*5];
        splineInput[3+i*9]=outputX[3+i*5];
        splineInput[4+i*9]=outputX[4+i*5];
        splineInput[5+i*9]=outputX[5+i*5];
        splineInput[6+i*9]=outputX[6+i*5];
        splineInput[7+i*9]=outputY[3+i*5];
        splineInput[8+i*9]=outputY[4+i*5];
        splineInput[9+i*9]=outputY[5+i*5];
        splineInput[10+i*9]=outputY[6+i*5];
    }
    logParameters(params->ctx,outputX,2+len*5,1);
    logParameters(params->ctx,outputY,2+len*5,1);
    logParameters(params->ctx,splineInput,2+outputX[0]*9,1);
    double t = params->start;
    while(t<params->start+params->length) {
        double colors[3];
        double complex _z = av_genutil_spline_2d(t,splineInput);
        double complex z = av_genutil_rotate(_z,params->rot);
        int x = floor(params->fx*creal(z) + params->w/2) + params->x;
        int y = floor(params->fy*cimag(z) + params->h/2) + params->y;
        av_genutil_get_color(params->cfunc, params->cp,t,params->cmod, params->is_rgb, colors);
        av_plot_form(params->form,x,y,colors[0],params->w,params->h,in->linesize[0],in->data[0]);
        av_plot_form(params->form,x,y,colors[1],params->w,params->h,in->linesize[1],in->data[1]);
        av_plot_form(params->form,x,y,colors[2],params->w,params->h,in->linesize[2],in->data[2]);
        t+=params->delta;
    }
    free(tarr);
    free(inputX);
    free(inputY);
    free(splineInput);
    free(outputX);
    free(outputY);
}

static void set_alpha(GenutilFuncParams *params, AVFrame *in, double t, PFUNC(f), double a,double alpha,int ofsX, int ofsY) {
    double complex _z = f(t,params->p);
    double complex z = av_genutil_rotate(_z,alpha);
    int x = floor(params->fx*creal(z) + params->w/2) + params->x + ofsX;
    int y = floor(params->fy*cimag(z) + params->h/2) + params->y + ofsY;
    if(params->is_packed_rgb) {
        av_plot_form_p(params->form,x,y,a,params->w,params->h,in->linesize[0],in->data[0]+params->rgba_map[A]);
    } else {
        av_plot_form(params->form,x,y,a,params->w,params->h,in->linesize[3],in->data[3]);
    }
}

static void curve(GenutilFuncParams *params, PFUNC(f),int n, AVFrame *in) {
    double count = params->count;
    double speed = params->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int j;
        double start = 0;
        int layers = params->layers;
        for(j=0;j<layers;j++) {
            double s = start + params->length*j/layers;
            double t = s + n*speed;
            set_alpha(params,in,t,f,1,alpha,0,0);
        }
    }
}

static void pcurve(GenutilFuncParams *params, PFUNC(f),int n, AVFrame *in) {
    double count = params->count;
    double speed = params->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int j;
        double start = 0;
        int layers = params->layers;
        int clen = params->length*speed;
        for(j=0;j<layers;j++) {
            double s = start + params->length*j/layers;
            double t = s + ((n%clen)+1)*speed;
            set_alpha(params,in,t,f,1,alpha,0,0);
        }
    }
}


static void curve2(GenutilFuncParams *params, PFUNC(f),int n, AVFrame *in) {
    double count = params->count;
    double speed = params->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int j;
        double start = 0;
        int layers = params->layers;
        for(j=0;j<layers;j++) {
            double s = start + params->length*j/layers;
            double t = s + n*speed;
            set_alpha(params,in,t,f,1,alpha,0,0);
            set_alpha(params,in,-t,f,1,alpha,0,0);
        }
    }
}

static void rnd_curve(GenutilFuncParams *params, PFUNC(f),int n, AVFrame *in) {
    double count = params->count;
    int k;
    for(k=0;k<count;k++) {
        double rnd1 = ((double)rand())/RAND_MAX;
        double rnd3 = ((double)rand())/RAND_MAX;
        double rnd4 = ((double)rand())/RAND_MAX;
        double alpha = rnd1 * 2 * M_PI;
        int j;
        double speed = params->speed * rnd4;
        double start = rnd3*params->length;
        int layers = params->layers;
        for(j=0;j<layers;j++) {
            double rnd2 = ((double)rand())/RAND_MAX;
            double s = start + rnd2*params->length*j/layers;
            double t = s + n*speed;
            set_alpha(params,in,t,f,1,alpha,0,0);
        }
    }
}


/******************** FFUNCS **************************/

typedef struct GenutilFunc {
    const char *name;
    void (*f)(GenutilFuncParams*, int,AVFrame *in);
} GenutilFunc;

void gzero(GenutilFuncParams *params, int n, AVFrame *in) {
}

static void genutil_hypo(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_hypocycloid,in);
}
static void genutil_lissg(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_lissajousG,in);
}

static void genutil_lissq(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_lissajousQ,in);
}

static void genutil_liss(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_lissajous,in);
}

static void genutil_leg(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_legendre,in);
}

static void genutil_cb(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_cb,in);
}

static void genutil_sinoid(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_sinusoid,in);
}

static void genutil_capri(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_capricornoid,in);
}

static void genutil_scara(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_scarabaeus,in);
}

static void genutil_tro(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_trochoid,in);
}

static void genutil_circ(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_circoid,in);
}

static void genutil_rhodo(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_rhodonea,in);
}

static void genutil_gielis(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_gielis,in);
}

static void genutil_super_rose(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_super_rose,in);
}

static void genutil_super_spiral(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_super_rose,in);
}

static void genutil_epi_spiral(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_epi_spiral,in);
}

static void genutil_spiral(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve2(params,av_genutil_spiral,in);
}

static void genutil_maclaurin(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_maclaurin,in);
}

static void genutil_harmonic(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_harmonic,in);
}

static void genutil_cornoid(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_cornoid,in);
}

static void genutil_line(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_line,in);
}

static void genutil_poly(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_poly,in);
}

static void genutil_spline(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_spline(params,in);
}

static void genutil_addoid(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_addoid,in);
}

static double *gen_poly(int len, int size, int (*f)(int,int,double*), int n, double *pm, int mode) {
    double *p = calloc(len*2+4,sizeof(double));
    int k=3;
    int j;
    int last = -1;
    int c = 0;
    double s6 = sqrt(3)/2;
    for(j=0;c<len;j++) {
        if(mode == 0) {
            int rnd = f(j,n,pm)%8;
            switch(rnd) {
                case 0: if(last!=4) {p[k] = p[k-2]; p[k+1] = p[k-1] - size;k+=2;last = rnd;c++;} break;
                case 1: if(last!=5) {p[k] = p[k-2] - size; p[k+1] = p[k-1] - size;k+=2;last = rnd;c++;} break;
                case 2: if(last!=6) {p[k] = p[k-2] - size; p[k+1] = p[k-1];k+=2;last = rnd;c++;} break;
                case 3: if(last!=7) {p[k] = p[k-2] - size; p[k+1] = p[k-1] + size;k+=2;last = rnd;c++;} break;
                case 4: if(last!=0) {p[k] = p[k-2]; p[k+1] = p[k-1] + size;k+=2;last = rnd;c++;} break;
                case 5: if(last!=1) {p[k] = p[k-2] + size; p[k+1] = p[k-1] + size;k+=2;last = rnd;c++;} break;
                case 6: if(last!=2) {p[k] = p[k-2] + size; p[k+1] = p[k-1];k+=2;last = rnd;c++;} break;
                case 7: if(last!=3) {p[k] = p[k-2] + size; p[k+1] = p[k-1] - size;k+=2;last = rnd;c++;} break;
            }
        } else if(mode == 1) {
            int rnd = f(j,n,pm)%4;
            switch(rnd) {
                case 0: if(last!=1) {p[k] = p[k-2]; p[k+1] = p[k-1] + size;k+=2;last = rnd;c++;} break;
                case 1: if(last!=0) {p[k] = p[k-2]; p[k+1] = p[k-1] - size;k+=2;last = rnd;c++;} break;
                case 2: if(last!=3) {p[k] = p[k-2] + size; p[k+1] = p[k-1];k+=2;last = rnd;c++;} break;
                case 3: if(last!=2) {p[k] = p[k-2] - size; p[k+1] = p[k-1];k+=2;last = rnd;c++;} break;
            }
        } else if(mode == 2) {
            int rnd = f(j,n,pm)%8;
            switch(rnd) {
                case 0: p[k] = p[k-2]; p[k+1] = p[k-1] - size;k+=2;c++; break;
                case 4: p[k] = p[k-2]; p[k+1] = p[k-1] + size;k+=2;c++; break;
                case 1: p[k] = p[k-2] - size; p[k+1] = p[k-1] - size;k+=2;c++; break;
                case 5: p[k] = p[k-2] + size; p[k+1] = p[k-1] + size;k+=2;c++; break;
                case 2: p[k] = p[k-2] - size; p[k+1] = p[k-1];k+=2;c++; break;
                case 6: p[k] = p[k-2] + size; p[k+1] = p[k-1];k+=2;c++; break;
                case 3: p[k] = p[k-2] - size; p[k+1] = p[k-1] + size;k+=2;c++; break;
                case 7: p[k] = p[k-2] + size; p[k+1] = p[k-1] - size;k+=2;c++; break;
            }
        } else if(mode == 3) {
            int rnd = f(j,n,pm)%4;
            switch(rnd) {
                case 0: p[k] = p[k-2]; p[k+1] = p[k-1] + size;k+=2;c++; break;
                case 1: p[k] = p[k-2]; p[k+1] = p[k-1] - size;k+=2;c++; break;
                case 2: p[k] = p[k-2] + size; p[k+1] = p[k-1];k+=2;c++; break;
                case 3: p[k] = p[k-2] - size; p[k+1] = p[k-1];k+=2;c++; break;
            }
        } else if(mode == 4) {
            int rnd = f(j,n,pm)%6;
            switch(rnd) {
                case 0: if(last!=4) {p[k] = p[k-2] - size; p[k+1] = p[k-1];k+=2;last = rnd;c++;} break;
                case 4: if(last!=0) {p[k] = p[k-2] + size; p[k+1] = p[k-1];k+=2;last = rnd;c++;} break;
                case 1: if(last!=5) {p[k] = p[k-2] - size*0.5; p[k+1] = p[k-1] - size*s6;k+=2;last = rnd;c++;} break;
                case 5: if(last!=1) {p[k] = p[k-2] + size*0.5; p[k+1] = p[k-1] + size*s6;k+=2;last = rnd;c++;} break;
                case 2: if(last!=3) {p[k] = p[k-2] - size*0.5; p[k+1] = p[k-1] + size*s6;k+=2;last = rnd;c++;} break;
                case 3: if(last!=2) {p[k] = p[k-2] + size*0.5; p[k+1] = p[k-1] - size*s6;k+=2;last = rnd;c++;} break;
             }
         } else if(mode == 5) {
            int rnd = f(j,n,pm)%12;
            switch(rnd) {
                case 0: if(last!=6) {p[k] = p[k-2]; p[k+1] = p[k-1] - size;k+=2;last = rnd;c++;} break;
                case 1: if(last!=7) {p[k] = p[k-2] - size*0.5; p[k+1] = p[k-1] - size*s6;k+=2;last = rnd;c++;} break;
                case 2: if(last!=8) {p[k] = p[k-2] - size*s6; p[k+1] = p[k-1] - size*0.5;k+=2;last = rnd;c++;} break;
                case 3: if(last!=9) {p[k] = p[k-2] - size; p[k+1] = p[k-1];k+=2;last = rnd;c++;} break;
                case 4: if(last!=10) {p[k] = p[k-2] - size*s6; p[k+1] = p[k-1] + size*0.5;k+=2;last = rnd;c++;} break;
                case 5: if(last!=11) {p[k] = p[k-2] - size*0.5; p[k+1] = p[k-1] + size*s6;k+=2;last = rnd;c++;} break;
                case 6: if(last!=0) {p[k] = p[k-2]; p[k+1] = p[k-1] + size;k+=2;last = rnd;c++;} break;
                case 7: if(last!=1) {p[k] = p[k-2] + size*0.5; p[k+1] = p[k-1] + size*s6;k+=2;last = rnd;c++;} break;
                case 8: if(last!=2) {p[k] = p[k-2] + size*s6; p[k+1] = p[k-1] + size*0.5;k+=2;last = rnd;c++;} break;
                case 9: if(last!=3) {p[k] = p[k-2] + size; p[k+1] = p[k-1];k+=2;last = rnd;c++;} break;
                case 10: if(last!=4) {p[k] = p[k-2] + size*s6; p[k+1] = p[k-1] - size*0.5;k+=2;last = rnd;c++;} break;
                case 11: if(last!=5) {p[k] = p[k-2] + size*0.5; p[k+1] = p[k-1] - size*s6;k+=2;last = rnd;c++;} break;
            }
        }
    }
    p[k] = p[k-2]; p[k+1] = p[k-1];
    p[0] = len;
    return p;
}

static int gp_rnd(int k, int n, double *p) {
    return rand();
}

static int gp_sin(int k, int n, double *p) {
    return (int)(floor(sin(n*k*p[0])*8));
}

static int gp_mod(int k, int n, double *p) {
    return (n%(k+1));
}

static int gp_mod1(int k, int n, double *p) {
    return ((n*k)%(k+1));
}

static int gp_q(int k, int n, double *p) {
    return k*(k*p[0]) *n*(n*p[1]);
}

static int gp_mul(int k, int n, double *p) {
    return ((int)floor(k*n*p[0]));
}

static int gp_ulam(int k, int n, double *p) {
    int j;
    int last = n;
    for(j=0;j<k;j++) {
        if(last%2==0) last = last/2;
        else last = 3*last+1;
    }
    return last;
}

static void interpolate_poly(double *p1, double *p2, double t) {
    t = (sin(t*M_PI-M_PI/2)+1)/2;
    for(int k=3; k<(p1[0]*2+2); k+=2) {
        p1[k] += t*(p2[k]-p1[k]);
        p1[k+1] += t*(p2[k+1]-p1[k+1]);
    }
}

static void genutil_rndpoly(GenutilFuncParams *params, int n, int (*f)(int,int,double*), AVFrame *in) {
    int len = floor(params->p[0]);
    int size = params->p[1];
    int interval = params->p[2];
    int pause = params->p[3];
    int steplen = interval+pause;
    int step = n/steplen;
    srand(step);
    double *p1 = gen_poly(len,size,f,step,&params->p[4],params->mode);
    int nn = n%steplen;
    if(nn >= pause) {
        srand(step+1);
        double *p2 = gen_poly(len,size,f,step+1,&params->p[4],params->mode);
        interpolate_poly(p1,p2,(double)(nn-pause)/(double)interval);
        free(p2);
    }
    int spline = params->spline;
    double *ps = params->p;
    params->p = p1;
    params->length = poly_length(p1);
    if(spline) {
      colored_spline(params,in);
    } else {
      colored_curve(params,av_genutil_poly_s,in);
    }
    params->p = ps;
    free(p1);
}

static void genutil_rndpoly_p(GenutilFuncParams *params, int n, int (*f)(int,int,double*), AVFrame *in) {
    int k;
    int count = params->count;
    double rot = params->rot;
    for(k=0;k<count;k++) {
        params->rot = rot + k*2*M_PI/count;
        genutil_rndpoly(params,n,f,in);
    }
}

static void genutil_poly_mul(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_rndpoly_p(params,n,gp_mul,in);
}

static void genutil_poly_mod(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_rndpoly_p(params,n,gp_mod,in);
}

static void genutil_poly_mod1(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_rndpoly_p(params,n,gp_mod1,in);
}

static void genutil_poly_sin(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_rndpoly_p(params,n,gp_sin,in);
}

static void genutil_poly_q(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_rndpoly_p(params,n,gp_q,in);
}

static void genutil_poly_rnd(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_rndpoly_p(params,n,gp_rnd,in);
}

static void genutil_poly_ulam(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_rndpoly_p(params,n,gp_ulam,in);
}

static double *gen_poly_d(int len, int size, double (*f)(int,int,double*), int n, double *pm) {
    double *p = calloc(len*2+4,sizeof(double));
    int k=3;
    int j;
    for(j=0;j<len;j++) {
        double phi = f(j,n,pm);
        p[k] = p[k-2] + (double)size * cos(phi); 
        p[k+1] = p[k-1] + (double)size * sin(phi);
        k+=2;
    }
    p[k] = p[k-2]; p[k+1] = p[k-1];
    p[0] = len;
    return p;
}


static void genutil_rndpoly_d(GenutilFuncParams *params, int n, double (*f)(int,int,double*), AVFrame *in) {
    int len = floor(params->p[0]);
    int size = params->p[1];
    int interval = params->p[2];
    int pause = params->p[3];
    int steplen = interval+pause;
    int step = n/steplen;
    srand(step);
    double *p1 = gen_poly_d(len,size,f,step,&params->p[4]);
    int nn = n%steplen;
    if(nn > pause) {
        srand(step+1);
        double *p2 = gen_poly_d(len,size,f,step+1,&params->p[4]);
        interpolate_poly(p1,p2,(double)(nn-pause)/(double)interval);
        free(p2);
    }
    double *ps = params->p;
    int spline = params->spline;
    params->p = p1;
    params->length = poly_length(p1);
    if(spline) {
      colored_spline(params,in);
    } else {
      colored_curve(params,av_genutil_poly_s,in);
    }
    params->p = ps;
    free(p1);
}

static void genutil_rndpoly_p_d(GenutilFuncParams *params, int n, double (*f)(int,int,double*), AVFrame *in) {
    int k;
    int count = params->count;
    double rot = params->rot;
    for(k=0;k<count;k++) {
        params->rot = rot + k*2*M_PI/count;
        genutil_rndpoly_d(params,n,f,in);
    }
}

static double dfrnd(int k, int n, double *pm) {
    return ((double)rand()/RAND_MAX)*2*M_PI; 
}

static void genutil_poly_drnd(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_rndpoly_p_d(params,n,dfrnd,in);
}

static double dfsin(int k, int n, double *pm) {
    return sin(k*n*pm[0])*pm[1]*M_PI;
}

static void genutil_poly_dsin(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_rndpoly_p_d(params,n,dfsin,in);
}

static double dfsink(int k, int n, double *pm) {
    return sin(k*pm[0])*n*pm[1]*M_PI;
}

static void genutil_poly_dsink(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_rndpoly_p_d(params,n,dfsink,in);
}

static double dfq(int k, int n, double *pm) {
    return fmod((double)k*k*n*pm[0], 2*M_PI);
}

static void genutil_poly_dq(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_rndpoly_p_d(params,n,dfq,in);
}

static double dfqs(int k, int n, double *pm) {
    return fmod((double)k*sin(k*pm[1])*n*pm[0], 2*M_PI);
}

static void genutil_poly_dqs(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_rndpoly_p_d(params,n,dfqs,in);
}

/******************************* hilbert *************************************/
static void move(uint8_t val, int *pos, uint8_t *buf) {
    buf[*pos]=val;
    (*pos) += 1;
}
enum {
  UP,
  LEFT,
  DOWN,
  RIGHT,
};
static void hilbert(int level,int direction, uint8_t *buf,int *k)
{
  if (level==1) {
    switch (direction) {
    case LEFT:
      move(RIGHT,k,buf);
      move(DOWN,k,buf);
      move(LEFT,k,buf);
      break;
    case RIGHT:
      move(LEFT,k,buf);
      move(UP,k,buf);
      move(RIGHT,k,buf);
      break;
    case UP:
      move(DOWN,k,buf);
      move(RIGHT,k,buf);
      move(UP,k,buf);
      break;
    case DOWN:
      move(UP,k,buf);
      move(LEFT,k,buf);
      move(DOWN,k,buf);
      break;
    } /* switch */
  } else {
    switch (direction) {
    case LEFT:
      hilbert(level-1,UP,buf,k);
      move(RIGHT,k,buf);
      hilbert(level-1,LEFT,buf,k);
      move(DOWN,k,buf);
      hilbert(level-1,LEFT,buf,k);
      move(LEFT,k,buf);
      hilbert(level-1,DOWN,buf,k);
      break;
    case RIGHT:
      hilbert(level-1,DOWN,buf,k);
      move(LEFT,k,buf);
      hilbert(level-1,RIGHT,buf,k);
      move(UP,k,buf);
      hilbert(level-1,RIGHT,buf,k);
      move(RIGHT,k,buf);
      hilbert(level-1,UP,buf,k);
      break;
    case UP:
      hilbert(level-1,LEFT,buf,k);
      move(DOWN,k,buf);
      hilbert(level-1,UP,buf,k);
      move(RIGHT,k,buf);
      hilbert(level-1,UP,buf,k);
      move(UP,k,buf);
      hilbert(level-1,RIGHT,buf,k);
      break;
    case DOWN:
      hilbert(level-1,RIGHT,buf,k);
      move(UP,k,buf);
      hilbert(level-1,DOWN,buf,k);
      move(LEFT,k,buf);
      hilbert(level-1,DOWN,buf,k);
      move(DOWN,k,buf);
      hilbert(level-1,LEFT,buf,k);
      break;
    } /* switch */
  } /* if */
}

static void genutil_hilbert(GenutilFuncParams *params, int n, AVFrame *in) {
    double *_p = params->p;
    int level = floor(_p[0]);
    int size = _p[1];
    double *p = calloc(8<<(level*2),sizeof(double));
    uint8_t *buf = calloc(3<<(level*2),sizeof(uint8_t));
    int len = 1<<(level*2);
    int pos = 0;
    hilbert(level,0,buf,&pos);
    int k,j;
    k=3;
    p[0] = len;
    for(j=0;j<len;j++) {
        switch(buf[j]) {
            case 0: p[k] = p[k-2]; p[k+1] = p[k-1] - size;k+=2; break;
            case 1: p[k] = p[k-2] - size; p[k+1] = p[k-1];k+=2; break;
            case 2: p[k] = p[k-2]; p[k+1] = p[k-1] + size;k+=2; break;
            case 3: p[k] = p[k-2] + size; p[k+1] = p[k-1];k+=2; break;
        }
    }
    params->p = p;
    colored_curve(params,av_genutil_poly,in);
    params->p = _p;
    free(buf);
    free(p);
}

/***** PLOT *********/
static void genutil_pkardioids(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_kardio,n,in);
}

static void genutil_plemg(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_lemniskateG,n,in);
}

static void genutil_plemb(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_lemniskateB,n,in);
}

static void genutil_pepi(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_epicycloid,n,in);
}

static void genutil_plissg(GenutilFuncParams *params, int n, AVFrame *in) {
    double t = 0;
    while(t<params->length) {
        int ofsX = params->x;
        int ofsY = params->y;
        set_alpha(params,in,t,av_genutil_lissajousG,1,0,ofsX,ofsY);
        t+=params->speed;
    }
}


static void genutil_plissq(GenutilFuncParams *params, int n, AVFrame *in) {
    double t = 0;
    while(t<params->length) {
        int ofsX = params->x;
        int ofsY = params->y;
        set_alpha(params,in,t,av_genutil_lissajousQ,1,0,ofsX,ofsY);
        t+=params->speed;
    }
}

static void genutil_plissp(GenutilFuncParams *params, int n, AVFrame *in) {
    double t = 0;
    while(t<params->length) {
        int ofsX = 0;
        int ofsY = 0;
        set_alpha(params,in,t,av_genutil_lissajous,1,0,ofsX,ofsY);
        t+=params->speed;
    }
}

static void genutil_pliss(GenutilFuncParams *params, int n, AVFrame *in) {
    double count = params->count;
    double speed = params->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int ofsX = 0;
        int ofsY = 0;
        int j;
        double start = 0;
        int layers = params->layers;
        for(j=0;j<layers;j++) {
            double s = start + params->length*j/layers;
            double t = s + n*speed;
            set_alpha(params,in,t,av_genutil_lissajous,1,alpha,ofsX,ofsY);
        }
    }
}

static void genutil_pleg(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_legendre,n,in);
}

static void genutil_phypo(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_hypocycloid,n,in);
}

static void genutil_prhodo(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_rhodonea,n,in);
}

static void genutil_pnodal(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_nodal,n,in);
}

static void genutil_ptalbot(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_talbot,n,in);
}

static void genutil_pfolium(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_folium,n,in);
}

static void genutil_pgielis(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_gielis,n,in);
}

static void genutil_psuper_spiral(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_super_spiral,n,in);
}

static void genutil_psuper_rose(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_super_rose,n,in);
}

static void genutil_pcb(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_cb,n,in);
}

static void genutil_pepi_spiral(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_epi_spiral,n,in);
}

static void genutil_pspiral(GenutilFuncParams *params, int n, AVFrame *in) {
    curve2(params,av_genutil_spiral,n,in);
}

static void genutil_pmaclaurin(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_maclaurin,n,in);
}

static void genutil_pharmonic(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_harmonic,n,in);
}

static void genutil_pcornoid(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_cornoid,n,in);
}

static void genutil_psinoid(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_sinusoid,n,in);
}

static void genutil_pcapri(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_capricornoid,n,in);
}

static void genutil_pscara(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_scarabaeus,n,in);
}

static void genutil_ptro(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_trochoid,n,in);
}

static void genutil_pcirco(GenutilFuncParams *params, int n, AVFrame *in) {
    colored_curve(params,av_genutil_circoid,in);
}

static void genutil_paddoid(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_addoid,n,in);
}

static void genutil_pline(GenutilFuncParams *params, int n, AVFrame *in) {
    pcurve(params,av_genutil_line,n,in);
}

static void genutil_ppoly(GenutilFuncParams *params, int n, AVFrame *in) {
    curve(params,av_genutil_poly,n,in);
}

static void genutil_ppolyu(GenutilFuncParams *params, int n, AVFrame *in) {
    rnd_curve(params,av_genutil_poly_u,n,in);
}

static void genutil_prndpoly_p(GenutilFuncParams *params, int (*f)(int,double), int n, AVFrame *in) {
    double *_p = params->p;
    int len = floor(_p[0]);
    int size = _p[1];
    double param = _p[2];
    double *p = calloc(len*2+4, sizeof(double));
    p[0] = len;
    int k=3;
    int j;
    int last = -1;
    for(j=0;j<len;j++) {
        //int rnd = rand() % 8;
        int rnd = f(j,param);
        switch(rnd) {
            case 0: if(last!=4) {p[k] = p[k-2]; p[k+1] = p[k-1] - size;k+=2;last = rnd;} break;
            case 4: if(last!=0) {p[k] = p[k-2]; p[k+1] = p[k-1] + size;k+=2;last = rnd;} break;
            case 1: if(last!=5) {p[k] = p[k-2] - size; p[k+1] = p[k-1] - size;k+=2;last = rnd;} break;
            case 5: if(last!=1) {p[k] = p[k-2] + size; p[k+1] = p[k-1] + size;k+=2;last = rnd;} break;
            case 2: if(last!=6) {p[k] = p[k-2] - size; p[k+1] = p[k-1];k+=2;last = rnd;} break;
            case 6: if(last!=2) {p[k] = p[k-2] + size; p[k+1] = p[k-1];k+=2;last = rnd;} break;
            case 3: if(last!=7) {p[k] = p[k-2] - size; p[k+1] = p[k-1] + size;k+=2;last = rnd;} break;
            case 7: if(last!=3) {p[k] = p[k-2] + size; p[k+1] = p[k-1] - size;k+=2;last = rnd;} break;
        }
    }
    //p[k] = p[k-2]; p[k+1] = p[k-1];
    params->p = p;
    curve(params,av_genutil_poly,n,in);
    free(p);
}
static int dfunc(int k, double p) {
    //if(k%2==0) return  ((((n*n)%(k+1))*k)%4)+4 ;
    //return (((n*n)%(k+1))*k)%4;
    //return (int)(floor(sin(n*k)*8))%8;
    //return (n%(k+1))%8;
    //return ((int)floor((k/8)*n*n/(n+1)))%8;
    //return ((n*k)%(k+1))%8;
    //return ((int)floor(k*1.01*n))%8;
    return ((int)floor(k*p))%8;
}
static void genutil_prndpoly(GenutilFuncParams *params,int n, AVFrame *in) {
     genutil_prndpoly_p(params,dfunc,n,in);
}

static void genutil_philbert(GenutilFuncParams *params, int n, AVFrame *in) {
    double *_p = params->p;
    int level = floor(_p[0]);
    int size = _p[1];
    double *p = calloc(4<<(level*2),sizeof(double));
    uint8_t *buf = calloc(2<<(level*2),sizeof(uint8_t));
    int len = 1<<(level*2);
    int pos = 0;
    hilbert(level,0,buf,&pos);
    int k,j;
    k=3;
    p[0] = len-1;
    for(j=0;j<len;j++) {
        switch(buf[j]) {
            case 0: p[k] = p[k-2]; p[k+1] = p[k-1] - size;k+=2; break;
            case 1: p[k] = p[k-2] - size; p[k+1] = p[k-1];k+=2; break;
            case 2: p[k] = p[k-2]; p[k+1] = p[k-1] + size;k+=2; break;
            case 3: p[k] = p[k-2] + size; p[k+1] = p[k-1];k+=2; break;
        }
    }
    params->p = p;
    curve(params,av_genutil_poly,n,in);
    params->p = _p;
    free(buf);
    free(p);
}

static void genutil_circs(GenutilFuncParams *params, int n, AVFrame *in) {
    double count = params->count;
    double speed = params->speed;
    int k;
    for(k=0;k<count;k++) {
        double alpha = (k/count) * 2 * M_PI;
        int ofsX = floor(cos(alpha) * params->p[0]);
        int ofsY = floor(sin(alpha) * params->p[1]);
        double start = 2*M_PI*k/count;
        int j;
        int layers = params->layers;
        for(j=0;j<layers;j++) {
            double s = start + params->length*j/layers;
            double t = s + n*speed;
            set_alpha(params,in,t,av_genutil_circ,1,0,ofsX,ofsY);
        }
    }
} 

// ============ AFFIN / FLAMES (later) =======/
static double complex affin(double complex z, double *p) {
    double x = p[0] * creal(z) + p[1] * cimag(z) + p[2];
    double y = p[3] * creal(z) + p[4] * cimag(z) + p[5];
    return x + I * y;
}

static void genutil_frac(GenutilFuncParams *params, int n, AVFrame *in) {
    double *p = params->p;
    int len = floor(p[0]);
    double complex start = p[1] + I * p[2];
    double p0 = p[3];
    double p1 = p[4];
    double p2 = p[5];
    double *w0 = &p[6];
    double *w1 = &p[12];
    double *w2 = &p[18];
    double *w3 = &p[24];
    int k;
    plot_point(params,in,0,start,1);
    double complex z = start;
    double complex z1;
    double r;
    for(k=0;k<len;k++) {
        r = (double)rand()/RAND_MAX;
        if(r<p0) z1 = affin(z,w0);
        else if(r<p1) z1 = affin(z,w1);
        else if(r<p2) z1 = affin(z,w2);
        else z1 = affin(z,w3);
        if(k>20)
            plot_point(params,in,0,creal(z1) - I * cimag(z1),1);
        z = z1;
    }
}



// ================= other maps (mira hopalong ... ) =================/
//

static void inc_pix(int x, int y, int w, int h,int *buf, int *max) {
    if(x>=0 && x<w && y>=0 && y<h) {
        buf[x + w * y]++;
	if(buf[x + w * y]>*max) *max = buf[x + w * y];
    }
}
static void inc_point(GenutilFuncParams *params, long double complex _z, int *buf, int*max) {
    double complex z = av_genutil_rotate(_z,params->rot);
    int x = floor(params->fx*creall(_z) + params->w/2) + params->x;
    int y = floor(params->fy*cimagl(_z) + params->h/2) + params->y;
    inc_pix(x,y,params->w,params->h,buf,max);
}


static void genutil_attr(GenutilFuncParams *params, AVFrame *in, long double complex (*f)(long double complex, double*)) {
    int *buf = calloc(params->w*params->h,sizeof(int));
    int max = 0;
    double *p = params->p;
    int len = floor(p[0]);
    long double complex start = p[1] + I * p[2];
    int k;
    long double complex z = start;
    long double complex z1;
    double colors[3];
    for(k=0;k<len;k++) {
        z1 = f(z,&p[5]);
        inc_point(params,creall(z1) - I * cimagl(z1),buf, &max);
        z = z1;
    }
    av_log(params->ctx, AV_LOG_INFO, "max: %d\n", max);
    int x,y;
    for(x=0;x<params->w;x++) {
        for(y=0;y<params->h;y++) {
	    int val = buf[x + params->w*y];
	    if(val>0) {
            av_genutil_get_color(params->cfunc, params->cp,val*p[3]+p[4],params->cmod, params->is_rgb, colors);
            av_plot_form(params->form,x,y,colors[0],params->w,params->h,in->linesize[0],in->data[0]);
            av_plot_form(params->form,x,y,colors[1],params->w,params->h,in->linesize[1],in->data[1]);
            av_plot_form(params->form,x,y,colors[2],params->w,params->h,in->linesize[2],in->data[2]);		 }

	}
    }
    free(buf);
}
static double miraf(long double x, double a) {
    return a*x-(1-a)*(2*x*x/(1+x*x));
}

static long double complex mira(long double complex z, double *p) {
    double a = p[0];
    double b = p[1];
    long double x = b * cimagl(z) + miraf(creall(z),a);
    long double y = -creall(z) + miraf(x,a);
    return x + I * y;
}

static void genutil_mira(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_attr(params,in,mira);
}

static long double complex kaneko1(long double complex z, double *p) {
    double a = p[0];
    double b = p[1];
    long double x = a * creall(z) + (1-a)*(1-b*cimagl(z)*cimagl(z));
    long double y = creall(z);
    return x + I * y;
}

static void genutil_kaneko1(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_attr(params,in,kaneko1);
}

static long double complex kaneko2(long double complex z, double *p) {
    double a = p[0];
    double b = p[1];
    long double x = a * creall(z) + (1-a)*(1-b*fabsl(cimagl(z)));
    long double y = creall(z);
    return x + I * y;
}

static void genutil_kaneko2(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_attr(params,in,kaneko2);
}

static long double complex gingerbreadman(long double complex z, double *p) {
    long double x = 1 - cimagl(z) + fabs(creall(z));
    long double y = creall(z);
    return x + I * y;
}

static void genutil_ginger(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_attr(params,in,gingerbreadman);
}

static long double complex hopalong(long double complex z, double *p) {
    long double a = p[0];
    long double b = p[1];
    long double c = p[2];
    long double x0 = creall(z);
    long double y0 = cimagl(z);
    long double x =  y0 - SIGN(x0)*sqrt(fabsl(b*x0-c));
    long double y = a - x0;
    return x + I * y;
}

static long double complex bedhead(long double complex z, double *p) {
    long double a = p[0];
    long double b = p[1];
    long double x0 = creall(z);
    long double y0 = cimagl(z);
    long double x = sinl(x0*y0/b)*y0+cosl(a*x0-y0);
    long double y = x0+sinl(y0)/b;
    return x + I * y;
}

static void genutil_hopa(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_attr(params,in,hopalong);
}

static void genutil_bedhead(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_attr(params,in,bedhead);
}

static double complex gexp(double complex z, double *p) {
    return cexp(z/cpow(p[0] + I * p[1],4));
}

static void genutil_exp(GenutilFuncParams *params, int n, AVFrame *in) {
    genutil_attr(params,in,gexp);
}
/***** lindenmeyer ***********/

static char *parseForm(GenutilFuncParams * params) {
    int iter = params->p[0];
    int k,i,j;
    GenutilGrowString *ss = av_genutil_growstring_new();
    av_genutil_growstring_append(ss,params->s[0]);
    for(k=0;k<iter;k++) {
        char *old = strdup(ss->str);
        av_genutil_growstring_reset(ss);
        size_t len = strlen(old);
        for(i=0;i<len;i++) {
            int found = 0;
            j=1;
            while(params->s[j]) {
                char *rule = params->s[j];
                j++;
                if(rule[0] == old[i]) {
                    found = 1;
                    av_genutil_growstring_append(ss,&rule[2]);
                }
            }
            if(!found) {
                char b[2];
                b[0] = old[i];
                b[1] = 0;
                av_genutil_growstring_append(ss,b);
            }
        }
        free(old);
    }
    char *ret = strdup(ss->str);
    av_genutil_growstring_free(ss);
    return ret;
}

static int genutil_count_brackets(char *str) {
    int k;
    int len = strlen(str);
    int max = 0;
    int current = 0;
    for(k=0;k<len;k++) {
       if(str[k] == '(') {
           current++;
           if(current > max) max = current;
       }
       if(str[k] == ')') current--;

    }
    return max;
}

static void genutil_interpreter(GenutilFuncParams *params,char *str, int b, double end, int cf, AVFrame *in) {
    double alpha = params->p[1];
    
    double *x = calloc(b*2,sizeof(double));
    double *y = calloc(b*2,sizeof(double));
    double *phi = calloc(b*2,sizeof(double));
    double r[4];
    x[0] = 0;
    y[0] = 0;
    phi[0] = 0;
    int k;
    int len = strlen(str);
    int br = 0;
    int col = 0;
    int count=0;
    double t = 0;
    for(k=0;k<len;k++) {
        switch(str[k]) {
            case '(':
            case '[':
                br++;
                x[br] = x[br-1];
                y[br] = y[br-1];
                phi[br] = phi[br-1];
                if(params->p[5] != 0)
                params->form = 7-br;
                break;
            case ')': br--;break;
            case '+': phi[br]+= (alpha*M_PI/180);break;
            case '-': phi[br]-= (alpha*M_PI/180);break;
            case '0': col = 0; break;
            case '1': col = 1; break;
            case '2': col = 2; break;
            case '3': col = 3; break;
            case '4': col = 4; break;
            case '5': col = 5; break;
            case '6': col = 6; break;
            case '7': col = 7; break;
            case '8': col = 8; break;
            case '9': col = 9; break;
            case 'a': params->form = 0; break;
            case 'b': params->form = 1; break;
            case 'c': params->form = 2; break;
            case 'd': params->form = 3; break;
            case 'e': params->form = 4; break;
            case 'f': params->form = 5; break;
            case 'g': params->form = 6; break;
            case 'F': 
              r[0] = x[br];
              r[1] = y[br];
              x[br] += sin(phi[br])*params->p[2];
              y[br] += cos(phi[br])*params->p[2];
              r[2] = x[br];
              r[3] = y[br];
              if(cf) {
                c_line_c(params,r,&t,params->delta,in);
              } else {
                c_line(params,r,params->colors[col],in);
              }
              count++;
              break;
            case 'G':
              r[0] = x[br];
              r[1] = y[br];
              x[br] += sin(phi[br])*params->p[3];
              y[br] += cos(phi[br])*params->p[3];
              r[2] = x[br];
              r[3] = y[br];
              if(cf) {
                c_line_c(params,r,&t,params->delta,in);
              } else {
                c_line(params,r,params->colors[col],in);
              }
              count++;
              break;
            default:
              break;
        }
        if(end>=0 && count>end) break;
    }
    free(x);
    free(y);
    free(phi);
}

static void genutil_lindenmeyer_pal(GenutilFuncParams *params, int n, AVFrame *in) {
    char *str = parseForm(params);
    av_log(params->ctx, AV_LOG_DEBUG, "str: %s\n", str);
    int max_brackets = genutil_count_brackets(str);
    genutil_interpreter(params,str,max_brackets,params->p[4],0,in);
    free(str);
}

static const char* pats[] = {"","-","--","---","+","++","+++","++++"};
static const int patlen[] = {0,1,2,3,1,2,3,4};

static const char* pats1[] = {"","-","--","---","----","-----","------","-------","+","++","+++","++++","+++++","++++++","+++++++","++++++++"};
static const int patlen1[] = {0,1,2,3,4,5,6,7,1,2,3,4,5,6,7,8};

static int rndf(int k, double *p) {
    return rand();
}

//static int testf(int k, double *p) {
//    return k*p[0];
//}

static int testsin(int k, double *p) {
    return ((sin(k*p[0])+1)/2) * 8;
}

static int testsin2(int k, double *p) {
    return ((sin(k*p[0])+1)/2) * 16;
}

static char *genrnd(int len, int (*f)(int,double*),double *p) {
    char *ret = calloc(len*6,sizeof(char));
    int pos = 0;
    ret[pos++] = 'F';
    for(int k=0;k<len;k++) {
        int rnd = f(k,p)%8;
        for(int j=0;j<patlen[rnd];j++) {
            ret[pos++]=pats[rnd][j];
        }
        ret[pos++]= rand()%8+48;;
        ret[pos++]='F';
    }
    return ret;
}

static char *genrnd2(int len, int (*f)(int,double*),double *p) {
    char *ret = calloc(len*12,sizeof(char));
    int pos = 0;
    ret[pos++] = 'F';
    for(int k=0;k<len;k++) {
        int rnd = f(k,p)%16;
        for(int j=0;j<patlen1[rnd];j++) {
            ret[pos++]=pats1[rnd][j];
        }
        ret[pos++]= rand()%8+48;;
        ret[pos++]='F';
    }
    return ret;
}

static void genutil_interpreter_cfunc(GenutilFuncParams *params,char *str, double end, AVFrame *in) {
    double alpha = params->p[1];
    
    double r[4];
    double x = 0;
    double y = 0;
    double phi = 0;
    int k;
    int len = strlen(str);
    int count=0;
    double t = 0;
    for(k=0;k<len;k++) {
        switch(str[k]) {
            case '+': phi+= (alpha*M_PI/180);break;
            case '-': phi-= (alpha*M_PI/180);break;
            case 'F': 
              r[0] = x;
              r[1] = y;
              x += sin(phi)*params->p[2];
              y += cos(phi)*params->p[2];
              r[2] = x;
              r[3] = y;
              c_line_c(params,r,&t,params->delta,in);
              count++;
              break;
            case 'G':
              r[0] = x;
              r[1] = y;
              x += sin(phi)*params->p[3];
              y += cos(phi)*params->p[3];
              r[2] = x;
              r[3] = y;
              c_line_c(params,r,&t,params->delta,in);
              count++;
              break;
            default:
              break;
        }
        if(end>=0 && count>end) break;
    }
}

static void genutil_lindenmeyer_t1(GenutilFuncParams *params, int n, AVFrame *in) {
    char *str = genrnd(floor(params->p[0]),testsin,&params->p[4]);
    av_log(params->ctx, AV_LOG_DEBUG, "str: %s\n", str);
    genutil_interpreter_cfunc(params,str,-1,in);
    free(str);
}

static void genutil_lindenmeyer_t2(GenutilFuncParams *params, int n, AVFrame *in) {
    char *str = genrnd2(floor(params->p[0]),testsin2,&params->p[4]);
    av_log(params->ctx, AV_LOG_DEBUG, "str: %s\n", str);
    genutil_interpreter_cfunc(params,str,-1,in);
    free(str);
}

static void genutil_lindenmeyer_rnd(GenutilFuncParams *params, int n, AVFrame *in) {
    srand(floor(fabs(params->p[4])));
    char *str = genrnd(floor(params->p[0]),rndf,NULL);
    av_log(params->ctx, AV_LOG_DEBUG, "str: %s\n", str);
    genutil_interpreter_cfunc(params,str,-1,in);
    free(str);
}

static void genutil_lindenmeyer_cfunc(GenutilFuncParams *params, int n, AVFrame *in) {
    char *str = parseForm(params);
    av_log(params->ctx, AV_LOG_DEBUG, "str: %s\n", str);
    int max_brackets = genutil_count_brackets(str);
    genutil_interpreter(params,str,max_brackets,params->p[4],1,in);
    free(str);
}

static GenutilFunc gfuncs[] = {
    {"hypo",genutil_hypo},
    {"lissg",genutil_lissg},
    {"lissq",genutil_lissq},
    {"liss",genutil_liss},
    {"leg",genutil_leg},
    {"cb",genutil_cb},
    {"sinoid",genutil_sinoid},
    {"capri",genutil_capri},
    {"scara",genutil_scara},
    {"tro",genutil_tro},
    {"circ",genutil_circ},
    {"rhodo",genutil_rhodo},
    {"superrose",genutil_super_rose},
    {"superspiral",genutil_super_spiral},
    {"epispiral",genutil_epi_spiral},
    {"spiral",genutil_spiral},
    {"gielis",genutil_gielis},
    {"maclaurin",genutil_maclaurin},
    {"harmonic",genutil_harmonic},
    {"cornoid",genutil_cornoid},
    {"line",genutil_line},
    {"poly",genutil_poly},
    {"spline",genutil_spline},
    {"polymul",genutil_poly_mul},
    {"polymod",genutil_poly_mod},
    {"polymod1",genutil_poly_mod1},
    {"polysin",genutil_poly_sin},
    {"polyq",genutil_poly_q},
    {"polyrnd",genutil_poly_rnd},
    {"polydsink",genutil_poly_dsink},
    {"polydsin",genutil_poly_dsin},
    {"polydrnd",genutil_poly_drnd},
    {"polydq",genutil_poly_dq},
    {"polydqs",genutil_poly_dqs},
    {"ulam",genutil_poly_ulam},
    {"frac",genutil_frac},
    {"mira",genutil_mira},
    {"kaneko1",genutil_kaneko1},
    {"kaneko2",genutil_kaneko2},
    {"ginger",genutil_ginger},
    {"hopa",genutil_hopa},
    {"bedhead",genutil_bedhead},
    {"exp",genutil_exp},
    {"hilbert",genutil_hilbert},
    {"lm",genutil_lindenmeyer_pal},
    {"lmc",genutil_lindenmeyer_cfunc},
    {"lmrnd",genutil_lindenmeyer_rnd},
    {"lmtest",genutil_lindenmeyer_t1},
    {"lmtest2",genutil_lindenmeyer_t2},
    {"addoid",genutil_addoid},
    {"phypo",genutil_phypo},
    {"pkardioids",genutil_pkardioids},
    {"plemg",genutil_plemg},
    {"plemb",genutil_plemb},
    {"pepi",genutil_pepi},
    {"pnodal",genutil_pnodal},
    {"ptalbot",genutil_ptalbot},
    {"pfolium",genutil_pfolium},
    {"plissg",genutil_plissg},
    {"plissq",genutil_plissq},
    {"plissp",genutil_plissp},
    {"pliss",genutil_pliss},
    {"pleg",genutil_pleg},
    {"pcb",genutil_pcb},
    {"psinoid",genutil_psinoid},
    {"pcapri",genutil_pcapri},
    {"pscara",genutil_pscara},
    {"ptro",genutil_ptro},
    {"pcirc",genutil_pcirco},
    {"pcircs",genutil_circs},
    {"prhodo",genutil_prhodo},
    {"psuperrose",genutil_psuper_rose},
    {"psuperspiral",genutil_psuper_spiral},
    {"pepispiral",genutil_pepi_spiral},
    {"pspiral",genutil_pspiral},
    {"pgielis",genutil_pgielis},
    {"pmaclaurin",genutil_pmaclaurin},
    {"pharmonic",genutil_pharmonic},
    {"pcornoid",genutil_pcornoid},
    {"pline",genutil_pline},
    {"ppoly",genutil_ppoly},
    {"ppolyu",genutil_ppolyu},
    {"prndpoly",genutil_prndpoly},
    {"paddoid",genutil_paddoid},
    {"philbert",genutil_philbert},
    {"zero",gzero},
    {NULL,NULL}
};

void (*av_genutil_get_ffunc(const char *name))(GenutilFuncParams*,int,AVFrame*) {
    int k=0;
    while(gfuncs[k].name) {
        if(!strcmp(name, gfuncs[k].name)) {
            return gfuncs[k].f;
        }
        k++;
    }
    return NULL;

}

void av_genutil_parse_ffunc(const char *nf, double *p, void (**f)(GenutilFuncParams*,int,AVFrame*)) {
    char *saveptr,*token;
    char *str = strdup(nf);
    int j;
    for (j = 0; ; j++, str = NULL) {
        token = av_strtok(str, " |", &saveptr);
        if (token == NULL)
            break;
        if(j==0) {
            *f = av_genutil_get_ffunc(token);
        } else {
            p[j-1] = atof(token);
        }
    }
    free(str);
}

/* ============================= NFUNCS *******************************************/

/*** spline ***/
static double p3(int n, double a,double b,double c,double d) {
    return a + b*n + c*n*n + d*n*n*n;
}

static double pcubic(int n, double *p) {
    int len = p[0];
    int end = p[1];
    int k,index;
    for(k=1;k<len;k++) {
        if(n<p[k*5+2]) {
            index = (k-1)*5 +2;
            return p3(n-p[index],p[index+1],p[index+2],p[index+3],p[index+4]);
        } 
    }
    if(n<end) {
        index = (len-1)*5+2;
        return p3(n-p[index],p[index+1],p[index+2],p[index+3],p[index+4]);
    } else {
        index = (len-1)*5+2;
        return p3(end-p[index],p[index+1],p[index+2],p[index+3],p[index+4]);
    }
}

static double pspline(int n, double *p) {
    int len = p[0];
    double *out = calloc(len*5+2, sizeof(double));
    av_genutil_spline(p,out);
    double ret =  pcubic(n,out);
    free(out);
    return ret;
}

static double npoly(int n, double *p) {
    double x = (double)n * p[0];
    return p[1] + p[2]*x + p[3]*x+x + p[4]*x*x*x + p[5]*x*x*x*x;
}

static double nsin(int n, double *p) {
    return (sin(p[0]*n - M_PI/2 + p[3]))*p[1]+p[2]; 
}

static double rsin(int n, double *p) {
    return (sin(p[0]*n + p[3]))*(p[1]==0.0?1:p[1])+p[2]; 
}

static double rcos(int n, double *p) {
    return (cos(p[0]*n + p[3]))*(p[1]==0.0?1:p[1])+p[2]; 
}

static double linsin(int n, double *p) {
    return p[0]*n+(sin(p[1]*n + p[3]))*p[2]; 
}

static double idn(int n, double *p) {
    return n;
}

static double nexp(int n, double *p) {
    return p[2] + p[1]*(exp(p[0]*n-p[3]));
}

static double nexpinv(int n, double *p) {
    return p[2] + p[1]/(exp(p[0]*n-p[3]));
}

static double nexp2(int n, double *p) {
    return p[2] + p[1]*(exp(p[0]*n*n-p[3]));
}

static double constant(int n, double *p) {
    return p[0];
}

static double sqn(int n, double *p) {
    double x = (double)n*p[0];
    return p[1] * (sqrt(x+0.25)-0.5);
}

static double ln(int n, double *p) {
    double x = (double)n*p[0];
    return p[1] * log(x+1);
}

static double lin(int n, double *p) {
    return p[0] + p[1]*n;
}

static double inv(int n, double *p) {
    double x = (double)n*p[0];
    return p[1]*(1/(p[2]*((double)x+1)));
}

static double pval(int n, double *p) {
    int k;
    if(n<p[0]) return 0;
    if(p[NMAXPARAMS-2]&&(n>p[NMAXPARAMS-2])) return p[NMAXPARAMS-1];
    for(k=0;k<NMAXPARAMS-2;k+=2) {
        if(k>0 && p[k]==0) return p[k-1]; 
        if(n>p[k] && n<=p[k+2])
            return p[k+1];
    }
    return 0;
}

static double pchain(int n, double *p) {
    int len = p[0];
    int section = n/len + 1;
    if(section > NMAXPARAMS) return 0;
    double x = (double)(n%len)/(double)len;
    return p[section] + (p[section+1]-p[section])*x;
}

static double gauss(int n, double *p) {
    double m = p[1]==0?1:p[1];
    double fac = p[2]==0?0.5:p[2];
    double x = (double)n*p[0];
    return m * exp(-fac * pow(x-p[3],2)); 
}

typedef struct NFunc {
    const char *name;
    double (*f)(int n, double *p);
} NFunc;

static NFunc nfuncs[] = {
    {"idn",idn},
    {"constant",constant},
    {"lin",lin},
    {"linsin",linsin},
    {"poly",npoly},
    {"pchain",pchain},
    {"pcubic",pcubic},
    {"pspline",pspline},
    {"pval",pval},
    {"sin",nsin},
    {"rsin",rsin},
    {"rcos",rcos},
    {"inv",inv},
    {"gauss",gauss},
    {"sq",sqn},
    {"ln",ln},
    {"exp",nexp},
    {"expinv",nexpinv},
    {"exp2",nexp2},
    {NULL,NULL}
};

double (*av_genutil_get_nfunc(const char *name))(int,double*) {
    int k=0;
    while(nfuncs[k].name) {
        if(!strcmp(name, nfuncs[k].name)) {
            return nfuncs[k].f;
        }
        k++;
    }
    return NULL;
}

void av_genutil_parse_nfunc(const char *nf, double *p, double (**f)(int,double*)) {
    char *saveptr,*token;
    char *str = strdup(nf);
    int j;
    for (j = 0; ; j++, str = NULL) {
        token = av_strtok(str, " |", &saveptr);
        if (token == NULL)
            break;
        if(j==0) {
            *f = av_genutil_get_nfunc(token);
        } else {
            p[j-1] = atof(token);
        }
    }
    free(str);
}


/* ============================= CFUNCS *******************************************/

static double cfsin(double t, double *p) {
    return (sin(p[0]*t+p[3])+1)*0.5*p[1]+p[2]; 
}

static double cfcos(double t, double *p) {
    return (cos(p[0]*t+p[3])+1)*0.5*p[1]+p[2]; 
}

static double cconstant(double t, double *p) {
    return p[0];
}

static double clin(double t, double *p) {
    return p[0] + p[1]*t;
}
static double cpchain(double t, double *p) {
    int len = p[0];
    int section = t/len + 1;
    if(section > CMAXPARAMS) return 0;
    double x = (double)(t-section*len)/(double)len;
    return p[section] + (p[section+1]-p[section])*x;
}

double czero(double t, double *p) {
    return 0;
}

static CFunc cfuncs[] = {
    {"sin",cfsin},
    {"cos",cfcos},
    {"lin",clin},
    {"constant",cconstant},
    {"pchain",cpchain},
    {NULL,NULL}
};

double (*av_genutil_get_cfunc(const char *name))(double,double*) {
    int k=0;
    while(cfuncs[k].name) {
        if(!strcmp(name, cfuncs[k].name)) {
            return cfuncs[k].f;
        }
        k++;
    }
    return NULL;
}

void av_genutil_parse_cfunc(const char *nf, double *p, double (**f)(double,double*)) {
    char *saveptr,*token;
    char *str = strdup(nf);
    int j;
    for (j = 0; ; j++, str = NULL) {
        token = av_strtok(str, " |", &saveptr);
        if (token == NULL)
            break;
        if(j==0) {
            *f = av_genutil_get_cfunc(token);
        } else {
            p[j-1] = atof(token);
        }
    }
    free(str);
}



static void get_hsv2rgb_params(double phi, double S, double V, uint8_t *h, double *p, double *q, double *t) {
    double H = phi<0?((180/M_PI)*phi)+360:(180/M_PI)*phi;
    *h = floor(H/60.0);
    double f = H/60.0 - *h;
    *p = V*(1-S);
    *q = V*(1-S*f);
    *t = V*(1-S*(1-f)); 
}

static double hsv_r(double phi, double S, double V) {
    double p,q,t;
    uint8_t h;
    get_hsv2rgb_params(phi, S, V, &h, &p, &q, &t);
    double retval = 0;
    switch(h) {
        case 0:
        case 5: retval = V; break;
        case 1: retval = q; break;
        case 2:
        case 3: retval = p; break;
        case 4: retval = t; break;        
    }
    return retval;
}

static double hsv_g(double phi, double S, double V) {
    double p,q,t;
    uint8_t h;
    get_hsv2rgb_params(phi, S, V, &h, &p, &q, &t);
    double retval = 0;
    switch(h) {
        case 1:
        case 2: retval = V; break;
        case 3: retval = q; break;
        case 4:
        case 5: retval = p; break;
        case 0: retval = t; break;        
    }
    return retval;
}

static double hsv_b(double phi, double S, double V) {
    double p,q,t;
    uint8_t h;
    get_hsv2rgb_params(phi, S, V, &h, &p, &q, &t);
    double retval = 0;
    switch(h) {
        case 3:
        case 4: retval = V; break;
        case 5: retval = q; break;
        case 0:
        case 1: retval = p; break;
        case 2: retval = t; break;        
    }
    return retval;
}

void av_genutil_get_color(double (*cfunc[3])(double,double*), double cp[3][10], double t, int cmod, int is_rgb, double *out) {
    double a,b,c;
    double r0,g0,b0;
    a = cfunc[0](t,cp[0]);
    b = cfunc[1](t,cp[1]);
    c = cfunc[2](t,cp[2]);
    switch(cmod) {
        case 0: // YUV
           if(is_rgb) {
                out[0] = a + 1.403 * (c-0.5);
                out[0] = out[0]>1?1:out[0];
                out[0] = out[0]<0?0:out[0];
                out[1] = a - 0.3455 * (b - 0.5) - (0.7169 * (c-0.5));
                out[1] = out[1]>1?1:out[1];
                out[1] = out[1]<0?0:out[1];
                out[2] = a + 1.7790 * (b - 0.5);
                out[2] = out[2]>1?1:out[2];
                out[2] = out[2]<0?0:out[2];
           } else {
               out[0] = a;
               out[1] = b;
               out[2] = c;
           }
           break;
        case 1: // RGB
           if(is_rgb) {
               out[0] = b;
               out[1] = c;
               out[2] = a;
           } else {
                out[0] = a * 0.299000  + b * 0.587000  + c * 0.114000;
                out[1] = a * -0.168736 + b * -0.331264 + c * 0.500000  + 0.5;
                out[2] = a * 0.500000  + b * -.418688  + c * -0.081312 + 0.5;  
           }
           break;
        case 2: // HSV
           if(a>=M_PI*2) a -= M_PI*2;
            r0 = hsv_r(a,b,c);
            g0 = hsv_g(a,b,c);
            b0 = hsv_b(a,b,c);
            if(is_rgb) {
               out[0] = g0;
               out[1] = b0;
               out[2] = r0;
            } else {
                out[0] = r0 * 0.299000  + g0 * 0.587000  + b0 * 0.114000;
                out[1] = r0 * -0.168736 + g0 * -0.331264 + b0 * 0.500000  + 0.5;
                out[2] = r0 * 0.500000  + g0 * -.418688  + b0 * -0.081312 + 0.5;  
            }

                
    }
}

double av_genutil_avg(double *p, int len) {
    double ret = 0;
    int k;
    for(k=0;k<len;k++) {
        ret += p[k];
    }
    return ret/(double)len;
}

void av_genutil_replace_char(char *str, char needle, char rep) {
    int k=0;
    while(str[k]) {
        if(str[k] == needle) str[k] = rep;
        k++;
    }
}

GenutilGrowString *av_genutil_growstring_new() {
    GenutilGrowString *ss = malloc(sizeof(GenutilGrowString));
    if(NULL == ss) goto fail;
    ss->str = malloc(SMARTY_SIZE_INIT);
    if(NULL == ss->str ) goto fail;
    ss->end = ss->str;
    *ss->end = '\0';
    ss->size = SMARTY_SIZE_INIT;
    return ss; 

    fail:
       free(ss);
       return NULL;
}

void av_genutil_growstring_append(GenutilGrowString *ss, char *suffix) {
    while(*suffix != '\0'){
      if(! ((ss->end - ss->str) < (ss->size-1)) ){
        size_t offset = ss->end - ss->str;
        char* newstr = realloc(ss->str, (ss->size * 2) );
        if(NULL == newstr){ ss->end--; break; }
        ss->str = newstr; 
        ss->size = ss->size * 2;
        ss->end = ss->str + offset; 
      }
      *(ss->end) = *suffix;
      suffix++;
      ss->end++;
    }
    *(ss->end) = '\0';
}

void av_genutil_growstring_reset(GenutilGrowString *ss) {
    ss->end = ss->str;
    *ss->end = '\0';
}

void av_genutil_growstring_free(GenutilGrowString *ss) {
    free(ss->str);
    free(ss);
}

