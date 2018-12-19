#include "genutil.h"



typedef struct NFunc {
    const char *name;
    double (*f)(int n, double *p);
} NFunc;

typedef struct CFunc {
    const char *name;
    double (*f)(double t, double *p);
} CFunc;

double complex av_genutil_rotate(double complex z, double phi) {
    return phi!=0.0?(cos(phi) + sin(phi)*I) * z:z;
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
    double x = cos(t) + a * sin(b*t) * sin(b*t);
    double y = sin(t) + c * sin(d*t) * sin(d*t);
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
    double x = (a*cos(2*t) + b*cos(t))*cos(t);
    double y = (a*cos(2*t) + b*cos(t))*sin(t);
    return x + I * y;
}


/* ============================= NFUNCS *******************************************/

static double npoly(int n, double *p) {
    double x = (double)n * p[0];
    return p[1] + p[2]*x + p[3]*x+x + p[4]*x*x*x + p[5]*x*x*x*x;
}

static double nsin(int n, double *p) {
    return (sin(p[0]*n - M_PI/2 + p[3])+1)*p[1]+p[2]; 
}

static double idn(int n, double *p) {
    return n;
}

static double constant(int n, double *p) {
    return p[0]?p[0]:1;
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
    return p[1]*(1/((double)x+1));
}

static double pchain(int n, double *p) {
    int len = p[0];
    int section = n/len + 1;
    if((section+1) >= 20) return 0;
    double x = (double)(n%len)/(double)len;
    return p[section] + (p[section+1]-p[section])*x;
}

static double gauss(int n, double *p) {
    double m = p[1]==0?1:p[1];
    double fac = p[2]==0?0.5:p[2];
    double x = (double)n*p[0];
    return m * exp(-fac * pow(x-p[3],2)); 
}

static NFunc nfuncs[] = {
    {"idn",idn},
    {"constant",constant},
    {"lin",lin},
    {"poly",npoly},
    {"pchain",pchain},
    {"sin",nsin},
    {"inv",inv},
    {"gauss",gauss},
    {"sq",sqn},
    {"ln",ln},
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

/* ============================= CFUNCS *******************************************/

static double cfsin(double t, double *p) {
    return (sin(p[0]*t)+1)*0.5*p[1]+p[2]; 
}

static double cfcos(double t, double *p) {
    return (cos(p[0]*t)+1)*0.5*p[1]+p[2]; 
}

static double cconstant(double t, double *p) {
    return p[0];
}

static double cpchain(double t, double *p) {
    int len = p[0];
    int section = t/len + 1;
    if((section+1) >= 10) return 0;
    double x = (double)(t-section*len)/(double)len;
    return p[section] + (p[section+1]-p[section])*x;
}

double czero(double t, double *p) {
    return 0;
}

static CFunc cfuncs[] = {
    {"sin",cfsin},
    {"cos",cfcos},
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

