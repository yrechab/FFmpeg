#include "plot.h"

static void av_plot_pix(int x, int y, double a, int w, int h,int linesize, uint8_t *buf) {
    if(x>=0 && x<w && y>=0 && y<h) {
        buf[x + linesize * y] = floor(a*255);
    }
}
static void av_plot_pix_p(int x, int y, double a, int w, int h,int linesize, uint8_t *buf) {
    if(x>=0 && x<w && y>=0 && y<h) {
        buf[x*4 + linesize * y] = floor(a*255);
    }
}
static void av_plot_pix_pc(int x, int y, uint8_t *color, int w, int h,int linesize, uint8_t *buf) {
    int k;
    if(x>=0 && x<w && y>=0 && y<h) {
        for(k=0;k<4;k++) {
            buf[x*4+k + linesize * y] = color[k];
        }
    }
}
#define PF(x) void (*x)(int,int,double,int,int,int,uint8_t*)

static uint8_t numbers[12][8] = {
   {
       0b00000000,
       0b01111110,
       0b01100110,
       0b01100110,
       0b01100110,
       0b01100110,
       0b01100110,
       0b01111110,
   },
   {
       0b00000000,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
   },
   {
       0b00000000,
       0b01111110,
       0b00000110,
       0b00000110,
       0b01111110,
       0b01100000,
       0b01100000,
       0b01111110,
   },
   {
       0b00000000,
       0b01111110,
       0b00000110,
       0b00000110,
       0b01111110,
       0b00000110,
       0b00000110,
       0b01111110,
   },
   {
       0b00000000,
       0b01100110,
       0b01100110,
       0b01100110,
       0b01111110,
       0b00000110,
       0b00000110,
       0b00000110,
   },
   {
       0b00000000,
       0b01111110,
       0b01100000,
       0b01100000,
       0b01111110,
       0b00000110,
       0b00000110,
       0b01111110,
   },
   {
       0b00000000,
       0b01111110,
       0b01100000,
       0b01100000,
       0b01111110,
       0b01100110,
       0b01100110,
       0b01111110,
   },
   {
       0b00000000,
       0b01111110,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
       0b00000110,
   },
   {
       0b00000000,
       0b01111110,
       0b01100110,
       0b01100110,
       0b01111110,
       0b01100110,
       0b01100110,
       0b01111110,
   },
   {
       0b00000000,
       0b01111110,
       0b01100110,
       0b01100110,
       0b01111110,
       0b00000110,
       0b00000110,
       0b01111110,
   },
   {
       0b00000000,
       0b00000000,
       0b00000000,
       0b00000000,
       0b00000000,
       0b00000000,
       0b00111000,
       0b00111000,
   },
   {
       0b00000000,
       0b00000000,
       0b00000000,
       0b00000000,
       0b01111110,
       0b00000000,
       0b00000000,
       0b00000000,
   },

};

static void av_plot_ch(char ch, int x, int y, int w, int h, int linesize, uint8_t *buf) {
    uint8_t *bmp;
    if(ch == '.') {
        bmp = numbers[10];
    } else if(ch == '-') {
        bmp = numbers[11];
    } else if(ch>=48 && ch<58) {
        bmp = numbers[ch-48];
    } else {
        return;
    }
    int k,j;
    for(k=0;k<8;k++) {
        uint8_t b = bmp[k];
        for(j=0;j<8;j++) {
            if( (1 << 7-j)&b) {
                av_plot_pix(x+j,y+k,1,w,h,linesize,buf);
            }
        }
    }
}

void av_plot_number(int x, int y, double z, int w, int h, int linesize, uint8_t *buf) {
    char s[50];
    sprintf(s,"%lf",z);
    char *ptr = s;
    int k = 0;
    while(*ptr) {
        av_plot_ch(*ptr,x+k,y,w,h,linesize,buf);
        ptr++;
        k+=9;
    }
}

void av_plot_int(int x, int y, int n, int w, int h, int linesize, uint8_t *buf) {
    char s[50];
    sprintf(s,"%d",n);
    char *ptr = s;
    int k = 0;
    while(*ptr) {
        av_plot_ch(*ptr,x+k,y,w,h,linesize,buf);
        ptr++;
        k+=9;
    }
}


static void av_plot_ch_p(char ch, uint8_t *color, int x, int y, int w, int h, int linesize, uint8_t *buf) {
    uint8_t *bmp;
    if(ch == '.') {
        bmp = numbers[10];
    } else if(ch == '-') {
        bmp = numbers[11];
    } else if(ch>=48 && ch<58) {
        bmp = numbers[ch-48];
    } else {
        return;
    }
    int k,j;
    for(k=0;k<8;k++) {
        uint8_t b = bmp[k];
        for(j=0;j<8;j++) {
            if( (1 << 7-j)&b) {
                av_plot_pix_pc(x+j,y+k,color,w,h,linesize,buf);
            }
        }
    }
}

void av_plot_number_p(int x, int y, double z, uint8_t *color, int w, int h, int linesize, uint8_t *buf) {
    char s[50];
    sprintf(s,"%lf",z);
    char *ptr = s;
    int k = 0;
    while(*ptr) {
        av_plot_ch_p(*ptr,color,x+k,y,w,h,linesize,buf);
        ptr++;
        k+=9;
    }
}

void av_plot_int_p(int x, int y, int n, uint8_t *color, int w, int h, int linesize, uint8_t *buf) {
    char s[50];
    sprintf(s,"%d",n);
    char *ptr = s;
    int k = 0;
    while(*ptr) {
        av_plot_ch_p(*ptr,color,x+k,y,w,h,linesize,buf);
        ptr++;
        k+=9;
    }
}



static void plot_form(PF(plotf),int form, int x, int y, double a, int w, int h,int linesize, uint8_t *buf) {
    int i,j;
    switch(form) {
        case 0: plotf(x,y,a,w,h,linesize,buf);
                break;
        case 1: plotf(x,y,a,w,h,linesize,buf);
                plotf(x+1,y,a,w,h,linesize,buf);
                plotf(x,y+1,a,w,h,linesize,buf);
                plotf(x+1,y+1,a,w,h,linesize,buf);
                break;
        case 2: plotf(x,y,a,w,h,linesize,buf);
                plotf(x+1,y,a,w,h,linesize,buf);
                break;
        case 3: plotf(x,y,a,w,h,linesize,buf);
                plotf(x,y+1,a,w,h,linesize,buf);
                break;
        case 4: plotf(x,y,a,w,h,linesize,buf);
                plotf(x-1,y,a,w,h,linesize,buf);
                plotf(x+1,y,a,w,h,linesize,buf);
                plotf(x,y-1,a,w,h,linesize,buf);
                plotf(x,y+1,a,w,h,linesize,buf);
                break;
        case 5: plotf(x,y,a,w,h,linesize,buf);
                plotf(x-1,y,a,w,h,linesize,buf);
                plotf(x+1,y,a,w,h,linesize,buf);
                plotf(x,y-1,a,w,h,linesize,buf);
                plotf(x,y+1,a,w,h,linesize,buf);
                plotf(x+1,y+1,a,w,h,linesize,buf);
                plotf(x-1,y-1,a,w,h,linesize,buf);
                break;
        case 6: plotf(x,y,a,w,h,linesize,buf);
                plotf(x+1,y,a,w,h,linesize,buf);
                plotf(x-1,y,a,w,h,linesize,buf);
                plotf(x+2,y,a,w,h,linesize,buf);
                plotf(x-2,y,a,w,h,linesize,buf);
                break;
        case 7: plotf(x,y,a,w,h,linesize,buf);
                plotf(x+1,y,a,w,h,linesize,buf);
                plotf(x-1,y,a,w,h,linesize,buf);
                plotf(x+2,y,a,w,h,linesize,buf);
                plotf(x-2,y,a,w,h,linesize,buf);
                plotf(x,y-1,a,w,h,linesize,buf);
                plotf(x+1,y-1,a,w,h,linesize,buf);
                plotf(x-1,y-1,a,w,h,linesize,buf);
                plotf(x,y-2,a,w,h,linesize,buf);
                plotf(x,y+1,a,w,h,linesize,buf);
                plotf(x+1,y+1,a,w,h,linesize,buf);
                plotf(x-1,y+1,a,w,h,linesize,buf);
                plotf(x,y+2,a,w,h,linesize,buf);
//  00X00
//  0XXX0
//  XXXXX
//  0XXX0
//  00X00
                break;
        case 8: plotf(x,y,a,w,h,linesize,buf);
                plotf(x+4,y,a,w,h,linesize,buf);
                plotf(x-4,y,a,w,h,linesize,buf);
                plotf(x+8,y,a,w,h,linesize,buf);
                plotf(x-8,y,a,w,h,linesize,buf);
                plotf(x+12,y,a,w,h,linesize,buf);
                plotf(x-12,y,a,w,h,linesize,buf);
                break;
        case 9: for(i=-2;i<=2;i++)
                    for(j=-2;j<=2;j++) 
                        plotf(x+i,y+j,a,w,h,linesize,buf);
                break;
        default:plotf(x,y,a,w,h,linesize,buf);
                break;
    }
}

void av_plot_form(int form, int x, int y, double a, int w, int h,int linesize, uint8_t *buf) {
    plot_form(av_plot_pix,form,x,y,a,w,h,linesize,buf);
}

void av_plot_form_p(int form, int x, int y, double a, int w, int h,int linesize, uint8_t *buf) {
    plot_form(av_plot_pix_p,form,x,y,a,w,h,linesize,buf);
}
