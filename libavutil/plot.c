#include "plot.h"

static void av_plot_pix(int x, int y, double a, int w, int h,int linesize, uint8_t *buf) {
    if(x>=0 && x<w && y>=0 && y<h) {
        buf[x + linesize * y] = floor(a*255);
    }
}

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

void av_plot_form(int form, int x, int y, double a, int w, int h,int linesize, uint8_t *buf) {
    switch(form) {
        case 0: av_plot_pix(x,y,a,w,h,linesize,buf);
                break;
        case 1: av_plot_pix(x,y,a,w,h,linesize,buf);
                av_plot_pix(x+1,y,a,w,h,linesize,buf);
                av_plot_pix(x,y+1,a,w,h,linesize,buf);
                av_plot_pix(x+1,y+1,a,w,h,linesize,buf);
                break;
        case 2: av_plot_pix(x,y,a,w,h,linesize,buf);
                av_plot_pix(x+1,y,a,w,h,linesize,buf);
                break;
        case 3: av_plot_pix(x,y,a,w,h,linesize,buf);
                av_plot_pix(x,y+1,a,w,h,linesize,buf);
                break;
        case 4: av_plot_pix(x,y,a,w,h,linesize,buf);
                av_plot_pix(x-1,y,a,w,h,linesize,buf);
                av_plot_pix(x+1,y,a,w,h,linesize,buf);
                av_plot_pix(x,y-1,a,w,h,linesize,buf);
                av_plot_pix(x,y+1,a,w,h,linesize,buf);
                break;
        case 5: av_plot_pix(x,y,a,w,h,linesize,buf);
                av_plot_pix(x-1,y,a,w,h,linesize,buf);
                av_plot_pix(x+1,y,a,w,h,linesize,buf);
                av_plot_pix(x,y-1,a,w,h,linesize,buf);
                av_plot_pix(x,y+1,a,w,h,linesize,buf);
                av_plot_pix(x+1,y+1,a,w,h,linesize,buf);
                av_plot_pix(x-1,y-1,a,w,h,linesize,buf);
                break;
        case 6: av_plot_pix(x,y,a,w,h,linesize,buf);
                av_plot_pix(x+1,y,a,w,h,linesize,buf);
                av_plot_pix(x-1,y,a,w,h,linesize,buf);
                av_plot_pix(x+2,y,a,w,h,linesize,buf);
                av_plot_pix(x-2,y,a,w,h,linesize,buf);
                break;
        case 7: av_plot_pix(x,y,a,w,h,linesize,buf);
                av_plot_pix(x+1,y,a,w,h,linesize,buf);
                av_plot_pix(x-1,y,a,w,h,linesize,buf);
                av_plot_pix(x+2,y,a,w,h,linesize,buf);
                av_plot_pix(x-2,y,a,w,h,linesize,buf);
                av_plot_pix(x,y-1,a,w,h,linesize,buf);
                av_plot_pix(x+1,y-1,a,w,h,linesize,buf);
                av_plot_pix(x-1,y-1,a,w,h,linesize,buf);
                av_plot_pix(x,y-2,a,w,h,linesize,buf);
                av_plot_pix(x,y+1,a,w,h,linesize,buf);
                av_plot_pix(x+1,y+1,a,w,h,linesize,buf);
                av_plot_pix(x-1,y+1,a,w,h,linesize,buf);
                av_plot_pix(x,y+2,a,w,h,linesize,buf);
//  00X00
//  0XXX0
//  XXXXX
//  0XXX0
//  00X00
                break;
        case 8: av_plot_pix(x,y,a,w,h,linesize,buf);
                av_plot_pix(x+4,y,a,w,h,linesize,buf);
                av_plot_pix(x-4,y,a,w,h,linesize,buf);
                av_plot_pix(x+8,y,a,w,h,linesize,buf);
                av_plot_pix(x-8,y,a,w,h,linesize,buf);
                av_plot_pix(x+12,y,a,w,h,linesize,buf);
                av_plot_pix(x-12,y,a,w,h,linesize,buf);
                break;
        default:av_plot_pix(x,y,a,w,h,linesize,buf);
                break;
    }
}


