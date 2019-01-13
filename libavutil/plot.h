#ifndef AVUTIL_PLOT_H
#define AVUTIL_PLOT_H

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

void av_plot_number(int x, int y, double z, int w, int h, int linesize, uint8_t *buf);
void av_plot_int(int x, int y, int n, int w, int h, int linesize, uint8_t *buf);
void av_plot_form(int form, int x, int y, double a, int w, int h,int linesize, uint8_t *buf);
void av_plot_form_p(int form, int x, int y, double a, int w, int h,int linesize, uint8_t *buf);
void av_plot_number_p(int x, int y, double z, uint8_t *color, int w, int h, int linesize, uint8_t *buf);
void av_plot_int_p(int x, int y, int n, uint8_t *color, int w, int h, int linesize, uint8_t *buf);

#endif
