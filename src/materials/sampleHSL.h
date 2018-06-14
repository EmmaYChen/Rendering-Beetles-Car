#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MATERIALS_SAMPLEHSL_H
#define PBRT_MATERIALS_SAMPLEHSL_H

#include "pbrt.h"
#include "geometry.h"

namespace pbrt {

/**
 * Converts an HSL color value to RGB. Conversion formula
 * adapted from http://en.wikipedia.org/wiki/HSL_color_space.
 * Assumes h is in the set [0,360], s, and l are contained in the set [0, 1] and
 * returns r, g, and b in the set [0, 1].
 *
 * @param   {number}  h       The hue
 * @param   {number}  s       The saturation
 * @param   {number}  l       The lightness
 * @return  {Array}           The RGB representation
 */
 class RGB{
 public:
    unsigned char R;
    unsigned char G;
    unsigned char B;
    RGB(unsigned char r, unsigned char g, unsigned char b){
        R = r;
        G = g;
        B = b;}
    bool Equals(RGB rgb){
        return (R == rgb.R) && (G == rgb.G) && (B == rgb.B);}
 };

 class HSL{
 public:
    int H;
    float S;
    float L;
    HSL(int h, float s, float l){
        H = h;
        S = s;
        L = l;}
    bool Equals(HSL hsl){
        return (H == hsl.H) && (S == hsl.S) && (L == hsl.L);}
 };

float HueToRGB(float v1, float v2, float vH) {
    if (vH < 0)
        vH += 1;
    if (vH > 1)
        vH -= 1;
    if ((6 * vH) < 1)
        return (v1 + (v2 - v1) * 6 * vH);
    if ((2 * vH) < 1)
        return v2;
    if ((3 * vH) < 2)
        return (v1 + (v2 - v1) * ((2.0f / 3) - vH) * 6);
    return v1;
 }

// Transform hsl to rgb
RGB HSLToRGB(HSL hsl) {
    unsigned char r = 0;
    unsigned char g = 0;
    unsigned char b = 0;
    if (hsl.S == 0)
    {
        r = g = b = (unsigned char)(hsl.L);
    }
    else
    {
        float v1, v2;
        float hue = (float)hsl.H / 360;
        v2 = (hsl.L < 0.5) ? (hsl.L * (1 + hsl.S)) : ((hsl.L + hsl.S) - (hsl.L * hsl.S));
        v1 = 2 * hsl.L - v2;
        r = (unsigned char)(HueToRGB(v1, v2, hue + (1.0f / 3)));
        g = (unsigned char)(HueToRGB(v1, v2, hue));
        b = (unsigned char)(HueToRGB(v1, v2, hue - (1.0f / 3)));
    }
    return RGB(r, g, b);
}

// sample from the color table
Point3f SampleFromHSL(){
    float s = 1;
    float l = 0.5;
    int H = (int)(6*(rand()% 60));
    // If we want to decrease flakes of one certain color
    // if (H<=30&& H>=0) return Point3f(1,1,1);
    RGB color_rgb = HSLToRGB(HSL(H,s,l));
    return Point3f(color_rgb.R,color_rgb.G,color_rgb.B);
}

// return the average color of the table
Point3f AverageColor(){
    Point3f average = Point3f(0,0,0);
    float s = 1;
    float l = 0.5;
    for (int i = 0;i<60;i++){
        int H = 6*i;
        RGB color_rgb = HSLToRGB(HSL(H,s,l));
        average += Point3f(color_rgb.R,color_rgb.G,color_rgb.B);
    }
    average = average/60;
    return average;}

}
#endif
