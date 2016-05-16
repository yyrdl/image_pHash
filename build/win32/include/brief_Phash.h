 

#ifndef _PHASH_H
#define _PHASH_H

 
#include <limits.h>
#include <math.h>
#include "dirent.h"
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#define __STDC_CONSTANT_MACROS

 
#define cimg_debug 0
#define cimg_display 0
#include "CImg.h"
using namespace cimg_library;

#if !defined(__GLIBC__) && !defined(_WIN32)
#include <sys/param.h>
#include <sys/sysctl.h>
#endif

using namespace std;
  
 

 
typedef unsigned long long ulong64;
typedef signed long long long64;
typedef unsigned char uint8_t;
typedef unsigned int uint32_t;
 

#ifdef __cplusplus
extern "C" {
#endif  

static CImg<float>* ph_dct_matrix(const int N);

CImg<float>* ph_dct_matrix(const int N){
    CImg<float> *ptr_matrix = new CImg<float>(N,N,1,1,1/sqrt((float)N));
    const float c1 = sqrt(2.0f/N); 
    for (int x=0;x<N;x++){
        for (int y=1;y<N;y++){
            *ptr_matrix->data(x,y) = c1*(float)cos((cimg::PI/2/N)*y*(2*x+1));
        }
    }
    return ptr_matrix;
}
 
CImg<float>* GetMHKernel(float alpha, float level){
    int sigma = (int)(4*pow((float)alpha,(float)level));
    static CImg<float> *pkernel = NULL;
    float xpos, ypos, A;
    if (!pkernel){
        pkernel = new CImg<float>(2*sigma+1,2*sigma+1,1,1,0);
        cimg_forXY(*pkernel,X,Y){
            xpos = pow(alpha,-level)*(X-sigma);
            ypos = pow(alpha,-level)*(Y-sigma);
            A = xpos*xpos + ypos*ypos;
            pkernel->atXY(X,Y) = (2-A)*exp(-A/2);
        }
    }
    return pkernel;
}

int ph_dct_imagehash(const char* file,ulong64 &hash){

    if (!file){
        return -1;
    }
    CImg<uint8_t> src;
    try {
        src.load(file);
		if(src._is_error_by_yyrdl){
			return -1;
		}
    } catch (CImgIOException ex){
        return -1;
    }
    CImg<float> meanfilter(7,7,1,1,1);
    CImg<float> img;
    if (src.spectrum() == 3){
        img = src.RGBtoYCbCr().channel(0).get_convolve(meanfilter);
    } else if (src.spectrum() == 4){
        int width = img.width();
        int height = img.height();
        int depth = img.depth();
        img = src.crop(0,0,0,0,width-1,height-1,depth-1,2).RGBtoYCbCr().channel(0).get_convolve(meanfilter);
    } else {
        img = src.channel(0).get_convolve(meanfilter);
    }

    img.resize(32,32);
    CImg<float> *C  = ph_dct_matrix(32);
    CImg<float> Ctransp = C->get_transpose();

    CImg<float> dctImage = (*C)*img*Ctransp;

    CImg<float> subsec = dctImage.crop(1,1,8,8).unroll('x');;

    float median = subsec.median();
    ulong64 one = 0x0000000000000001;
    hash = 0x0000000000000000;
    for (int i=0;i< 64;i++){
        float current = subsec(i);
        if (current > median)
            hash |= one;
        one = one << 1;
    }

    delete C;
    return 0;
}
 
int ph_hamming_distance(const ulong64 hash1,const ulong64 hash2){
    ulong64 x = hash1^hash2;
    const ulong64 m1  = 0x5555555555555555ULL;
    const ulong64 m2  = 0x3333333333333333ULL;
    const ulong64 h01 = 0x0101010101010101ULL;
    const ulong64 m4  = 0x0f0f0f0f0f0f0f0fULL;
    x -= (x >> 1) & m1;
    x = (x & m2) + ((x >> 2) & m2);
    x = (x + (x >> 4)) & m4;
    return (x * h01)>>56;
}
  
int ph_bitcount8(uint8_t val){
    int num = 0;
    while (val){
        ++num;
        val &= val - 1;
    }
    return num;
}

 
double ph_hammingdistance2(uint8_t *hashA, int lenA, uint8_t *hashB, int lenB){
    if (lenA != lenB){
        return -1.0;
    }
    if ((hashA == NULL) || (hashB == NULL) || (lenA <= 0)){
        return -1.0;
    }
    double dist = 0;
    uint8_t D = 0;
    for (int i=0;i<lenA;i++){
        D = hashA[i]^hashB[i];
        dist = dist + (double)ph_bitcount8(D);
    }
    double bits = (double)lenA*8;
    return dist/bits;

}
 
#ifdef __cplusplus
}
#endif
#endif


 
