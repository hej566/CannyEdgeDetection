#ifndef FONCTIONDEMO_H
#define FONCTIONDEMO_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define SQUARE(X) ((X)*(X))
#define MAX(i,j)  ((i)>(j)?(i):(j))
#define MIN(i,j)  ((i)>(j)?(i):(j))

#define NBCHAR 200

#define FFT   1
#define IFFT -1
#define FFT2D 2

#define GREY_LEVEL 255
#define PI 3.141592654

#define WHITE 255
#define BLACK 0


float*  fmatrix_allocate_1d(int);
float** fmatrix_allocate_2d(int,int);
void    free_fmatrix_1d(float*);
void    free_fmatrix_2d(float**);
float** LoadImagePgm(char*,int*,int*);
void    SaveImagePgm(char*,float**,int,int);

void    fourn(float*,unsigned long*,int,int);
void    FFTDD(float**,float**,int,int);
void    IFFTDD(float**,float**,int,int);
void    Mod(float**,float**,float**,int,int);
void    Mult(float**,float,int,int);
void    Recal(float**,int,int);
void    MultMatrix(float**,float**,float**,float**,float**,float**,int,int);
void    SquareMatrix(float**,float**,float**,float**,int,int);
float   funcgauss2D(float x,float y,float var);
void    compute_histo(float** mat,int lgth,int wdth,float* hist);

void    ComputeGuss(float** filter,float sigma,int length,int width);
void    Convolution(float** filter, float** ImgR, int length, int width);
void    Grad_X(float** ImgR, float** Temp_X, int length, int width);
void    Grad_Y(float** ImgR, float** Temp_Y, int length, int width);
void    Grad(float** Temp_X, float** Temp_Y, float** MagGrad, float** AngGrad, int length, int width);
void    Appro(float** AngGrad, int length,int width);
void    Suppe(float** MagGrad, float** AngGrad, float** SupGrad, int length,int width);
void    Hyste(int tau_L, int tau_H, float** SupGrad,float** HystGrad,float** AngGrad,int length, int width);


#endif
