#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Function.h"

#define NAME_IMG_IN  "photograph"
#define NAME_IMG_CANNY "canny"

int main(int argc,int** argv)
{
  
  //Define variables
  int i,j,k,l;
  int length,width;
  float sum;
  float tau_L;
  float tau_H;
  float p_H;
  float sigma;
  float** filter;
  float** ImgR;
  float** Temp_X;
  float** Temp_Y;
  float** MagGrad;
  float** AngGrad;
  float** SupGrad;
  float** HystGrad;
  float*  hist;

  //Allocate memories
  ImgR=LoadImagePgm(NAME_IMG_IN,&length,&width);
  hist=fmatrix_allocate_1d(GREY_LEVEL);
  filter=fmatrix_allocate_2d(length,width);
  Temp_X=fmatrix_allocate_2d(length,width);
  Temp_Y=fmatrix_allocate_2d(length,width);
  MagGrad=fmatrix_allocate_2d(length,width);
  AngGrad=fmatrix_allocate_2d(length,width);
  SupGrad=fmatrix_allocate_2d(length,width);
  HystGrad=fmatrix_allocate_2d(length,width);
  
  //Initialize arrays
  for(i=0;i<length;i++){
     for(j=0;j<width;j++){
        filter[i][j]=0.0;
        Temp_X[i][j]=0.0;
        Temp_Y[i][j]=0.0;
	MagGrad[i][j]=0.0;
	AngGrad[i][j]=0.0;
 	SupGrad[i][j]=0.0;
    	HystGrad[i][j]=0.0;
     }
  }
	
  //Get value from keyboard
  printf("Input the value of tau_L: ");
  scanf("%f",&tau_L);
  printf("Input the value of tau_H: ");
  scanf("%f",&tau_H);
  printf("Input the sigma of filter Gaussien: ");
  scanf("%f",&sigma);
  
  //Implement Canny algo
  ComputeGuss(filter,sigma,length,width);
  Convolution(filter,ImgR,length,width);
  Grad_X(ImgR,Temp_X,length,width);
  Grad_Y(ImgR,Temp_Y,length,width);
  Grad(Temp_X, Temp_Y, MagGrad, AngGrad, length, width);
  Recal(MagGrad,length,width);
  Appro(AngGrad,length,width);
  Suppe(MagGrad,AngGrad,SupGrad,length,width);
  Recal(SupGrad,length,width);
  Hyste(tau_L,tau_H,SupGrad,HystGrad,AngGrad,length,width);
 
  SaveImagePgm(NAME_IMG_CANNY,HystGrad,length,width);
  system("display canny.pgm&");

  printf("\n Ending ... \n\n\n");
  return 0; 	 
}
