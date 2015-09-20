#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "Function.h"


float* fmatrix_allocate_1d(int hsize)
{
  float* matrix;

  matrix=(float*)malloc(sizeof(float)*hsize); 
  if (matrix==NULL) printf("problem of allocation");

  return matrix; 
}


float** fmatrix_allocate_2d(int vsize,int hsize)
{
  int i;
  float** matrix;
  float *imptr;

  matrix=(float**)malloc(sizeof(float*)*vsize);
  if (matrix==NULL) printf("problem of allocation");

  imptr=(float*)malloc(sizeof(float)*hsize*vsize);
  if (imptr==NULL) printf("problem of allocation");
 
  for(i=0;i<vsize;i++,imptr+=hsize) matrix[i]=imptr;
  return matrix;
}


void free_fmatrix_1d(float* pmat)
 { 
  free(pmat); 
 }


void free_fmatrix_2d(float** pmat)
 { 
  free(pmat[0]);
  free(pmat);
 }


float** LoadImagePgm(char* name,int *length,int *width)
{
  int i,j,k;
  unsigned char var;
  char buff[NBCHAR];
  float** mat;

  char stringTmp1[NBCHAR],stringTmp2[NBCHAR],stringTmp3[NBCHAR];
 
  int ta1,ta2,ta3;
  FILE *fic;

  /*-----load image name.pgm-----*/
  strcpy(buff,name);
  strcat(buff,".pgm");
  printf("---> Open %s",buff);

  /*----open file----*/
  fic=fopen(buff,"r");
  if (fic==NULL)
    { printf("\n- error happen %s  -\n",buff);
      exit(-1); }

  /*--get image info--*/
  fgets(stringTmp1,100,fic);
  fgets(stringTmp2,100,fic);
  fscanf(fic,"%d %d",&ta1,&ta2);
  fscanf(fic,"%d\n",&ta3);

  /*--display image info--*/
  printf("\n\n--info--");
  printf("\n----------");
  printf("\n%s%s%d %d \n%d\n",stringTmp1,stringTmp2,ta1,ta2,ta3);

  *length=ta1;
  *width=ta2;
  mat=fmatrix_allocate_2d(*length,*width);
   
  /*--load image matrix--*/
     for(i=0;i<*length;i++)
      for(j=0;j<*width;j++)  
        { fread(&var,1,1,fic);
          mat[i][j]=var; }

  fclose(fic);

  return(mat);
}


void SaveImagePgm(char* name,float** mat,int length,int width)
 {
  int i,j,k;
  char buff[NBCHAR];
  FILE* fic;
  time_t tm;

  /*--extension--*/
  strcpy(buff,name);
  strcat(buff,".pgm");

  /*--open file--*/
  fic=fopen(buff,"w");
    if (fic==NULL) 
        { printf(" Problem of save image%s",buff); 
          exit(-1); }
  printf("\n Save %s in format pgm\n",name);

  /*--save info--*/
  fprintf(fic,"P5");
  if (ctime(&tm)==NULL) fprintf(fic,"\n#\n");
  else fprintf(fic,"\n# IMG Module, %s",ctime(&tm));
  fprintf(fic,"%d %d",width,length);
  fprintf(fic,"\n255\n");

  /*--save data--*/
   for(i=0;i<length;i++)
      for(j=0;j<width;j++) 
        fprintf(fic,"%c",(char)mat[i][j]);
   
  /*--close file--*/
   fclose(fic); 
 } 

/*--------------*/
/* FOURIER -----*/
/*--------------*/
/*------------------------------------------------*/
/*  FOURN ----------------------------------------*/
/*------------------------------------------------*/
void fourn(float data[], unsigned long nn[], int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}


/*----------------------------------------------------------*/
/* FFTDD                                                    */
/*----------------------------------------------------------*/
void FFTDD(float** mtxR,float** mtxI,int lgth, int wdth)
{
 int i,j;
 int posx,posy;

 float* data;
 float* ImgFreqR;
 float* ImgFreqI;
 unsigned long* nn;

 /*allocate memory*/
 data=(float*)malloc(sizeof(float)*(2*wdth*lgth)+1);
 ImgFreqR=(float*)malloc(sizeof(float)*(wdth*lgth));
 ImgFreqI=(float*)malloc(sizeof(float)*(wdth*lgth));
 nn=(unsigned long*)malloc(sizeof(unsigned long)*(FFT2D+1)); 

 /*fill nn*/
 nn[1]=lgth; nn[2]=wdth;

 /*fill data*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
   { data[2*(i*lgth+j)+1]=mtxR[i][j];
     data[2*(i*lgth+j)+2]=mtxI[i][j]; }

 /*FFTDD*/
 fourn(data,nn,FFT2D,FFT);

 /*fill data*/
 for(i=0;i<(wdth*lgth);i++)
  { ImgFreqR[i]=data[(2*i)+1];
    ImgFreqI[i]=data[(2*i)+2];  }

 /*matrix conversion*/
 for(i=0;i<(wdth*lgth);i++)
  { posy=(int)(i/wdth);
    posx=(int)(i%wdth);

    mtxR[posy][posx]=ImgFreqR[i];///(wdth*lgth);  
    mtxI[posy][posx]=ImgFreqI[i];}///(wdth*lgth); }

 /*free memory*/
 free(data);
 free(ImgFreqR);
 free(ImgFreqI);
 free(nn);
}


/*----------------------------------------------------------*/
/* IFFTDD                                                   */
/*----------------------------------------------------------*/
void IFFTDD(float** mtxR,float**  mtxI,int lgth,int wdth)
{
 int i,j;
 int posx,posy;

 float* data;
 float* ImgFreqR;
 float* ImgFreqI;
 unsigned long* nn;

 /*allocate memory*/
 data=(float*)malloc(sizeof(float)*(2*wdth*lgth)+1);
 ImgFreqR=(float*)malloc(sizeof(float)*(wdth*lgth));
 ImgFreqI=(float*)malloc(sizeof(float)*(wdth*lgth));
 nn=(unsigned long*)malloc(sizeof(unsigned long)*(FFT2D+1));

 /*fill nn*/
 nn[1]=lgth; nn[2]=wdth;

 /*fill data*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
   { data[2*(i*lgth+j)+1]=mtxR[i][j];
     data[2*(i*lgth+j)+2]=mtxI[i][j]; }

 /*FFTDD*/
 fourn(data,nn,FFT2D,IFFT);

 /*fill data*/
 for(i=0;i<(wdth*lgth);i++)
  { ImgFreqR[i]=data[(2*i)+1];
    ImgFreqI[i]=data[(2*i)+2]; }

 /*Matrix conversion*/
 for(i=0;i<(wdth*lgth);i++)
  { posy=(int)(i/wdth);
    posx=(int)(i%wdth);

   mtxR[posy][posx]=ImgFreqR[i]/(wdth*lgth);  
   mtxI[posy][posx]=ImgFreqI[i]/(wdth*lgth); }

 /*free memory*/
 free(data);
 free(ImgFreqR);
 free(ImgFreqI);
 free(nn);
}

void Mod(float** matM,float** matR,float** matI,int lgth,int wdth)
{
 int i,j;

 /*Calcul module*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
 matM[i][j]=sqrt((matR[i][j]*matR[i][j])+(matI[i][j]*matI[i][j]));
}


void Mult(float** mat,float coef,int lgth,int wdth)
{
 int i,j;

 /*multiplication*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
    { mat[i][j]*=coef;
      if (mat[i][j]>GREY_LEVEL) mat[i][j]=GREY_LEVEL; }
}


void Recal(float** mat,int lgth,int wdth)
{
 int i,j;
 float max,min;

 /*Initialization*/
 min=mat[0][0];

 /*look for min*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    if (mat[i][j]<min) min=mat[i][j];

 /*plus min*/
   for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    mat[i][j]-=min;

   max=mat[0][0];
 /*look for max*/
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
    if (mat[i][j]>max) max=mat[i][j];

 /*matrix caliberation*/
 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   mat[i][j]*=(GREY_LEVEL/max);      
}


void MultMatrix(float** matRout,float** matIout,float** mat1Rin,float** mat1Iin,
float** mat2Rin,float** mat2Iin,int lgth,int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   { matRout[i][j]=mat1Rin[i][j]*mat2Rin[i][j]-mat1Iin[i][j]*mat2Iin[i][j];
     matIout[i][j]=mat1Rin[i][j]*mat2Iin[i][j]+mat2Rin[i][j]*mat1Iin[i][j]; }
}


void SquareMatrix(float** matRout,float** matIout,float** matRin,float** matIin,int lgth,int wdth)
{
 int i,j;

 for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
   { matRout[i][j]=SQUARE(matRin[i][j])-SQUARE(matIin[i][j]);
     matIout[i][j]=2*matRin[i][j]*matIin[i][j]; }
}


float funcgauss2D(float x,float y,float var)
{
 return(exp(-(SQUARE(x)+SQUARE(y))/(2*var))/(2.0*PI*var));
}


void compute_histo(float** mat,int lgth,int wdth,float* hist)
{
 int i,j;

  for(i=0;i<=GREY_LEVEL;i++) hist[i]=0.0;

  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++)
    if ((mat[i][j]>=0)&&(mat[i][j]<=GREY_LEVEL))
       hist[(int)(mat[i][j])]++;
  
  for(i=0;i<=GREY_LEVEL;i++)  hist[i]/=(wdth*lgth);
}


void ComputeGuss(float** filter,float sigma,int length,int width)
{
	int i,j;
	float sum=0;
	for(i=0;i<length;i++){
		for(j=0;j<width;j++){
			filter[i][j]=funcgauss2D((float)(i-length/2),(float)(j-width/2),SQUARE(sigma));	
		}	
	}	
}

void Convolution(float** filter, float** ImgR, int length, int width)
{
	int i,j;
	float temp1,temp2;
	float** ImgI;
	float** filterI;
	float** OutI;
	float** OutR;

        filterI=fmatrix_allocate_2d(length,width);
	ImgI=fmatrix_allocate_2d(length,width);
	OutI=fmatrix_allocate_2d(length,width);	
	OutR=fmatrix_allocate_2d(length,width);

        FFTDD(ImgR,ImgI,length,width);
	FFTDD(filter,filterI,length,width);
	MultMatrix(OutR,OutI,filter,filterI,ImgR,ImgI,length,width);
	IFFTDD(OutR,OutI,length,width);
	
	for(i=0;i<length/2;i++){
		for(j=0;j<width/2;j++){
                        temp1=OutR[i][j];
			OutR[i][j]=OutR[i+length/2][j+width/2];
			OutR[i+length/2][j+width/2]=temp1;
			temp2=OutR[length/2+i][j];
			OutR[length/2+i][j]=OutR[i][width/2+j];
			OutR[i][width/2+j]=temp2;				
		}	
	}
    
	for(i=0;i<length;i++){
		for(j=0;j<width;j++){  
			ImgR[i][j]=OutR[i][j];		
		}	
	}
}

void Grad_X(float** ImgR, float** Temp_X, int length, int width)
{
	int i,j;
	for(i=0;i<length-1;i++){
		for(j=0;j<width;j++){
			
			 Temp_X[i][j]=ImgR[i][j]-ImgR[i+1][j];					
		}	
	}
}

void Grad_Y(float** ImgR, float** Temp_Y, int length, int width)
{
	int i,j;
	for(i=0;i<length;i++){
		for(j=0;j<width-1;j++){
			
			 Temp_Y[i][j]=ImgR[i][j]-ImgR[i][j+1];					
		}	
	}
}

void Grad(float** Temp_X, float** Temp_Y, float** MagGrad, float** AngGrad, int length, int width)
{
	int i,j;
	for(i=0;i<length;i++){
		for(j=0;j<width;j++){
			MagGrad[i][j]=sqrt(SQUARE(Temp_X[i][j])+SQUARE(Temp_Y[i][j]));	
			AngGrad[i][j]=atan2(-Temp_X[i][j], Temp_Y[i][j]);	
		}
	}

}

void Appro(float** AngGrad, int length,int width)
{
	int i,j;
	for(i=0;i<length;i++){
		for(j=0;j<width;j++){
			if(((AngGrad[i][j]+PI>=-PI/8)&&(AngGrad[i][j]+PI<=PI/8))||((AngGrad[i][j]>=-PI/8)&&(AngGrad[i][j]<=PI/8)))
			   AngGrad[i][j]=0;
			else if(((AngGrad[i][j]+PI>=PI/8)&&(AngGrad[i][j]+PI<=PI*3/8))||((AngGrad[i][j]>=PI/8)&&(AngGrad[i][j]<=PI*3/8)))
			   AngGrad[i][j]=45;
			else if(((AngGrad[i][j]+PI>=PI*3/8)&&(AngGrad[i][j]+PI<=PI*5/8))||((AngGrad[i][j]>=PI*3/8)&&(AngGrad[i][j]<=PI*5/8)))
			   AngGrad[i][j]=90;
			else if(((AngGrad[i][j]+PI>=PI*5/8)&&(AngGrad[i][j]+PI<=PI*7/8))||((AngGrad[i][j]>=PI*5/8)&&(AngGrad[i][j]<=PI*7/8)))
			   AngGrad[i][j]=135;
			else if (AngGrad[i][j]>=PI*7/8)
			   AngGrad[i][j]=0;
		}
	}
}

void Suppe(float** MagGrad, float** AngGrad, float** SupGrad, int length,int width)
{
	int i,j;
	for(i=1;i<length-1;i++){
		for(j=1;j<width-1;j++){
			if (AngGrad[i][j]==0){
				if((MagGrad[i][j]<MagGrad[i][j+1])||(MagGrad[i][j]<MagGrad[i][j-1])){
					SupGrad[i][j]=0;
				}
				else{
					SupGrad[i][j]=MagGrad[i][j];				
				}
			}
    			
			if (AngGrad[i][j]==45){
				if((MagGrad[i][j]<MagGrad[i-1][j+1])||(MagGrad[i][j]<MagGrad[i+1][j-1])){
					SupGrad[i][j]=0;
				}
				else{
					SupGrad[i][j]=MagGrad[i][j];				
				}
			}
			
			if (AngGrad[i][j]==90){
				if((MagGrad[i][j]<MagGrad[i-1][j])||(MagGrad[i][j]<MagGrad[i+1][j])){
					SupGrad[i][j]=0;
				}
				else{
					SupGrad[i][j]=MagGrad[i][j];				
				}
			}
			
			if (AngGrad[i][j]==135){
				if((MagGrad[i][j]<MagGrad[i-1][j-1])||(MagGrad[i][j]<MagGrad[i+1][j+1])){
					SupGrad[i][j]=0;
				}
				else{
					SupGrad[i][j]=MagGrad[i][j];				
				}
			}		
		}
	}
}

void Hyste(int tau_L, int tau_H, float** SupGrad,float** HystGrad,float** AngGrad,int length, int width)
{
	int i,j;

	for(i=1;i<length-1;i++){
		for(j=1;j<width-1;j++){
			if (SupGrad[i][j]>tau_H){
				HystGrad[i][j]=255;				
			}
                        if (SupGrad[i][j]<tau_L){
				HystGrad[i][j]=0.0;
			}
			else{
				if((AngGrad[i][j]==0)&&(i-1>0)){
					if ((SupGrad[i-1][j]==255)||(SupGrad[i+1][j]==255)){
						HystGrad[i][j]=255;					
					}				
				}
				if((AngGrad[i][j]==45)&&(i-1>0)&&(j-1>0)){
					if ((SupGrad[i-1][j-1]==255)||(SupGrad[i+1][j+1]==255)){
						HystGrad[i][j]=255;					
					}				
				}
				if((AngGrad[i][j]==90)&&(j-1>0)){
					if ((SupGrad[i][j-1]==255)||(SupGrad[i][j+1]==255)){
						HystGrad[i][j]=255;					
					}				
				}
				if((AngGrad[i][j]==135)&&(j-1>0)&&(i-1>0)){
					if ((SupGrad[i+1][j-1]==255)||(SupGrad[i-1][j+1]==255)){
						HystGrad[i][j]=255;					
					}				
				}										
			}	
		}
	}
}



