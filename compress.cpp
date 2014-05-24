#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv\cv.h>
#include <math.h>
#include <iostream>
#include <conio.h>
#include <cstdlib>

#define PI 3.14
#define N 8

using namespace std;
using namespace cv;




int zigzag[] = {
		0, 1, 8, 16, 9, 2, 3, 10,
		17, 24, 32, 25, 18, 11, 4, 5,
		12, 19, 26, 33, 40, 48, 41, 34,
		27, 20, 13, 6, 7, 14, 21, 28,
		35, 42, 49, 56, 57, 50, 43, 36,
		29, 22, 15, 23, 30, 37, 44, 51,
		58, 59, 52, 45, 38, 31, 39, 46,
		53, 60, 61, 54, 47, 55, 62, 63};


unsigned char yq[8][8] = {
		{16,11,	10,	16,	24,	40,	51,	61},
		{12,12,	14,	19,	26,	58,	60,	55},
		{14,13,	16,	24,	40,	57,	69,	56},
		{14,17,	22,	29,	51,	87,	80,	62},
		{18,22,	37,	56,	68,	109,103,77},
		{24,35,	55,	64,	81,	104,113,92},
		{49,64,	78,	87,	103,121,120,101},
		{72,92,	95,	98,	112,100,103,99}};


class rgbmac
{
public:
	rgbmac(){};
	unsigned char R;unsigned char G;unsigned char B;
};

class ycbcrmac
{
public:
	ycbcrmac(){};
	unsigned char Y;unsigned char Cb;unsigned char Cr;
};

class yuy
{
public:
	yuy();
	unsigned char *y;unsigned char *cb;unsigned char *cr;
};
yuy::yuy()
{
	y = (unsigned char*)malloc(256*sizeof(char));
	cb = (unsigned char*)malloc(64*sizeof(char));
	cr = (unsigned char*)malloc(64*sizeof(char));
}


class ysplit
{
public:
	ysplit();
	unsigned char *y1;unsigned char *y2;unsigned char *y3;unsigned char *y4;
};

ysplit::ysplit()
{
	y1 = (unsigned char*)malloc(64*sizeof(unsigned char));y2 = (unsigned char*)malloc(64*sizeof(unsigned char));
	y3 = (unsigned char*)malloc(64*sizeof(unsigned char));y4 = (unsigned char*)malloc(64*sizeof(unsigned char));
}

class ysplit2
{
public:
	ysplit2();
	int *y1;int *y2;int *y3;int *y4;
};

ysplit2::ysplit2()
{
	y1 = (int*)malloc(64*sizeof(int));y2 = (int*)malloc(64*sizeof(int));
	y3 = (int*)malloc(64*sizeof(int));y4 = (int*)malloc(64*sizeof(int));
}


class DCT
{
public:
	DCT();
	int *y;int *cb;int *cr;            //DCT can be negative
};
DCT::DCT()
{
	y = (int*)malloc(256*sizeof(int));
	cb = (int*)malloc(64*sizeof(int));
    cr = (int*)malloc(64*sizeof(int));
}

class remacro
{
public:
	remacro();
	char *y1;char *y2;char *y3;
	char *y4;char *cb;char *cr;
};

remacro::remacro()
{
	y1 = (char*)malloc(64*sizeof(char));y2 = (char*)malloc(64*sizeof(char));
	y3 = (char*)malloc(64*sizeof(char));y4 = (char*)malloc(64*sizeof(char));
	cb = (char*)malloc(64*sizeof(char));cr = (char*)malloc(64*sizeof(char));
}


class Rune
{
public:
	Rune::Rune(){
		int run=0;
	char level = 0;
	int last = 0;
	}
	int run;
	char level;
	int last;
};




int yuy_420(ycbcrmac*,yuy);								//YCrCb to YUY420                 
void yblock(yuy,ysplit);								//split YUY into four 8*8 parts
void dct_conv(DCT,ysplit,ysplit2);						//find DCT of Y component and for a ful 16*16 macroblock
void dct_min(unsigned char*,int**);						//find Discrete Cosine Transform(DCT)
void sim_quant(DCT,ysplit2);							//quantize the 8*8 macroblocks
void reorder(ysplit2,remacro,DCT);						//reorder the block to group zeros
Rune* runlevel(char*);									//run level encoding
void inv_sim_quant(DCT,ysplit2);						//reverse quantise the values
void inv_dct_min(int*,unsigned char**);					//inverse DCT
void idct_conv(unsigned char**,ysplit);					//converge 8*8 blocks of Y values
void inv_yuy(yuy,ycbcrmac*);							//covert from YUY420 to YCrCb

//void print_out(Rune*,FILE*);							//print runlevel encoded tuples to file

//// Huffman encoding has not been implemented///////////////////////////////



int main(int argc,char**argv)
{

	
	//FILE *out;
	//out = fopen("comprendenew.txt","w");
	
	//Rune *ans1,*ans2,*ans3,*ans4,*ans5,*ans6;           // object runlevel encoding
	//remacro remac;									  //


	ysplit yb;			 // stores 4 parts of Y 16*16 macroblock
	ysplit2 yb2;
	DCT dir;             // stores DCT of macroblock
	rgbmac rgb[256];     // stores RGB values of pixel
	ycbcrmac ycbcr[256]; // stores YCbCr values of pixel
	yuy yuymac;          // stores YUY 4:2:0 values of macroblock
	
	int row,col,wr;
	
	////////////////////////////////// Load the image//////////////////////////////////
	if( argc != 2)
	{
     cout <<" Usage: display_image ImageToLoadAndDisplay" << endl;
     return -1;
    }
    Mat image;
    image = imread(argv[1], CV_LOAD_IMAGE_COLOR);   // Read the file
    if(! image.data )                              // Check for invalid input
    {
        cout <<  "Could not open or find the image" << std::endl ;
        return -1;
    }

	row = (image.rows)/16;
	col = (image.cols)/16;
	Mat im = image(Range(0,col*16),Range(0,row*16));
	Mat image2;
	image2= Mat::zeros(col*16,row*16,CV_8UC3);


	///////////////////////////////////////////////////////////////////////////////////
	
	Mat imt;
	cv::cvtColor(im,imt,CV_BGR2YCrCb,0);
	CV_Assert(im.depth() != sizeof(uchar));


	for(int by=0;by<row;by++)          
	{
		for(int bx=0;bx<col;bx++)      
		{
			Mat_<Vec3b> _imt = imt;
			for(int i=0;i<16;i++)
			{
				for(int j=0;j<16;j++)
				{
					ycbcr[i*16 + j].Y  = _imt(i+(by*16),j+(bx*16))[0];
					ycbcr[i*16 + j].Cr = _imt(i+(by*16),j+(bx*16))[1];
					ycbcr[i*16 + j].Cb = _imt(i+(by*16),j+(bx*16))[2];

				}
			}
			
			wr = yuy_420(ycbcr,yuymac);			//convert to YUY420
			
			
			yblock(yuymac,yb);					//		split Y[256] into 4 parts
			dct_min(yuymac.cb,&dir.cb);			////
			dct_min(yuymac.cr,&dir.cr);			////	DCT of chrominance componenets stored in DCT dir object
			dct_conv(dir,yb,yb2);				//		DCT of 4* Y components stored in yb2 object /as a whole 16*16 in dir object
			sim_quant(dir,yb2);					//		quantize values of luminance and chrominance
			
			///////////////////////////////////////////////////////



			//carries out run level encoding,reodeing, prints them to file
			/*
			reorder(yb2,remac,dir);		//reorder for runlength encoding
			ans1 = runlevel((char*)yb2.y1);ans2 = runlevel((char*)yb2.y2);ans3 = runlevel((char*)yb2.y3);ans4 = runlevel((char*)yb2.y4);
			ans5 = runlevel((char*)dir.cb);ans6 = runlevel((char*)dir.cr);
			
			fprintf(out,"{");
			print_out(ans1,out);fprintf(out,".");
			print_out(ans2,out);fprintf(out,".");
			print_out(ans3,out);fprintf(out,".");
			print_out(ans4,out);fprintf(out,".");
			print_out(ans5,out);fprintf(out,".");
			print_out(ans6,out);fprintf(out,".");
			fprintf(out,"}");
			*/

			//     Decompression starts here

			inv_sim_quant(dir,yb2);
			inv_dct_min(yb2.y1,&yb.y1);			//		inverse DCT of 8*8 Y blocks
			inv_dct_min(yb2.y2,&yb.y2);			//
			inv_dct_min(yb2.y3,&yb.y3);			//
			inv_dct_min(yb2.y4,&yb.y4);			//

			inv_dct_min(dir.cb,&yuymac.cb);		//		inverse DCT of 8*8 chrominance blocks
			inv_dct_min(dir.cr,&yuymac.cr);		//

			idct_conv(&yuymac.y,yb);			//		converge 8*8 values to form yuy macroblock object
			inv_yuy(yuymac,ycbcr);				//		convert YUY420 values to YCrCb
		
			Mat_<Vec3b> _image2 = image2;
			for(int i=0;i<16;i++)
			{
				for(int j=0;j<16;j++)
				{
					_image2(i+(by*16),j+(bx*16))[0] = ycbcr[i*16+j].Y;
					_image2(i+(by*16),j+(bx*16))[1] = ycbcr[i*16+j].Cr;
					_image2(i+(by*16),j+(bx*16))[2] = ycbcr[i*16+j].Cb;

					
				}
			}
			cout<<"done till here"<<" block id "<<(by*(row) + bx)<<endl;
						
			
		}
	}
	//fclose(out);
	
	// show result of lossy compression
	Mat im3;
	cv::cvtColor(image2,im3,CV_YCrCb2BGR,0);
	namedWindow("Display",CV_WINDOW_NORMAL|CV_WINDOW_KEEPRATIO);
	imshow("Display",im3);
	imwrite("C:/Users/sushrut141/Pictures/not_img/out_img.tif",im3);
	waitKey(0);
	return 0;
	
	}




int yuy_420(ycbcrmac *ycbcr,yuy yuymac)
{
	int pix,r=0;
	for(int i=0;i<16;i++)
	{
		for(int j=0;j<16;j++)
		{
			pix = (16*i)+j;
			if (!(i&1)&&!(j&1))
			{
				yuymac.y[pix] =ycbcr[pix].Y;
				yuymac.cb[r] = ycbcr[pix].Cb;
				yuymac.cr[r] = ycbcr[pix].Cr;
				r++;
			}
			else
			{
				yuymac.y[pix] =ycbcr[pix].Y;
			}
		}
	}
	return r;
}

void yblock(yuy yuymac,ysplit yb)
{
	int bx,by,idx;
	for(bx=0;bx<2;bx++)
	{
		for(by=0;by<2;by++)
		{
			idx = bx*2 + by;
			switch(idx)
			{
			case 0:
				for(int i=0;i<8;i++)
				{
					for(int j=0;j<8;j++)
					{
						yb.y1[i*8+j] = yuymac.y[i*16+j];
					}
				}break;
				
			case 1:
				for(int i=0;i<8;i++)
				{
					for(int j=0;j<8;j++)
					{
						yb.y2[i*8+j] = yuymac.y[(8 + j) +(i*16)];
					}
				}break;
				
			case 2:
				for(int i=0;i<8;i++)
				{
					for(int j=0;j<8;j++)
					{
						yb.y3[i*8+j] = yuymac.y[(i+8)*16 + j];
					}
				}break;
				
			case 3:
				for(int i=0;i<8;i++)
				{
					for(int j=0;j<8;j++)
					{
						yb.y4[i*8+j] = yuymac.y[(i+8)*16+j+8];
					}
				}break;
			}
		}
	}
}

void dct_min(unsigned char *f,int**F)
{
	int u,v,i,j;
	double a[8];
	double sum,coef;
	a[0] = sqrt(1.0/8.0);
	for(int i=1;i<8;i++)
	{
		a[i] = sqrt(2.0/8.0);;
	}
	
	for ( u = 0; u < N; ++u ) 
	{
		for ( v = 0; v < N; ++v ) 
		{
			sum = 0.0;
			for ( i = 0; i < N; ++i ) 
			{
				for ( j = 0; j < N; ++j )
				{
					coef = cos((2*i+1)*u*PI/(2*N))*cos((2*j+1)*v*PI/(2*N));
					sum += ((double)(*(f+i*N+j))) * coef; //f[i][j] * coef
				} 
					(*F)[u*N + v] = (int)((a[u] * a[v] * sum)+0.5);
			} 
		} 
	} 

}


void dct_conv(DCT dir,ysplit yb,ysplit2 yb2)
{
	int bx,by,loc;
	dct_min(yb.y1,&yb2.y1);
	dct_min(yb.y2,&yb2.y2);
	dct_min(yb.y3,&yb2.y3);
	dct_min(yb.y4,&yb2.y4);
	for(bx=0;bx<2;bx++)
	{
		for(by=0;by<2;by++)
		{
			loc = bx*2 + by;
			switch(loc)
			{
			case 0:
				for(int i=0;i<8;i++)
				{
					for(int j=0;j<8;j++)
					{
						dir.y[i*8 + j] = yb2.y1[i*8 + j];
					}
				}break;
			case 1:
				for(int i=0;i<8;i++)
				{
					for(int j=0;j<8;j++)
					{
						dir.y[i*16 + j+8] = yb2.y2[i*8 + j];
					}
				}break;
			case 2:
				for(int i=0;i<8;i++)
				{
					for(int j=0;j<8;j++)
					{
						dir.y[(i+8)*16 + j] = yb2.y3[i*8 + j];
					}
				}break;
			case 3:
				for(int i=0;i<8;i++)
				{
					for(int j=0;j<8;j++)
					{
						dir.y[(i+8)*16 + j+8] = yb2.y4[i*8 + j];
					}
				}break;
			}
		}
	}

}      

void sim_quant(DCT dir,ysplit2 yb2)
{
	for(int i=0;i<8;i++)
	{
		for(int j=0;j<8;j++)
		{
			yb2.y1[i*8+j] = yb2.y1[i*8+j]/yq[i][j]; yb2.y2[i*8+j] = yb2.y2[i*8+j]/yq[i][j];
			yb2.y3[i*8+j] = yb2.y3[i*8+j]/yq[i][j]; yb2.y4[i*8+j] = yb2.y4[i*8+j]/yq[i][j];

			dir.cb[i*8+j] = dir.cb[i*8+j]/yq[i][j];dir.cr[i*8+j] = dir.cr[i*8+j]/yq[i][j];

		}
	}
}


void reorder ( ysplit2 yb2, remacro remac,DCT dir)
{
	int k, i1, j1;
	k = 0;
	for ( int i = 0; i < 8; i++ )
	{
		for ( int j = 0; j < 8; j++ )
		{
			i1 = zigzag[k] / 8;
			j1 = zigzag[k] % 8;
			remac.y1[i*8+j] = yb2.y1[i1*8+j1];remac.y2[i*8+j] = yb2.y2[i1*8+j1];
			remac.y3[i*8+j] = yb2.y3[i1*8+j1];remac.y4[i*8+j] = yb2.y4[i1*8+j1];
			remac.cb[i*8+j] = dir.cb[i1*8+j1];remac.cr[i*8+j] = dir.cr[i1*8+j1];
			k++;
		}
	}
}

Rune* runlevel(char *y1)
{
	Rune *run2;
	Rune *runs = (Rune*)malloc(64*sizeof(Rune));
	int k=0;
	unsigned char len=0;
	for(int i=0;i<64;i++)
	{
		if(y1[i]==0)
		{
			len++;
		}
		else
		{
		runs[k].run = len;
		runs[k].level = y1[i];
		runs[k].last = 0;
		k++;len = 0;}
	}
		if(k==0)
		{
			run2 = (Rune*)malloc(sizeof(Rune));
			run2[k].run = 63;
		run2[k].level = 0;
		run2[k].last = 1;
		free(runs);
		return run2;
		}
		else
		{
		//cout<<"it is"<<(int)runs[k-1].last<<" "<<(int)runs[k-2].last<<"///"<<endl;
		runs[k-1].last = 1;
			
	run2 = (Rune*)malloc((k)*sizeof(Rune));
	for(int i = 0;i<k;i++)
	{
		run2[i] = runs[i];
	}
	free(runs);
	return run2;}
}


void inv_sim_quant(DCT dir,ysplit2 yb2)
{
	for(int i=0;i<8;i++)
	{
		for(int j=0;j<8;j++)
		{
			yb2.y1[i*8+j] = yb2.y1[i*8+j]*(yq[i][j]); yb2.y2[i*8+j] = yb2.y2[i*8+j]*(yq[i][j]);
			yb2.y3[i*8+j] = yb2.y3[i*8+j]*(yq[i][j]); yb2.y4[i*8+j] = yb2.y4[i*8+j]*(yq[i][j]);

			dir.cb[i*8+j] = dir.cb[i*8+j]*(yq[i][j]);dir.cr[i*8+j] = dir.cr[i*8+j]*yq[i][j];

		}
	}
}

void inv_dct_min(int *F,unsigned char**f)
{

	int u,v,i,j;
	double a[8];
	double sum,coef;
	a[0] = sqrt(1.0/8.0);
	for(int i=1;i<8;i++)
	{
		a[i] = sqrt(2.0/8.0);;
	}
	
	for ( i = 0; i < 8; ++i ) 
	{
		for ( j = 0; j < 8; ++j ) 
		{
			sum = 0.0;
			for ( u = 0; u < 8; ++u ) 
			{
				for ( v = 0; v < 8; ++v )
				{
					coef = cos((2*i+1)*u*PI/(2*8))*cos((2*j+1)*v*PI/(2*8));
					sum += ((double)(*(F+(u*8)+v))) * coef*a[u] * a[v] ; //f[i][j] * coef
				} 
				(*f)[i*8 + j] = (unsigned char)(sum+0.5);
			} 
		} 
	} 
}


void idct_conv(unsigned char**y,ysplit yb)
{
	int loc;
	for(int by=0;by<2;by++)
	{
		for(int bx=0;bx<2;bx++)
		{
			loc = by*2 + bx;
			switch(loc)
			{
			case 0:
				for(int i=0;i<8;i++)
				{
					for(int j=0;j<8;j++)
					{
						(*y)[i*16 + j] = yb.y1[i*8 + j];
					}
				}break;
			case 1:
				for(int i=0;i<8;i++)
				{
					for(int j=0;j<8;j++)
					{
						(*y)[i*16 + (j+8)] = yb.y2[i*8 + j];
					}
				}break;
			case 2:
				for(int i=0;i<8;i++)
				{
					for(int j=0;j<8;j++)
					{
						(*y)[(i+8)*16 + j] = yb.y3[i*8 + j];
					}
				}break;
			case 3:
				for(int i=0;i<8;i++)
				{
					for(int j=0;j<8;j++)
					{
						(*y)[(i+8)*16 + (j+8)] = yb.y4[i*8 + j];
					}
				}break;
			}
		}
	}
		
}


void inv_yuy(yuy yb,ycbcrmac *ycbcr)
{
	int k = 0;
	for(int i=0;i<256;i++)
	{
		ycbcr[i].Y = yb.y[i];
	}
	for(int i=0;i<16;i++)
	{
		for(int j=0;j<16;j++)
		{
			
			if((!(i&1)&&!(j&1))==1)
			{
				
				ycbcr[i*16 + j].Cb = yb.cb[k];
				ycbcr[i*16 + (j+1)].Cb = yb.cb[k];
				ycbcr[(i+1)*16 + j].Cb = yb.cb[k];
				ycbcr[(i+1)*16 + j+1].Cb = yb.cb[k];

				ycbcr[i*16 + j].Cr = yb.cr[k];
				ycbcr[i*16 + (j+1)].Cr = yb.cr[k];
				ycbcr[(i+1)*16 + j].Cr = yb.cr[k];
				ycbcr[(i+1)*16 + j+1].Cr = yb.cr[k];
				k++;
			}
		}
	}
	

}

void print_out(Rune* a,FILE*out)
{
	for(int sc=0;;sc++)
			{
				fprintf(out,"%d",a[sc].run);
				fprintf(out,",");
				fprintf(out,"%d",a[sc].level);
				fprintf(out,",");
				fprintf(out,"%d",a[sc].last);
				fprintf(out,"|");
				if(a[sc].last==1)
					break;
			}
}
