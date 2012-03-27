#include <stdio.h>
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"
#include "cvaux.h"
#include "skinFilter.h"


int main( int argc, char ** argv)
{
	
	
	IplImage * src = cvLoadImage(argv[1], 1);
	// mask imgage
	IplImage * dst_luv = cvCreateImage ( cvGetSize(src), src->depth, 1);
	IplImage * dst_new = cvCreateImage ( cvGetSize(src), src->depth, 1);
	IplImage * dst_rgb = cvCreateImage ( cvGetSize(src), src->depth, 1);
	IplImage * dst_hsv = cvCreateImage ( cvGetSize(src), src->depth, 1);
	IplImage * dst_CbCr = cvCreateImage ( cvGetSize(src), src->depth, 1);

	// pre-processing of Illumination Compensation 
	//    using grey-world

	//grey_World ( src);
		
	// initial
	cvZero (dst_luv);
	cvZero (dst_new);
	cvZero (dst_rgb);
	cvZero (dst_hsv);
	cvZero (dst_CbCr);
	// skin Filter eg smooth file 
	SkinFilter_LUV ( src, dst_luv);
	SkinFilter_NEW ( src, dst_new);
	SkinFilter_HSV ( src, dst_hsv);
	SkinFilter_RGB ( src, dst_rgb);
	SkinFilter_Cbcr( src, dst_CbCr);	
	
	// morphology
	morPho ( dst_luv, dst_luv);
	morPho ( dst_new, dst_new);
	morPho ( dst_CbCr, dst_CbCr);

	cvSaveImage("luv.jpg", dst_luv, NULL);
    cvSaveImage("Cbcr.jpg", dst_CbCr, NULL);
    cvSaveImage("new.jpg", dst_new, NULL);

	

/*
	imgShow ( "src", src);
	imgShow ( "luv", dst_luv);
	imgShow ( "new", dst_new);
	imgShow ( "rgb", dst_rgb);
	imgShow ( "hsv", dst_hsv);
	imgShow ( "Cbcr", dst_CbCr);
*/
//	cvWaitKey (0);
	imgRelease ( "src", src);
    imgRelease ( "luv", dst_luv);
    imgRelease ( "new", dst_new);
    imgRelease ( "rgb", dst_rgb);
    imgRelease ( "hsv", dst_hsv);
	imgRelease ( "Cbcr", dst_CbCr);
	
	return 0;
}
