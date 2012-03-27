#include "cv.h"
#include "highgui.h"
#include <stdio.h>
int inRange ( float tmin, float tmax, float data);
void SkinFilter_LUV(IplImage* inputImg, IplImage * pmask);
void SkinFilter_NEW (IplImage *inputImg, IplImage * pmask);
CvScalar BGR2New ( CvScalar pix);
void CvtBgr2New (IplImage *src, IplImage *dst);


CvScalar BGR2Hsv ( CvScalar pix);
void CvtBgr2Hsv (IplImage *src, IplImage *dst);
void SkinFilter_HSV (IplImage *inputImg, IplImage * pmask);

CvScalar BGR2Hist ( CvScalar pix);
void CvtBgr2Hist (IplImage *src, IplImage *dst);
void SkinFilter_RGB (IplImage *inputImg, IplImage * pmask);

CvScalar BGR2Cbcr ( CvScalar pix);
void CvtBgr2Cbcr (IplImage *src, IplImage *dst);
void SkinFilter_Cbcr (IplImage *inputImg, IplImage * pmask);

void Illum_Comp_YCrCb ( IplImage *src, IplImage *dst);

void grey_World ( IplImage *src);

void imgShow ( char *name, IplImage *img);
void imgRelease ( char *name, IplImage *img);

void morPho ( IplImage *src, IplImage *dst);
