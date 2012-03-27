/*
	skin color segmentation module
			
	convert RGB iamge to YUV YCbcr new respectively using the 
		human face color threshold

	copyright France Telecom HCI
	
	author: jiwei  buptjiwei@gmail.com

	2011.3.29  the first edition 

*/
#include "skinFilter.h"
#define MAX(a,b) ((a) < (b) ? (b) : (a))

// range in threshold
int inRange ( float tmin, float tmax, float data)
{
	int flag;
	if ( data >= tmin && data <= tmax)
                return 1;
        else
                return 0;

}

// comp in three 
double mymin ( double a, double b, double c)
{
	double min = a;
	if ( min > b) min = b;
	if ( min > c) min = c;
	return min;
}


double mymax ( double a, double b, double c)
{
	double max = a;
	if ( max < b) max = b;
	if ( max < c) max = c;
	return max;
}



// RGB ==> LUV
void SkinFilter_LUV(IplImage* inputImg, IplImage * pmask)
{
	
	IplImage * imgSmth = cvCreateImage ( cvGetSize( inputImg), inputImg->depth, inputImg->nChannels);
	IplImage * imgLuv = cvCreateImage ( cvGetSize( inputImg), inputImg->depth, inputImg->nChannels);
	cvSmooth ( inputImg, imgSmth, CV_BLUR, 3, 3, 0, 0);
	cvCvtColor ( imgSmth, imgLuv, CV_RGB2Luv);
/*
	int winWid = src->width/4;
        int winHei = src->height/4;
        int xstep = src->width/8;
        int ystep =  src->height/8;
*/
	int winWid = 1;
        int winHei = 1;
        int xstep = 1;
        int ystep =  1;
        int x=0, y=0;
	CvScalar meanimgLuv;
	CvScalar temp;
        temp.val[0] = 255;
        float thesUmin = 71.0, thesUmax = 104.0;
        float thesVmin = 67.0, thesVmax = 140.0;

	for ( y=0; y<(imgLuv->height - winHei); y+=ystep)
	{
		for ( x=0; x<(imgLuv->width - winWid); x+=xstep)
		{
			cvSetImageROI ( imgLuv, cvRect(x, y, winWid, winHei));
			cvAvgSdv ( imgLuv, &meanimgLuv, NULL, NULL);
			//printf (" %f %f\n", meanimgLuv.val[1], meanimgLuv.val[2]);
			if ( inRange ( thesUmin, thesUmax, meanimgLuv.val[1]) && inRange(thesVmin, thesVmax, meanimgLuv.val[2]))
			{
                                cvSetImageROI(pmask, cvRect( x, y, winWid, winHei));
                                cvSet ( pmask, temp, NULL);
                                cvResetImageROI ( pmask);
                        }

                        cvResetImageROI (imgLuv);

		}
	}
	cvReleaseImage ( &imgSmth);
        cvReleaseImage ( &imgLuv);


}

CvScalar BGR2New ( CvScalar pix)
{
	double res = 0.0;
	double res1 = 0.0;
	CvScalar ret;
	ret.val[0] = 0.0;
	res = 0.140209042551032500 * pix.val[0] + 0.587043074451121360* pix.val[1]+ 0.298936021293775390*pix.val[2];
	res1 = MAX ( pix.val[0], pix.val[1]);
	//ret.val[0] = ((res - res1<1.0)? (res-res1):0.9999);
	ret.val[0] = res - res1;
	return ret;
}

// like cvCvtColor RGB --> new

void CvtBgr2New (IplImage *src, IplImage *dst)
{
	IplImage * img32F = cvCreateImage ( cvGetSize( src), IPL_DEPTH_32F, src->nChannels);
	double scale = 1.0 / 255.0;
	cvConvertScale ( src, img32F, scale, 0);
	CvScalar pix;
	CvScalar value;
	
	int x,y;

	for ( y = 0; y < src->height; y++)
	{
		for ( x = 0; x < src->width; x++)
		{
			pix = cvGet2D ( img32F, y, x);
			value = BGR2New (pix);
			cvSet2D ( dst, y, x, value);
		}
	}

	cvReleaseImage ( &img32F);


}


// RGB --> new
void SkinFilter_NEW (IplImage *inputImg, IplImage * pmask)
{
	IplImage * imgSmth = cvCreateImage ( cvGetSize( inputImg), inputImg->depth, inputImg->nChannels);
        IplImage * imgNew = cvCreateImage ( cvGetSize( inputImg), IPL_DEPTH_32F, 1);
        cvSmooth ( inputImg, imgSmth, CV_BLUR, 3, 3, 0, 0);
        //cvCvtColor ( imgSmth, imgLuv, CV_RGB2Luv);
	CvtBgr2New ( imgSmth, imgNew);
/*
        int winWid = src->width/4;
        int winHei = src->height/4;
        int xstep = src->width/8;
        int ystep =  src->height/8;
*/

        int winWid = 1;
        int winHei = 1;
        int xstep = 1;
        int ystep =  1;
        int x=0, y=0;
        CvScalar meanimgNew;
        CvScalar temp;
        temp.val[0] = 255;
	double thresmin = 0.0251;
	double thresmax = 0.1177;
        for ( y=0; y<(imgNew->height - winHei); y+=ystep)
        {
                for ( x=0; x<(imgNew->width - winWid); x+=xstep)
                {
                        cvSetImageROI ( imgNew, cvRect(x, y, winWid, winHei));
                        cvAvgSdv ( imgNew, &meanimgNew, NULL, NULL);
                        //printf (" %f %f\n", meanimgLuv.val[1], meanimgLuv.val[2]);
                        if ( inRange (thresmin, thresmax, meanimgNew.val[0]))
                        {
                                cvSetImageROI(pmask, cvRect( x, y, winWid, winHei));
                                cvSet ( pmask, temp, NULL);
                                cvResetImageROI ( pmask);
                        }

                        cvResetImageROI (imgNew);

                }
        }

	cvReleaseImage ( &imgSmth);
	cvReleaseImage ( &imgNew);
}


// cvCvt BGR-->HSV
CvScalar BGR2Hsv ( CvScalar pix)
{
	CvScalar ret;
	ret.val[0] = 0.0;
	int B,G,R;
	double H1, H, S, V;
	B = pix.val[0];
	G = pix.val[1];
	R = pix.val[2];
	H1 = acos( 0.5*((R-G)+(R-B)) / sqrt( (double)(R-G)*(R-G) + (R-B)*(G-B) ) ) *180.0/CV_PI;
	if (B<=G) { H = H1;}
	else { H = 360 - H1;}
	S = (double)( mymax(R,G,B) - mymin(R,G,B) )/ mymax(R,G,B);
	V = (double)( mymax(R,G,B) ) / 255;
	ret = cvScalar ( H, S, V, 0);
	return ret;
	
}

// cvt bgr--> hsv
void CvtBgr2Hsv (IplImage *src, IplImage *dst)
{
        IplImage * img32F = cvCreateImage ( cvGetSize( src), IPL_DEPTH_32F, src->nChannels);
        double scale = 1.0 / 255.0;
        cvConvertScale ( src, img32F, scale, 0);
        CvScalar pix;
        CvScalar value;

        int x,y;

        for ( x=0; x < src->width; x++)
        {
                for (y=0; y < src->height; y++)
                {
                        pix = cvGet2D ( img32F, y, x);
                        value = BGR2Hsv ( pix);
                        cvSet2D ( dst, y, x, value);

                }
        }

	cvReleaseImage ( &img32F);
}


// RGB --> HSV (not work)

void SkinFilter_HSV (IplImage *inputImg, IplImage * pmask)
{
	IplImage *imgSmth = cvCreateImage ( cvGetSize (inputImg), inputImg->depth, inputImg->nChannels);
	//IplImage *imgHsv = cvCreateImage ( cvGetSize (inputImg), IPL_DEPTH_32F, inputImg->nChannels);
	IplImage *imgHsv = cvCreateImage ( cvGetSize (inputImg), inputImg->depth, inputImg->nChannels);
	cvSmooth ( inputImg, imgSmth, CV_BLUR, 3, 3, 0, 0);
	//cvCvtColor ( imgSmth, imgHsv, CV_BGR2HSV);
	//CvtBgr2Hsv ( imgSmth, imgHsv);
	cvCvtColor ( imgSmth, imgHsv, CV_BGR2HSV);
	int winWid = 1;
        int winHei = 1;
        int xstep = 1;
        int ystep =  1;
        int x=0, y=0;
	CvScalar meanimgHsv;
	CvScalar temp;
	temp.val[0] = 255;
	double thesHmin = 0.00, thesHmax = 50.0;
	double thesSmin = 0.20, thesSmax = 0.68;
	double thesVmin = 0.35, thesVmax = 1.00;
	
	for ( y=0; y<(imgHsv->height - winHei); y+=ystep)
        {
                for ( x=0; x<(imgHsv->width - winWid); x+=xstep)
                {
                        cvSetImageROI ( imgHsv, cvRect(x, y, winWid, winHei));
                        cvAvgSdv ( imgHsv, &meanimgHsv, NULL, NULL);
                        //printf (" %f %f\n", meanimgLuv.val[1], meanimgLuv.val[2]);
                        if ( inRange (thesHmin, thesHmax, meanimgHsv.val[0]) && inRange ( thesSmin, thesSmax, meanimgHsv.val[1]) && inRange ( thesVmin, thesVmax, meanimgHsv.val[2]))
                        {
                                cvSetImageROI(pmask, cvRect( x, y, winWid, winHei));
                                cvSet ( pmask, temp, NULL);
                                cvResetImageROI ( pmask);
                        }

                        cvResetImageROI (imgHsv);

                }
        }

	cvReleaseImage ( &imgSmth);
	cvReleaseImage ( &imgHsv);
	
	

}


// BGR histogram
CvScalar BGR2Hist ( CvScalar pix)
{
	CvScalar ret;
	int B, G, R, sum;
	double nB, nG, nR;
	B = pix.val[0];
	G = pix.val[1];
	R = pix.val[2];

	sum = B + G + R;
	nB = (double) B / sum;
	nG = (double) G / sum;
	nR = (double) R / sum;

	ret.val[0] = nB;
	ret.val[1] = nG;
	ret.val[2] = nR;

	return ret;
}

// cvt
void CvtBgr2Hist ( IplImage *src, IplImage *dst)
{
	IplImage * img32F = cvCreateImage ( cvGetSize (src), IPL_DEPTH_32F, src->nChannels);
	double scale = 1.0 / 255;
	cvConvertScale ( src, img32F, scale, 0);
	CvScalar pix, value;
	int x, y;
	for ( y=0; y<src->height; y++)
	{
		for ( x=0; x<src->width; x++)
		{
			pix = cvGet2D ( img32F, y, x);
			value = BGR2Hist ( pix);
			cvSet2D ( dst, y, x, value);
		} 
	}

	cvReleaseImage ( &img32F);
}

void SkinFilter_RGB (IplImage *inputImg, IplImage * pmask)
{
	IplImage * imgSmth = cvCreateImage ( cvGetSize( inputImg), inputImg->depth, inputImg->nChannels);
	IplImage * imgHist = cvCreateImage ( cvGetSize( inputImg), IPL_DEPTH_32F, inputImg->nChannels);
	cvSmooth ( inputImg, imgSmth, CV_BLUR, 3, 3, 0, 0);
	CvtBgr2Hist ( imgSmth, imgHist);
	int winWid = 1;
        int winHei = 1;
        int xstep = 1;
        int ystep =  1;
        int x=0, y=0;
	CvScalar meanimgHist;
	CvScalar temp;
        temp.val[0] = 255;
        float thesGmin = 0.28, thesGmax = 0.363;
        float thesRmin = 0.36, thesRmax = 0.465;
	for ( y=0; y<(imgHist->height - winHei); y+=ystep)
        {
                for ( x=0; x<(imgHist->width - winWid); x+=xstep)
                {
                        cvSetImageROI ( imgHist, cvRect(x, y, winWid, winHei));
                        cvAvgSdv ( imgHist, &meanimgHist, NULL, NULL);
                        //printf (" %f %f\n", meanimgLuv.val[1], meanimgLuv.val[2]);
                        if ( inRange (thesGmin, thesGmax, meanimgHist.val[1]) && inRange ( thesRmin, thesRmax, meanimgHist.val[2]))
                        {
                                cvSetImageROI(pmask, cvRect( x, y, winWid, winHei));
                                cvSet ( pmask, temp, NULL);
                                cvResetImageROI ( pmask);
                        }

                        cvResetImageROI (imgHist);

                }
        }

	cvReleaseImage ( &imgSmth);
	cvReleaseImage ( &imgHist);
		
}


// BGR --> Cbcr
void SkinFilter_Cbcr(IplImage* inputImg, IplImage * pmask)
{

        IplImage * imgSmth = cvCreateImage ( cvGetSize( inputImg), inputImg->depth, inputImg->nChannels);
        IplImage * imgYCbcr = cvCreateImage ( cvGetSize( inputImg), inputImg->depth, inputImg->nChannels);
	IplImage * imgIllumYCb =  cvCreateImage ( cvGetSize( inputImg), inputImg->depth, inputImg->nChannels);
        cvSmooth ( inputImg, imgSmth, CV_BLUR, 3, 3, 0, 0);
        cvCvtColor ( imgSmth, imgYCbcr, CV_BGR2YCrCb);
	
	CvScalar avgYcbcr;
	avgYcbcr = cvAvg ( imgYCbcr, NULL);
	if ( avgYcbcr.val[0] <60.0 || avgYcbcr.val[0] > 125.0)
	{
		//printf("YCbCr!\n");
	//	Illum_Comp_YCrCb ( imgYCbcr, imgYCbcr);
		//printf("Illuminated compensation!\n");
	}
/*
        int winWid = src->width/4;
        int winHei = src->height/4;
        int xstep = src->width/8;
        int ystep =  src->height/8;
*/
        int winWid = 1;
        int winHei = 1;
        int xstep = 1;
        int ystep =  1;
        int x=0, y=0;
        CvScalar meanimgYCbcr;
        CvScalar temp;
        temp.val[0] = 255;
        double thesCbmin = 77.0, thesCbmax = 127.0;
        double thesCrmin = 133.0, thesCrmax = 173.0;

        for ( y=0; y<(imgYCbcr->height - winHei); y+=ystep)
        {
                for ( x=0; x<(imgYCbcr->width - winWid); x+=xstep)
                {
                        cvSetImageROI ( imgYCbcr, cvRect(x, y, winWid, winHei));
                        cvAvgSdv ( imgYCbcr, &meanimgYCbcr, NULL, NULL);
                        //printf (" %f %f\n", meanimgLuv.val[1], meanimgLuv.val[2]);
                        if ( inRange ( thesCbmin, thesCbmax, meanimgYCbcr.val[2]) && inRange(thesCrmin, thesCrmax, meanimgYCbcr.val[1]))
                        {
                                cvSetImageROI(pmask, cvRect( x, y, winWid, winHei));
                                cvSet ( pmask, temp, NULL);
                                cvResetImageROI ( pmask);
                        }

                        cvResetImageROI (imgYCbcr);

                }
        }
        cvReleaseImage ( &imgSmth);
        cvReleaseImage ( &imgYCbcr);


}


// Illumination Compensation just for YCrcb
void Illum_Comp_YCrCb ( IplImage *src, IplImage *dst)
{
        IplImage * img_Y = cvCreateImage ( cvGetSize (src), src->depth, 1);
        IplImage * img_Cb = cvCreateImage ( cvGetSize (src), src->depth, 1);
        IplImage * img_Cr = cvCreateImage ( cvGetSize (src), src->depth, 1);

        cvSplit ( src, img_Y, img_Cb, img_Cr, NULL);

//        cvNamedWindow ("y", 1);
//        cvShowImage ( "y", img_Y);
        int x,y, ix,iy, locx,locy, fx, fy, smx, smy;
        int winW = 10, winH = 10;
        CvScalar meanY;
        ix = src->width / winW;
        iy = src->height / winH;
        fx = src->width % winW;
        fy = src->height % winH;
        locx=0;
        locy=0;
        double minY;
        if ( fx) {smx = ix + 1;}
        else {smx = ix;}
        if ( fy) {smy = iy + 1;}
        else {smy = iy;}
//      printf( "%d %d| %d %d\n",smx, smy, ix, iy);
        IplImage * maskY = cvCreateImage ( cvSize(smx, smy),IPL_DEPTH_8U,1);
        cvZero (maskY);


        for ( x=0; x<ix; x++)
        {
                locy=0;
                for ( y=0; y<iy; y++)
                {
//                      printf( "%d %d %d %d \n",x,y,locx,locy );
                        cvSetImageROI ( img_Y, cvRect(locx, locy, winW, winH));
                        //meanY = cvAvg ( img_Y, NULL);
                        cvMinMaxLoc ( img_Y, &minY, NULL, NULL, NULL, NULL);
                        meanY.val[0] = minY;
                        cvSet2D ( maskY, y, x, meanY);
                        cvResetImageROI ( img_Y);
                        locy = locy + winH;
                }
                locx = locx + winW;
        }

        IplImage *resizeY =  cvCreateImage ( cvGetSize ( src), IPL_DEPTH_8U, 1);
        cvResize ( maskY, resizeY, CV_INTER_LINEAR);

        IplImage * subRes =  cvCreateImage ( cvGetSize ( src), IPL_DEPTH_8U, 1);
        cvSub ( img_Y, resizeY, subRes, NULL);
/*
        int ri, rj;
        CvScalar pixsub;
        for ( rj=0; rj<subRes->height; rj++)
        {
                for ( ri=0; ri<subRes->width; ri++)
                {
                        pixsub = cvGet2D ( subRes, rj, ri);
                        if ( pixsub.val[0] < 0.000000)
                        {
                                pixsub.val[0] = pixsub.val[0] * (-1.0);
                                cvSet2D ( subRes, rj, ri, pixsub);
                        }
                }
        }
*/
        // EqualizeHist
        IplImage * histY = cvCreateImage ( cvGetSize ( src), IPL_DEPTH_8U, 1);


        cvEqualizeHist ( subRes, histY);

        cvMerge ( histY, img_Cb, img_Cr, NULL, dst);


//       cvNamedWindow ("ysub", 1);
//     cvShowImage ("ysub", subRes);
//      cvNamedWindow ("ymask", 1);
//     cvShowImage ("ymask", resizeY);
    cvNamedWindow ("yhist", 1);
    cvShowImage ("yhist", histY);

        cvReleaseImage ( &histY);
        cvReleaseImage ( &subRes);
        cvReleaseImage ( &resizeY);
        cvReleaseImage ( &maskY);
        cvReleaseImage ( &img_Y);
        cvReleaseImage ( &img_Cb);
        cvReleaseImage ( &img_Cr);
//        cvReleaseImage ( &histY);

        //cvWaitKey(0);
}

// Illumination compensation using grey world algorithm
void grey_World ( IplImage *src)
{
	CvScalar meanRgb;
        meanRgb = cvAvg ( src, NULL);

        double avgGray;
        avgGray = ( meanRgb.val[0] +  meanRgb.val[1] +  meanRgb.val[2]) / 3.0;

        IplImage * dst_r = cvCreateImage ( cvGetSize (src), IPL_DEPTH_32F, 1);
        IplImage * dst_g = cvCreateImage ( cvGetSize (src), IPL_DEPTH_32F, 1);
        IplImage * dst_b = cvCreateImage ( cvGetSize (src), IPL_DEPTH_32F, 1);
        IplImage * dst32F = cvCreateImage ( cvGetSize (src), IPL_DEPTH_32F, 3);
        IplImage * merge = cvCreateImage ( cvGetSize (src), IPL_DEPTH_32F, 3);



        double scaler, scaleg, scaleb;
        cvConvertScale ( src, dst32F, 1.0, 0.0);

        scaler = meanRgb.val[0] / avgGray;
        scaleg = meanRgb.val[1] / avgGray;
        scaleb = meanRgb.val[2] / avgGray;

        cvSplit ( dst32F, dst_r, dst_g, dst_b, NULL);
	
	if ( scaler > 0.0000000000001)
        {
                cvConvertScale ( dst_r, dst_r, scaler, 0);
        }

        if ( scaleg > 0.0000000000001)
        {
                cvConvertScale ( dst_g, dst_g, scaleg, 0);
        }

        if ( scaleb > 0.0000000000001)
        {
                cvConvertScale ( dst_b, dst_b, scaleb, 0);
        }

        double max_r, max_g, max_b, factor;
        double max_r1, max_g1, max_b1;;
//        max_r = imgMax ( dst_r);
//        max_g = imgMax ( dst_g);
//        max_b = imgMax ( dst_b);

	cvMinMaxLoc ( dst_r, NULL, &max_r, NULL, NULL, NULL);
        cvMinMaxLoc ( dst_g, NULL, &max_g, NULL, NULL, NULL);
        cvMinMaxLoc ( dst_b, NULL, &max_b, NULL, NULL, NULL);


	factor = max_r;
        if ( factor < max_g) factor = max_g;
        if ( factor < max_b) factor = max_b;

        double avgScale;
        avgScale = 1.0 / factor;

        cvMerge ( dst_r, dst_g, dst_b, NULL, merge);

	if ( factor > 1)
        {
                cvConvertScale ( merge, merge, avgScale, 0.0);
        }

	IplImage * mer_r = cvCreateImage ( cvGetSize (src), IPL_DEPTH_32F, 1);
        IplImage * mer_g = cvCreateImage ( cvGetSize (src), IPL_DEPTH_32F, 1);
        IplImage * mer_b = cvCreateImage ( cvGetSize (src), IPL_DEPTH_32F, 1);


	cvSplit ( merge, mer_r, mer_g, mer_b, NULL);
        cvNormalize ( mer_r, mer_r, 0.0, 255.0, CV_MINMAX, NULL);
        cvNormalize ( mer_g, mer_g, 0.0, 255.0, CV_MINMAX, NULL);
        cvNormalize ( mer_b, mer_b, 0.0, 255.0, CV_MINMAX, NULL);
        cvMerge ( mer_r, mer_g, mer_b, NULL, dst32F);
        cvConvertScale ( dst32F, src, 1.0, 0.0);

	cvReleaseImage ( &mer_r);
	cvReleaseImage ( &mer_g);
	cvReleaseImage ( &mer_b);
        cvReleaseImage ( &dst32F);
        cvReleaseImage ( &merge);
        cvReleaseImage ( &dst_r);
        cvReleaseImage ( &dst_g);
        cvReleaseImage ( &dst_b);
        //cvReleaseImage ( &);





	
}





void imgShow ( char *name, IplImage *img)
{
	cvNamedWindow ( name, 1);
	cvShowImage ( name, img);
//	cvWaitKey (0);
//	cvDestroyWindow ( name);
//	cvReleaseImage ( &img);
}

void imgRelease ( char *name, IplImage *img)
{
	cvDestroyWindow ( name);
        cvReleaseImage ( &img);
	
}

// morphology the last step, Close --> Open

void morPho ( IplImage *src, IplImage *dst)
{
	int nSize = 3;
	IplConvKernel *structElem1 = cvCreateStructuringElementEx ( nSize +2,
	nSize,
	floor ( ( nSize + 2) / 2.0),
	floor ( nSize / 2.0),
	CV_SHAPE_ELLIPSE,
	NULL);

	cvMorphologyEx ( src, dst, NULL, structElem1, CV_MOP_CLOSE, 1);
	

	nSize = 11;
	IplConvKernel *structElem2 = cvCreateStructuringElementEx ( nSize,
	nSize,
	floor ( nSize / 2.0),
	floor ( nSize / 2.0),
	CV_SHAPE_ELLIPSE,
	NULL);
	cvMorphologyEx ( src, dst, NULL, structElem2, CV_MOP_OPEN, 1);
}
