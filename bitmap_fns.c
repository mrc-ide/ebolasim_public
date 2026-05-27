/* (c) 2004-13 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
	All rights reserved. Copying and distribution prohibited without prior permission. */

#include "model_constants.h"
#include "model_structs.h"
#include "model_fns.h"
#include "model_globals_ext.h"
#include "binio.h"

#define PREVSCALE 50
#define IMM_PREVSCALE 5
#define BWCOLS 58
#define PREVCOLS 87

void CaptureBitmap(int ns, int tp)
{
	int i, j, x, y, f, mi;
	static double mx;
	static int fst = 1;
	double prev, age;

	mi = (int)(P.bwidth * P.bheight);
	if (fst)
	{
		fst = 0;
		mx = 0;
		for (i = 0; i < mi; i++) bmi[i] = 0;
		for (i = 0; i < P.N; i++)
		{
			x = ((int)(Households[Hosts[i].hh].loc_x * P.scalex)) - P.bminx;
			y = ((int)(Households[Hosts[i].hh].loc_y * P.scaley)) - P.bminy;
			if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
			{
				j = y * bmh->width + x;
				if ((j < bmh->imagesize) && (j >= 0))
				{
					bmi[j]++;
					if (bmi[j] > mx) mx = (double)bmi[j];
				}
			}
		}
		mx = log(1.001 * mx);
		for (i = 0; i < P.NMC; i++)
		{
			f = 0;
			if ((i / P.nmch == (i + 1) / P.nmch) && ((Mcells[i].country != Mcells[i + 1].country) || ((P.DoAdunitBoundaryOutput) && ((AdUnits[Mcells[i].adunit].id % P.AdunitLevel1Mask) / P.AdunitBitmapDivisor != (AdUnits[Mcells[i + 1].adunit].id % P.AdunitLevel1Mask) / P.AdunitBitmapDivisor))))
			{
				f = 1;
				//fprintf(stderr,"mcells %i, %i\n",i,i+1);
				//fprintf(stderr,"mcells countries %i, %i\n",Mcells[i].country,Mcells[i+1].country);
				//fprintf(stderr,"mcells adunits %i, %i\n",Mcells[i].adunit,Mcells[i+1].adunit);
				//fprintf(stderr,"mcells pop %i, %i\n",Mcells[i].n,Mcells[i+1].n);
			}
			else if ((i > 0) && (i / P.nmch == (i - 1) / P.nmch) && (Mcells[i].country != Mcells[i - 1].country)) f = 1;
			else if ((i < P.NMC - P.nmch) && ((Mcells[i].country != Mcells[i + P.nmch].country) || ((P.DoAdunitBoundaryOutput) && ((AdUnits[Mcells[i].adunit].id % P.AdunitLevel1Mask) / P.AdunitBitmapDivisor != (AdUnits[Mcells[i + P.nmch].adunit].id % P.AdunitLevel1Mask) / P.AdunitBitmapDivisor)))) f = 1;
			else if ((i > P.nmch) && (Mcells[i].country != Mcells[i - P.nmch].country)) f = 1;
			if (f)
			{
				x = (int)(P.mcwidth * (((double)(i / P.nmch)) + 0.5) * P.scalex) - P.bminx;
				y = (int)(P.mcheight * (((double)(i % P.nmch)) + 0.5) * P.scaley) - P.bminy;
				if ((x >= 0) && (x < P.bwidth) && (y >= 0) && (y < P.bheight))
				{
					j = y * bmh->width + x;
					if ((j < bmh->imagesize) && (j >= 0)) bmi[j] = -1;
				}
			}
		}
		for (i = 0; i < P.bwidth / 2; i++)
		{
			prev = floor(1.99999 * ((double)i) * PREVCOLS / ((double)P.bwidth) * 2);
			f = BWCOLS + ((int)prev);
			for (j = 0; j < 10; j++)
			{
				bm[(j + P.bheight + 5) * bmh->width + P.bwidth / 4 + i] = f;
			}
		}
	}
#pragma omp parallel for private(i,j,prev) schedule(static,5000)
	for (i = 0; i < mi; i++)
	{
		if (bmi[i] == -1)
			bm[i] = BWCOLS - 1; /* black for country boundary */
		else if ((bmi4[i] == 0) && (bmi2[i] == 0) && (bmi3[i] > 0))
			bm[i] = (unsigned char)(3 * BWCOLS + BWCOLS * log(bmi3[i]) / mx);  /* green for recovered */
		else if (bmi2[i] > 0)
			bm[i] = (unsigned char)(BWCOLS + BWCOLS * log(bmi2[i]) / mx); /* red for infected */
		else if (bmi4[i] > 0)
			bm[i] = (unsigned char)(2 * BWCOLS + BWCOLS * log(bmi[i]) / mx); /* blue for treated */
		else if (bmi[i] > 0)
			bm[i] = (unsigned char)(BWCOLS * log(bmi[i]) / mx); /* grey for just people */
		else
			bm[i] = 0;
	}
}


//void CaptureBitmap(int ns,int tp)
//{
//	int i,j,x,y,f,mi;
//	static double mx;
//	static int fst=1;
//	double prev,age;
//
//	mi=(int) (P.bwidth*P.bheight);
//	if(fst)
//		{
//		fst=0;
//		mx=0;
//		for(i=0;i<mi;i++) bmi[i]=0;
//		for(i=0;i<P.N;i++)
//			{
//			x=((int) (Households[Hosts[i].hh].loc_x*P.scalex))-P.bminx;
//			y=((int) (Households[Hosts[i].hh].loc_y*P.scaley))-P.bminy;
//			if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
//				{
//				j=y*bmh->width+x;
//				if((j<bmh->imagesize)&&(j>=0))
//					{
//					bmi[j]++;
//					if(bmi[j]>mx) mx=(double) bmi[j];
//					}
//				}
//			}
//		mx=log(1.001*mx);
//		for(i=0;i<P.NMC;i++)
//			{
//			f=0;
//			if((i/P.nmch==(i+1)/P.nmch)&&((Mcells[i].country!=Mcells[i+1].country)||((P.DoAdunitBoundaryOutput)&&((AdUnits[Mcells[i].adunit].id%P.AdunitLevel1Mask)/P.AdunitBitmapDivisor !=(AdUnits[Mcells[i+1].adunit].id%P.AdunitLevel1Mask)/P.AdunitBitmapDivisor)))) f=1;
//			else if((i>0)&&(i/P.nmch==(i-1)/P.nmch)&&(Mcells[i].country!=Mcells[i-1].country)) f=1;
//			else if((i<P.NMC-P.nmch)&&((Mcells[i].country!=Mcells[i+P.nmch].country)||((P.DoAdunitBoundaryOutput)&&((AdUnits[Mcells[i].adunit].id%P.AdunitLevel1Mask)/P.AdunitBitmapDivisor!=(AdUnits[Mcells[i+P.nmch].adunit].id%P.AdunitLevel1Mask)/P.AdunitBitmapDivisor)))) f=1;
//			else if((i>P.nmch)&&(Mcells[i].country!=Mcells[i-P.nmch].country)) f=1;
//			if(f)
//				{
//				x=(int) (P.mcwidth*(((double) (i/P.nmch))+0.5)*P.scalex)-P.bminx;
//				y=(int) (P.mcheight*(((double) (i%P.nmch))+0.5)*P.scaley)-P.bminy;
//				if((x>=0)&&(x<P.bwidth)&&(y>=0)&&(y<P.bheight))
//					{
//					j=y*bmh->width+x;
//					if((j<bmh->imagesize)&&(j>=0)) bmi[j]=-1;
//					}
//				}
//			}
//		for(i=0;i<P.bwidth/2;i++)
//			{
//			prev=floor(1.99999*((double) i)*PREVCOLS/((double) P.bwidth)*2);
//			f=BWCOLS+((int) prev);
//			for(j=0;j<10;j++)
//				{
//				bm[(j+P.bheight+5)*bmh->width+P.bwidth/4+i]=f;
//				}
//			}
//		}
//#pragma omp parallel for private(i,j,prev) schedule(static,5000)
//	for(i=0;i<mi;i++)
//		{
//		if(bmi[i]==-1)
//			bm[i]=BWCOLS-1; /* black for country boundary */
//		else if((bmi4[i]==0)&&(bmi2[i]==0)&&(bmi3[i]>0))
//			bm[i]=(unsigned char) (3*BWCOLS+BWCOLS*log(bmi3[i])/mx);  /* green for recovered */
//		else if(bmi2[i]>0)
//			bm[i]=(unsigned char) (BWCOLS+BWCOLS*log(bmi2[i])/mx); /* red for infected */
//		else if(bmi4[i]>0)
//			bm[i]=(unsigned char) (2*BWCOLS+BWCOLS*log(bmi[i])/mx); /* blue for treated */
//		else if(bmi[i]>0)
//			bm[i]=(unsigned char) (BWCOLS*log(bmi[i])/mx); /* grey for just people */
//		}
//}

void CaptureMeanBitmap(int ns)
{
	int i, j, mi;
	static double mx;
	static int fst = 1;

	mi = (int)bmh->imagesize;
	if (fst)
	{
		fst = 0;
		mx = 0;
		for (j = 0; j < bmh->imagesize; j++)
			if (bmi[j] > mx) mx = bmi[j];
		mx = log(1.001 * mx);
	}
#pragma omp parallel for private(i) schedule(static,500) //added i to private
	for (i = 0; i < mi; i++)
	{
		bmi2[i] = ceil(TSMean[ns].bmi2[i] / ((float)P.NRactual));
		bmi3[i] = ceil(TSMean[ns].bmi3[i] / ((float)P.NRactual));
		bmi4[i] = ceil(TSMean[ns].bmi4[i] / ((float)P.NRactual));
	}
#pragma omp parallel for private(i) schedule(static,500) //added i to private
	for (i = 0; i < mi; i++)
	{
		if (bmi[i] == -1)
			bm[i] = BWCOLS - 1; /* black for country boundary */
		else if ((bmi4[i] == 0) && (bmi2[i] == 0) && (bmi3[i] > 0))
			bm[i] = (unsigned char)(3 * BWCOLS + BWCOLS * log(bmi3[i]) / mx);  /* green for recovered */
		else if (bmi2[i] > 0)
			bm[i] = (unsigned char)(BWCOLS + BWCOLS * log(bmi2[i]) / mx); /* red for infected */
		else if (bmi4[i] > 0)
			bm[i] = (unsigned char)(2 * BWCOLS + BWCOLS * log(bmi4[i]) / mx); /* blue for treated */
		else if (bmi[i] > 0)
			bm[i] = (unsigned char)(BWCOLS * log(bmi[i]) / mx); /* grey */
	}
}


void OutputBitmap(double t2, int tp)
{
	FILE* dat;
	char buf[1024], OutF[1024];
	int i, j;
	static int cn1 = 0, cn2 = 0, cn3 = 0, cn4 = 0;
	size_t a;

	if (tp == 0)
	{
		j = cn1;
		cn1++;
		sprintf(OutF, "%s", OutFile);
	}
	else if (tp == 1)
	{
		j = cn2;
		cn2++;
		sprintf(OutF, "Mean.%s", OutFile);
	}
	else if (tp == 2)
	{
		j = cn3;
		cn3++;
		sprintf(OutF, "Min.%s", OutFile);
	}
	else if (tp == 3)
	{
		j = cn4;
		cn4++;
		sprintf(OutF, "Max.%s", OutFile);
	}

#ifdef IMAGE_MAGICK
	using namespace Magick;
	fprintf(stderr, "\noutputing ImageMagick stuff");
	sprintf(buf, "%s.bmp", OutF);
	if (!(dat = fopen(buf, "wb"))) ERR_CRITICAL("Unable to open bitmap file\n");
	fprintf(dat, "BM");
	//fwrite_big((void *) &bmf,sizeof(unsigned char),(sizeof(bitmap_header)/sizeof(unsigned char))+bmh->imagesize,dat);
	fwrite_big((void*)bmf, sizeof(bitmap_header), 1, dat);
	for (i = 0; i < bmh->imagesize; i++) fputc(bm[i], dat);
	fclose(dat);
	Image bmap(buf);
	sprintf(buf, "%s.%d.png", OutF, j);
	ColorRGB white(1.0, 1.0, 1.0);
	bmap.transparent(white);
	bmap.write(buf);
#elif defined(WIN32_BM)	
	//Windows specific bitmap manipulation code - could be recoded using LIBGD or another unix graphics library
	using namespace Gdiplus;

	HBITMAP hbm;
	wchar_t wbuf[1024];
	static UINT palsize;
	static ColorPalette* palette;

	//Add new bitmap to AVI
	if ((P.OutputBitmap == 1) && (tp == 0)) AddAviFrame(avi, bmpdib, (unsigned char*)(&bmh->palette[0][0]));

	//This transfers HBITMAP to GDI+ Bitmap object
	Bitmap* gdip_bmp = Bitmap::FromHBITMAP(bmpdib, NULL);
	//Now change White in palette (first entry) to be transparent
	if ((cn1 == 1) && (tp == 0))
	{
		palsize = gdip_bmp->GetPaletteSize();
		palette = (ColorPalette*)malloc(palsize);
	}
	i = gdip_bmp->GetPalette(palette, palsize);
	palette->Flags = PaletteFlagsHasAlpha;
	palette->Entries[0] = 0x00ffffff; // Transparent white 
	gdip_bmp->SetPalette(palette);
	//Now save as png
	sprintf(buf, "%s.%05i.png", OutF, j + 1); //sprintf(buf,"%s.ge\\%s.%05i.png",OutFileBase,OutF,j+1);
	mbstowcs_s(&a, wbuf, strlen(buf) + 1, buf, _TRUNCATE);
	gdip_bmp->Save(wbuf, &encoderClsid, NULL);
	delete gdip_bmp;
#else
	sprintf(buf, "%s.bmp", OutF);
	if (!(dat = fopen(buf, "wb"))) ERR_CRITICAL("Unable to open bitmap file\n");
	fprintf(dat, "BM");
	fwrite_big((void*)bmf, sizeof(unsigned char), sizeof(bitmap_header) / sizeof(unsigned char) + bmh->imagesize, dat);
	fclose(dat);
#endif
}

void InitBMHead()
{
	int i, j, k, k2, R, B, G;
	double x, f;

	fprintf(stderr, "Initialising bitmap\n");
	k = 4 * P.bwidth * P.bheight2;
	k2 = sizeof(bitmap_header) / sizeof(unsigned char);

	if (!(bmf = (unsigned char*)malloc((k + k2) * sizeof(unsigned char))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");
	bm = &(bmf[k2]);
	bmp = &(bmf[12]);
	bmh = (bitmap_header*)bmf;
	bmh->spare = 0;
	bmh->boffset = 2 + sizeof(bitmap_header);
	bmh->headersize = 40;
	bmh->width = P.bwidth;
	bmh->height = P.bheight2;
	bmh->PlanesAndBitspp = 8 * 65536 + 1;
	bmh->compr = 0;
	bmh->imagesize = bmh->width * bmh->height;
	bmh->filesize = 2 + bmh->imagesize + ((unsigned int)sizeof(bitmap_header));
	bmh->hres = bmh->vres = (int)(bmh->width * 10);
	bmh->colours = BWCOLS * 4;
	bmh->impcol = 0;
	for (i = 0; i < 256; i++)
		bmh->palette[i][3] = 0;
	for (j = 0; j < BWCOLS; j++)
	{
		B = 255 - 255 * j / (BWCOLS - 1);
		bmh->palette[j][0] = bmh->palette[j][1] = bmh->palette[j][2] = (unsigned char)B;
		bmh->palette[BWCOLS + j][0] = 0;
		bmh->palette[BWCOLS + j][1] = 0;
		bmh->palette[BWCOLS + j][2] = (unsigned char)B;
		bmh->palette[2 * BWCOLS + j][0] = (unsigned char)B;
		bmh->palette[2 * BWCOLS + j][1] = 0;
		bmh->palette[2 * BWCOLS + j][2] = 0;
		bmh->palette[3 * BWCOLS + j][0] = 0;
		bmh->palette[3 * BWCOLS + j][1] = (unsigned char)B;
		bmh->palette[3 * BWCOLS + j][2] = 0;
	}
	if (!(bmi = (float*)malloc(bmh->imagesize * sizeof(float))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");
	if (!(bmi2 = (float*)malloc(bmh->imagesize * sizeof(float))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");
	if (!(bmi3 = (float*)malloc(bmh->imagesize * sizeof(float))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");
	if (!(bmi4 = (float*)malloc(bmh->imagesize * sizeof(float))))
		ERR_CRITICAL("Unable to allocate storage for bitmap\n");

#ifdef WIN32_BM
	bmpdib = CreateDIBSection(GetDC(NULL), (BITMAPINFO*)bmp, DIB_RGB_COLORS, (void**)&bm, NULL, NULL);
	Gdiplus::GdiplusStartupInput gdiplusStartupInput;
	Gdiplus::GdiplusStartup(&m_gdiplusToken, &gdiplusStartupInput, NULL);

	using namespace Gdiplus;
	UINT  num = 0;          // number of image encoders
	UINT  size = 0;         // size of the image encoder array in bytes

	ImageCodecInfo* pImageCodecInfo = NULL;
	GetImageEncodersSize(&num, &size);
	pImageCodecInfo = (ImageCodecInfo*)(malloc(size));
	GetImageEncoders(num, size, pImageCodecInfo);
	for (UINT j = 0; j < num; ++j)
	{
		if (wcscmp(pImageCodecInfo[j].MimeType, L"image/png") == 0)
		{
			encoderClsid = pImageCodecInfo[j].Clsid;
			j = num;
		}
	}
	free(pImageCodecInfo);
	char buf[1024];
	sprintf(buf, "%s.ge", OutFileBase);
	if (!(CreateDirectory(buf, NULL))) fprintf(stderr, "Unable to create directory %s\n", buf);

#endif


}


void HSB2RGB(double h1, double s1, double v1, int* r, int* g, int* b)
{
	int i;
	double  h, s, v, p, q, t;
	double f;

	h = h1; s = s1; v = v1;
	if (s == 0)
	{
		*r = *g = *b = v * 255;  /* If sat=0 the color is grey as determined by value */
		return;
	}
	else
	{
		h *= 6;
		if (h >= 6)  h = 0;
		i = (int)h;
		f = h - floor(h);
		p = 255 * v * (1 - s);
		q = 255 * v * (1 - s * f);
		t = 255 * v * (1 - (s * (1 - f)));
		switch (i)
		{
		case 0: *r = v * 255 + 0.5;     /* Add 0.5 to all values in order to round */
			*g = t + 0.5;         /* instead of truncating */
			*b = p + 0.5;
			return;
		case 1: *r = q + 0.5;
			*g = v * 255 + 0.5;
			*b = p + 0.5;
			return;
		case 2: *r = p + 0.5;
			*g = v * 255 + 0.5;
			*b = t + 0.5;
			return;
		case 3: *r = p + 0.5;
			*g = q + 0.5;
			*b = v * 255 + 0.5;
			return;
		case 4: *r = t + 0.5;
			*g = p + 0.5;
			*b = v * 255 + 0.5;
			return;
		case 5: *r = v * 255 + 0.5;
			*g = p + 0.5;
			*b = q + 0.5;
			return;
		}
	}

}
