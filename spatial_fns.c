/* (c) 2004-13 Neil Ferguson, Imperial College London (neil.ferguson@imperial.ac.uk)
	All rights reserved. Copying and distribution prohibited without prior permission. */

#include "model_constants.h"
#include "model_structs.h"
#include "model_fns.h"
#include "model_globals_ext.h"
#include "binio.h"


double ExpKernel(double r2)
{
	return exp(-sqrt(r2) / P.KernelScale);
}

double PowerKernel(double r2)
{
	double t;

	t = -P.KernelShape * log(sqrt(r2) / P.KernelScale + 1);

	return (t < -690) ? 0 : exp(t);
}

double PowerKernelB(double r2)
{
	double t;

	t = 0.5 * P.KernelShape * log(r2 / (P.KernelScale * P.KernelScale));

	return (t > 690) ? 0 : (1 / (exp(t) + 1));
}

double PowerKernelUS(double r2)
{
	double t;

	t = log(sqrt(r2) / P.KernelScale + 1);

	return (t < -690) ? 0 : (exp(-P.KernelShape * t) + P.KernelP3 * exp(-P.KernelP4 * t)) / (1 + P.KernelP3);
}


double GaussianKernel(double r2)
{
	return exp(-r2 / (P.KernelScale * P.KernelScale));
}

double StepKernel(double r2)
{
	return (r2 > P.KernelScale * P.KernelScale) ? 0 : 1;
}

double PowerExpKernel(double r2)
{
	double d, t, t2;

	d = sqrt(r2);
	t = -P.KernelShape * log(d / P.KernelScale + 1);

	return (t < -690) ? 0 : exp(t - pow(d / P.KernelP3, P.KernelP4));
}

double numKernel(double r2)
{
	double t, s;

	t = r2 / P.KernelDelta;
	if (t > NKR)
	{
		fprintf(stderr, "** %lg  %lg  %lg**\n", r2, P.KernelDelta, t);
		ERR_CRITICAL("r too large in NumKernel\n");
	}
	s = t * NK_HR;
	if (s < NKR)
	{
		t = s - floor(s);
		t = (1 - t) * nKernelHR[(int)s] + t * nKernelHR[(int)(s + 1)];
	}
	else
	{
		s = t - floor(t);
		t = (1 - s) * nKernel[(int)t] + s * nKernel[(int)(t + 1)];
	}
	return t;
}


double dist2UTM(double x1, double y1, double x2, double y2)
{
	double x, y, cy1, cy2, yt, xi, yi;

	x = fabs(x1 - x2) / 2;
	y = fabs(y1 - y2) / 2;
	xi = floor(x);
	yi = floor(y);
	x -= xi;
	y -= yi;
	x = (1 - x) * sinx[(int)xi] + x * sinx[((int)xi) + 1];
	x = x * x;
	y = (1 - y) * sinx[(int)yi] + y * sinx[((int)yi) + 1];
	y = y * y;
	yt = fabs(y1 + P.SpatialBoundingBox[1]);
	yi = floor(yt);
	cy1 = yt - yi;
	cy1 = (1 - cy1) * cosx[((int)yi)] + cy1 * cosx[((int)yi) + 1];
	yt = fabs(y2 + P.SpatialBoundingBox[1]);
	yi = floor(yt);
	cy2 = yt - yi;
	cy2 = (1 - cy2) * cosx[((int)yi)] + cy2 * cosx[((int)yi) + 1];
	x = fabs(1000 * (y + x * cy1 * cy2));
	xi = floor(x);
	x -= xi;
	y = (1 - x) * asin2sqx[((int)xi)] + x * asin2sqx[((int)xi) + 1];
	return 4 * EARTHRADIUS * EARTHRADIUS * y;
}

double dist2(person* a, person* b)
{
	double x, y;

	if (P.DoUTM_coords)
		return dist2UTM(Households[a->hh].loc_x, Households[a->hh].loc_y, Households[b->hh].loc_x, Households[b->hh].loc_y);
	else
	{
		x = fabs(Households[a->hh].loc_x - Households[b->hh].loc_x);
		y = fabs(Households[a->hh].loc_y - Households[b->hh].loc_y);
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.width * 0.5) x = P.width - x;
			if (y > P.height * 0.5) y = P.height - y;
		}
		return x * x + y * y;
	}
}

double dist2_cc(cell* a, cell* b)
{
	double x, y;
	int l, m;

	l = (int)(a - Cells);
	m = (int)(b - Cells);
	if (P.DoUTM_coords)
		return dist2UTM(P.cwidth * fabs((double)(l / P.nch)), P.cheight * fabs((double)(l % P.nch)),
			P.cwidth * fabs((double)(m / P.nch)), P.cheight * fabs((double)(m % P.nch)));
	else
	{
		x = P.cwidth * fabs((double)(l / P.nch - m / P.nch));
		y = P.cheight * fabs((double)(l % P.nch - m % P.nch));
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.width * 0.5) x = P.width - x;
			if (y > P.height * 0.5) y = P.height - y;
		}
		return x * x + y * y;
	}
}

double dist2_cc_min(cell* a, cell* b)
{
	double x, y;
	int l, m, i, j;

	l = (int)(a - Cells);
	m = (int)(b - Cells);
	i = l; j = m;
	if (P.DoUTM_coords)
	{
		if (P.cwidth * ((double)abs(m / P.nch - l / P.nch)) > PI)
		{
			if (m / P.nch > l / P.nch)
				j += P.nch;
			else if (m / P.nch < l / P.nch)
				i += P.nch;
		}
		else
		{
			if (m / P.nch > l / P.nch)
				i += P.nch;
			else if (m / P.nch < l / P.nch)
				j += P.nch;
		}
		if (m % P.nch > l % P.nch)
			i++;
		else if (m % P.nch < l % P.nch)
			j++;
		return dist2UTM(P.cwidth * fabs((double)(i / P.nch)), P.cheight * fabs((double)(i % P.nch)),
			P.cwidth * fabs((double)(j / P.nch)), P.cheight * fabs((double)(j % P.nch)));
	}
	else
	{
		if ((P.DoPeriodicBoundaries) && (P.cwidth * ((double)abs(m / P.nch - l / P.nch)) > P.width * 0.5))
		{
			if (m / P.nch > l / P.nch)
				j += P.nch;
			else if (m / P.nch < l / P.nch)
				i += P.nch;
		}
		else
		{
			if (m / P.nch > l / P.nch)
				i += P.nch;
			else if (m / P.nch < l / P.nch)
				j += P.nch;
		}
		if ((P.DoPeriodicBoundaries) && (P.height * ((double)abs(m % P.nch - l % P.nch)) > P.height * 0.5))
		{
			if (m % P.nch > l % P.nch)
				j++;
			else if (m % P.nch < l % P.nch)
				i++;
		}
		else
		{
			if (m % P.nch > l % P.nch)
				i++;
			else if (m % P.nch < l % P.nch)
				j++;
		}
		x = P.cwidth * fabs((double)(i / P.nch - j / P.nch));
		y = P.cheight * fabs((double)(i % P.nch - j % P.nch));
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.width * 0.5) x = P.width - x;
			if (y > P.height * 0.5) y = P.height - y;
		}
		return x * x + y * y;
	}
}


double dist2_mm(microcell* a, microcell* b)
{
	double x, y;
	int l, m;

	l = (int)(a - Mcells);
	m = (int)(b - Mcells);
	if (P.DoUTM_coords)
		return dist2UTM(P.mcwidth * fabs((double)(l / P.nmch)), P.mcheight * fabs((double)(l % P.nmch)),
			P.mcwidth * fabs((double)(m / P.nmch)), P.mcheight * fabs((double)(m % P.nmch)));
	else
	{
		x = P.mcwidth * fabs((double)(l / P.nmch - m / P.nmch));
		y = P.mcheight * fabs((double)(l % P.nmch - m % P.nmch));
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.width * 0.5) x = P.width - x;
			if (y > P.height * 0.5) y = P.height - y;
		}
		return x * x + y * y;
	}
}


//added function to find minimum distance
double dist2_mm_min(microcell* a, microcell* b)
{
	double x, y;
	int l, m, i, j;

	l = (int)(a - Mcells);
	m = (int)(b - Mcells);
	i = l; j = m;
	if (P.DoUTM_coords)
	{
		if (P.mcwidth * ((double)abs(m / P.nmch - l / P.nmch)) > PI)
		{
			if (m / P.nmch > l / P.nmch)
				j += P.nmch;
			else if (m / P.nmch < l / P.nmch)
				i += P.nmch;
		}
		else
		{
			if (m / P.nmch > l / P.nmch)
				i += P.nmch;
			else if (m / P.nmch < l / P.nmch)
				j += P.nmch;
		}
		if (m % P.nmch > l % P.nmch)
			i++;
		else if (m % P.nmch < l % P.nmch)
			j++;
		return dist2UTM(P.mcwidth * fabs((double)(i / P.nmch)), P.mcheight * fabs((double)(i % P.nmch)),
			P.mcwidth * fabs((double)(j / P.nmch)), P.mcheight * fabs((double)(j % P.nmch)));
	}
	else
	{
		if ((P.DoPeriodicBoundaries) && (P.mcwidth * ((double)abs(m / P.nmch - l / P.nmch)) > P.width * 0.5))
		{
			if (m / P.nmch > l / P.nmch)
				j += P.nmch;
			else if (m / P.nmch < l / P.nmch)
				i += P.nmch;
		}
		else
		{
			if (m / P.nmch > l / P.nmch)
				i += P.nmch;
			else if (m / P.nmch < l / P.nmch)
				j += P.nmch;
		}
		if ((P.DoPeriodicBoundaries) && (P.height * ((double)abs(m % P.nmch - l % P.nmch)) > P.height * 0.5))
		{
			if (m % P.nmch > l % P.nmch)
				j++;
			else if (m % P.nmch < l % P.nmch)
				i++;
		}
		else
		{
			if (m % P.nmch > l % P.nmch)
				i++;
			else if (m % P.nmch < l % P.nmch)
				j++;
		}
		x = P.mcwidth * fabs((double)(i / P.nmch - j / P.nmch));
		y = P.mcheight * fabs((double)(i % P.nmch - j % P.nmch));
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.width * 0.5) x = P.width - x;
			if (y > P.height * 0.5) y = P.height - y;
		}
		return x * x + y * y;
	}
}

double dist2_raw(double ax, double ay, double bx, double by)
{
	double x, y;

	if (P.DoUTM_coords)
		return dist2UTM(ax, ay, bx, by);
	else
	{
		x = fabs(ax - bx);
		y = fabs(ay - by);
		if (P.DoPeriodicBoundaries)
		{
			if (x > P.width * 0.5) x = P.width - x;
			if (y > P.height * 0.5) y = P.height - y;
		}
		return x * x + y * y;
	}
}
