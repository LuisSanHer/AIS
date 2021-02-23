/*
 * rand.c
 *
 *  Created on: 24/01/2011
 *      Author: Saúl Zapotecas
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rand.h"

double oldrand[55];
int jrand;

/* Get seed number for random and start it up */
void randomize(double seed)
{
	int j1;
	for (j1 = 0; j1 <= 54; j1++)
	{
		oldrand[j1] = 0.0;
	}
	jrand = 0;
	warmup_random(seed);
	return;
}

/* Get randomize off and running */
void warmup_random(double seed)
{
	int j1, ii;
	double new_random, prev_random;
	oldrand[54] = seed;
	new_random = 0.000000001;
	prev_random = seed;
	for (j1 = 1; j1 <= 54; j1++)
	{
		ii = (21 * j1) % 54;
		oldrand[ii] = new_random;
		new_random = prev_random - new_random;
		if (new_random < 0.0)
		{
			new_random += 1.0;
		}
		prev_random = oldrand[ii];
	}
	advance_random();
	advance_random();
	advance_random();
	jrand = 0;
	return;
}

/* Create next batch of 55 random numbers */
void advance_random()
{
	int j1;
	double new_random;
	for (j1 = 0; j1 < 24; j1++)
	{
		new_random = oldrand[j1] - oldrand[j1 + 31];
		if (new_random < 0.0)
		{
			new_random = new_random + 1.0;
		}
		oldrand[j1] = new_random;
	}
	for (j1 = 24; j1 < 55; j1++)
	{
		new_random = oldrand[j1] - oldrand[j1 - 24];
		if (new_random < 0.0)
		{
			new_random = new_random + 1.0;
		}
		oldrand[j1] = new_random;
	}
}

/* Fetch a single random number between 0.0 and 1.0 */
double randomperc()
{
	jrand++;
	if (jrand >= 55)
	{
		jrand = 1;
		advance_random();
	}
	return ((double) oldrand[jrand]);
}

/* Fetch a single random integer between low and high including the bounds */
int rndint(int low, int high)
{
	int res;
	if (low >= high)
	{
		res = low;
	}
	else
	{
		res = low + (randomperc() * (high - low + 1));
		if (res > high)
		{
			res = high;
		}
	}
	return (res);
}

/* Fetch a single random real number between low and high including the bounds */
double rndreal(double low, double high)
{
	return (low + (high - low) * randomperc());
}

double box_muller(double m, double s) /* normal random variate generator */
{
	/* mean m, standard deviation s */
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last) /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do
		{
			x1 = 2.0 * randomperc() - 1.0;
			x2 = 2.0 * randomperc() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while (w >= 1.0);

		w = sqrt((-2.0 * log(w)) / w);
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}
	return (m + y1 * s);
}

double Gauss(double sigma)
{
	/* System generated locals */
	double ret_val;
	/* Local variables */
	static double u, x, y, u0, u1, u2;
	/* 	SIGMA	--> standard deviation */
	/* L1: */
	u = randomperc();
	u0 = randomperc();
	if (u >= .919544406)
	{
		goto L2;
	}
	x = (u0 + u * .825339283) * 2.40375766 - 2.11402808;
	goto L10;
	L2: if (u < .965487131)
	{
		goto L4;
	}
	L3: u1 = randomperc();
	y = sqrt(4.46911474 - log(u1) * 2.);
	u2 = randomperc();
	if (y * u2 > 2.11402808)
	{
		goto L3;
	}
	goto L9;
	L4: if (u < .949990709)
	{
		goto L6;
	}
	L5: u1 = randomperc();
	y = u1 * .273629336 + 1.84039875;
	u2 = randomperc();
	if (exp(y * -.5 * y) * .39894228 - .443299126 + y * .209694057
			< u2 * .0427025816)
	{
		goto L5;
	}
	goto L9;
	L6: if (u < .925852334)
	{
		goto L8;
	}
	L7: u1 = randomperc();
	y = u1 * 1.55066917 + .289729574;
	u2 = randomperc();
	if (exp(y * -.5 * y) * .39894228 - .443299126 + y * .209694057
			< u2 * .0159745227)
	{
		goto L7;
	}
	goto L9;
	L8: u1 = randomperc();
	y = u1 * .289729574;
	u2 = randomperc();
	if (exp(y * -.5 * y) * .39894228 - .382544556 < u2 * .0163977244)
	{
		goto L8;
	}
	L9: x = y;
	if (u0 >= .5)
	{
		x = -y;
	}
	L10: ret_val = sigma * x;
	return ret_val;
}

double N(double m, double sigma)
{
	return box_muller(m, sigma);
	return (m + Gauss(sigma));
}

void shuffle(int *a, int n)
{
	int *b = (int*) malloc(sizeof(int) * n);
	int i, rnd, tmp;
	int *aa = (int*) a;
	memcpy(b, a, n * sizeof(int));
	for (i = 0; i < n; i++)
	{
		rnd = rndint(i, n - 1);
		tmp = b[rnd];
		b[rnd] = b[i];
		aa[i] = tmp;
	}
	free(b);
}
