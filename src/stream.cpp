/*-----------------------------------------------------------------------*/
/* Program: STREAM                                                       */
/* Revision: $Id: stream.c,v 5.10 2013/01/17 16:01:06 mccalpin Exp mccalpin $ */
/* Original code developed by John D. McCalpin                           */
/* Programmers: John D. McCalpin                                         */
/*              Joe R. Zagar                                             */
/*                                                                       */
/* This program measures memory transfer rates in MB/s for simple        */
/* computational kernels coded in C.                                     */
/*-----------------------------------------------------------------------*/
/* Copyright 1991-2013: John D. McCalpin                                 */
/*-----------------------------------------------------------------------*/
/* License:                                                              */
/*  1. You are free to use this program and/or to redistribute           */
/*     this program.                                                     */
/*  2. You are free to modify this program for your own use,             */
/*     including commercial use, subject to the publication              */
/*     restrictions in item 3.                                           */
/*  3. You are free to publish results obtained from running this        */
/*     program, or from works that you derive from this program,         */
/*     with the following limitations:                                   */
/*     3a. In order to be referred to as "STREAM benchmark results",     */
/*         published results must be in conformance to the STREAM        */
/*         Run Rules, (briefly reviewed below) published at              */
/*         http://www.cs.virginia.edu/stream/ref.html                    */
/*         and incorporated herein by reference.                         */
/*         As the copyright holder, John McCalpin retains the            */
/*         right to determine conformity with the Run Rules.             */
/*     3b. Results based on modified source code or on runs not in       */
/*         accordance with the STREAM Run Rules must be clearly          */
/*         labelled whenever they are published.  Examples of            */
/*         proper labelling include:                                     */
/*           "tuned STREAM benchmark results"                            */
/*           "based on a variant of the STREAM benchmark code"           */
/*         Other comparable, clear, and reasonable labelling is          */
/*         acceptable.                                                   */
/*     3c. Submission of results to the STREAM benchmark web site        */
/*         is encouraged, but not required.                              */
/*  4. Use of this program or creation of derived works based on this    */
/*     program constitutes acceptance of these licensing restrictions.   */
/*  5. Absolutely no warranty is expressed or implied.                   */
/*-----------------------------------------------------------------------*/
#include "stream.h"

#include <float.h>
#include <stdlib.h>
#include <omp.h>
#include <array>


#define STREAM_ARRAY_SIZE 10000000

#define NTIMES 10

#define OFFSET 0

#define HLINE "-------------------------------------------------------------\n"

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

#ifndef STREAM_TYPE
#define STREAM_TYPE double
#endif

static std::array<STREAM_TYPE,STREAM_ARRAY_SIZE + OFFSET>  a,b,c;

static std::array<double, 4> avgtime = {0},
							maxtime = {0}, 
							mintime = {FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX};
/*
static char label[4][12] = {"Copy:      ", "Scale:     ",
						 "Add:       ", "Triad:     "};
*/

static std::array<double, 4> bytes = {
	2 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE,
	2 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE,
	3 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE,
	3 * sizeof(STREAM_TYPE) * STREAM_ARRAY_SIZE};

extern double mysecond();
extern void checkSTREAMresults();

#ifdef _OPENMP
extern int omp_get_num_threads();
#endif

/*

*/


/* void s_setup(FILE* file)
{
	// --- SETUP --- determine precision and check timing --- 

	fprintf(file, HLINE);
	fprintf(file, "STREAM version $Revision: 5.10 $\n");
	fprintf(file, HLINE);
	int BytesPerWord = sizeof(STREAM_TYPE);
	fprintf(file, "This system uses %d bytes per array element.\n", BytesPerWord);

	fprintf(file, HLINE);

	fprintf(file, "Array size = %llu (elements), Offset = %d (elements)\n", (unsigned long long)STREAM_ARRAY_SIZE, OFFSET);
	fprintf(file, "Memory per array = %.1f MiB (= %.1f GiB).\n",
			BytesPerWord * ((double)STREAM_ARRAY_SIZE / 1024.0 / 1024.0),
			BytesPerWord * ((double)STREAM_ARRAY_SIZE / 1024.0 / 1024.0 / 1024.0));
	fprintf(file, "Total memory required = %.1f MiB (= %.1f GiB).\n",
			(3.0 * BytesPerWord) * ((double)STREAM_ARRAY_SIZE / 1024.0 / 1024.),
			(3.0 * BytesPerWord) * ((double)STREAM_ARRAY_SIZE / 1024.0 / 1024. / 1024.));
	fprintf(file, "Each kernel will be executed %d times.\n", NTIMES);
	fprintf(file, " The *best* time for each kernel (excluding the first iteration)\n");
	fprintf(file, " will be used to compute the reported bandwidth.\n");
	fprintf(file, HLINE);
	fflush(file);
} */

void s_main(std::array<double,4> result)
{
	int quantum, checktick();
	int k;
	ssize_t j;
	STREAM_TYPE scalar;
	double t, times[4][NTIMES];

#ifdef _OPENMP
	printf(HLINE);
#pragma omp parallel
	{
#pragma omp master
		{
			omp_display_affinity("Thread Affinity: %0.3L %.8n %.15{thread_affinity} %.12H");
			k = omp_get_num_threads();
			printf("Number of Threads requested = %i\n", k);
		}
	}
#endif

#ifdef _OPENMP
	k = 0;
#pragma omp parallel
#pragma omp critical
	k++;
	printf("Number of Threads counted = %i\n", k);
#endif

	/* Get initial value for system clock. */
#pragma omp parallel for
	for (j = 0; j < STREAM_ARRAY_SIZE; j++)
	{
		a[j] = 1.0;
		b[j] = 2.0;
		c[j] = 0.0;
	}

	// printf(HLINE);

	if ((quantum = checktick()) < 1)
	{
		quantum = 1;
	}

	t = mysecond();
#pragma omp parallel for
	for (j = 0; j < STREAM_ARRAY_SIZE; j++)
		a[j] = 2.0E0 * a[j];
	t = 1.0E6 * (mysecond() - t);

	if ((t / quantum) <= 20)
	{
		printf("Increase size of the arrays\n");
		return;
	}

	/*	--- MAIN LOOP --- repeat test cases NTIMES times --- */

	scalar = 3.0;
	for (k = 0; k < NTIMES; k++)
	{
		times[0][k] = mysecond();
// Copy
#pragma omp parallel for
		for (j = 0; j < STREAM_ARRAY_SIZE; j++)
			c[j] = a[j];

		times[0][k] = mysecond() - times[0][k];

		times[1][k] = mysecond();
// Scale
#pragma omp parallel for
		for (j = 0; j < STREAM_ARRAY_SIZE; j++)
			b[j] = scalar * c[j];

		times[1][k] = mysecond() - times[1][k];

		times[2][k] = mysecond();
// Add
#pragma omp parallel for
		for (j = 0; j < STREAM_ARRAY_SIZE; j++)
			c[j] = a[j] + b[j];

		times[2][k] = mysecond() - times[2][k];

		times[3][k] = mysecond();
// Triad
#pragma omp parallel for
		for (j = 0; j < STREAM_ARRAY_SIZE; j++)
			a[j] = b[j] + scalar * c[j];

		times[3][k] = mysecond() - times[3][k];
	}

	/*	--- SUMMARY --- */

	for (k = 1; k < NTIMES; k++) /* note -- skip first iteration */
	{
		for (j = 0; j < 4; j++)
		{
			avgtime[j] = avgtime[j] + times[j][k];
			mintime[j] = MIN(mintime[j], times[j][k]);
			maxtime[j] = MAX(maxtime[j], times[j][k]);
		}
	}

	// printf("Function    Best Rate MB/s  Avg time     Min time     Max time\n");
	for (j = 0; j < 4; j++)
	{
		avgtime[j] = avgtime[j] / (double)(NTIMES - 1);
		result[j] = 1.0E-06 * bytes[j] / mintime[j];

				printf("%s%12.1f  %11.6f  %11.6f  %11.6f\n", "id",//Need Id
					   1.0E-06 * bytes[j] / mintime[j],
					   avgtime[j],
					   mintime[j],
					   maxtime[j]);
		
	}
	// printf(HLINE);

	/* --- Check Results --- */
	checkSTREAMresults();
}

#define M 20

int checktick()
{
	int i, minDelta, Delta;
	double t1, t2, timesfound[M];

	/*  Collect a sequence of M unique time values from the system. */

	for (i = 0; i < M; i++)
	{
		t1 = mysecond();
		while (((t2 = mysecond()) - t1) < 1.0E-6)
			;
		timesfound[i] = t1 = t2;
	}

	/*
	 * Determine the minimum difference between these M values.
	 * This result will be our estimate (in microseconds) for the
	 * clock granularity.
	 */

	minDelta = 1000000;
	for (i = 1; i < M; i++)
	{
		Delta = (int)(1.0E6 * (timesfound[i] - timesfound[i - 1]));
		minDelta = MIN(minDelta, MAX(Delta, 0));
	}

	return (minDelta);
}

/* A gettimeofday routine to give access to the wall
   clock timer on most UNIX-like systems.  */

#include <sys/time.h>

double mysecond()
{
	struct timeval tp;
	struct timezone tzp;

	gettimeofday(&tp, &tzp);
	return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

#ifndef abs
#define abs(a) ((a) >= 0 ? (a) : -(a))
#endif
void checkSTREAMresults()
{
	STREAM_TYPE aj, bj, cj, scalar;
	STREAM_TYPE aSumErr, bSumErr, cSumErr;
	STREAM_TYPE aAvgErr, bAvgErr, cAvgErr;
	double epsilon;
	ssize_t j;
	int k, ierr, err;

	/* reproduce initialization */
	aj = 1.0;
	bj = 2.0;
	cj = 0.0;
	/* a[] is modified during timing check */
	aj = 2.0E0 * aj;
	/* now execute timing loop */
	scalar = 3.0;
	for (k = 0; k < NTIMES; k++)
	{
		cj = aj;
		bj = scalar * cj;
		cj = aj + bj;
		aj = bj + scalar * cj;
	}

	/* accumulate deltas between observed and expected results */
	aSumErr = 0.0;
	bSumErr = 0.0;
	cSumErr = 0.0;
	for (j = 0; j < STREAM_ARRAY_SIZE; j++)
	{
		aSumErr += abs(a[j] - aj);
		bSumErr += abs(b[j] - bj);
		cSumErr += abs(c[j] - cj);
		// if (j == 417) printf("Index 417: c[j]: %f, cj: %f\n",c[j],cj);	// MCCALPIN
	}
	aAvgErr = aSumErr / (STREAM_TYPE)STREAM_ARRAY_SIZE;
	bAvgErr = bSumErr / (STREAM_TYPE)STREAM_ARRAY_SIZE;
	cAvgErr = cSumErr / (STREAM_TYPE)STREAM_ARRAY_SIZE;

	if (sizeof(STREAM_TYPE) == 4)
	{
		epsilon = 1.e-6;
	}
	else if (sizeof(STREAM_TYPE) == 8)
	{
		epsilon = 1.e-13;
	}
	else
	{
		printf("WEIRD: sizeof(STREAM_TYPE) = %lu\n", sizeof(STREAM_TYPE));
		epsilon = 1.e-6;
	}

	err = 0;
	if (abs(aAvgErr / aj) > epsilon)
	{
		err++;
		printf("Failed Validation on array a[], AvgRelAbsErr > epsilon (%e)\n", epsilon);
		printf("     Expected Value: %e, AvgAbsErr: %e, AvgRelAbsErr: %e\n", aj, aAvgErr, abs(aAvgErr) / aj);
		ierr = 0;
		for (j = 0; j < STREAM_ARRAY_SIZE; j++)
		{
			if (abs(a[j] / aj - 1.0) > epsilon)
			{
				ierr++;
#ifdef VERBOSE
				if (ierr < 10)
				{
					printf("         array a: index: %ld, expected: %e, observed: %e, relative error: %e\n",
						   j, aj, a[j], abs((aj - a[j]) / aAvgErr));
				}
#endif
			}
		}
		printf("     For array a[], %d errors were found.\n", ierr);
	}
	if (abs(bAvgErr / bj) > epsilon)
	{
		err++;
		printf("Failed Validation on array b[], AvgRelAbsErr > epsilon (%e)\n", epsilon);
		printf("     Expected Value: %e, AvgAbsErr: %e, AvgRelAbsErr: %e\n", bj, bAvgErr, abs(bAvgErr) / bj);
		printf("     AvgRelAbsErr > Epsilon (%e)\n", epsilon);
		ierr = 0;
		for (j = 0; j < STREAM_ARRAY_SIZE; j++)
		{
			if (abs(b[j] / bj - 1.0) > epsilon)
			{
				ierr++;
#ifdef VERBOSE
				if (ierr < 10)
				{
					printf("         array b: index: %ld, expected: %e, observed: %e, relative error: %e\n",
						   j, bj, b[j], abs((bj - b[j]) / bAvgErr));
				}
#endif
			}
		}
		printf("     For array b[], %d errors were found.\n", ierr);
	}
	if (abs(cAvgErr / cj) > epsilon)
	{
		err++;
		printf("Failed Validation on array c[], AvgRelAbsErr > epsilon (%e)\n", epsilon);
		printf("     Expected Value: %e, AvgAbsErr: %e, AvgRelAbsErr: %e\n", cj, cAvgErr, abs(cAvgErr) / cj);
		printf("     AvgRelAbsErr > Epsilon (%e)\n", epsilon);
		ierr = 0;
		for (j = 0; j < STREAM_ARRAY_SIZE; j++)
		{
			if (abs(c[j] / cj - 1.0) > epsilon)
			{
				ierr++;
#ifdef VERBOSE
				if (ierr < 10)
				{
					printf("         array c: index: %ld, expected: %e, observed: %e, relative error: %e\n",
						   j, cj, c[j], abs((cj - c[j]) / cAvgErr));
				}
#endif
			}
		}
		printf("     For array c[], %d errors were found.\n", ierr);
	}
	if (err == 0)
	{
		printf("Solution Validates: avg error less than %e on all three arrays\n", epsilon);
	}
#ifdef VERBOSE
	// printf ("Results Validation Verbose Results: \n");
	// printf ("    Expected a(1), b(1), c(1): %f %f %f \n",aj,bj,cj);
	// printf ("    Observed a(1), b(1), c(1): %f %f %f \n",a[1],b[1],c[1]);
	printf("    Rel Errors on a, b, c:     %e %e %e \n", abs(aAvgErr / aj), abs(bAvgErr / bj), abs(cAvgErr / cj));
#endif
}