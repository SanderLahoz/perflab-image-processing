/*******************************************************************
 * 
 * driver.c - Driver program for CS:APP Performance Lab
 * 
 * In kernels.c, students generate an arbitrary number of rotate and
 * smooth test functions, which they then register with the driver
 * program using the add_rotate_function() and add_smooth_function()
 * functions.
 * 
 * The driver program runs and measures the registered test functions
 * and reports their performance.
 * 
 * Copyright (c) 2002, R. Bryant and D. O'Hallaron, All rights
 * reserved.  May not be used, modified, or copied without permission.
 *
 ********************************************************************/

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include "fcyc.h"
#include "defs.h"
#include "config.h"

/* Structure that identifies the student */
extern student_t student; 

/* Keep track of a number of different test functions */
#define MAX_BENCHMARKS 100
#define DIM_CNT 4

/* Misc constants */
#define BSIZE 64     /* cache block size in bytes */ // W: s/32/64 `cat /proc/cpuinfo`
#define MAX_DIM 8192 /* 1024 + 256 */ // W: bigger arrays now, removed the unexplained +256
#define ODD_DIM 96   /* not a power of 2 */

/* Margin of error for blend calculation */
#define EPSILON 5 // a very generous error margin

/* fast versions of min and max */
#define min(a,b) (a < b ? a : b)
#define max(a,b) (a > b ? a : b)

/* This struct characterizes the results for one benchmark test */
typedef struct {
    lab_test_func tfunct; /* The test function */
    double cpes[DIM_CNT]; /* One CPE result for each dimension */
    char *description;    /* ASCII description of the test function */
    unsigned short valid; /* The function is tested if this is non zero */
} bench_t;

/* The range of image dimensions that we will be testing */
static int test_dim_rotate[] = {64, 128, 256, 512};
static int test_dim_blend[] = {64, 128, 256, 512};
static int test_dim_smooth[] = {64, 128, 256, 512};

/* Baseline CPEs (see config.h) */
static double rotate_baseline_cpes[] = {R64, R128, R256, R512};
static double blend_baseline_cpes[] = {B64, B128, B256, B512};
static double smooth_baseline_cpes[] = {S64, S128, S256, S512};

/* These hold the results for all benchmarks */
static bench_t benchmarks_rotate[MAX_BENCHMARKS];
static bench_t benchmarks_rotate_t[MAX_BENCHMARKS];
static bench_t benchmarks_blend[MAX_BENCHMARKS];
static bench_t benchmarks_blend_v[MAX_BENCHMARKS];
static bench_t benchmarks_smooth[MAX_BENCHMARKS];
static bench_t benchmarks_smooth_n[MAX_BENCHMARKS];

/* These give the sizes of the above lists */
static int rotate_benchmark_count = 0;
static int rotate_t_benchmark_count = 0;
static int blend_benchmark_count = 0;
static int blend_v_benchmark_count = 0;
static int smooth_benchmark_count = 0;
static int smooth_n_benchmark_count = 0;

/* 
 * An image is a dimxdim matrix of pixels stored in a 1D array.  The
 * data array holds three images (the input original, a copy of the original, 
 * and the output result array. There is also an additional BSIZE bytes
 * of padding for alignment to cache block boundaries.
 */
static pixel data[(3*MAX_DIM*MAX_DIM) + (BSIZE/sizeof(pixel))];

/* Various image pointers */
static pixel *orig = NULL;         /* original image */
static pixel *copy_of_orig = NULL; /* copy of original for checking result */
static pixel *result = NULL;       /* result image */

static int benchmark_only = 0;

/* Keep track of the best blend, rotate and smooth score for grading */
static double rotate_maxmean = 0.0;
static char *rotate_maxmean_desc = NULL;

static double rotate_t_maxmean = 0.0;
static char *rotate_t_maxmean_desc = NULL;

static double blend_maxmean = 0.0;
static char *blend_maxmean_desc = NULL;

static double blend_v_maxmean = 0.0;
static char *blend_v_maxmean_desc = NULL;

static double smooth_maxmean = 0.0;
static char *smooth_maxmean_desc = NULL;

static double smooth_n_maxmean = 0.0;
static char *smooth_n_maxmean_desc = NULL;

/* background color for blend functions */
pixel bgc; // non-static to make this a global variable accessible from kernel.c
static pixel copy_of_bgc; // backup

/******************** Functions begin *************************/

void add_rotate_function(lab_test_func f, char *description) 
{
    benchmarks_rotate[rotate_benchmark_count].tfunct = f;
    benchmarks_rotate[rotate_benchmark_count].description = description;
    benchmarks_rotate[rotate_benchmark_count].valid = 0;
    rotate_benchmark_count++;
}

void add_rotate_t_function(lab_test_func f, char *description) 
{
    benchmarks_rotate_t[rotate_t_benchmark_count].tfunct = f;
    benchmarks_rotate_t[rotate_t_benchmark_count].description = description;
    benchmarks_rotate_t[rotate_t_benchmark_count].valid = 0;
    rotate_t_benchmark_count++;
}

void add_blend_function(lab_test_func f, char *description) 
{
    benchmarks_blend[blend_benchmark_count].tfunct = f;
    benchmarks_blend[blend_benchmark_count].description = description;
    benchmarks_blend[blend_benchmark_count].valid = 0;
    blend_benchmark_count++;
}

void add_blend_v_function(lab_test_func f, char *description) 
{
    benchmarks_blend_v[blend_v_benchmark_count].tfunct = f;
    benchmarks_blend_v[blend_v_benchmark_count].description = description;
    benchmarks_blend_v[blend_v_benchmark_count].valid = 0;
    blend_v_benchmark_count++;
}

void add_smooth_function(lab_test_func f, char *description) 
{
    benchmarks_smooth[smooth_benchmark_count].tfunct = f;
    benchmarks_smooth[smooth_benchmark_count].description = description;
    benchmarks_smooth[smooth_benchmark_count].valid = 0;  
    smooth_benchmark_count++;
}

void add_smooth_n_function(lab_test_func f, char *description) 
{
    benchmarks_smooth_n[smooth_n_benchmark_count].tfunct = f;
    benchmarks_smooth_n[smooth_n_benchmark_count].description = description;
    benchmarks_smooth_n[smooth_n_benchmark_count].valid = 0;  
    smooth_n_benchmark_count++;
}

/* 
 * random_in_interval - Returns random integer in interval [low, high) 
 */
static int random_in_interval(int low, int high) 
{
    int size = high - low;
    return (rand()% size) + low;
}

/*
 * create - creates a dimxdim image aligned to a BSIZE byte boundary
 */
static void create(int dim)
{
    int i, j;

    /* Align the images to BSIZE byte boundaries */
    char *orig_b = (char *) data;
    //while ((unsigned)orig_b % BSIZE)
    while ((uintptr_t)orig_b % BSIZE) // W: changed typecast to avoid a compiler warning
	orig_b++;
    orig = (pixel *) orig_b;

    result = orig + dim*dim;
    copy_of_orig = result + dim*dim;

    for (i = 0; i < dim; i++) {
	for (j = 0; j < dim; j++) {
	    /* Original image initialized to random colors */
	    orig[RIDX(i,j,dim)].red = random_in_interval(0, 65536);
	    orig[RIDX(i,j,dim)].green = random_in_interval(0, 65536);
	    orig[RIDX(i,j,dim)].blue = random_in_interval(0, 65536);
	    orig[RIDX(i,j,dim)].alpha = random_in_interval(0, 65536);

	    /* Copy of original image for checking result */
	    copy_of_orig[RIDX(i,j,dim)].red = orig[RIDX(i,j,dim)].red;
	    copy_of_orig[RIDX(i,j,dim)].green = orig[RIDX(i,j,dim)].green;
	    copy_of_orig[RIDX(i,j,dim)].blue = orig[RIDX(i,j,dim)].blue;
	    copy_of_orig[RIDX(i,j,dim)].alpha = orig[RIDX(i,j,dim)].alpha;

	    /* Result image initialized to all black */
	    result[RIDX(i,j,dim)].red = 0;
	    result[RIDX(i,j,dim)].green = 0;
	    result[RIDX(i,j,dim)].blue = 0;
	    result[RIDX(i,j,dim)].alpha = 0; // fully transparent pixel
	}
    }
    /* original background color initialized to random colors */
    bgc.red = random_in_interval(0, 65536);
    bgc.green = random_in_interval(0, 65536);
    bgc.blue = random_in_interval(0, 65536);
    bgc.alpha = USHRT_MAX;

    /* copy of original background color for checking results */
    copy_of_bgc.red = bgc.red;
    copy_of_bgc.green = bgc.green;
    copy_of_bgc.blue = bgc.blue;
    copy_of_bgc.alpha = bgc.alpha;
    
    return;
}


/* 
 * compare_pixels - Returns 1 if the two arguments don't have same RGB
 *    values, 0 o.w.  
 */
static int compare_pixels(pixel p1, pixel p2) 
{
    return 
	(p1.red != p2.red) || 
	(p1.green != p2.green) || 
	(p1.blue != p2.blue) ||
	(p1.alpha != p2.alpha);
}

static int compare_pixels_epsilon(pixel p1, pixel p2, unsigned short epsilon) 
{
    return 
	( abs( (int)p1.red   - (int)p2.red   ) > epsilon ) || 
	( abs( (int)p1.green - (int)p2.green ) > epsilon ) || 
	( abs( (int)p1.blue  - (int)p2.blue  ) > epsilon ) || 
	( abs( (int)p1.alpha - (int)p2.alpha ) > epsilon );
}

/* Make sure the orig array is unchanged */
static int check_orig(int dim) 
{
    int i, j;

    for (i = 0; i < dim; i++) 
	for (j = 0; j < dim; j++) 
	    if (compare_pixels(orig[RIDX(i,j,dim)], copy_of_orig[RIDX(i,j,dim)])) {
		printf("Error: Original image has been changed!\n");
		return 1;
	    }
    
    return 0;
}

static pixel check_blended_pixel(int dim, int i, int j, pixel *src) {
    pixel result;
    float a = ( (float)(src[RIDX(i,j,dim)].alpha) ) / USHRT_MAX;

    result.red   = (a * src[RIDX(i,j,dim)].red  ) + ( (1 - a) * bgc.red  );
    result.green = (a * src[RIDX(i,j,dim)].green) + ( (1 - a) * bgc.green);
    result.blue  = (a * src[RIDX(i,j,dim)].blue) +  ( (1 - a) * bgc.blue );
    result.alpha = USHRT_MAX; // opaque
 
    return result;
}

/* 
 * check_blend - Make sure the blend function actually works.  The
 * orig array should not have been tampered with!  
 */
static int check_blend(int dim) {
    int err = 0;
    int i, j;
    int badi = 0;
    int badj = 0;
    pixel right, wrong;

    /* return 1 if original image has been changed */
    if (check_orig(dim)) 
	return 1; 

    /* return 1 if the background color has been changed */
    if ( compare_pixels ( bgc, copy_of_bgc ) ) {
      printf("Error: Original background color has been changed!\n");
      return 1;
    }	    
    
    for (i = 0; i < dim; i++) {
	for (j = 0; j < dim; j++) {
	    pixel blended = check_blended_pixel(dim, i, j, orig);
	    if (compare_pixels_epsilon(result[RIDX(i,j,dim)], blended, EPSILON)) {
		err++;
		badi = i;
		badj = j;
		wrong = result[RIDX(i,j,dim)];
		right = blended;
	    }
	}
    }

    if (err) {
	printf("\n");
	printf("ERROR: Dimension=%d, %d errors\n", dim, err);    
	printf("E.g., \n");
	printf("You have dst[%d][%d].{red,green,blue,alpha} = {%d,%d,%d,%d}\n",
	       badi, badj, wrong.red, wrong.green, wrong.blue, wrong.alpha);
	printf("It should be dst[%d][%d].{red,green,blue,alpha} = {%d,%d,%d,%d}\n",
	       badi, badj, right.red, right.green, right.blue, right.alpha);
    }

    return err;
}

/* 
 * check_rotate - Make sure the rotate actually works. 
 * The orig array should not have been tampered with! 
 */
static int check_rotate(int dim) 
{
    int err = 0;
    int i, j;
    int badi = 0;
    int badj = 0;
    pixel orig_bad, res_bad;

    /* return 1 if the original image has been  changed */
    if (check_orig(dim)) 
	return 1; 

    for (i = 0; i < dim; i++) 
	for (j = 0; j < dim; j++) 
	    if (compare_pixels(orig[RIDX(i,j,dim)], 
			       result[RIDX(dim-1-j,i,dim)])) {
		err++;
		badi = i;
		badj = j;
		orig_bad = orig[RIDX(i,j,dim)];
		res_bad = result[RIDX(dim-1-j,i,dim)];
	    }

    if (err) {
	printf("\n");
	printf("ERROR: Dimension=%d, %d errors\n", dim, err);    
	printf("E.g., The following two pixels should have equal value:\n");
	printf("src[%d][%d].{red,green,blue,alpha} = {%d,%d,%d,%d}\n",
	       badi, badj, orig_bad.red, orig_bad.green, orig_bad.blue, orig_bad.alpha);
	printf("dst[%d][%d].{red,green,blue,alpha} = {%d,%d,%d,%d}\n",
	       (dim-1-badj), badi, res_bad.red, res_bad.green, res_bad.blue, res_bad.alpha);
    }

    return err;
}

static pixel check_average(int dim, int i, int j, pixel *src) {
    pixel result;
    int num = 0;
    int ii, jj;
    int sum0, sum1, sum2, sum3;
    int top_left_i, top_left_j;
    int bottom_right_i, bottom_right_j;

    top_left_i = max(i-1, 0);
    top_left_j = max(j-1, 0);
    bottom_right_i = min(i+1, dim-1); 
    bottom_right_j = min(j+1, dim-1);

    sum0 = sum1 = sum2 = sum3 = 0;
    for(ii=top_left_i; ii <= bottom_right_i; ii++) {
	for(jj=top_left_j; jj <= bottom_right_j; jj++) {
	    num++;
	    sum0 += (int) src[RIDX(ii,jj,dim)].red;
	    sum1 += (int) src[RIDX(ii,jj,dim)].green;
	    sum2 += (int) src[RIDX(ii,jj,dim)].blue;
	    sum3 += (int) src[RIDX(ii,jj,dim)].alpha;
	}
    }
    result.red = (unsigned short) (sum0/num);
    result.green = (unsigned short) (sum1/num);
    result.blue = (unsigned short) (sum2/num);
    result.alpha = (unsigned short) (sum3/num);
 
    return result;
}


/* 
 * check_smooth - Make sure the smooth function actually works.  The
 * orig array should not have been tampered with!  
 */
static int check_smooth(int dim) {
    int err = 0;
    int i, j;
    int badi = 0;
    int badj = 0;
    pixel right, wrong;

    /* return 1 if original image has been changed */
    if (check_orig(dim)) 
	return 1; 

    for (i = 0; i < dim; i++) {
	for (j = 0; j < dim; j++) {
	    pixel smoothed = check_average(dim, i, j, orig);
	    if (compare_pixels(result[RIDX(i,j,dim)], smoothed)) {
		err++;
		badi = i;
		badj = j;
		wrong = result[RIDX(i,j,dim)];
		right = smoothed;
	    }
	}
    }

    if (err) {
	printf("\n");
	printf("ERROR: Dimension=%d, %d errors\n", dim, err);    
	printf("E.g., \n");
	printf("You have dst[%d][%d].{red,green,blue,alpha} = {%d,%d,%d,%d}\n",
	       badi, badj, wrong.red, wrong.green, wrong.blue, wrong.alpha);
	printf("It should be dst[%d][%d].{red,green,blue,alpha} = {%d,%d,%d,%d}\n",
	       badi, badj, right.red, right.green, right.blue, right.alpha);
    }

    return err;
}


void func_wrapper(void *arglist[]) 
{
    pixel *src, *dst;
    int mydim;
    lab_test_func f;

    f = (lab_test_func) arglist[0];
    mydim = *((int *) arglist[1]);
    src = (pixel *) arglist[2];
    dst = (pixel *) arglist[3];

    (*f)(mydim, src, dst);

    return;
}

void run_rotate_benchmark(int idx, int dim) 
{
    benchmarks_rotate[idx].tfunct(dim, orig, result);
}

void test_rotate(int bench_index) 
{
    int i;
    int test_num;
    char *description = benchmarks_rotate[bench_index].description;
  
    for (test_num = 0; test_num < DIM_CNT; test_num++) {
	int dim;

	/* Check for odd dimension */
	create(ODD_DIM);
	run_rotate_benchmark(bench_index, ODD_DIM);

	if ( ! benchmark_only ) {
	    if (check_rotate(ODD_DIM)) {
		printf("Benchmark \"%s\" failed correctness check for dimension %d.\n\n",
		       benchmarks_rotate[bench_index].description, ODD_DIM);
		return;
	    }
	}

	/* Create a test image of the required dimension */
	dim = test_dim_rotate[test_num];
	create(dim);
#ifdef DEBUG
	printf("DEBUG: Running benchmark \"%s\"\n", benchmarks_rotate[bench_index].description);
#endif

	/* Check that the code works */
	run_rotate_benchmark(bench_index, dim);
	if ( ! benchmark_only ) {
	    if (check_rotate(dim)) {
		printf("Benchmark \"%s\" failed correctness check for dimension %d.\n",
		       benchmarks_rotate[bench_index].description, dim);
		return;
	    }
	}

	/* Measure CPE */
	{
	    double num_cycles, cpe;
	    int tmpdim = dim;
	    void *arglist[4];
	    double dimension = (double) dim;
	    double work = dimension*dimension;
#ifdef DEBUG
	    printf("DEBUG: dimension=%.1f\n",dimension);
	    printf("DEBUG: work=%.1f\n",work);
#endif
	    arglist[0] = (void *) benchmarks_rotate[bench_index].tfunct;
	    arglist[1] = (void *) &tmpdim;
	    arglist[2] = (void *) orig;
	    arglist[3] = (void *) result;

	    create(dim);
	    num_cycles = fcyc_v((test_funct_v)&func_wrapper, arglist); 
	    cpe = num_cycles/work;
	    benchmarks_rotate[bench_index].cpes[test_num] = cpe;
	}
    }

    /* 
     * Print results as a table 
     */
    printf("Rotate, version \"%s\":\n", description);
    printf("Dim\t");
    for (i = 0; i < DIM_CNT; i++)
	printf("\t%d", test_dim_rotate[i]);
    printf("\tMean\n");
  
    printf("Your CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", benchmarks_rotate[bench_index].cpes[i]);
    }
    printf("\n");

    printf("Baseline CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", rotate_baseline_cpes[i]);
    }
    printf("\n");

    /* Compute Speedup */
    {
	double prod, ratio, mean;
	prod = 1.0; /* Geometric mean */
	printf("Speedup\t");
	for (i = 0; i < DIM_CNT; i++) {
	    if (benchmarks_rotate[bench_index].cpes[i] > 0.0) {
		ratio = rotate_baseline_cpes[i]/
		    benchmarks_rotate[bench_index].cpes[i];
	    }
	    else {
		printf("Fatal Error: Non-positive CPE value...\n");
		exit(EXIT_FAILURE);
	    }
	    prod *= ratio;
	    printf("\t%.1f", ratio);
	}

	/* Geometric mean */
	mean = pow(prod, 1.0/(double) DIM_CNT);
	printf("\t%.1f", mean);
	printf("\n\n");
	if (mean > rotate_maxmean) {
	    rotate_maxmean = mean;
	    rotate_maxmean_desc = benchmarks_rotate[bench_index].description;
	}
    }


#ifdef DEBUG
    fflush(stdout);
#endif
    return;  
}

void run_rotate_t_benchmark(int idx, int dim) 
{
    benchmarks_rotate_t[idx].tfunct(dim, orig, result);
}

void test_rotate_t(int bench_index) 
{
    int i;
    int test_num;
    char *description = benchmarks_rotate_t[bench_index].description;
  
    for (test_num = 0; test_num < DIM_CNT; test_num++) {
	int dim;

	/* Check for odd dimension */
	create(ODD_DIM);
	run_rotate_t_benchmark(bench_index, ODD_DIM);

	if ( ! benchmark_only ) {
	    if (check_rotate(ODD_DIM)) {
		printf("Benchmark \"%s\" failed correctness check for dimension %d.\n\n",
		       benchmarks_rotate_t[bench_index].description, ODD_DIM);
		return;
	    }
	}

	/* Create a test image of the required dimension */
	dim = test_dim_rotate[test_num];
	create(dim);
#ifdef DEBUG
	printf("DEBUG: Running benchmark \"%s\"\n", benchmarks_rotate_t[bench_index].description);
#endif

	/* Check that the code works */
	run_rotate_t_benchmark(bench_index, dim);
	if ( ! benchmark_only ) {
	    if (check_rotate(dim)) {
		printf("Benchmark \"%s\" failed correctness check for dimension %d.\n",
		       benchmarks_rotate_t[bench_index].description, dim);
		return;
	    }
	}

	/* Measure CPE */
	{
	    double num_cycles, cpe;
	    int tmpdim = dim;
	    void *arglist[4];
	    double dimension = (double) dim;
	    double work = dimension*dimension;
#ifdef DEBUG
	    printf("DEBUG: dimension=%.1f\n",dimension);
	    printf("DEBUG: work=%.1f\n",work);
#endif
	    arglist[0] = (void *) benchmarks_rotate_t[bench_index].tfunct;
	    arglist[1] = (void *) &tmpdim;
	    arglist[2] = (void *) orig;
	    arglist[3] = (void *) result;

	    create(dim);
	    num_cycles = fcyc_v((test_funct_v)&func_wrapper, arglist); 
	    cpe = num_cycles/work;
	    benchmarks_rotate_t[bench_index].cpes[test_num] = cpe;
	}
    }

    /* 
     * Print results as a table 
     */
    printf("Rotate_T, version \"%s\":\n", description);
    printf("Dim\t");
    for (i = 0; i < DIM_CNT; i++)
	printf("\t%d", test_dim_rotate[i]);
    printf("\tMean\n");
  
    printf("Your CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", benchmarks_rotate_t[bench_index].cpes[i]);
    }
    printf("\n");

    printf("Baseline CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", rotate_baseline_cpes[i]);
    }
    printf("\n");

    /* Compute Speedup */
    {
	double prod, ratio, mean;
	prod = 1.0; /* Geometric mean */
	printf("Speedup\t");
	for (i = 0; i < DIM_CNT; i++) {
	    if (benchmarks_rotate_t[bench_index].cpes[i] > 0.0) {
		ratio = rotate_baseline_cpes[i]/
		    benchmarks_rotate_t[bench_index].cpes[i];
	    }
	    else {
		printf("Fatal Error: Non-positive CPE value...\n");
		exit(EXIT_FAILURE);
	    }
	    prod *= ratio;
	    printf("\t%.1f", ratio);
	}

	/* Geometric mean */
	mean = pow(prod, 1.0/(double) DIM_CNT);
	printf("\t%.1f", mean);
	printf("\n\n");
	if (mean > rotate_t_maxmean) {
	    rotate_t_maxmean = mean;
	    rotate_t_maxmean_desc = benchmarks_rotate_t[bench_index].description;
	}
    }


#ifdef DEBUG
    fflush(stdout);
#endif
    return;  
}

void run_blend_benchmark(int idx, int dim) 
{
    benchmarks_blend[idx].tfunct(dim, orig, result);
}

void test_blend(int bench_index) 
{
    int i;
    int test_num;
    char *description = benchmarks_blend[bench_index].description;
  
    for (test_num = 0; test_num < DIM_CNT; test_num++) {
	int dim;

	/* Check for odd dimension */
	create(ODD_DIM);
	run_blend_benchmark(bench_index, ODD_DIM);

	if ( ! benchmark_only ) {
	    if (check_blend(ODD_DIM)) {
		printf("Benchmark \"%s\" failed correctness check for dimension %d.\n\n",
		       benchmarks_blend[bench_index].description, ODD_DIM);
		return;
	    }
	}

	/* Create a test image of the required dimension */
	dim = test_dim_blend[test_num];
	create(dim);
#ifdef DEBUG
	printf("DEBUG: Running benchmark \"%s\"\n", benchmarks_blend[bench_index].description);
#endif

	/* Check that the code works */
	run_blend_benchmark(bench_index, dim);
	if ( ! benchmark_only ) {
	    if (check_blend(dim)) {
		printf("Benchmark \"%s\" failed correctness check for dimension %d.\n",
		       benchmarks_blend[bench_index].description, dim);
		return;
	    }
	}

	/* Measure CPE */
	{
	    double num_cycles, cpe;
	    int tmpdim = dim;
	    void *arglist[4];
	    double dimension = (double) dim;
	    double work = dimension*dimension;
#ifdef DEBUG
	    printf("DEBUG: dimension=%.1f\n",dimension);
	    printf("DEBUG: work=%.1f\n",work);
#endif
	    arglist[0] = (void *) benchmarks_blend[bench_index].tfunct;
	    arglist[1] = (void *) &tmpdim;
	    arglist[2] = (void *) orig;
	    arglist[3] = (void *) result;

	    create(dim);
	    num_cycles = fcyc_v((test_funct_v)&func_wrapper, arglist); 
	    cpe = num_cycles/work;
	    benchmarks_blend[bench_index].cpes[test_num] = cpe;
	}
    }

    /* 
     * Print results as a table 
     */
    printf("Blend,  version \"%s\":\n", description);
    printf("Dim\t");
    for (i = 0; i < DIM_CNT; i++)
	printf("\t%d", test_dim_blend[i]);
    printf("\tMean\n");
  
    printf("Your CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", benchmarks_blend[bench_index].cpes[i]);
    }
    printf("\n");

    printf("Baseline CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", blend_baseline_cpes[i]);
    }
    printf("\n");

    /* Compute Speedup */
    {
	double prod, ratio, mean;
	prod = 1.0; /* Geometric mean */
	printf("Speedup\t");
	for (i = 0; i < DIM_CNT; i++) {
	    if (benchmarks_blend[bench_index].cpes[i] > 0.0) {
		ratio = blend_baseline_cpes[i]/
		    benchmarks_blend[bench_index].cpes[i];
	    }
	    else {
		printf("Fatal Error: Non-positive CPE value...\n");
		exit(EXIT_FAILURE);
	    }
	    prod *= ratio;
	    printf("\t%.1f", ratio);
	}

	/* Geometric mean */
	mean = pow(prod, 1.0/(double) DIM_CNT);
	printf("\t%.1f", mean);
	printf("\n\n");
	if (mean > blend_maxmean) {
	    blend_maxmean = mean;
	    blend_maxmean_desc = benchmarks_blend[bench_index].description;
	}
    }


#ifdef DEBUG
    fflush(stdout);
#endif
    return;  
}

void run_blend_v_benchmark(int idx, int dim) 
{
    benchmarks_blend_v[idx].tfunct(dim, orig, result);
}

void test_blend_v(int bench_index) 
{
    int i;
    int test_num;
    char *description = benchmarks_blend_v[bench_index].description;
  
    for (test_num = 0; test_num < DIM_CNT; test_num++) {
	int dim;

	/* Check for odd dimension */
	create(ODD_DIM);
	run_blend_v_benchmark(bench_index, ODD_DIM);

	if ( ! benchmark_only ) {
	    if (check_blend(ODD_DIM)) {
		printf("Benchmark \"%s\" failed correctness check for dimension %d.\n\n",
		       benchmarks_blend_v[bench_index].description, ODD_DIM);
		return;
	    }
	}

	/* Create a test image of the required dimension */
	dim = test_dim_blend[test_num];
	create(dim);
#ifdef DEBUG
	printf("DEBUG: Running benchmark \"%s\"\n", benchmarks_blend_v[bench_index].description);
#endif

	/* Check that the code works */
	run_blend_v_benchmark(bench_index, dim);
	if ( ! benchmark_only ) {
	    if (check_blend(dim)) {
		printf("Benchmark \"%s\" failed correctness check for dimension %d.\n",
		       benchmarks_blend_v[bench_index].description, dim);
		return;
	    }
	}

	/* Measure CPE */
	{
	    double num_cycles, cpe;
	    int tmpdim = dim;
	    void *arglist[4];
	    double dimension = (double) dim;
	    double work = dimension*dimension;
#ifdef DEBUG
	    printf("DEBUG: dimension=%.1f\n",dimension);
	    printf("DEBUG: work=%.1f\n",work);
#endif
	    arglist[0] = (void *) benchmarks_blend_v[bench_index].tfunct;
	    arglist[1] = (void *) &tmpdim;
	    arglist[2] = (void *) orig;
	    arglist[3] = (void *) result;

	    create(dim);
	    num_cycles = fcyc_v((test_funct_v)&func_wrapper, arglist); 
	    cpe = num_cycles/work;
	    benchmarks_blend_v[bench_index].cpes[test_num] = cpe;
	}
    }

    /* 
     * Print results as a table 
     */
    printf("Blend_V,  version \"%s\":\n", description);
    printf("Dim\t");
    for (i = 0; i < DIM_CNT; i++)
	printf("\t%d", test_dim_blend[i]);
    printf("\tMean\n");
  
    printf("Your CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", benchmarks_blend_v[bench_index].cpes[i]);
    }
    printf("\n");

    printf("Baseline CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", blend_baseline_cpes[i]);
    }
    printf("\n");

    /* Compute Speedup */
    {
	double prod, ratio, mean;
	prod = 1.0; /* Geometric mean */
	printf("Speedup\t");
	for (i = 0; i < DIM_CNT; i++) {
	    if (benchmarks_blend_v[bench_index].cpes[i] > 0.0) {
		ratio = blend_baseline_cpes[i]/
		    benchmarks_blend_v[bench_index].cpes[i];
	    }
	    else {
		printf("Fatal Error: Non-positive CPE value...\n");
		exit(EXIT_FAILURE);
	    }
	    prod *= ratio;
	    printf("\t%.1f", ratio);
	}

	/* Geometric mean */
	mean = pow(prod, 1.0/(double) DIM_CNT);
	printf("\t%.1f", mean);
	printf("\n\n");
	if (mean > blend_v_maxmean) {
	    blend_v_maxmean = mean;
	    blend_v_maxmean_desc = benchmarks_blend_v[bench_index].description;
	}
    }


#ifdef DEBUG
    fflush(stdout);
#endif
    return;  
}

void run_smooth_benchmark(int idx, int dim) 
{
    benchmarks_smooth[idx].tfunct(dim, orig, result);
}

void test_smooth(int bench_index) 
{
    int i;
    int test_num;
    char *description = benchmarks_smooth[bench_index].description;
  
    for(test_num=0; test_num < DIM_CNT; test_num++) {
	int dim;

	/* Check correctness for odd (non power of two dimensions */
	create(ODD_DIM);
	run_smooth_benchmark(bench_index, ODD_DIM);
	if ( ! benchmark_only ) {
	    if (check_smooth(ODD_DIM)) {
		printf("Benchmark \"%s\" failed correctness check for dimension %d.\n",
		       benchmarks_smooth[bench_index].description, ODD_DIM);
		return;
	    }
	}

	/* Create a test image of the required dimension */
	dim = test_dim_smooth[test_num];
	create(dim);

#ifdef DEBUG
	printf("DEBUG: Running benchmark \"%s\"\n", benchmarks_smooth[bench_index].description);
#endif
	/* Check that the code works */
	run_smooth_benchmark(bench_index, dim);
	if ( ! benchmark_only ) {
	    if (check_smooth(dim)) {
		printf("Benchmark \"%s\" failed correctness check for dimension %d.\n",
		       benchmarks_smooth[bench_index].description, dim);
		return;
	    }
	}

	/* Measure CPE */
	{
	    double num_cycles, cpe;
	    int tmpdim = dim;
	    void *arglist[4];
	    double dimension = (double) dim;
	    double work = dimension*dimension;
#ifdef DEBUG
	    printf("DEBUG: dimension=%.1f\n",dimension);
	    printf("DEBUG: work=%.1f\n",work);
#endif
	    arglist[0] = (void *) benchmarks_smooth[bench_index].tfunct;
	    arglist[1] = (void *) &tmpdim;
	    arglist[2] = (void *) orig;
	    arglist[3] = (void *) result;
        
	    create(dim);
	    num_cycles = fcyc_v((test_funct_v)&func_wrapper, arglist); 
	    cpe = num_cycles/work;
	    benchmarks_smooth[bench_index].cpes[test_num] = cpe;
	}
    }

    /* Print results as a table */
    printf("Smooth: Version = %s:\n", description);
    printf("Dim\t");
    for (i = 0; i < DIM_CNT; i++)
	printf("\t%d", test_dim_smooth[i]);
    printf("\tMean\n");
  
    printf("Your CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", benchmarks_smooth[bench_index].cpes[i]);
    }
    printf("\n");

    printf("Baseline CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", smooth_baseline_cpes[i]);
    }
    printf("\n");

    /* Compute speedup */
    {
	double prod, ratio, mean;
	prod = 1.0; /* Geometric mean */
	printf("Speedup\t");
	for (i = 0; i < DIM_CNT; i++) {
	    if (benchmarks_smooth[bench_index].cpes[i] > 0.0) {
		ratio = smooth_baseline_cpes[i]/
		    benchmarks_smooth[bench_index].cpes[i];
	    }
	    else {
		printf("Fatal Error: Non-positive CPE value...\n");
		exit(EXIT_FAILURE);
	    }
	    prod *= ratio;
	    printf("\t%.1f", ratio);
	}
	/* Geometric mean */
	mean = pow(prod, 1.0/(double) DIM_CNT);
	printf("\t%.1f", mean);
	printf("\n\n");
	if (mean > smooth_maxmean) {
	    smooth_maxmean = mean;
	    smooth_maxmean_desc = benchmarks_smooth[bench_index].description;
	}
    }

    return;  
}

void run_smooth_n_benchmark(int idx, int dim) 
{
    benchmarks_smooth_n[idx].tfunct(dim, orig, result);
}

void test_smooth_n(int bench_index) 
{
    int i;
    int test_num;
    char *description = benchmarks_smooth_n[bench_index].description;
  
    for(test_num=0; test_num < DIM_CNT; test_num++) {
	int dim;

	/* Check correctness for odd (non power of two dimensions */
	create(ODD_DIM);
	run_smooth_n_benchmark(bench_index, ODD_DIM);
	if ( ! benchmark_only ) {
	    if (check_smooth(ODD_DIM)) {
		printf("Benchmark \"%s\" failed correctness check for dimension %d.\n",
		       benchmarks_smooth_n[bench_index].description, ODD_DIM);
		return;
	    }
	}

	/* Create a test image of the required dimension */
	dim = test_dim_smooth[test_num];
	create(dim);

#ifdef DEBUG
	printf("DEBUG: Running benchmark \"%s\"\n", benchmarks_smooth_n[bench_index].description);
#endif
	/* Check that the code works */
	run_smooth_n_benchmark(bench_index, dim);
	if ( ! benchmark_only ) {
	    if (check_smooth(dim)) {
		printf("Benchmark \"%s\" failed correctness check for dimension %d.\n",
		       benchmarks_smooth_n[bench_index].description, dim);
		return;
	    }
	}

	/* Measure CPE */
	{
	    double num_cycles, cpe;
	    int tmpdim = dim;
	    void *arglist[4];
	    double dimension = (double) dim;
	    double work = dimension*dimension;
#ifdef DEBUG
	    printf("DEBUG: dimension=%.1f\n",dimension);
	    printf("DEBUG: work=%.1f\n",work);
#endif
	    arglist[0] = (void *) benchmarks_smooth_n[bench_index].tfunct;
	    arglist[1] = (void *) &tmpdim;
	    arglist[2] = (void *) orig;
	    arglist[3] = (void *) result;
        
	    create(dim);
	    num_cycles = fcyc_v((test_funct_v)&func_wrapper, arglist); 
	    cpe = num_cycles/work;
	    benchmarks_smooth_n[bench_index].cpes[test_num] = cpe;
	}
    }

    /* Print results as a table */
    printf("Smooth: Version = %s:\n", description);
    printf("Dim\t");
    for (i = 0; i < DIM_CNT; i++)
	printf("\t%d", test_dim_smooth[i]);
    printf("\tMean\n");
  
    printf("Your CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", benchmarks_smooth_n[bench_index].cpes[i]);
    }
    printf("\n");

    printf("Baseline CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", smooth_baseline_cpes[i]);
    }
    printf("\n");

    /* Compute speedup */
    {
	double prod, ratio, mean;
	prod = 1.0; /* Geometric mean */
	printf("Speedup\t");
	for (i = 0; i < DIM_CNT; i++) {
	    if (benchmarks_smooth_n[bench_index].cpes[i] > 0.0) {
		ratio = smooth_baseline_cpes[i]/
		    benchmarks_smooth_n[bench_index].cpes[i];
	    }
	    else {
		printf("Fatal Error: Non-positive CPE value...\n");
		exit(EXIT_FAILURE);
	    }
	    prod *= ratio;
	    printf("\t%.1f", ratio);
	}
	/* Geometric mean */
	mean = pow(prod, 1.0/(double) DIM_CNT);
	printf("\t%.1f", mean);
	printf("\n\n");
	if (mean > smooth_n_maxmean) {
	    smooth_n_maxmean = mean;
	    smooth_n_maxmean_desc = benchmarks_smooth_n[bench_index].description;
	}
    }

    return;  
}


void usage(char *progname) 
{
    fprintf(stderr, "Usage: %s [-hqg] [-f <func_file>] [-d <dump_file>]\n", progname);    
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h         Print this message\n");
    fprintf(stderr, "  -l         Larger images\n");
    fprintf(stderr, "  -o <oper>  Benchmark only the specified operation <oper> (rotate, rotate_t, smooth, smooth_a, blend, blend_t)\n");
    fprintf(stderr, "  -b         Benchmark only, skip validity check (use only if you know your functions are correct)\n");
    fprintf(stderr, "  -g         Autograder mode: checks only blend(), rotate() and smooth()\n");
    fprintf(stderr, "  -f <file>  Get test function names from dump file <file>\n");
    fprintf(stderr, "  -d <file>  Emit a dump file <file> for later use with -f\n");
    fprintf(stderr, "  -q         Quit after dumping (use with -d )\n");
    exit(EXIT_FAILURE);
}



int main(int argc, char *argv[])
{
    int i;
    int quit_after_dump = 0;
    int skip_studentname_check = 0;
    int autograder = 0;
    int seed = 1729;
    char c = '0';
    char *bench_func_file = NULL;
    char *func_dump_file = NULL;
    int bench_rotate = 0;
    int bench_rotate_t = 0;
    int bench_blend = 0;
    int bench_blend_v = 0;
    int bench_smooth = 0;
    int bench_smooth_n = 0;
    int large_input = 0;

    /* register all the defined functions */
    register_rotate_functions();
    register_rotate_t_functions();
    register_blend_functions();
    register_blend_v_functions();
    register_smooth_functions();
    register_smooth_n_functions();

    /* parse command line args */
    while ((c = getopt(argc, argv, "lo:tgqf:d:s:h")) != -1)
	switch (c) {

	case 'l':
	    large_input = 1;
	    break;

	/* benchmark only - skip validation checks.
	 * saves maybe 1/6 of the driver's execution time.
	 * only do this if you are certain that your is correct! */
	case 'b':
	    benchmark_only = 1;
	    break;
	    
	/* benchmark only the specified operation.*/
        case 'o':
	    if ( strcmp(optarg, "rotate") == 0 ) {
		bench_rotate = 1;
	    } else
	    if ( strcmp(optarg, "rotate_t") == 0 ) {
		bench_rotate_t = 1;
	    } else
	    if ( strcmp(optarg, "blend") == 0 ) {
		bench_blend = 1;
	    } else
	    if ( strcmp(optarg, "blend_v") == 0 ) {
		bench_blend_v = 1;
	    } else
	    if ( strcmp(optarg, "smooth") == 0 ) {
		bench_smooth = 1;
	    } else
	    if ( strcmp(optarg, "smooth_n") == 0 ) {
		bench_smooth_n = 1;
	    } else {
		usage(argv[0]);
	    }
	    break;
	    
	case 't': /* skip student name check (hidden flag) */
	    skip_studentname_check = 1;
	    break;

	case 's': /* seed for random number generator (hidden flag) */
	    seed = atoi(optarg);
	    break;

	case 'g': /* autograder mode (checks only rotate() and smooth()) */
	    autograder = 1;
	    break;

	case 'q':
	    quit_after_dump = 1;
	    break;

	case 'f': /* get names of benchmark functions from this file */
	    bench_func_file = strdup(optarg);
	    break;

	case 'd': /* dump names of benchmark functions to this file */
	    func_dump_file = strdup(optarg);
	    {
		int i;
		FILE *fp = fopen(func_dump_file, "w");	

		if (fp == NULL) {
		    printf("Can't open file %s\n",func_dump_file);
		    exit(-5);
		}

		if ( bench_rotate ) {
		    for(i = 0; i < rotate_benchmark_count; i++) {
			fprintf(fp, "R:%s\n", benchmarks_rotate[i].description); 
		    }
		}
		if ( bench_rotate_t ) {
		    for(i = 0; i < rotate_t_benchmark_count; i++) {
			fprintf(fp, "T:%s\n", benchmarks_rotate_t[i].description); 
		    }
		}
		if ( bench_blend ) {
		    for(i = 0; i < blend_benchmark_count; i++) {
			fprintf(fp, "B:%s\n", benchmarks_blend[i].description); 
		    }
		}
		if ( bench_blend_v ) {
		    for(i = 0; i < blend_v_benchmark_count; i++) {
			fprintf(fp, "V:%s\n", benchmarks_blend_v[i].description); 
		    }
		}
		if ( bench_smooth ) {
		    for(i = 0; i < smooth_benchmark_count; i++) {
			fprintf(fp, "S:%s\n", benchmarks_smooth[i].description); 
		    }
		}
		if ( bench_smooth_n ) {
		    for(i = 0; i < smooth_n_benchmark_count; i++) {
			fprintf(fp, "N:%s\n", benchmarks_smooth_n[i].description); 
		    }
		}
		fclose(fp);
	    }
	    break;

	case 'h': /* print help message */
	    usage(argv[0]);

	default: /* unrecognized argument */
	    usage(argv[0]);
	}

    if (quit_after_dump) 
	exit(EXIT_SUCCESS);

    /* Print student info */
    if (!skip_studentname_check) {
	if (strcmp("bovik", student.alias) == 0) {
	    printf("%s: Please fill in the student struct in kernels.c.\n", argv[0]);
	    exit(1);
	}
	printf("Name:  %s\n", student.name);
	printf("Alias: %s\n", student.alias);
	printf("Email: %s\n", student.email);
	printf("\n");
    }

    if ( large_input ) {
	printf("Large images; this will take a while (esp. if you benchmark many functions).\n\n");
	test_dim_blend[0]  = 1024;
	test_dim_blend[1]  = 2048;
	test_dim_blend[2]  = 4096;
	test_dim_blend[3]  = 8192;
	test_dim_rotate[0] = 1024;
	test_dim_rotate[1] = 2048;
	test_dim_rotate[2] = 4096;
	test_dim_rotate[3] = 8192;
	test_dim_smooth[0] = 1024;
	test_dim_smooth[1] = 2048;
	test_dim_smooth[2] = 4096;
	test_dim_smooth[3] = 8192;
	blend_baseline_cpes[0]  = B1024;
	blend_baseline_cpes[1]  = B2048;
	blend_baseline_cpes[2]  = B4096;
	blend_baseline_cpes[3]  = B8192;
	rotate_baseline_cpes[0] = R1024;
	rotate_baseline_cpes[1] = R2048;
	rotate_baseline_cpes[2] = R4096;
	rotate_baseline_cpes[3] = R8192;
	smooth_baseline_cpes[0] = S1024;
	smooth_baseline_cpes[1] = S2048;
	smooth_baseline_cpes[2] = S4096;
	smooth_baseline_cpes[3] = S8192;
    }

    srand(seed);

    /* 
     * If we are running in autograder mode, we will only test
     * the rotate(), rotate_t(), blend(), blend_v(), smooth(), and smooth_n() functions.
     */
    if (autograder) {
	if ( ! bench_rotate && ! bench_rotate_t && ! bench_blend && ! bench_blend_v && ! bench_smooth && ! bench_smooth_n ) {
	    bench_rotate = 1;
	    bench_rotate_t = 1;
	    bench_blend = 1;
	    bench_blend_v = 1;
	    bench_smooth = 1;
	    bench_smooth_n = 1;
	}
	
	if ( bench_rotate ) {
	    rotate_benchmark_count = 1;
	    benchmarks_rotate[0].tfunct = rotate;
	    benchmarks_rotate[0].description = "rotate() function";
	    benchmarks_rotate[0].valid = 1;
	}
	
	if ( bench_rotate_t ) {
	    rotate_t_benchmark_count = 1;
	    benchmarks_rotate_t[0].tfunct = rotate_t;
	    benchmarks_rotate_t[0].description = "rotate_t() function";
	    benchmarks_rotate_t[0].valid = 1;
	}	    
	if ( bench_blend ) {
	    blend_benchmark_count = 1;
	    benchmarks_blend[0].tfunct = blend;
	    benchmarks_blend[0].description = "blend() function";
	    benchmarks_blend[0].valid = 1;
	}
	
	if ( bench_blend_v ) {
	    blend_v_benchmark_count = 1;
	    benchmarks_blend_v[0].tfunct = blend_v;
	    benchmarks_blend_v[0].description = "blend_v() function";
	    benchmarks_blend_v[0].valid = 1;
	}
	
	if ( bench_smooth ) {
	    smooth_benchmark_count = 1;
	    benchmarks_smooth[0].tfunct = smooth;
	    benchmarks_smooth[0].description = "smooth() function";
	    benchmarks_smooth[0].valid = 1;
	}
	
	if ( bench_smooth_n ) {
	    smooth_n_benchmark_count = 1;
	    benchmarks_smooth_n[0].tfunct = smooth_n;
	    benchmarks_smooth_n[0].description = "smooth_n() function";
	    benchmarks_smooth_n[0].valid = 1;
	}
    }

    /* 
     * If the user specified a file name using -f, then use
     * the file to determine the versions of rotate and smooth to test
     */
    else if (bench_func_file != NULL) {
	char flag;
	char func_line[256];
	FILE *fp = fopen(bench_func_file, "r");

	if (fp == NULL) {
	    printf("Can't open file %s\n",bench_func_file);
	    exit(-5);
	}
    
	while(func_line == fgets(func_line, 256, fp)) {
	    char *func_name = func_line;
	    char **strptr = &func_name;
	    char *token = strsep(strptr, ":");
	    flag = token[0];
	    func_name = strsep(strptr, "\n");
#ifdef DEBUG
	    printf("Function Description is %s\n",func_name);
#endif

	    if (flag == 'R') {
		for(i=0; i<rotate_benchmark_count; i++) {
		    if (strcmp(benchmarks_rotate[i].description, func_name) == 0) {
			bench_rotate = 1;
			benchmarks_rotate[i].valid = 1;
		    }
		}
	    }
	    else
	    if (flag == 'T') {
		for(i=0; i<rotate_t_benchmark_count; i++) {
		    if (strcmp(benchmarks_rotate_t[i].description, func_name) == 0) {
			bench_rotate_t = 1;
			benchmarks_rotate_t[i].valid = 1;
		    }
		}
	    }
	    else
	    if (flag == 'B') {
		for(i=0; i<blend_benchmark_count; i++) {
		    if (strcmp(benchmarks_blend[i].description, func_name) == 0) {
			bench_blend = 1;
			benchmarks_blend[i].valid = 1;
		    }
		}
	    }
	    else
	    if (flag == 'V') {
		for(i=0; i<blend_v_benchmark_count; i++) {
		    if (strcmp(benchmarks_blend_v[i].description, func_name) == 0) {
			bench_blend_v = 1;
			benchmarks_blend_v[i].valid = 1;
		    }
		}
	    }
	    else
	    if (flag == 'S') {
		for(i=0; i<smooth_benchmark_count; i++) {
		    if (strcmp(benchmarks_smooth[i].description, func_name) == 0) {
			bench_smooth = 1;
			benchmarks_smooth[i].valid = 1;
		    }
		}
	    }      
	    else
	    if (flag == 'N') {
		for(i=0; i<smooth_n_benchmark_count; i++) {
		    if (strcmp(benchmarks_smooth_n[i].description, func_name) == 0) {
			bench_smooth_n = 1;
			benchmarks_smooth_n[i].valid = 1;
		    }
		}
	    }      
	}

	fclose(fp);
    }

    /* 
     * If the user didn't specify a dump file using -f, then 
     * test all of the functions
     */
    else { /* set all valid flags to 1 */
	if ( ! bench_rotate && ! bench_rotate_t && ! bench_blend && ! bench_blend_v && ! bench_smooth && ! bench_smooth_n ) {
	    bench_rotate = 1;
	    bench_rotate_t = 1;
	    bench_blend = 1;
	    bench_blend_v = 1;
	    bench_smooth = 1;
	    bench_smooth_n = 1;
	}

	if ( bench_rotate ) {
	    for (i = 0; i < rotate_benchmark_count; i++)
		benchmarks_rotate[i].valid = 1;
	}
	if ( bench_rotate_t ) {
	    for (i = 0; i < rotate_t_benchmark_count; i++)
		benchmarks_rotate_t[i].valid = 1;
	}
	if ( bench_blend ) {
	    for (i = 0; i < blend_benchmark_count; i++)
		benchmarks_blend[i].valid = 1;
	}
	if ( bench_blend_v ) {
	    for (i = 0; i < blend_v_benchmark_count; i++)
		benchmarks_blend_v[i].valid = 1;
	}
	if ( bench_smooth ) {
	    for (i = 0; i < smooth_benchmark_count; i++)
		benchmarks_smooth[i].valid = 1;
	}
	if ( bench_smooth_n ) {
	    for (i = 0; i < smooth_n_benchmark_count; i++)
		benchmarks_smooth_n[i].valid = 1;
	}
    }

    /* Set measurement (fcyc) parameters */
    set_fcyc_cache_size(1 << 26); /* 70 mB cache size */
    set_fcyc_clear_cache(1); /* clear the cache before each measurement */
    set_fcyc_compensate(1); /* try to compensate for timer overhead */

    if ( bench_rotate ) {
	printf("Benchmarking rotate...\n");
	for (i = 0; i < rotate_benchmark_count; i++) {
	    if (benchmarks_rotate[i].valid)
		test_rotate(i);
	}
    }
    if ( bench_rotate_t ) {
	printf("Benchmarking rotate_t...\n");
	for (i = 0; i < rotate_t_benchmark_count; i++) {
	    if (benchmarks_rotate_t[i].valid)
		test_rotate_t(i);
	}
    }
    if ( bench_blend ) {
	printf("Benchmarking blend...\n");
	for (i = 0; i < blend_benchmark_count; i++) {
	    if (benchmarks_blend[i].valid)
		test_blend(i);
	}
    }
    if ( bench_blend_v ) {
	printf("Benchmarking blend_v...\n");
	for (i = 0; i < blend_v_benchmark_count; i++) {
	    if (benchmarks_blend_v[i].valid)
		test_blend_v(i);
	}
    }
    if ( bench_smooth ) {
	printf("Benchmarking smooth...\n");
	for (i = 0; i < smooth_benchmark_count; i++) {
	    if (benchmarks_smooth[i].valid)
		test_smooth(i);
	}
    }
    if ( bench_smooth_n ) {
	printf("Benchmarking smooth_n...\n");
	for (i = 0; i < smooth_n_benchmark_count; i++) {
	    if (benchmarks_smooth_n[i].valid)
		test_smooth_n(i);
	}
    }

    if (autograder) {
	printf("bestscores:%.1f:%.1f:%.1f:%.1f:%.1f:\n", rotate_maxmean, rotate_t_maxmean, blend_maxmean, blend_v_maxmean, smooth_maxmean);
    }
    else {
	printf("Summary of Your Best Scores:\n");
	if ( bench_rotate ) {
	    printf("  Rotate:   %3.1f (%s)\n", rotate_maxmean, rotate_maxmean_desc);
	}
	if ( bench_rotate_t ) {
	    printf("  Rotate_T: %3.1f (%s)\n", rotate_t_maxmean, rotate_t_maxmean_desc);
	}
	if ( bench_blend ) {
	    printf("  Blend:    %3.1f (%s)\n", blend_maxmean, blend_maxmean_desc);
	}
	if ( bench_blend_v ) {
	    printf("  Blend_V:  %3.1f (%s)\n", blend_v_maxmean, blend_v_maxmean_desc);
	}
	if ( bench_smooth ) {
	    printf("  Smooth:   %3.1f (%s)\n", smooth_maxmean, smooth_maxmean_desc);
	}
	if ( bench_smooth_n ) {
	    printf("  Smooth_N: %3.1f (%s)\n", smooth_n_maxmean, smooth_n_maxmean_desc);
	}
    }

    printf("\n");
    return 0;
}
