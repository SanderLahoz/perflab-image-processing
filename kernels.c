/********************************************************
 * Kernels to be optimized for the OS&C prflab.
 * Acknowledgment: This lab is an extended version of the
 * CS:APP Performance Lab
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "defs.h"
#include "smooth.h" // helper functions for naive_smooth
#include "blend.h"  // helper functions for naive_blend

/*
 * Please fill in the following struct
 */
student_t student = {
    "salc",                     /* ITU alias */
    "Sander Lahoz Christensen", /* Full name */
    "salc@itu.dk",              /* Email address */
};

/******************************************************************************
 * ROTATE KERNEL
 *****************************************************************************/

// Your different versions of the rotate kernel go here
char my_rotate_descr[] = "single_loop_rotate: No nested loops";
void my_rotate(int dim, pixel *src, pixel *dst)
{
    int k, i, j;
    int iterations = dim * dim;
    for (k = 0; k < iterations; k++)
    {
        i = k / dim;
        j = k - i * dim;
        dst[RIDX(dim - 1 - j, i, dim)] = src[RIDX(i, j, dim)];
    }
}

char my_rotate_descr_v2[] = "less_multiplication: remove multiplication in inner loop";
void my_rotate_v2(int dim, pixel *src, pixel *dst)
{
    int i, j;
    for (i = 0; i < dim; i++)
    {
        pixel *s = src + i * dim;
        pixel *d = dst + (dim - 1) * dim + i;
        for (j = 0; j < dim; j++)
        {
            *d = *s;
            s++;
            d -= dim;
        }
    }
}

char my_rotate_descr_v3[] = "inter_changed_loops: changing the loop order (dst stride 1)";
void my_rotate_v3(int dim, pixel *src, pixel *dst)
{
    int i, j;
    for (j = 0; j < dim; j++)
    {
        pixel *d = dst + (dim - 1 - j) * dim;
        pixel *s = src + j;
        for (i = 0; i < dim; i++)
        {
            *d = *s;
            d++;
            s += dim;
        }
    }
}

#define B 32
char my_rotate_descr_v4[] = "cache_blocking: BxB cache blocking";
void my_rotate_v4(int dim, pixel *src, pixel *dst)
{
    int ii, jj, i, j;
    for (ii = 0; ii < dim; ii += B)
        for (jj = 0; jj < dim; jj += B)
            for (i = ii; i < ii + B; i++)
                for (j = jj; j < jj + B; j++)
                    dst[RIDX(dim - 1 - j, i, dim)] = src[RIDX(i, j, dim)];
}

char my_rotate_descr_v5[] = "blocking_optimized: cache blocking + pointer walking";
void my_rotate_v5(int dim, pixel *src, pixel *dst)
{
    int ii, jj, i, j;
    for (ii = 0; ii < dim; ii += B)
        for (jj = 0; jj < dim; jj += B)
        {
            pixel *s = src + ii * dim + jj;
            pixel *d = dst + (dim - 1 - jj) * dim + ii;
            for (i = 0; i < B; i++)
            {
                pixel *sp = s;
                pixel *dp = d;
                for (j = 0; j < B; j++)
                {
                    *dp = *sp;
                    sp++;
                    dp -= dim;
                }
                s += dim;
                d++;
            }
        }
}

char my_rotate_descr_v6[] = "unrolling_8_times: blocking + pointers + 8x unrolling";
void my_rotate_v6(int dim, pixel *src, pixel *dst)
{
    int ii, jj, i, j;
    int d1 = dim, d2 = 2 * dim, d3 = 3 * dim, d4 = 4 * dim;
    int d5 = 5 * dim, d6 = 6 * dim, d7 = 7 * dim, d8 = 8 * dim;

    for (ii = 0; ii < dim; ii += B)
        for (jj = 0; jj < dim; jj += B)
        {
            pixel *s = src + ii * dim + jj;
            pixel *d = dst + (dim - 1 - jj) * dim + ii;
            for (i = 0; i < B; i++)
            {
                pixel *sp = s;
                pixel *dp = d;
                for (j = 0; j < B; j += 8)
                {
                    dp[0] = sp[0];
                    dp[-d1] = sp[1];
                    dp[-d2] = sp[2];
                    dp[-d3] = sp[3];
                    dp[-d4] = sp[4];
                    dp[-d5] = sp[5];
                    dp[-d6] = sp[6];
                    dp[-d7] = sp[7];
                    sp += 8;
                    dp -= d8;
                }
                s += dim;
                d++;
            }
        }
}

#define TB 32

/* Rotates the source row-bands [ii_start, ii_end). Each band is TB rows. */
static void rotate_bands(int dim, pixel *src, pixel *dst, int ii_start, int ii_end)
{
    pixel buf[TB * TB]; /* buf[j*TB + i] = src(ii+i, jj+j) */
    int ii, jj, i, j;

    for (ii = ii_start; ii < ii_end; ii += TB)
        for (jj = 0; jj < dim; jj += TB)
        {
            /* Phase 1: sequential reads from src, transpose into buf */
            for (i = 0; i < TB; i++)
            {
                pixel *sp = src + (ii + i) * dim + jj;
                pixel *bp = buf + i;
                for (j = 0; j < TB; j += 8)
                {
                    bp[0 * TB] = sp[0];
                    bp[1 * TB] = sp[1];
                    bp[2 * TB] = sp[2];
                    bp[3 * TB] = sp[3];
                    bp[4 * TB] = sp[4];
                    bp[5 * TB] = sp[5];
                    bp[6 * TB] = sp[6];
                    bp[7 * TB] = sp[7];
                    sp += 8;
                    bp += 8 * TB;
                }
            }
            /* Phase 2: write whole 256-byte row segments to dst */
            pixel *d = dst + (dim - 1 - jj) * dim + ii;
            for (j = 0; j < TB; j++)
            {
                memcpy(d, buf + j * TB, TB * sizeof(pixel));
                d -= dim;
            }
        }
}

char rotate_buf_descr[] = "rotate_buf: tile staged through contiguous buffer";
void rotate_buf(int dim, pixel *src, pixel *dst)
{
    rotate_bands(dim, src, dst, 0, dim);
}

/*
 * naive_rotate - The naive baseline version of rotate
 */
/* stride pattern, visualization (we recommend that you draw this for your functions):
    dst         src
    3 7 B F     0 1 2 3
    2 6 A E     4 5 6 7
    1 5 9 D     8 9 A B
    0 4 8 C     C D E F
 */
char naive_rotate_descr[] = "naive_rotate: Naive baseline implementation";
void naive_rotate(int dim, pixel *src, pixel *dst)
{
    int i, j;

    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            dst[RIDX(dim - 1 - j, i, dim)] = src[RIDX(i, j, dim)];
}

/*
 * rotate - Your current working version of rotate
 * IMPORTANT: This is the version you will be graded on
 */
char rotate_descr[] = "rotate: Current working version";
void rotate(int dim, pixel *src, pixel *dst)
{
    rotate_buf(dim, src, dst);
}

/*
 * register_rotate_functions - Register all of your different versions
 *     of the rotate kernel with the driver by calling the
 *     add_rotate_function() for each test function.
 */
void register_rotate_functions()
{
    add_rotate_function(&rotate, rotate_descr);
    /* ... Register additional test functions here */
    add_rotate_function(&my_rotate, my_rotate_descr);
    add_rotate_function(&my_rotate_v2, my_rotate_descr_v2);
    add_rotate_function(&my_rotate_v3, my_rotate_descr_v3);
    add_rotate_function(&my_rotate_v4, my_rotate_descr_v4);
    add_rotate_function(&my_rotate_v5, my_rotate_descr_v5);
    add_rotate_function(&my_rotate_v6, my_rotate_descr_v6);
    add_rotate_function(&rotate_buf, rotate_buf_descr);
}

/******************************************************************************
 * ROTATE_T KERNEL
 *****************************************************************************/

// Your different versions of the rotate_t kernel go here
// (i.e. rotate with multi-threading)

#define NTHREADS 4 // Should match the amount of CPU cores.
typedef struct
{
    int dim;
    int row_start; /* inclusive, multiple of B */
    int row_end;   /* exclusive, multiple of B */
    pixel *src;
    pixel *dst;
} rot_args_t;

static void *rotate_worker(void *arg)
{
    rot_args_t *a = (rot_args_t *)arg;
    int dim = a->dim;
    pixel *src = a->src, *dst = a->dst;
    int ii, jj, i, j;
    int d1 = dim, d2 = 2 * dim, d3 = 3 * dim, d4 = 4 * dim;
    int d5 = 5 * dim, d6 = 6 * dim, d7 = 7 * dim, d8 = 8 * dim;

    for (ii = a->row_start; ii < a->row_end; ii += B)
        for (jj = 0; jj < dim; jj += B)
        {
            pixel *s = src + ii * dim + jj;
            pixel *d = dst + (dim - 1 - jj) * dim + ii;
            for (i = 0; i < B; i++)
            {
                pixel *sp = s;
                pixel *dp = d;
                for (j = 0; j < B; j += 8)
                {
                    dp[0] = sp[0];
                    dp[-d1] = sp[1];
                    dp[-d2] = sp[2];
                    dp[-d3] = sp[3];
                    dp[-d4] = sp[4];
                    dp[-d5] = sp[5];
                    dp[-d6] = sp[6];
                    dp[-d7] = sp[7];
                    sp += 8;
                    dp -= d8;
                }
                s += dim;
                d++;
            }
        }
    return NULL;
}

char my_rotate_t_descr[] = "rotate_t: pthreads, row-range split, blocked + unrolled";
void my_rotate_t(int dim, pixel *src, pixel *dst)
{
    pthread_t tid[NTHREADS];
    rot_args_t args[NTHREADS];
    int t;
    int blocks = dim / B; /* number of B-row bands */
    int per = blocks / NTHREADS;
    int extra = blocks % NTHREADS;
    int band = 0;

    for (t = 0; t < NTHREADS; t++)
    {
        int n = per + (t < extra ? 1 : 0);
        args[t].dim = dim;
        args[t].src = src;
        args[t].dst = dst;
        args[t].row_start = band * B;
        args[t].row_end = (band + n) * B;
        band += n;
        pthread_create(&tid[t], NULL, rotate_worker, &args[t]);
    }
    for (t = 0; t < NTHREADS; t++)
        pthread_join(tid[t], NULL);
}

/*
 * rotate_t - Your current working version of rotate_t
 * IMPORTANT: This is the version you will be graded on
 */
char rotate_t_descr[] = "rotate_t: Current working version";
void rotate_t(int dim, pixel *src, pixel *dst)
{
    naive_rotate(dim, src, dst);
}

/*********************************************************************
 * register_rotate_t_functions - Register all of your different versions
 *     of the rotate_t kernel with the driver by calling the
 *     add_rotate_t_function() for each test function. When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.
 *********************************************************************/

void register_rotate_t_functions()
{
    add_rotate_t_function(&rotate_t, rotate_t_descr);
    /* ... Register additional test functions here */
    add_rotate_t_function(&my_rotate_t, my_rotate_t_descr);
    add_rotate_t_function(&rotate_buf, rotate_buf_descr);
}

/******************************************************************************
 * SMOOTH KERNEL
 *****************************************************************************/

// Your different versions of the smooth kernel go here

/*
 * naive_smooth - The naive baseline version of smooth
 */
char naive_smooth_descr[] = "naive_smooth: Naive baseline implementation";
void naive_smooth(int dim, pixel *src, pixel *dst)
{
    int i, j;

    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            dst[RIDX(i, j, dim)] = avg(dim, i, j, src); // `avg` defined in smooth.c
}

char smooth_descr[] = "smooth: Current working version";
void smooth(int dim, pixel *src, pixel *dst)
{
    naive_smooth(dim, src, dst);
}

/*
 * register_smooth_functions - Register all of your different versions
 *     of the smooth kernel with the driver by calling the
 *     add_smooth_function() for each test function.
 */

void register_smooth_functions()
{
    add_smooth_function(&smooth, smooth_descr);
    /* ... Register additional test functions here */
}

/******************************************************************************
 * SMOOTH_N KERNEL
 *****************************************************************************/

// Your different versions of the smooth_n kernel go here
// (i.e. where anything goes, including multithreading and SIMD).

char smooth_n_descr[] = "smooth_n: Current working version";
void smooth_n(int dim, pixel *src, pixel *dst)
{
    naive_smooth(dim, src, dst);
}

/*
 * register_smooth_n_functions - Register all of your different versions
 *     of the smooth_n kernel with the driver by calling the
 *     add_smooth_n_function() for each test function.
 */
void register_smooth_n_functions()
{
    add_smooth_n_function(&smooth_n, smooth_n_descr);
    /* ... Register additional test functions here */
}

/******************************************************************************
 * BLEND KERNEL
 *****************************************************************************/

// Your different versions of the blend kernel go here.

char naive_blend_descr[] = "naive_blend: Naive baseline implementation";
void naive_blend(int dim, pixel *src, pixel *dst) // reads global variable `pixel bgc`
{
    int i, j;

    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++)
            blend_pixel(&src[RIDX(i, j, dim)], &dst[RIDX(i, j, dim)], &bgc); // `blend_pixel` defined in blend.c
}

char blend_descr[] = "blend: Current working version";
void blend(int dim, pixel *src, pixel *dst)
{
    naive_blend(dim, src, dst);
}

/*
 * register_blend_functions - Register all of your different versions
 *     of the blend kernel with the driver by calling the
 *     add_blend_function() for each test function.
 */
void register_blend_functions()
{
    add_blend_function(&blend, blend_descr);
    /* ... Register additional test functions here */
}

/******************************************************************************
 * BLEND_V KERNEL
 *****************************************************************************/

// Your different versions of the blend_v kernel go here
// (i.e. with vectorization, aka. SIMD).

char blend_v_descr[] = "blend_v: Current working version";
void blend_v(int dim, pixel *src, pixel *dst)
{
    naive_blend(dim, src, dst);
}

/*
 * register_blend_v_functions - Register all of your different versions
 *     of the blend_v kernel with the driver by calling the
 *     add_blend_function() for each test function.
 */
void register_blend_v_functions()
{
    add_blend_v_function(&blend_v, blend_v_descr);
    /* ... Register additional test functions here */
}
