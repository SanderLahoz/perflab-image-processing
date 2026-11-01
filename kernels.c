/********************************************************
 * Kernels to be optimized for the OS&C prflab.
 * Acknowledgment: This lab is an extended version of the
 * CS:APP Performance Lab
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
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
    naive_rotate(dim, src, dst);
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
}

/******************************************************************************
 * ROTATE_T KERNEL
 *****************************************************************************/

// Your different versions of the rotate_t kernel go here
// (i.e. rotate with multi-threading)

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
