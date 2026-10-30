
#ifndef _SMOOTH_H_
#define _SMOOTH_H_

/* A struct used to compute averaged pixel value */
typedef struct {
    int red;
    int green;
    int blue;
    int alpha;
    int num;
} pixel_sum;

/* Compute min and max of two integers, respectively */
int min(int a, int b);
int max(int a, int b);

/* initialize_pixel_sum - Initializes all fields of sum to 0 */
void initialize_pixel_sum(pixel_sum *sum);

/* accumulate_sum - Accumulates field values of p in corresponding 
 * fields of sum */
void accumulate_sum(pixel_sum *sum, pixel p);

/* assign_sum_to_pixel - Computes averaged pixel value in current_pixel */
void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum);

/* avg - Returns averaged pixel value at (i,j) */
pixel avg(int dim, int i, int j, pixel *src);

#endif
