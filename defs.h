/*
 * driver.h - Various definitions for the Performance Lab.
 * 
 * DO NOT MODIFY ANYTHING IN THIS FILE
 */
#ifndef _DEFS_H_
#define _DEFS_H_

#include <stdlib.h>

#define RIDX(i,j,n) ((i)*(n)+(j))

typedef struct {
    char *alias, *name, *email;
} student_t;

extern student_t student;

typedef struct {
    unsigned short red;
    unsigned short green;
    unsigned short blue;
    unsigned short alpha;
} pixel;

extern pixel bgc;

typedef void (*lab_test_func) (int, pixel*, pixel*);

void rotate(int, pixel *, pixel *);
void rotate_t(int, pixel *, pixel *);
void blend(int, pixel *, pixel *);
void blend_v(int, pixel *, pixel *);
void smooth(int, pixel *, pixel *);
void smooth_n(int, pixel *, pixel *);

void register_rotate_functions(void);
void register_rotate_t_functions(void);
void register_blend_functions(void);
void register_blend_v_functions(void);
void register_smooth_functions(void);
void register_smooth_n_functions(void);
void add_rotate_function(lab_test_func, char*);
void add_rotate_t_function(lab_test_func, char*);
void add_blend_function(lab_test_func, char*);
void add_blend_v_function(lab_test_func, char*);
void add_smooth_function(lab_test_func, char*);
void add_smooth_n_function(lab_test_func, char*);

#endif /* _DEFS_H_ */
