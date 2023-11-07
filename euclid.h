#ifndef EUCLID_H
#define EUCLID_H
#include "basics.h"

typedef struct
{
	void *inc[2];
	char flag;
	int addr;
} point;

typedef struct
{
	point *center;
	point *radius;
	int addr;
} circle;

typedef struct
{
	point *a;
	point *b;
	int addr;
} line;

typedef struct
{
	array_voidstar *points;
	array_voidstar *circles;
	array_voidstar *lines;
	array_char history;
} sc_constr;



void point_init_rooted(point *p, double *x, double *y);
void point_init_ll(point *p, line *a, line *b);
void point_init_lc(point *p, line *a, circle *b, char lr);
void point_init_cl(point *p, circle *b, line *a, char lr);
void point_init_cc(point *p, circle *a, circle *b, char lr);
void point_init(point *p, void *a, void *b, char case_);
void point_coords(point *p, double *x, double *y);

void circle_init(circle *c, point *center, point *radius);
void circle_coords(circle *c, double *cx, double *cy, double *r);


void line_init(line *l, point *a, point *b);
void line_coords(line *l, double *ax, double *ay, double *bx, double *by);

void sc_constr_init(sc_constr *c);
void add2sc_constr(sc_constr *c, void *a, char mode);
void sc_constr_undo(sc_constr *c);
char sc_constr_check_inter_ll(sc_constr *c, int l1, int l2);
char sc_constr_check_inter_lc(sc_constr *c, int l_, int c_);
char sc_constr_check_inter_cc(sc_constr *c, int c1, int c2);
void add_point_sc_constr_ll(sc_constr *c, int l1, int l2);
void add_point_sc_constr_lc(sc_constr *c, int l_, int c_, char lr);
void add_point_sc_constr_cc(sc_constr *c, int c1, int c2, char lr);
void add_line_sc_constr_pp(sc_constr *c, int a_, int b_);
void add_circle_sc_constr_pp(sc_constr *c, int c_, int r_);
void free_sc_constr(sc_constr *c);

void sc_constr_point_coords_rem(sc_constr *c, array_double *xs, array_double *ys, array_char *covered, int i);
void sc_constr_point_coords_rem_exp(sc_constr *c, array_double *xs, array_double *ys, array_char *covered, int i, double *x_, double *y_);
void sc_constr_points_coords(sc_constr *c, array_double *xs, array_double *ys);
char sc_constr_check_singular_lc(sc_constr *sc, int l, int c);
void sc_constr_print_point_coords(sc_constr *sc);
void sc_constr_print_line_coords(sc_constr *sc);
void sc_constr_print_circle_coords(sc_constr *sc);

char line_line_intersection_exp(double a1x, double a1y, double b1x, double b1y, double a2x, double a2y, double b2x, double b2y, double *x, double *y);
char line_circle_intersection_exp(double ax, double ay, double bx, double by, double cx, double cy, double rx, double ry, double *x1, double *y1, double *x2, double *y2);
char circle_circle_intersection_exp(double c1x, double c1y, double r1x, double r1y, double c2x, double c2y, double r2x, double r2y, double *x1, double *y1, double *x2, double *y2);

#endif
