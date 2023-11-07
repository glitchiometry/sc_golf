// RESUME: avoid duplicate points/lines/circles

#include "euclid.h"

#define LR_MASK 8
#define COMPLEX_RCODE 1
#define PARALLEL_RCODE 2
#define EXCEP_LC_RCODE 3
#define COMPLEX_FLAG -31
#define INF_FLAG -15
const double HALF_SQRT3 = 0.5 * sqrt(3.0);

void point_init_rooted(point *p, double *x, double *y)
{
	(*p).flag = 1;
	(*p).inc[0] = (void *) x;
	(*p).inc[1] = (void *) y;
}

void point_init_ll(point *p, line *a, line *b)
{
	point_init(p, (void *) a, (void *) b, 0);
}

void point_init_lc(point *p, line *a, circle *b, char lr)
{
	point_init(p, (void *) a, (void *) b, 2 | (lr << 3));
	//(*p).flag = 2 | (lr << 3);
	//(*p).inc[0] = (void *) a;
	//(*p).inc[1] = (void *) b;
}

void point_init_cl(point *p, circle *b, line *a, char lr)
{
	point_init_lc(p, a, b, lr);
	//point_init(p, (void *) b, (void *) a, 4 | (lr << 3));
	//(*p).flag = 4 | (lr << 3);
	//(*p).inc[0] = (void *) a;
	//(*p).inc[1] = (void *) b;
}

void point_init_cc(point *p, circle *a, circle *b, char lr)
{
	point_init(p, (void *) a, (void *) b, 6 | (lr << 3));
}

void point_init(point *p, void *a, void *b, char case_)
{
	char case_disp[8];
	int opp = 7;
	for (int i = 0; i < 8; i++) 
	{
		case_disp[opp] = (case_ >> i) & 1 ? '1' : '0';
		opp -= 1;
	}
	//printf("Initializing point with case %s = %d\n", case_disp, case_);
	(*p).flag = case_;
	(*p).inc[0] = a;
	(*p).inc[1] = b;
}

char line_line_intersection_exp(double l1a_x, double l1a_y, double l1b_x, double l1b_y, 
				double l2a_x, double l2a_y, double l2b_x, double l2b_y, double *x_, double *y_)
{
	l1b_x -= l1a_x;
	l1b_y -= l1a_y;
	l2b_x -= l2a_x;
	l2b_y -= l2a_y;
	l2a_x -= l1a_x;
	l2a_y -= l1a_y;
	double denom = l1b_y * l2b_x - l1b_x * l2b_y;
	if (denom == 0)
	{
		printf("Parallel lines (intersect at infinity)\n");
		printf("(Feature not yet supported)\n");
		(*x_) = 1e99;
		(*y_) = 1e99;
		return PARALLEL_RCODE;
	}
	double ta = (l2b_x * l2a_y - l2b_y * l2a_x) / denom;
	(*x_) = l1a_x + l1b_x * ta;
	(*y_) = l1a_y + l1b_y * ta;
	//;;;printf("Line-line intersection: %g %g\n", (*x_), (*y_));
	return 0;
}

char line_circle_intersection_exp(double lx1, double ly1, double lx2, double ly2, double ccx, double ccy, double crx, double cry, double *x1, double *y1, double *x2, double *y2)
{
	//;;;printf("\tLine-circle intersection (exp): (%g %g: %g %g), <%g, %g: %g, %g>\n", lx1, ly1, lx2, ly2, ccx, ccy, crx, cry);
	double drx = crx - ccx, dry = cry - ccy, vx = lx2 - lx1, vy = ly2 - ly1;
	double rsq = drx * drx + dry * dry, t, vsq = vx * vx + vy * vy, discr;

	double ly1mcy = ly1 - ccy, lx1mcx = lx1 - ccx;
	double vdotl1mc = vx * lx1mcx + vy * ly1mcy, l1mcsq = lx1mcx * lx1mcx + ly1mcy * ly1mcy;
	discr = vdotl1mc * vdotl1mc - vsq * (l1mcsq - rsq);
	if (discr >= 0)
	{
		discr = sqrt(4 * discr);
		double denom = 1. / (2 * vsq);
		//double m2b = -4 * vdotl1mc;
		double mb = -2 * vdotl1mc;
		double t1 = (mb + discr) * denom;
		double t2 = (mb - discr) * denom;
		(*x1) = lx1 + t1 * vx;
		(*y1) = ly1 + t1 * vy;
		(*x2) = lx1 + t2 * vx;
		(*y2) = ly1 + t2 * vy;
		//;;;printf("Line-circle intersection (generic): (%g, %g), (%g, %g)\n", (*x1), (*y1), (*x2), (*y2));
		return 0;
	}
	else
	{
		//;;;printf("Line and circle failed to intersect: line:(%g %g: %g %g), circle (%g %g: %g %g)\n", lx1, ly1, lx2, ly2, ccx, ccy, crx, cry);
		return COMPLEX_RCODE;
	}
}

char line_circle_intersection(line *l, circle *c, double *x1, double *y1, double *x2, double *y2)
{
	double lx1, ly1, lx2, ly2, ccx, ccy, crx, cry;
	line_coords(l, &lx1, &ly1, &lx2, &ly2);
	point_coords((*c).center, &ccx, &ccy);
	point_coords((*c).radius, &crx, &cry);
	if (!((*l).a == (*c).center) && !((*l).b == (*c).center)) {}
	else
	{
		double drx = crx - ccx, dry = cry - ccy, vx = lx2 - lx1, vy = ly2 - ly1;
		double rsq = drx * drx + dry * dry, vsq = vx * vx + vy * vy;

		double ratio = sqrt(rsq / vsq);
		//;;;printf("Exceptional line-circle intersection: line (%g, %g; %g, %g), circle (%g, %g; %g, %g, r = %g)\n", lx1, ly1, lx2, ly2, ccx, ccy, crx, cry, sqrt(rsq));
		// r * (b - a) / || b - a ||
		(*x1) = ratio * vx;
		(*y1) = ratio * vy;
		(*x2) = -(*x1) + ccx;
		(*y2) = -(*y1) + ccy;
		(*x1) += ccx;
		(*y1) += ccy;
		//;;;printf("\tIntersection: (%g %g), (%g %g)\n", (*x1), (*y1), (*x2), (*y2));
		return EXCEP_LC_RCODE;
	}
	return line_circle_intersection_exp(lx1, ly1, lx2, ly2, ccx, ccy, crx, cry, x1, y1, x2, y2);
}

char circle_circle_intersection_exp(double c1x, double c1y, double r1x, double r1y, double c2x, double c2y, double r2x, double r2y, double *x1, double *y1, double *x2, double *y2)
{
	//;;;printf("Circle circle intersection: %g %g %g %g, %g %g %g %g\n", c1x, c1y, r1x, r1y, c2x, c2y, r2x, r2y);
	double c1xmc2x, c1ymc2y, c1mc2sq, r2sqmr1sq;
	c1xmc2x = c1x - c2x;
	c1ymc2y = c1y - c2y;
	double r1sq, r2sq;
	double alpha, c1xmc2xsq, c1ymc2ysq, discr;
	r1x -= c1x;
	r1y -= c1y;
	r2x -= c2x;
	r2y -= c2y;
	r1sq = r1x * r1x + r1y * r1y;
	r2sq = r2x * r2x + r2y * r2y;
	r2sqmr1sq = r2sq - r1sq;
	c1xmc2xsq = c1xmc2x * c1xmc2x;
	c1ymc2ysq = c1ymc2y * c1ymc2y;
	c1mc2sq = c1xmc2xsq + c1ymc2ysq;
	alpha = r2sqmr1sq + c1mc2sq;
	if (c1xmc2xsq >= c1ymc2ysq)
	{
		double C = 0.25 * alpha * alpha - r2sq * c1xmc2xsq;
		double B = -c1ymc2y * alpha;
		double A = c1mc2sq;
		discr = B * B - 4 * A * C;
		if (discr >= 0)
		{
			discr = sqrt(discr);
			double denom = 1. / (2 * A);
			(*y1) = (-B + discr) * denom;
			(*y2) = (-B - discr) * denom;
			double inv2c1xmc2x = 1. / (2 * c1xmc2x), tc1ymc2y = 2 * c1ymc2y;
			(*x1) = (alpha - tc1ymc2y * (*y1)) * inv2c1xmc2x;
			(*x2) = (alpha - tc1ymc2y * (*y2)) * inv2c1xmc2x;
		}
		else
		{
			printf("No overlap detected\n");
			return COMPLEX_RCODE;
		}
	}
	else
	{
		double C = 0.25 * alpha * alpha - r2sq * c1ymc2ysq;
		double B = -c1xmc2x * alpha;
		double A = c1mc2sq;
		discr = B * B - 4 * A * C;
		if (discr >= 0)
		{
			discr = sqrt(discr);
			double denom = 1. / (2 * A);
			(*x1) = (-B + discr) * denom;
			(*x2) = (-B - discr) * denom;
			double inv2c1ymc2y = 1. / (2 * c1ymc2y), tc1xmc2x = 2 * c1xmc2x;
			(*y1) = (alpha - tc1xmc2x * (*x1)) * inv2c1ymc2y;
			(*y2) = (alpha - tc1xmc2x * (*x2)) * inv2c1ymc2y;
		}
		else
		{
			return COMPLEX_RCODE;
		}
	}
	// Swap intersection points to match convention for lr_flag calculation:
	// Test: x1 x (c1 - c2) 
	double cp_test = (*x1) * c1ymc2y - (*y1) * c1xmc2x;
	if (cp_test > 0) {}
	else
	{
		double aux = (*x1);
		(*x1) = (*x2);
		(*x2) = aux;
		aux = (*y1);
		(*y1) = (*y2);
		(*y2) = aux;
	}
	(*x1) += c2x;
	(*y1) += c2y;
	(*x2) += c2x;
	(*y2) += c2y;
	//;;;printf("Circle-circle intersection (generic): (%g, %g), (%g, %g)\n", (*x1), (*y1), (*x2), (*y2));
	return 0;
}

char circle_circle_intersection(circle *c1, circle *c2, double *x1, double *y1, double *x2, double *y2)
{
		
	double c1x, c1y, c2x, c2y; 
	double r1x, r1y, r2x, r2y;
	point_coords((*c1).center, &c1x, &c1y);
	point_coords((*c2).center, &c2x, &c2y);
	point_coords((*c1).radius, &r1x, &r1y);
	point_coords((*c2).radius, &r2x, &r2y);
	if ((*c1).center != (*c2).radius || (*c1).radius != (*c2).center) {}
	else
	{
		r1x -= c1x;
		r1y -= c1y;
		r2x -= c2x;
		r2y -= c2y;
		//;;;printf("Exceptional circle-circle intersection: circle 1: (%g, %g; %g), circle 2: (%g, %g; %g)", c1x, c1y, sqrt(r1x * r1x + r1y * r1y), c2x, c2y, sqrt(r2x * r2x + r2y * r2y));
		double c1xmc2x, c1ymc2y, c1mc2sq, r2sqmr1sq;
		c1xmc2x = c1x - c2x;
		c1ymc2y = c1y - c2y;
		double mdptx = 0.5 * (c1x + c2x);
		double mdpty = 0.5 * (c1y + c2y);
		double delx = -HALF_SQRT3 * c1ymc2y;
		double dely = HALF_SQRT3 * c1xmc2x; // delx * c1ymc2y - dely * c1xmc2x
		(*x1) = mdptx - delx;
		(*y1) = mdpty - dely;
		(*x2) = mdptx + delx;
		(*y2) = mdpty + dely;
		//;;;printf("\tIntersection: (%g, %g), (%g, %g)\n", (*x1), (*y1), (*x2), (*y2));
		return 0;
	}
	return circle_circle_intersection_exp(c1x, c1y, r1x, r1y, c2x, c2y, r2x, r2y, 
			x1, y1, x2, y2);
}

void point_coords_ll(point *p, double *x, double *y)
{
	// RESUME
	line * l1 = (line *) (*p).inc[0];
	line * l2 = (line *) (*p).inc[1];
	double a1x, a1y, b1x, b1y, a2x, a2y, b2x, b2y;
	line_coords(l1, &a1x, &a1y, &b1x, &b1y);
	line_coords(l2, &a2x, &a2y, &b2x, &b2y);
	char status = line_line_intersection_exp(a1x, a1y, b1x, b1y, a2x, a2y, b2x, b2y, x, y);
	/*
	b1x -= a1x;
	b1y -= a1y;
	b2x -= a2x;
	b2y -= a2y;
	// (a1x, a1y) + (b1x, b1y) t1 = (a2x, a2y) + (b2x, b2y) t2
	// b1 t1 - b2 t2 = a2 - a1
	a2x -= a1x;
	a2y -= a1y;
	// b1x t1 - b2x t2 = a2x
	// 	t1 = (a2x + b2x t2) / b1x
	// b1y (a2x + b2x t2) - b2y b1x t2 = a2y b1x
	//	t2 = (a2y b1x - a2x b1y) / (b1y b2x - b2y b1x)
	// 	t1 = (a1y b2x - a1x b2y) / (b2y b1x - b1y b2x)
	double t = (a1y * b2x - a1x * b2y) / (b2y * b1x - b1y * b2x);
	(*x) = a1x + b1x * t;
	(*y) = a1y + b1y * t;
	*/
}

void point_coords(point *p, double *x, double *y)
{
	if ((*p).flag == COMPLEX_FLAG || (*p).flag == INF_FLAG)
	{
		printf("Setting point to infinite coords\n");
		(*x) = 1e99;
		(*y) = 1e99;
		return;
	}
	if ((*p).flag == 1)
	{
		(*x) = *((double *) (*p).inc[0]);
		(*y) = *((double *) (*p).inc[1]);
		return;
	}
	char mode_test = ((*p).flag >> 1) & 3;
	if (mode_test) 
	{
		double x1, y1, x2, y2;
		if (mode_test == 1 || mode_test == 2)
		{
			line *l_;
			circle *c_;
			switch (mode_test)
			{
				case 1:
					l_ = (line *) (*p).inc[1];
					c_ = (circle *) (*p).inc[0];
					break;
				case 2:
					l_ = (line *) (*p).inc[0];
					c_ = (circle *) (*p).inc[1];
					break;
			}
			char status = line_circle_intersection(l_, c_, &x1, &y1, &x2, &y2);
			if (status == 0)
			{
				double ax, ay, bx, by;
				line_coords(l_, &ax, &ay, &bx, &by);
				bx -= ax;
				by -= ay;
				double x2mx1 = x2 - x1;
				double y2my1 = y2 - y1;
				double test_quant = x2mx1 * bx + y2my1 * by;
				char fb_state = (*p).flag & LR_MASK;
				char aligned_ = test_quant >= 0;
				if (fb_state && aligned_ || !(fb_state || aligned_))
				{
					(*x) = x1;
					(*y) = y1;
				}
				else if (!fb_state && aligned_ || fb_state && !aligned_)
				{
					(*x) = x2;
					(*y) = y2;
				}
				return;
			}
			else if (status == EXCEP_LC_RCODE)
			{
				if ((*p).flag & LR_MASK)
				{
					(*x) = x2;
					(*y) = y2;
				}
				else
				{
					(*x) = x1;
					(*y) = y1;
				}
				return;
			}
			else if (status == COMPLEX_RCODE)
			{
				(*p).flag = COMPLEX_FLAG;
				return;
			}
		}
		else if (mode_test == 3)
		{
			char mode_ = (*p).flag & LR_MASK;
			double circ_x, circ_y, rad_x, rad_y;
			double cir_x, cir_y, ra_x, ra_y, rad_, rad__;
			char status;
			circle *c_;
			c_ = (circle *) (*p).inc[0];
			//;;circle *c__;
			//;;c__ = (circle *) (*p).inc[1];
			status = circle_circle_intersection(c_, (circle *) (*p).inc[1], 
					&x1, &y1, &x2, &y2);
			point_coords((*c_).center, &circ_x, &circ_y);
			//;;point_coords((*c__).center, &cir_x, &cir_y);
			//;;circle_coords(c_, &circ_x, &circ_y, &rad_);
			//;;circle_coords(c__, &cir_x, &cir_y, &rad__);
			if (status == 0) {}
			else
			{
				if (status == COMPLEX_RCODE) 
				{
					(*p).flag = COMPLEX_FLAG;
					return;
				}
			}
			//double ave_x, ave_y, loc_x1 = x1 - circ_x, loc_y1 = y1 - circ_y;
			//ave_x = x1 + x2 - 2 * circ_x;
			//ave_y = y1 + y2 - 2 * circ_y;
			//double cross_1;
			//cross_1 = ave_x * loc_y1 - ave_y * loc_x1;
			//;;;printf("\t lr flag: %d\n", mode_ & 1);

			//;;printf("CC point %d, flag %d: (%g, %g; %g) -> (%g, %g; %g), ave = %g %g, loc = %g %g\n", (*p).addr, (*p).flag, circ_x, circ_y, rad_, cir_x, cir_y, rad__, ave_x, ave_y, loc_x1, loc_y1);
			/*if ((cross_1 <= 0) && ((*p).flag & LR_MASK) || !(cross_1 <= 0 || (*p).flag & LR_MASK))
			{
				(*x) = x1;
				(*y) = y1;
			}
			else
			{
				(*x) = x2;
				(*y) = y2;
			}*/
			char lr_flag = ((*p).flag >> 3) & 1; // RESUME: LR_
			if (lr_flag)
			{
				(*x) = x2;
				(*y) = y2;
			}
			else
			{
				(*x) = x1;
				(*y) = y1;
			}
			return;
		}
	}
	else
	{
		// Line-line intersection
		line *a = (line *) (*p).inc[0];
		line *b = (line *) (*p).inc[1];
		double a1x, a2x, a1y, a2y, b1x, b2x, b1y, b2y;
		line_coords(a, &a1x, &a1y, &a2x, &a2y);
		line_coords(b, &b1x, &b1y, &b2x, &b2y);
		char status = line_line_intersection_exp(a1x, a1y, a2x, a2y, b1x, b1y, b2x, b2y, x, y);
		if (status == 0) {}
		else
		{
			(*p).flag |= INF_FLAG;
		}
	}
}

void circle_init(circle *c, point *center, point *radius)
{
	(*c).center = center;
	(*c).radius = radius;
}

void circle_coords(circle *c, double *cx, double *cy, double *r)
{
	double rx, ry;
	point_coords((*c).center, cx, cy);
	point_coords((*c).radius, &rx, &ry);
	rx -= (*cx);
	ry -= (*cy);
	(*r) = sqrt(rx * rx + ry * ry);
}


void line_init(line *l, point *a, point *b)
{
	(*l).a = a;
	(*l).b = b;
}

void line_coords(line *l, double *ax, double *ay, double *bx, double *by)
{
	point_coords((*l).a, ax, ay);
	point_coords((*l).b, bx, by);
}


void sc_constr_init(sc_constr *c)
{
	(*c).points = (array_voidstar *) calloc(1, sizeof(array_voidstar));
	(*c).lines = (array_voidstar *) calloc(1, sizeof(array_voidstar));
	(*c).circles = (array_voidstar *) calloc(1, sizeof(array_voidstar));
	array_voidstar_init((*c).points, 1);
	array_voidstar_init((*c).circles, 1);
	array_voidstar_init((*c).lines, 1);
	array_char_init(&((*c).history), 1);
}

void add2sc_constr(sc_constr *c, void *a, char mode)
{
	add2array_char(&((*c).history), mode);
	switch (mode)
	{
		case 'p':
			(*((point *) a)).addr = (*(*c).points).len;
			add2array_voidstar((*c).points, a);
			break;
		case 'l':
			(*((line *) a)).addr = (*(*c).lines).len;
			add2array_voidstar((*c).lines, a);
			break;
		case 'c':
			(*((circle *) a)).addr = (*(*c).circles).len;
			add2array_voidstar((*c).circles, a);
			break;
	}
}

void sc_constr_undo(sc_constr *c)
{
	int new_len = (*c).history.len - 1;
	char last_mode = (*c).history.e[new_len];
	remove_array_char(&((*c).history), new_len);
	switch (last_mode)
	{
		int last_addr;
		case 'p':
			last_addr = (*(*c).points).len - 1;
			free((point *) (*(*c).points).e[last_addr]);
			(*(*c).points).e[last_addr] = NULL;
			remove_array_voidstar((*c).points, last_addr, NULL);
			break;
		case 'l':
			last_addr = (*(*c).lines).len - 1;
			free((line *) (*(*c).lines).e[last_addr]);
			(*(*c).lines).e[last_addr] = NULL;
			remove_array_voidstar((*c).lines, last_addr, NULL);
			break;
		case 'c':
			last_addr = (*(*c).circles).len - 1;
			free((circle *) (*(*c).circles).e[last_addr]);
			(*(*c).circles).e[last_addr] = NULL;
			remove_array_voidstar((*c).circles, last_addr, NULL);
			break;
	}
}

void free_sc_constr(sc_constr *c)
{
	free_array_voidstar((*c).points, NULL);
	free_array_voidstar((*c).lines, NULL);
	free_array_voidstar((*c).circles, NULL);
	free((*c).points);
	free((*c).lines);
	free((*c).circles);
	free_array_char(&((*c).history));
}



void sc_constr_line_coords_rem(sc_constr *c, array_double *xs, array_double *ys, array_char *covered, int i)
{
	line *line_i = (line *) (*(*c).lines).e[i];
	int a_i = (*(*line_i).a).addr;
	int b_i = (*(*line_i).b).addr;
	if (!(*covered).e[a_i]) sc_constr_point_coords_rem(c, xs, ys, covered, a_i);
	if (!(*covered).e[b_i]) sc_constr_point_coords_rem(c, xs, ys, covered, b_i);
}

void sc_constr_circle_coords_rem(sc_constr *c, array_double *xs, array_double *ys, array_char *covered, int i)
{
	circle *c_i = (circle *) (*(*c).circles).e[i];
	int pc = (*(*c_i).center).addr;
	int pr = (*(*c_i).radius).addr;
	if (!(*covered).e[pc]) sc_constr_point_coords_rem(c, xs, ys, covered, pc);
	if (!(*covered).e[pr]) sc_constr_point_coords_rem(c, xs, ys, covered, pr);
}

void sc_constr_point_coords_rem(sc_constr *c, array_double *xs, array_double *ys, array_char *covered, int i)
{
	if (covered != NULL && (*covered).e[i]) return;

	point *p = (point *) (*(*c).points).e[i];
	//;;;;printf("Determining coordinates associated with point %d with flag %d\n", i, (*p).flag);
	if ((*p).flag == 1)
	{
		(*xs).e[i] = *((double *) (*p).inc[0]);
		(*ys).e[i] = *((double *) (*p).inc[1]);
		//printf("Setting coordinates for rooted point %d at %g %g\n", i, (*xs).e[i], (*ys).e[i]);
	}
	else
	{
		//;;;printf("Intersection point\n");
		int mode = (*p).flag >> 1;
		if (!(mode & 3))
		{
			// line line intersection
			line *line_1 = (line *) (*p).inc[0];
			line *line_2 = (line *) (*p).inc[1];
			int l1 = (*line_1).addr;
			int l2 = (*line_2).addr;
			if (covered != NULL) 
			{
				sc_constr_line_coords_rem(c, xs, ys, covered, l1);
				sc_constr_line_coords_rem(c, xs, ys, covered, l2);
			}
			int l1a = (*(*line_1).a).addr;
			int l1b = (*(*line_1).b).addr;
			int l2a = (*(*line_2).a).addr;
			int l2b = (*(*line_2).b).addr;
			double l1a_x = (*xs).e[l1a], l1a_y = (*ys).e[l1a], l1b_x = (*xs).e[l1b], l1b_y = (*ys).e[l1b], l2a_x = (*xs).e[l2a], l2a_y = (*ys).e[l2a], l2b_x = (*xs).e[l2b], l2b_y = (*ys).e[l2b];
			char status = line_line_intersection_exp(l1a_x, l1a_y, l1b_x, l1b_y, l2a_x, l2a_y, l2b_x, l2b_y, &((*xs).e[i]), &((*ys).e[i]));
			//;;;printf("Line-line intersection (exp): lines %d and %d -> (%g, %g)\n", l1, l2, (*xs).e[i], (*ys).e[i]);
		}
		else if (((mode & 3) == 1) || ((mode & 3) == 2))
		{
			//;;;printf("\tLine-circle intersection\n");
			// Line circle || circle line
			line *l_;
			circle *c_;
			// RESUME: Check this!
			if (mode & 3 == 1)
			{
				l_ = (line *) (*p).inc[1];
				c_ = (circle *) (*p).inc[0];
			}
			else
			{
				l_ = (line *) (*p).inc[0];
				c_ = (circle *) (*p).inc[1];
			}
			if (covered != NULL)
			{
				sc_constr_line_coords_rem(c, xs, ys, covered, (*l_).addr);
				sc_constr_circle_coords_rem(c, xs, ys, covered, (*c_).addr);
			}
			int l_a_i = (*(*l_).a).addr;
			int l_b_i = (*(*l_).b).addr;
			int c_c_i = (*(*c_).center).addr;
			int c_r_i = (*(*c_).radius).addr;
			double x1, y1, x2, y2;
			char status;
			double ax = (*xs).e[l_a_i], ay = (*ys).e[l_a_i], 
			       bx = (*xs).e[l_b_i], by = (*ys).e[l_b_i], 
			       ccx = (*xs).e[c_c_i], ccy = (*ys).e[c_c_i], 
			       crx = (*xs).e[c_r_i], cry = (*ys).e[c_r_i];
			status = line_circle_intersection_exp(
				ax, ay, bx, by, ccx, ccy, crx, cry,	
				&x1, &y1, &x2, &y2);
			//;; ;;printf("Line-circle intersection (exp): line %d, circle %d, (%g, %g), (%g, %g)\n", (*l_).addr, (*c_).addr, x1, y1, x2, y2);
			double x2mx1 = x2 - x1, y2my1 = y2 - y1;
			bx -= ax;
			by -= ay;
			double test_eval = bx * x2mx1 + by * y2my1;
			if ((*p).flag & LR_MASK)
			{
				if (test_eval > 0)
				{
					(*xs).e[i] = x1;
					(*ys).e[i] = y1;
				}
				else
				{
					(*xs).e[i] = x2;
					(*ys).e[i] = y2;
				}
			}
			else
			{
				if (test_eval > 0)
				{
					(*xs).e[i] = x2;
					(*ys).e[i] = y2;
				}
				else
				{
					(*xs).e[i] = x1;
					(*ys).e[i] = y1;
				}
			}
		}	
		else if (mode & 3 == 3)
		{
			// circle-circle
			//;;;printf("\tCircle-circle intersection\n");
			circle *c1 = (circle *) (*p).inc[0];
			circle *c2 = (circle *) (*p).inc[1];
			if (covered != NULL)
			{
				sc_constr_circle_coords_rem(c, xs, ys, covered, (*c1).addr);
				sc_constr_circle_coords_rem(c, xs, ys, covered, (*c2).addr);
			}
			int c1c_i = (*(*c1).center).addr;
			int c1r_i = (*(*c1).radius).addr;
			int c2c_i = (*(*c2).center).addr;
			int c2r_i = (*(*c2).radius).addr;
			double x1, y1, x2, y2;
			//;; ;;printf("Circle coordinates: (%g %g: %g %g), (%g %g: %g %g)\n", (*xs).e[c1c_i], (*ys).e[c1c_i], (*xs).e[c1r_i], (*ys).e[c1r_i], (*xs).e[c2c_i], (*ys).e[c2c_i], (*xs).e[c2r_i], (*ys).e[c2r_i]);
			circle_circle_intersection_exp(
				(*xs).e[c1c_i], (*ys).e[c1c_i], 
				(*xs).e[c1r_i], (*ys).e[c1r_i], 
				(*xs).e[c2c_i], (*ys).e[c2c_i], 
				(*xs).e[c2r_i], (*ys).e[c2r_i], 
				&x1, &y1, &x2, &y2);
			//;; ;;printf("Circle circle intersection (exp): circles %d and %d, (%g, %g), (%g, %g)\n", (*c1).addr, (*c2).addr, x1, y1, x2, y2);
			char lr_flag = ((*p).flag >> 3) & 1; // RESUME: LR_
			/*double x1ref = x1 - (*xs).e[c1c_i];
			double y1ref = y1 - (*ys).e[c1c_i];
			double x2ref = x2 - (*xs).e[c1c_i];
			double y2ref = y2 - (*ys).e[c1c_i];*/
			if (lr_flag)
			{
				(*xs).e[i] = x2;
				(*ys).e[i] = y2;
			}
			else
			{
				(*xs).e[i] = x1;
				(*ys).e[i] = y1;
			}
			/*
			ave_x = x1ref + x2ref;
			ave_y = y1ref + y2ref;
			double test_eval = ave_x * y1ref - ave_y * x1ref;
			char test_bool = test_eval > 0;
			if (lr_flag)
			{
				if (test_bool)
				{
					(*xs).e[i] = x2;
					(*ys).e[i] = y2;
				}
				else
				{
					(*xs).e[i] = x1;
					(*ys).e[i] = y1;
				}
			}
			else
			{
				if (test_bool)
				{
					(*xs).e[i] = x1;
					(*ys).e[i] = y1;
				}
				else
				{
					(*xs).e[i] = x2;
					(*ys).e[i] = y2;
				}
			}
			*/
		}
	}
	if (covered != NULL) (*covered).e[i] = 1;
}

void sc_constr_points_coords(sc_constr *c, array_double *xs, array_double *ys)
{
	/*
	array_char covered;
	array_char_init(&covered, (*(*c).points).len);
	covered.len = (*(*c).points).len;
	for (int i = 0; i < (*(*c).points).len; i++)
	{
		covered.e[i] = 0;
	}
	*/
	for (int i = 0; i < (*(*c).points).len; i++)
	{
		//sc_constr_point_coords_rem(c, xs, ys, &covered, i);
		sc_constr_point_coords_rem(c, xs, ys, NULL, i);
	}
	//free_array_char(&covered);
}

char sc_constr_check_singular_lc(sc_constr *sc, int l, int c)
{
	// check if 'l' passes through the center of 'c' (in general it might also be worth checking 
	// 	whether 'l' is tangent to 'c')
	// Also, in the future it would be ideal to use a polynomial API for these calculations
	line *l_ = (line *) (*(*sc).lines).e[l];
	circle *c_ = (circle *) (*(*sc).circles).e[c];
	if ((*l_).a == (*c_).center || (*l_).b == (*c_).center) return 1;
	else
	{
		double ax, ay, bx, by, cx, cy, rx, ry;
		line_coords(l_, &ax, &ay, &bx, &by);
		point_coords((*c_).center, &cx, &cy);
		bx -= ax;
		by -= ay;
		rx = cx - ax;
		ry = cy - ay;
		if (rx * by - ry * bx == 0) return 1;
	}
	return 0;
}

void add_point_sc_constr_ll(sc_constr *c, int l1, int l2)
{
	point *p = (point *) calloc(1, sizeof(point));
	(*p).inc[0] = (*(*c).lines).e[l1];
	(*p).inc[1] = (*(*c).lines).e[l2];
	(*p).flag = 0;
	add2sc_constr(c, (void *) p, 'p');	
}

void add_point_sc_constr_lc(sc_constr *c, int l_, int c_, char lr)
{
	point *p = (point *) calloc(1, sizeof(point));
	(*p).inc[0] = (*(*c).lines).e[l_];
	(*p).inc[1] = (*(*c).circles).e[c_];
	(*p).flag = 4 | (lr << 3);
	add2sc_constr(c, (void *) p, 'p');
}

void add_point_sc_constr_cc(sc_constr *c, int c1, int c2, char lr)
{
	point *p = (point *) calloc(1, sizeof(point));
	(*p).inc[0] = (*(*c).circles).e[c1];
	(*p).inc[1] = (*(*c).circles).e[c2];
	(*p).flag = 6 | (lr << 3);
	add2sc_constr(c, (void *) p, 'p');
}

void add_line_sc_constr_pp(sc_constr *c, int a_, int b_)
{
	line *l_ = (line *) calloc(1, sizeof(line));
	(*l_).a = (point *) (*(*c).points).e[a_];
	(*l_).b = (point *) (*(*c).points).e[b_];
	add2sc_constr(c, (void *) l_, 'l');
}

void add_circle_sc_constr_pp(sc_constr *c, int c_, int r_)
{
	circle *crc = (circle *) calloc(1, sizeof(circle));
	(*crc).center = (point *) (*(*c).points).e[c_];
	(*crc).radius = (point *) (*(*c).points).e[r_];
	add2sc_constr(c, (void *) crc, 'c');
	double cx, cy, rx, ry;
	point_coords((*crc).center, &cx, &cy);
	point_coords((*crc).radius, &rx, &ry);
	//printf("New circle from %d to %d: center = %g %g, radial point = %g %g\n", c_, r_, cx, cy, rx, ry);
}

char sc_constr_check_inter_ll(sc_constr *c, int l1, int l2)
{
	line *ln1 = (line *) (*(*c).lines).e[l1];
	line *ln2 = (line *) (*(*c).lines).e[l2];
	double a1x, a1y, b1x, b1y, a2x, a2y, b2x, b2y;
	point_coords((*ln1).a, &a1x, &a1y);
	point_coords((*ln1).b, &b1x, &b1y);
	point_coords((*ln2).a, &a2x, &a2y);
	point_coords((*ln2).b, &b2x, &b2y);
	b1x -= a1x;
	b1y -= a1y;
	b2x -= a2x;
	b2y -= a2y;
	return b1x * b2y - b1y * b2x != 0;
}

char sc_constr_check_inter_lc(sc_constr *c, int l_, int c_)
{
	line *ln_ = (line *) (*(*c).lines).e[l_];
	circle *crc = (circle *) (*(*c).circles).e[c_];
	double ax, ay, bx, by, ccx, ccy, crx, cry;
	point_coords((*ln_).a, &ax, &ay);
	point_coords((*ln_).b, &bx, &by);
	point_coords((*crc).center, &ccx, &ccy);
	point_coords((*crc).radius, &crx, &cry);
	bx -= ax;
	by -= ay;
	ax -= ccx;
	ay -= ccy;
	crx -= ccx;
	cry -= ccy;
	double rsq = crx * crx + cry * cry, asq = ax * ax + ay * ay;
	double adotv = ax * bx + ay * by;
	double vsq = bx * bx + by * by;
	double B = 2 * adotv;
	double discr = B * B - 4 * vsq * (asq - rsq);
	return discr >= 0;
	// (ax + vx t - ccx)^2 + (ay + vy t - ccy)^2 = rsq
	// a'^2 + 2 ax' vx t + 2 ay' vy t + v^2 t^2 = rsq
}

char sc_constr_check_inter_cc(sc_constr *c, int c1, int c2)
{
	circle *c1_ = (circle *) (*(*c).circles).e[c1];
	circle *c2_ = (circle *) (*(*c).circles).e[c2];
	double c1x, c1y, r1x, r1y, c2x, c2y, r2x, r2y;
	point_coords((*c1_).center, &c1x, &c1y);
	point_coords((*c1_).radius, &r1x, &r1y);
	point_coords((*c2_).center, &c2x, &c2y);
	point_coords((*c2_).radius, &r2x, &r2y);
	r1x -= c1x;
	r1y -= c1y;
	r2x -= c2x;
	r2y -= c2y;
	double r1sq = r1x * r1x + r1y * r1y, r2sq = r2x * r2x + r2y * r2y;
	c1x -= c2x;
	c1y -= c2y;
	double c1xsq = c1x * c1x, c1ysq = c1y * c1y, c1sq = c1xsq + c1ysq, r2sqmr1sq = r2sq - r1sq;
	double B = -2 * r2sqmr1sq * c1y * c1x, C = 0.25 * r2sqmr1sq * r2sqmr1sq;
	if (c1xsq > c1ysq)
	{

		// x = (r2sq - r1sq - 2 y c1y) / (2 c1x)
		// y^2 c1^2 - 2 y (r2sqmr1sq) c1y c1x + (r2sqmr1sq)^2/(4) = r2sq c1x^2
		C -= c1xsq * r2sq;
		return (B * B - 4 * c1sq * C) > 0;
	}
	else
	{
		C -= c1ysq * r2sq;
		return (B * B - 4 * c1sq * C) > 0;
		// y = (r2sq - r1sq - 2 x c1x) / (2 c1y)
		// x^2 c1^2 - 2 x c1y c1x (r2sqmr1sq) + (r2sqmr1sq)^2 / (4) = r2sq c1y^2
		
	}
	// (x - c1x)^2 + (y - c1y)^2 = r1sq
	// (x)^2 + (y)^2 = r2sq
	// -2 x c1x - 2 y c1y = r1sq - r2sq
	
}

void sc_constr_print_point_coords(sc_constr *sc)
{
	printf("Point coordinates: (diagnostic) \n");
	for (int i = 0; i < (*(*sc).points).len; i++)
	{
		point *p = (point *) (*(*sc).points).e[i];
		double px, py;
		point_coords(p, &px, &py);
		printf("%d: %g %g\n", i, px, py);
	}
}

void sc_constr_print_circle_coords(sc_constr *sc)
{
	printf("Circle coordinates: (diagnostic)\n");
	for (int i = 0; i < (*(*sc).circles).len; i++)
	{
		circle *c_ = (circle *) (*(*sc).circles).e[i];
		double cx, cy, rx, ry;
		point_coords((*c_).center, &cx, &cy);
		point_coords((*c_).radius, &rx, &ry);
		rx -= cx;
		ry -= cy;
		printf("%d: %g %g: %g %g r = %g\n", i, cx, cy, rx, ry, sqrt(rx * rx + ry * ry));
	}
}

void sc_constr_print_line_coords(sc_constr *sc)
{
	printf("Line coordinates: (diagnostic)\n");
	for (int i = 0; i < (*(*sc).lines).len; i++)
	{
		line *l_ = (line *) (*(*sc).lines).e[i];
		double ax, ay, bx, by;
		point_coords((*l_).a, &ax, &ay);
		point_coords((*l_).b, &bx, &by);
		printf("%d: %g %g: %g %g\n", i, ax, ay, bx, by);
	}
}
