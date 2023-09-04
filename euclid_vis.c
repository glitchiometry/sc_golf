#include "euclid_vis.h"

#define TWO_PI 6.283185307179586
#define NPTS_CIRCLE 100
#define POINT_TEMPLATE_LEN 8
int point_template_x[POINT_TEMPLATE_LEN] = {2, 1, -1, -2, -2, -1, 1, 2};
int point_template_y[POINT_TEMPLATE_LEN] = {1, 2, 2, 1, -1, -2, -2, -1};

char ival_contains(double *xbnds, double x)
{
	return xbnds[0] <= x && xbnds[1] >= x;
}

char rect_contains(double *xbnds, double *ybnds, double x, double y)
{
	return ival_contains(xbnds, x) && ival_contains(ybnds, y);
}

void sc_constr_interface_init(sc_constr_interface *scci, sc_constr *sc, double *xbnds, double *ybnds, int scr_len_x, int scr_len_y)
{
	printf("Initializing sc_constr_interface over rectangle [%g, %g]x[%g, %g]\n", xbnds[0], xbnds[1], ybnds[0], ybnds[1]);
	(*scci).sc = sc;
	(*scci).xbnds[0] = xbnds[0];
	(*scci).xbnds[1] = xbnds[1];
	(*scci).ybnds[0] = ybnds[0];
	(*scci).ybnds[1] = ybnds[1];
	(*scci).screen_len_x = scr_len_x;
	(*scci).screen_len_y = scr_len_y;
	double inv_wid_x = (*scci).screen_len_x * (1. / ((*scci).xbnds[1] - (*scci).xbnds[0]));
	double inv_wid_y = (*scci).screen_len_y * (1. / ((*scci).ybnds[1] - (*scci).ybnds[0]));
	array_int_init(&((*scci).points_x_int), (*(*sc).points).len);
	array_int_init(&((*scci).points_y_int), (*(*sc).points).len);
	(*scci).points_x_int.len = (*(*sc).points).len;
	(*scci).points_y_int.len = (*(*sc).points).len;
	array_double_init(&((*scci).points_x), (*(*sc).points).len);
	array_double_init(&((*scci).points_y), (*(*sc).points).len);
	(*scci).points_x.len = (*(*sc).points).len;
	(*scci).points_y.len = (*(*sc).points).len;
	array_char_init(&((*scci).points_active), (*(*sc).points).len);
	(*scci).points_active.len = (*(*sc).points).len;
	array_int_init(&((*scci).active_points), 1);
	array_voidstar_init(&((*scci).circ_data), (*(*sc).circles).len);
	array_int_init(&((*scci).active_circles), 1);
	array_voidstar_init(&((*scci).line_data), (*(*sc).lines).len);
	array_int_init(&((*scci).active_lines), 1);
	sc_constr_points_coords(sc, &((*scci).points_x), &((*scci).points_y));
	// Add discretized points
	for (int i = 0; i < (*(*sc).points).len; i++)
	{
		if (rect_contains((*scci).xbnds, (*scci).ybnds, (*scci).points_x.e[i], (*scci).points_y.e[i]))
		{
			(*scci).points_active.e[i] = 1;
			add2array_int(&((*scci).active_points), i);
		}
		else
		{
			(*scci).points_active.e[i] = 0;
		}
		(*scci).points_x_int.e[i] = (int) (((*scci).points_x.e[i] - (*scci).xbnds[0]) * inv_wid_x);
		(*scci).points_y_int.e[i] = (int) (((*scci).points_y.e[i] - (*scci).ybnds[0]) * inv_wid_y);
		printf("Adding discretized point: %g %g -> %d %d\n", (*scci).points_x.e[i], (*scci).points_y.e[i], (*scci).points_x_int.e[i], (*scci).points_y_int.e[i]);
	}
	// Add line data
	for (int i = 0; i < (*(*sc).lines).len; i++)
	{
		line *l_i = (line *) (*(*sc).lines).e[i];
		double ax, ay, bx, by;
		int x_upper_int, x_lower_int, y_upper_int, y_lower_int;
		//	void line_rectangle_intersection_int(double ax, double ay, double bx, double by, sc_constr_interface *scci, double inv_wid_x, double inv_wid_y, int *xx, int *yy, int *x__, int *y__)

		line_coords(l_i, &ax, &ay, &bx, &by);
		// Determine intersections on top and bottom boundaries
		int x_a_int, y_a_int, x_b_int, y_b_int;
		char status;
		line_rectangle_intersection_int(ax, ay, bx, by, scci, inv_wid_x, inv_wid_y, &x_a_int, &y_a_int, &x_b_int, &y_b_int, &status);
		line_render_data *l_ = (line_render_data *) calloc(1, sizeof(line_render_data));
		// Make sure that lines are oriented consistently (in case prudence requires adding a graphical representation.)
		(*l_).a.x = x_a_int;
		(*l_).a.y = y_a_int;
		(*l_).b.x = x_b_int;
		(*l_).b.y = y_b_int;
		(*l_).vis = status == 1;
		add2array_voidstar(&((*scci).line_data), (void *) l_);
		if ((*l_).vis)
		{
			add2array_int(&((*scci).active_lines), i);
		}
	}
	// Add circle data
	for (int i = 0; i < (*(*sc).circles).len; i++)
	{
		circle *c_i = (circle *) (*(*sc).circles).e[i];
		double cx, cy, cr;
		circle_coords(c_i, &cx, &cy, &cr);
		circle_render_data *cdat = (circle_render_data *) calloc(1, sizeof(circle_render_data));
		circle_render_data_init(cdat, scci, cx, cy, cr);
		add2array_voidstar(&((*scci).circ_data), (void *) cdat);
		if ((*cdat).vis)
		{
			add2array_int(&((*scci).active_circles), i);
		}
	}
}

void add_point_sc_constr_interface(sc_constr_interface *scci, int point_addr)
{
	double inv_wid_x = (*scci).screen_len_x / ((*scci).xbnds[1] - (*scci).xbnds[0]);
	double inv_wid_y = (*scci).screen_len_y / ((*scci).ybnds[1] - (*scci).ybnds[0]);
	if (point_addr < (*(*(*scci).sc).points).len && -1 < point_addr)
	{
		add_mem_array_double_until(&((*scci).points_x), point_addr);
		add_mem_array_double_until(&((*scci).points_y), point_addr);
		add_mem_array_int_until(&((*scci).points_x_int), point_addr);
		add_mem_array_int_until(&((*scci).points_y_int), point_addr);
		add_mem_array_char_until(&((*scci).points_active), point_addr);
		if (point_addr < (*scci).points_x.len) {}
		else
		{
			int new_len = point_addr + 1;
			(*scci).points_x.len = new_len;
			(*scci).points_y.len = new_len;
			(*scci).points_x_int.len = new_len;
			(*scci).points_y_int.len = new_len;
			(*scci).points_active.len = new_len;
		}
		int x_, y_;
		sc_constr_point_coords_rem((*scci).sc, &((*scci).points_x), &((*scci).points_y), NULL, point_addr);
		x_ = (int) (((*scci).points_x.e[point_addr] - (*scci).xbnds[0]) * inv_wid_x);
		y_ = (int) (((*scci).points_y.e[point_addr] - (*scci).ybnds[0]) * inv_wid_y);
		(*scci).points_x_int.e[point_addr] = x_;
		(*scci).points_y_int.e[point_addr] = y_;
		if (rect_contains((*scci).xbnds, (*scci).ybnds, (*scci).points_x.e[point_addr], (*scci).points_y.e[point_addr]))
		{
			(*scci).points_active.e[point_addr] = 1;
			add2array_int(&((*scci).active_points), point_addr);
		}
		else
		{
			(*scci).points_active.e[point_addr] = 0;
		}
		
	}
}

void add_line_sc_constr_interface(sc_constr_interface *scci, int line_addr)
{
	line *l_i;
	if (line_addr < (*(*(*scci).sc).lines).len) l_i = (line *) (*(*(*scci).sc).lines).e[line_addr];
	else 
	{
		printf("Error: attempting to access non-existent line address %d of %d\n", line_addr, (*(*(*scci).sc).lines).len);
		exit(EXIT_FAILURE);
	}
	double ax, ay, bx, by;
	int x_upper_int, x_lower_int, y_upper_int, y_lower_int;
	//	void line_rectangle_intersection_int(double ax, double ay, double bx, double by, sc_constr_interface *scci, double inv_wid_x, double inv_wid_y, int *xx, int *yy, int *x__, int *y__)
	int a_i = (*(*l_i).a).addr;
	int b_i = (*(*l_i).b).addr;
	ax = (*scci).points_x.e[a_i];
	ay = (*scci).points_y.e[a_i];
	bx = (*scci).points_x.e[b_i];
	by = (*scci).points_y.e[b_i];
	printf("Initializing line render data: %g %g, %g %g\n", ax, ay, bx, by);
	// Determine intersections on top and bottom boundaries
	int x_a_int, y_a_int, x_b_int, y_b_int;
	double inv_wid_x = (*scci).screen_len_x / ((*scci).xbnds[1] - (*scci).xbnds[0]);
	double inv_wid_y = (*scci).screen_len_y / ((*scci).ybnds[1] - (*scci).ybnds[0]);
	char status;
	line_rectangle_intersection_int(ax, ay, bx, by, scci, inv_wid_x, inv_wid_y, &x_a_int, &y_a_int, &x_b_int, &y_b_int, &status);
	line_render_data *l_ = (line_render_data *) calloc(1, sizeof(line_render_data));
	printf("Setting line render data coordinates to %d %d, %d %d\n", x_a_int, y_a_int, x_b_int, y_b_int);
	(*l_).a.x = x_a_int;
	(*l_).a.y = y_a_int;
	(*l_).b.x = x_b_int;
	(*l_).b.y = y_b_int;
	(*l_).vis = status > 0;
	add2array_voidstar(&((*scci).line_data), (void *) l_);
	if ((*l_).vis)
	{
		add2array_int(&((*scci).active_lines), line_addr);
	}
}

	

void free_sc_constr_interface(sc_constr_interface *scci)
{
	for (int i = 0; i < (*scci).circ_data.len; i++)
	{
		circle_render_data *cdat = (circle_render_data *) (*scci).circ_data.e[i];
		free_circle_render_data(cdat);
	}
	free_array_double(&((*scci).points_x));
	free_array_double(&((*scci).points_y));
	free_array_int(&((*scci).points_x_int));
	free_array_int(&((*scci).points_y_int));
	free_array_voidstar(&((*scci).circ_data), NULL);
	free_array_voidstar(&((*scci).line_data), NULL);
	free_array_char(&((*scci).points_active));
	free_array_int(&((*scci).active_points));
	free_array_int(&((*scci).active_lines));
	free_array_int(&((*scci).active_circles));
}

void sc_constr_interface_resize(sc_constr_interface *scci, sc_constr *sc, double *xbnds, double *ybnds, int scr_len_x, int scr_len_y)
{
	free_sc_constr_interface(scci);

	sc_constr_interface_init(scci, sc, xbnds, ybnds, scr_len_x, scr_len_y);
}

void add_circle_sc_constr_interface(sc_constr_interface *scci, int circle_addr)
{
	if (circle_addr < (*(*(*scci).sc).circles).len) 
	{
		circle *c_;
		c_ = (circle *) (*(*(*scci).sc).circles).e[circle_addr];
		int c_i = (*(*c_).center).addr;
		int r_i = (*(*c_).radius).addr;
		double cx = (*scci).points_x.e[c_i], cy = (*scci).points_y.e[c_i], rx = (*scci).points_x.e[r_i], ry = (*scci).points_y.e[r_i];
		rx -= cx;
		ry -= cy;
		double r = sqrt(rx * rx + ry * ry);
		circle_render_data *cdat_ = (circle_render_data *) calloc(1, sizeof(circle_render_data));
		circle_render_data_init(cdat_, scci, cx, cy, r);
		if ((*cdat_).vis)
		{
			add2array_int(&((*scci).active_circles), circle_addr);
		}
		add2array_voidstar(&((*scci).circ_data), (void *) cdat_);
	}
	else
	{
		printf("Attempting to access a circle that is not yet registered in associated sc_constr instance\n");
		exit(EXIT_FAILURE);
	}
}

void circle_render_data_init(circle_render_data *cdata, sc_constr_interface *scci, double cx, double cy, double r)
{
	double xmax = cx + r;
	double xmin = cx - r;
	double ymax = cy + r;
	double ymin = cy - r;
	if (xmax > (*scci).xbnds[0] && xmin < (*scci).xbnds[1] && ymax > (*scci).ybnds[0] && ymin < (*scci).ybnds[1]) {}
	else
	{
		(*cdata).bdry = NULL;
		(*cdata).len = 0;
		(*cdata).vis = 0;
		return;
	}
	(*cdata).vis = 0;
	int xs[NPTS_CIRCLE];
	int ys[NPTS_CIRCLE];
	int len = 0;
	double dtheta = TWO_PI / (NPTS_CIRCLE - 1);
	double theta = dtheta;
	char resume_flag = 0;
	double inv_wid_x = (*scci).screen_len_x * (1.0 / ((*scci).xbnds[1] - (*scci).xbnds[0]));
	double inv_wid_y = (*scci).screen_len_y * (1.0 / ((*scci).ybnds[1] - (*scci).ybnds[0]));
	for (int i = 0; i < NPTS_CIRCLE; i++)
	{
		double x = r * cos(theta) + cx;
		double y = r * sin(theta) + cy;
		if ((*cdata).vis) {}
		else
		{
			if (x < (*scci).xbnds[0] || x > (*scci).xbnds[1] || y < (*scci).ybnds[0] || y > (*scci).ybnds[1]) {}
			else
			{
				(*cdata).vis = 1;
			}
		}
		/*
		char test_xmin = x < (*scci).xbnds[0] ? 1 : 0;
		char test_xmax = x > (*scci).xbnds[1] ? 1 : 0;
		char test_ymin = y < (*scci).ybnds[0] ? 1 : 0;
		char test_ymax = y > (*scci).ybnds[1] ? 1 : 0;
		char next_resume_flag = test_xmin | test_xmax << 1 | test_ymin << 2 | test_ymax << 3;
		if (!next_resume_flag)
		{
			if (!resume_flag)
			{}
			else
			{
				double last_x = r * cos(theta - dtheta) + cx, last_y = r * sin(theta - dtheta) + cy;
				int *aux_x;
				int *aux_y;
				segment_rectangle_intersection_int(x, y, 0, last_x, last_y, resume_flag, scci, inv_wid_x, inv_wid_y, &(xs[len]), &(ys[len]), aux_x, aux_y);
				len += 1;
				resume_flag = 0;
			}
			xs[len] = (int) ((x - (*scci).xbnds[0]) * inv_wid_x);
			ys[len] = (int) ((y - (*scci).ybnds[0]) * inv_wid_y);
			len += 1;
		}
		else
		{
			if (!resume_flag)
			{
				// Determine the boundary intersection
				double last_x = cx + r * cos(theta - dtheta), last_y = cy + r * sin(theta - dtheta);
				int *aux_x;
				int *aux_y;
				segment_rectangle_intersection_int(x, y, next_resume_flag, last_x, last_y, 0, scci, inv_wid_x, inv_wid_y, &(xs[len]), &(ys[len]), aux_x, aux_y);

				len += 1;
			}
			resume_flag = next_resume_flag;
		}
		*/
		xs[len] = (int) ((x - (*scci).xbnds[0]) * inv_wid_x);
		ys[len] = (int) ((y - (*scci).ybnds[0]) * inv_wid_y);
		len += 1;
		theta += dtheta;
	}
	// Allocate SDL points
	(*cdata).bdry = (SDL_Point *) calloc(len, sizeof(SDL_Point));
	for (int i = 0; i < len; i++)
	{
		(*cdata).bdry[i].x = xs[i];
		(*cdata).bdry[i].y = ys[i];
	}
	(*cdata).center.x = (int) ((cx - (*scci).xbnds[0]) * inv_wid_x);
	(*cdata).center.y = (int) ((cy - (*scci).ybnds[0]) * inv_wid_y);
	(*cdata).len = len;
}

void circle_render_data_init_test(circle_render_data *cdata, sc_constr_interface *scci, double cx, double cy, double r)
{
	(*cdata).center.x = cx;
	(*cdata).center.y = cy;
	double xmin = cx - r;
	double xmax = cx + r;
	double ymin = cy - r;
	double ymax = cy + r;
	char ovl_state = 0;
	ovl_state |= (xmin < (*scci).xbnds[0]);
	ovl_state |= (xmax > (*scci).xbnds[1]) << 1;
	ovl_state |= (ymin < (*scci).ybnds[0]) << 2;
	ovl_state |= (ymax > (*scci).ybnds[1]) << 3;
	if (ovl_state)
	{
		if (ovl_state & 1 && ovl_state & 4)
		{
			// check if lower left corner is within the circle
			double delx = (*scci).xbnds[0] - cx;
			double dely = (*scci).ybnds[0] - cy;
			double delsq = delx * delx + dely * dely;
		}
		if (ovl_state & 1 && ovl_state & 8)
		{
			// check if upper left corner is within the circle
		}
		if (ovl_state & 2 && ovl_state & 4)
		{
			// check lower right corner
		}
		if (ovl_state & 2 && ovl_state & 8)
		{
			// check upper right corner
		}
	}
}

void free_circle_render_data(circle_render_data *cdata)
{
	free((*cdata).bdry);
}

void line_rectangle_intersection_int(double ax, double ay, double bx, double by, sc_constr_interface *scci, double inv_wid_x, double inv_wid_y, int *xx, int *yy, int *x__, int *y__, char *status)
{
	//;;;printf("Line-rectangle intersection: points (%g, %g), (%g, %g), rectangle [%g, %g]x[%g, %g]\n", ax, ay, bx, by, (*scci).xbnds[0], (*scci).xbnds[1], (*scci).ybnds[0], (*scci).ybnds[1]);
	bx -= ax;
	by -= ay;
	double t_upper = ((*scci).ybnds[1] - ay) / by;
	double t_lower = ((*scci).ybnds[0] - ay) / by;
	double t_right = ((*scci).xbnds[1] - ax) / bx;
	double t_left = ((*scci).xbnds[0] - ax) / bx;
	double x_upper = ax + t_upper * bx, x_lower = ax + t_lower * bx, y_right = ay + t_right * by, y_left = ay + t_left * by;
	int *xs[2] = {xx, x__};
	int *ys[2] = {yy, y__};
	double ts[2];
	char count = 0;
	if (ival_contains((*scci).xbnds, x_upper))
	{
		(*xs[count]) = (int) ((x_upper - (*scci).xbnds[0]) * inv_wid_x);
		(*ys[count]) = (*scci).screen_len_y;
		printf("\t %d %d\n", *(xs[count]), *(ys[count]));
		ts[count] = t_upper;
		count += 1;
	}
	if (ival_contains((*scci).xbnds, x_lower))
	{
		(*xs[count]) = (int) ((x_lower - (*scci).xbnds[0]) * inv_wid_x);
		(*ys[count]) = 0;
		printf("\t %d %d\n", *(xs[count]), *(ys[count]));
		ts[count] = t_lower;
		count += 1;
	}
	if (ival_contains((*scci).ybnds, y_right))
	{
		(*xs[count]) = (*scci).screen_len_x;
		(*ys[count]) = (int) ((y_right - (*scci).ybnds[0]) * inv_wid_y);
		printf("\t %d %d\n", *(xs[count]), *(ys[count]));
		ts[count] = t_right;
		count += 1;
	}
	if (ival_contains((*scci).ybnds, y_left))
	{
		(*xs[count]) = 0;
		(*ys[count]) = (int) ((y_left - (*scci).ybnds[0]) * inv_wid_y);
		printf("\t %d %d\n", *(xs[count]), *(ys[count]));
		ts[count] = t_left;
	}
	if (ts[0] <= ts[1]) {}
	else
	{
		int aux = (*xx);
		(*xx) = (*x__);
		(*x__) = aux;
		aux = (*yy);
		(*yy) = (*y__);
		(*y__) = aux;
	}
	(*status) = count;
}

// FIX THIS!
void segment_rectangle_intersection_int(double ax, double ay, char bdry_a, double bx, double by, char bdry_b, sc_constr_interface *scci, double inv_wid_x, double inv_wid_y, int *xx, int *yy, int *x__, int *y__)
{
	bx -= ax;
	by -= ay;
	if (!(bdry_a && bdry_b))
	{
		x__ = y__ = NULL;
		bdry_b = bdry_a ? bdry_a : bdry_b;
		// Determine the boundary intersection
		if (bdry_b & 4)
		{
			// Overlap with lower boundary
			double x_ = ax + bx * ((*scci).ybnds[0] - ay) / by;
			if (ival_contains((*scci).xbnds, x_))
			{
				(*xx) = (int) ((x_ - (*scci).xbnds[0]) * inv_wid_x);
				(*yy) = 0;
			}
			else
			{
				double y_;
				if (x_ > (*scci).xbnds[1])
				{
					y_ = ay + by * ((*scci).xbnds[1] - ax) / bx;
					(*xx) = (*scci).screen_len_x;
					(*yy) = (int) ((y_ - (*scci).ybnds[0]) * inv_wid_y);
				}
				else
				{
					y_ = ay + by * ((*scci).xbnds[0] - ax) / bx;
					(*xx) = 0;
					(*yy) = (int) ((y_ - (*scci).ybnds[0]) * inv_wid_y);
				}
			}
		}
		else if (bdry_b & 8)
		{
			double x_ = ax + bx * ((*scci).ybnds[1] - ay) / by;
			if (ival_contains((*scci).xbnds, x_))
			{
				(*xx) = (int) ((x_ - (*scci).xbnds[0]) * inv_wid_x);
				(*yy) = (*scci).screen_len_y;
			}
			else
			{
				double y_;
				if (x_ > (*scci).xbnds[1])
				{
					y_ = ay + by * ((*scci).xbnds[1] - ax) / bx;
					(*xx) = (*scci).screen_len_x;
					(*yy) = (int) ((y_ - (*scci).ybnds[1]) * inv_wid_y);
				}
				else
				{
					y_ = ay + by * ((*scci).xbnds[0] - ax) / bx;
					(*xx) = 0;
					(*yy) = (int) ((y_ - (*scci).ybnds[1]) * inv_wid_y);
				}
			}
		}
		else if (bdry_b & 1)
		{
			double y_ = ay + by * ((*scci).xbnds[0] - ax) / bx;
			if (ival_contains((*scci).ybnds, y_))
			{
				(*yy) = (int) ((y_ - (*scci).ybnds[0]) * inv_wid_y);
				(*xx) = 0;
			}
			else
			{
				printf("This message should not appear...Okay you can leave now...No really, the corner cases were covered in the last part...\n");
			}
		}
		else if (bdry_b & 2)
		{
			double y_ = ay + by * ((*scci).xbnds[1] - ax) / bx;
			if (ival_contains((*scci).ybnds, y_))
			{
				(*yy) = (int) ((y_ - (*scci).ybnds[1]) * inv_wid_y);
				(*xx) = (*scci).screen_len_x;
			}
			else
			{
				printf("This message should not appear...Okay you can leave now...No really, the corner cases were covered in the last part...\n");
			}
		}
	}
	else
	{
		// This case shouldn't appear unless the program samples line segments (instead of points.)
		// 	The line segment might cross a corner of the rectangle, depending on the values of 
		//	bdry_a and bdry_b
		// (ax, ay) + (bx, by) t
		//	ay + by t = y1, t = 
		double t_upper = ((*scci).ybnds[1] - ay) / by;
		double t_lower = ((*scci).ybnds[0] - ay) / by;
		double t_right = ((*scci).xbnds[1] - ax) / bx;
		double t_left = ((*scci).xbnds[0] - ax) / bx;
		double x_upper, x_lower, y_right, y_left;
		int *xs[2] = {xx, x__};
		int *ys[2] = {yy, y__};
		int count = 0;
		if (0 <= t_upper && t_upper <= 1.0)
		{
			x_upper = ax + t_upper * bx;
			if ((*scci).xbnds[0] <= x_upper && x_upper <= (*scci).xbnds[1])
			{
				(*xs[count]) = (int) ((x_upper - (*scci).xbnds[0]) * inv_wid_x);
				(*ys[count]) = (*scci).screen_len_y;
				count += 1;
			}
		}
		if (0 <= t_lower && t_lower <= 1.0)
		{
			x_lower = ax + t_lower * bx;
			if ((*scci).xbnds[0] <= x_lower && x_lower <= (*scci).xbnds[1])
			{
				(*xs[count]) = (int) ((x_lower - (*scci).xbnds[0]) * inv_wid_x);
				(*ys[count]) = 0;
				count += 1;
			}
		}
		if (0 <= t_right && t_right <= 1.0)
		{
			y_right = ay + t_right * by;
			if ((*scci).ybnds[0] <= y_right && y_right <= (*scci).ybnds[1])
			{
				(*xs[count]) = (*scci).screen_len_x;
				(*ys[count]) = (int) ((y_right - (*scci).ybnds[0]) * inv_wid_y);
				count += 1;
			}
		}
		if (0 <= t_left && t_left <= 1.0)
		{
			y_left = ay + t_left * by;
			if ((*scci).ybnds[0] <= y_left && y_left <= (*scci).ybnds[1])
			{
				(*xs[count]) = 0;
				(*ys[count]) = (int) ((y_right - (*scci).ybnds[0]) * inv_wid_y);
			}
		}
	}
}

void render_circle(circle_render_data *cdata, SDL_Renderer *rndrr)
{
	SDL_RenderDrawLines(rndrr, (*cdata).bdry, (*cdata).len);
}

void render_line(line_render_data *ldata, SDL_Renderer *rndrr)
{
	SDL_RenderDrawLine(rndrr, (*ldata).a.x, (*ldata).a.y, (*ldata).b.x, (*ldata).b.y);
}

void render_point(int x, int y, SDL_Renderer *rndrr)
{
	SDL_Point marker[POINT_TEMPLATE_LEN];
	for (int i = 0; i < POINT_TEMPLATE_LEN; i++)
	{
		marker[i].x = point_template_x[i] + x;
		marker[i].y = point_template_y[i] + y;
	}
	SDL_RenderDrawLines(rndrr, marker, POINT_TEMPLATE_LEN);
}

void render_sc_constr(sc_constr_interface *scci, SDL_Renderer *rndrr)
{
	for (int i = 0; i < (*scci).active_points.len; i++)
	{
		int pi = (*scci).active_points.e[i];
		render_point((*scci).points_x_int.e[pi], (*scci).points_y_int.e[pi], rndrr);
	}
	for (int i = 0; i < (*scci).active_lines.len; i++)
	{
		int li = (*scci).active_lines.e[i];
		render_line((line_render_data *) (*scci).line_data.e[li], rndrr);
	}
	for (int i = 0; i < (*scci).active_circles.len; i++)
	{
		int ci = (*scci).active_circles.e[i];
		render_circle((circle_render_data *) (*scci).circ_data.e[ci], rndrr);
	}
}

void render_sc_constr_hlts(sc_constr_interface *scci, array_char *hltd_points, array_char *hltd_lines, array_char *hltd_circles, SDL_Renderer *rndrr)
{
	for (int i = 0; i < (*scci).active_points.len; i++)
	{
		int pi = (*scci).active_points.e[i];
		if (!(*hltd_points).e[pi]) render_point((*scci).points_x_int.e[pi], (*scci).points_y_int.e[pi], rndrr);
		else
		{
			SDL_SetRenderDrawColor(rndrr, 255, 0, 0, 0);
			render_point((*scci).points_x_int.e[pi], (*scci).points_y_int.e[pi], rndrr);
			SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 0);
		}
	}
	for (int i = 0; i < (*scci).active_lines.len; i++)
	{
		int li = (*scci).active_lines.e[i];
		if (!(*hltd_lines).e[li]) render_line((line_render_data *) (*scci).line_data.e[li], rndrr);
		else
		{
			SDL_SetRenderDrawColor(rndrr, 0, 255, 0, 0);
			render_line((line_render_data *) (*scci).line_data.e[li], rndrr);
			SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 0);
		}
	}
	for (int i = 0; i < (*scci).active_circles.len; i++)
	{
		int ci = (*scci).active_circles.e[i];
		if (!(*hltd_circles).e[ci]) render_circle((circle_render_data *) (*scci).circ_data.e[ci], rndrr);
		else
		{
			SDL_SetRenderDrawColor(rndrr, 0, 0, 255, 0);
			render_circle((circle_render_data *) (*scci).circ_data.e[ci], rndrr);
			SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 0);
		}
	}
}

void sc_constr_interface_remove_last_point(sc_constr_interface *scci)
{
	remove_array_int(&((*scci).points_x_int), (*scci).points_x_int.len - 1);
	remove_array_int(&((*scci).points_y_int), (*scci).points_x_int.len);
	remove_array_double(&((*scci).points_x), (*scci).points_x_int.len);
	remove_array_double(&((*scci).points_y), (*scci).points_x_int.len);
	if ((*scci).points_active.e[(*scci).points_x.len])
	{
		int api = (*scci).active_points.len;
		char found = 0;
		while (api > 0)
		{
			api -= 1;
			if ((*scci).active_points.e[api] == (*scci).points_x.len)
			{
				remove_array_int(&((*scci).active_points), api);
				found = 1;
				break;
			}
		}
		if (found) {}
		else
		{
			printf("Something weird happened in sc_constr_interface_remove_last_point! Fix it ASAP\n");
			exit(EXIT_FAILURE);
		}
	}
	remove_array_char(&((*scci).points_active), (*scci).points_x.len);
}

void sc_constr_interface_remove_last_line(sc_constr_interface *scci)
{
	int last_ = (*scci).line_data.len - 1;
	line_render_data *l_ = (line_render_data *) (*scci).line_data.e[last_];
	if ((*l_).vis)
	{
		int ali = (*scci).active_lines.len;
		char found = 0;
		while (ali > 0)
		{
			ali -= 1;
			if ((*scci).active_lines.e[ali] == last_)
			{
				found = 1;
				remove_array_int(&((*scci).active_lines), ali);
				break;
			}
		}
		if (found) {}
		else
		{
			printf("Something weird happened in sc_constr_interface_remove_last_line.\n");
			exit(EXIT_FAILURE);
		}
	}
	remove_array_voidstar(&((*scci).line_data), last_, NULL);
}

void sc_constr_interface_remove_last_circle(sc_constr_interface *scci)
{
	int last_index = (*scci).circ_data.len - 1;
	circle_render_data *crd = (circle_render_data *) (*scci).circ_data.e[last_index];
	if ((*crd).vis)
	{
		int aci = (*scci).active_circles.len;
		char found = 0;
		while (aci > 0)
		{
			aci -= 1;
			if ((*scci).active_circles.e[aci] == last_index)
			{
				found = 1;
				remove_array_int(&((*scci).active_circles), aci);
				break;
			}
		}
		if (found) {}
		else
		{
			printf("Something weird happened in sc_constr_interface_remove_last_circle\n");
			exit(EXIT_FAILURE);
		}
	}
	free_circle_render_data((circle_render_data *) (*scci).circ_data.e[last_index]);
	remove_array_voidstar(&((*scci).circ_data), (*scci).circ_data.len - 1, NULL);
}

