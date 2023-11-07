#include "SDL.h"
#include "SDL_events.h"
#include "SDL_keyboard.h"
#include "euclid_vis.h"
#include "euclid.h"
#include "system.h"
#include "time.h"
#define epsilon 0.1
#define epsq 0.01
#define N_T_PTS 3
#define MAX_ROOTED_PTS 100
#define MAX_N_HOLES 100
#define DIGIT_WIDTH 35

typedef struct
{
	Uint8 *kbstate;
	SDL_Event *e;	
} scg_data;

int n_holes = 0;
int n_starting = 0;
char completed = 0;
const Uint8 *kbstate;
int n_rooted_pts = 0;
double *rooted_x;
double *rooted_y;

double hole_wid = epsilon;
double hole_widsq = epsilon * epsilon;

int SCR_LEN_X = 960;
int SCR_LEN_Y = 540;
int default_text_pos_y;
int vertex_mark_x[16] = {3, 4, 3, 2, 1, 0, -1, -2, -3, -4, -3, -2, -1, 0, 1, 2};
int vertex_mark_y[16] = {-1, 0, 1, 2, 3, 4, 3, 2, 1, 0, -1, -2, -3, -4, -3, -2};
const char *vdriver;
SDL_Texture *tex; 
double rnd();
void closest_point(int x, int y, array_int *xs, array_int *ys, int *i_);
void closest_point_excluding(int x, int y, array_int *xs, array_int *ys, int *i_, int *excl, int len);
void closest_circle(double x, double y, sc_constr_interface *scci, int *i_);
void closest_line(double x, double y, sc_constr_interface *scci, int *i_);
void render_vertex(int x, int y, SDL_Renderer *rndrr);
void render_hole(int x, int y, SDL_Renderer *rndrr);
void render_vertex_highlighted(int x, int y, SDL_Renderer *rndrr);

void test_sc_constr_points();
const int cutoff_dist_px = 100;
int cutoff_distsq;

void print_state_vars();
void reset_state_vars();
void set_cutoff_distsq();
void set_conv_factors();
void add_rooted_point(double x, double y);
void clear_hltd();


double px_wid = 0.01;
double wid_x = 0.01;
double wid_y = 0.01;
double inv_wid_x;
double inv_wid_y;
double mouse_x, mouse_y;
char mouse_reset = 1;
int i_, j_;
double _x1, _y1, _x2, _y2;
int _x1_, _y1_, _x2_, _y2_;
double zm_xbnds[2], zm_ybnds[2];
int prosp_x, prosp_y;
double lr_x[2], lr_y[2];

// Consider collecting state variables into a structure, or 
// 	integrating them implicitly in the program structure 
// 		(e.g. with distinct loops for finding points of intersection, 
// 		adding curves, etc.)
// Subprograms/loops
void welcome_loop();
void welcome_render_step();
void set_n_holes_loop();
void set_number_loop(char *imname, int *n, int n_min, int n_max);
void set_double_loop(char *imname, double *val, double _min_, double _max_);
void set_number_render_step(SDL_Texture *tex, char *digits, int len);
void main_loop();
void add_point_loop();
void select_lr_loop();
char select_curve_loop(int *i, char *select_mode);
int select_line_loop();
int select_circle_loop();
void add_curve_loop();
void add_line_loop();
void add_circle_loop();
void ctrl_loop();
void zoom_loop();
void high_score_loop();
void help_loop();
void render_step();

// State variables (to be integrated implicitly into the program flow)
char sc_init = 0;
char intersection_mode = 0;
char select_lr_flag = 0;
char lr_case;
char select_mode = 0;
int select_curve = -1;
char apm_bit = 0;
char add_point = 0;
char n_lines_point = 0;
char n_circles_point = 0;
char add_line = 0;
char add_circle = 0;
char ctrl_mode = 0;
char zoom_mode = 0;

// Straight-edge and compass construction
sc_constr sc;
// The total number of points, lines, circles, and undos used in construction
int tally[4] = {0, 0, 0, 0};
// The coordinates of target points
double *t_xs;
double *t_ys;
// Whether a point in the construction is within the prescribed distance of each target point
char *t_score;
// Total energy of a point particle configuration where target points have negative charge 
// 	and constructed points have positive charge.
double electroscore;

// Straightedge-compass interface to display
sc_constr_interface scci;
sc_constr_interface scci_t;
array_char hltd_points;
array_char hltd_lines;
array_char hltd_circles;
SDL_Window *win;
SDL_Renderer *rndrr;

double compute_electroscore();

int filter_events(void *data, SDL_Event *e)
{
	if ((*e).type == SDL_MOUSEMOTION) return 0;
	else return 1;
}

void random_rect(double *xbnds, double *ybnds, double *x, double *y);

SDL_Texture *update_SDL_texture(char *img_name, SDL_Renderer *rndrr);

double _distsq_(double _x1, double _y1, double _x2, double _y2)
{
	double delxsq = _x2 - _x1;
	double delysq = _y2 - _y1;
	delxsq *= delxsq;
	delysq *= delysq;
	return delxsq + delysq;
}

void print_t_score()
{
	int tally_total = tally[0] + tally[1] + tally[2] + tally[3] + 10;
	electroscore = compute_electroscore();
	printf("Scores: %d points, %d lines, %d circles, %d undo ops, electroscore = %g\t", tally[0], tally[1], tally[2], tally[3], electroscore);
	for (int i = 0; i < n_holes; i++)
	{
		printf("%d, ", t_score[i]);
	}
	printf("\n");
}

void exit_failure();
char query_enter(const Uint8 *kbstate)
{
	return kbstate[SDL_SCANCODE_RETURN] == 1 || kbstate[SDL_SCANCODE_RETURN2] == 1 || kbstate[SDL_SCANCODE_EQUALS] == 1;
}
char query_ctrl(const Uint8 *kbstate)
{
	return kbstate[SDL_SCANCODE_LCTRL] == 1 || kbstate[SDL_SCANCODE_RCTRL] == 1;
}

void update_t_scores(int point_addr)
{
	double epsilon_sq = epsilon * epsilon;
	completed = 1;
	for (int iii = 0; iii < n_holes; iii++)
	{
		double delsq_iii = _distsq_(scci.points_x.e[point_addr], scci.points_y.e[point_addr], t_xs[iii], t_ys[iii]);
		t_score[iii] = delsq_iii < epsilon_sq;
		completed = completed && t_score[iii];
	}
}

void compute_t_score_i(int i)
{
	if (!t_score[i])
	{
		double epsilon_sq = epsilon * epsilon;
		for (int ii = 0; ii < (*(sc.points)).len; ii++)
		{
			double distsq_i_ii = _distsq_(scci.points_x.e[ii], scci.points_y.e[ii], t_xs[i], t_ys[i]);
			if (distsq_i_ii < epsilon_sq)
			{
				t_score[i] = 1;
				break;
			}
		}
	}
}

void exit_program();

int main(int argc, char *argv[])
{
	// Initialize global variables
	default_text_pos_y = SCR_LEN_Y - 50;
	srand(time(NULL));
	array_char_init(&hltd_points, 1);
	array_char_init(&hltd_lines, 1);
	array_char_init(&hltd_circles, 1);
	
	SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS | SDL_INIT_TIMER);
	// Get the name of the graphics driver:
	int N_vdrivers = SDL_GetNumVideoDrivers();
	if (N_vdrivers > 0) vdriver = SDL_GetVideoDriver(0);
	else
	{
		exit_failure();
	}
	SDL_VideoInit(vdriver);
	win = SDL_CreateWindow("", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SCR_LEN_X, SCR_LEN_Y, SDL_WINDOW_SHOWN);
	rndrr = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
	kbstate = SDL_GetKeyboardState(NULL);
	welcome_loop();
	// exit_program(win, rndrr);
	return 0;
} // END MAIN

void closest_point(int x, int y, array_int *xs, array_int *ys, int *i_)
{
	return closest_point_excluding(x, y, xs, ys, i_, NULL, 0);
}

void closest_point_excluding(int x, int y, array_int *xs, array_int *ys, int *i_, int *excl, int len)
{
	long unsigned int min_delsq = -1;
	int j_excl = 0;
	(*i_) = -1;
	for (int i = 0; i < (*xs).len; i++)
	{
		if (j_excl < len && i != excl[j_excl]) {}
		else if (j_excl < len)
		{
			j_excl += 1;
			continue;
		}
		int delx, dely, delsq;
		delx = (*xs).e[i] - x;
		dely = (*ys).e[i] - y;
		delsq = delx * delx + dely * dely;
		if (delsq < min_delsq && delsq < cutoff_distsq)
		{
			(*i_) = i;
			min_delsq = delsq;
		}
	}
}

void closest_line(double x, double y, sc_constr_interface *scci, int *i_)
{
	double min_distsq = 1e99;
	for (int i = 0; i < (*scci).line_data.len; i++)
	{
		line * l_ = (line *) (*((*((*scci).sc)).lines)).e[i];
		double ax, ay, bx, by;
		int a_i = (*(*l_).a).addr;
		int b_i = (*(*l_).b).addr;
		ax = (*scci).points_x.e[a_i];
		bx = (*scci).points_x.e[b_i];
		ay = (*scci).points_y.e[a_i];
		by = (*scci).points_y.e[b_i];
		bx -= ax;
		by -= ay;
		ax -= x;
		ay -= y;
		// (ax + bx t, ay + by t)
		//	ax^2 + 2ax bx t + bx^2 t^2 + ay^2 + 2ay by t + by^2 t^2
		//	ax bx + bx^2 t + ay by + by^2 t = 0
		//	(ax bx + ay by) + b^2 t = 0
		//	 t = -a.b / b^2
		double bsq = bx * bx + by * by, adotb = ax * bx + ay * by, asq = ax * ax + ay * ay;
		double t = -adotb / bsq;
		double distsq = asq + (2 * adotb + bsq * t) * t;
		if (distsq < min_distsq)
		{
			(*i_) = i;
			min_distsq = distsq;
		}
	}
}

void closest_circle(double x, double y, sc_constr_interface *scci, int *i_)
{
	double min_distsq = 1e99;
	for (int i = 0; i < (*scci).circ_data.len; i++)
	{
		circle * c_ = (circle *) (*((*((*scci).sc)).circles)).e[i];
		double cx, cy, rx, ry;
		int c_i = (*(*c_).center).addr;
		int r_i = (*(*c_).radius).addr;
		cx = (*scci).points_x.e[c_i];
		cy = (*scci).points_y.e[c_i];
		rx = (*scci).points_x.e[r_i];
		ry = (*scci).points_y.e[r_i];
		rx -= cx;
		ry -= cy;
		rx = sqrt(rx * rx + ry * ry);
		double delx = x - cx, dely = y - cy, delsq, distsq;
		delsq = delx * delx + dely * dely;
		distsq = sqrt(delsq) - rx;
		distsq *= distsq;
		if (distsq < min_distsq)
		{
			(*i_) = i;
			min_distsq = distsq;
		}
	}
}

void render_vertex(int x, int y, SDL_Renderer *rndrr)
{
	for (int i = 0; i < 16; i++)
	{
		SDL_RenderDrawPoint(rndrr, x + vertex_mark_x[i], y + vertex_mark_y[i]);
	}
}

void render_vertex_highlighted(int x, int y, SDL_Renderer *rndrr)
{
	//printf("Rendering highlighted vertex at %d %d\n", x, y);
	for (int i = 0; i < 16; i++)
	{
		int x_ = x + vertex_mark_x[i];
		int y_ = y + vertex_mark_y[i];
		SDL_RenderDrawPoint(rndrr, x_, y_);
		SDL_RenderDrawPoint(rndrr, x_ + vertex_mark_x[i], y_ + vertex_mark_y[i]);
	}
	//printf("(done)\n");
}

void add_rooted_point(double x, double y)
{
	//printf("Adding rooted point at %g %g\n", x, y);
	rooted_x[n_rooted_pts] = x;
	rooted_y[n_rooted_pts] = y;
	add2array_char(&hltd_points, 0);
	point *p = (point *) calloc(1, sizeof(point));
	(*p).inc[0] = (void *) &rooted_x[n_rooted_pts];
	(*p).inc[1] = (void *) &rooted_y[n_rooted_pts];
	(*p).flag = 1;
	add2sc_constr(&sc, (void *) p, 'p');
	n_rooted_pts += 1;
	//printf("(done)\n");
}

void print_state_vars()
{
	printf("mouse: %g %g, control mode: %d (z = %d), add_point: (%d, %d, %d; %d %d %d), add_line: %d, add_circle: %d \n", mouse_x, mouse_y, ctrl_mode, zoom_mode, add_point, intersection_mode, apm_bit, select_mode, select_curve, select_lr_flag, add_line, add_circle);
}

void clear_hltd()
{
	for (int i = 0; i < hltd_points.len; i++) hltd_points.e[i] = 0;
	for (int i = 0; i < hltd_lines.len; i++) hltd_lines.e[i] = 0;
	for (int i = 0; i < hltd_circles.len; i++) hltd_circles.e[i] = 0;
}

void reset_state_vars()
{
	add_point = intersection_mode = apm_bit = select_mode = select_lr_flag = add_line = add_circle = ctrl_mode = n_lines_point = n_circles_point = 0;
	select_curve = -1;
	zoom_mode = 0;
	clear_hltd();
}

void set_cutoff_distsq()
{
	cutoff_distsq = cutoff_dist_px * cutoff_dist_px;
}

void set_conv_factors()
{
	wid_x = (scci.xbnds[1] - scci.xbnds[0]) / scci.screen_len_x;
	wid_y = (scci.ybnds[1] - scci.ybnds[0]) / scci.screen_len_y;
	px_wid = wid_x;
	inv_wid_x = scci.screen_len_x / (scci.xbnds[1] - scci.xbnds[0]);
	inv_wid_y = scci.screen_len_y / (scci.ybnds[1] - scci.ybnds[0]);
	//printf("pixel widths: %g %g\n", wid_x, wid_y);
}

void random_rect(double *xbnds, double *ybnds, double *x, double *y)
{
	(*x) = rnd() * (xbnds[1] - xbnds[0]) + xbnds[0];
	(*y) = rnd() * (ybnds[1] - ybnds[0]) + ybnds[0];
}

double rnd()
{
	return ((double) rand()) / RAND_MAX;
}

SDL_Texture *update_SDL_texture(char *img_name, SDL_Renderer *rndrr)
{
	SDL_Surface *surf = SDL_LoadBMP(img_name);
	SDL_Texture *new_tex = SDL_CreateTextureFromSurface(rndrr, surf);
	SDL_FreeSurface(surf);
	return new_tex;
}

void main_loop()
{
	// Reset variables
	reset_state_vars();
	tex = update_SDL_texture("baseline_c.bmp", rndrr);
	//const Uint8 *kbstate = SDL_GetKeyboardState(NULL);
	SDL_SetEventFilter(filter_events, NULL);
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0)
	{
		render_step();
		if (e.type == SDL_MOUSEBUTTONUP)
		{
			mouse_reset = 1;
		}
		if (e.type == SDL_QUIT)
		{
			exit_program();
		}
		if (e.type == SDL_MOUSEBUTTONDOWN && mouse_reset)
		{
			mouse_x = e.button.x * wid_x + scci.xbnds[0];
			mouse_y = e.button.y * wid_y + scci.ybnds[0];
			print_state_vars();
			print_t_score();
			mouse_reset = 0;
		}
		if (query_ctrl(kbstate))
		{
			ctrl_mode = 1;
			ctrl_loop();
			ctrl_mode = 0;
		}
		if (kbstate[SDL_SCANCODE_C] == 1)
		{
			add_circle_loop();
		}
		if (kbstate[SDL_SCANCODE_L] == 1)
		{
			add_line_loop();
		}
		if (kbstate[SDL_SCANCODE_P] == 1)
		{
			add_point_loop();
		}
		if (kbstate[SDL_SCANCODE_Q] == 1)
		{
			exit_program();
		}
	}
}

void render_step()
{
	// Try to limit background textures to a rectangle (of smallest possible area)
	SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 0);
	SDL_RenderClear(rndrr);
	int tex_w, tex_h, tex_acc;
	Uint32 tex_fmt;
	SDL_QueryTexture(tex, &tex_fmt, &tex_acc, &tex_w, &tex_h);
	SDL_Rect tex_rect;
	tex_rect.x = (SCR_LEN_X - tex_w) / 2;
	tex_rect.y = SCR_LEN_Y - tex_h;
	tex_rect.y = tex_rect.y < default_text_pos_y? tex_rect.y : default_text_pos_y;
	tex_rect.w = tex_w;
	tex_rect.h = tex_h;
	SDL_RenderCopy(rndrr, tex, NULL, &tex_rect);
	// Draw background
	SDL_SetRenderDrawColor(rndrr, 230, 230, 230, 0);
	// ...
	// Draw points, lines, circles
	// render_sc_constr(&scci, rndrr);
	// Render highlighted features
	// ...
	render_sc_constr_hlts(&scci, &hltd_points, &hltd_lines, &hltd_circles, rndrr);
	if (select_lr_flag)
	{
		render_vertex_highlighted(prosp_x, prosp_y, rndrr);
	}
	SDL_SetRenderDrawColor(rndrr, 230, 0, 230, 0);
	for (int i = 0; i < n_holes; i++)
	{
		int x_i = (int) ((t_xs[i] - scci.xbnds[0]) * inv_wid_x);
		int y_i = (int) ((t_ys[i] - scci.ybnds[0]) * inv_wid_y);
		render_vertex(x_i, y_i, rndrr);
	}
	SDL_SetRenderDrawColor(rndrr, 230, 230, 230, 0);
	SDL_RenderPresent(rndrr);
}

void add_point_loop()
{
	//printf("add_point_loop\n");
	if ((*(sc.lines)).len + (*(sc.circles)).len > 1) {}
	else
	{
		printf("Unable to add point with fewer than two curves\n");
		return;
	}
	// Update history
	// Initialize state vars (and possibly screen data)
	apm_bit = 0;
	intersection_mode = 0;
	n_lines_point = 0;
	n_circles_point = 0;
	int selection[2] = {-1, -1};
	char sel_mode[2] = {-1, -1};
	char sel_index = 0;
	while (1)
	{
		if (sel_index < 2)
		{
			char status = select_curve_loop(&selection[sel_index], &sel_mode[sel_index]);
			if (selection[sel_index] < 0)
			{
				if (sel_index) sel_index = 0;
				else main_loop();
			}
			else sel_index += 1;
		}
		else 
		{
			i_ = selection[0];
			j_ = selection[1];
			if (sel_mode[0] == 'c' || sel_mode[1] == 'c')
			{
				intersection_mode = sel_mode[0] == 'c' | ((sel_mode[1] == 'c') << 1);
				// Set lr_case
				select_lr_loop();
				if (lr_case < 0) 
				{
					sel_index -= 1;
					if (sel_mode[sel_index] == 'c') hltd_circles.e[selection[sel_index]] = 0;
					else hltd_lines.e[selection[sel_index]] = 0;
				}
				else break;
			}
			else break;
		}
	}
	if (sel_mode[0] == 'c') hltd_circles.e[i_] = 0;
	else hltd_lines.e[i_] = 0;
	if (sel_mode[1] == 'c') hltd_circles.e[j_] = 0;
	else hltd_lines.e[j_] = 0;
	intersection_mode = sel_mode[0] == 'c' | ((sel_mode[1] == 'c') << 1);
	int point_addr = (*(sc.points)).len;
	if (sel_mode[1] == 'c' || sel_mode[0] == 'c')
	{
		//select_lr_loop();
		if (sel_mode[0] == 'c')
		{
			if (sel_mode[1] == 'c') add_point_sc_constr_cc(&sc, i_, j_, lr_case);
			else add_point_sc_constr_lc(&sc, j_, i_, lr_case);
		}
		else
		{
			if (sel_mode[1] == 'c') add_point_sc_constr_lc(&sc, i_, j_, lr_case);
		}
	}
	else add_point_sc_constr_ll(&sc, i_, j_);
	add_point_sc_constr_interface(&scci, point_addr);
	add2array_char(&hltd_points, 0);
	tally[0] += 1;
	//printf("Adding point at %g %g\n", scci.points_x.e[point_addr], scci.points_y.e[point_addr]);
	update_t_scores(point_addr);
	i_ = -1;
	j_ = -1;
	main_loop();
}

void select_lr_loop()
{
	// Move the next few lines into the loop
	select_lr_flag = 1;
	lr_case = 0;
	// Determine the prospective point(s) and correspondence with lr_flag
	if (intersection_mode == 1 || intersection_mode == 2)
	{
		line *l;
		circle *c;
		if (intersection_mode == 1)
		{
			l = (line *) (*(sc.lines)).e[j_];
			c = (circle *) (*(sc.circles)).e[i_];
		}
		else
		{
			l = (line *) (*(sc.lines)).e[i_];
			c = (circle *) (*(sc.circles)).e[j_];
		}
		double cx, cy, rx, ry, ax, bx, ay, by;
		int a_ = (*(*l).a).addr, b_ = (*(*l).b).addr, c_ = (*(*c).center).addr, r_ = (*(*c).radius).addr;
		ax = scci.points_x.e[a_];
		ay = scci.points_y.e[a_];
		bx = scci.points_x.e[b_];
		by = scci.points_y.e[b_];
		cx = scci.points_x.e[c_];
		cy = scci.points_y.e[c_];
		rx = scci.points_x.e[r_];
		ry = scci.points_y.e[r_];
		line_circle_intersection_exp(ax, ay, bx, by, cx, cy, rx, ry, &lr_x[0], &lr_y[0], &lr_x[1], &lr_y[1]);
	}
	else if (intersection_mode == 3)
	{
		circle *c1 = (circle *) (*(sc.circles)).e[i_];
		circle *c2 = (circle *) (*(sc.circles)).e[j_];
		int c1i = (*(*c1).center).addr, r1i = (*(*c1).radius).addr, c2i = (*(*c2).center).addr, r2i = (*(*c2).radius).addr;
		double c1x, c1y, r1x, r1y, c2x, c2y, r2x, r2y;
		c1x = scci.points_x.e[c1i];
		c1y = scci.points_y.e[c1i];
		r1x = scci.points_x.e[r1i];
		r1y = scci.points_y.e[r1i];
		c2x = scci.points_x.e[c2i];
		c2y = scci.points_y.e[c2i];
		r2x = scci.points_x.e[r2i];
		r2y = scci.points_y.e[r2i];
		circle_circle_intersection_exp(c1x, c1y, r1x, r1y, c2x, c2y, r2x, r2y, &lr_x[0], &lr_y[0], &lr_x[1], &lr_y[1]);
	}
	prosp_x = (int) ((lr_x[lr_case] - scci.xbnds[0]) * inv_wid_x);
	prosp_y = (int) ((lr_y[lr_case] - scci.ybnds[0]) * inv_wid_y);
	tex = update_SDL_texture("addpoint_3_c.bmp", rndrr);
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0)
	{
		render_step();
		if (e.type == SDL_KEYDOWN)
		{
			if (kbstate[SDL_SCANCODE_DOWN] == 1 || kbstate[SDL_SCANCODE_UP] == 1)
			{
				lr_case = !lr_case;
				prosp_x = (int) ((lr_x[lr_case] - scci.xbnds[0]) * inv_wid_x);
				prosp_y = (int) ((lr_y[lr_case] - scci.ybnds[0]) * inv_wid_y);
			}
			if (query_enter(kbstate)) break;
			if (kbstate[SDL_SCANCODE_Q] == 1)
			{
				exit_program();
			}
			if (kbstate[SDL_SCANCODE_ESCAPE] == 1)
			{
				select_lr_flag = 0;
				lr_case = -1;
				return;
			}
		}	
		else if (e.type == SDL_QUIT)
		{
			exit_program();
		}
	}
	select_lr_flag = 0;
	return;
}

char select_curve_loop(int *i, char *select_mode)
{
	// Update display...
	tex = update_SDL_texture("addpoint_0_c.bmp", rndrr);
	// RESUME
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0)
	{
		// Update display
		render_step(); // RESUME: define this!
		// Parse input
		if (e.type == SDL_MOUSEBUTTONDOWN)
		{
			// Determine the curve closest to the mouse button
		}
		if (e.type == SDL_QUIT) exit_program();
		if (e.type == SDL_KEYDOWN)
		{
			if (kbstate[SDL_SCANCODE_L] == 1)
			{
				(*i) = select_line_loop();
				(*select_mode) = 'l';
				break;
			}
			else if (kbstate[SDL_SCANCODE_C] == 1)
			{
				(*i) = select_circle_loop();
				(*select_mode) = 'c';
				break;
			}
			if (kbstate[SDL_SCANCODE_Q] == 1) exit_program();
			if (kbstate[SDL_SCANCODE_ESCAPE] == 1) 
			{
				(*i) = -1;
				(*select_mode) = 'n';
				return -1;
			}
		}
	}
}


int select_line_loop()
{
	if ((*(sc.lines)).len > 0) {}
	else return -1;
	int li = 0;
	while (hltd_lines.e[li] == 1 && li < hltd_lines.len) li += 1;
	if (li < hltd_lines.len) {}
	else return -1;
	tex = update_SDL_texture("addpoint_1_c.bmp", rndrr);
	hltd_lines.e[li] = 1;
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0)
	{
		render_step();
		if (e.type == SDL_KEYDOWN)
		{
			if (kbstate[SDL_SCANCODE_Q] == 1) exit_program();
			if (kbstate[SDL_SCANCODE_UP] == 1)
			{
				hltd_lines.e[li] = 0;
				do
				{
					li = (li + 1) % (scci.line_data).len;
				} while (hltd_lines.e[li] == 1);
				hltd_lines.e[li] = 1;
			}
			if (kbstate[SDL_SCANCODE_DOWN] == 1)
			{
				hltd_lines.e[li] = 0;
				do
				{
					li -= 1;
					li = li > -1 ? li : (scci.line_data).len - 1;
				} while (hltd_lines.e[li] == 1);
				hltd_lines.e[li] = 1;
			}
			if (query_enter(kbstate))
			{
				return li;
			}
			if (kbstate[SDL_SCANCODE_Q] == 1) exit_program();
			if (kbstate[SDL_SCANCODE_ESCAPE] == 1) 
			{
				hltd_lines.e[li] = 0;
				return -1;
			}
		}
		else if (e.type == SDL_MOUSEBUTTONDOWN)
		{
			int x = e.button.x, y = e.button.y;
			double x_double = scci.xbnds[0] + x * wid_x;
			double y_double = scci.ybnds[0] + y * wid_y;
			hltd_lines.e[li] = 0;
			int li0 = li;
			closest_line(x_double, y_double, &scci, &li);
			if (hltd_lines.e[li] == 1) 
			{
				li = li0;
				printf("Singular intersections not supported.\n");
			}
			hltd_lines.e[li] = 1;
		}
		if (e.type == SDL_QUIT)
		{
			exit_program();
		}
	} // END WHILE
} // END SELECT LINE LOOP

int select_circle_loop()
{
	if (hltd_circles.len > 0) {}
	else return -1;
	int ci = 0;
	while (hltd_circles.e[ci] == 1 && ci < hltd_circles.len) ci += 1;
	if (ci < hltd_circles.len) {}
	else return -1;
	hltd_circles.e[ci] = 1;
	tex = update_SDL_texture("addpoint_2_c.bmp", rndrr);
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0)
	{
		render_step();
		if (e.type == SDL_KEYDOWN)
		{
			if (kbstate[SDL_SCANCODE_UP] == 1)
			{
				hltd_circles.e[ci] = 0;
				do
				{
					ci = (ci + 1) % (scci.circ_data).len;
				} while (hltd_circles.e[ci] == 1);
				hltd_circles.e[ci] = 1;
			}
			if (kbstate[SDL_SCANCODE_DOWN] == 1)
			{
				hltd_circles.e[ci] = 0;
				do
				{
					ci -= 1;
					ci = ci > -1 ? ci : (scci.circ_data).len - 1;
				} while (hltd_circles.e[ci] == 1);
				hltd_circles.e[ci] = 1;
			}
			if (query_enter(kbstate))
			{
				return ci;
			}
			if (kbstate[SDL_SCANCODE_ESCAPE] == 1)
			{
				hltd_circles.e[ci] = 0;
				return -1;
			}
			if (kbstate[SDL_SCANCODE_Q] == 1) exit_program();
		}
		else if (e.type == SDL_MOUSEBUTTONDOWN)
		{
			int x = e.button.x, y = e.button.y;
			double x_double = scci.xbnds[0] + x * wid_x;
			double y_double = scci.ybnds[0] + y * wid_y;
			hltd_circles.e[ci] = 0;
			int ci0 = ci;
			closest_circle(x_double, y_double, &scci, &ci);
			if (hltd_circles.e[ci] == 1)
			{
				ci = ci0; 
			}
			hltd_circles.e[ci] = 1;
		}
		if (e.type == SDL_QUIT) exit_program();
	}
}

void select_points_loop(int *cpi)
{
	if (hltd_points.len > 1) {}
	else 
	{
		printf("Not enough points to define a line...\n");
		cpi[0] = -1;
		cpi[1] = -1;
		return;
	}
	int exc = -1;
	cpi[0] = 0;
	cpi[1] = 0;
	char cpi_index = 0;
	array_int *xs = &(scci.points_x_int);
	array_int *ys = &(scci.points_y_int);
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0 && cpi_index < 2)
	{
		render_step();
		if (e.type == SDL_MOUSEBUTTONDOWN)
		{
			int cx = e.button.x;
			int cy = e.button.y;
			hltd_points.e[cpi[cpi_index]] = 0;
			closest_point_excluding(cx, cy, xs, ys, &cpi[cpi_index], &exc, 1);
			if (cpi[cpi_index] > -1)
			{
				hltd_points.e[cpi[cpi_index]] = 1;
				exc = cpi[cpi_index];
				cpi_index += 1;
			}
		}
		if (e.type == SDL_KEYDOWN)
		{
			if (kbstate[SDL_SCANCODE_Q] == 1)
			{
				exit_program();
			}
			if (kbstate[SDL_SCANCODE_UP] == 1)
			{
				hltd_points.e[cpi[cpi_index]] = 0;
				do
				{
					cpi[cpi_index] = (cpi[cpi_index] + 1) % (*(sc.points)).len;
				} while (hltd_points.e[cpi[cpi_index]] != 0);
				hltd_points.e[cpi[cpi_index]] = 1;
			}
			if (kbstate[SDL_SCANCODE_DOWN] == 1)
			{
				hltd_points.e[cpi[cpi_index]] = 0;
				do
				{
					cpi[cpi_index] -= 1;
					cpi[cpi_index] = cpi[cpi_index] > -1 ? cpi[cpi_index] : (*(sc.points)).len - 1;
				} while (hltd_points.e[cpi[cpi_index]] != 0);
				hltd_points.e[cpi[cpi_index]] = 1;
			}
			if (query_enter(kbstate))
			{
				exc = cpi[cpi_index];
				cpi_index += 1;
			}
			if (kbstate[SDL_SCANCODE_ESCAPE] == 1)
			{
				main_loop();
			}
		}
		if (e.type == SDL_QUIT) exit_program();
	}
	hltd_points.e[cpi[0]] = 0;
	hltd_points.e[cpi[1]] = 0;
}

void add_line_loop()
{
	if ((*(sc.points)).len > 0) {}
	else return;
	tex = update_SDL_texture("addline_c.bmp", rndrr);
	int cpi[2];
	select_points_loop(cpi);
	double x_i_ = scci.points_x.e[cpi[0]];
	double y_i_ = scci.points_y.e[cpi[0]];
	double x_j_ = scci.points_x.e[cpi[1]];
	double y_j_ = scci.points_y.e[cpi[1]];
	// printf(" from %d to %d, (%g, %g) to (%g, %g)\n", cpi[0], cpi[1], x_i_, y_i_, x_j_, y_j_);
	int line_addr = (*(sc.lines)).len;
	add_line_sc_constr_pp(&sc, cpi[0], cpi[1]);
	add_line_sc_constr_interface(&scci, line_addr);
	int x_a, y_a, x_b, y_b;
	line_render_data *lrd = (line_render_data *) scci.line_data.e[line_addr];
	x_a = (*lrd).a.x;
	y_a = (*lrd).a.y;
	x_b = (*lrd).b.x;
	y_b = (*lrd).b.y;
	// printf("Line end points: (%g %g), (%g %g), vis = %d\n", scci.xbnds[0] + x_a * wid_x, scci.ybnds[0] + y_a * wid_y, scci.xbnds[0] + x_b * wid_x, scci.ybnds[0] + y_b * wid_y, (*lrd).vis);
	add2array_char(&hltd_lines, 0);
	tally[1] += 1;
	main_loop();
}

void add_circle_loop()
{
	if ((*(sc.points)).len > 0) {}
	else return;
	tex = update_SDL_texture("addcircle_c.bmp", rndrr);
	int cpi[2];
	select_points_loop(cpi);
	double x_i_ = scci.points_x.e[cpi[0]];
	double y_i_ = scci.points_y.e[cpi[0]];
	double x_j_ = scci.points_x.e[cpi[1]];
	double y_j_ = scci.points_y.e[cpi[1]];
	// printf(" from %d to %d, (%g, %g) to (%g, %g)\n", cpi[0], cpi[1], x_i_, y_i_, x_j_, y_j_);
	int circle_addr = (*(sc.circles)).len;
	add_circle_sc_constr_pp(&sc, cpi[0], cpi[1]);
	add_circle_sc_constr_interface(&scci, circle_addr);
	add2array_char(&hltd_circles, 0);
	tally[1] += 1;
	main_loop();
}

void ctrl_loop()
{
	tex = update_SDL_texture("ctrl_mode_c.bmp", rndrr);
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0)
	{
		render_step();
		if (e.type == SDL_KEYDOWN)
		{
			if (kbstate[SDL_SCANCODE_Q] == 1)
			{
				exit_program();
			}
			if (kbstate[SDL_SCANCODE_ESCAPE] == 1)
			{
				ctrl_mode = 0;
				main_loop();
			}
			if (kbstate[SDL_SCANCODE_Z] == 1)
			{
				zoom_loop();
				zoom_mode = 1;
			}
			if (kbstate[SDL_SCANCODE_U] == 1 && sc.history.len > 0 && e.type == SDL_KEYDOWN)
			{
				tally[3] += 1;
				// Undo the last operation
				char last_op = sc.history.e[sc.history.len - 1];
				//printf("Undoing operation %c\n", last_op);
				if (last_op == 'p')
				{
					if ((*(sc.points)).len == n_rooted_pts) 
					{
						last_op = '\0';
					}
					else
					{
						// Check if the point is rooted
						point *lp = (point *) (*sc.points).e[(*(sc.points)).len - 1];
						if ((*lp).flag == 1 && n_rooted_pts > 0)
						{
							n_rooted_pts -= 1;
						}
						// Remove the last point from scci
						sc_constr_interface_remove_last_point(&scci);
					}
				}
				else if (last_op == 'c')
				{
					//
					sc_constr_interface_remove_last_circle(&scci);
				}
				else if (last_op == 'l')
				{
					//
					sc_constr_interface_remove_last_line(&scci);
				}
				if (last_op == 'p' || last_op == 'c' || last_op == 'l') sc_constr_undo(&sc);
			}
		}
		if (e.type == SDL_QUIT) exit_program();
	}
}

void zoom_loop()
{
	char pressed = 0;
	tex = update_SDL_texture("zoom_mode_c.bmp", rndrr);
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0)
	{
		render_step();
		if (e.type == SDL_MOUSEBUTTONDOWN && !pressed)
		{
			_x1_ = e.button.x;
			_y1_ = e.button.y;
			double x__ = scci.xbnds[0] + _x1_ * wid_x;
			double y__ = scci.ybnds[0] + _y1_ * wid_y;
			double hl_x = (scci.xbnds[1] - scci.xbnds[0]) * 0.25;
			double hl_y = (scci.ybnds[1] - scci.ybnds[0]) * 0.25;
			zm_xbnds[0] = x__ - hl_x;
			zm_xbnds[1] = x__ + hl_x;
			zm_ybnds[0] = y__ - hl_y;
			zm_ybnds[1] = y__ + hl_y;
			//printf("Zooming field of view to [%g, %g]x[%g, %g] centered at %g %g\n", zm_xbnds[0], zm_xbnds[1], zm_ybnds[0], zm_ybnds[1], x__, y__);
			sc_constr_interface_resize(&scci, &sc, zm_xbnds, zm_ybnds, scci.screen_len_x, scci.screen_len_y);
			set_conv_factors();
			pressed = 1;
		}
		if (e.type == SDL_MOUSEBUTTONUP)
		{
			pressed = 0;
		}
		if (e.type == SDL_KEYDOWN)
		{
			if (kbstate[SDL_SCANCODE_MINUS] && !pressed)
			{
				pressed = 1;
				double hl_x = (scci.xbnds[1] - scci.xbnds[0]) * 0.5;
				double hl_y = (scci.ybnds[1] - scci.ybnds[0]) * 0.5;
				zm_xbnds[0] = scci.xbnds[0] - hl_x;
				zm_xbnds[1] = scci.xbnds[1] + hl_x;
				zm_ybnds[0] = scci.ybnds[0] - hl_y;
				zm_ybnds[1] = scci.ybnds[1] + hl_y;
				//printf("Zooming out to [%g, %g]x[%g, %g]\n", zm_xbnds[0], zm_xbnds[1], zm_ybnds[0], zm_ybnds[1]);
				sc_constr_interface_resize(&scci, &sc, zm_xbnds, zm_ybnds, scci.screen_len_x, scci.screen_len_y);
				set_conv_factors();
			}
			else if (kbstate[SDL_SCANCODE_LCTRL] == 1 || kbstate[SDL_SCANCODE_RCTRL] == 1 || kbstate[SDL_SCANCODE_ESCAPE] == 1)
			{
				ctrl_loop();
			}
			else if (kbstate[SDL_SCANCODE_Q] == 1) exit_program();
		}
		if (e.type == SDL_KEYUP) pressed = 0;
		if (e.type == SDL_QUIT) exit_program();
	}
}

void high_score_loop();
void help_loop();

void free_globals()
{
	if (n_holes > 0)
	{
		free(t_xs);
		free(t_ys);
		free(t_score);
	}
	if (n_starting > 0)
	{
		free(rooted_x);
		free(rooted_y);
	}
	if (sc_init)
	{
		free_sc_constr(&sc);
		free_sc_constr_interface(&scci);
	}
	free_array_char(&hltd_points);
	free_array_char(&hltd_lines);
	free_array_char(&hltd_circles);
	SDL_DestroyRenderer(rndrr);
	SDL_DestroyWindow(win);
}

void exit_program()
{
	free_globals();
	SDL_VideoQuit();
	SDL_Quit();
	exit(0);
}

void exit_failure()
{
	free_globals();
	SDL_VideoQuit();
	SDL_Quit();
	exit(EXIT_FAILURE);
}

double compute_electroscore()
{
	double e_score = 0;
	for (int i = 0; i < n_holes; i++)
	{
		for (int ii = 0; ii < (*(sc.points)).len; ii++)
		{
			double delsq_i_ii = _distsq_(scci.points_x.e[ii], scci.points_y.e[ii], t_xs[i], t_ys[i]);
			e_score -= 1. / sqrt(delsq_i_ii);
		}
	}
	for (int i = 0; i < (*(sc.points)).len; i++)
	{
		for (int ii = 0; ii < i; ii++)
		{
			double delsq_i_ii = _distsq_(scci.points_x.e[ii], scci.points_y.e[ii], scci.points_x.e[i], scci.points_y.e[i]);
			e_score += 1. / sqrt(delsq_i_ii);
		}
	}
	return e_score;
}

void init_t_points()
{
	double xybnds[2] = {SCR_LEN_X * px_wid, SCR_LEN_Y * px_wid};
	double xbnds_[2] = {0, xybnds[0]};
	double ybnds_[2] = {0, xybnds[1]};
	t_xs = (double *) calloc(n_holes, sizeof(double));
	t_ys = (double *) calloc(n_holes, sizeof(double));
	t_score = (char *) calloc(n_holes, sizeof(char));
	for (int i = 0; i < n_holes; i++)
	{
		random_rect(xbnds_, ybnds_, &t_xs[i], &t_ys[i]);
		t_score[i] = 0;
	}	
}

void set_number_digits(int n, array_char *n_digits)
{
	(*n_digits).len = 0;
	while (n > 0)
	{
		char dig = n % 10;
		n /= 10;
		add2array_char(n_digits, dig);
	}
}

void decrement_digits(array_char *n_digits, int pt)
{
	if (pt > -1 && pt < (*n_digits).len)
	{
		int last_elem = (*n_digits).len - 1;
		if ((*n_digits).e[pt] > 1)
		{
			(*n_digits).e[pt] -= 1;
		}
		else if (pt < last_elem && ((*n_digits).e[pt] > 0))
		{
			(*n_digits).e[pt] = 0;
		}
		else if (pt == last_elem) // Then the last element must equal 1
		{
			do
			{
				remove_array_char(n_digits, last_elem);
				last_elem -= 1;
			}
			while ((*n_digits).e[last_elem] == 0);
		}
		else
		{
			(*n_digits).e[pt] = 9;
			decrement_digits(n_digits, pt + 1);
		}
	}
}

void increment_digits(array_char *n_digits, int pt)
{
	if (pt < (*n_digits).len)
	{
		if ((*n_digits).e[pt] < 9) (*n_digits).e[pt] += 1;
		else
		{
			(*n_digits).e[pt] = 0;
			increment_digits(n_digits, pt + 1);
		}
	}
	else
	{
		add2array_char(n_digits, 1);
	}
}

void set_number_loop(char *imname, int *n, int n_min, int n_max)
{
	(*n) = (*n) >= n_min ? (*n) : n_min;
	(*n) = (*n) <= n_max ? (*n) : n_max;
	SDL_Event e;
	char pressed = 0;
	SDL_Texture *ltex = update_SDL_texture(imname, rndrr);
	array_char n_digits;
	array_char_init(&n_digits, 2);
	set_number_digits((*n), &n_digits);
	int last_elem = n_digits.len - 1;
	while (SDL_WaitEvent(&e) >= 0)
	{
		set_number_render_step(ltex, n_digits.e, n_digits.len);
		if (e.type == SDL_KEYDOWN)
		{
			if (kbstate[SDL_SCANCODE_DOWN] == 1 && !pressed)
			{
				if ((*n) > n_min) 
				{
					(*n) -= 1;
					decrement_digits(&n_digits, 0);
				}
			}
			else if (kbstate[SDL_SCANCODE_UP] == 1 && !pressed) 
			{
				if ((*n) < n_max)
				{
					increment_digits(&n_digits, 0);
					(*n) += 1;
				}
			}
			pressed = 1;
			if (query_enter(kbstate))
			{
				free_array_char(&n_digits);
				return;
			}
			if (kbstate[SDL_SCANCODE_Q] == 1)
			{
				free_array_char(&n_digits);
				exit_program();
			}
		}
		if (e.type == SDL_KEYUP) pressed = 0;
		if (e.type == SDL_QUIT) 
		{
			free_array_char(&n_digits);
			exit_program();
		}
	}
}

void set_double_step(SDL_Texture *msg, char *digits1, int len1, char *digits2, int len2)
{
	
}

void set_number_render_step(SDL_Texture *msg, char *digits, int len)
{
	SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 0);
	SDL_RenderClear(rndrr);
	int tex_acc, tex_w, tex_h;
	Uint32 tex_fmt;
	SDL_QueryTexture(msg, &tex_fmt, &tex_acc, &tex_w, &tex_h);
	SDL_Rect tex_rect;
	tex_rect.x = (SCR_LEN_X - tex_w) / 2;
	tex_rect.y = (SCR_LEN_Y - tex_h) / 2;
	tex_rect.w = tex_w;
	tex_rect.h = tex_h;
	SDL_RenderCopy(rndrr, msg, NULL, &tex_rect);
	// Display n_holes either to the right of the message or below
	int base_x = tex_rect.x + (tex_w - len * DIGIT_WIDTH) / 2;
	int base_y = tex_rect.y + tex_h;
	for (int i = len - 1; i > -1; i--)
	{
		int digit_w, digit_h, digit_acc;
		char digit_name[256];
		sprintf(digit_name, "digit_%d_c.bmp", digits[i]);
		SDL_Texture *ltex = update_SDL_texture(digit_name, rndrr);
		Uint32 fmt;
		SDL_QueryTexture(ltex, &fmt, &digit_acc, &digit_w, &digit_h);
		SDL_Rect digit_rect;
		digit_rect.x = base_x;
		digit_rect.y = base_y;
		digit_rect.w = digit_w;
		digit_rect.h = digit_h;
		// Define rectangle for ith digit
		// Copy ith digit to renderer within specified rectangle
		SDL_RenderCopy(rndrr, ltex, NULL, &digit_rect);
		base_x += digit_w;
	}
	SDL_RenderPresent(rndrr);
}


void welcome_loop()
{
	tex = update_SDL_texture("scgolf_main_c.bmp", rndrr);
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0)
	{
		welcome_render_step();
		if (e.type == SDL_KEYDOWN || e.type == SDL_MOUSEBUTTONDOWN)
		{
			n_holes = 1;
			set_number_loop("select_n_holes_c.bmp", &n_holes, 1, MAX_N_HOLES);
			//printf("Number of holes set to %d\n", n_holes);
			n_starting = 2;
			set_number_loop("select_n_starting_points_c.bmp", &n_starting, 2, MAX_ROOTED_PTS);
			//printf("Number of starting points set to %d\n", n_starting);
			// Initialize positions of holes, rooted points, and sc structures
			sc_constr_init(&sc);
			init_t_points();
			double xybnds[2] = {SCR_LEN_X * px_wid, SCR_LEN_Y * px_wid};
			double xbnds_[2] = {0, xybnds[0]};
			double ybnds_[2] = {0, xybnds[1]};
			rooted_x = (double *) calloc(n_starting, sizeof(double));
			rooted_y = (double *) calloc(n_starting, sizeof(double));
			for (int i = 0; i < n_starting; i++)
			{
				double x_, y_;
				random_rect(xbnds_, ybnds_, &x_, &y_);
				add_rooted_point(x_, y_);
			}
			for (int i = 0; i < n_holes; i++)
			{
				random_rect(xbnds_, ybnds_, &t_xs[i], &t_ys[i]);
				t_score[i] = 0;
			}
			sc_constr_interface_init(&scci, &sc, xbnds_, ybnds_, SCR_LEN_X, SCR_LEN_Y);
			sc_init = 1;
			for (int i = 0; i < n_holes; i++)
			{
				compute_t_score_i(i);
			}
			set_cutoff_distsq();
			set_conv_factors();
			wid_x = (scci.xbnds[1] - scci.xbnds[0]) / scci.screen_len_x;
			wid_y = (scci.ybnds[1] - scci.ybnds[0]) / scci.screen_len_y;
			main_loop();
		}
	}
}

void welcome_render_step()
{
	SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 0);
	SDL_RenderClear(rndrr);
	int tex_acc, tex_w, tex_h;
	Uint32 tex_fmt;
	SDL_QueryTexture(tex, &tex_fmt, &tex_acc, &tex_w, &tex_h);
	SDL_Rect tex_rect;
	tex_rect.x = (SCR_LEN_X - tex_w) / 2;
	tex_rect.y = (SCR_LEN_Y - tex_h) / 2;
	tex_rect.w = tex_w;
	tex_rect.h = tex_h;
	SDL_RenderCopy(rndrr, tex, NULL, &tex_rect);
	SDL_RenderPresent(rndrr);
}

void save_game(char *ofprefix)
{
	mkdir_s(ofprefix);
	char ofname[256];
	sprintf(ofname, "%s/sc_constr.dat", ofprefix);
	FILE *ofile = fopen(ofname, "w");
	if (ofile != NULL)
	{
		// Print sequence of SC operations
		int ip, il, ic;
	       	ip = il = ic = 0;
		for (int i = 0; i < sc.history.len; i++)
		{
			point *a;
			point *b;
			switch (sc.history.e[i])
			{
				case 'p':
					point *p_ = (point *) (*(sc.points)).e[ip];
					if ((*p_).flag == 1)
					{
						fprintf(ofile, "p.%d.%d\n", (*p_).flag, ip);
					}
					else
					{
						int inc0, inc1;
						char mode = (*p_).flag >> 1;
						char imode = mode & 3;
						if (imode & 1)
						{
							circle *c0 = (circle *) (*p_).inc[0];
							inc0 = (*c0).addr;
						}
						else
						{
							line *l0 = (line *) (*p_).inc[0];
							inc0 = (*l0).addr;
						}
						if (imode & 2)
						{
							circle *c1 = (circle *) (*p_).inc[1];
							inc1 = (*c1).addr;
						}
						else
						{
							line *l1 = (line *) (*p_).inc[1];
							inc1 = (*l1).addr;
						}
						fprintf(ofile, "p.%d.%d.%d\n", (*p_).flag, inc0, inc1);
					}
					ip += 1;
					break;
				case 'l':
					line *l_ = (line *) (*(sc.lines)).e[il];
					a = (point *) (*l_).a;
					b = (point *) (*l_).b;
					fprintf(ofile, "l.%d.%d\n", (*a).addr, (*b).addr);
					il += 1;
					break;
				case 'c':
					circle *c_ = (circle *) (*(sc.circles)).e[ic];
					a = (point *) (*c_).center;
					b = (point *) (*c_).radius;
					fprintf(ofile, "c.%d.%d\n", (*a).addr, (*b).addr);
					ic += 1;
					break;
			}
		}
		fclose(ofile);
	}
	sprintf(ofname, "%s/rooted_pts.dat", ofprefix);
	ofile = fopen(ofname, "w");
	if (ofile != NULL)
	{
		for (int i = 0; i < n_starting; i++)
		{
			fprintf(ofile, "%g %g\n", rooted_x[i], rooted_y[i]);
		}
		fclose(ofile);
	}
	sprintf(ofname, "%s/holes.dat", ofprefix);
	ofile = fopen(ofname, "w");
	if (ofile != NULL)
	{
		for (int i = 0; i < n_holes; i++)
		{
			fprintf(ofile, "%g %g\n", t_xs[i], t_ys[i]);
		}
		fclose(ofile);
	}
	sprintf(ofname, "scores.dat");
	ofile = fopen(ofname, "a");
	if (ofile != NULL)
	{
		char found = 0, complete;
		char linebuf[256];
		int n_pts, n_lines, n_circles;
		fpos_t last_pos;
		double escore;
		fgetpos(ofile, &last_pos);
		while (fscanf(ofile, "%s.%d.%d.%d.%d.%lg", linebuf, &n_pts, &n_lines, &n_circles, &complete, &escore) != EOF)
		{
			if (strcmp(linebuf, ofprefix) == 0)
			{
				found = 1;
				break;
			}
			fgetpos(ofile, &last_pos);
		}
		fsetpos(ofile, &last_pos);
		sprintf(linebuf, "%s.%d.%d.%d.%d.%g", ofprefix, scci.points_x.len, (*(sc.lines)).len, (*(sc.circles)).len, completed, compute_electroscore());
		for (int i = 8; i < 256; i++) 
		{
			if (linebuf[i] == '\0') linebuf[i] = '@';
		}
		fprintf(ofile, "%s\n", linebuf);
		fclose(ofile);
	}
}

void test_sc_constr_points()
{
	for (int i = 0; i < (*(sc.points)).len; i++)
	{
		double x, y;
		point_coords((point *) (*(sc.points)).e[i], &x, &y);
		printf("%g %g %d %d \n", x, y, i, (*((point *) (*(sc.points)).e[i])).addr);
	}
}

