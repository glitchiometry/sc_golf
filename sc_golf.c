#include "SDL.h"
#include "SDL_events.h"
#include "SDL_keyboard.h"
#include "euclid_vis.h"
#include "euclid.h"
#include "system.h"
#include "time.h"
#define TOL_SQ 1e-14
#define epsilon 0.05
#define epsq 1.0
#define N_T_PTS 3
#define MAX_ROOTED_PTS 100
#define MAX_N_HOLES 100
#define MAX_COUNT_RELAX 1000000
#define DIGIT_WIDTH 36
#define DIGIT_HEIGHT 40

char shift_map[128];

void menu_loop();
void menu_render_step();

void init_shift_map();

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
double *rooted_x = NULL;
double *rooted_y = NULL;

double hole_wid = epsilon;
double hole_widsq = epsilon * epsilon;
double hole_wid_incr = 0.01;

aarray_char save_prog_msg;
int save_prog_msg_len;
aarray_char save_transcript_msg;
int save_transcript_msg_len;

int SCR_LEN_X = 960;
int SCR_LEN_Y = 540;
int default_text_pos_y;
int vertex_mark_x[16] = {3, 4, 3, 2, 1, 0, -1, -2, -3, -4, -3, -2, -1, 0, 1, 2};
int vertex_mark_y[16] = {-1, 0, 1, 2, 3, 4, 3, 2, 1, 0, -1, -2, -3, -4, -3, -2};
int digit_pixel_width[10];
int digit_pixel_height[10];
int digit_ascii_addr[10];
char ascii_format_offset_v[128];
char ascii_pixel_width[128];
char ascii_pixel_height[128];
int minus_pixel_width;
int minus_pixel_height;
const char *vdriver;
SDL_Texture *tex;

clock_t proc_time = 0;
int blink_time = 50000;

void enter_string_render_step(char *buf, int len, double esrs, int base_x, int base_y);
int enter_string_loop(char *buf, int *len, double scale, int base_x, int base_y);
int string_pixel_len(char *str, int len);
double point_electroscore(int point_addr);
int load_game(char *ofprefix);
void load_game_loop();
void load_game_render_step(aarray_char *opts, aarray_char *dates, int s);
void save_game(char *ofprefix, int len);
void save_game_loop();
void save_game_render_step();
void init_sentence(aarray_char *s, char *msg);

void exit_program();
void set_pixel_dimensions();
void render_image_box(SDL_Texture *img, int base_x, int base_y, int *tex_wd, int *tex_ht, double scale);
void render_string(char *str, int str_len, int base_x, int base_y, char fb, double scale);
void render_sentence(aarray_char *s, int base_x, int base_y, char fb, double scale);

double rnd();
void closest_point(int x, int y, array_int *xs, array_int *ys, int *i_);
void closest_point_excluding(int x, int y, array_int *xs, array_int *ys, int *i_, int *excl, int len);
void closest_circle(double x, double y, sc_constr_interface *scci, int *i_);
void closest_line(double x, double y, sc_constr_interface *scci, int *i_);
void render_vertex(int x, int y, SDL_Renderer *rndrr);
void render_electroscore(SDL_Renderer *rndrr, int base_x, int base_y, int prec);
void render_hole(int x, int y, SDL_Renderer *rndrr);
void render_vertex_highlighted(int x, int y, SDL_Renderer *rndrr);

void test_sc_constr_points();
const int cutoff_dist_px = 100;
int cutoff_distsq;
void resize_t_data(double *xbnds, double *ybnds);
void print_t_points();
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
void render_digit(char j, int base_x, int base_y);
void render_integer(int n, int *base_x, int *base_y, char fb);
void render_ascii(char c, int base_x, int base_y, double scale);
void render_double(double n, int n_dec, int base_x, int base_y, char fb);
void finalscore_loop();
void finalscore_render_step();
void welcome_loop();
void welcome_render_step();
int choose_list_loop(char *ttl_msg, aarray_char *opt_list, int base_x, int base_y, double scale);
void choose_list_render_step(char *ttl_msg, aarray_char *opt_list, int base_x, int base_y, double scale, int s);
void set_n_holes_loop();
void set_hole_width_loop();
void set_hole_width_render_step();
void set_integer_loop(char *imname, int *n, int n_min, int n_max, int base_x, int base_y);
void set_double_loop(char *imname, double *val, double incr, double _min_, double _max_, int base_x, int base_y);
void set_double_render_step(char *msg, char *digits, int len, int dec_pt_pos, int prec, int base_x, int base_y);

void set_integer_render_step(SDL_Texture *tex, char *digits, int len, int base_x, int base_y);
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
char continue_flag = 0;
char save_game_flag = 0;
char main_loop_init = 0;
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
double *t_xs = NULL;
double *t_ys = NULL;
circle_render_data *t_data;
// Whether a point in the construction is within the prescribed distance of each target point
char *t_score = NULL;
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
void query__exit_save(SDL_Event e);
char query_exit(SDL_Event e);

char query_enter(const Uint8 *kbstate)
{
	return kbstate[SDL_SCANCODE_RETURN] == 1 || kbstate[SDL_SCANCODE_RETURN2] == 1 || kbstate[SDL_SCANCODE_EQUALS] == 1;
}
char query_escape(const Uint8 *kbstate)
{
	return kbstate[SDL_SCANCODE_ESCAPE] == 1;
}
char query_ctrl(const Uint8 *kbstate)
{
	return kbstate[SDL_SCANCODE_LCTRL] == 1 || kbstate[SDL_SCANCODE_RCTRL] == 1;
}

void update_t_scores(int point_addr)
{
	completed = 1;
	for (int iii = 0; iii < n_holes; iii++)
	{
		double delsq_iii = _distsq_(scci.points_x.e[point_addr], scci.points_y.e[point_addr], t_xs[iii], t_ys[iii]);
		t_score[iii] = t_score[iii] || (delsq_iii < hole_widsq);
		completed = completed && t_score[iii];
	}
}

double point_electroscore(int point_addr)
{
	double lscore = 0;
	for (int i = 0; i < n_holes; i++)
	{
		double delsq = _distsq_(scci.points_x.e[point_addr], scci.points_y.e[point_addr], t_xs[i], t_ys[i]);
		lscore -= 1.0 / sqrt(delsq);
	}
	for (int i = 0; i < scci.points_x.len; i++)
	{
		if (i == point_addr) continue;
		double delsq = _distsq_(scci.points_x.e[point_addr], scci.points_y.e[point_addr], scci.points_x.e[i], scci.points_y.e[i]);
		lscore += 1.0 / sqrt(delsq);
	}
	return lscore;
}

void compute_t_score_i(int i)
{
	if (!t_score[i])
	{
		for (int ii = 0; ii < (*(sc.points)).len; ii++)
		{
			double distsq_i_ii = _distsq_(scci.points_x.e[ii], scci.points_y.e[ii], t_xs[i], t_ys[i]);
			if (distsq_i_ii < hole_widsq)
			{
				t_score[i] = 1;
				break;
			}
		}
	}
}
void init_sentence(aarray_char *s, char *msg);
void exit_program();

int main(int argc, char *argv[])
{
	// Initialize global variables
	set_cutoff_distsq();
	default_text_pos_y = SCR_LEN_Y - 50;
	srand(time(NULL));
	aarray_char_init(&save_prog_msg, 1);
	aarray_char_init(&save_transcript_msg, 1);
	char sp_msg[] = "Would you like to save your progress? (Y/N)";
	init_sentence(&save_prog_msg, sp_msg);
	save_prog_msg_len = strlen(sp_msg);
	char st_msg[] = "Would you like to save a transcript of your game? (Y/N)";
	init_sentence(&save_transcript_msg, st_msg);
	save_transcript_msg_len = strlen(st_msg);
	init_shift_map();
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
	set_pixel_dimensions();
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
		if (j_excl < len)
		{
			if (i != excl[j_excl]) {}
			else
			{
				j_excl += 1;
				continue;
			}
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

void render_target_vertex(int i, SDL_Renderer *rndrr)
{
	int x = t_data[i].center.x;
	int y = t_data[i].center.y;
	for (int i = 0; i < 16; i++)
	{
		SDL_RenderDrawPoint(rndrr, x + vertex_mark_x[i], y + vertex_mark_y[i]);
	}
	render_circle(&(t_data[i]), rndrr);
}

void render_electroscore(SDL_Renderer *rndrr, int base_x, int base_y, int prec)
{
	int p10 = 1;
	for (int i = 0; i < prec; i++) p10 *= 10;
	int approx_score = (int) (p10 * electroscore);
	SDL_Texture *ltex = update_SDL_texture("electroscore_c.bmp", rndrr);
	int ltw, lth, lt_acc;
	Uint32 ltfmt;
	SDL_QueryTexture(ltex, &ltfmt, &lt_acc, &ltw, &lth);
	SDL_Rect lt_rect;
	lt_rect.x = base_x;
	lt_rect.y = base_y;
	lt_rect.w = ltw;
	lt_rect.h = lth;
	SDL_RenderCopy(rndrr, ltex, NULL, &lt_rect);
	base_x += ltw + 15;
	render_double(electroscore, prec, base_x, base_y, 0);
	/*
	// Render first two digits
	SDL_Rect tex_rect;
	//int base_x = (19 * SCR_LEN_X) / 20;
	//int base_y = SCR_LEN_Y / 20;
	for (int i = 0; i < prec; i++)
	{
		char digit_i = approx_score % 10;
		approx_score /= 10;
		int digit_w, digit_h, digit_acc;
		char digit_name[256];
		sprintf(digit_name, "ascii_letters/ascii_%d_c.bmp", digit_ascii_addr[digit_i]);
		SDL_Texture *ltex = update_SDL_texture(digit_name, rndrr);
		Uint32 fmt;
		SDL_QueryTexture(ltex, &fmt, &digit_acc, &digit_w, &digit_h);
		SDL_Rect digit_rect;
		base_x -= digit_w;
		digit_rect.x = base_x;
		digit_rect.y = base_y;
		digit_rect.w = digit_w;
		digit_rect.h = digit_h;
		// Define rectangle for ith digit
		// Copy ith digit to renderer within specified rectangle
		SDL_RenderCopy(rndrr, ltex, NULL, &digit_rect);
	}
	// Draw decimal point
	tex_rect.x = base_x - 10;
	tex_rect.y = base_y + 40;
	tex_rect.w = 5;
	tex_rect.h = 5;
	SDL_SetRenderDrawColor(rndrr, 0, 0, 0, SDL_ALPHA_OPAQUE);
	SDL_RenderFillRect(rndrr, &tex_rect);
	SDL_SetRenderDrawColor(rndrr, 230, 230, 230, SDL_ALPHA_OPAQUE);
	// Draw digits before decimal
	do
	{
		char digit_i = approx_score % 10;
		approx_score /= 10;
		int digit_w, digit_h, digit_acc;
		char digit_name[256];
		sprintf(digit_name, "ascii_letters/ascii_%d_c.bmp", digit_ascii_addr[digit_i]);
		SDL_Texture *ltex = update_SDL_texture(digit_name, rndrr);
		Uint32 fmt;
		SDL_QueryTexture(ltex, &fmt, &digit_acc, &digit_w, &digit_h);
		SDL_Rect digit_rect;
		base_x -= digit_w;
		digit_rect.x = base_x;
		digit_rect.y = base_y;
		digit_rect.w = digit_w;
		digit_rect.h = digit_h;
		// Define rectangle for ith digit
		// Copy ith digit to renderer within specified rectangle
		SDL_RenderCopy(rndrr, ltex, NULL, &digit_rect);

	} while (approx_score != 0 && base_x > 0);
	*/
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
	printf("mouse: %g %g\n", mouse_x, mouse_y);
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
	inv_wid_x = 1. / wid_x;
	inv_wid_y = 1. / wid_y;
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
	main_loop_init = 1;
	if (completed && !continue_flag)
	{
		finalscore_loop();
	}
	// Reset variables
	reset_state_vars();
	tex = update_SDL_texture("baseline_c.bmp", rndrr);
	//const Uint8 *kbstate = SDL_GetKeyboardState(NULL);
	electroscore = compute_electroscore();
	SDL_SetEventFilter(filter_events, NULL);
	SDL_Event e;
	render_step();
	while (SDL_WaitEvent(&e) >= 0)
	{
		query__exit_save(e);
		if (e.type == SDL_MOUSEBUTTONUP)
		{
			mouse_reset = 1;
		}
		if (e.type == SDL_MOUSEBUTTONDOWN && mouse_reset)
		{
			mouse_x = e.button.x * wid_x + scci.xbnds[0];
			mouse_y = e.button.y * wid_y + scci.ybnds[0];
			print_state_vars();
			print_t_score();
			mouse_reset = 0;
		}
		if (query_escape(kbstate))
		{
			save_game_loop();
			welcome_loop();
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
	render_electroscore(rndrr, SCR_LEN_X / 10, SCR_LEN_Y / 10, 2);
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
	SDL_SetRenderDrawColor(rndrr, 230, 0, 230, SDL_ALPHA_OPAQUE);
	for (int i = 0; i < n_holes; i++)
	{

		if (t_data[i].vis == 1)
		{
			render_target_vertex(i, rndrr);
		}
		//render_vertex(x_i, y_i, rndrr);
	}
	SDL_SetRenderDrawColor(rndrr, 230, 230, 230, SDL_ALPHA_OPAQUE);
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
	// First, determine if the point is conspicuously close to an existing point
	// 	(or equals an existing point to within numerical precision)
	
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
	double px, py;
	sc_constr_point_coords_rem_exp(&sc, &(scci.points_x), &(scci.points_y), NULL, point_addr, &px, &py);
	char valid_pt = 1;
	for (int i = 0; i < scci.points_x.len; i++)
	{
		double delx = scci.points_x.e[i] - px;
		double dely = scci.points_y.e[i] - py;
		double delsq = delx * delx + dely * dely;
		if (delsq > TOL_SQ) {}
		else
		{
			valid_pt = 0;
			break;
		}
	}
	if (valid_pt)
	{
		add_point_sc_constr_interface(&scci, point_addr);
		add2array_char(&hltd_points, 0);
		tally[0] += 1;
		//printf("Adding point at %g %g\n", scci.points_x.e[point_addr], scci.points_y.e[point_addr]);
		update_t_scores(point_addr);
		electroscore += point_electroscore(point_addr);
	}
	else
	{
		sc_constr_undo(&sc);
	}
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
		query__exit_save(e);
		if (e.type == SDL_KEYDOWN)
		{
			if (kbstate[SDL_SCANCODE_DOWN] == 1 || kbstate[SDL_SCANCODE_UP] == 1)
			{
				lr_case = !lr_case;
				prosp_x = (int) ((lr_x[lr_case] - scci.xbnds[0]) * inv_wid_x);
				prosp_y = (int) ((lr_y[lr_case] - scci.ybnds[0]) * inv_wid_y);
			}
			if (query_enter(kbstate)) break;
			if (kbstate[SDL_SCANCODE_ESCAPE] == 1)
			{
				select_lr_flag = 0;
				lr_case = -1;
				return;
			}
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
		query__exit_save(e);
		// Parse input
		if (e.type == SDL_MOUSEBUTTONDOWN)
		{
			// Determine the curve closest to the mouse button
		}
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
		query__exit_save(e);
		if (e.type == SDL_KEYDOWN)
		{
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
		query__exit_save(e);
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
		query__exit_save(e);
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
	tally[2] += 1;
	main_loop();
}

void ctrl_loop()
{
	tex = update_SDL_texture("ctrl_mode_c.bmp", rndrr);
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0)
	{
		render_step();
		query__exit_save(e);
		if (e.type == SDL_KEYDOWN)
		{
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
				if (last_op == 'p' || last_op == 'c' || last_op == 'l') 
				{
					tally[3] += 1;
					sc_constr_undo(&sc);
				}
			}
		}
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
		query__exit_save(e);
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
			resize_t_data(zm_xbnds, zm_ybnds);
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
				resize_t_data(scci.xbnds, scci.ybnds);
			}
			else if (kbstate[SDL_SCANCODE_LCTRL] == 1 || kbstate[SDL_SCANCODE_RCTRL] == 1 || kbstate[SDL_SCANCODE_ESCAPE] == 1)
			{
				ctrl_loop();
			}
		}
		if (e.type == SDL_KEYUP) pressed = 0;
	}
}

void high_score_loop();
void help_loop();

void free_globals()
{
	if (t_xs != NULL)
	{
		free(t_xs);
		free(t_ys);
		free(t_score);
		for (int i = 0; i < n_holes; i++)
		{
			free_circle_render_data(&(t_data[i]));
		}
		free(t_data);
	}
	if (rooted_x != NULL)
	{
		free(rooted_x);
		free(rooted_y);
	}
	if (sc_init)
	{
		free_sc_constr(&sc);
		free_sc_constr_interface(&scci);
	}
	free_aarray_char(&save_prog_msg);
	free_aarray_char(&save_transcript_msg);
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

void resize_t_data(double *xbnds, double *ybnds)
{
	for (int i = 0; i < n_holes; i++)
	{
		free_circle_render_data(&(t_data[i]));
		circle_render_data_init_exp(&(t_data[i]), xbnds, ybnds, SCR_LEN_X, SCR_LEN_Y, t_xs[i], t_ys[i], hole_wid);
	}
}

char check_t_overlap(double x, double y, int n_set)
{
	for (int i = 0; i < n_set; i++)
	{
		double delsq = _distsq_(x, y, t_xs[i], t_ys[i]);
		if (delsq < hole_widsq) return 1;
	}
	return 0;
}

void relax_t_points()
{
	double f_x[n_holes];
	double f_y[n_holes];
	double cutoffsq = 4 * hole_widsq;
	int count = 0;
	double fsq = 1.0;
	while ((count < MAX_COUNT_RELAX) && (fsq > 1e-5))
	{
		count += 1;
		for (int i = 0; i < n_holes; i++)
		{
			f_x[i] = 0;
			f_y[i] = 0;
		}
		for (int i = 0; i < n_holes; i++)
		{
			for (int ii = 0; ii < i; ii++)
			{
				double dx = t_xs[i] - t_xs[ii], dy = t_ys[i] - t_ys[ii], dxsq;
				dxsq = dx * dx + dy * dy;
				if (dxsq < cutoffsq)
				{
					dxsq = 1.0 / dxsq;
					dx *= dxsq;
					dy *= dxsq;
					f_x[i] += dx;
					f_y[i] += dy;
					f_x[ii] -= dx;
					f_y[ii] -= dy;
				}
			}
			for (int ii = 0; ii < n_rooted_pts; ii++)
			{
				double dx = t_xs[i] - rooted_x[ii], dy = t_ys[i] - rooted_y[ii], dxsq;
				dxsq = dx * dx + dy * dy;
				if (dxsq < hole_widsq)
				{
					dxsq = 1.0 / dxsq;
					dx *= dxsq;
					dy *= dxsq;
					f_x[i] += dx;
					f_y[i] += dy;
				}
			}
		}
		fsq = 0;
		for (int i = 0; i < n_holes; i++)
		{
			t_xs[i] += f_x[i];
			t_ys[i] += f_y[i];
			fsq += f_x[i] * f_x[i] + f_y[i] * f_y[i];
		}
	}
}

void init_t_points()
{
	double xybnds[2] = {SCR_LEN_X * px_wid, SCR_LEN_Y * px_wid};
	double txbnds_[2] = {0.15 * xybnds[0], 0.85 * xybnds[0]};
	double tybnds_[2] = {0.15 * xybnds[1], 0.85 * xybnds[1]};
	double xbnds_[2] = {0, xybnds[0]};
	double ybnds_[2] = {0, xybnds[1]};
	t_xs = (double *) calloc(n_holes, sizeof(double));
	t_ys = (double *) calloc(n_holes, sizeof(double));
	t_data = (circle_render_data *) calloc(n_holes, sizeof(circle_render_data));
	t_score = (char *) calloc(n_holes, sizeof(char));
	for (int i = 0; i < n_holes; i++)
	{
		random_rect(txbnds_, tybnds_, &t_xs[i], &t_ys[i]);
	}
	relax_t_points();
	for (int i = 0; i < n_holes; i++)
	{
		circle_render_data_init_exp(&(t_data[i]), xbnds_, ybnds_, SCR_LEN_X, SCR_LEN_Y, t_xs[i], t_ys[i], hole_wid);
		t_score[i] = 0;
	}
}

void set_integer_digits(int n, array_char *n_digits)
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

void set_integer_loop(char *imname, int *n, int n_min, int n_max, int base_x, int base_y)
{
	(*n) = (*n) >= n_min ? (*n) : n_min;
	(*n) = (*n) <= n_max ? (*n) : n_max;
	SDL_Event e;
	char pressed = 0;
	SDL_Texture *ltex = update_SDL_texture(imname, rndrr);
	array_char n_digits;
	array_char_init(&n_digits, 2);
	set_integer_digits((*n), &n_digits);
	int last_elem = n_digits.len - 1;
	int ltex_wid, ltex_ht, ltex_acc;
	Uint32 ltex_fmt;
	SDL_QueryTexture(ltex, &ltex_fmt, &ltex_acc, &ltex_wid, &ltex_ht);
	/*int base_x = SCR_LEN_X - ltex_wid, base_y = SCR_LEN_Y - ltex_ht;
	base_x /= 2;
	base_y /= 2;*/
	while (SDL_WaitEvent(&e) >= 0)
	{
		set_integer_render_step(ltex, n_digits.e, n_digits.len, base_x, base_y); 
		char exit_flag = query_exit(e);
		if (exit_flag)
		{
			free_array_char(&n_digits);
			if (main_loop_init) save_game_loop();
			exit_program();
		}
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
		}
		if (e.type == SDL_KEYUP) pressed = 0;
	}
}

void set_double_render_step(char *msg, char *digits, int len, int dec_pt_pos, int prec, int base_x, int base_y)
{
	SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 0);
	SDL_RenderClear(rndrr);
	render_string(msg, strlen(msg), base_x, base_y, 0, 1.0);
	int tex_w, tex_h;
	tex_w = string_pixel_len(msg, strlen(msg));
	tex_h = DIGIT_HEIGHT;
	base_x += tex_w + 20;
	int i = len;
	int prec_lim = dec_pt_pos - prec - 1;
	if (prec_lim < 0)
	{
		printf("Precision warning: not enough digits stored past decimal point for desired precision\n");
	}
	if (len <= dec_pt_pos)
	{
		i = dec_pt_pos;
		render_ascii('0', base_x, base_y, 1.0);
		base_x += ascii_pixel_width['0'];
		render_ascii('.', base_x, base_y + ascii_pixel_height['4'], 1.0);
		base_x += ascii_pixel_width['.'];
		while (i > len)
		{
			i -= 1;
			render_ascii('0', base_x, base_y, 1.0);
			base_x += ascii_pixel_width['0'];
		}
	}
	else
	{
		i = len;
	}
	int penultimate = prec_lim + 1;
	i -= 1;
	while (i > penultimate)
	{
		render_digit(digits[i], base_x, base_y);
		base_x += ascii_pixel_width[digit_ascii_addr[digits[i]]];
		if (i == dec_pt_pos)
		{
			render_ascii('.', base_x, base_y + ascii_pixel_height['4'], 1.0);
			base_x += ascii_pixel_width['.'];
		}
		i -= 1;
	}
	if (digits[prec_lim] > 4)
	{
		render_digit(digits[penultimate] + 1, base_x, base_y);
	}
	else
	{
		render_digit(digits[penultimate], base_x, base_y);
	}
	SDL_RenderPresent(rndrr);
}

void render_digit(char j, int base_x, int base_y)
{
	char imname[256];
	int i_j = digit_ascii_addr[j];
	sprintf(imname, "ascii_letters/ascii_%d_c.bmp", i_j);
	SDL_Texture *lt = update_SDL_texture(imname, rndrr);
	SDL_Rect lt_rect;
	lt_rect.x = base_x;
	lt_rect.y = base_y;
	lt_rect.w = ascii_pixel_width[i_j];
	lt_rect.h = ascii_pixel_height[i_j];
	SDL_RenderCopy(rndrr, lt, NULL, &lt_rect);
}

void render_double(double n, int n_dec, int base_x, int base_y, char fb)
{
	int p10 = 1;
	for (int i = 0; i < n_dec; i++) p10 *= 10;
	double n_ = n >= 0 ? n : -n;
	int n0 = (int) n_;
	int frac_part = (int) ((n_ - n0) * p10 * 10);
	int last_digit = frac_part % 10;
	frac_part /= 10;
	frac_part = last_digit < 5 ? frac_part : frac_part + 1;
	if (!fb) 
	{
		if (n < 0)
		{
			render_ascii('-', base_x, base_y + DIGIT_HEIGHT / 2, 1.0);
			base_x += ascii_pixel_width['-'];
		}
		render_integer(n0, &base_x, &base_y, fb);
		render_ascii('.', base_x, base_y + ascii_pixel_height['4'], 1.0);
		base_x += ascii_pixel_width['.'];
		while (frac_part < p10)
		{
			p10 /= 10;
			render_digit(0, base_x, base_y);
			base_x += ascii_pixel_width['0'];
		}
		render_integer(frac_part, &base_x, &base_y, fb);
	}
	else
	{
		render_integer(frac_part, &base_x, &base_y, fb);
		while (frac_part < p10)
		{
			p10 /= 10;
			base_x -= ascii_pixel_width['0'];
			render_digit(0, base_x, base_y);
		}
		render_ascii('.', base_x, base_y + DIGIT_HEIGHT, 1.0);
		base_x -= ascii_pixel_width['.'];
		render_integer(n0, &base_x, &base_y, fb);
	}

}

void render_integer(int n, int *base_x, int *base_y, char fb)
{
	int len = 0, n_ = n >= 0 ? n : -n;
	char digits[64];
	int total_len = 0;
	if (n_ == 0)
	{
		len = 1;
		digits[0] = 0;
	}
	if (!fb)
	{
		while (n_ > 0)
		{
			digits[len] = n_ % 10;
			int incr = ascii_pixel_width[digit_ascii_addr[digits[len]]];
			total_len += incr;
			n_ /= 10;
			len += 1;
		}
		if (n < 0)
		{
			render_ascii('-', (*base_x) - ascii_pixel_width['-'], (*base_y), 1.0);
		}
		do
		{
			len -= 1;
			render_digit(digits[len], (*base_x), (*base_y));
			int incr = ascii_pixel_width[digit_ascii_addr[digits[len]]];
			(*base_x) += incr;
		} while (len > 0);
	}
	else
	{
		while (n_ > 0)
		{
			char di = n_ % 10;
			n_ /= 10;
			(*base_x) -= ascii_pixel_width[digit_ascii_addr[di]];
			render_digit(di, (*base_x), (*base_y));
		}
		if (n < 0)
		{
			render_ascii('-', (*base_x) - ascii_pixel_width['-'], (*base_y), 1.0);
		}
	}
}

void set_integer_render_step(SDL_Texture *msg, char *digits, int len, int base_x, int base_y)
{
	SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 0);
	SDL_RenderClear(rndrr);
	int tex_acc, tex_w, tex_h;
	Uint32 tex_fmt;
	SDL_QueryTexture(msg, &tex_fmt, &tex_acc, &tex_w, &tex_h);
	SDL_Rect tex_rect;
	tex_rect.x = base_x;
	tex_rect.y = base_y;
	//tex_rect.x = (SCR_LEN_X - tex_w) / 2;
	//tex_rect.y = (SCR_LEN_Y - tex_h) / 2;
	tex_rect.w = tex_w;
	tex_rect.h = tex_h;
	SDL_RenderCopy(rndrr, msg, NULL, &tex_rect);
	// Display n_holes either to the right of the message or below
	base_x = tex_rect.x + (tex_w - len * DIGIT_WIDTH) / 2;
	base_y = tex_rect.y + tex_h;
	int i = len;
	do
	{
		i -= 1;
		int digit_w, digit_h, digit_acc;
		char digit_name[256];
		sprintf(digit_name, "ascii_letters/ascii_%d_c.bmp", digit_ascii_addr[digits[i]]);
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
	} while (i > 0);
	SDL_RenderPresent(rndrr);
}


void finalscore_loop()
{
	// Prepare the SDL_Texture
	SDL_Event e;
	finalscore_render_step();
	while (SDL_WaitEvent(&e) >= 0)
	{
		//finalscore_render_step();
		if (e.type == SDL_KEYDOWN)
		{
			if (kbstate[SDL_SCANCODE_C] == 1)
			{
				continue_flag = 1;
				main_loop();
			}
			if (kbstate[SDL_SCANCODE_M] == 1)
			{
				save_game_loop();
				welcome_loop();
			}
		}
		query__exit_save(e);
	}
}

void render_image_box(SDL_Texture *img, int base_x, int base_y, int *tex_wd, int *tex_ht, double scale)
{
	int tex_acc;
	Uint32 tex_fmt;
	SDL_QueryTexture(img, &tex_fmt, &tex_acc, tex_wd, tex_ht);
	SDL_Rect tex_window;
	tex_window.x = base_x;
	tex_window.y = base_y;
	tex_window.w = (int) (scale * (*tex_wd));
	tex_window.h = (int) (scale * (*tex_ht));
	SDL_RenderCopy(rndrr, img, NULL, &tex_window);
}

// RESUME: Check if a new rectangle is needed for each field/image
void finalscore_render_step()
{
	SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 0);
	SDL_RenderClear(rndrr);
	SDL_SetRenderDrawColor(rndrr, 230, 230, 230, 255);
	int tex_wd, tex_ht, tex_acc;
	Uint32 tex_fmt;

	//SDL_Rect tex_window;
	SDL_Texture *title_msg = update_SDL_texture("score_c.bmp", rndrr);
	SDL_QueryTexture(title_msg, &tex_fmt, &tex_acc, &tex_wd, &tex_ht);
	int base_x = (SCR_LEN_X - tex_wd) / 2;
	int base_y = tex_ht / 2;
	render_image_box(title_msg, base_x, base_y, &tex_wd, &tex_ht, 1.0);
	// Render options
	char continue_msg[] = "Continue (C)";
	int cmsg_len = string_pixel_len(continue_msg, strlen(continue_msg));
	base_x = (SCR_LEN_X) - cmsg_len;
	base_y = SCR_LEN_Y - DIGIT_HEIGHT - 15;
	render_string(continue_msg, strlen(continue_msg), base_x, base_y, 0, 0.7);
	base_x = 10;
	char exit_msg[] = "Quit (Q)";
	render_string(exit_msg, strlen(exit_msg), base_x, base_y, 0, 0.7);
	char newgame_msg[] = "Menu (M)";
	base_x = (SCR_LEN_X - strlen(newgame_msg) * DIGIT_WIDTH) / 2;
	render_string(newgame_msg, strlen(newgame_msg), base_x, base_y, 0, 0.7);
	/*
	SDL_QueryTexture(title_msg, &tex_fmt, &tex_acc, &tex_wd, &tex_ht);
	tex_window.x = (SCR_LEN_X - tex_wd) / 2;
	tex_window.y = tex_ht / 2;
	tex_window.w = tex_wd;
	tex_window.h = tex_ht;
	SDL_RenderCopy(rndrr, title_msg, NULL, &tex_window);
	*/
	SDL_Texture *points_img = update_SDL_texture("point_count_c.bmp", rndrr);
	base_y = tex_ht << 1;
	base_x = SCR_LEN_X / 20;
       	render_image_box(points_img, base_x, base_y, &tex_wd, &tex_ht, 1.0);	
	/*SDL_QueryTexture(score_name, &tex_fmt, &tex_acc, &tex_wd, &tex_ht);
	tex_window.x = base_x;
	tex_window.y = base_y;
	tex_window.w = tex_wd;
	tex_window.h = tex_ht;
	SDL_RenderCopy(rndrr, score_name, NULL, &tex_window);
	*/
	// Render the associated number
	int aux_bx = base_x + tex_wd + 10, aux_by = base_y;
	render_integer((*(sc.points)).len, &aux_bx, &aux_by, 0);
	base_y += (tex_ht * 7) / 5;
	SDL_Texture *lines_img = update_SDL_texture("line_count_c.bmp", rndrr);
	render_image_box(lines_img, base_x, base_y, &tex_wd, &tex_ht, 1.0);
	aux_bx = base_x + tex_wd + 10;
	aux_by = base_y;
	render_integer((*(sc.lines)).len, &aux_bx, &aux_by, 0);
	/*SDL_QueryTexture(score_name, &tex_fmt, &tex_acc, &tex_wd, &tex_ht);
	SDL_Rect tex_rect;
	tex_rect.x = (SCR_LEN_X - tex_wd) / 2;
	tex_rect.y = (SCR_LEN_Y - tex_ht) / 2;
	tex_rect.w = tex_wd;
	tex_rect.h = tex_ht;
	SDL_RenderCopy(rndrr, tex, NULL, &tex_rect);
	aux_bx = base_x + (7 * tex_wd) / 5;
	aux_by = base_y;
	render_integer((*(sc.lines)).len, &aux_bx, &aux_by, 0);*/
	base_y += (tex_ht * 7) / 5;
	SDL_Texture *circles_img = update_SDL_texture("circle_count_c.bmp", rndrr);
	render_image_box(circles_img, base_x, base_y, &tex_wd, &tex_ht, 1.0);
	/*SDL_QueryTexture(score_name, &tex_fmt, &tex_acc, &tex_wd, &tex_ht);
	tex_window.x = base_x;
	tex_window.y = base_y;
	tex_window.w = tex_wd;
	tex_window.h = tex_ht;
	SDL_RenderCopy(rndrr, score_name, NULL, &tex_window);*/
	aux_bx = base_x + tex_wd + 10;
	aux_by = base_y;
	render_integer((*(sc.circles)).len, &aux_bx, &aux_by, 0);
	base_y += (tex_ht * 7) / 5;
	SDL_Texture *escore_img = update_SDL_texture("electroscore_c.bmp", rndrr);
	render_image_box(escore_img, base_x, base_y, &tex_wd, &tex_ht, 1.0);
	/*SDL_QueryTexture(score_name, &tex_fmt, &tex_acc, &tex_wd, &tex_ht);
	tex_window.x = base_x;
	tex_window.y = base_y;
	tex_window.w = tex_wd;
	tex_window.h = tex_ht;
	SDL_RenderCopy(rndrr, score_name, NULL, &tex_window);*/
	aux_bx = base_x + tex_wd + 10;
	aux_by = base_y;
	render_double(electroscore, 2, aux_bx, aux_by, 0);
	SDL_RenderPresent(rndrr);
}

void welcome_loop()
{
	n_holes = 0;
	n_starting = 0;
	completed = 0;
	n_rooted_pts = 0;
	main_loop_init = 0;
	if (rooted_x != NULL)
	{
		free(rooted_x);
		free(rooted_y);
		rooted_x = NULL;
		rooted_y = NULL;
	}
	if (t_xs != NULL)
	{
		free(t_xs);
		free(t_ys);
		free(t_data);
		free(t_score);
		t_xs = NULL;
		t_ys = NULL;
		t_data = NULL;
		t_score = NULL;
	}
	if (sc_init)
	{
		sc_init = 0;
		free_sc_constr(&sc);
		free_sc_constr_interface(&scci);
	}
	hltd_points.len = 0;
	hltd_lines.len = 0;
	hltd_circles.len = 0;
	hole_wid = epsilon;
	hole_widsq = epsilon * epsilon;
	hole_wid_incr = 0.01;
	tex = update_SDL_texture("scgolf_main_c.bmp", rndrr);
	welcome_render_step();
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0)
	{
		//welcome_render_step();
		if (query_exit(e))
		{
			exit_program();
		}
		if (e.type == SDL_KEYDOWN || e.type == SDL_MOUSEBUTTONDOWN)
		{
			if (kbstate[SDL_SCANCODE_N] == 1)
			{
				n_holes = 1;
				SDL_Texture *ltex = update_SDL_texture("select_n_holes_c.bmp", rndrr);
				int ltex_wid, ltex_ht, ltex_acc;
				Uint32 ltex_fmt;
				SDL_QueryTexture(ltex, &ltex_fmt, &ltex_acc, &ltex_wid, &ltex_ht);
				int base_x = SCR_LEN_X - ltex_wid, base_y = SCR_LEN_Y - ltex_ht;
				base_x /= 2;
				base_y /= 2;
				set_integer_loop("select_n_holes_c.bmp", &n_holes, 1, MAX_N_HOLES, base_x, base_y);
				//printf("Number of holes set to %d\n", n_holes);
				n_starting = 2;
				ltex = update_SDL_texture("select_n_starting_points_c.bmp", rndrr);
				SDL_QueryTexture(ltex, &ltex_fmt, &ltex_acc, &ltex_wid, &ltex_ht);
				base_x = (SCR_LEN_X - ltex_wid) / 2;
				base_y = (SCR_LEN_Y - ltex_ht) / 2;
				set_integer_loop("select_n_starting_points_c.bmp", &n_starting, 2, MAX_ROOTED_PTS, base_x, base_y);
				//printf("Number of starting points set to %d\n", n_starting);
				// Initialize positions of holes, rooted points, and sc structures
				sc_constr_init(&sc);
				double xybnds[2] = {SCR_LEN_X * px_wid, SCR_LEN_Y * px_wid};
				double xbnds_[2] = {0, xybnds[0]};
				double ybnds_[2] = {0, xybnds[1]};
				double rxbnds[2] = {0.15 * xybnds[0], 0.85 * xybnds[0]};
				double rybnds[2] = {0.15 * xybnds[1], 0.85 * xybnds[1]};
				rooted_x = (double *) calloc(n_starting, sizeof(double));
				rooted_y = (double *) calloc(n_starting, sizeof(double));
				for (int i = 0; i < n_starting; i++)
				{
					double x_, y_;
					random_rect(rxbnds, rybnds, &x_, &y_);
					add_rooted_point(x_, y_);
				}
				
				sc_constr_interface_init(&scci, &sc, xbnds_, ybnds_, SCR_LEN_X, SCR_LEN_Y);
				sc_init = 1;
				set_hole_width_loop();
				init_t_points();
				for (int i = 0; i < n_holes; i++)
				{
					t_score[i] = 0;
				}
				for (int i = 0; i < n_holes; i++)
				{
					compute_t_score_i(i);
				}
				//set_cutoff_distsq();
				set_conv_factors();
				main_loop();
			}
			else if (kbstate[SDL_SCANCODE_L] == 1)
			{
				load_game_loop();
			}
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
	tex_rect.y = SCR_LEN_Y / 20;
	tex_rect.w = tex_w;
	tex_rect.h = tex_h;
	SDL_RenderCopy(rndrr, tex, NULL, &tex_rect);
	int base_y = tex_rect.y + tex_h + 20;
	render_string("Load game (L)", 13, SCR_LEN_X / 5, base_y, 0, 1.0);
	render_string("New game (N)", 12, SCR_LEN_X / 5, base_y + DIGIT_HEIGHT + 30, 0, 1.0);
	render_string("Quit (Q)", 8, SCR_LEN_X / 5, base_y + 2 * DIGIT_HEIGHT + 60, 0, 1.0);
	SDL_RenderPresent(rndrr);
}

char valid_alphanumeric(char c)
{
	return c == 45 || c == 43 || c == 46 || (47 < c && c < 58) || (64 < c && c < 91) || c == 95 || (96 < c && c < 123);
}

void enter_string_render_step(char *buf, int len, double esrs, int base_x, int base_y)
{
	SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 0);
	SDL_RenderClear(rndrr);
	SDL_SetRenderDrawColor(rndrr, 230, 230, 230, SDL_ALPHA_OPAQUE);
	render_string(buf, len, base_x, base_y, 0, esrs);
	if (proc_time < blink_time)
	{
		int str_wid = (int) (string_pixel_len(buf, len) * esrs);
		int height = (int) (DIGIT_HEIGHT * esrs);
		render_ascii('_', base_x + str_wid, base_y + height, esrs);
	}
}

int enter_string_loop(char *buf, int *len, double scale, int base_x, int base_y)
{
	enter_string_render_step(buf, (*len), scale, base_x, base_y);
	SDL_RenderPresent(rndrr);
	char pressed = 0;
	SDL_Event e;
	while (1)
	{
		proc_time += 1;
		proc_time %= (blink_time << 1);
		if (proc_time == 0 || proc_time == blink_time)
		{
			enter_string_render_step(buf, (*len), scale, base_x, base_y);
			SDL_RenderPresent(rndrr);
		}
		if (SDL_PollEvent(&e) >= 0)
		{
			
			if (query_escape(kbstate))
			{
				return -1;
			}
			if (e.type == SDL_QUIT)
			{
				exit_program();
			}
			if (e.type == SDL_KEYUP)
			{
				pressed = 0;
			}
			if (e.type == SDL_KEYDOWN)
			{
				if (query_enter(kbstate))
				{
					return 0;
				}
				// Check if key is a valid alpha-numeric character
				char c = e.key.keysym.sym;
				if (valid_alphanumeric(c) && !pressed)
				{
					if (kbstate[SDL_SCANCODE_RSHIFT] == 1 || kbstate[SDL_SCANCODE_LSHIFT] == 1)
					{
						printf("%c -> %c\n", c, shift_map[c]);
						c = shift_map[c] != -1 ? shift_map[c] : c;
					}
					pressed = 1;
					buf[(*len)] = c;
					(*len) += 1;
					enter_string_render_step(buf, (*len), scale, base_x, base_y);
					SDL_RenderPresent(rndrr);
				}
				if (kbstate[SDL_SCANCODE_BACKSPACE] == 1 && (*len) > 0 && ! pressed)
				{
					pressed = 1;
					(*len) -= 1;
					buf[(*len)] = '\0';
					enter_string_render_step(buf, (*len), scale, base_x, base_y);
					SDL_RenderPresent(rndrr);
				}
			}
		}
	}
}

void save_game_render_step()
{
	SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 0);
	SDL_RenderClear(rndrr);
	if (!completed) 
	{
		render_sentence(&save_prog_msg, (SCR_LEN_X - save_prog_msg_len * DIGIT_WIDTH) / 2, (SCR_LEN_Y * 3) / 20, 0, 0.8);
	}
	else 
	{
		render_sentence(&save_transcript_msg, (SCR_LEN_X - save_transcript_msg_len * DIGIT_WIDTH) / 2, (SCR_LEN_Y * 3) / 20, 0, 0.8);
	}
	SDL_RenderPresent(rndrr);
}

void save_game_loop()
{
	save_game_render_step();
	save_game_flag = 1;
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0)
	{
		if (e.type == SDL_KEYDOWN)
		{
			if (query_exit(e))
			{
				exit_program();
			}
			if (kbstate[SDL_SCANCODE_Y] == 1)
			{
				int buf_len = 0;
				char buf[256];
				int status = enter_string_loop(buf, &buf_len, 0.8, ( 3 * SCR_LEN_X) / 10, (3 * SCR_LEN_Y) / 10); 
				if ((status > -1) && (buf_len > 0)) save_game(buf, buf_len);
				save_game_flag = 0;
				return;
			}
			if (query_escape(kbstate) || kbstate[SDL_SCANCODE_N] == 1)
			{
				save_game_flag = 0;
				return;
			}
		}
	}
}

void delete_saved_game(char *ofprefix)
{
	// Update the scores file
	// Remove the associated folder
}

void save_game(char *ofprefix, int len)
{
	char aux_prefix[256];
	for (int i = 0; i < len; i++) aux_prefix[i] = ofprefix[i];
	for (int i = len; i < 256; i++) aux_prefix[i] = '\0';
	mkdir_s(aux_prefix);
	char ofname[512];
	time_t mdy = time(NULL);
	char *mdyhr = ctime(&mdy);
	printf("Save game: %s %d\n", aux_prefix, len);
	sprintf(ofname, "%s/sc_constr.dat", aux_prefix);
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
	sprintf(ofname, "%s/rooted_pts.dat", aux_prefix);
	ofile = fopen(ofname, "w");
	if (ofile != NULL)
	{
		fprintf(ofile, "%d\n", n_rooted_pts);
		for (int i = 0; i < n_starting; i++)
		{
			fprintf(ofile, "%g %g\n", rooted_x[i], rooted_y[i]);
		}
		fclose(ofile);
	}
	sprintf(ofname, "%s/holes.dat", aux_prefix);
	ofile = fopen(ofname, "w");
	if (ofile != NULL)
	{
		fprintf(ofile, "%d %g\n", n_holes, hole_wid);
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
		char linebuf[512];
		int n_pts, n_lines, n_circles, aux_complete;
		fpos_t last_pos;
		double escore;
		time_t tdata;
		fgetpos(ofile, &last_pos);
		while (fscanf(ofile, "%s.%d.%d.%d.%d.%lg.%ld.", linebuf, &n_pts, &n_lines, &n_circles, &aux_complete, &escore, &tdata) != EOF)
		{
			if (strcmp(linebuf, aux_prefix) == 0)
			{
				found = 1;
				break;
			}
			fgetpos(ofile, &last_pos);
		}
		fsetpos(ofile, &last_pos);
		sprintf(linebuf, "%s.%d.%d.%d.%d.%g.%ld.", aux_prefix, scci.points_x.len, (*(sc.lines)).len, (*(sc.circles)).len, completed, compute_electroscore(), mdy);
		int len = strlen(linebuf);
		for (int i = len; i < 256; i++) 
		{
			linebuf[i] = '@';
		}
		for (int i = 256; i < 512; i++) linebuf[i] = '\0';
		fprintf(ofile, "%s\n", linebuf);
		fclose(ofile);
	}
}

int load_game(char *fprefix)
{
	if (sc_init) {}
	else
	{
		sc_constr_init(&sc);
	}
	char buf[256];
	sprintf(buf, "%s/rooted_pts.dat", fprefix);
	FILE *ifile = fopen(buf, "r");
	if (ifile != NULL)
	{
		if (fscanf(ifile, "%d\n", &n_starting) == 1)
		{
			rooted_x = (double *) calloc(n_starting, sizeof(double));
			rooted_y = (double *) calloc(n_starting, sizeof(double));
			n_rooted_pts = 0;
			double x, y;
			while (fscanf(ifile, "%lg %lg\n", &x, &y) != EOF)
			{
				add_rooted_point(x, y);
			}
		}
		fclose(ifile);
	}
	double xbnds_[2] = {0, SCR_LEN_X * px_wid};
	double ybnds_[2] = {0, SCR_LEN_Y * px_wid};
	sc_constr_interface_init(&scci, &sc, xbnds_, ybnds_, SCR_LEN_X, SCR_LEN_Y);
	sc_init = 1;
	sprintf(buf, "%s/holes.dat", fprefix);
	ifile = fopen(buf, "r");
	if (ifile != NULL)
	{
		if (fscanf(ifile, "%d %lg", &n_holes, &hole_wid) == 2)
		{
			hole_widsq = hole_wid * hole_wid;
			t_xs = (double *) calloc(n_holes, sizeof(double));
			t_ys = (double *) calloc(n_holes, sizeof(double));
			t_data = (circle_render_data *) calloc(n_holes, sizeof(circle_render_data));
			t_score = (char *) calloc(n_holes, sizeof(char));
			for (int i = 0; i < n_holes; i++)
			{
				int status = fscanf(ifile, "%lg %lg", &(t_xs[i]), &(t_ys[i]));
				if (status == 2) {}
				else break;
			}
			for (int i = 0; i < n_holes; i++)
			{
				circle_render_data_init_exp(&(t_data[i]), xbnds_, ybnds_, SCR_LEN_X, SCR_LEN_Y, t_xs[i], t_ys[i], hole_wid);
				t_score[i] = 0;
			}
		}
		fclose(ifile);
	}
	sprintf(buf, "%s/sc_constr.dat", fprefix);
	ifile = fopen(buf, "r");
	if (ifile != NULL)
	{
		char line_buf[256];
		char mode;
		int inc1, inc2, inc3, n_args;
		//while ((n_args = fscanf(ifile, "%c.%d.%d.%d\n", &mode, &inc1, &inc2, &inc3)) != EOF)
		while (fscanf(ifile, "%s\n", line_buf) != EOF)
		{
			n_args = sscanf(line_buf, "%c.%d.%d.%d", &mode, &inc1, &inc2, &inc3);
			switch (mode)
			{
				case 'p':
					int pflag = inc1;
					if (n_args == 3)
					{
						if (inc1 == 1) {}
						else
						{
							printf("Something weird happened!\n");
							exit(EXIT_FAILURE);
						}
					}
					if (n_args == 4)
					{
						int point_addr = (*(sc.points)).len;
						if (pflag == 0)
						{
							add_point_sc_constr_ll(&sc, inc2, inc3);
						}
						else
						{
							char lr = (pflag >> 3) & 1;
							char imode = pflag >> 1;
							if (imode & 3 == 1)
							{
								add_point_sc_constr_lc(&sc, inc3, inc2, lr);
							}
							if ((imode & 3) == 2)
							{
								add_point_sc_constr_lc(&sc, inc2, inc3, lr);
							}
							if ((imode & 3) == 3)
							{
								add_point_sc_constr_cc(&sc, inc2, inc3, lr);
							}
						}
						add_point_sc_constr_interface(&scci, point_addr);
						double px, py;
						point_coords((point *) (*(sc.points)).e[point_addr], &px, &py);

					}
					add2array_char(&hltd_points, 0);
					break;
				case 'l':
					int line_addr = (*(sc.lines)).len;
					add_line_sc_constr_pp(&sc, inc1, inc2);
					add_line_sc_constr_interface(&scci, line_addr);
					add2array_char(&hltd_lines, 0);
					break;
				case 'c':
					int circ_addr = (*(sc.circles)).len;
					double cx, cy, rx, ry;
					point_coords((point *) (*(sc.points)).e[inc1], &cx, &cy);
					point_coords((point *) (*(sc.points)).e[inc2], &rx, &ry);
					add_circle_sc_constr_pp(&sc, inc1, inc2);
					add_circle_sc_constr_interface(&scci, circ_addr);
					add2array_char(&hltd_circles, 0);
					// RESUME: make sure that hlted arrays are updated as well,
					// 	and that scci is initialized.
					break;
			}
		}
		fclose(ifile);
	}
	return 0;
}

void load_game_loop()
{
	aarray_char save_games;
	aarray_char dates;
	aarray_char_init(&dates, 1);
	aarray_char_init(&save_games, 1);
	FILE *ifile = fopen("scores.dat", "r");
	if (ifile != NULL)
	{
		char line_buf[512];
		while (fscanf(ifile, "%s", line_buf) != EOF)
		{
			int i_ = save_games.len;
			extend_aarray_char(&save_games);
			int i = 0;
			while (line_buf[i] != '.')
			{
				add2array_char(&(save_games.e[i_]), line_buf[i]);
				i += 1;
			}
			int i0 = 0;
			while (line_buf[i] != '@' && line_buf[i] != '\0')
			{
				line_buf[i0] = line_buf[i];
			       	i += 1;
				i0 += 1;
			}
			char sgname[256];
			int n_pts, n_lines, n_circles, aux_complete;
		       	double escore;
			time_t tdata;
			sscanf(line_buf, ".%d.%d.%d.%d.%lg.%ld.", &n_pts, &n_lines, &n_circles, &aux_complete, &escore, &tdata);
			char *mdy_ = ctime(&tdata);
			extend_aarray_char(&dates);
			for (int ii = 0; ii < strlen(mdy_); ii++)
			{
				add2array_char(&(dates.e[i_]), mdy_[ii]);
			}
		}
		fclose(ifile);
	}
	// int selection = 0;
	int base_x = SCR_LEN_X / 10;
	int base_y = SCR_LEN_Y / 15;
	int selection = 0;
	load_game_render_step(&save_games, &dates, selection);
	int pressed = 0;
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0)
	{
		if (query_exit(e))
		{
			free_aarray_char(&save_games);
			free_aarray_char(&dates);
			exit_program();
		}
		if (e.type == SDL_KEYDOWN)
		{
			if (query_escape(kbstate))
			{
				free_aarray_char(&save_games);
				free_aarray_char(&dates);
				welcome_loop();
			}
			if (query_enter(kbstate))
			{
				break;
			}
			if (!pressed)
			{
				pressed = 1;
				if (kbstate[SDL_SCANCODE_UP] == 1 && selection > 0)
				{
					selection -= 1;
					load_game_render_step(&save_games, &dates, selection);
				}
				if (kbstate[SDL_SCANCODE_DOWN] == 1 && selection < save_games.len - 1)
				{
					selection += 1;
					load_game_render_step(&save_games, &dates, selection);
				}
			}
		}
		if (e.type == SDL_KEYUP)
		{
			pressed = 0;
		}
	}
	char status = load_game(save_games.e[selection].e);
	if (status == 0)
	{
		free_aarray_char(&save_games);
		free_aarray_char(&dates);
		for (int i = 0; i < (*(sc.points)).len; i++)
		{
			update_t_scores(i);
		}
		electroscore = compute_electroscore();
		main_loop();
	}
	else if (sc_init)
	{
		printf("Unable to load %s\n", save_games.e[selection].e);
		free_sc_constr(&sc);
		free_sc_constr_interface(&scci);
	}
	/*
	 * int selection = choose_list_loop("Load game:", &save_games, base_x, base_y, 0.8);
	if (selection > -1)
	{
		char status = load_game(save_games.e[selection].e);
		if (status == 0)
		{
			free_aarray_char(&save_games);
			for (int i = 0; i < (*(sc.points)).len; i++)
			{
				update_t_scores(i);
			}
			electroscore = compute_electroscore();
			main_loop();
		}
		else if (sc_init)
		{
			printf("Unable to load %s\n", save_games.e[selection].e);
			free_sc_constr(&sc);
			free_sc_constr_interface(&scci);
		}
	}
	*/
}

void load_game_render_step(aarray_char *sgs, aarray_char *dates, int selection)
{
	double scale = 0.7;
	SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 0);
	SDL_RenderClear(rndrr);
	int base_x = (SCR_LEN_X / 20);
	int opp_base_x = (19 * SCR_LEN_X) / 20;
	int base_y = (SCR_LEN_Y / 20);
	SDL_SetRenderDrawColor(rndrr, 230, 230, 230, SDL_ALPHA_OPAQUE);
	render_string("Load Game:", 10, (SCR_LEN_X - 10 * DIGIT_WIDTH) / 2, SCR_LEN_Y / 10, 0, 1.0);
	base_y = SCR_LEN_Y / 10 + DIGIT_HEIGHT + 10;
	int incr_y = (int) (DIGIT_HEIGHT * scale) + 20;
	for (int i = 0; i < (*sgs).len; i++)
	{
		render_string((*sgs).e[i].e, (*sgs).e[i].len, base_x, base_y, 0, scale);
		render_string((*dates).e[i].e, (*dates).e[i].len, opp_base_x, base_y, 1, scale);
		if (i == selection)
		{
			SDL_Point box[5];
			box[0].x = base_x - 20; box[0].y = base_y - 5;
			box[1].x = opp_base_x + 20; box[1].y = box[0].y;
			box[2].x = box[1].x; box[2].y = base_y + incr_y - 5;
			box[3].x = box[0].x;
			box[3].y = box[2].y;
			box[4].x = box[0].x;
			box[4].y = box[0].y;
			SDL_SetRenderDrawColor(rndrr, 100, 100, 230, 255);
			SDL_RenderDrawLines(rndrr, box, 5);
		}
		base_y += incr_y;
	}
	SDL_RenderPresent(rndrr);
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

void set_pixel_dimensions()
{
	for (char i = 0; i < 33; i++)
	{
		ascii_pixel_width[i] = -1;
		ascii_pixel_height[i] = -1;
	}
	for (char i = 33; i < 127; i++)
	{
		Uint32 tex_fmt;
		int tex_acc;
		char imname[256];
		sprintf(imname, "ascii_letters/ascii_%d_c.bmp", i);
		SDL_Texture *lt = update_SDL_texture(imname, rndrr);
		int apw, aph;
		SDL_QueryTexture(lt, &tex_fmt, &tex_acc, &apw, &aph);
		ascii_pixel_width[i] = (char) apw;
		ascii_pixel_height[i] = (char) aph;
		ascii_format_offset_v[i] = DIGIT_HEIGHT - aph;
	}
	ascii_pixel_width[' '] = DIGIT_WIDTH;
	ascii_pixel_height[' '] = ascii_pixel_height['0'];
	ascii_format_offset_v['g'] = ascii_format_offset_v['a'];
	ascii_format_offset_v['y'] = ascii_format_offset_v['u'];
	ascii_format_offset_v['q'] = ascii_format_offset_v['g'];
	ascii_format_offset_v['j'] = ascii_format_offset_v['i'];
	ascii_format_offset_v['p'] = ascii_format_offset_v['k'];
	ascii_format_offset_v['-'] = DIGIT_HEIGHT / 2;
	ascii_format_offset_v['_'] = ascii_format_offset_v['.'];
	ascii_format_offset_v['\''] = (9 * DIGIT_HEIGHT) / 10;
	ascii_format_offset_v['"'] = ascii_format_offset_v['\''];
	ascii_format_offset_v['='] = DIGIT_HEIGHT / 2;
	ascii_format_offset_v['*'] = (7 * DIGIT_HEIGHT) / 10;
	ascii_format_offset_v['`'] = ascii_format_offset_v['\''];
	for (char i = 0; i < 10; i++)
	{
		char buf[32];
		sprintf(buf, "%d", i);
		digit_ascii_addr[i] = (int) buf[0];
	}
	char imname[256];
	sprintf(imname, "ascii_letters/ascii_%d_c.bmp", '-');
	SDL_Texture *lt = update_SDL_texture(imname, rndrr);
	Uint32 tex_fmt;
	int tex_acc;
	SDL_QueryTexture(lt, &tex_fmt, &tex_acc, &minus_pixel_width, &minus_pixel_height);
}



void render_ascii(char c, int base_x, int base_y, double scale)
{
	char imname[256];
	sprintf(imname, "ascii_letters/ascii_%d_c.bmp", c);
	SDL_Texture *lt = update_SDL_texture(imname, rndrr);
	int lt_w, lt_h, lt_acc;
	Uint32 lt_fmt;
	SDL_QueryTexture(lt, &lt_fmt, &lt_acc, &lt_w, &lt_h);
	SDL_Rect lt_rect;
	lt_rect.x = base_x;
	lt_rect.y = base_y;
	lt_rect.w = (int) (scale * lt_w);
	lt_rect.h = (int) (scale * lt_h);
	SDL_RenderCopy(rndrr, lt, NULL, &lt_rect);
}

void set_digits(int n, char *digits, int *len)
{
	(*len) = 0;
	while (n > 0)
	{
		digits[(*len)] = n % 10;
		n /= 10;
		(*len) += 1;
	}
}

void set_double_loop(char *msg, double *val, double incr, double _min_, double _max_, int base_x, int base_y)
{
	SDL_Event e;
	int prec = 2;
	int dec_pt_pos = 3;
	int p10 = 1000;
	int len = 0;
	double prox_val = (*val) * p10;
	double prox_incr = incr * p10;
	double prox_max_ = p10 * _max_;
	double prox_min_ = p10 * _min_;
	int delay = 5;
	int n0 = (int) prox_val;
	char digits[64];
	set_digits(n0, digits, &len);
	char pressed = 0;
	int count;
	set_double_render_step(msg, digits, len, dec_pt_pos, prec, base_x, base_y);
	while (SDL_WaitEvent(&e) >= 0)
	{
		query__exit_save(e);
		if (e.type == SDL_KEYDOWN)
		{
			if (pressed == 0)
			{
				if (kbstate[SDL_SCANCODE_UP] == 1) 
				{
					pressed = 1;
				}
				if (kbstate[SDL_SCANCODE_DOWN] == 1) 
				{
					pressed = -1;
				}
			}
			if (query_enter(kbstate))
			{
				(*val) = prox_val / p10;
				return;
			}
			else if (query_escape(kbstate))
			{
				return;
			}
		}
		if (e.type == SDL_KEYUP && pressed != 0)
		{
			if ((pressed == 1) && (prox_val < prox_max_)) 
			{
				prox_val += prox_incr;
			}
		       	if ((pressed == -1) && (prox_val > prox_min_)) 
			{
				prox_val -= prox_incr;
			}
			n0 = (int) prox_val;
			set_digits(n0, digits, &len);
			set_double_render_step(msg, digits, len, dec_pt_pos, prec, base_x, base_y);
			pressed = 0;
			count = 0;
		}
		if (pressed != 0)
		{
			if ((kbstate[SDL_SCANCODE_UP] == 1) && (prox_val < prox_max_)) 
			{
				count += 1;
			}
			if ((kbstate[SDL_SCANCODE_DOWN] == 1) && (prox_val > prox_min_)) 
			{
				count += 1;
			}
			if (count >= delay)
			{
				if ((prox_val > prox_min_) && (pressed == -1)) prox_val -= prox_incr;
				if ((prox_val < prox_max_) && (pressed == 1)) prox_val += prox_incr;
				count = 0;
				n0 = (int) prox_val;
				set_digits(n0, digits, &len);
				set_double_render_step(msg, digits, len, dec_pt_pos, prec, base_x, base_y);
			}
		}
		
	}
}

void render_sentence(aarray_char *s, int base_x, int base_y, char fb, double scale)
{
	base_x = base_x >= 0 ? base_x : 0;
	int base_x0 = base_x;
	for (int i = 0; i < (*s).len; i++)
	{
		int word_len = 0;
		for (int ii = 0; ii < (*s).e[i].len; ii++) word_len += (int) (scale * ascii_pixel_width[(*s).e[i].e[ii]]);
		if (base_x + word_len < SCR_LEN_X) {}
		else
		{
			base_x = base_x0;
			base_y += DIGIT_HEIGHT;
		}
		render_string((*s).e[i].e, (*s).e[i].len, base_x, base_y, 0, scale);
		base_x += word_len + (int) (15 * scale);
	}
}

void render_string(char *str, int str_len, int base_x, int base_y, char fb, double scale)
{
	base_x = base_x > 0 ? base_x : 0;
	base_y = base_y > 0 ? base_y : 0;
	int base_x0 = base_x;
	if (fb == 0)
	{
		for (int i = 0; i < str_len; i++)
		{
			if (str[i] != ' ')
			{
				int lwd, lht, lacc;
				Uint32 lfmt;
				char imname[256];
				sprintf(imname, "ascii_letters/ascii_%d_c.bmp", str[i]);
				SDL_Texture *img = update_SDL_texture(imname, rndrr);
				SDL_QueryTexture(img, &lfmt, &lacc, &lwd, &lht);
				render_image_box(img, base_x, base_y + (int) (scale * ascii_format_offset_v[str[i]]), &lwd, &lht, scale); 
				base_x += (int) (scale * lwd);
			}
			else
			{
				base_x += (int) (scale * DIGIT_WIDTH);
			}
			if (base_x > SCR_LEN_X)
			{
				base_x = base_x0;
				base_y += (int) (scale * DIGIT_HEIGHT) + 1;
			}
		}
	}
	else
	{
		int i = str_len;
		do 
		{
			i -= 1;
			if (str[i] != ' ')
			{
				int lwd, lht, lacc;
				Uint32 lfmt;
				char imname[256];
				sprintf(imname, "ascii_letters/ascii_%d_c.bmp", str[i]);
				SDL_Texture *img = update_SDL_texture(imname, rndrr);
				SDL_QueryTexture(img, &lfmt, &lacc, &lwd, &lht);
				base_x -= (int) (scale * lwd);
				if (base_x < 0)
				{
					base_x = base_x0 - (int) (scale * lwd);
					base_y += (int) (scale * DIGIT_HEIGHT);
				}
				render_image_box(img, base_x, base_y + (int) (scale * ascii_format_offset_v[str[i]]), &lwd, &lht, scale);
			}
			else
			{
				base_x -= (int) (scale * DIGIT_WIDTH);
				if (base_x < 0)
				{
					base_x = base_x0 - (int) (scale * DIGIT_WIDTH);
					base_y += (int) (scale * DIGIT_HEIGHT);
				}

			}
		} while (i > 0);
	}
}

void set_hole_width_render_step()
{
	SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 0);
	SDL_RenderClear(rndrr);
	SDL_SetRenderDrawColor(rndrr, 230, 230, 230, SDL_ALPHA_OPAQUE);
	render_string("Hole radius:", 10, SCR_LEN_X / 2 - 300, SCR_LEN_Y / 2 - 40, 0, 1.0);
	render_double(hole_wid, 2, (SCR_LEN_X + 300)/ 2, SCR_LEN_Y / 2 + DIGIT_HEIGHT, 0);
	SDL_RenderPresent(rndrr);
}

void set_hole_width_loop()
{
	char state = 0;
	char count = 0;
	char delay = 10;
	set_double_loop("Hole radius:", &hole_wid, hole_wid_incr, 0.01, 10.0, SCR_LEN_X / 2 - 300, SCR_LEN_Y / 2 - 40);
	return;
	/*
	set_hole_width_render_step();
	SDL_Event e;
	while (SDL_WaitEvent(&e) >= 0)
	{
		query__exit_save(e);
		if (e.type == SDL_KEYDOWN)
		{
			if (kbstate[SDL_SCANCODE_UP] == 1 && (state == 0))
			{
				state = 1;
			}
			if (kbstate[SDL_SCANCODE_DOWN] == 1 && (state == 0))
			{
				state = -1;
			}
			if (query_enter(kbstate))
			{
				hole_widsq = hole_wid * hole_wid;
				return;
			}
		}
		if (e.type == SDL_KEYUP)
		{
			if (state != 0)
			{
				hole_wid += state == 1 ? hole_wid_incr : -hole_wid_incr;
				count = 0;
				set_hole_width_render_step();
			}
			if (kbstate[SDL_SCANCODE_UP] == 0 && kbstate[SDL_SCANCODE_DOWN] == 0) state = 0;
		}
		if (state != 0)
		{
			count += 1;
			if (count > delay)
			{
				hole_wid += state == 1 ? hole_wid_incr : -hole_wid_incr;
				count = 0;
				set_hole_width_render_step();
			}
		}
	}
	*/
}

char query_exit(SDL_Event e)
{
	return e.type == SDL_QUIT || (e.type == SDL_KEYDOWN && kbstate[SDL_SCANCODE_Q] == 1);
}

void query__exit_save(SDL_Event e)
{
	char exit_flag = query_exit(e);
	if (exit_flag)
	{
		if (main_loop_init) save_game_loop();
		exit_program();
	}
}

void init_sentence(aarray_char *s, char *msg)
{
	int i = 0;
	while (1)
	{
		while (msg[i] != '\0' && msg[i] == ' ' || msg[i] == '\t')
		{
			i += 1;
		}
		if (msg[i] != '\0')
		{
			int wrd_addr = (*s).len;
			extend_aarray_char(s);
			while (msg[i] != '\0' && msg[i] != ' ' && msg[i] != '\t')
			{
				add2array_char(&((*s).e[wrd_addr]), msg[i]);
				i += 1;
			}
		}
		else break;
	}
}

void init_shift_map()
{
	for (int i = 0; i < 128; i++) shift_map[i] = -1;
	shift_map['-'] = '_';
	shift_map['\''] = '"';
	shift_map['='] = '+';
	for (int i = 97; i < 123; i++)
	{
		shift_map[i] = i - 32;
	}
}

// TBD


void menu_loop()
{

}

void menu_render_step()
{

}

int string_pixel_len(char *str, int len)
{
	int plen = 0;
	for (int i = 0; i < len; i++)
	{
		plen += ascii_pixel_width[str[i]];
	}
	return plen;
}

int choose_list_loop(char *ttl_msg, aarray_char *opt_list, int base_x, int base_y, double scale)
{
	int s = 0;
	choose_list_render_step(ttl_msg, opt_list, base_x, base_y, scale, s);
	SDL_Event e;
	char pressed = 0;
	int last_elem = (*opt_list).len - 1;
	while (SDL_WaitEvent(&e) >= 0)
	{
		if (query_exit(e))
		{
			if (!save_game_flag && main_loop_init) save_game_loop();
			exit_program();
		}
		if (e.type == SDL_KEYDOWN)
		{
			if (query_escape(kbstate))
			{
				return -1;
			}
			if (query_enter(kbstate))
			{
				return s;
			}
			if (!pressed)
			{
				pressed = 1;
				if (kbstate[SDL_SCANCODE_DOWN] == 1 && s < last_elem)
				{
					s += 1;
					choose_list_render_step(ttl_msg, opt_list, base_x, base_y, scale, s);
				}
				if (kbstate[SDL_SCANCODE_UP] == 1 && s > 0)
				{
					s -= 1;
					choose_list_render_step(ttl_msg, opt_list, base_x, base_y, scale, s);
				}
			}
		}
		if (e.type == SDL_KEYUP) pressed = 0;
	}
}

void choose_list_render_step(char *ttl_msg, aarray_char *opt_list, int base_x, int base_y, double scale, int s) 
{
	SDL_SetRenderDrawColor(rndrr, 0, 0, 0, SDL_ALPHA_OPAQUE);
	SDL_RenderClear(rndrr);
	SDL_SetRenderDrawColor(rndrr, 230, 230, 230, SDL_ALPHA_OPAQUE);
	// Render title message
	int msg_w, msg_h;
	int ttl_msg_len = strlen(ttl_msg);
	msg_w = string_pixel_len(ttl_msg, ttl_msg_len);
	msg_h = DIGIT_HEIGHT;
	render_string(ttl_msg, strlen(ttl_msg), base_x, base_y, 0, 1.0);
	base_y += msg_h + 10;
	SDL_Point sbox[5];
	int incr_y = (int) (scale * (DIGIT_HEIGHT + 10)) + 1;
	for (int i = 0; i < (*opt_list).len; i++)
	{
		if (i == s)
		{
			sbox[0].x = base_x; 
			sbox[0].y = base_y;
			sbox[1].x = base_x + (int) (scale * string_pixel_len((*opt_list).e[i].e, (*opt_list).e[i].len)); 
			sbox[1].y = base_y;
			sbox[2].x = sbox[1].x; sbox[2].y = base_y + incr_y - 3;
			sbox[3].x = sbox[0].x; sbox[3].y = sbox[2].y;
			sbox[4].x = sbox[0].x; sbox[4].y = sbox[0].y;
		}
		render_string((*opt_list).e[i].e, (*opt_list).e[i].len, base_x, base_y, 0, scale);
		base_y += incr_y;
	}
	SDL_SetRenderDrawColor(rndrr, 100, 100, 230, SDL_ALPHA_OPAQUE);
	SDL_RenderDrawLines(rndrr, sbox, 5);
	SDL_SetRenderDrawColor(rndrr, 230, 230, 230, SDL_ALPHA_OPAQUE);
	SDL_RenderPresent(rndrr);
}

