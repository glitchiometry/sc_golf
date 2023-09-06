#include "SDL.h"
#include "SDL_events.h"
#include "SDL_keyboard.h"
#include "euclid_vis.h"
#include "euclid.h"

#define N_T_PTS 3
#define MAX_ROOTED_PTS 2

int n_rooted_pts = 0;
double rooted_x[MAX_ROOTED_PTS];
double rooted_y[MAX_ROOTED_PTS];
int SCR_LEN_X = 960;
int SCR_LEN_Y = 540;
int vertex_mark_x[16] = {3, 4, 3, 2, 1, 0, -1, -2, -3, -4, -3, -2, -1, 0, 1, 2};
int vertex_mark_y[16] = {-1, 0, 1, 2, 3, 4, 3, 2, 1, 0, -1, -2, -3, -4, -3, -2};
const char *vdriver;

double rnd();
void closest_point(int x, int y, array_int *xs, array_int *ys, int *i_);
void closest_point_excluding(int x, int y, array_int *xs, array_int *ys, int *i_, int *excl, int len);
void closest_circle(double x, double y, sc_constr_interface *scci, int *i_);
void closest_line(double x, double y, sc_constr_interface *scci, int *i_);
void render_vertex(int x, int y, SDL_Renderer *rndrr);
void render_vertex_highlighted(int x, int y, SDL_Renderer *rndrr);

const int cutoff_dist_px = 100;
int cutoff_distsq;

void print_state_vars();
void reset_state_vars();
void set_cutoff_distsq();
void set_conv_factors();
void add_rooted_point(double x, double y);

double px_wid = 0.01;
double wid_x;
double wid_y;
double inv_wid_x;
double inv_wid_y;

char add_point_mode = 0;
char select_lr_flag = 0;
char select_mode = 0;
int select_curve = -1;
char apm_bit = 0;
char add_point = 0;
char add_line = 0;
char add_circle = 0;
char ctrl_mode = 0;
char zoom_mode = 0;

sc_constr sc;
int tally[4] = {0, 0, 0, 0};
double t_xs[N_T_PTS];
double t_ys[N_T_PTS];
double t_score[N_T_PTS];
sc_constr_interface scci;
sc_constr_interface scci_t;
array_char hltd_points;
array_char hltd_lines;
array_char hltd_circles;

void random_rect(double *xbnds, double *ybnds, double *x, double *y);

SDL_Texture *update_SDL_texture(char *img_name);

double _distsq_(double x1, double y1, double x2, double y2)
{
	double delxsq = x2 - x1;
	double delysq = y2 - y1;
	delxsq *= delxsq;
	delysq *= delysq;
	return delxsq + delysq;
}

void print_t_score()
{
	int tally_total = tally[0] + tally[1] + tally[2] + tally[3];
	printf("Scores: ");
	for (int i = 0; i < N_T_PTS; i++)
	{
		printf("%g, ", tally_total * sqrt(sqrt(t_score[i])));
	}
	printf("\n");
}

int main(int argc, char *argv[])
{
	// Initialize global variables
	double xybnds[2] = {SCR_LEN_X * px_wid, SCR_LEN_Y * px_wid};
	double xbnds_[2] = {0, xybnds[0]};
	double ybnds_[2] = {0, xybnds[1]};
	sc_constr_init(&sc);
	for (int i = 0; i < MAX_ROOTED_PTS; i++)
	{
		double x_, y_;
		random_rect(xbnds_, ybnds_, &x_, &y_);
		add_rooted_point(x_, y_);
	}
	for (int i = 0; i < N_T_PTS; i++)
	{
		random_rect(xbnds_, ybnds_, &t_xs[i], &t_ys[i]);
		t_score[i] = 1e99;
		for (int ii = 0; ii < MAX_ROOTED_PTS; ii++)
		{
			double distsq_i_ii = _distsq_(rooted_x[ii], rooted_y[ii], t_xs[i], t_ys[i]);
			t_score[i] = t_score[i] < distsq_i_ii ? t_score[i] : distsq_i_ii;
		}
	}
	sc_constr_interface_init(&scci, &sc, xbnds_, ybnds_, SCR_LEN_X, SCR_LEN_Y);
	set_cutoff_distsq();
	set_conv_factors();

	array_char_init(&hltd_points, 1);
	array_char_init(&hltd_lines, 1);
	array_char_init(&hltd_circles, 1);
	wid_x = (scci.xbnds[1] - scci.xbnds[0]) / scci.screen_len_x;
	wid_y = (scci.ybnds[1] - scci.ybnds[0]) / scci.screen_len_y;
	SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS | SDL_INIT_TIMER);
	// Get the name of the graphics driver:
	int N_vdrivers = SDL_GetNumVideoDrivers();
	if (N_vdrivers > 0) vdriver = SDL_GetVideoDriver(0);
	else
	{
		free_sc_constr(&sc);
		free_sc_constr_interface(&scci);
		free_array_char(&hltd_points);
		free_array_char(&hltd_lines);
		free_array_char(&hltd_circles);
		SDL_Quit();
		exit(EXIT_FAILURE);
	}
	SDL_VideoInit(vdriver);
	SDL_Window *win = SDL_CreateWindow("", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, SCR_LEN_X, SCR_LEN_Y, SDL_WINDOW_SHOWN);
	SDL_Renderer *rndrr = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
	SDL_Texture *tex = update_SDL_texture("baseline.bmp", rndrr);

	int i_, j_;
	double x1, y1, x2, y2;
	int x1_, y1_, x2_, y2_;
	double zm_xbnds[2], zm_ybnds[2];
	int prosp_x, prosp_y;
	double lr_x[2], lr_y[2];
	char lr_case;
	const Uint8 *kbstate = SDL_GetKeyboardState(NULL);
	while (1)
	{
		SDL_Event e;
		if (SDL_WaitEvent(&e)) {
			if (e.type == SDL_QUIT) break;
			if (e.type == SDL_MOUSEBUTTONDOWN)
			{
				print_state_vars();
				print_t_score();
			}
			if (kbstate[SDL_SCANCODE_ESCAPE] == 1)
			{
				reset_state_vars();
				tex = update_SDL_texture("baseline.bmp", rndrr);
			}
			if (kbstate[SDL_SCANCODE_LCTRL] == 1 || kbstate[SDL_SCANCODE_RCTRL] == 1 && !ctrl_mode)
			{
				ctrl_mode = 1;
			}
			if (ctrl_mode)
			{
				if (kbstate[SDL_SCANCODE_ESCAPE] == 1)
				{
					ctrl_mode = 0;
					continue;
				}
				if (kbstate[SDL_SCANCODE_Z] == 1)
				{
					zoom_mode = 1;
				}
				if (zoom_mode)
				{
					if (zoom_mode == 1 && e.type == SDL_MOUSEBUTTONDOWN)
					{
						x1_ = e.button.x;
						y1_ = e.button.y;
						double x__ = scci.xbnds[0] + x1_ * wid_x;
						double y__ = scci.ybnds[0] + y1_ * wid_y;
						double hl_x = (scci.xbnds[1] - scci.xbnds[0]) * 0.25;
						double hl_y = (scci.ybnds[1] - scci.ybnds[0]) * 0.25;
						zm_xbnds[0] = x__ - hl_x;
						zm_xbnds[1] = x__ + hl_x;
						zm_ybnds[0] = y__ - hl_y;
						zm_ybnds[1] = y__ + hl_y;
						printf("Zooming field of view to [%g, %g]x[%g, %g] centered at %g %g\n", zm_xbnds[0], zm_xbnds[1], zm_ybnds[0], zm_ybnds[1], x__, y__);
						sc_constr_interface_resize(&scci, &sc, zm_xbnds, zm_ybnds, scci.screen_len_x, scci.screen_len_y);
						set_conv_factors();
						zoom_mode = 0;
					}
					else if (zoom_mode == 1 && kbstate[SDL_SCANCODE_MINUS])
					{
						double hl_x = (scci.xbnds[1] - scci.xbnds[0]) * 0.5;
						double hl_y = (scci.ybnds[1] - scci.ybnds[0]) * 0.5;
						zm_xbnds[0] = scci.xbnds[0] - hl_x;
						zm_xbnds[1] = scci.xbnds[1] + hl_x;
						zm_ybnds[0] = scci.ybnds[0] - hl_y;
						zm_ybnds[1] = scci.ybnds[1] + hl_y;
						printf("Zooming out to [%g, %g]x[%g, %g]\n", zm_xbnds[0], zm_xbnds[1], zm_ybnds[0], zm_ybnds[1]);
						sc_constr_interface_resize(&scci, &sc, zm_xbnds, zm_ybnds, scci.screen_len_x, scci.screen_len_y);
						set_conv_factors();
						zoom_mode = 0;
					}
				}
				if (kbstate[SDL_SCANCODE_U] == 1 && sc.history.len > 0 && e.type == SDL_KEYDOWN)
				{
					tally[3] += 1;
					// Undo the last operation
					char last_op = sc.history.e[sc.history.len - 1];
					printf("Undoing operation %c\n", last_op);
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
			if (select_mode != 0)
			{
				if (select_mode == 1 || select_mode == 2 || select_mode == 3) {}
				else
				{
					printf("This message should not appear! Fix it ASAP!\n");
					exit(EXIT_FAILURE);
				}
				if (e.type == SDL_KEYDOWN)
				{
					if (kbstate[SDL_SCANCODE_UP] == 1)
					{
						if (select_mode == 1) 
						{
							hltd_lines.e[select_curve] = 0;
							do
							{
								select_curve = (select_curve + 1) % (scci.line_data).len;
							} while (apm_bit == 2 && (add_point_mode & 1) == 0 && select_curve == i_);
							hltd_lines.e[select_curve] = 1;
						}
						else if (select_mode == 2)
						{
							hltd_circles.e[select_curve] = 0;
							do
							{
								select_curve = (select_curve + 1) % (scci.circ_data).len;
							} while (apm_bit == 2 && (add_point_mode & 1) && select_curve == i_);
							hltd_circles.e[select_curve] = 1;
						}
						else
						{
							select_lr_flag = !select_lr_flag;
							lr_case = !lr_case;
							prosp_x = (int) ((lr_x[lr_case] - scci.xbnds[0]) * inv_wid_x);
							prosp_y = (int) ((lr_y[lr_case] - scci.ybnds[0]) * inv_wid_y);
						}
					}
					if (kbstate[SDL_SCANCODE_DOWN] == 1)
					{
						if (select_mode == 1) 
						{
							hltd_lines.e[select_curve] = 0;
							do
							{
								select_curve -= 1;
								select_curve = select_curve > -1 ? select_curve : (scci.line_data).len - 1;
							} while (apm_bit == 2 && (add_point_mode & 1) == 0 && select_curve == i_);
							hltd_lines.e[select_curve] = 1;
						}
						else if (select_mode == 2) 
						{
							hltd_circles.e[select_curve] = 0;
							do
							{
								select_curve -= 1;
								select_curve = select_curve > -1 ? select_curve : (scci.circ_data).len - 1;
							} while ((apm_bit == 2) && (add_point_mode & 1) == 0 && (select_curve == i_));
							hltd_circles.e[select_curve] = 1;
						}
						else
						{
							select_lr_flag = !select_lr_flag;
							lr_case = !lr_case;
							prosp_x = (int) ((lr_x[lr_case] - scci.xbnds[0]) * inv_wid_x);
							prosp_y = (int) ((lr_y[lr_case] - scci.ybnds[0]) * inv_wid_y);
						}
					}
					if (kbstate[SDL_SCANCODE_RETURN] == 1 || kbstate[SDL_SCANCODE_RETURN2] == 1)
					{
						if (select_mode == 3)
						{
							// Finalize the lr flag

							select_mode = 0;
						}
						else
						{
							if (apm_bit == 2)
							{
								printf("Setting second curve to %d\n", j_);
								j_ = select_curve;
								select_curve = -1;
								if (add_point_mode != 0)
								{
									select_lr_flag = 0;
									select_mode = 3;
									// Determine the prospective point(s)
									double x, y;
									if (add_point_mode == 1 || add_point_mode == 2)
									{
										line *l;
										circle *c;
										if (add_point_mode == 1)
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
										double x2mx1 = lr_x[1] - lr_x[0], y2my1 = lr_y[1] - lr_y[0];
										bx -= ax;
										by -= ay;
										double vdotdiff = bx * x2mx1 + by * y2my1;
										char test_vdiff = vdotdiff > 0;
										if (test_vdiff && !select_lr_flag || !test_vdiff && select_lr_flag)
										{
											x = lr_x[1];
											y = lr_y[1];
											lr_case = 1;
										}
										else if (test_vdiff && select_lr_flag || !test_vdiff && !select_lr_flag)
										{
											x = lr_x[0];
											y = lr_y[0];
											lr_case = 0;
										}
									}
									else if (add_point_mode == 3)
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
										double x1ref = lr_x[0] - c1x, y1ref = lr_y[0] - c1y, x2ref = lr_x[1] - c1x, y2ref = lr_y[1] - c1y;
										double x2px1 = x2ref + x1ref;
										double y2py1 = y2ref + y1ref;
										double cross_p = x2px1 * y1ref - x1ref * y2py1;
										char cross_test = cross_p > 0;
										if (cross_test && select_lr_flag || !cross_test && !select_lr_flag)
										{
											x = lr_x[1];
											y = lr_y[1];
											lr_case = 1;
										}
										else
										{
											x = lr_x[0];
											y = lr_y[0];
											lr_case = 0;
										}
									}
									prosp_x = (int) ((x - scci.xbnds[0]) * inv_wid_x);
									prosp_y = (int) ((y - scci.ybnds[0]) * inv_wid_y);
								}
								else 
								{
									select_mode = 0;
									select_curve = -1;
								}
							}
							else if (apm_bit == 1)
							{
								printf("Setting first curve to %d\n", i_);
								i_ = select_curve;
								select_mode = 0;
							}
							else
							{
								printf("This message should not appear! Fix it quickly! (apm_bit = %d)\n", apm_bit);
							}
							
						}
					}
				}
				else if (e.type == SDL_MOUSEBUTTONDOWN)
				{
					int x = e.button.x, y = e.button.y;
					double x_double = scci.xbnds[0] + x * px_wid;
					double y_double = scci.ybnds[0] + y * px_wid;
					if (select_mode == 1)
					{
						hltd_lines.e[select_curve] = 0;
						int orig_select_curve = select_curve;
						closest_line(x_double, y_double, &scci, &select_curve);
						if (apm_bit == 2 && add_point_mode & 1 == 0 && select_curve == i_) 
						{
							select_curve = orig_select_curve;
							printf("Singular intersections not supported.\n");
						}
						hltd_lines.e[select_curve] = 1;
					}
					else if (select_mode == 2)
					{
						hltd_circles.e[select_curve] = 0;
						int orig_select_curve = select_curve;
						closest_circle(x_double, y_double, &scci, &select_curve);
						if (apm_bit == 2 && add_point_mode & 1 && select_curve == i_)
						{
							select_curve = orig_select_curve;
						}
						hltd_circles.e[select_curve] = 1;
					}
				}
			}
			else
			{
				if (add_point)
				{
					if (apm_bit == 2)
					{
						if (e.type == SDL_KEYDOWN)
						{
							if (kbstate[SDL_SCANCODE_RETURN] == 1 || kbstate[SDL_SCANCODE_RETURN2] == 1)
							{
								// Add point depending on add_point_mode
								int point_addr = scci.points_x.len;
								point *p = (point *) calloc(1, sizeof(point));
								if (add_point_mode == 0)
								{
									// line-line
									add_point_sc_constr_ll(&sc, i_, j_);
									hltd_lines.e[i_] = 0;
									hltd_lines.e[j_] = 0;
								}
								else if (add_point_mode == 1)
								{
									// circle-line
									add_point_sc_constr_lc(&sc, j_, i_, select_lr_flag);
									hltd_lines.e[j_] = 0;
									hltd_circles.e[i_] = 0;
								}
								else if (add_point_mode == 2)
								{
									// line-circle
									add_point_sc_constr_lc(&sc, i_, j_, select_lr_flag);
									hltd_lines.e[i_] = 0;
									hltd_circles.e[j_] = 0;
								}
								else if (add_point_mode == 3)
								{
									// circle-circle
									hltd_circles.e[j_] = 0;
									hltd_circles.e[i_] = 0;
									add_point_sc_constr_cc(&sc, i_, j_, select_lr_flag);
								}
								add_point_sc_constr_interface(&scci, point_addr);
								add2array_char(&hltd_points, 0);
								tally[0] += 1;
								printf("Adding point at %g %g\n", scci.points_x.e[point_addr], scci.points_y.e[point_addr]);
								for (int iii = 0; iii < N_T_PTS; iii++)
								{
									double delsq_iii = _distsq_(scci.points_x.e[point_addr], scci.points_y.e[point_addr], t_xs[iii], t_ys[iii]);
									t_score[iii] = t_score[iii] < delsq_iii ? t_score[iii] : delsq_iii;
								}
								apm_bit = 0;
								add_point = 0;
								add_point_mode = 0;
								tex = update_SDL_texture("baseline.bmp", rndrr);
							}
							else if (kbstate[SDL_SCANCODE_ESCAPE] == 1)
							{
								reset_state_vars();
								tex = update_SDL_texture("baseline.bmp", rndrr);
							}
						}
					}
					else
					{
						if (e.type == SDL_KEYDOWN)
						{
							if (kbstate[SDL_SCANCODE_L] == 1) 
							{
								// Select line
								select_mode = 1;
								select_curve = 0;
								apm_bit += 1;
								hltd_lines.e[select_curve] = 1;
								tex = update_SDL_texture("addpoint_1.bmp", rndrr);
							}
							else if (kbstate[SDL_SCANCODE_C] == 1)
							{
								// Select circle
								select_mode = 2;
								select_curve = 0;
								add_point_mode |= (1 << apm_bit);
								apm_bit += 1;
								hltd_circles.e[select_curve] = 1;
								tex = update_SDL_texture("addpoint_2.bmp", rndrr);
							}
							else if (kbstate[SDL_SCANCODE_ESCAPE])
							{
								reset_state_vars();
								tex = update_SDL_texture("baseline.bmp", rndrr);
							}
						}
						else if (e.type == SDL_MOUSEBUTTONDOWN && apm_bit == 0 
								&& n_rooted_pts < MAX_ROOTED_PTS)
						{
							if (n_rooted_pts < MAX_ROOTED_PTS)
							{
								/*int x = e.button.x;
								int y = e.button.y;
								double x_ = scci.xbnds[0] + x * px_wid;
								double y_ = scci.ybnds[0] + y * px_wid;
								rooted_x[n_rooted_pts] = x_;
								rooted_y[n_rooted_pts] = y_;
								*/
								add2array_char(&hltd_points, 0);
								int addr_p = scci.points_x.len;
								point *p = (point *) calloc(1, sizeof(point));
								(*p).flag = 1;
								(*p).inc[0] = (void *) &rooted_x[n_rooted_pts];
								(*p).inc[1] = (void *) &rooted_y[n_rooted_pts];
								printf("Defining new point with coordinates %g %g\n", *((double *) (*p).inc[0]), *((double *) (*p).inc[1]));
								add2sc_constr(&sc, (void *) p, 'p');
								add_point_sc_constr_interface(&scci, addr_p);
								add_point = 0;
								n_rooted_pts += 1;
								tex = update_SDL_texture("baseline.bmp", rndrr);
							}
							else
							{
								printf("Unable to define more than %d rooted points\n", MAX_ROOTED_PTS);
							}
						}
					}
				}
				else
				{
					if (e.type == SDL_MOUSEBUTTONDOWN)
					{
						int x = e.button.x;
						int y = e.button.y;
						int cp_i;
						array_int *xs = &(scci.points_x_int);
						array_int *ys = &(scci.points_y_int);
						closest_point(x, y, xs, ys, &cp_i);
						printf("Mouse: %g %g (closest point: %d)\n", wid_x * x + scci.xbnds[0], y * wid_y + scci.ybnds[0], cp_i);

						if (kbstate[SDL_SCANCODE_L] == 1 && !add_line && (*xs).len > 1)
						{
							closest_point(x, y, xs, ys, &i_);
							if (i_ > -1) 
							{
								printf("Adding line...");
								add_line = 1;
								tex = update_SDL_texture("addline.bmp", rndrr);
							}
							else printf("No points in immediate vicinity\n");
						}
						else if (add_line)
						{
							closest_point_excluding(x, y, xs, ys, &j_, &i_, 1);
							if (j_ > -1)
							{
								add_line = 0;
								printf(" from %d to %d\n", i_, j_);
								int line_addr = (*(sc.lines)).len;
								add_line_sc_constr_pp(&sc, i_, j_);
								add_line_sc_constr_interface(&scci, line_addr);
								int x_a, y_a, x_b, y_b;
								line_render_data *lrd = (line_render_data *) scci.line_data.e[line_addr];
								x_a = (*lrd).a.x;
								y_a = (*lrd).a.y;
								x_b = (*lrd).b.x;
								y_b = (*lrd).b.y;
								printf("Line end points: (%g %g), (%g %g), vis = %d\n", scci.xbnds[0] + x_a * wid_x, scci.ybnds[0] + y_a * wid_y, scci.xbnds[0] + x_b * wid_x, scci.ybnds[0] + y_b * wid_y, (*lrd).vis);
								add2array_char(&hltd_lines, 0);
								tally[1] += 1;
								tex = update_SDL_texture("baseline.bmp", rndrr);
							}
						}
						else if (kbstate[SDL_SCANCODE_C] == 1 && !add_circle && (*xs).len > 1)
						{
							printf("Adding circle...");
							closest_point(x, y, xs, ys, &i_);
							if (i_ > -1) 
							{
								add_circle = 1;
								tex = update_SDL_texture("addcircle.bmp", rndrr);
							}
							else
							{
								printf("No points in immediate vicinity\n");
							}
						}
						else if (add_circle)
						{
							closest_point_excluding(x, y, xs, ys, &j_, &i_, 1);
							if (j_ > -1)
							{
								add_circle = 0;
								printf(" from %d to %d (circle index %d)\n", i_, j_, (*(sc.circles)).len);
								int circle_index = (*(sc.circles)).len;
								add_circle_sc_constr_pp(&sc, i_, j_);
								add_circle_sc_constr_interface(&scci, circle_index);
								add2array_char(&hltd_circles, 0);
								tally[2] += 1;
								tex = update_SDL_texture("baseline.bmp", rndrr);
							}
						}
					}
					else if (e.type == SDL_KEYDOWN)
					{
						if (kbstate[SDL_SCANCODE_P])
						{
							tex = update_SDL_texture("addpoint_0.bmp");
							add_point = 1;
							apm_bit = 0;
							add_point_mode = 0;
						}
					}
				}
			}
		}
		SDL_SetRenderDrawColor(rndrr, 0, 0, 0, 255);
		SDL_RenderClear(rndrr);
		SDL_RenderCopy(rndrr, tex, NULL, NULL);
		// Draw background
		SDL_SetRenderDrawColor(rndrr, 230, 230, 230, 0);
		// ...
		// Draw points, lines, circles
		// render_sc_constr(&scci, rndrr);
		// Render highlighted features
		// ...
		render_sc_constr_hlts(&scci, &hltd_points, &hltd_lines, &hltd_circles, rndrr);
		if (select_mode == 3)
		{
			render_vertex_highlighted(prosp_x, prosp_y, rndrr);
		}
		SDL_SetRenderDrawColor(rndrr, 230, 0, 230, 0);
		for (int i = 0; i < N_T_PTS; i++)
		{
			int x_i = (int) ((t_xs[i] - scci.xbnds[0]) * inv_wid_x);
			int y_i = (int) ((t_ys[i] - scci.ybnds[0]) * inv_wid_y);
			render_vertex(x_i, y_i, rndrr);
		}
		SDL_SetRenderDrawColor(rndrr, 230, 230, 230, 0);
		SDL_RenderPresent(rndrr);
	}
	free_array_char(&hltd_points);
	free_array_char(&hltd_lines);
	free_array_char(&hltd_circles);
	free_sc_constr(&sc);
	free_sc_constr_interface(&scci);
	SDL_DestroyRenderer(rndrr);
	SDL_DestroyWindow(win);
	SDL_VideoQuit();
	SDL_Quit();
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
	for (int i = 0; i < 16; i++)
	{
		int x_ = x + vertex_mark_x[i];
		int y_ = y + vertex_mark_y[i];
		SDL_RenderDrawPoint(rndrr, x_, y_);
		SDL_RenderDrawPoint(rndrr, x_ + vertex_mark_x[i], y_ + vertex_mark_y[i]);
	}
}

void add_rooted_point(double x, double y)
{
	printf("Adding rooted point at %g %g\n", x, y);
	rooted_x[n_rooted_pts] = x;
	rooted_y[n_rooted_pts] = y;

	point *p = (point *) calloc(1, sizeof(point));
	(*p).inc[0] = (void *) &rooted_x[n_rooted_pts];
	(*p).inc[1] = (void *) &rooted_y[n_rooted_pts];
	(*p).flag = 1;
	add2sc_constr(&sc, (void *) p, 'p');
	n_rooted_pts += 1;
}

void print_state_vars()
{
	printf("control mode: %d (z = %d), add_point: (%d, %d, %d; %d %d %d), add_line: %d, add_circle: %d \n", ctrl_mode, zoom_mode, add_point, add_point_mode, apm_bit, select_mode, select_curve, select_lr_flag, add_line, add_circle);
}

void reset_state_vars()
{
	add_point = add_point_mode = apm_bit = select_mode = select_lr_flag = add_line = add_circle = ctrl_mode = 0;
	select_curve = -1;
	zoom_mode = 0;
	for (int i = 0; i < hltd_points.len; i++)
	{
		hltd_points.e[i] = 0;
	}
	for (int i = 0; i < hltd_lines.len; i++) hltd_lines.e[i] = 0;
	for (int i = 0; i < hltd_circles.len; i++) hltd_circles.e[i] = 0;
}

void set_cutoff_distsq()
{
	cutoff_distsq = cutoff_dist_px * cutoff_dist_px;
}

void set_conv_factors()
{
	wid_x = (scci.xbnds[1] - scci.xbnds[0]) / scci.screen_len_x;
	wid_y = (scci.ybnds[1] - scci.ybnds[0]) / scci.screen_len_y;
	inv_wid_x = scci.screen_len_x / (scci.xbnds[1] - scci.xbnds[0]);
	inv_wid_y = scci.screen_len_y / (scci.ybnds[1] - scci.ybnds[0]);
	printf("pixel widths: %g %g\n", wid_x, wid_y);
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

