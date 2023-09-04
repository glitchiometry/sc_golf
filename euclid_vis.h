#ifndef EUCLID_VIS_H
#define EUCLID_VIS_H
#include "euclid.h"
#include "SDL_video.h"
#include "SDL_rect.h"
#include "SDL_render.h"

typedef struct
{
	sc_constr *sc;
	double xbnds[2];
	double ybnds[2];
	int screen_len_x;
	int screen_len_y;
	array_double points_x;
	array_double points_y;
	array_int points_x_int;
	array_int points_y_int;
	array_char points_active;
	array_int active_points;
	array_voidstar circ_data;
	array_int active_circles;
	array_voidstar line_data;
	array_int active_lines;
} sc_constr_interface;

void sc_constr_interface_init(sc_constr_interface *scci, sc_constr *sc, double *xbnds, double *ybnds, int scr_len_x, int scr_len_y);

void free_sc_constr_interface(sc_constr_interface *scci);

void sc_constr_interface_resize(sc_constr_interface *scci, sc_constr *sc, double *xbnds, double *ybnds, int scr_len_x, int scr_len_y);
void add_point_sc_constr_interface(sc_constr_interface *scci, int point_addr);
void add_line_sc_constr_interface(sc_constr_interface *scci, int line_addr);
void add_circle_sc_constr_interface(sc_constr_interface *scci, int circ_addr);
typedef struct
{
	SDL_Point *bdry;
	SDL_Point center;
	int len;
	char vis;
} circle_render_data;

void circle_render_data_init(circle_render_data *cdata, sc_constr_interface *scci, double cx, double cy, double r);
void free_circle_render_data(circle_render_data *cdata);

typedef struct
{
	SDL_Point a;
	SDL_Point b;
	char vis;
} line_render_data;

void line_render_data_init(line_render_data *ldata, sc_constr_interface *scci, double ax, double ay, double bx, double by);
void line_rectangle_intersection_int(double ax, double ay, double bx, double by, sc_constr_interface *scci, double inv_wid_x, double inv_wid_y, int *xx, int *yy, int *x__, int *y__, char *status);

void segment_rectangle_intersection_int(double ax, double ay, char bdry_a, double bx, double by, char bdry_b, sc_constr_interface *scci, double inv_wid_x, double inv_wid_y, int *x_, int *y_, int *x__, int *y__);

void render_circle(circle_render_data *cdata, SDL_Renderer *rndrr);
void render_line(line_render_data *ldata, SDL_Renderer *rndrr);
void render_point(int x, int y, SDL_Renderer *rndrr);
void render_sc_constr(sc_constr_interface *scci, SDL_Renderer *rndrr);

void render_sc_constr_hlts(sc_constr_interface *scci, array_char *hltd_pts, array_char *hltd_lines, array_char *hltd_circles, SDL_Renderer *rndrr);
void sc_constr_interface_remove_last_point(sc_constr_interface *scci);
void sc_constr_interface_remove_last_line(sc_constr_interface *scci);
void sc_constr_interface_remove_last_circle(sc_constr_interface *scci);
#endif
