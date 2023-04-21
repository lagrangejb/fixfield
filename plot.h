/*
 *  plot.h
 *  ringdesign
 *
 *  Copyright 2009 Kyoto University. All rights reserved.
 *
 */

#ifndef PLOT 
#define PLOT

#include "constants.h"
#include "common.h"


static void aspect_radial(char *filename, struct Cell *cell);
static void aspect_straight(char *filename, struct Cell *cell);
static void aspect_realbend(char *filename, struct Cell *cell);
static void aspect_purebend(char *filename, struct Cell *cell);
static void aspect_spiral(char *filename, struct Cell *cell);
static void aspect_tilt(char *filename, struct Cell *cell);
static void aspect_straighttilt(char *filename, struct Cell *cell);
static void change_coord_vffa_sect(double *x, double *y, double r, double th, int comp_nb, struct Cell *cell);
static void aspect_sect_vffa(char *filename, struct Cell *cell);
static void aspect_rect_vffa(char *filename, struct Cell *cell);
static void aspect_boun(char *filename, struct Cell *cell);
static void aspect_fringe_en_radial(char *filename, struct Cell *cell);
static void aspect_fringe_ex_radial(char *filename, struct Cell *cell);
static void aspect_fringe_en_straight(char *filename, struct Cell *cell);
static void aspect_fringe_ex_straight(char *filename, struct Cell *cell);
static void aspect_fringe_en_spiral(char *filename, struct Cell *cell);
static void aspect_fringe_ex_spiral(char *filename, struct Cell *cell);
static void aspect_fringe_en_tilt(char *filename, struct Cell *cell);
static void aspect_fringe_ex_tilt(char *filename, struct Cell *cell);
static void aspect_fringe_en_straighttilt(char *filename, struct Cell *cell);
static void aspect_fringe_ex_straighttilt(char *filename, struct Cell *cell);

#endif
