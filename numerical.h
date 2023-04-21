/*
 *  numerical.h
 *  ringdesign
 *
 *  Copyright 2009 Kyoto University. All rights reserved.
 *
 */

#ifndef NUMERICAL 
#define NUMERICAL

#include "constants.h"
#include "common.h"
#include "nrutil.h"

static void write_trackout(FILE *wfile, struct Particle *part, double bprobe[]);
static int arewedone(struct Particle *part, struct Cell *cell);
static void moveback_totheboun(struct Particle *part, struct Boundary *theboun, void (*derivs)(struct Particle*, double[], double[], struct Cell*,
	void(*)(double,double,double,double*,double*,double*,struct Cell*,int)), struct Cell *cell,
	void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
static int rk4(struct Particle *part, double dfds[], double h, void (*derivs)(struct Particle*, double[], double[], struct Cell*,
	void(*)(double,double,double,double*,double*,double*,struct Cell*,int)), struct Cell* cell,
	void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));


#endif
