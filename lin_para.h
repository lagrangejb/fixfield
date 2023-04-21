/*
 *  lin_para.h
 *  ringdesign
 *
 *  Copyright 2010 Kyoto University. All rights reserved.
 *
 */

#ifndef LIN_PARA 
#define LIN_PARA

#include "constants.h"
#include "common.h"

static double halfintegertunex(struct Particle *reference, struct Lattice *latt, double ampxprime);
static double halfintegertunez(struct Particle *reference, struct Lattice *latt, double ampzprime);
//static void write_efbposition(char *filename, double s0, double x_clo, struct Lattice *latt);


#endif
