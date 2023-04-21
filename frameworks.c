/*
 *  frameworks.c
 *  ringdesign
 *
 *  Copyright 2009 Kyoto University. All rights reserved.
 *
 */

#include "frameworks.h"

//initialize array fstart
extern void getready_part(struct Particle *part, struct Cell *cell)
{
	//change part->fwk, assuming that we pass from part->fwr to cell->framework by a translation of cell->deltar along X 
	part->x += cell->deltar;	
	part->fwk = cell->framework;
	
	if(fabs(part->y) > 2*TINYLENGTH) {
		printf("!Warning in ini_fstart: y = %le is not 0 at the entrance of an Cell! is it what you really want?\n", part->y);
		printf("!!!ERROR in ini_fstart: y = %le is not 0 at the entrance of an Cell!!!\n ==>particle is marked as 'LOST_ERR1'\n", fabs(part->y));
		part->status = LOST_ERR1;
	}
	part->y = 0;
}

//update Particle part from values in fstart[]
extern void update_part(int awd_answer, struct Particle *part, struct Cell *cell)
{
	double shouldbeone;
	struct Framework fwk_exitface;
	struct Boundary exitboun;
	
	if(awd_answer == CUP) exitboun = cell->instru;
	else if(awd_answer == YES) exitboun = cell->boun;
	else errorstop("!!!ERROR in update_part: unexpected value of awd_answer");
	
	if(fabs(part->y*cos(exitboun.thmax) - part->x*sin(exitboun.thmax) - exitboun.ymax) > 2*TINYLENGTH) {
		printf("!!!ERROR in update_part: particle is not on the Cell boundary (y = %le)!!!\n ==>particle is marked as 'LOST_ERR2'\n", fabs(part->y*cos(exitboun.thmax) - part->x*sin(exitboun.thmax) - exitboun.ymax));
		part->status = LOST_ERR2;
	}
	
	shouldbeone = sqrt(part->ux*part->ux + part->uy*part->uy + part->uz*part->uz);
	//if(fabs(shouldbeone - 1) > TINYDIMLESS) printf("\n!Warning in update_part: fabs(ux*ux+uy*uy+uz*uz - 1) = %le, while you expect it is exactly zero... this will limit the computation accuracy\n", fabs(shouldbeone - 1));
	part->ux = part->ux/shouldbeone;
	part->uy = part->uy/shouldbeone;
	part->uz = part->uz/shouldbeone;
	
	//set part->fwk on the Cell exit face
	if(awd_answer == CUP) find_fwk_instru(&fwk_exitface, cell);	
	if(awd_answer == YES) find_fwk_exitface(&fwk_exitface, cell);	
	changefwk_part(part, &fwk_exitface);
}

//change framework: point, local to global
extern void fwk_pt_loctoglob(double *x, double *y, struct Framework *framework)
{
	double x_glob, y_glob;
	
	x_glob = framework->xc + *x*cos(framework->ae) - *y*sin(framework->ae);
	y_glob = framework->yc + *y*cos(framework->ae) + *x*sin(framework->ae);
	
	*x = x_glob;
	*y = y_glob;
}

//change framework: vector, local to global
extern void fwk_vect_loctoglob(double *x, double *y, struct Framework *framework)
{
	double x_glob, y_glob;
	
	x_glob = *x*cos(framework->ae) - *y*sin(framework->ae);
	y_glob = *y*cos(framework->ae) + *x*sin(framework->ae);
	
	*x = x_glob;
	*y = y_glob;
}

//change framework: point, global to local
extern void fwk_pt_globtoloc(double *x, double *y, struct Framework *framework)
{
	double x_loc, y_loc;
	
	x_loc = (*x - framework->xc)*cos(framework->ae) + (*y - framework->yc)*sin(framework->ae);
	y_loc = (*y - framework->yc)*cos(framework->ae) - (*x - framework->xc)*sin(framework->ae);
	
	*x = x_loc;
	*y = y_loc;
}

//change framework: point, global to local
extern void fwk_vect_globtoloc(double *x, double *y, struct Framework *framework)
{
	double x_loc, y_loc;
	
	x_loc = *x*cos(framework->ae) + *y*sin(framework->ae);
	y_loc = *y*cos(framework->ae) - *x*sin(framework->ae);
	
	*x = x_loc;
	*y = y_loc;
}

//change framework: point, local1 to local2
extern void fwk_pt_fwktofwk(double *x, double *y, struct Framework *fwk1, struct Framework *fwk2)
{
	fwk_pt_loctoglob(x, y, fwk1);
	fwk_pt_globtoloc(x, y, fwk2);
}

//change framework: vector, local1 to local2
extern void fwk_vect_fwktofwk(double *x, double *y, struct Framework *fwk1, struct Framework *fwk2)
{
	fwk_vect_loctoglob(x, y, fwk1);
	fwk_vect_globtoloc(x, y, fwk2);
}

//change the framework in which particles coordinates are expressed
extern void changefwk_part(struct Particle *part, struct Framework *newfwk)
{
	fwk_pt_fwktofwk(&(part->x), &(part->y), &(part->fwk), newfwk);
	fwk_vect_fwktofwk(&(part->ux), &(part->uy), &(part->fwk), newfwk);
	part->fwk = *newfwk;
}

extern void find_fwk_exitface(struct Framework *fwk_exitface, struct Cell *cell)
{
	// check for error in Boundary
	if(cell->boun.thmax != 0 && cell->boun.ymax != 0) errorstop("!!!ERROR in find_fwk_exitface: both boundary conditions are != 0 !!!\n");
	
	*fwk_exitface = cell->framework;
	
	fwk_exitface->ae += cell->boun.thmax;
	fwk_exitface->xc += -cell->boun.ymax*sin(cell->framework.ae);
	fwk_exitface->yc += cell->boun.ymax*cos(cell->framework.ae);
}

extern void find_fwk_instru(struct Framework *fwk_exitface, struct Cell *cell)
{
	// check for error in Boundary
	if(cell->instru.thmax != 0 && cell->instru.ymax != 0) errorstop("!!!ERROR in find_fwk_instru: both boundary conditions are != 0 !!!\n");
	
	*fwk_exitface = cell->framework;
	
	fwk_exitface->ae += cell->instru.thmax;
	fwk_exitface->xc += -cell->instru.ymax*sin(cell->framework.ae);
	fwk_exitface->yc += cell->instru.ymax*cos(cell->framework.ae);
}

//move the framework of every Cell of the latt so that what was the exit face of the latt becomes its entrance face
extern void latt_move_exit_to_entrance(struct Lattice *latt)
{
	int j;
	struct Framework fwk_exitface;
	
	find_fwk_exitface(&fwk_exitface, &(latt->cell[latt->nbcell - 1]));
	
	for(j = 0; j < latt->nbcell; j++) {
		shitf_fwk(&fwk_exitface, -1*latt->cell[j].deltar);
		latt->cell[j].framework = fwk_exitface;
		find_fwk_exitface(&fwk_exitface, &(latt->cell[j]));
	}
}

//shift framework origin of deltax its x axis (deltax > = < 0)
extern void shitf_fwk(struct Framework *fwk, double deltax)
{
	fwk->xc += deltax*cos(fwk->ae);
	fwk->yc += deltax*sin(fwk->ae);
}