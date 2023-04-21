/*
 *  m_alloc.c
 *  ringdesign
 *
 *  Copyright 2009 Kyoto University. All rights reserved.
 *
 */

#include "m_alloc.h"

//----------------------------------------------------------------------//
//------------------------ Particle and Beam ---------------------------//
//----------------------------------------------------------------------//
extern void alloc_beam(struct Beam *beam, int npart)
{
	beam->npart = npart;
	if((beam->part = calloc(npart, sizeof(struct Particle))) == NULL) errorstop("memory allocation failure in alloc_beam");
}

extern void free_beam(struct Beam *beam)
{
	free(beam->part);
}

extern void copy_beam(struct Beam *copy, struct Beam *ini_beam)
{
	int i;
	
	alloc_beam(copy, ini_beam->npart);
	for(i = 0; i < ini_beam->npart; i++) copy->part[i] = ini_beam->part[i];
}

//----------------------------------------------------------------------//
//----------- Various Cells of struct Lattice ------------------//
//----------------------------------------------------------------------//
extern struct Cell *alloccell(int nbcell)
{
	struct Cell *cell;
	
	//allocate memory and display error if failed
	if((cell = calloc(nbcell, sizeof(struct Cell))) == NULL) errorstop("allocate failure in alloccell");
	return cell;
}

extern double **allocmatrix(int n1, int n2)
{
	int i;
	double **m;
	
	//1st dimension, allocate memory and display error if failed
	if((m = calloc(n1, sizeof(double*))) == NULL) errorstop("allocate failure 1 in allocmatrix");
	//2nd dimension, allocate memory and display error if failed
	for(i = 0; i < n1; i++) {
		if((m[i] = calloc((int)(n2), sizeof(double))) == NULL) errorstop("allocate failure 2 in allocmatrix");
	}
	return m;
}

extern struct Map allocmap(int nnode_r, int nnode_th, int nnode_z)
{
	int ir, ith;
	struct Map map;
	
	if ((map.node = calloc(nnode_r, sizeof(struct Mapnode**))) == NULL) errorstop("allocation failure 1 in allocmap");
	for(ir = 0; ir < nnode_r; ir++) {
		if ((map.node[ir] = calloc(nnode_th, sizeof(struct Mapnode*))) == NULL) errorstop("allocation failure 2 in allocmap");
		for(ith = 0; ith < nnode_th; ith++) {
			if ((map.node[ir][ith] = calloc(nnode_z, sizeof(struct Mapnode))) == NULL) errorstop("allocation failure 3 in allocmap");
		}
	}
	return map;
}

extern struct Cavity *alloccav(int periodicity)
{
	struct Cavity *cav;
	
	if((cav = calloc(periodicity, sizeof(struct Cavity))) == NULL) errorstop("allocate failure in alloccav");
	
	return cav;
}

extern void free_map(struct Map *map, int nnode_r, int nnode_th)
{
	int ir, ith;
	
	for(ir = 0; ir < nnode_r; ir++) {
		for(ith = 0; ith < nnode_th; ith++) {
			free(map->node[ir][ith]);
		}
		free(map->node[ir]);
	}
	free(map->node);
}

//free memory allocated in the various Cells of a struct Lattice
extern void free_latt(struct Lattice *latt)
{
	int i1, i2;
	
	if(latt->nbcell == 0) return;
	
	for(i1 = 0; i1 < latt->nbcell; i1++) {
		//if the Cell is a field map
		if(test_cell_map(&(latt->cell[i1])) == YES) free_map(&(latt->cell[i1].map), latt->cell[i1].map.nnodes[0], latt->cell[i1].map.nnodes[1]);
		
		//if the Cell is a cavity
		else if(strcmp(latt->cell[i1].keyword, "rf-thingap") == 0) free(latt->cell[i1].cav);
		
		//if the Cell is a thin-dipole
		else if(strcmp(latt->cell[i1].keyword, "thin-hdipole") == 0) {
			free(latt->cell[i1].mpara[0]);
			free(latt->cell[i1].mpara);
		}
		
		//else (except in case of a drift and collimator)
		else if(strcmp(latt->cell[i1].keyword, "drift") != 0 && strcmp(latt->cell[i1].keyword, "collimator") != 0 && strcmp(latt->cell[i1].keyword, "float-faraday") != 0){
			for(i2 = 0; i2 < latt->cell[i1].nbcomp; i2++) {
				free(latt->cell[i1].mpara[i2]);
				free(latt->cell[i1].efben[i2]);
				free(latt->cell[i1].efbex[i2]);
				free(latt->cell[i1].alierror[i2]);
			}
			free(latt->cell[i1].mpara);
			free(latt->cell[i1].efben);
			free(latt->cell[i1].efbex);
			free(latt->cell[i1].alierror);
		}
		//if the Cell is a drift or a collimator, there is nothing to free
	}
	free(latt->cell);
}


//copy a Lattice, with cell stored at a diffent memory address
extern void copy_latt(struct Lattice *copy, struct Lattice *ini_latt)
{
	int i;
	
	*copy = *ini_latt;
	//allocate memory for copy.cell
	copy->cell = alloccell(copy->nbcell);
	//copy all ini_latt->cell[i]
	for(i = 0; i < copy->nbcell; i++) {
		strcpy(copy->cell[i].keyword, ini_latt->cell[i].keyword);
		copy->cell[i] = ini_latt->cell[i];
	}
	//Warning: do not forget to free the allocated memory!
}

extern void copy_latt_1period(struct Lattice *reference_latt, struct Lattice *copy_latt)
{
	int i, j, nbcell, n1, n2;
	struct Framework fwk;
	
	
	//initialize copy_latt
	nbcell = reference_latt->nbcell*reference_latt->periodicity;
	copy_latt->periodicity = 1;
	copy_latt->nbcell = nbcell;
	copy_latt->cellnumber = 0;
	fwk = reference_latt->cell[0].framework;
	
	//allocate memory for err_latt->cell
	copy_latt->cell = alloccell(nbcell);
	for(j = 0; j < reference_latt->periodicity; j++) {
		for(i = 0; i < reference_latt->nbcell; i++) {
			
			copy_latt->cell[i+reference_latt->nbcell*j] = reference_latt->cell[i];
			shitf_fwk(&fwk, -1.*reference_latt->cell[i].deltar);
			copy_latt->cell[i+reference_latt->nbcell*j].framework = fwk;
			find_fwk_exitface(&fwk, &(copy_latt->cell[i+reference_latt->nbcell*j]));
			
			if(test_cell_map(&(reference_latt->cell[i])) == YES) {
				copy_cell_map(&(copy_latt->cell[i+reference_latt->nbcell*j]), &(reference_latt->cell[i]));
				//errorstop("Error in copy_latt_1period: memory allocation for a copy of a field map is not yet implemented");
			}
			
			else if(strcmp(reference_latt->cell[i].keyword, "rf-thingap") == 0) {
				//allocate memory for the err_latt->cell[i+reference_latt->nbcell*j].cav
				copy_latt->cell[i+reference_latt->nbcell*j].cav = alloccav(1);
				//initialize err_latt->cell[i+reference_latt->nbcell*j].cav
				copy_latt->cell[i+reference_latt->nbcell*j].cav[0] = reference_latt->cell[i].cav[j];
			}
			
			else if(strcmp(reference_latt->cell[i].keyword, "drift") == 0 || strcmp(reference_latt->cell[i].keyword, "collimator") == 0); //do nothing
			
			else {	
				//allocate memory and copy mpara, efben and efbex
				copy_latt->cell[i+reference_latt->nbcell*j].mpara	= allocmatrix(reference_latt->cell[i].nbcomp, SIZE_MPARA);
				for(n1 = 0; n1 < reference_latt->cell[i].nbcomp; n1++) {
					for(n2 = 0; n2 < SIZE_MPARA; n2++) {
						copy_latt->cell[i+reference_latt->nbcell*j].mpara[n1][n2] = reference_latt->cell[i].mpara[n1][n2];
					}
				}
				copy_latt->cell[i+reference_latt->nbcell*j].efben	= allocmatrix(reference_latt->cell[i].nbcomp, SIZE_EFB);
				for(n1 = 0; n1 < reference_latt->cell[i].nbcomp; n1++) {
					for(n2 = 0; n2 < SIZE_EFB; n2++) {
						copy_latt->cell[i+reference_latt->nbcell*j].efben[n1][n2] = reference_latt->cell[i].efben[n1][n2];
					}
				}
				copy_latt->cell[i+reference_latt->nbcell*j].efbex	= allocmatrix(reference_latt->cell[i].nbcomp, SIZE_EFB);
				for(n1 = 0; n1 < reference_latt->cell[i].nbcomp; n1++) {
					for(n2 = 0; n2 < SIZE_EFB; n2++) {
						copy_latt->cell[i+reference_latt->nbcell*j].efbex[n1][n2] = reference_latt->cell[i].efbex[n1][n2];
					}
				}
			}
		}
	}
}

//copy a Lattice, with cell stored at a different memory address with less cells
extern void copy_latt_n_first_cells(struct Lattice *copy, struct Lattice *ini_latt, int nb_cell)
{
	int i;
	
	*copy = *ini_latt;
	copy->nbcell = nb_cell;
	//allocate memory for copy.cell
	copy->cell = alloccell(copy->nbcell);
	//copy all ini_latt->cell[i]
	for(i = 0; i < copy->nbcell; i++) {
		strcpy(copy->cell[i].keyword, ini_latt->cell[i].keyword);
		copy->cell[i] = ini_latt->cell[i];
		if(test_cell_map(&(ini_latt->cell[i])) == YES) {
			copy_cell_map(&(copy->cell[i]), &(ini_latt->cell[i]));
		}
		
	}
	//Warning: do not forget to free the allocated memory!
}

//warning: not debugged
//copy cell with field map with memory allocation
extern void copy_cell_map(struct Cell *copycell, struct Cell *inicell)
{
	int i,j,k,n;
	
	printf("\n\n\n\n\twarning: copy_cell_map not debugged\n\n\n\n\n\n");
	copycell = inicell;
	copycell->map = allocmap(inicell->map.nnodes[0], inicell->map.nnodes[1], inicell->map.nnodes[2]);
	for(i=0;i<3;i++) copycell->map.nnodes[i] = inicell->map.nnodes[i];
	copycell->map.sym = inicell->map.sym;
	for(i=0;i<3;i++) copycell->map.stepsize[i] = inicell->map.stepsize[i];
	for(i=0;i<6;i++) copycell->map.mapdim[i] = inicell->map.mapdim[i];
	for(i=0;i<inicell->map.nnodes[0];i++) {
		for(j=0;j<inicell->map.nnodes[1];j++) {
			for(k=0;k<inicell->map.nnodes[2];k++) {
				for(n=0;n<3;n++) {
					copycell->map.node[i][j][k].coord[n] = inicell->map.node[i][j][k].coord[n];
					copycell->map.node[i][j][k].b[n] = inicell->map.node[i][j][k].b[n];
				}
			}
		}
	}
}

extern void gene_alignerror_latt(struct Lattice *reference_latt, struct Lattice *err_latt)
{
	int i;
	
	copy_latt_1period(reference_latt, err_latt);
	for(i=0;i<err_latt->nbcell;i++) {
		if(test_cell_map(&(err_latt->cell[i])) == YES) errorstop("alignment error in field maps not implemented!\n");
		err_latt->cell[i].alierror	= allocmatrix(err_latt->cell[i].nbcomp, SIZE_ALIERR);
	}
}
