/*
 *  bricabrac.c
 *  ringdesign
 *
 *  Copyright 2009 Kyoto University. All rights reserved.
 *
 */


#include "bricabrac.h"



// Integration of Bz along a straight line in a cell
extern void int_bz_liney_cell(double *intbdl, double *intbdlabs, double x, double ymin, double ymax, int nstep, struct Cell *cell,
							  void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double ystep, y, z, bprobe[4];
	
	//ystep = (ymax - ymin)/nstep; 
	ystep = comp_step(ymin, ymax, nstep);
	*intbdl = 0;		//int{bz.dl}
	*intbdlabs = 0;		//int{|bz|.dl}
	z = 0;				// allways in the midplane
	
	for(i = 0; i < nstep + 1; i++) {
		y = ymin + i*ystep;
		if(get_bfield(x, y, z, &bprobe[1], &bprobe[2], &bprobe[3], cell, add_contribution_comp) != ALIVE) errorstop("int_bz_liney_cell: get_bfield != ALIVE\n");
		*intbdl += bprobe[3]*ystep;
		*intbdlabs += fabs(bprobe[3])*ystep;
	}
}

// Integration of Bz along a radius in a cell
extern void int_bz_radius_cell(double *intbdl, double *intbdlabs, double r, double thmin, double thmax, int nstep, struct Cell *cell,
							  void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double thstep, th, z, bprobe[4], x, y;
	
	//thstep = (thmax - thmin)/nstep; 
	thstep = comp_step(thmin, thmax, nstep);
	*intbdl = 0;		//int{bz.dl}
	*intbdlabs = 0;		//int{|bz|.dl}
	z = 0;				// allways in the midplane
	
	for(i = 0; i < nstep + 1; i++) {
		th = thmin + i*thstep;
		x = r*cos(th);
		y = r*sin(th);
		if(get_bfield(x, y, z, &bprobe[1], &bprobe[2], &bprobe[3], cell, add_contribution_comp) != ALIVE) errorstop("int_bz_radius_cell: get_bfield != ALIVE\n");
		*intbdl += bprobe[3]*r*thstep;
		*intbdlabs += fabs(bprobe[3])*r*thstep;
	}
}

// Integration of Bz along a radius in a cell, and min and max value of Bz
extern void int_bz_max_min_radius_cell(double *intbdl, double *intbdlabs, double *bmax, double *bmin, double r, double thmin, double thmax, int nstep, struct Cell *cell,
							  void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double thstep, th, z, bprobe[4], x, y;
	
	*bmin = 1.e9;
	*bmax = -1.e9;
	thstep = comp_step(thmin, thmax, nstep);
	*intbdl = 0;		//int{bz.dl}
	*intbdlabs = 0;		//int{|bz|.dl}
	z = 0;				// always in the midplane
	
	for(i = 0; i < nstep + 1; i++) {
		th = thmin + i*thstep;
		x = r*cos(th);
		y = r*sin(th);
		if(get_bfield(x, y, z, &bprobe[1], &bprobe[2], &bprobe[3], cell, add_contribution_comp) != ALIVE) errorstop("int_bz_radius_cell: get_bfield != ALIVE\n");
		*intbdl += bprobe[3]*r*thstep;
		*intbdlabs += fabs(bprobe[3])*r*thstep;
		if(*bmin>bprobe[3]) *bmin = bprobe[3];
		if(*bmax<bprobe[3]) *bmax = bprobe[3];
	}
}

// Computation of m from Bz*l in a cell
extern void m_loc_bl_cell(char *outfilename, double ymin, double ymax, int nstep_y, double xmin, double xmax, int nstep_x, int nxref, double mref, struct Cell *cell,
						  void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double xstep, m[nstep_x+1], x[nstep_x+1], intbdl[nstep_x+1], intbdlabs[nstep_x+1], errorm[nstep_x+1], errorbdl[nstep_x+1], intbdldes[nstep_x+1];
	FILE *wfile;
	
	//calcul de BL
	//xstep = (xmax - xmin)/nstep_x;
	xstep = comp_step(xmin, xmax, nstep_x);
	for(i = 0; i < nstep_x + 1; i++) {
		x[i] = xmin + i*xstep;
		int_bz_liney_cell(&(intbdl[i]), &(intbdlabs[i]), x[i], ymin, ymax, nstep_y, cell, add_contribution_comp);
	}
	//calcul de m_loc
	for(i = 0; i < nstep_x + 1; i++) {
		if(i != nxref) {
			m[i] = log(intbdl[i]/intbdl[nxref])/(x[i] - x[nxref]);
			errorm[i] = (m[i] - mref)/mref;
		}
	}
	if(nxref != nstep_x) m[nxref] = (m[nxref - 1] + m[nxref + 1])/2;
	else m[nxref] = m[nxref - 1];
	//plot the desired Bl
	for(i = 0; i < nstep_x + 1; i++) {
		if(i != nxref) {
			intbdldes[i] = intbdl[nxref]*exp(mref*(x[i] - x[nxref]));
			errorbdl[i] = (fabs(intbdl[i]) - fabs(intbdldes[i]))/fabs(intbdldes[i]);
		}
	}
	intbdldes[nxref] = intbdl[nxref];
	errorbdl[nxref] = 0;
	errorm[nxref] = (m[nxref] - mref)/mref;
	//data writing
	wfile = fopen(outfilename,"w");
	for(i = 0; i < nstep_x + 1; i++) {
		fprintf(wfile, "%lf  %lf  %lf  %le  %le  %le\n", x[i], intbdl[i], m[i], errorm[i], intbdldes[i], errorbdl[i]);
	}
	fclose(wfile);
}

// Computation of k from Bz*l in a cell
extern void k_loc_bl_cell(char *outfilename, double thmin, double thmax, int nstep_th, double rmin, double rmax, int nstep_r, int nrref, double kref, struct Cell *cell,
						  void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double rstep, k[nstep_r+1], r[nstep_r+1], intbdl[nstep_r+1], intbdlabs[nstep_r+1], errork[nstep_r+1], errorbdl[nstep_r+1], intbdldes[nstep_r+1];
	FILE *wfile;
	
	//calculation of BL
	rstep = (rmax - rmin)/nstep_r;
	printf("rstep = %le\n",rstep);
	for(i = 0; i < nstep_r + 1; i++) {
		r[i] = rmin + i*rstep;
		printf("r[%i] = %lf\t",i, r[i]);
		int_bz_radius_cell(&(intbdl[i]), &(intbdlabs[i]), r[i], thmin, thmax, nstep_th, cell, add_contribution_comp);
		printf("intbdl[%i] = %le\n", i, intbdl[i]);
	}

	//calculation of k_loc
	for(i = 0; i < nstep_r + 1; i++) {
		if(i != nrref) {
			k[i] = log(intbdl[i]/intbdl[nrref])/log(r[i]/r[nrref])-1.;
			errork[i] = (k[i] - kref)/kref;
		}
	}
	if(nrref != nstep_r-1) k[nrref] = (k[nrref - 1] + k[nrref + 1])/2;
	else k[nrref] = k[nrref - 1];

	//data writing
	wfile = fopen(outfilename,"w");
	for(i = 0; i < nstep_r + 1; i++) fprintf(wfile, "%lf  %lf  %lf  %le  %le  %le\n", r[i], intbdl[i], k[i], errork[i], intbdldes[i], errorbdl[i]);
	fclose(wfile);
}

extern double comp_k_value_local(double r, double deltar, double bdeltar, double b) 
{
	return r*(bdeltar-b)/(b*deltar);
}

extern void k_loc_cell2(char *outfilename, double thmin, double thmax, int nstep_th, double rmin, double rmax, int nstep_r, double deltar, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double rstep, k_cell, k_d, k_f, k_f_cent, k_d_cent, r, intbdl, intbdlabs, bmin, bmax, intbdl_dr, intbdlabs_dr, bmin_dr, bmax_dr;
	FILE *wfile;
	
	//calculation of BL
	rstep = comp_step(rmin, rmax, nstep_r);
	//printf("rstep = %le\n",rstep);
	wfile = fopen(outfilename,"w");
	fprintf(wfile, "r, k_cell, k_f, k_d, k_f_cent, k_d_cent, intbdl, intbdl_abs, bmin, bmax, intbdl_dr, intbdl_abs_dr, bmin_dr, bmax_dr\n");
	
	for(i = 0; i < nstep_r + 1; i++) {
		r = rmin + i*rstep;
		int_bz_max_min_radius_cell(&(intbdl), &(intbdlabs), &bmax, &bmin, r, thmin, thmax, nstep_th, cell, add_contribution_comp);
		int_bz_max_min_radius_cell(&(intbdl_dr), &(intbdlabs_dr), &bmax_dr, &bmin_dr, r+deltar, thmin, thmax, nstep_th, cell, add_contribution_comp);
		k_cell = comp_k_value_local(r, deltar, intbdl_dr, intbdl);
		k_d = comp_k_value_local(r, deltar, (intbdlabs_dr+intbdl_dr)/2., (intbdlabs+intbdl)/2.); //positive part of the integral (D magnet)
		k_f = comp_k_value_local(r, deltar, (intbdlabs_dr-intbdl_dr)/2., (intbdlabs-intbdl)/2.); //negative part of the integral (F magnet)
		if (bmin<0 && bmin_dr<0) k_f_cent = comp_k_value_local(r, deltar, bmin_dr, bmin); 
		else k_f_cent = 0;
		if (bmax>0 && bmax_dr>0) k_d_cent = comp_k_value_local(r, deltar, bmax_dr, bmax); 
		else k_d_cent = 0;
		fprintf(wfile, "%lf	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le\n", r, k_cell, k_f, k_d, k_f_cent, k_d_cent, intbdl, intbdlabs, bmin, bmax, intbdl_dr, intbdlabs_dr, bmin_dr, bmax_dr);
	}
	fclose(wfile);
}


// Computation of m from Bz*l in a lattice
extern void m_loc_bl_latt(char *outfilename, int nstep_y, double xmin, double xmax, int nstep_x, int nxref, double mref, struct Lattice *latt,
						  void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i, j;
	double xstep, m[nstep_x+1], x[nstep_x+1], intbdl[nstep_x+1], intbdlabs[nstep_x+1], errorm[nstep_x+1], errorbdl[nstep_x+1], intbdldes[nstep_x+1], temp_intbdl, temp_intbdlabs;
	FILE *wfile;
	//calcul de BL
	xstep = (xmax - xmin)/nstep_x;
	for(i = 0; i < nstep_x + 1; i++) {
		x[i] = xmin + i*xstep;
		intbdl[i] = 0;
		intbdlabs[i] = 0;
		//printf("x[i] = %lf\n", x[i]);
		for(j = 0; j < latt->nbcell; j++) {
			if(strcmp(latt->cell[j].keyword, "drift") == 0) ;
			else {
				int_bz_liney_cell(&(temp_intbdl), &(temp_intbdlabs), x[i], 0., latt->cell[j].boun.ymax, nstep_y, &(latt->cell[j]), add_contribution_comp);
				intbdl[i] += temp_intbdl;
				intbdlabs[i] += temp_intbdlabs;
			}
		}
	}
	//calcul de m_loc
	for(i = 0; i < nstep_x + 1; i++) {
		if(i != nxref) {
			m[i] = log(intbdlabs[i]/intbdlabs[nxref])/(x[i] - x[nxref]);
			errorm[i] = (m[i] - mref)/mref;
		}
	}
	if(nxref != nstep_x) m[nxref] = (m[nxref - 1] + m[nxref + 1])/2;
	else m[nxref] = m[nxref - 1];
	//plot the desired Bl
	for(i = 0; i < nstep_x + 1; i++) {
		if(i != nxref) {
			intbdldes[i] = intbdlabs[nxref]*exp(mref*(x[i] - x[nxref]));
			errorbdl[i] = (intbdlabs[i] - intbdldes[i])/intbdldes[i];
		}
	}
	intbdldes[nxref] = intbdlabs[nxref];
	errorbdl[nxref] = 0;
	errorm[nxref] = (m[nxref] - mref)/mref;
	//data writing
	wfile = fopen(outfilename,"w");
	for(i = 0; i < nstep_x + 1; i++) fprintf(wfile, "%lf  %lf  %lf  %le  %le  %le\n", x[i], intbdlabs[i], m[i], errorm[i], intbdldes[i], errorbdl[i]);
	fclose(wfile);
}

// Computation of m from Bz in a cell (no integration)
extern void m_loc_y_cell(char *outfilename, double y, double xmin, double xmax, int nstep_x, int nxref, double mref, struct Cell *cell,
						 void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i, accept_flag;
	double xstep, m[nstep_x+1], x[nstep_x+1], bx, by, bz[nstep_x+1], bzabs[nstep_x+1], error[nstep_x+1];
	FILE *wfile;
	
	//calcul de B
	xstep = (xmax - xmin)/nstep_x;
	for(i = 0; i < nstep_x + 1; i++) {
		x[i] = xmin + i*xstep;
		accept_flag = get_bfield(x[i], y, 0, &bx, &by, &bz[i], cell, add_contribution_comp);
		if(accept_flag != ALIVE) errorstop("m_loc_y_cell: accept_flag != ALIVE\n");
		bzabs[i] = fabs(bz[i]);
	}
	
	//calcul de m
	for(i = 0; i < nstep_x + 1; i++) {
		if(i != nxref) {
			m[i] = log(bz[i]/bz[nxref])/(x[i] - x[nxref]);
			error[i] = (m[i] - mref)/mref;
		}
		
	}
	m[nxref] = (m[nxref - 1] + m[nxref + 1])/2;
	error[nxref] = (m[nxref] - mref)/mref;
	
	
	//data writing
	wfile = fopen(outfilename,"w");
	for(i = 0; i < nstep_x + 1; i++) {
		fprintf(wfile, "%lf  %lf  %lf  %le\n", x[i], bz[i], m[i], error[i]);
	}
	fclose(wfile);
}

//EXAMPLE of error study on B0
extern void scan_err_b0(char *outfilename7, char *outfilename3, char *outfile_final, struct Lattice *reference_latt, struct Particle *reference7, struct Particle *reference3, long *idum_err, int nblatt_err)
{
	int n;
	double rms, deltax7, deltaxp7, deltax3, deltaxp3, deltaz7, deltaz3, uhp_test7, uhp_test3, deltazp7, deltazp3, rmsdeltax7, rmsdeltaxp7, rmsdeltax3, rmsdeltaxp3, rmsdeltaz7, rmsdeltaz3, rmsdeltazp7, rmsdeltazp3;
	struct Particle test_part7, test_part3, part_ref7, part_ref3;
	struct Lattice err_latt;
	
	copy_latt_1period(reference_latt, &err_latt);
	FILE *rfile7;
	rfile7 = fopen(outfilename7,"w");
	FILE *rfile3;
	rfile3 = fopen(outfilename3,"w");
	FILE *wfile;
	wfile = fopen(outfile_final,"w");
	
	part_ref7 = *reference7;
	part_ref7.hat = -2;
	part_ref3 = *reference3;
	part_ref3.hat = -2;
	part_oneturn(&part_ref7, reference_latt,NULL);
	part_oneturn(&part_ref3, reference_latt,NULL);
	
	for(rms = 1.e-3; rms < 1.e-1; rms *= 10) {
		printf("rms = %.8f\n", rms);
		rmsdeltax7 = 0;
		rmsdeltaxp7 = 0;
		rmsdeltaz7 = 0;
		rmsdeltazp7 = 0;
		rmsdeltax3 = 0;
		rmsdeltaxp3 = 0;
		rmsdeltaz3 = 0;
		rmsdeltazp3 = 0;
		for(n = 0; n < nblatt_err; n++) {
			test_part7 = *reference7;
			test_part3 = *reference3;
			gene_b0_error_latt(&err_latt, idum_err, rms);
			//printf("n = %i\n", n);
			//find_closed_orbite_x(&(test_part7), &(test_part7.x), 1.e-3, &err_latt, NO);
			part_oneturn(&(test_part7), &err_latt,NULL);
			//printf("test_part11.x = %lf\n", test_part7.x);
			//find_closed_orbite_x(&(test_part3), &(test_part3.x), 1.e-3, &err_latt, NO);
			part_oneturn(&(test_part3), &err_latt,NULL);
			deltax7 = fabs(test_part7.x - part_ref7.x);
			deltax3 = fabs(test_part3.x - part_ref3.x);
			deltaxp7 = fabs(atan_ratio(test_part7.ux, test_part7.uy) - atan_ratio(part_ref7.ux, part_ref7.uy));
			deltaxp3 = fabs(atan_ratio(test_part3.ux, test_part3.uy) - atan_ratio(part_ref3.ux, part_ref3.uy));
			deltaz7 = fabs(test_part7.z);
			deltaz3 = fabs(test_part3.z);
			uhp_test7 = sqrt(pow(test_part7.ux, 2) + pow(test_part7.uy, 2));
			uhp_test3 = sqrt(pow(test_part3.ux, 2) + pow(test_part3.uy, 2));
			deltazp7 = fabs(atan_ratio(test_part7.uz, uhp_test7));
			deltazp3 = fabs(atan_ratio(test_part3.uz, uhp_test3));
			
			fprintf(rfile7, "%le  \t  %li  \t  %lf  \t  %lf  \t  %lf  \t  %lf\n", rms, *idum_err, deltax7, deltaxp7, deltaz7, deltazp7);
			fprintf(rfile3, "%le  \t  %li  \t  %lf  \t  %lf  \t  %lf  \t  %lf\n", rms, *idum_err, deltax3, deltaxp3, deltaz3, deltazp3);
			rmsdeltax7 += deltax7*deltax7;
			rmsdeltaxp7 += deltaxp7*deltaxp7;
			rmsdeltaz7 += deltaz7*deltaz7;
			rmsdeltazp7 += deltazp7*deltazp7;
			rmsdeltax3 += deltax3*deltax3;
			rmsdeltaxp3 += deltaxp3*deltaxp3;
			rmsdeltaz3 += deltaz3*deltaz3;
			rmsdeltazp3 += deltazp3*deltazp3;
			//	printf("rmsdeltax7 = %le\n", rmsdeltax7);
			//	printf("rmsdeltaxp7 = %le\n", rmsdeltaxp7);
			//	printf("rmsdeltaz7 = %le\n", rmsdeltaz7);
			//	printf("rmsdeltazp7 = %le\n", rmsdeltazp7);
		}
		rmsdeltax7 = sqrt(rmsdeltax7/nblatt_err);
		rmsdeltaxp7 = sqrt(rmsdeltaxp7/nblatt_err);
		rmsdeltaz7 = sqrt(rmsdeltaz7/nblatt_err);
		rmsdeltazp7 = sqrt(rmsdeltazp7/nblatt_err);
		rmsdeltax3 = sqrt(rmsdeltax3/nblatt_err);
		rmsdeltaxp3 = sqrt(rmsdeltaxp3/nblatt_err);
		rmsdeltaz3 = sqrt(rmsdeltaz3/nblatt_err);
		rmsdeltazp3 = sqrt(rmsdeltazp3/nblatt_err);
		//	printf("\n\nfinal:\n");
		//	printf("rmsdeltax7 = %le\n", rmsdeltax7);
		//	printf("rmsdeltaxp7 = %le\n", rmsdeltaxp7);
		//	printf("rmsdeltaz7 = %le\n", rmsdeltaz7);
		//	printf("rmsdeltazp7 = %le\n", rmsdeltazp7);
		fprintf(wfile, "%le  \t  %le  \t  %le  \t  %le  \t  %le  \t  %le  \t  %le  \t  %le  \t  %le\n", rms, rmsdeltax7, rmsdeltaxp7, rmsdeltaz7, rmsdeltazp7, rmsdeltax3, rmsdeltaxp3, rmsdeltaz3, rmsdeltazp3);
	}
	
	fclose(rfile7);
	fclose(rfile3);
	fclose(wfile);
	free_latt(&err_latt);
}

// generate an gaussian error on B0 for each component of the lattice 
extern void gene_b0_error_latt(struct Lattice *err_latt, long *idum, double rms_percentage_error)
{
	int i, n1;
	double a;
	
	for(i = 0; i < err_latt->nbcell; i++) {
		if(test_cell_map(&(err_latt->cell[i])) == NO && strcmp(err_latt->cell[i].keyword, "rf-thingap") != 0  && strcmp(err_latt->cell[i].keyword, "drift") != 0  && strcmp(err_latt->cell[i].keyword, "collimator") != 0) {
			for(n1 = 0; n1 < err_latt->cell[i].nbcomp; n1++) {
				a = gasdev(idum);
				a = a*rms_percentage_error;
				err_latt->cell[i].mpara[n1][2] = err_latt->cell[i].mpara[n1][2] * (1. + a);
			}
		}
	}
}

//EXAMPLE of alignment error study 
extern void scan_error(char *outfilename7, char *outfilename3, char *outfile_final, struct Lattice *reference_latt, struct Particle *reference7, struct Particle *reference3, long *idum_err, int nblatt_err)
{
	int i, n;
	double rms, deltax7, deltaxp7, deltax3, deltaxp3, deltaz7, deltaz3, uhp_test7, uhp_test3, deltazp7, deltazp3, rmsdeltax7, rmsdeltaxp7, rmsdeltax3, rmsdeltaxp3, rmsdeltaz7, rmsdeltaz3, rmsdeltazp7, rmsdeltazp3;
	struct Particle test_part7, test_part3, part_ref7, part_ref3;
	struct Lattice err_latt;
	
	copy_latt_1period(reference_latt, &err_latt);
	for(i = 0; i < reference_latt->nbcell; i++) {
		if(strcmp(reference_latt->cell[i].keyword, "rf-thingap") == 0 || strcmp(reference_latt->cell[i].keyword, "drift") == 0 || strcmp(reference_latt->cell[i].keyword, "collimator") == 0); //do nothing
		else err_latt.cell[i].alierror = allocmatrix(reference_latt->cell[i].nbcomp, SIZE_ALIERR);
	}
	
	FILE *rfile7;
	rfile7 = fopen(outfilename7,"w");
	FILE *rfile3;
	rfile3 = fopen(outfilename3,"w");
	FILE *wfile;
	wfile = fopen(outfile_final,"w");
	
	part_ref7 = *reference7;
	part_ref7.hat = -2;
	part_ref3 = *reference3;
	part_ref3.hat = -2;
	part_oneturn(&part_ref7, reference_latt,NULL);
	part_oneturn(&part_ref3, reference_latt,NULL);
	//printf("part_ref7.x = %lf\n", part_ref7.x);
	
	for(rms = 1.e-6; rms < 1.e-1; rms *= 10) {
		printf("rms = %.8f\n", rms);
		rmsdeltax7 = 0;
		rmsdeltaxp7 = 0;
		rmsdeltaz7 = 0;
		rmsdeltazp7 = 0;
		rmsdeltax3 = 0;
		rmsdeltaxp3 = 0;
		rmsdeltaz3 = 0;
		rmsdeltazp3 = 0;
		for(n = 0; n < nblatt_err; n++) {
			test_part7 = *reference7;
			test_part3 = *reference3;
			gene_alignerror_latt_old(reference_latt, &err_latt, idum_err, rms, 0., 2.175);
			//printf("n = %i\n", n);
			part_oneturn(&(test_part7), &err_latt,NULL);
			//printf("test_part11.x = %lf\n", test_part7.x);
			part_oneturn(&(test_part3), &err_latt,NULL);
			deltax7 = fabs(test_part7.x - part_ref7.x);
			deltax3 = fabs(test_part3.x - part_ref3.x);
			deltaxp7 = fabs(atan_ratio(test_part7.ux, test_part7.uy) - atan_ratio(part_ref7.ux, part_ref7.uy));
			deltaxp3 = fabs(atan_ratio(test_part3.ux, test_part3.uy) - atan_ratio(part_ref3.ux, part_ref3.uy));
			deltaz7 = fabs(test_part7.z);
			deltaz3 = fabs(test_part3.z);
			uhp_test7 = sqrt(pow(test_part7.ux, 2) + pow(test_part7.uy, 2));
			uhp_test3 = sqrt(pow(test_part3.ux, 2) + pow(test_part3.uy, 2));
			deltazp7 = fabs(atan_ratio(test_part7.uz, uhp_test7));
			deltazp3 = fabs(atan_ratio(test_part3.uz, uhp_test3));
			
			fprintf(rfile7, "%le  \t  %li  \t  %lf  \t  %lf  \t  %lf  \t  %lf\n", rms, *idum_err, deltax7, deltaxp7, deltaz7, deltazp7);
			fprintf(rfile3, "%le  \t  %li  \t  %lf  \t  %lf  \t  %lf  \t  %lf\n", rms, *idum_err, deltax3, deltaxp3, deltaz3, deltazp3);
			rmsdeltax7 += deltax7*deltax7;
			rmsdeltaxp7 += deltaxp7*deltaxp7;
			rmsdeltaz7 += deltaz7*deltaz7;
			rmsdeltazp7 += deltazp7*deltazp7;
			rmsdeltax3 += deltax3*deltax3;
			rmsdeltaxp3 += deltaxp3*deltaxp3;
			rmsdeltaz3 += deltaz3*deltaz3;
			rmsdeltazp3 += deltazp3*deltazp3;
			//	printf("rmsdeltax7 = %le\n", rmsdeltax7);
			//	printf("rmsdeltaxp7 = %le\n", rmsdeltaxp7);
			//	printf("rmsdeltaz7 = %le\n", rmsdeltaz7);
			//	printf("rmsdeltazp7 = %le\n", rmsdeltazp7);
		}
		rmsdeltax7 = sqrt(rmsdeltax7/nblatt_err);
		rmsdeltaxp7 = sqrt(rmsdeltaxp7/nblatt_err);
		rmsdeltaz7 = sqrt(rmsdeltaz7/nblatt_err);
		rmsdeltazp7 = sqrt(rmsdeltazp7/nblatt_err);
		rmsdeltax3 = sqrt(rmsdeltax3/nblatt_err);
		rmsdeltaxp3 = sqrt(rmsdeltaxp3/nblatt_err);
		rmsdeltaz3 = sqrt(rmsdeltaz3/nblatt_err);
		rmsdeltazp3 = sqrt(rmsdeltazp3/nblatt_err);
		//	printf("\n\nfinal:\n");
		//	printf("rmsdeltax7 = %le\n", rmsdeltax7);
		//	printf("rmsdeltaxp7 = %le\n", rmsdeltaxp7);
		//	printf("rmsdeltaz7 = %le\n", rmsdeltaz7);
		//	printf("rmsdeltazp7 = %le\n", rmsdeltazp7);
		fprintf(wfile, "%le  \t  %le  \t  %le  \t  %le  \t  %le  \t  %le  \t  %le  \t  %le  \t  %le\n", rms, rmsdeltax7, rmsdeltaxp7, rmsdeltaz7, rmsdeltazp7, rmsdeltax3, rmsdeltaxp3, rmsdeltaz3, rmsdeltazp3);
	}
	fclose(rfile7);
	fclose(rfile3);
	fclose(wfile);
	free_latt(&err_latt);
}

// adjust the fields of the element to keep the same FD ratio, with the particle on the wanted closed orbit (r0)
// ! only works with a 1-element lattice !
// works also with straights : r0 = x0 and k = novero.
extern int adjustclosedorbit(struct Lattice *latt, struct Particle *part, double r0, double k, double *b00, double *b01)
{
	int i, j;
	double rnew, nux, nuz, betax, alphax, betaz, alphaz, b0origin0, b0origin1;

	b0origin0 = latt->cell[0].mpara[0][2];
	b0origin1 = latt->cell[0].mpara[1][2];
	
	for (i = 0; i < 10; i++) {
		//if (find_closed_orbite_xxp(part, &rnew, &part->ux, &part->uy, 1.e-5, latt, NO) == TRUE) {
		if (find_closed_orbite_x(part, &rnew, 5.e-6, latt, NO) == TRUE) {
			part->x = rnew;
			tune_calc_matrix(part, &nux, &nuz, &betax, &alphax, &betaz, &alphaz, 1.e-4, 1.e-5, 1.e-4, 1.e-5, latt, part_cross_latt, NO, NULL);
			if (fabs(rnew - r0) < 1.e-7) {
				printf("\n\n\t\t\tb0adjust = %lf and %lf [T]\n\n\n", latt->cell[0].mpara[0][2], latt->cell[0].mpara[1][2]);
				*b00 = latt->cell[0].mpara[0][2];
				*b01 = latt->cell[0].mpara[1][2];
				return TRUE;
			}
			//for(j = 0; j < latt->cell[0].nbcomp; j++) printf("latt->cell[0].mpara[%i][2] = %lf\n", j, latt->cell[0].mpara[j][2]);
			//printf("rnew = %lf, r0 = %lf, k = %lf, exp(k*(rnew - r0)) = %lf\n", rnew, r0, k, exp(k*(rnew - r0)));
			if(strcmp(latt->cell[0].keyword, "ffag-r-he") == 0 ||
			   strcmp(latt->cell[0].keyword, "ffag-r-lin") == 0 ||
			   strcmp(latt->cell[0].keyword, "ffag-r-enge") == 0 ||
			   strcmp(latt->cell[0].keyword, "ffag-spi-he") == 0 ||
			   strcmp(latt->cell[0].keyword, "ffag-spi-lin") == 0 ||
			   strcmp(latt->cell[0].keyword, "ffag-spi-enge") == 0 ||
			   strcmp(latt->cell[0].keyword, "ffag-spi-fullenge") == 0)  for(j = 0; j < latt->cell[0].nbcomp; j++) latt->cell[0].mpara[j][2] = latt->cell[0].mpara[j][2]*pow(rnew/r0, k);
			
			else if(strcmp(latt->cell[0].keyword, "ffag-s-lin") == 0 ||
					strcmp(latt->cell[0].keyword, "ffag-s-he") == 0 ||
					strcmp(latt->cell[0].keyword, "ffag-s-enge") == 0) for(j = 0; j < latt->cell[0].nbcomp; j++) latt->cell[0].mpara[j][2] = latt->cell[0].mpara[j][2]*exp(k*(rnew - r0));
			else printf("\n \nproblem in adjustclosedorbitbend, cannot read keyword!\n");
			//for(j = 0; j < latt->cell[0].nbcomp; j++) printf("latt->cell[0].mpara[%i][2] = %lf\n", j, latt->cell[0].mpara[j][2]); 
		}
		else {
			printf("\n \nproblem in adjustclosedorbitbend, closed orbit not found!\n");
			latt->cell[0].mpara[0][2] = b0origin0;
			latt->cell[0].mpara[1][2] = b0origin1;
			return FALSE;
		}
	}
	printf("\n \nin adjustclosedorbitbend precision not achieved: eps = %le, increase the number of turns or decrease demanded precision\n \n", fabs(rnew - r0));
	printf("\n\n\t\t\tb0adjust = %lf and %lf [T]\n\n\n", latt->cell[0].mpara[0][2], latt->cell[0].mpara[1][2]);
	latt->cell[0].mpara[0][2] = b0origin0;
	latt->cell[0].mpara[1][2] = b0origin1;
	return FALSE;
}

extern void iteration_k2(char *outfilename, struct Lattice *latt, struct Particle *reference_part, double eps_clo, double k2step, double k2max, double xref)
{
	int i, n, m, accept_flag, flag;
	double qx1, qz1, betax, alphax, betaz, alphaz, qx2, qz2, brhomin, brhomax, dqx, dqz;
	struct Particle test_part;
	FILE *wfile;
	wfile = fopen(outfilename,"w");
	
	do {
		test_part = *reference_part;
		test_part.hat = -2;
		flag = 0;
		//change k2
		for(n=0;n<latt->cell[0].nbcomp;n++) latt->cell[0].mpara[n][4] += k2step;
		printf("\n\nk2 = %lf\n", latt->cell[0].mpara[0][4]);
		//find energy which closed orbit = xref
		accept_flag = find_closed_orbite_xxp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), eps_clo, latt, NO);
		if(accept_flag != TRUE) errorstop("!!!ERROR in iteration_k2 (1), closed orbit not found");
		
		if(test_part.x>xref) {
			brhomax = test_part.brho;
			m = 0;
			for(;;) {
				//printf("\n");
				test_part.brho *=0.98;
				accept_flag = find_closed_orbite_xxp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), eps_clo, latt, NO);
				if(accept_flag != TRUE) errorstop("!!!ERROR in iteration_k2 (2), closed orbit not found");
				if(test_part.x<xref) {
					brhomin = test_part.brho;
					//printf("brhomin = %lf, brhomax = %lf\n", brhomin, brhomax);
					break;
				}
				m++;
				if(m>5) errorstop("!!!ben alors ca marche ou pas ?");
			}
		}
		else if(test_part.x<xref) {
			brhomin = test_part.brho;
			m = 0;
			for(;;) {
				test_part.brho *=1.02;
				accept_flag = find_closed_orbite_xxp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), eps_clo, latt, NO);
				if(accept_flag != TRUE) errorstop("!!!ERROR in iteration_k2 (2), closed orbit not found");
				if(test_part.x>xref) {
					brhomax = test_part.brho;
					//printf("brhomin = %lf, brhomax = %lf\n", brhomin, brhomax);
					break;
				}
				m++;
				if(m>5) errorstop("!!!ben alors ca marche ou pas (2)?");
			}
		}
		else if(test_part.x == xref) flag = 1;
		i = 0;
		for(;;) {
			//printf("entree dans la boucle brhomin-brhomax\n");
			if(flag == 1) break;
			test_part.brho = (brhomax+brhomin)/2.;
			//printf("brhomin = %lf, brhomax = %lf, test_part.brho = %lf\n", brhomin, brhomax, test_part.brho);
			accept_flag = find_closed_orbite_xxp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), eps_clo, latt, NO);
			if(accept_flag != TRUE) errorstop("!!!ERROR in iteration_k2 (3), closed orbit not found");
			if(fabs(test_part.x-xref)<1.e-3) break;
			
			if(test_part.x>xref) brhomax = test_part.brho;
			else if(test_part.x<xref) brhomin = test_part.brho;
			
			if(i == 20) errorstop("!!!ERROR in iteration_k2 (4), energy not found");
			i++;
		}
		printf("\n brho1 = %lf\n", test_part.brho);
		//first tune
		tune_calc_matrix(&test_part, &qx1, &qz1, &betax, &alphax, &betaz, &alphaz, 1.e-5, 1.e-5, 1.e-5, 1.e-5, latt, part_cross_latt, YES, NULL);
		
		
		*reference_part = test_part;
		
		//change energy
		test_part.brho *=1.2; 
		
		accept_flag = find_closed_orbite_xxp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), eps_clo, latt, NO);
		if(accept_flag != TRUE) errorstop("!!!ERROR in iteration_k2 (5), closed orbit not found");
		
		//second tune
		tune_calc_matrix(&test_part, &qx2, &qz2, &betax, &alphax, &betaz, &alphaz, 1.e-5, 1.e-5, 1.e-5, 1.e-5, latt, part_cross_latt, YES, NULL);
		dqx = qx1-qx2;
		dqz = qz1-qz2;
		printf("dqx = %le\n", dqx);
		printf("dqz = %le\n", dqz);
		fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n", latt->cell[0].mpara[0][4], reference_part->brho, test_part.brho, reference_part->x, test_part.x, qx1, qz1, qx2, qz2, dqx, dqz);
		
		
	} while(latt->cell[0].mpara[0][4]<k2max);
	
	fprintf(wfile,"\n");
	fclose(wfile);
	
}

//works only for one cell lattice straight scaling FFAG
extern int noangleatend(struct Lattice *latt, struct Particle *reference_part, int print_option) 
{
	int i;
	double b0origin; 
	struct Particle test_part;
	
	b0origin = latt->cell[0].mpara[2][2];
	
	for (i = 0; i < 1000; i++) {
		test_part = *reference_part;
		test_part.hat = -2;	// no acceleration, no output
		part_cross_latt(&test_part, latt,NULL);
		if (fabs(test_part.ux) < 1.e-9) {
			if(print_option == YES) printf("\n\n\t\t\tb0adjust = %lf [T]\n\n\n", latt->cell[0].mpara[2][2]);
			return TRUE;
		}
		latt->cell[0].mpara[2][2] -= (test_part.ux)*10;
		if(print_option == YES) printf("ux=%le, b0=%lf\n", test_part.ux, latt->cell[0].mpara[2][2]);
	}
	if(print_option == YES) printf("\n \nin noangleatend precision not achieved: ux = %le, increase the number of turns or decrease demanded precision\n \n", fabs(test_part.ux));
	if(print_option == YES) printf("\n\n\t\t\tb0adjust = %lf [T]\n\n\n", latt->cell[0].mpara[2][2]);

	latt->cell[0].mpara[2][2] = b0origin;
	return FALSE;
}

//works only for straight scaling FFAG, with a 1-cell lattice (matching cell, starting in the middle of the magnet)
extern void scan_alpha_zero(char *outfilename, struct Lattice *latt, struct Particle *reference_part)
{
	int i,j,k,accept_flag;
	double betax,alphax,betax0,alphax0,betaz,alphaz,betaz0,alphaz0, magic_number, length_step, middle_step, bd_step, store_bf0, store_bd0, store_ymax, store_centerd;
	struct Particle test_part;
	FILE *wfile;
	wfile = fopen(outfilename,"w");
	fprintf(wfile, "noangle	betax\t	alphax\t	betaz\t	alphaz\t	Bf0\t	Bd0\t	Centerd\t	cell_length\n");
	
	store_bf0 = latt->cell[0].mpara[1][2];
	store_bd0 = latt->cell[0].mpara[2][2];
	store_ymax = latt->cell[0].boun.ymax;
	store_centerd = latt->cell[0].mpara[2][0];
	
	get_periodic_twiss(&betax0, &alphax0, &betaz0, &alphaz0, 1.e-4, 1.e-4, 1.e-4, 1.e-4, reference_part, latt,0);
	
	magic_number = 0.01;
	length_step = 0.01;
	middle_step = 0.01;
	bd_step = 0.001;
	
	latt->cell[0].efben[2][0] = -0.1;
	latt->cell[0].efbex[2][0] = 0.1;
	//for linear fringe field only:
	//latt->cell[0].efben[2][1] = 0.099;
	//latt->cell[0].efbex[2][1] = 0.099;
	printf("\n");
	for(i=0;i<12;i++) {
		latt->cell[0].boun.ymax = 4.9 + i*length_step; //change the cell length
		for(j=0;j<15;j++) {
			latt->cell[0].mpara[2][0] = 3.65 + j*middle_step; //change the position of half D magnet
			if(latt->cell[0].mpara[2][0] + 1.2 < latt->cell[0].boun.ymax) { //check that the cell length match the position of the half D magnet
				CLRSCR();
				//printf("%c[2K\r", 27);
				printf("i = %i, j = %i\t",i,j);
				fflush(stdout);
				for(k=0;k<40;k++) {
					//printf("k = %i\n",k);
					latt->cell[0].mpara[1][2] = -0.78 + k*bd_step; //change the Bf0 field
					test_part = *reference_part;
					test_part.hat = -2;
					accept_flag = noangleatend(latt, &test_part, NO);
					if(accept_flag != TRUE) printf("accept_flag!=TRUE\n");
					else {
						get_twissx_atinstru(&betax,&alphax,betax0,alphax0,1.e-4,1.e-4,&test_part,latt);
						get_twissz_atinstru(&betaz,&alphaz,betaz0,alphaz0,1.e-4,1.e-4,&test_part,latt);
						if(fabs(alphax) < magic_number && fabs(alphaz) < magic_number) fprintf(wfile, "%i	%lf	%lf	%lf	%lf	%le	%le	%lf	%lf\n", accept_flag,betax,alphax,betaz,alphaz, latt->cell[0].mpara[1][2], latt->cell[0].mpara[2][2], latt->cell[0].mpara[2][0], latt->cell[0].boun.ymax);
					}
				}
			}
		}
	}
	printf("\n\n");
	fprintf(wfile,"\n");
	fclose(wfile);
	
	latt->cell[0].mpara[1][2] = store_bf0;
	latt->cell[0].mpara[2][2] = store_bd0;
	latt->cell[0].boun.ymax = store_ymax;
	latt->cell[0].mpara[2][0] = store_centerd;
	latt->cell[0].efben[2][0] = -0.2;
	latt->cell[0].efbex[2][0] = 0.2;
}

extern void scan_alpha_zero2(char *outfilename, struct Lattice *latt, struct Particle *reference_part)
{
	int i,j,k,accept_flag;
	double betax,alphax,betax0,alphax0,betaz,alphaz,betaz0,alphaz0, magic_number, length_step, middle_step, bd_step, store_bf0, store_bd0, store_ymax, store_centerd;
	struct Particle test_part;
	struct Lattice rev_latt;
	FILE *wfile;
	wfile = fopen(outfilename,"w");
	fprintf(wfile, "noangle	betax\t	alphax\t	betaz\t	alphaz\t	Bf0\t	Bd0\t	Centerd\t	cell_length\n");
	
	store_bf0 = latt->cell[0].mpara[1][2];
	store_bd0 = latt->cell[0].mpara[2][2];
	store_ymax = latt->cell[0].boun.ymax;
	store_centerd = latt->cell[0].mpara[2][0];
	
	get_periodic_twiss(&betax0, &alphax0, &betaz0, &alphaz0, 1.e-4, 1.e-4, 1.e-4, 1.e-4, reference_part, latt,0);
	printf("bx=%lf, ax=%lf, bz=%lf, az=%lf\n", betax0, alphax0, betaz0, alphaz0);
	//betax0 = 18.304;
	//alphax0 = 0.;
	//betaz0 = 21.781;
	//alphaz0 = 0.;
	magic_number = 0.01;
	length_step = 0.02;
	middle_step = 0.02;
	bd_step = 0.01;
	
	latt->cell[0].efben[2][0] = -0.1;
	latt->cell[0].efbex[2][0] = 0.1;
	//for linear fringe field only:
	//latt->cell[0].efben[2][1] = 0.099;
	//latt->cell[0].efbex[2][1] = 0.099;
	
	for(i=0;i<20;i++) {
		printf("\ni = %i\n",i);
		latt->cell[0].boun.ymax = 6.2 + i*length_step; //change the cell length
		for(j=0;j<30;j++) {
			latt->cell[0].mpara[0][0] = 1.2 + j*middle_step; //change the position of half D magnet
			if(latt->cell[0].mpara[0][0] < latt->cell[0].mpara[1][0] - 0.6) { //check that the cell length match the position of the half D magnet
				printf("j = %i\t",j);
				fflush(stdout);
				for(k=0;k<15;k++) {
					//printf("k = %i\n",k);
					latt->cell[0].mpara[1][2] = -4.5 + k*bd_step; //change the Bf0 field
					test_part = *reference_part;
					test_part.hat = -2;
					
					reverse_latt(&rev_latt, latt);
					accept_flag = noangleatend(&rev_latt, &test_part, NO); //adjust the half D magnet B0
					//if(accept_flag != TRUE) errorstop("!!!ERROR in scan_alpha_zero, no-angle closed orbit not found\n");
					if(accept_flag != TRUE) printf("accept_flag!=TRUE\n");
					else {
						//printf("in!\n");
						latt->cell[0].mpara[0][2] = rev_latt.cell[0].mpara[2][2];
						get_twissx_atinstru(&betax,&alphax,betax0,alphax0,1.e-4, 1.e-4,&test_part,latt);
						get_twissz_atinstru(&betaz,&alphaz,betaz0,alphaz0,1.e-4,1.e-4,&test_part,latt);
						if(fabs(alphax) < magic_number && fabs(alphaz) < magic_number) {
							//printf("x=%lf,z=%lf\n", alphax, alphaz);
							//printf("OK!\n");
							fprintf(wfile, "%i	%lf	%lf	%lf	%lf	%le	%le	%lf	%lf\n", accept_flag,betax,alphax,betaz,alphaz,latt->cell[0].mpara[1][2],rev_latt.cell[0].mpara[0][2],latt->cell[0].mpara[2][0],latt->cell[0].boun.ymax);
						}
					}
				}
			}
		}
	}//
	printf("\n\n");
	fprintf(wfile,"\n");
	fclose(wfile);
	free_latt(&(rev_latt));
	
	latt->cell[0].mpara[1][2] = store_bf0;
	latt->cell[0].mpara[2][2] = store_bd0;
	latt->cell[0].boun.ymax = store_ymax;
	latt->cell[0].mpara[2][0] = store_centerd;
	latt->cell[0].efben[2][0] = -0.2;
	latt->cell[0].efbex[2][0] = 0.2;
}

extern int adjust_b_maxangle(struct Lattice *latt, struct Particle *reference_part, double auth_angle, double precision, int print_option)
{
	int i;
	double b0origin, qx, qz, betax, alphax, betaz, alphaz;
	struct Particle part;
	
	b0origin = latt->cell[0].mpara[1][2];
	
	
	for (i = 0; i < 100; i++) {
		max_angle = 0.;
		part = *reference_part;
		part.hat = -2;	// no acceleration, no output
		if (find_closed_orbite_xxp(&part, &(part.x), &(part.ux), &(part.uy), 1.e-8, latt, NO) == TRUE) {
			part.hat = 3;
			part_cross_latt(&part,latt,NULL);
			part.hat = -2;
			if(max_angle < auth_angle + precision && max_angle> auth_angle - precision) {
				if(print_option == YES)printf("\n\n\t\t\tb0adjust = %lf [T]\n\t\t\tmax_angle = %lf [mrad]\n\n", latt->cell[0].mpara[1][2], max_angle*1000.);
				//*reference_part = part;
				tune_calc_matrix(&part, &qx, &qz, &betax, &alphax, &betaz, &alphaz, 1.e-4, 1.e-4, 1.e-4, 1.e-4, latt, part_cross_latt, YES, NULL);
				return TRUE;
			}
			latt->cell[0].mpara[1][2] -= (max_angle - auth_angle)*2.;
			if(print_option == YES) {
				CLRSCR();
				printf("#%i, angle = %le [rad], b0 = %le", i, max_angle, latt->cell[0].mpara[1][2]);
				fflush(stdout);
			}
		}
		else {
			printf("\n \nproblem in adjust_b_maxangle, closed orbit not found!\n");
			//latt->cell[0].mpara[1][2] = b0origin;
			return FALSE;
		}
	}
	if(print_option == YES) printf("\n \nin noangleatend precision not achieved: ux = %le, increase the number of turns or decrease demanded precision\n \n", fabs(part.ux));
	if(print_option == YES) printf("\n\n\t\t\tb0adjust = %lf [T]\n\t\tmax_angle = %lf [mrad]\n\n", latt->cell[0].mpara[2][2], max_angle*1000.);
	return FALSE;
}

extern void gene_ellibunch_twiss_x(struct Beam *beam, double emitx, double betax, double alphax, int nbparts)
//it is assumed here that beam->part[0] has already been initialized, and this will be taken as the center of the bunch of particle that is going to be loaded in beam.
{
	int i, n;
	double t, p, x_ctr, xprime_ctr, z,
		x, xprime, y, zprime, px, py, pz,
		tx;
	struct Particle central_part;
	
	//test errors and warnings
	if(nbparts <= 0) errorstop("!!!ERROR in gene_ellibunch:\n nsteps[1] <= 0");
	if(beam->part[0].uy == 0) errorstop("!!!ERROR in gene_ellibunch: beam->part[0].uy must not be 0!!!\n");
	
	//save beam->part[0] in central part and free the memory already alocated for beam
	central_part = beam->part[0];
	free_beam(beam);
	alloc_beam(beam, nbparts+1);
	beam->part[0] = central_part;
	
	//compute central particle clue parameters (total momentum, x, xprime, z, zprime)
	t = beam->part[0].t;
	p = fabs(beam->part[0].brho*beam->part[0].q);
	x_ctr = beam->part[0].x;
	xprime_ctr = atan(beam->part[0].ux/beam->part[0].uy);
	z = beam->part[0].z;
	zprime = atan(beam->part[0].uz/sqrt(pow(beam->part[0].ux, 2) + pow(beam->part[0].uy, 2)));
	
	//generate a 2D beam, along hor. ellipse centered on central particle
	for(i = 0; i < nbparts; i++) {
		tx = TWOPI*i/(nbparts);
		
		x = x_ctr + cos(tx)*sqrt(betax*emitx);
		xprime	= xprime_ctr + (sin(tx)*sqrt(betax*emitx)-alphax*(x-x_ctr))/betax;
		y = 0;
		px = p*sin(xprime)*cos(zprime);
		py = p*cos(xprime)*cos(zprime);
		pz = p*sin(zprime);
		
		n = 1 + i;
		
		gene_singlepart(&(beam->part[n]),beam->part[0].m0,beam->part[0].q, x, y, z, px, py, pz, t, NO);
	}
}

extern void gene_ellibunch_twiss_z(struct Beam *beam, double emitz, double betaz, double alphaz, int nbparts)
//it is assumed here that beam->part[0] has already been initialized, and this will be taken as the center of the bunch of particle that is going to be loaded in beam.
{
	int i, n;
	double t, p, z_ctr, zprime_ctr, z,
		x, xprime, y, zprime, px, py, pz,
		tz;
	struct Particle central_part;
	
	//test errors and warnings
	if(nbparts <= 0) errorstop("!!!ERROR in gene_ellibunch:\n nsteps[1] <= 0");
	if(beam->part[0].uy == 0) errorstop("!!!ERROR in gene_ellibunch: beam->part[0].uy must not be 0!!!\n");
	
	//save beam->part[0] in central part and free the memory already alocated for beam
	central_part = beam->part[0];
	free_beam(beam);
	alloc_beam(beam, nbparts+1);
	beam->part[0] = central_part;
	
	//compute central particle clue parameters (total momentum, x, xprime, z, zprime)
	t = beam->part[0].t;
	p = fabs(beam->part[0].brho*beam->part[0].q);
	x = beam->part[0].x;
	xprime = atan(beam->part[0].ux/beam->part[0].uy);
	z_ctr = beam->part[0].z;
	zprime_ctr = atan(beam->part[0].uz/sqrt(pow(beam->part[0].ux, 2) + pow(beam->part[0].uy, 2)));
	
	//generate a 2D beam, along vert. ellipse centered on central particle
	for(i = 0; i < nbparts; i++) {
		tz = TWOPI*i/(nbparts);
		z = z_ctr + cos(tz)*sqrt(betaz*emitz);
		zprime	= zprime_ctr + (sin(tz)*sqrt(betaz*emitz)-alphaz*(z-z_ctr))/betaz;
		
		y = 0;
		px = p*sin(xprime)*cos(zprime);
		py = p*cos(xprime)*cos(zprime);
		pz = p*sin(zprime);
		n = 1 + i;
		gene_singlepart(&(beam->part[n]),beam->part[0].m0,beam->part[0].q, x, y, z, px, py, pz, t, NO);
	}
}

// emit_x=rms_emit_x
extern void gene_4D_beam_gasdev(struct Beam *beam, double emit_x, double twiss_betax, double twiss_alphax, double emit_z, double twiss_betaz, double twiss_alphaz, double factor_rms_tot, int nparts, long seed)
{
	int i;
	long idum;
	double t_ctr, p_ctr, x_ctr, xprime_ctr, z_ctr, zprime_ctr,
	t, p, x, xprime, y, z, zprime, px, py, pz,
	xn, xpn, zn, zpn;
	struct Particle central_part;
	
	if(debug == YES) printf("gene_4D_beam_gasdev: nparts = %i\n", nparts);
	
	//test errors and warnings
	if(nparts <= 0) errorstop("!!!ERROR in gene_4D_beam_gasdev:\n nparts <= 0");
	if(beam->part[0].uy == 0) errorstop("!!!ERROR in gene_4D_beam_gasdev: beam->part[0].uy = 0!!!\n");
	if(twiss_betax <= 0) errorstop("!!!ERROR in gene_4D_beam_gasdev: twiss_betax <= 0");
	if(twiss_betaz <= 0) errorstop("!!!ERROR in gene_4D_beam_gasdev: twiss_betaz <= 0");
	
	//save beam->part[0] in central_part and free the memory already allocated for beam, then allocate new memory
	central_part = beam->part[0];
	free_beam(beam);
	alloc_beam(beam, nparts+1);
	beam->part[0] = central_part;
	
	//compute central particle clue parameters (total momentum, x, xprime, z, zprime)
	t_ctr = beam->part[0].t;
	p_ctr = fabs(beam->part[0].brho*beam->part[0].q);
	x_ctr = beam->part[0].x;
	xprime_ctr = atan(beam->part[0].ux/beam->part[0].uy);
	z_ctr = beam->part[0].z;
	zprime_ctr = atan(beam->part[0].uz/sqrt(pow(beam->part[0].ux, 2) + pow(beam->part[0].uy, 2)));
	
	//initialise idum (for the random generator)
	if(seed >=0) idum = -seed;
	else idum = seed;
	
	
	for(i = 1; i < nparts+1 ; i++) {
		//gaussian normalized transverse coordinates
		for(;;) {
			xn = gasdev(&idum);
			xpn = gasdev(&idum);
			if(fabs(xn) < factor_rms_tot && fabs(xpn) < factor_rms_tot) break;
		}
		for(;;) {
			zn = gasdev(&idum);
			zpn = gasdev(&idum);
			if(fabs(zn) < factor_rms_tot && fabs(zpn) < factor_rms_tot) break;
		}
		
		//transform normalized transverse coordinates (xn,xpn,zn,zpn) into real transverse coordinates (x,xp,z,zp)
		x = xn*sqrt(emit_x*twiss_betax);
		xprime	= (xpn*sqrt(emit_x*twiss_betax)-twiss_alphax*(x))/twiss_betax;
		z = zn*sqrt(emit_z*twiss_betaz);
		zprime	= (zpn*sqrt(emit_z*twiss_betaz)-twiss_alphaz*(z))/twiss_betaz;
		
		t = t_ctr;
		p = p_ctr;
		
		x += x_ctr;
		xprime	+= xprime_ctr;
		
		z += z_ctr;
		zprime	+= zprime_ctr;
		
		y = 0;
		px = p*sin(xprime)*cos(zprime);
		py = p*cos(xprime)*cos(zprime);
		pz = p*sin(zprime);
		
		//generate a particle at (x,y,z) with momentum (px,py,pz) and time t.
		gene_singlepart(&(beam->part[i]),beam->part[0].m0,beam->part[0].q, x, y, z, px, py, pz, t, NO);
	}// and it is a loop, it does it again nparts times...
}

//compute average radius from trackout (check trackout output form before use!!)
extern void compute_av_radius(char *trackout, char *outfilename)
{
	int i;
	double s, x, y, z, bx, by, bz, ux, uy, uz, brho, t, r, rtot, brhoini;
	FILE *wfile;
	wfile = fopen(outfilename,"w");
	FILE *rfile = NULL;
	rfile = fopen(trackout, "r");
	if(rfile == NULL) errorstop("!!!ERROR in compute_av_radius: cannot open inputfile!!!");
	brhoini = 0;
	i=0;
	for(;;) {
		fscanf(rfile, "%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le", &s, &x, &y, &z, &bx, &by, &bz, &ux, &uy, &uz, &brho, &t, &r);
		//printf("%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le\n", s, x, y, z, bx, by, bz, ux, uy, uz, brho, t, r);
		if(i>2000) {
			fprintf(wfile, "%le  %le\n", brhoini, rtot);
			break;
		}
		if(fabs(brho-brhoini)> TINYDIMLESS) {
			rtot = rtot/((double)i);
			fprintf(wfile, "%le  %le\n", brhoini, rtot);
			printf("%le  %le", brhoini, rtot);
			rtot = 0;
			i=0;
		}
		else printf("s=%le\n", s);
		rtot += r; 
		i++;
		brhoini = brho;
		newline(rfile);
	}
	fclose(rfile);
	fclose(wfile);
}

//ith derivative of simple Enge denominator u=(1+ee)*(1+es) (with only c1!=0) with respect to theta
extern double deriv_u_dt(int order, double rlce, double rlcs, double ee, double es) {
	return (ee*pow(-rlce,order)+es*pow(rlcs,order)+ee*es*pow((rlcs-rlce),order));
}

//ith derivative of simple Enge denominator u=(1+ee)*(1+es) (with only c1!=0) with respect to r
extern double deriv_u_dr(int order, double r, double ee, double es, double rlcte, double rlcts) {
	double a,b,c;
	a=rec_mul_minus(rlcte,order-1);
	//b=pow(-1.0,order)*rec_mul_plus(rlcts,order-1);
	b=rec_mul_minus(-rlcts,order-1);
	c=rec_mul_minus((rlcte-rlcts),order-1);
	return (a*ee+b*es+c*ee*es)/(pow(r,order));
}

// parameters of iterative derivative of Enge denominator u=(1+ee)*(1+es) (with only c1!=0) with respect to theta (diudti=ee*ai+es*bi+ee*es*ci)
extern void derudt(double rlce, double rlcs, double ai1, double bi1, double ci1, double *ai, double *bi, double *ci) {
	*ai = -rlce*ai1;
	*bi = rlcs*bi1;
	*ci = (rlcs-rlce)*ci1;
}

// parameters of iterative derivative of Enge denominator u=(1+ee)*(1+es) (with only c1!=0) with respect to r (diudri=(ee*ai+es*bi+ee*es*ci)/pow(r,i))
extern void derudr(int order, double rlcte, double rlcts, double ai1, double bi1, double ci1, double *ai, double *bi, double *ci)
{
	*ai = ai1*(rlcte-order+1);
	*bi = bi1*(-rlcts-order+1);
	*ci = ci1*(rlcte-rlcts-order+1);
}

extern void compute_field_spi_enge(struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), 
	double r0, double rmax, int nbr, double th0, double thmax, int nbth, double z0, double zmax, int nbz)
{
	int i,j,k;
	double bx,by,bz,x,y,z,r,th, zstep,rstep,thstep;
	FILE *wfile;
	
	wfile = fopen("data/track.dat","w");
	
	write_enge = YES;
	//zstep = (zmax-z0)/nbz;
	zstep = comp_step(z0, zmax, nbz);
	//rstep = (rmax-r0)/nbr;
	rstep = comp_step(r0, rmax, nbr);
	//thstep = (thmax-th0)/nbth;
	thstep = comp_step(th0, thmax, nbth);
	
	if(nbz>1) nbz+=1; 
	if(nbr>1) nbr+=1; 
	if(nbth>1) nbth+=1; 
	
	for(i=0;i<nbz;i++) {
		printf("%i	%i\n",i, nbz);
		z = z0 + i*zstep;
		for(j=0;j<nbr;j++) {
			r = r0 + j*rstep;
			for(k=0;k<nbth;k++) {
				th = k*thstep;
				x = r*cos(th);
				y = r*sin(th);
				if(get_bfield(x, y, z, &bx, &by, &bz, cell,add_contribution_comp)!=ALIVE) errorstop("bfield!=ALIVE");
				fprintf(wfile,"%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n", r, th, x, y, z, bx, by, bz);
			}
		}
	}
	write_enge = NO;
	fclose(wfile);
}

extern void compute_fieldmap(char *filename, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int),
							double r0, double rmax, int nbr, double th0, double thmax, int nbth, double z0, double zmax, int nbz)
{
	int i,j,k, circ=YES;
	double x,y,z,bx,by,bz,br,bth,r,th, rstep,thstep,zstep;
	
	FILE *wfile;

	wfile = fopen(filename,"w");
	fprintf(wfile,"%i %i %i\n",nbr,nbth,nbz);
	fprintf(wfile,"(%i, %i, %i) nodes\n",nbr,nbth,nbz);
	fprintf(wfile,"from (%le, %le, %le) to (%le, %le, %le)\n",r0, th0, z0, rmax, thmax*180./PI, zmax);
	
	//if(nbz>1) zstep = (zmax-z0)/(nbz-1);
	//else zstep = 0.;
	zstep = comp_step(z0, zmax, nbz);
	//if(nbr>1) rstep = (rmax-r0)/(nbr-1);
	//else rstep = 0.;
	rstep = comp_step(r0, rmax, nbr);
	//if(nbth>1) thstep = (thmax-th0)/(nbth-1);
	//else thstep = 0.;
	thstep = comp_step(th0, thmax, nbth);
	fprintf(wfile,"cylindrical grid with step (%le [m], %le [deg], %le [m])\n", rstep,thstep*180./PI,zstep);
	fprintf(wfile, "r [m], theta [deg], z [m], Br [T], Btheta, Bz\n");
	
	printf("rstep=%le, thstep=%le, zstep = %le\n", rstep,thstep,zstep);
	
	for(i=0;i<nbz;i++) {
		z = z0 + i*zstep;
		for(j=0;j<nbr;j++) {
			r = r0 + j*rstep;
			for(k=0;k<nbth;k++) {
				th = k*thstep;
				if(circ == YES) {
					x = r*cos(th);
					y = r*sin(th);
				}
				else {
					x = r;
					y = th;
				}
				//printf("%le, %le, %le\n", x, y, z);
				if(get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp)!=ALIVE) errorstop("bfield not alive");
				//printf("%lf, %lf, %lf	%le, %le, %le\n", x, y, z, bx, by, bz);
				
				//x = r*cos(th-thmax/2.);
				//y = r*sin(th-thmax/2.);
				//fprintf(wfile,"%lf	%lf	%lf	%lf	%lf	%le	%le	%le\n", x, y, z, r, th*180./PI, bx, by, bz);
				//fprintf(wfile,"%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", x, y, z, r, th, bx, by, bz);
				//fprintf(wfile,"%lf	%lf	%lf	%lf	%lf	%lf\n", x*100, z*100, y*100, bx*10000, bz*10000, by*10000);
				br = bx*cos(th) + by*sin(th);
				bth = by*cos(th) - bx*sin(th);
				fprintf(wfile,"%lf	%lf	%lf	%lf	%lf	%lf\n", r, th*180./PI, z, br, bth, bz);
				//fprintf(wfile,"%lf	%lf	%lf	%lf	%lf	%lf\n", r, th*180./PI, z, bx, by, bz);
			}
		}
	}
	fclose(wfile);
}

extern void analysis_field_spi_enge(char *rfilename, char *outfilename, int order, int nb_elements)
{
	int i,ic,max_b;
	double b[8], b_sum[8];
	double r,t,z;
	FILE *rfile = NULL;
	FILE *wfile;
	
	rfile = fopen(rfilename, "r");
	if(rfile==NULL) errorstop("cannot open rfile!");
	wfile = fopen(outfilename, "w");
	b_sum[0]=0.;
	for(i=0;i<8;i++) b_sum[i]=0.0;
	max_b = floor(order/2.0)+2;
	while(!feof(rfile)) {
		if(order == 12) fscanf(rfile,"%i\t%le\t%le\t%le\t %le\t%le\t%le\t%le\t%le\t%le\t%le\t%le",&ic,&r,&t,&z,&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7],&b[0]);
		else if(order == 11 || order == 10) fscanf(rfile,"%i\t%le\t%le\t%le\t %le\t%le\t%le\t%le\t%le\t%le\t%le",&ic,&r,&t,&z,&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[0]);
		else if(order == 9 || order == 8) fscanf(rfile,"%i\t%le\t%le\t%le\t %le\t%le\t%le\t%le\t%le\t%le",&ic,&r,&t,&z,&b[1],&b[2],&b[3],&b[4],&b[5],&b[0]);
		else if(order == 7 || order == 6) fscanf(rfile,"%i\t%le\t%le\t%le\t %le\t%le\t%le\t%le\t%le",&ic,&r,&t,&z,&b[1],&b[2],&b[3],&b[4],&b[0]);
		else if(order == 5 || order == 4) fscanf(rfile,"%i\t%le\t%le\t%le\t %le\t%le\t%le\t%le",&ic,&r,&t,&z,&b[1],&b[2],&b[3],&b[0]);
		else if(order == 3 || order == 2) fscanf(rfile,"%i\t%le\t%le\t%le\t %le\t%le\t%le",&ic,&r,&t,&z,&b[1],&b[2],&b[0]);
		else if(order == 1) fscanf(rfile,"%i\t%le\t%le\t%le\t %le\t%le",&ic,&r,&t,&z,&b[1],&b[0]);
		//printf("%i\t%le\t%le\t%le\t %le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",ic,r,t,z,b[1],b[2],b[3],b[4],b[5],b[6],b[7],b[0]);
		if(ic==0) for(i=0;i<max_b;i++) b_sum[i]=0.0;
		
		for(i=0;i<max_b;i++) b_sum[i]+=b[i];
		//printf("%i\t%le\t%le\t%le\t %le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",ic,r,t,z,b_sum[1],b_sum[2],b_sum[3],b_sum[4],b_sum[5],b_sum[6],b_sum[7],b_sum[0]);
		//printf("%i\t%le\t%le\t%le\t",ic,r,t,z);
		//for(i=1;i<max_b;i++) printf("%le\t",b[i]);
		//printf("%le\n",b[0]);
		if(ic==nb_elements-1) {
			fprintf(wfile,"%le\t%le\t%le",r,t,z);
			//printf("%i\t%le\t%le\t%le",ic,r,t,z);
			for(i=0;i<max_b;i++) {
				fprintf(wfile,"\t%le",b_sum[i]);
				//printf("\t%le",b_sum[i]);
			}
			fprintf(wfile,"\n");
			//printf("\n");
		}
	}
	fclose(rfile);
	fclose(wfile);
}

extern int maxwell_test(double x, double y, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int),
	double *divb, double *curlbx, double *curlby, double *curlbz, int doyouglob)//, double *dbzdx, double *dbxdz, double *dbzdy, double *dbydz, double *dbxdy, double *dbydx)
{
	double bx, by, bz, bxpdx, bxpdy, bxpdz, bypdy, bypdx, bypdz, bzpdz, bzpdx, bzpdy, dglob_to_loc_x, dglob_to_loc_y;
	double bxmdx, bxmdy, bxmdz, bymdy, bymdx, bymdz, bzmdz, bzmdx, bzmdy;
	double dbxdx,dbydy,dbzdz,dbzdx,dbxdz,dbzdy,dbydz,dbxdy,dbydx;
	
	get_bfield(x, y, z, &bx, &by, &bz, cell,add_contribution_comp);
	//printf("bx=%le, by=%le, bz=%le\n",bx, by, bz);
	dglob_to_loc_x = dx;
	dglob_to_loc_y = 0.;
	if(doyouglob==YES) fwk_vect_globtoloc(&dglob_to_loc_x, &dglob_to_loc_y, &(cell->framework));
	get_bfield(x+dglob_to_loc_x, y+dglob_to_loc_y, z, &bxpdx, &bypdx, &bzpdx, cell, add_contribution_comp);
	get_bfield(x-dglob_to_loc_x, y-dglob_to_loc_y, z, &bxmdx, &bymdx, &bzmdx, cell, add_contribution_comp);
	//printf("bxpdx=%le, bypdx=%le, bzpdx=%le\n",bxpdx, bypdx, bzpdx);
	dglob_to_loc_x = 0;
	dglob_to_loc_y = dy;
	if(doyouglob==YES) fwk_vect_globtoloc(&dglob_to_loc_x, &dglob_to_loc_y, &(cell->framework));
	get_bfield(x+dglob_to_loc_x, y+dglob_to_loc_y, z, &bxpdy, &bypdy, &bzpdy, cell, add_contribution_comp);
	get_bfield(x-dglob_to_loc_x, y-dglob_to_loc_y, z, &bxmdy, &bymdy, &bzmdy, cell, add_contribution_comp);
	//printf("bxpdy=%le, bypdy=%le, bzpdy=%le\n",bxpdy, bypdy, bzpdy);
	get_bfield(x, y, z+dz, &bxpdz, &bypdz, &bzpdz, cell, add_contribution_comp);
	get_bfield(x, y, z-dz, &bxmdz, &bymdz, &bzmdz, cell, add_contribution_comp);
	//printf("bxpdz=%le, bypdz=%le, bzpdz=%le\n",bxpdz, bypdz, bzpdz);
	
	//*divb = (bxpdx-bx)/dx + (bypdy-by)/dy + (bzpdz-bz)/dz;
	//*curlbx = (bzpdx - bz)/dx - (bxpdz - bx)/dz;
	//*curlby = (bzpdy - bz)/dy - (bypdz - by)/dz;
	//*curlbz = (bxpdy - bx)/dy - (bypdx - by)/dx;
	//dbxdx = (bxpdx-bx)/dx;
	//dbydy = (bypdy-by)/dy;
	//dbzdz = (bzpdz-bz)/dz;
	dbxdx = (bxpdx-bxmdx)/(2*dx); //second order cancelled
	dbydy = (bypdy-bymdy)/(2*dy);
	dbzdz = (bzpdz-bzmdz)/(2*dz);
	*divb = dbxdx + dbydy + dbzdz;
	//printf("%le + %le + %le = %le\n", dbxdx, dbydy, dbzdz, *divb);
	//dbzdx = (bzpdx - bz)/dx;
	//dbxdz = (bxpdz - bx)/dz;
	dbzdx = (bzpdx - bzmdx)/(2*dx);
	dbxdz = (bxpdz - bxmdz)/(2*dz);
	*curlbx = dbzdx - dbxdz;
	//printf("%le - %le = %le\n", dbzdx, dbxdz, *curlbx);
	//dbzdy = (bzpdy - bz)/dy;
	//dbydz = (bypdz - by)/dz;
	dbzdy = (bzpdy - bzmdy)/(2*dy);
	dbydz = (bypdz - bymdz)/(2*dz);
	*curlby = dbzdy - dbydz;
	//printf("%le - %le = %le\n", dbzdy, dbydz, *curlby);
	//dbxdy = (bxpdy - bx)/dy;
	//dbydx = (bypdx - by)/dx;
	dbxdy = (bxpdy - bxmdy)/(2*dy);
	dbydx = (bypdx - bymdx)/(2*dx);
	*curlbz = dbxdy - dbydx;
	//printf("%le - %le = %le\n", dbxdy, dbydx, *curlbz);
	return TRUE;
}

extern void maxwell_test_txtfile(char *textfile, double x, double y, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	double divb, curlbx, curlby, curlbz;//, dbzdx, dbxdz, dbzdy, dbydz, dbxdy, dbydx;
	double r=sqrt(x*x+y*y),th=atan_ratio(y, x);
	FILE *wfile;
	
	wfile = fopen(textfile, "a");
	if(maxwell_test(x, y, z, dx, dy, dz, cell, add_contribution_comp, &divb, &curlbx, &curlby, &curlbz,NO)!=TRUE) {
		divb = 0;
		curlbx = 0;
		curlby = 0;
		curlbz = 0;
	}
	//, &dbzdx, &dbxdz, &dbzdy, &dbydz, &dbxdy, &dbydx);
	//fprintf(wfile, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",r,th,z, divb, curlbx, dbzdx, dbxdz, curlby, dbydz, dbzdy, curlbz, dbydx, dbxdy); 
	if(cell->boun.thmax != 0) fprintf(wfile, "%le\t%le\t%le\t%le\t%le\t%le\t%le\n",r,th,z, divb, curlbx, curlby, curlbz); 
	else fprintf(wfile, "%le\t%le\t%le\t%le\t%le\t%le\t%le\n",x,y,z, divb, curlbx, curlby, curlbz); 
	fclose(wfile);
}

extern void maxwell_test_str_txtfile(char *textfile, double x, double y, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	double divb, curlbx, curlby, curlbz;//, dbzdx, dbxdz, dbzdy, dbydz, dbxdy, dbydx;
	FILE *wfile;
	
	wfile = fopen(textfile, "a");
	if(maxwell_test(x, y, z, dx, dy, dz, cell, add_contribution_comp, &divb, &curlbx, &curlby, &curlbz,NO)!=TRUE) { //, &dbzdx, &dbxdz, &dbzdy, &dbydz, &dbxdy, &dbydx);
		divb = 0;
		curlbx = 0;
		curlby = 0;
		curlbz = 0;
	}
	
	//fprintf(wfile, "%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\t%le\n",r,th,z, divb, curlbx, dbzdx, dbxdz, curlby, dbydz, dbzdy, curlbz, dbydx, dbxdy); 
	fprintf(wfile, "%le\t%le\t%le\t%le\t%le\t%le\t%le\n",x,y,z, divb, curlbx, curlby, curlbz); 
	fclose(wfile);
}

extern void maxwell_test_dxeffect(char *textfile, double x, double y, double z, double dxmin, double dxmax, int nbstep, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double dx, dxstep;
	
	emptyfile(textfile);
	
	dxstep = (dxmax-dxmin)/(nbstep);
	
	for(i=0;i<nbstep+1;i++) {
		dx = dxmin+dxstep*i;
		maxwell_test_txtfile(textfile, x, y, z, dx, dx, dx, cell, add_contribution_comp);	
	}
}

extern void maxwell_test_across_cell(char *textfile, double r, double z, double dx, int th_step, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double th, stepth,x,y;
	
	emptyfile(textfile);
	stepth = cell->boun.thmax/th_step;
	for(i=0;i<th_step;i++) {
		th = i*stepth;
		x = r*cos(th);
		y = r*sin(th);
		maxwell_test_txtfile(textfile, x, y, z, dx, dx, dx, cell, add_contribution_comp);
	}
}

extern double maxwell_test_max(double x, double y, double z, double dx, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), int doyouglob)
{
	double divb=0, curlbx=0, curlby=0, curlbz=0, outmax=0.0;
	
	if(maxwell_test(x, y, z, dx, dx, dx, cell, add_contribution_comp, &divb, &curlbx, &curlby, &curlbz, doyouglob)!=TRUE) {
		divb = 0;
		curlbx = 0;
		curlby = 0;
		curlbz = 0;
	}
	printf("div,curlx,y,z=(%le,%le,%le,%le)\n",divb, curlbx, curlby, curlbz);
	outmax = MAX(outmax,fabs(divb));
	outmax = MAX(outmax,fabs(curlbx));
	outmax = MAX(outmax,fabs(curlby));
	outmax = MAX(outmax,fabs(curlbz));
	
	return outmax;
}

extern double maxwell_test_across_cell_max(double r, double z, double dx, int th_step, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	//double divb, curlbx, curlby, curlbz;//, dbzdx, dbxdz, dbzdy, dbydz, dbxdy, dbydx;
	double th, stepth,x,y, test_max, outmax=0.0;
		
	stepth = cell->boun.thmax/th_step;
	
	for(i=0;i<th_step;i++) {
		th = i*stepth;
		x = r*cos(th);
		y = r*sin(th);
		test_max = maxwell_test_max(x, y, z, dx, cell, add_contribution_comp,NO);
		outmax = MAX(outmax,test_max);
	}
	return outmax;
}

extern void maxwell_test_map_radius(char *txt_prefix, double r0, double rmax, int r_step, double th0, double thmax, int th_step, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), int doyoucirc)
{
	char name[300];
	int i,j;
	double r, th, step_r, step_th, x, y, divb, curlbx, curlby, curlbz;
	FILE *wfilediv, *wfilex, *wfiley, *wfilez;
	
	sprintf(name,"%s_div.dat", txt_prefix);
	wfilediv = fopen(name, "w");
	sprintf(name,"%s_curlx.dat", txt_prefix);
	wfilex = fopen(name, "w");
	sprintf(name,"%s_curly.dat", txt_prefix);
	wfiley = fopen(name, "w");
	sprintf(name,"%s_curlz.dat", txt_prefix);
	wfilez = fopen(name, "w");
	
	step_r = comp_step(r0, rmax, r_step);
	step_th = comp_step(th0, thmax, th_step);
	printf("r0=%le, rmax=%le, r_step=%le\n", r0, rmax, step_r);
	printf("th0=%le, thmax=%le, th_step=%le\n", th0, thmax, step_th);
	
	
	for(j=0;j<th_step;j++) {
		th = th0 + j*step_th;
		for(i=0;i<r_step;i++) {
			r=r0 + i*step_r;
			if(doyoucirc == YES) {
				x = r*cos(th);
				y = r*sin(th);
			}
			else {
				x = r;
				y = th;
			}
			printf("j=%i,i=%i,x=%le,y=%le\n",j,i,x,y);
			if(maxwell_test(x, y, z, dx, dy, dz, cell, add_contribution_comp, &divb, &curlbx, &curlby, &curlbz, YES)!=TRUE) { //, &dbzdx, &dbxdz, &dbzdy, &dbydz, &dbxdy, &dbydx);
				divb = 0;
				curlbx = 0;
				curlby = 0;
				curlbz = 0;
			}
			
			//printf("loc (%le,%le)\t", x,y);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			//printf("glob (%le,%le)\n", x,y);
			fprintf(wfilediv, "%le	%le	%le\n", y, x, divb);
			fprintf(wfilex, "%le	%le	%le\n", y, x, curlbx);
			fprintf(wfiley, "%le	%le	%le\n", y, x, curlby);
			fprintf(wfilez, "%le	%le	%le\n", y, x, curlbz);
		}
		fprintf(wfilediv, "\n");
		fprintf(wfilex, "\n");
		fprintf(wfiley, "\n");
		fprintf(wfilez, "\n");
	}
	fclose(wfilediv);
	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
}

extern double convergence_radius(double dx, int th_step, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i, n, origin_order;
	double max1,max2,z, z0=0., z1=2*cell->collim.zmax;
	if(strcmp(cell->keyword, "ffag-spi-enge") != 0 && strcmp(cell->keyword, "ffag-spi-fullenge") != 0) errorstop("convergence_radius works only with ffag-spi-enge and ffag-spi-fullenge cells\n");
	
	origin_order = cell->efben[0][2];
	
	for(n=0;n<50;n++) {
		z=(z1+z0)/2.;
		for(i=0;i<cell->nbcomp;i++) cell->efben[i][2] = 2;
		max1 = maxwell_test_across_cell_max(cell->mpara[0][1], z, dx, th_step, cell, add_contribution_comp);
		for(i=0;i<cell->nbcomp;i++) cell->efben[i][2] = 12;
		max2 = maxwell_test_across_cell_max(cell->mpara[0][1], z, dx, th_step, cell, add_contribution_comp);
		if(max2 > max1) z1=z;
		else z0=z;
		if((z1-z0)<0.001) {
			for(i=0;i<cell->nbcomp;i++) cell->efben[i][2] = origin_order;
			if(fabs(z1-cell->collim.zmax)<TINYDIMLESS) return z1;
			else return z0;
		}
		printf("try %i, z0=%le,z1=%le,max1=%le,max2=%le\n",n,z0,z1,max1,max2);
	}
	printf("too many tries!\n");
	for(i=0;i<cell->nbcomp;i++) cell->efben[i][2] = origin_order;
	return 0;
}

extern void compute_convergence_limit_vffa(struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))//, int doyouglob)
{
	int i,l;
	double x_0,y_0,angle,r,x,y,xtest,xmax,xmin,dx=1.e-4,prec=1.e-5;
	//double max1,max2;
	double bx2,by2,bz2,bx10,by10,bz10;
	struct Cell temp_cell;
	
	temp_cell.deltar = 0; //default value
	temp_cell.nbcomp = 1; //default value
	temp_cell.instrutype = NO; //default value
	temp_cell.framework = cell->framework; 	//set cell.framework
	temp_cell.mpara = allocmatrix(temp_cell.nbcomp, SIZE_MPARA);
	temp_cell.efben = allocmatrix(temp_cell.nbcomp, SIZE_EFB);
	temp_cell.efbex = allocmatrix(temp_cell.nbcomp, SIZE_EFB);
	temp_cell.stepsize = cell->stepsize;
	temp_cell.deltar = cell->deltar;
	temp_cell.collim.rmin = cell->collim.rmin;
	temp_cell.collim.rmax = cell->collim.rmax;
	temp_cell.collim.zmin = cell->collim.zmin;
	temp_cell.collim.zmax = cell->collim.zmax;
	temp_cell.boun.thmax = cell->boun.thmax;
	temp_cell.boun.ymax = cell->boun.ymax;
	for(i=0;i<cell->nbcomp;i++) {
		for(l=0;l<SIZE_MPARA;l++) temp_cell.mpara[0][l] = cell->mpara[i][l];
		for(l=0;l<SIZE_EFB;l++) temp_cell.efben[0][l] = cell->efben[i][l];
		for(l=0;l<SIZE_EFB;l++) temp_cell.efbex[0][l] = cell->efbex[i][l];
		temp_cell.mpara[0][6] = temp_cell.collim.rmax;
		x_0 = temp_cell.mpara[0][1]*cos(temp_cell.mpara[0][0]);
		y_0 = temp_cell.mpara[0][1]*sin(temp_cell.mpara[0][0]);
		//printf("x0=%le,y0=%le,r0=%le,theta0=%le\n",x_0,y_0, temp_cell.mpara[0][1],temp_cell.mpara[0][0]);
		angle = temp_cell.mpara[0][0] + temp_cell.mpara[0][5];
		//find xstart
		xmax = sqrt(cell->collim.rmax*cell->collim.rmax - temp_cell.efbex[0][0]*temp_cell.efbex[0][0]) - temp_cell.mpara[0][1];
		xmin = temp_cell.mpara[0][1] - sqrt(cell->collim.rmin*cell->collim.rmin - temp_cell.efbex[0][0]*temp_cell.efbex[0][0]);
		//printf("xmin = %le, xmax = %le\n", xmin, xmax);
		if(xmin>xmax) {
			xmax = xmin;
			temp_cell.collim.rmax = 2*temp_cell.mpara[0][1] - temp_cell.collim.rmin;
		}
		else temp_cell.collim.rmax = cell->collim.rmax;
		//printf("collim: %le\n", temp_cell.collim.rmax);
		for(l=0;l<50;l++) {
			y = y_0 + temp_cell.efbex[0][0]*cos(angle) + xmax*sin(angle);
			x = x_0 - temp_cell.efbex[0][0]*sin(angle) + xmax*cos(angle);
			r = sqrt(x*x+y*y);
			//printf("l=%i,dif=%le, r=%le, xmax=%le\n",l,temp_cell.collim.rmax-r,r,xmax);
			if(fabs(temp_cell.collim.rmax-r)<prec) break;
			else xmax += 0.5*(temp_cell.collim.rmax - r);
		}
		temp_cell.collim.rmax += 0.1;
		xmin = 0.0;
		xtest = xmax;
		printf("xmax = %le\n",xmax);
		for(l=0;l<50;l++) {
			x = x_0 + xtest*cos(angle) - temp_cell.efbex[0][0]*sin(angle);
			y = y_0 + xtest*sin(angle) + temp_cell.efbex[0][0]*cos(angle);
			//printf("comp %i, try %i, x=%le, y=%le, r=%le, xtest = %le\n", i,l, x,y,sqrt(x*x+y*y),xtest);
			temp_cell.efben[0][2] = 2;
			//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, temp_cell.mpara[0][3], &bx2, &by2, &bz2, &temp_cell, add_contribution_comp);
			temp_cell.efben[0][2] = 10;
			//max2=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, temp_cell.mpara[0][3], &bx10, &by10, &bz10, &temp_cell, add_contribution_comp);
			if(fabs(bx10-bx2)>fabs(bx2) || fabs(by10-by2)>fabs(by2) || fabs(bz10-bz2)>fabs(bz2)) xmax=xtest;
			//if(max2 > max1) xmax=xtest;
			else xmin=xtest;
			//printf("\txmin=%le,xmax=%le,max1=%le,max2=%le\n",xmin,xmax,max1,max2);
			if((xmax-xmin)<10.*dx) break;
			xtest = (xmax+xmin)/2.;
		}
		cell->mpara[i][6] = xmin;
		printf("component %i, conv_lim = %le\n",i,xmin);
	}
	free(temp_cell.mpara[0]);
	free(temp_cell.efben[0]);
	free(temp_cell.efbex[0]);
	free(temp_cell.mpara);
	free(temp_cell.efben);
	free(temp_cell.efbex);
}

extern void read_fieldmap_boundaries(char *file, char *fileout)
{
	int j,flag=0;
	char buf[MAX_CHARINLINE];
	double x,y,z,bx,by,bz;
	FILE *rfile=NULL;
	FILE *wfile;
	
	wfile = fopen(fileout, "w");
	rfile = fopen(file,"r");
	if(rfile==NULL) errorstop("cannot open file");
	
	for(j = 0; j < 8; j++) newline(rfile);
	
	while(!feof(rfile)) {
	    if(fgets(buf, MAX_CHARINLINE-1, rfile) == NULL) break;
	    buf[MAX_CHARINLINE-1]='\0';

	    if(sscanf(buf, "%lf  %lf  %lf  %lf  %lf  %lf", &x,&z,&y,&bx,&bz,&by) != 6) {
	    	printf("why?\n");
			continue;
	    }
		if(fabs(z-0.25)<TINYDIMLESS) {
		//if(z==0 & flag==0) {
			printf("%lf  %lf  %lf  %lf  %lf  %lf\n", x,z,y,bx,bz,by);
			fprintf(wfile,"%le	%le	%le\n",x,y,z);
			//flag=1;
		}
		if(fabs(z)>TINYDIMLESS) flag=0;
	}
	fclose(rfile);
	fclose(wfile);
	
}

extern void remove_useless_lines_eritmap(char * file, char *fileout)
{
	int n=0;
	char buf[MAX_CHARINLINE];
	double x,y,z,bx,by,bz;
	double z_bef,y_bef;
	FILE *rfile=NULL;
	FILE *wfile;
	
	wfile = fopen(fileout, "w");
	rfile = fopen(file,"r");
	if(rfile==NULL) errorstop("cannot open file");
	
	while(!feof(rfile)) {
		n++;
	    if(fgets(buf, MAX_CHARINLINE-1, rfile) == NULL) break;
	    buf[MAX_CHARINLINE-1]='\0';
		
	    if(sscanf(buf, "%lf  %lf  %lf  %lf  %lf  %lf", &x,&z,&y,&bx,&bz,&by) != 6) {
	    	fprintf(wfile, "%s",buf);
			continue;
		}
		if(fabs(z)>TINYDIMLESS) {
		//if(y_bef<0 && y<0 && fabs(z-z_bef)>0.3 && fabs(z)<TINYDIMLESS) {
			printf("%i  %lf  %lf  %lf  %lf  %lf  %lf\n", n,x,z,y,bx,bz,by);
			continue;
		}
		else fprintf(wfile, "%s",buf);
		z_bef = z;
		y_bef = y;
	}
	fclose(rfile);
	fclose(wfile);
}

extern void compare_field(char *file_mod, char *file_map, struct Cell *cell_mod, struct Cell *cell_map, void(*add_contribution_comp_mod)(double,double,double,double*,double*,double*,struct Cell*,int), 
	void(*add_contribution_comp_map)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i,nb_th;
	double x,y,z,bx,by,bz;
	double r,th,th_step;
	
	FILE *wfile_mod;
	FILE *wfile_map;
	
	wfile_mod = fopen(file_mod,"w");
	wfile_map = fopen(file_map,"w");
	
	z = 0.05;
	r = 2.3;
	nb_th = 1800;
	th_step = cell_mod->boun.thmax/(nb_th);
	
	for(i=0;i<nb_th+1;i++) {
		th = i*th_step;
		x = r*cos(th);
		y = r*sin(th);
		printf("r=%lf, th=%le\n",r,th);
		if(get_bfield(x, y, z, &bx, &by, &bz, cell_mod, add_contribution_comp_mod)==ALIVE) {
			fprintf(wfile_mod,"%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", x, y, z, r, th, bx, by, bz);
		}
		else errorstop("bfield_mod not alive");
		if(get_bfield(x, y, z, &bx, &by, &bz, cell_map, add_contribution_comp_map)==ALIVE) {
			fprintf(wfile_map,"%lf	%lf	%lf	%lf	%lf	%lf	%lf	%lf\n", x, y, z, r, th, bx, by, bz);
		}
		else errorstop("bfield_map not alive");
	}
	fclose(wfile_mod);
	fclose(wfile_map);
}

extern void gnuplot_fit_fringe_field(char *txtfilename1, char *xcolumn1, char *ycolumn1, struct Cell *cell, double r, double a, char *with1, 
				 char *with2, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	char fitname[40];
	FILE *gp;
	
	sprintf(fitname ,"data/%s-fit.log",psfilename);
	remove(fitname);
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	gp=popen("/usr/local/bin/gnuplot -persist","w");
	fprintf(gp," set tics font 'times,30.'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'times,30.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -5,0 font 'times,30.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		//fprintf(gp," set terminal postscript color\n");
		fprintf(gp," set terminal postscript eps enhanced\n");
		fprintf(gp," set output 'output/%s.eps' \n", psfilename);
		fprintf(gp,"set fit logfile 'data/%s-fit.log'\n", psfilename);
		fprintf(gp," f(x) = ");
		
		fprintf(gp,"b*(%lf/%lf)**%lf/((1+exp(a*%lf/%lf*(%lf + %lf - x)))*(1+exp(a*%lf/%lf*(x - %lf - %lf)))) + ", 
		r, cell->mpara[0][1], cell->mpara[0][3], 
		cell->mpara[0][1], cell->efben[0][0], cell->mpara[0][0], cell->efben[0][0], 
		cell->mpara[0][1], cell->efbex[0][0], cell->mpara[0][0], cell->efbex[0][0]);
		fprintf(gp,"c*(%lf/%lf)**%lf/((1+exp(a*%lf/%lf*(%lf + %lf - x)))*(1+exp(a*%lf/%lf*(x - %lf - %lf)))) + ", 
		r, cell->mpara[1][1], cell->mpara[1][3], 
		cell->mpara[1][1], cell->efben[1][0], cell->mpara[1][0], cell->efben[1][0], 
		cell->mpara[1][1], cell->efbex[1][0], cell->mpara[1][0], cell->efbex[1][0]);
		fprintf(gp,"b*(%lf/%lf)**%lf/((1+exp(a*%lf/%lf*(%lf + %lf - x)))*(1+exp(a*%lf/%lf*(x - %lf - %lf))))", 
		r, cell->mpara[2][1], cell->mpara[2][3], 
		cell->mpara[2][1], cell->efben[2][0], cell->mpara[2][0], cell->efben[2][0], 
		cell->mpara[2][1], cell->efbex[2][0], cell->mpara[2][0], cell->efbex[2][0]);
		
		
		
		//for(ic=0;ic<cell->nbcomp;ic++) {
		//	fprintf(gp,"%le*(%le/%le)**%le/((1+exp(a*(%le + %le - x)))*(1+exp(a*(x - %le - %le)))", 
		//	cell->mpara[ic][2], r, cell->mpara[ic][1], cell->mpara[ic][3], cell->mpara[ic][0], cell->efben[ic][0], cell->mpara[ic][0], cell->efbex[ic][0]);
		//	if(ic<cell->nbcomp-1) fprintf(gp," + ");
		//}
		fprintf(gp,"\n");
		//fprintf(gp," f(x) = a*exp(-(x-b)*(x-b)/c/c)\n");
		fprintf(gp,"a = %le\n", a);
		fprintf(gp,"b = -1\n");
		fprintf(gp,"c = 1\n");
		
		//fprintf(gp,"b = %le\n", b);
		//fprintf(gp,"c = %le\n", c);
		
		fprintf(gp,"fit f(x) '%s' using %s:%s via a,b,c\n", txtfilename1, xcolumn1, ycolumn1);
		//fprintf(gp,"fit f(x) '%s' using %s:%s via a\n", txtfilename1, xcolumn1, ycolumn1);
		fprintf(gp,"plot '%s' using %s:%s with %s, f(x) with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, with2);
		//fprintf(gp," set output\n");
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using %s:%s with %s, f(x) with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, with2);
	pclose(gp);
}

// change field map for zgoubi (all coordinates are growing in the field map)
extern void shift_theta_field_map(char *readname, char *writename, int nsteps_r, int nsteps_th, int nsteps_z)
{
	char buf[MAX_CHARINLINE];
	int j,ir,ith,iz;
	double a, b, c, d, e, f;
	struct Map map;
	FILE *mapfile=NULL;
	FILE *wfile;
	
	map = allocmap(nsteps_r, nsteps_th, nsteps_z);
	
	mapfile = fopen(readname, "r");
	if(mapfile==NULL) errorstop("cannot open field map");
	
	wfile = fopen(writename, "w");
	
	for(j = 0; j < 8; j++) {
		//newline(mapfile);
		if(fgets(buf, MAX_CHARINLINE-1, mapfile) == NULL) errorstop("strange in shift_theta_field_map");
		fprintf(wfile,"%s", buf);
	}
	
	for(ir = 0; ir < nsteps_r; ir++) {
		for(iz = 0; iz < nsteps_z; iz++) {
			for(ith = 0; ith < nsteps_th; ith++) {
				fscanf(mapfile,"%le %le %le %le %le %le", &a, &b, &c, &d, &e, &f);
				printf("%le %le %le %le %le %le\n", a, b, c, d, e, f);
				map.node[ir][nsteps_th-1-ith][iz].coord[0] = a;
				map.node[ir][nsteps_th-1-ith][iz].coord[1] = b;
				map.node[ir][nsteps_th-1-ith][iz].coord[2] = c;
				map.node[ir][nsteps_th-1-ith][iz].b[0]	= d;
				map.node[ir][nsteps_th-1-ith][iz].b[1]	= e;
				map.node[ir][nsteps_th-1-ith][iz].b[2]	= f;
			}
		}
	}
	
	for(ir = 0; ir < nsteps_r; ir++) {
		for(iz = 0; iz < nsteps_z; iz++) {
			for(ith = 0; ith < nsteps_th; ith++) {
				fprintf(wfile,"%le	%le	%le	%le	%le	%le\n", 
				map.node[ir][ith][iz].coord[0], map.node[ir][ith][iz].coord[1], map.node[ir][ith][iz].coord[2],
				map.node[ir][ith][iz].b[0], map.node[ir][ith][iz].b[1], map.node[ir][ith][iz].b[2]);
			}
		}
	}
	fclose(mapfile);
	fclose(wfile);
	free_map(&map, nsteps_r, nsteps_th);
}

extern void test_symmetry_map(struct Map *map, int iz)
{
	char name_filex[200],name_filey[200],name_filez[200], outx[200],outy[200],outz[200];
	int ir,ith,nr,nth,nz;
	double x,y,r,th,z,bx,by,bz,bxs,bys,bzs,difx,dify,difz;
	double bminx=100.,bminy=100.,bminz=100.;
	double bmaxx=-100.,bmaxy=-100.,bmaxz=-100.;
	FILE *wfilex;
	FILE *wfiley;
	FILE *wfilez;
	sprintf(name_filex, "data/dif_map_z%i_x.dat",iz);
	sprintf(name_filey, "data/dif_map_z%i_y.dat",iz);
	sprintf(name_filez, "data/dif_map_z%i_z.dat",iz);
	sprintf(outx, "output/dif_map_z%i_x.eps",iz);
	sprintf(outy, "output/dif_map_z%i_y.eps",iz);
	sprintf(outz, "output/dif_map_z%i_z.eps",iz);
	wfilex = fopen(name_filex, "w");
	wfiley = fopen(name_filey, "w");
	wfilez = fopen(name_filez, "w");
	
	nr = map->nnodes[0];
	nth = map->nnodes[1];
	nz = map->nnodes[2];
	if(iz<0 || iz>nz-1) errorstop("vertical coordinate outside map!!");
	printf("symmetry test of map at z = %le\n",map->node[0][0][iz].coord[2]);
	for(ir=0;ir<nr;ir++) {
		for(ith=0;ith<nth/2;ith++) {
			r = map->node[ir][ith][iz].coord[0];
			th = map->node[ir][ith][iz].coord[1];
			z = map->node[ir][ith][iz].coord[2];
			x = r*cos(th);
			y = r*sin(th);
			bx = map->node[ir][ith][iz].b[0];
			by = map->node[ir][ith][iz].b[1];
			bz = map->node[ir][ith][iz].b[2];
			bxs = map->node[ir][nth-1-ith][iz].b[0];
			bys = map->node[ir][nth-1-ith][iz].b[1];
			bzs = map->node[ir][nth-1-ith][iz].b[2];
			difx = (bxs - bx)/bx;
			dify = (bys + by)/by;
			difz = (bzs - bz)/bz;
			//difx = (bxs - bx);
			//dify = (bys + by);
			//difz = (bzs - bz);
			if(fabs(bx)<0.001) difx = 0;
			if(fabs(by)<0.001) dify = 0;
			if(fabs(bz)<0.001) difz = 0;
			
			if(difx<bminx) bminx = difx;
			if(difx>bmaxx) bmaxx = difx;
			if(dify<bminy) bminy = dify;
			if(dify>bmaxy) bmaxy = dify;
			if(difz<bminz) bminz = difz;
			if(difz>bmaxz) bmaxz = difz;
			fprintf(wfilex, "%le	%le	%le\n",x,y,difx);
			fprintf(wfiley, "%le	%le	%le\n",x,y,dify);
			fprintf(wfilez, "%le	%le	%le\n",x,y,difz);
		}
		fprintf(wfilex, "\n");
		fprintf(wfiley, "\n");
		fprintf(wfilez, "\n");
	}
	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
	printf("min/max(x,y,z)=(%le/%le, %le/%le, %le/%le)\n", bminx,bmaxx,bminy,bmaxy,bminz,bmaxz);
	easyplot3dmap(name_filex, "x [m]", "y [m]", "b [T]", NULL, NULL, outx, "size ratio -1\nset xtics 0.2", bminx,bmaxx);
	easyplot3dmap(name_filey, "x [m]", "y [m]", "b [T]", NULL, NULL, outy, "size ratio -1\nset xtics 0.2", bminy,bmaxy);
	easyplot3dmap(name_filez, "x [m]", "y [m]", "b [T]", NULL, NULL, outz, "size ratio -1\nset xtics 0.2", bminz,bmaxz);
}

extern void check_diff_map_coordinates(struct Map *map_tosca, struct Map *map_mod)
{
	int ir,ith,iz,nr,nth;
	double diffr,diffth,r_tosca,r_mod,th_tosca,th_mod,x,y;
	double diffmin_r, diffmax_r, diffmin_th, diffmax_th;
	FILE *wfiler;
	FILE *wfileth;
	
	wfiler = fopen("data/dif_map_coord_r.dat", "w");
	wfileth = fopen("data/dif_map_coord_th.dat", "w");
	
	nr = map_tosca->nnodes[0];
	nth = map_tosca->nnodes[1];
	iz = 20;
	printf("coordinates comparison test of maps at z = %le\n",map_tosca->node[0][0][iz].coord[2]);
	for(ir=0;ir<nr;ir++) {
		for(ith=0;ith<nth;ith++) {
			r_tosca = map_tosca->node[ir][ith][iz].coord[0];
			th_tosca = map_tosca->node[ir][ith][iz].coord[1];
			r_mod = map_mod->node[ir][ith][iz].coord[0];
			th_mod = map_mod->node[ir][ith][iz].coord[1];
			x = r_tosca*cos(th_tosca);
			y = r_tosca*sin(th_tosca);
			diffr = r_tosca - r_mod;
			diffth = th_tosca - th_mod;
			if(diffr<diffmin_r) diffmin_r = diffr;
			if(diffr>diffmax_r) diffmax_r = diffr;
			if(diffth<diffmin_th) diffmin_th = diffth;
			if(diffth>diffmax_th) diffmax_th = diffth;
			fprintf(wfiler, "%le	%le	%le\n", x,y,diffr);
			fprintf(wfileth, "%le	%le	%le\n",x,y,diffth);
		}
		fprintf(wfiler, "\n");
		fprintf(wfileth, "\n");
	}
	fclose(wfiler);
	fclose(wfileth);
	printf("min/max(r,th)=(%le/%le, %le/%le)\n", diffmin_r,diffmax_r,diffmin_th,diffmax_th);
	easyplot3dmap("data/dif_map_coord_r.dat", "x [m]", "y [m]", NULL, NULL, NULL, "output/dif_map_coord_r.eps", "size ratio -1\nset xtics 0.2", diffmin_r,diffmax_r);
	easyplot3dmap("data/dif_map_coord_th.dat", "x [m]", "y [m]", NULL, NULL, NULL, "output/dif_map_coord_th.eps", "size ratio -1\nset xtics 0.2", diffmin_th,diffmax_th);
	
}

extern void test_maxwell_cylmap(struct Map *map, int ir, int iz, char *textfile)
{
	int ith;
	double r,th,z,stepr,stepth,stepz,br,bth,bz,brpdr,bthpdr,bzpdr,brpdth,bthpdth,bzpdth,brpdz,bthpdz,bzpdz,divb,curlbr,curlbth,curlbz;
	FILE *wfile;
	
	wfile = fopen(textfile, "w");
	stepr = map->stepsize[0];
	stepth = map->stepsize[1];
	stepz = map->stepsize[2];
	r = map->mapdim[0] + ir*stepr;
	z = map->mapdim[4] + iz*stepz;
	for(ith=0;ith<map->nnodes[1]-1;ith++) {
		th = ith*stepth;
		printf("r,th,z=(%le,%le,%le)\n",r,th,z);
		br  = map->node[ir][ith][iz].b[0]*cos(th) + map->node[ir][ith][iz].b[1]*sin(th);
		bth = map->node[ir][ith][iz].b[1]*cos(th) - map->node[ir][ith][iz].b[0]*sin(th);
		bz  = map->node[ir][ith][iz].b[2];
		brpdr  = map->node[ir+1][ith][iz].b[0]*cos(th) + map->node[ir+1][ith][iz].b[1]*sin(th);
		bthpdr = map->node[ir+1][ith][iz].b[1]*cos(th) - map->node[ir+1][ith][iz].b[0]*sin(th);
		bzpdr  = map->node[ir+1][ith][iz].b[2];
		brpdth  = map->node[ir][ith+1][iz].b[0]*cos(th) + map->node[ir][ith+1][iz].b[1]*sin(th);
		bthpdth = map->node[ir][ith+1][iz].b[1]*cos(th) - map->node[ir][ith+1][iz].b[0]*sin(th);
		bzpdth  = map->node[ir][ith+1][iz].b[2];
		brpdz  = map->node[ir][ith][iz+1].b[0]*cos(th) + map->node[ir][ith][iz+1].b[1]*sin(th);
		bthpdz = map->node[ir][ith][iz+1].b[1]*cos(th) - map->node[ir][ith][iz+1].b[0]*sin(th);
		bzpdz  = map->node[ir][ith][iz+1].b[2];
		//drbrdr = r*(brpdr - br)/stepr + br;dbthdth = (bthpdth - bth)/stepth;dbzdz = (bzpdz - bz)/stepz;divb = drbrdr/r + dbthdth/r + dbzdz;
		divb = (brpdr - br)/stepr + br/r + (bthpdth - bth)/stepth + (bzpdz - bz)/stepz;
		//dbzdth = (bzpdth - bz)/stepth;dbthdz = (bthpdz - bth)/stepz;curlbr = dbzdth/r - dbthdz;
		curlbr = (bzpdth - bz)/stepth/r - (bthpdz - bth)/stepz;
		//dbrdz = (brpdz - br)/stepz;dbzdr = (bzpdr - bz)/stepr;curlbth = dbrdz - dbzdr;
		curlbth = (brpdz - br)/stepz - (bzpdr - bz)/stepr;
		//drbthdr = r*(bthpdth - bth)/stepth + bth;dbrdth = (brpdth - br)/stepth;curlbz = (drbthdr - dbrdth)/r;
		curlbz = (bthpdth - bth)/stepth + (bth - (brpdth - br)/stepth)/r;
		fprintf(wfile, "%le	%le	%le	%le	%le	%le	%le\n",r,th*180./PI,z,divb,curlbr,curlbth, curlbz);
	}
	fclose(wfile);
}

extern void compute_map(char *filename, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double unitlength, double bscale,
	int loop_order, int line_order, int header_skip, double rstart, double rend, int nbr, double rshift, double thstart, double thend, int nbth, double thshift, double zstart, double zend, int nbz, double zshift)
{
	int i,i1,i2,i3,n1,n2,n3,ir,ith,iz;
	double x,y,z,bx,by,bz,r,th, rstep,thstep,zstep,coord[4],btemp[4];
	
	FILE *wfile;
	wfile = fopen(filename,"w");
	fprintf(wfile,"%i     %i     %i\n",nbr,nbth,nbz);
	for(i=0;i<header_skip-1;i++) fprintf(wfile,"\n"); //jump header_skip -1 lines of header

	//if(nbz>1) zstep = (zend-zstart)/(nbz-1);
	//else zstep = 0.;
	zstep = comp_step(zstart, zend, nbz-1);
	//if(nbr>1) rstep = (rend-rstart)/(nbr-1);
	//else rstep = 0.;
	rstep = comp_step(rstart, rend, nbr-1);
	//if(nbth>1) thstep = (thend-thstart)/(nbth-1);
	//else thstep = 0.;
	thstep = comp_step(thstart, thend, nbth-1);
	
	if(loop_order/100%10 == 1) n1=nbr;
	else if(loop_order/100%10 == 2) n1=nbth;
	else n1 = nbz;
	if(loop_order/10%10 == 1) n2=nbr;
	else if(loop_order/10%10 == 2) n2=nbth;
	else n2 = nbz;
	if(loop_order%10 == 1) n3=nbr;
	else if(loop_order%10 == 2) n3=nbth;
	else n3 = nbz;
	
	for(i1 = 0; i1 < n1; i1++) {
		if(loop_order/100%10 == 1) ir=i1;
		else if(loop_order/100%10 == 2) ith=i1;
		else iz = i1;
		for(i2 = 0; i2 < n2; i2++) {
			if(loop_order/10%10 == 1) ir=i2;
			else if(loop_order/10%10 == 2) ith=i2;
			else iz = i2;
			for(i3 = 0; i3 < n3; i3++) {
				if(loop_order%10 == 1) ir=i3;
				else if(loop_order%10 == 2) ith=i3;
				else iz = i3;
				r = rstart + ir*rstep;
				th = thstart + ith*thstep;
				z = zstart + iz*zstep;
				x = r*cos(th);
				y = r*sin(th);
				//printf("r=%lf, th=%le, z=%le\n",r,th,z);
				if(get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp)==ALIVE) {
					x = (r+rshift)*cos(th+thshift);
					y = (r+rshift)*sin(th+thshift);
					z += zshift;
					if(line_order/100%10 == 1) {
						coord[1] = x;
						btemp[1] = bx;
					}
					else if(line_order/100%10 == 2) {
						coord[1] = y;
						btemp[1] = by;
					}
					else {
						coord[1] = z;
						btemp[1] = bz;
					}
					if(line_order/10%10 == 1) {
						coord[2] = x;
						btemp[2] = bx;
					}
					else if(line_order/10%10 == 2) {
						coord[2] = y;
						btemp[2] = by;
					}
					else {
						coord[2] = z;
						btemp[2] = bz;
					}
					if(line_order%10 == 1) {
						coord[3] = x;
						btemp[3] = bx;
					}
					else if(line_order%10 == 2) {
						coord[3] = y;
						btemp[3] = by;
					}
					else {
						coord[3] = z;
						btemp[3] = bz;
					}
					fprintf(wfile,"%lf	%lf	%lf	%lf	%lf	%lf\n", coord[1]*unitlength, coord[2]*unitlength, coord[3]*unitlength, btemp[1]/bscale, btemp[2]/bscale, btemp[3]/bscale);
				}
				else errorstop("bfield not alive");
			}
		}
	}
	fclose(wfile);
}

// works only with a single cell lattice with a 3D map
extern void create_sym_map_latt(struct Lattice *latt_ori, struct Lattice *latt_sym)
{
	int i, ir,ith,iz,ith2,nsteps_r,nsteps_th,nsteps_z;
	double th2;
	
	latt_sym->periodicity = latt_ori->periodicity;
	//latt_sym->nbcell = 		latt_ori->nbcell;
	latt_sym->nbcell = 1;
	latt_sym->cell = alloccell(latt_sym->nbcell);
	
	latt_sym->cell[0].nbcomp = 		latt_ori->cell[0].nbcomp;
	latt_sym->cell[0].stepsize = 	latt_ori->cell[0].stepsize;
	strcpy(latt_sym->cell[0].keyword,latt_ori->cell[0].keyword);
	latt_sym->cell[0].instrutype = 	latt_ori->cell[0].instrutype;
	latt_sym->cell[0].framework = 	latt_ori->cell[0].framework;
	latt_sym->cell[0].deltar = 		latt_ori->cell[0].deltar;
	latt_sym->cell[0].collim.rmin = latt_ori->cell[0].collim.rmin;
	latt_sym->cell[0].collim.rmax = latt_ori->cell[0].collim.rmax;
	latt_sym->cell[0].collim.zmin = latt_ori->cell[0].collim.zmin;
	latt_sym->cell[0].collim.zmax = latt_ori->cell[0].collim.zmax;
	latt_sym->cell[0].boun.thmax = 	latt_ori->cell[0].boun.thmax;
	latt_sym->cell[0].boun.ymax = 	latt_ori->cell[0].boun.ymax;
	
	nsteps_r = latt_ori->cell[0].map.nnodes[0];
	nsteps_th= latt_ori->cell[0].map.nnodes[1];
	nsteps_z = latt_ori->cell[0].map.nnodes[2];
	
	latt_sym->cell[0].map = allocmap(nsteps_r, nsteps_th, nsteps_z);
	latt_sym->cell[0].map.sym = latt_ori->cell[0].map.sym;
	
	for(i=0;i<3;i++) latt_sym->cell[0].map.nnodes[i] = 	latt_ori->cell[0].map.nnodes[i];
	for(i=0;i<3;i++) latt_sym->cell[0].map.stepsize[i]= latt_ori->cell[0].map.stepsize[i];
	for(i=0;i<6;i++) latt_sym->cell[0].map.mapdim[i] = 	latt_ori->cell[0].map.mapdim[i];
	
	for(ir=0;ir<nsteps_r;ir++) {
		for(iz=0;iz<nsteps_z;iz++) {
			for(ith=0;ith<(nsteps_th-1)/2+1;ith++) {
				ith2 = nsteps_th-1-ith;
				th2 = latt_sym->cell[0].map.stepsize[1]*ith2;
				for(i=0;i<3;i++) {
					latt_sym->cell[0].map.node[ir][ith][iz].coord[i]= 	latt_ori->cell[0].map.node[ir][ith][iz].coord[i]; 
					latt_sym->cell[0].map.node[ir][ith][iz].b[i]	= 	latt_ori->cell[0].map.node[ir][ith][iz].b[i];
				}
				latt_sym->cell[0].map.node[ir][ith2][iz].coord[0] = 	latt_ori->cell[0].map.node[ir][ith][iz].coord[0]; 
				//latt_sym->cell[0].map.node[ir][ith2][iz].coord[1] = th2; 
				latt_sym->cell[0].map.node[ir][ith2][iz].coord[1] = 	latt_ori->cell[0].map.node[ir][ith2][iz].coord[1];
				latt_sym->cell[0].map.node[ir][ith2][iz].coord[2] = 	latt_ori->cell[0].map.node[ir][ith][iz].coord[0]; 
				latt_sym->cell[0].map.node[ir][ith2][iz].b[0] = 		latt_ori->cell[0].map.node[ir][ith][iz].b[0];
				latt_sym->cell[0].map.node[ir][ith2][iz].b[1] = 		-latt_ori->cell[0].map.node[ir][ith][iz].b[1];
				latt_sym->cell[0].map.node[ir][ith2][iz].b[2] = 		latt_ori->cell[0].map.node[ir][ith][iz].b[2];
			}
		}
	}
}

extern void compare_map_mapfile(struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), char *mapfile, double bscale, double coord_scale, char *outfile)
{
	int ir,it,iz,flag_bfield;
	double x,y,z,bx,by,bz,bx_file,by_file,bz_file,difbx,difby,difbz,r,th;
	FILE *rfile=NULL;
	FILE *wfile;
	rfile = fopen(mapfile, "r");
	if(rfile==NULL) errorstop("cannot open map file\n");
	for(ir=0;ir<8;ir++) newline(rfile);
	
	wfile = fopen(outfile,"w");
	
	for(ir=0;ir<cell->map.nnodes[0];ir++) {
		for(iz=0;iz<cell->map.nnodes[2];iz++) {
			for(it=0;it<cell->map.nnodes[1];it++) {
				fscanf(rfile, "%le	%le	%le	%le	%le	%le", &x,&z,&y,&bx_file,&bz_file,&by_file);
				x/=coord_scale;
				y/=coord_scale;
				r=sqrt(x*x+y*y);
				th = atan_ratio(y, x);
				x=r*cos(th+15.*PI/180.); //only for main ring map!!!
				y=r*sin(th+15.*PI/180.); //only for main ring map!!!
				z/=coord_scale;
				if(fabs(y)<TINYLENGTH*10.) ;
				else {
					//printf("avant:ir=%i, it=%i, iz=%i\n",ir,it,iz);
					flag_bfield = get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp);
					if(flag_bfield!=ALIVE) errorstop("bfield not alive");
					difbx = bx - bx_file/bscale;
					difby = by - by_file/bscale;
					difbz = bz - bz_file/bscale;
					//if(fabs(difbx)>TINYDIMLESS) printf("not good for difbx=%le at x=%le,y=%le,z=%le\n",bx,x,y,z);
					//if(fabs(difby)>TINYDIMLESS) printf("not good for difby=%le at x=%le,y=%le,z=%le\n",by,x,y,z);
					//if(fabs(difbz)>TINYDIMLESS) printf("not good for difbz=%le at x=%le,y=%le,z=%le\n",bz,x,y,z);
					if(fabs(difbx)>1.e-9 || fabs(difby)>1.e-9 || fabs(difbz)>1.e-9) fprintf(wfile,"%le	%le	%le	%le	%le	%le\n",x,y,z,difbx,difby,difbz);
				}
			}
		}
	}
	fclose(rfile);
	fclose(wfile);
}

extern void filter_fieldmap(char *fieldmap, char *output)
{
	int i,nblines,jumphead =8;
	double x,y,z,bx,by,bz;
	FILE *rfile=NULL;
	FILE *wfile;
	nblines = get_nb_lines_file(fieldmap);
	rfile=fopen(fieldmap, "r");
	if(rfile==NULL) errorstop("cannot open fieldmap");
	for(i=0;i<jumphead;i++) newline(rfile);
	wfile = fopen(output,"w");
	
	for(i=0;i<nblines-jumphead;i++) {
		fscanf(rfile,"%le	%le	%le	%le	%le	%le",&x,&z,&y,&bx,&bz,&by);
		printf("%le	%le	%le	%le	%le	%le\n",x,z,y,bx,bz,by);
		//if(bx<0 && by<0 && bz<0) {
			//if(bx>-2.e-2 && bx<-0.5e-2 && by>-5.5e-2 && by<-3.e-2 && bz>-6.0 && bz<-3.5) {
			if(bz>-6.0 && bz<-3.5) {
					fprintf(wfile,"%le	%le	%le	%le	%le	%le\n",x,z,y,bx,bz,by);
			}
			//}
	}
	
	fclose(rfile);
	fclose(wfile);
}

extern void da_study_fets(struct Lattice *latt, struct Particle *part_ref, char *fileout, double angle)
{
	double qx,qz,beta_x,beta_z,alpha_x,alpha_z, ex,ez, axmax, azmax, x_offset,z_offset,gamma_x,gamma_z;
	FILE *wfile;
	wfile = fopen(fileout,"a");
	
	tune_calc_matrix(part_ref, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-4, 1.e-4, 1.e-4, 1.e-4, latt, part_cross_latt, YES, NULL);
	
	gamma_x = (1+alpha_x*alpha_x)/beta_x;
	gamma_z = (1+alpha_z*alpha_z)/beta_z;
	x_offset = sqrt(10.e-6/gamma_x);
	z_offset = sqrt(10.e-6/gamma_z);
	printf("offset x,z: %le, %le\n",x_offset, z_offset);
	acceptancex_auto(NULL, part_ref, latt, 2500, 0.001, &axmax, 0.001, z_offset); 
	acceptancez_auto(NULL, part_ref, latt, 2500, 0.001, &azmax, 0.001, x_offset);
	
	ex = axmax*axmax*gamma_x;
	ez = azmax*azmax*gamma_z;
	fprintf(wfile,"%lf	%le	%le	%le	%le	%le	%le\n", angle, qx,qz,ex,ez,axmax,azmax);
	fclose(wfile);
}

extern void scan_k2(char *outfile, struct Lattice *latt, struct Particle *part_ref, double k2_start, double k2_end, int nbstep_k2, double eps_clo, double percentage_mom, double amp)
{
	int i,n;
	double k2, k2_step, qx,qz,qx2,qz2,beta_x,alpha_x,beta_z,alpha_z,dqx,dqz;
	struct Particle test_part, test_part2;
	FILE *wfile;
	wfile =fopen(outfile,"w");
	
	//k2_step = (k2_end - k2_start)/nbstep_k2;
	k2_step = comp_step(k2_start, k2_end, nbstep_k2);

	for(i=0;i<nbstep_k2;i++) {
		test_part = *part_ref;
		test_part.hat = -2;
		test_part2 = *part_ref;
		test_part2.hat = -2;
		test_part2.brho *=percentage_mom;
		k2 = k2_start+k2_step*i;
		printf("k2=%le\n",k2);
		for(n=0;n<latt->cell[0].nbcomp;n++) latt->cell[0].mpara[n][4] = k2;
		if(find_closed_orbite_xxp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), eps_clo, latt, YES) == FALSE) errorstop("closed orbit not found\n");
		if(find_closed_orbite_xxp(&test_part2, &(test_part2.x), &(test_part2.ux), &(test_part2.uy), eps_clo, latt, YES) == FALSE) errorstop("closed orbit not found\n");
		tune_calc_matrix(&test_part, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, amp, amp, amp, amp, latt, part_cross_latt, YES, NULL);
		tune_calc_matrix(&test_part2, &qx2, &qz2, &beta_x, &alpha_x, &beta_z, &alpha_z, amp, amp, amp, amp, latt, part_cross_latt, YES, NULL);
		dqx = qx2-qx;
		dqz = qz2-qz;
		fprintf(wfile, "%le	%le	%le	%le	%le	%le	%le\n",dqx,dqz, k2, qx,qz,qx2,qz2);
	}
	fclose(wfile);
}//*/

extern void compare_diffield_added_separate(char *file_bx, char *file_by, char *file_bz, double r0, double rmax, int r_step, double th0, double thmax, int th_step, double z, struct Lattice *latt1, struct Lattice *latt2, void(*add_contribution_comp1)(double,double,double,double*,double*,double*,struct Cell*,int), 
	void(*add_contribution_comp2)(double,double,double,double*,double*,double*,struct Cell*,int), double *bxmin, double *bxmax, double *bymin, double *bymax, double *bzmin, double *bzmax)
{
	int i,j;
	double x,y,bx1,by1,bz1, bx2,by2,bz2, difx,dify,difz;
	double r, th, step_r, step_th;
	double minx1,maxx1,miny1,maxy1,minz1,maxz1, minx2,maxx2,miny2,maxy2,minz2,maxz2;
	
	FILE *wfile_bx;
	FILE *wfile_by;
	FILE *wfile_bz;
	wfile_bx = fopen(file_bx,"w");
	wfile_by = fopen(file_by,"w");
	wfile_bz = fopen(file_bz,"w");
	
	FILE *wfile_x1;
	FILE *wfile_y1;
	FILE *wfile_z1;
	FILE *wfile_x2;
	FILE *wfile_y2;
	FILE *wfile_z2;
	wfile_x1 = fopen("data/fieldbx_1.dat","w");
	wfile_y1 = fopen("data/fieldby_1.dat","w");
	wfile_z1 = fopen("data/fieldbz_1.dat","w");
	wfile_x2 = fopen("data/fieldbx_2.dat","w");
	wfile_y2 = fopen("data/fieldby_2.dat","w");
	wfile_z2 = fopen("data/fieldbz_2.dat","w");
	
	
	*bxmin = 1.e9;
	*bxmax = -1.e9;
	*bymin = 1.e9;
	*bymax = -1.e9;
	*bzmin = 1.e9;
	*bzmax = -1.e9;
	
	minx1 = 1.e9;
	maxx1 = -1.e9;
	miny1 = 1.e9;
	maxy1 = -1.e9;
	minz1 = 1.e9;
	maxz1 = -1.e9;
	minx2 = 1.e9;
	maxx2 = -1.e9;
	miny2 = 1.e9;
	maxy2 = -1.e9;
	minz2 = 1.e9;
	maxz2 = -1.e9;
		
	
	//if(r_step == 1) step_r = 0.;
	//else step_r = (rmax-r0)/(r_step);
	step_r = comp_step(r0, rmax, r_step);
	//if(th_step == 1) step_th = 0.;
	//else step_th = (thmax-th0)/(th_step);
	step_th = comp_step(th0, thmax, th_step);
	//printf("step_r = %le, step_th = %le\n", step_r, step_th);
	
	for(j=0;j<th_step;j++) {
		th = th0 + j*step_th;
		for(i=0;i<r_step;i++) {
			r=r0 + i*step_r;
			x = r*cos(th);
			y = r*sin(th);
			//printf("r=%lf, th=%le, thmax=%le\n",r,th, latt2->cell[0].boun.thmax);
			//printf("x=%le,y=%le\n",x,y);
			if(get_bfield(x, y, z, &bx1, &by1, &bz1, &(latt1->cell[0]), add_contribution_comp1)!=ALIVE) errorstop("bfield of latt1 not alive");
			if(th<latt2->cell[0].boun.thmax) {
				//printf("2: x=%le,y=%le\n",x,y);
				if(get_bfield(x, y, z, &bx2, &by2, &bz2, &(latt2->cell[0]), add_contribution_comp2)!=ALIVE) errorstop("bfield of latt2, cell0 not alive");
			}
			else {
				fwk_pt_globtoloc(&x, &y, &(latt2->cell[1].framework));
				//printf("2: x=%le,y=%le\n",x,y);
				if(get_bfield(x, y, z, &bx2, &by2, &bz2, &(latt2->cell[1]), add_contribution_comp2)!=ALIVE) errorstop("bfield of latt2, cell1 not alive");
				fwk_pt_loctoglob(&x, &y, &(latt2->cell[1].framework));
			}
			difx = bx1-bx2;
			*bxmin = MIN(*bxmin,difx);
			*bxmax = MAX(*bxmax,difx);
			dify = by1-by2;
			*bymin = MIN(*bymin,dify);
			*bymax = MAX(*bymax,dify);
			difz = bz1-bz2;
			*bzmin = MIN(*bzmin,difz);
			*bzmax = MAX(*bzmax,difz);
			fprintf(wfile_bx,"%le	%le	%le\n", y, x, difx);
			fprintf(wfile_by,"%le	%le	%le\n", y, x, dify);
			fprintf(wfile_bz,"%le	%le	%le\n", y, x, difz);
			
			minx1 = MIN(minx1,bx1);
			maxx1 = MAX(maxx1,bx1);
			miny1 = MIN(miny1,by1);
			maxy1 = MAX(maxy1,by1);
			minz1 = MIN(minz1,bz1);
			maxz1 = MAX(maxz1,bz1);
			minx2 = MIN(minx2,bx2);
			maxx2 = MAX(maxx2,bx2);
			miny2 = MIN(miny2,by2);
			maxy2 = MAX(maxy2,by2);
			minz2 = MIN(minz2,bz2);
			maxz2 = MAX(maxz2,bz2);
			fprintf(wfile_x1,"%le	%le	%le\n", y, x, bx1);
			fprintf(wfile_y1,"%le	%le	%le\n", y, x, by1);
			fprintf(wfile_z1,"%le	%le	%le\n", y, x, bz1);
			fprintf(wfile_x2,"%le	%le	%le\n", y, x, bx2);
			fprintf(wfile_y2,"%le	%le	%le\n", y, x, by2);
			fprintf(wfile_z2,"%le	%le	%le\n", y, x, bz2);
		}
		fprintf(wfile_bx,"\n");
		fprintf(wfile_by,"\n");
		fprintf(wfile_bz,"\n");
		
		fprintf(wfile_x1,"\n");
		fprintf(wfile_y1,"\n");
		fprintf(wfile_z1,"\n");
		fprintf(wfile_x2,"\n");
		fprintf(wfile_y2,"\n");
		fprintf(wfile_z2,"\n");
	}
	fclose(wfile_bx);
	fclose(wfile_by);
	fclose(wfile_bz);
	
	fclose(wfile_x1);
	fclose(wfile_y1);
	fclose(wfile_z1);
	fclose(wfile_x2);
	fclose(wfile_y2);
	fclose(wfile_z2);
	easyplot3dmap("data/fieldbx_1.dat", "x [m]", "y [m]", "bx [T]", NULL, NULL, "output/fieldbx_1.eps", "size ratio -1\nset grid", minx1,maxx1);
	easyplot3dmap("data/fieldby_1.dat", "x [m]", "y [m]", "by [T]", NULL, NULL, "output/fieldby_1.eps", "size ratio -1\nset grid", miny1,maxy1);
	easyplot3dmap("data/fieldbz_1.dat", "x [m]", "y [m]", "bz [T]", NULL, NULL, "output/fieldbz_1.eps", "size ratio -1\nset grid", minz1,maxz1);
	easyplot3dmap("data/fieldbx_2.dat", "x [m]", "y [m]", "bx [T]", NULL, NULL, "output/fieldbx_2.eps", "size ratio -1\nset grid", minx2,maxx2);
	easyplot3dmap("data/fieldby_2.dat", "x [m]", "y [m]", "by [T]", NULL, NULL, "output/fieldby_2.eps", "size ratio -1\nset grid", miny2,maxy2);
	easyplot3dmap("data/fieldbz_2.dat", "x [m]", "y [m]", "bz [T]", NULL, NULL, "output/fieldbz_2.eps", "size ratio -1\nset grid", minz2,maxz2);
}

extern double fn_co(struct Lattice *latt, struct Particle *part, double x[], double betax, double betaz)
{
	double xi[4], xf[4];
	struct Particle part_temp;
	part_temp = *part;
	part_temp.x  = x[0];
	part_temp.ux = x[1];
	part_temp.z  = x[2];
	part_temp.uz = x[3];
	xi[0] = part_temp.x;
	xi[1] = part_temp.ux;
	xi[2] = part_temp.z;
	xi[3] = part_temp.uz;
	//printf("%le, %le, %le, %le\n",xi[0],xi[1],xi[2],xi[3]);
	part_cross_latt(&part_temp, latt, NULL);
	if(part_temp.status != ALIVE) errorstop("part killed while crossing latt");
	xf[0] = part_temp.x;
	xf[1] = part_temp.ux;
	xf[2] = part_temp.z;
	xf[3] = part_temp.uz;
	return (xf[0]-xi[0])*(xf[0]-xi[0]) + betax*(xf[1]-xi[1])*(xf[1]-xi[1]) + (xf[2]-xi[2])*(xf[2]-xi[2]) + betaz*(xf[3]-xi[3])*(xf[3]-xi[3]);
}

extern void find_co_nelmin(struct Lattice *latt, struct Particle *part, double betax, double betaz, double fn_co(struct Lattice*, struct Particle*, double x[], double, double), 
	int n, double xmin[], double *ynewlo, double reqmin, double step[], int konvge, int kcount, int *icount, int *numres, int *ifault)
{
	int i, ihi, ilo, j, jcount, l, nn;
	double ccoeff = 0.5, del, dn, dnn, ecoeff = 2.0, eps = 1.0, rcoeff = 1.0, rq, x, *y, y2star, ylo, ystar, z, *p, *p2star, *pbar, *pstar;
	
	//Check the input parameters.
	if(reqmin <= 0.0 || n < 1 || konvge < 1) {
		*ifault = 1;
		return;
	}
	
	double start[n];
	struct Particle part_temp;
	start[0] = part->x;
	start[1] = part->ux;
	start[2] = part->z;
	start[3] = part->uz;
	
	p = (double *) malloc(n * (n + 1) * sizeof(double));
	pstar = (double *) malloc(n * sizeof(double));
	p2star = (double *) malloc(n * sizeof(double));
	pbar = (double *) malloc(n * sizeof(double));
	y = (double *) malloc((n + 1) * sizeof(double));
	
	*icount = 0;
	*numres = 0;
	
	jcount = konvge; 
	dn = (double) (n);
	nn = n + 1;
	dnn = (double) (nn);
	del = 1.0;
	rq = reqmin * dn;
	
	//Initial or restarted loop.
	
	for(;;) {
		for(i = 0; i < n; i++) p[i+n*n] = start[i];
		//y[n] = fn(start);
		part_temp = *part;
		y[n] = fn_co(latt, &part_temp, start, betax, betaz);
		*icount = *icount + 1;
	
		for(j = 0; j < n; j++) {
			x = start[j];
			start[j] = start[j] + step[j] * del;
			for(i = 0; i < n; i++) p[i+j*n] = start[i];
			//y[n] = fn(start);
			part_temp = *part;
			y[j] = fn_co(latt, &part_temp, start, betax, betaz);
			*icount = *icount + 1;
			start[j] = x;
		}
		//The simplex construction is complete.
		//Find highest and lowest Y values.	YNEWLO = Y(IHI) indicates
		//the vertex of the simplex to be replaced.
		ylo = y[0];
		ilo = 0;
		for(i = 1; i < nn; i++) {
			if(y[i] < ylo) {
				ylo = y[i];
				ilo = i;
			}
		}
		//CLRSCR();
		//printf("try # %i, ylo = %le, x=%le, ux=%le, z=%le, uz=%le\n", *icount, ylo, part_temp.x, part_temp.ux, part_temp.z, part_temp.uz);
		//fflush(stdout);
		//Inner loop.
		for(;;) {
			if(kcount <= *icount) break;
			*ynewlo = y[0];
			ihi = 0;
			for(i = 1; i < nn; i++) {
				if(*ynewlo < y[i]) {
					*ynewlo = y[i];
					ihi = i;
				}
			}
			//CLRSCR();
			//printf("try # %i, ynewlo = %le, x=%le, ux=%le, z=%le, uz=%le\n", *icount, *ynewlo, part_temp.x, part_temp.ux, part_temp.z, part_temp.uz);
			//printf("try # %i, ynewlo = %le\n", *icount, *ynewlo);
			//fflush(stdout);
			//printf("ynewlo=%le\n", *ynewlo);
			//Calculate PBAR, the centroid of the simplex vertices
			//excepting the vertex with Y value YNEWLO.
			for(i = 0; i < n; i++) {
				z = 0.0;
				for(j = 0; j < nn; j++) z = z + p[i+j*n];
				z = z - p[i+ihi*n];	
				pbar[i] = z / dn;
			}
			//Reflection through the centroid.
			for(i = 0; i < n; i++) pstar[i] = pbar[i] + rcoeff * (pbar[i] - p[i+ihi*n]);
			//ystar = fn(pstar);
			part_temp = *part;
			ystar = fn_co(latt, &part_temp, pstar, betax, betaz);
			*icount = *icount + 1;
			//Successful reflection, so extension.
			if(ystar < ylo) {
				for(i = 0; i < n; i++) p2star[i] = pbar[i] + ecoeff * (pstar[i] - pbar[i]);
				//y2star = fn(p2star);
				part_temp = *part;
				y2star = fn_co(latt, &part_temp, p2star, betax, betaz);
				*icount = *icount + 1;
				//Check extension.
				if(ystar < y2star) {
					for(i = 0; i < n; i++) p[i+ihi*n] = pstar[i];
					y[ihi] = ystar;
				}
				//Retain extension or contraction.
				else {
					for(i = 0; i < n; i++) p[i+ihi*n] = p2star[i];
					y[ihi] = y2star;
				}
			}
			//No extension.
			else {
				l = 0;
				for(i = 0; i < nn; i++) if(ystar < y[i]) l = l + 1;
	
				if(l > 1) {
					for(i = 0; i < n; i++) p[i+ihi*n] = pstar[i];
					y[ihi] = ystar;
				}
				//Contraction on the Y(IHI) side of the centroid.
				else if(l == 0) {
					for(i = 0; i < n; i++) p2star[i] = pbar[i] + ccoeff * (p[i+ihi*n] - pbar[i]);
					//y2star = fn(p2star);
					part_temp = *part;
					y2star = fn_co(latt, &part_temp, p2star, betax, betaz);
					*icount = *icount + 1;
					//Contract the whole simplex.
					if(y[ihi] < y2star) {
						for(j = 0; j < nn; j++) {
							for(i = 0; i < n; i++) {
								p[i+j*n] = (p[i+j*n] + p[i+ilo*n]) * 0.5;
								xmin[i] = p[i+j*n];
							}
							//y[j] = fn(xmin);
							part_temp = *part;
							y[j] = fn_co(latt, &part_temp, xmin, betax, betaz);
							*icount = *icount + 1;
						}
						ylo = y[0];
						ilo = 0;
						for(i = 1; i < nn; i++) {
							if(y[i] < ylo) {
								ylo = y[i];
								ilo = i;
							}
						}
						continue;
					}
					//Retain contraction.
					else {
						for(i = 0; i < n; i++) p[i+ihi*n] = p2star[i];
						y[ihi] = y2star;
					}
				}
				//Contraction on the reflection side of the centroid.
				else if(l == 1) {
					for(i = 0; i < n; i++) p2star[i] = pbar[i] + ccoeff * (pstar[i] - pbar[i]);
					//y2star = fn(p2star);
					part_temp = *part;
					y2star = fn_co(latt, &part_temp, p2star, betax, betaz);
					*icount = *icount + 1;
					//Retain reflection?
					if(y2star <= ystar) {
						for(i = 0; i < n; i++) p[i+ihi*n] = p2star[i];
						y[ihi] = y2star;
					}
					else {
						for(i = 0; i < n; i++) p[i+ihi*n] = pstar[i];
						y[ihi] = ystar;
					}
				}
			}
			//Check if YLO improved.
			if(y[ihi] < ylo) {
				ylo = y[ihi];
				ilo = ihi;
			}
			jcount = jcount - 1;
	
			if(0 < jcount) continue;
			//Check to see if minimum reached.
			if(*icount <= kcount) {
				jcount = konvge;
				z = 0.0;
				for(i = 0; i < nn; i++) z = z + y[i];
				x = z / dnn;
	
				z = 0.0;
				for(i = 0; i < nn; i++) z = z + pow(y[i] - x, 2);
				if(z <= rq) break;
			}
		}
		//Factorial tests to check that YNEWLO is a local minimum.
		for(i = 0; i < n; i++) xmin[i] = p[i+ilo*n];
		*ynewlo = y[ilo];
		//CLRSCR();
		//printf("trys # %i, eps = %le\n", *icount, *ynewlo);
		//printf("trys # %i, ynewlo = %le, x=%le, ux=%le, z=%le, uz=%le", *icount, *ynewlo, part_temp.x, part_temp.ux, part_temp.z, part_temp.uz);
		//fflush(stdout);
	
		if(kcount < *icount) {
			*ifault = 2;
			break;
		}
		*ifault = 0;
		for(i = 0; i < n; i++) {
			del = step[i] * eps;
			xmin[i] = xmin[i] + del;
			//z = fn(xmin);
			part_temp = *part;
			z = fn_co(latt, &part_temp, xmin, betax, betaz);
			*icount = *icount + 1;
			if(z < *ynewlo) {
				*ifault = 2;
				break;
			}
			xmin[i] = xmin[i] - del - del;
			//z = fn(xmin);
			part_temp = *part;
			z = fn_co(latt, &part_temp, xmin, betax, betaz);
			*icount = *icount + 1;
			if(z < *ynewlo) {
				*ifault = 2;
				break;
			}
			xmin[i] = xmin[i] + del;
		}
		if(*ifault == 0) break;
		//Restart the procedure.
		for(i = 0; i < n; i++) start[i] = xmin[i];
		del = eps;
		*numres = *numres + 1;
	}
	free(p);
	free(pstar);
	free(p2star);
	free(pbar);
	free(y);
	return;
}

extern int put_on_co_nelmin(struct Lattice *latt, struct Particle *part, double reqmin, double stepx, double stepux, double stepz, double stepuz, int doyouprintf)
{
	int icount, numres, ifault;
	double xco[4], ynewlo, step[4];
	int n=4;
	int konvge=1;
	int kcount=1e6;
	double betax=1.0;//10.0;
	double betaz=1.0;//10.0;
	step[0] = stepx;
	step[1] = stepux;
	step[2] = stepz;
	step[3] = stepuz;
	
	find_co_nelmin(latt, part, betax, betaz, fn_co, n, xco, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault);
	if(ifault==0) {
		part->x  = xco[0];
		part->ux = xco[1];
		part->z  = xco[2];
		part->uz = xco[3];
	}
	if(doyouprintf==YES) {
		if(reqmin>ynewlo) printf("\nClosed orbit found: x_clo = %.8f [m], xprime_clo = %le [deg], \nz_clo = %lf [m], zprime_clo = %lf [deg], eps = %le < %le\n", part->x, atan(part->ux/part->uy)*180./(PI), part->z, atan(part->uz/part->uy)*180./(PI), ynewlo, reqmin);
		else printf("\nCLOSED ORBIT NOT FOUND: x = %.8f [m], xprime = %le [deg], \nz = %lf [m], zprime = %lf [deg], eps = %le > %le\n", part->x, atan(part->ux/part->uy)*180./(PI), part->z, atan(part->uz/part->uy)*180./(PI), ynewlo, reqmin);
		printf("nb of iterations=%i, nb of restarts=%i, ifault=%i\n", icount, numres, ifault);
		
	}
	if(ifault==0) return TRUE;
	else return FALSE;
}

extern void maxwell_test_vffa(char *file_input, struct Lattice *latt, double interpol_order, double r0, double rmax, int r_step, double th0, double thmax, int th_step, double z, double dx, double dy, double dz, int doyoucirc)
{
	char name_div[300], name_curlx[300], name_curly[300], name_curlz[300], name_eps[300];
	int i,j;
	double th0_temp, xmin, xmax, ymin, ymax, valmin, valmax;
	
	sprintf(name_div,  "%s_div.dat"  , file_input);
	sprintf(name_curlx,"%s_curlx.dat", file_input);
	sprintf(name_curly,"%s_curly.dat", file_input);
	sprintf(name_curlz,"%s_curlz.dat", file_input);
	emptyfile(name_div);
	emptyfile(name_curlx);
	emptyfile(name_curly);
	emptyfile(name_curlz);
	
	for(i=0;i<latt->nbcell;i++) for(j=0;j<latt->cell[i].nbcomp;j++) latt->cell[i].efben[j][2] = interpol_order;
	for(i=0;i<latt->nbcell;i++) {
		if(fabs(th0)<TINYLENGTH && i>0) th0_temp = dx;
		else th0_temp = th0;
		maxwell_test_map_radius(file_input, r0, rmax, r_step, th0_temp, thmax, th_step, z, dx,dy,dz, &(latt->cell[i]), add_field_comp_VFFA_rect_enge, doyoucirc);
	}
	sprintf(name_eps, "%s_div.eps", file_input);
	find_boundaries_heated_map_file(name_div, &xmin, &xmax, &ymin, &ymax, &valmin, &valmax);
	easyplot3dmap(name_div,   "x [m]", "y [m]", "div B [T/m]"     , "[0:2.]", "[9.7:10.2]", name_eps, "size ratio -1\nset grid", valmin, valmax);
	sprintf(name_eps, "%s_curlx.eps", file_input);
	find_boundaries_heated_map_file(name_curlx, &xmin, &xmax, &ymin, &ymax, &valmin, &valmax);
	easyplot3dmap(name_curlx, "x [m]", "y [m]", "(curl B)_x [T/m]", "[0:2.]", "[9.7:10.2]", name_eps, "size ratio -1\nset grid", valmin, valmax);
	sprintf(name_eps, "%s_curly.eps", file_input);
	find_boundaries_heated_map_file(name_curly, &xmin, &xmax, &ymin, &ymax, &valmin, &valmax);
	easyplot3dmap(name_curly, "x [m]", "y [m]", "(curl B)_y [T/m]", "[0:2.]", "[9.7:10.2]", name_eps, "size ratio -1\nset grid", valmin, valmax);
	sprintf(name_eps, "%s_curlz.eps", file_input);
	find_boundaries_heated_map_file(name_curlz, &xmin, &xmax, &ymin, &ymax, &valmin, &valmax);
	easyplot3dmap(name_curlz, "x [m]", "y [m]", "(curl B)_z [T/m]", "[0:2.]", "[9.7:10.2]", name_eps, "size ratio -1\nset grid", valmin, valmax);
	
}

extern void brute_force_search(struct Lattice *latt, struct Particle *reference, double xmin, double xmax, int nbx, double zmin, double zmax, int nbz, double uxmin, double uxmax, int nbux, double uzmin, double uzmax, int nbuz)
{
	int i,j,k,l;
	double xstep, uxstep, zstep, uzstep, x,z, ux,uz;
	struct Particle part;
	FILE *wfile;
	wfile=fopen("data/bruteforce.dat","w");
	
	xstep = (xmax-xmin)/(nbx);
	zstep = (zmax-zmin)/(nbz);
	uxstep = (uxmax-uxmin)/(nbux);
	uzstep = (uzmax-uzmin)/(nbuz);
	
	for(i=0;i<nbx;i++) {
		x = xmin + i*xstep;
		for(j=0;j<nbz;j++) {
			z = zmin + j*zstep;
			printf("i=%i, j=%i\n",i,j);
			for(k=0;k<nbux;k++) {
				ux = uxmin + k*uxstep;
				for(l=0;l<nbuz;l++) {
					uz = uzmin + l*uzstep;
					//printf("i=%i, j=%i, k=%i, l=%i\n",i,j,k,l);
					//printf("x=%le, z=%le, ux=%le, uz=%le\n",x,z,ux,uz);
					part = *reference;
					part.x = x;
					part.z = z;
					part.ux = ux;
					part.uz = uz;
					part.uy = sqrt(1.0 - (part.ux*part.ux + part.uz*part.uz));
					//part_cross_latt(&part,latt,NULL);
					part_oneturn(&part,latt,NULL);
					if(part.status==ALIVE) fprintf(wfile, "x=%le, z=%le, ux=%le, uz=%le\n",x,z,ux,uz);
				}
			}
		}
	}
	fclose(wfile);
}

extern void compare_field_point(char *file_dif_prefix, double x_1, double y_1, double z_1, struct Cell *cell1, void(*add_contribution_comp1)(double,double,double,double*,double*,double*,struct Cell*,int), 
double x_2, double y_2, double z_2, struct Cell *cell2, void(*add_contribution_comp2)(double,double,double,double*,double*,double*,struct Cell*,int), int doyouradial)
{
	char namex[300], namey[300], namez[300]; 
	double b1x,b1y,b1z,b2x,b2y,b2z;
	FILE *wfilex;
	FILE *wfiley;
	FILE *wfilez;
	
	sprintf(namex,  "%s_bx.dat", file_dif_prefix);
	sprintf(namey,  "%s_by.dat", file_dif_prefix);
	sprintf(namez,  "%s_bz.dat", file_dif_prefix);
	wfilex = fopen(namex, "a");
	wfiley = fopen(namey, "a");
	wfilez = fopen(namez, "a");
	
	if(get_bfield(x_1, y_1, z_1, &b1x, &b1y, &b1z, cell1, add_contribution_comp1)!=ALIVE) errorstop("field in cell1 not alive");
	if(get_bfield(x_2, y_2, z_2, &b2x, &b2y, &b2z, cell2, add_contribution_comp2)!=ALIVE) errorstop("field in cell2 not alive");
	
	fprintf(wfilex,"%lf	%lf	%le\n", y_1, x_1, b1x-b2x);
	fprintf(wfiley,"%lf	%lf	%le\n", y_1, x_1, b1y-b2y);
	fprintf(wfilez,"%lf	%lf	%le\n", y_1, x_1, b1z-b2z);
	
	if(doyouradial==YES) {
		char namer[300], nameth[300];
		double th, b1r, b1th, b2r, b2th;
		FILE *wfiler;
		FILE *wfileth;
		sprintf(namer, "%s_br.dat", file_dif_prefix);
		sprintf(nameth, "%s_bth.dat", file_dif_prefix);
		wfiler = fopen(namer, "a");
		wfileth = fopen(nameth, "a");
		th = atan_ratio(y_1, x_1);
		b1r = b1x*cos(th) + b1y*sin(th);
		b1th = -b1x*sin(th) + b1y*cos(th);
		b2r = b2x*cos(th) + b2y*sin(th);
		b2th = -b2x*sin(th) + b2y*cos(th);
		fprintf(wfiler, "%lf	%lf	%le\n", y_1, x_1, b1r-b2r);
		fprintf(wfileth,"%lf	%lf	%le\n", y_1, x_1, b1th-b2th);
		fclose(wfiler);
		fclose(wfileth);
	}

	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
}

extern void compare_field_cell_heatmap(char *file_dif_prefix, struct Cell *cell1, void(*add_contribution_comp1)(double,double,double,double*,double*,double*,struct Cell*,int), struct Cell *cell2, void(*add_contribution_comp2)(double,double,double,double*,double*,double*,struct Cell*,int), 
  double r0, double rmax, int nbr, double th0, double thmax, int nbth, double z0, double zmax, int nbz, int doyouradial)
{
	char namex[300], namey[300], namez[300], namer[300], nameth[300]; 
	int i=0;
	int j=0,k=0;
	double x,y,z,r,th, rstep,thstep,zstep;
	double b1x,b1y,b1z,b2x,b2y,b2z;
	FILE *wfilex;
	FILE *wfiley;
	FILE *wfiler;
	FILE *wfileth;
	FILE *wfilez;
	
	//doyou_long_boun_sym = NO;
	
	if(nbr <= 1 && nbz <= 1) errorstop("map not possible, only 1 dimension is scanned!");
	if(nbr != 1 && nbz != 1) errorstop("map not possible, 3 dimensions scanned");
	if(nbth <= 1) errorstop("map not possible, longitudinal dimension not scanned");
	
	sprintf(namex,  "%s_bx.dat", file_dif_prefix);
	sprintf(namey,  "%s_by.dat", file_dif_prefix);
	sprintf(namez,  "%s_bz.dat", file_dif_prefix);
	wfilex = fopen(namex, "w");
	wfiley = fopen(namey, "w");
	wfilez = fopen(namez, "w");
	if(doyouradial==YES) {
		sprintf(namer, "%s_br.dat", file_dif_prefix);
		sprintf(nameth, "%s_bth.dat", file_dif_prefix);
		wfiler = fopen(namer, "w");
		wfileth = fopen(nameth, "w");
	}
	
	zstep = comp_step(z0, zmax, nbz);
	rstep = comp_step(r0, rmax, nbr);
	thstep = comp_step(th0, thmax, nbth);
	
	//if(cell1->boun.thmax!=0) printf("rstep=%le[m], thstep=%le[deg], zstep=%le[m]\n",rstep,thstep*180./PI,zstep);
	if(doyouradial==YES) printf("rstep=%le[m], thstep=%le[deg], zstep=%le[m]\n",rstep,thstep*180./PI,zstep);
	else printf("xstep=%le[m], ystep=%le[m], zstep=%le[m]\n",rstep,thstep,zstep);
	
	for(j=0;j<nbr;j++) {
		r = r0 + j*rstep;
		for(k=0;k<nbth;k++) {
			th = th0+k*thstep;
			//th = thmax - k*thstep;
			for(i=0;i<nbz;i++) {
				z = z0 + i*zstep;
				//x = r*cos(th);
				//y = r*sin(th);
				x = r;
				y = th;
				//printf("j=%i, k=%i, (x,y,z): (%le,%le,%le)\n", j, k, x, y, z);
				if(get_bfield(x, y, z, &b1x, &b1y, &b1z, cell1, add_contribution_comp1)!=ALIVE) {
					//b1x=0; b1y=0; b1z=0;
					errorstop("field in cell1 not alive");
				}
				//printf("1: (%le,%le,%le),",b1x,b1y,b1z);
				if(get_bfield(x, y, z, &b2x, &b2y, &b2z, cell2, add_contribution_comp2)!=ALIVE) {
					//b2x=0; b2y=0; b2z=0;
					errorstop("field in cell2 not alive");
				}
				//printf(" 2: (%le,%le,%le)\n\n\n\n\n",b2x,b2y,b2z);
				if(nbz==1) {
					fprintf(wfilex,"%lf	%lf	%le\n", y, x, b1x*1.873670e-01-b2x);
					fprintf(wfiley,"%lf	%lf	%le\n", y, x, b1y*1.873670e-01-b2y);
					fprintf(wfilez,"%lf	%lf	%le\n", y, x, b1z*1.873670e-01-b2z);
				}
				else if(nbr==1) {
					fprintf(wfilex,"%lf	%lf	%le\n", y, z, b1x*1.873670e-01-b2x);
					fprintf(wfiley,"%lf	%lf	%le\n", y, z, b1y*1.873670e-01-b2y);
					fprintf(wfilez,"%lf	%lf	%le\n", y, z, b1z*1.873670e-01-b2z);
				}
				else errorstop("wrong, nbx and nbz");
				
				if(doyouradial==YES) {
					double th, b1r, b1th, b2r, b2th;
					th = atan_ratio(y, x);
					b1r = b1x*cos(th) + b1y*sin(th);
					b1th = -b1x*sin(th) + b1y*cos(th);
					b2r = b2x*cos(th) + b2y*sin(th);
					b2th = -b2x*sin(th) + b2y*cos(th);
					fprintf(wfiler, "%lf	%lf	%le\n", y, x, b1r-b2r);
					fprintf(wfileth,"%lf	%lf	%le\n", y, x, b1th-b2th);
				}
			}
			if(nbr==1) {
				fprintf(wfilex,"\n"); // to build file to plot easyplot3dmap
				fprintf(wfiley,"\n"); // to build file to plot easyplot3dmap
				fprintf(wfilez,"\n"); // to build file to plot easyplot3dma
			}
			
		}
		fprintf(wfilex,"\n"); // to build file to plot easyplot3dmap
		fprintf(wfiley,"\n"); // to build file to plot easyplot3dmap
		fprintf(wfilez,"\n"); // to build file to plot easyplot3dmap
		if(doyouradial==YES) {
			fprintf(wfiler,"\n"); // to build file to plot easyplot3dmap
			fprintf(wfileth,"\n"); // to build file to plot easyplot3dmap
		}
	}
	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
	if(doyouradial==YES) {
		fclose(wfiler);
		fclose(wfileth);
	}
}

extern void find_boundaries_2d_file(char *file, double *xmin, double *xmax, double *ymin, double *ymax)
{
	int i,nblines;
	double x, y;
	FILE *rfile=NULL;
	
	nblines = get_nb_lines_file(file);
	
	rfile = fopen(file, "r");
	if(rfile==NULL) errorstop("cannot open rfile\n");
	
	*xmin = 1.e9;
	*xmax = -1.e9;
	*ymin = 1.e9;
	*ymax = -1.e9;
	for(i=0;i<nblines;i++) {
		fscanf(rfile, "%le	%le", &x,&y);
		if(x==x && y==y) {
			if(x<*xmin) *xmin = x;
			if(x>*xmax) *xmax = x;
			if(y<*ymin) *ymin = y;
			if(y>*ymax) *ymax = y;
		}
	}
	fclose(rfile);
	printf("in %s:\n xrange:[%lf:%lf], yrange:[%lf:%lf]\n", file, *xmin, *xmax, *ymin, *ymax);
}

extern void find_boundaries_heated_map_file(char *file, double *xmin, double *xmax, double *ymin, double *ymax, double *val_min, double *val_max)
{
	int i,nblines;
	double x, y, val, x_val_min, y_val_min, x_val_max, y_val_max;
	FILE *rfile=NULL;
	
	nblines = get_nb_lines_file(file);
	
	rfile = fopen(file, "r");
	if(rfile==NULL) errorstop("cannot open rfile\n");
	
	*xmin = 1.e9;
	*xmax = -1.e9;
	*ymin = 1.e9;
	*ymax = -1.e9;
	*val_min = 1.e9;
	*val_max = -1.e9;
	for(i=0;i<nblines;i++) {
		fscanf(rfile, "%le	%le	%le", &x,&y,&val);
		if(x==x && y==y && val==val) {
			if(x<*xmin) *xmin = x;
			if(x>*xmax) *xmax = x;
			if(y<*ymin) *ymin = y;
			if(y>*ymax) *ymax = y;
			if(val<*val_min) {
				*val_min = val;
				x_val_min = x;
				y_val_min = y;
			}
			if(val>*val_max) {
				*val_max = val;
				x_val_max = x;
				y_val_max = y;
			}
		}
	}
	fclose(rfile);
	printf("in %s:\n \tmin value %le at (%le,%le)\n \tmax value %le at (%le,%le)\n",file,*val_min,x_val_min,y_val_min, *val_max, x_val_max, y_val_max);
	printf(" xrange:[%lf:%lf], yrange:[%lf:%lf]\n", *xmin, *xmax, *ymin, *ymax);
	
}

extern void find_boundaries_heated_map_file2(char *file, char *xrange, char *yrange, double *val_min, double *val_max)
{
	double xmin, xmax, ymin, ymax;
	
	find_boundaries_heated_map_file(file, &xmin, &xmax, &ymin, &ymax, val_min, val_max);
	sprintf(xrange,"[%lf:%lf]", xmin, xmax);
	sprintf(yrange,"[%lf:%lf]", ymin, ymax);
	//printf("xrange:[%lf:%lf]", xmin, xmax);
	//printf("yrange:[%lf:%lf]", ymin, ymax);
}

// compute field map only in the convergence radius of cell2
extern void compute_fieldmap_special(struct Cell *cell2, char *filename, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int),
							double r0, double rmax, int nbr, double th0, double thmax, int nbth, double z0, double zmax, int nbz)
{
	char namex[300], namey[300], namez[300], outx[300], outy[300], outz[300];
	int i,j,k, ic, circ=YES;
	double x,y,z,bx,by,bz,r,th, rstep,thstep,zstep;
	double xmin,xmax,ymin,ymax,bmin,bmax;
	double x_0, y_0, newx;
	//FILE *wfile;
	FILE *wfilex, *wfiley, *wfilez;

//	wfile = fopen(filename,"w");
//	fprintf(wfile,"%i     %i     %i\n\n\n\n\n\n\n\n",nbr,nbz,nbth);
	sprintf(namex,"%s_bx.dat", filename);
	wfilex=fopen(namex, "w");
	sprintf(namey,"%s_by.dat", filename);
	wfiley=fopen(namey, "w");
	sprintf(namez,"%s_bz.dat", filename);
	wfilez=fopen(namez, "w");

	zstep = comp_step(z0, zmax, nbz);
	rstep = comp_step(r0, rmax, nbr);
	thstep = comp_step(th0, thmax, nbth);
	
	printf("rstep=%le, thstep=%le, zstep = %le\n", rstep,thstep,zstep);
	
	for(i=0;i<nbz+1;i++) {
		z = z0 + i*zstep;
		for(j=0;j<nbr+1;j++) {
			r = r0 + j*rstep;
			for(k=0;k<nbth+1;k++) {
				th = k*thstep;
				if(circ == YES) {
					x = r*cos(th);
					y = r*sin(th);
				}
				else {
					x = r;
					y = th;
				}
				if(th < 6*PI/180.) ic = 0;
				else if(th < 18*PI/180.) ic = 1;
				else ic = 2;
				x_0 = cell2->mpara[ic][1]*cos(cell2->mpara[ic][0]);
				y_0 = cell2->mpara[ic][1]*sin(cell2->mpara[ic][0]);
				newx = (y-y_0)*sin(cell2->mpara[ic][0]+cell2->mpara[ic][5]) + (x-x_0)*cos(cell2->mpara[ic][0]+cell2->mpara[ic][5]);
				if(fabs(newx)<cell2->mpara[ic][6]) {
					if(get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp)!=ALIVE) errorstop("bfield not alive");
				}
				else {
					bx = 0;
					by = 0;
					bz = 0;
				}
				
				fprintf(wfilex,"%le\t%le\t%le\n", y, x, bx);
				fprintf(wfiley,"%le\t%le\t%le\n", y, x, by);
				fprintf(wfilez,"%le\t%le\t%le\n", y, x, bz);
			}
			fprintf(wfilex,"\n"); // to build file to plot easyplot3dmap
			fprintf(wfiley,"\n"); // to build file to plot easyplot3dmap
			fprintf(wfilez,"\n"); // to build file to plot easyplot3dmap
		}
	}
	//fclose(wfile);
	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
	
	sprintf(outx,"%s_bx.eps", filename);
	sprintf(outy,"%s_by.eps", filename);
	sprintf(outz,"%s_bz.eps", filename);
	
	//easyplot3dmap(namex, "x [m]", "y [m]", "[0:2]", NULL, NULL, outx, "size ratio -1\nset grid\nset ytics 0.2", bxmin,bxmax);
	//easyplot3dmap(namey, "x [m]", "y [m]", "[0:2]", NULL, NULL, outy, "size ratio -1\nset grid\nset ytics 0.2", bymin,bymax);
	//easyplot3dmap(namez, "x [m]", "y [m]", "[0:2]", NULL, NULL, outz, "size ratio -1\nset grid\nset ytics 0.2", bzmin,bzmax);
	find_boundaries_heated_map_file(namex, &xmin, &xmax, &ymin, &ymax, &bmin, &bmax);
	easyplot3dmap(namex, "x [m]", "y [m]", NULL, NULL, NULL, outx, "size ratio -1\nset grid",bmin, bmax);// bxmin,bxmax);
	find_boundaries_heated_map_file(namey, &xmin, &xmax, &ymin, &ymax, &bmin, &bmax);
	easyplot3dmap(namey, "x [m]", "y [m]", NULL, NULL, NULL, outy, "size ratio -1\nset grid",bmin, bmax);// bymin,bymax);
	find_boundaries_heated_map_file(namez, &xmin, &xmax, &ymin, &ymax, &bmin, &bmax);
	easyplot3dmap(namez, "x [m]", "y [m]", NULL, NULL, NULL, outz, "size ratio -1\nset grid",bmin, bmax);// bzmin,bzmax);
	
	//easyplot3dmap_2d(namex, "data/aspect.dat", "using 2:1 with lines lt 1 lc 0 lw 2", "x [m]", "y [m]", NULL, NULL, outx, "size ratio -1\nset grid", bxmin,bxmax);
	//easyplot3dmap_2d(namey, "data/aspect.dat", "using 2:1 with lines lt 1 lc 0 lw 2", "x [m]", "y [m]", NULL, NULL, outy, "size ratio -1\nset grid", bymin,bymax);
	//easyplot3dmap_2d(namez, "data/aspect.dat", "using 2:1 with lines lt 1 lc 0 lw 2", "x [m]", "y [m]", NULL, NULL, outz, "size ratio -1\nset grid", bzmin,bzmax);
}

extern void file_vffa_btemp(char *trackout)
{
	char buf[MAX_CHARINLINE];
	double x_temp,y_temp,z_temp,bxtemp,bytemp,bztemp;
	double x_track,y_track,z_track,s,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10;
	FILE *wfile;
	FILE *rfile_temp=NULL;
	FILE *rfile_track=NULL;
	
	printf("enter file_vffa_btemp for %s", trackout);
	rfile_temp = fopen("data/field_vffa_loc.dat", "r");
	if(rfile_temp==NULL) errorstop("cannot open rfile_temp");
	rfile_track = fopen(trackout, "r");
	if(rfile_track==NULL) errorstop("cannot open rfile_track");
	
	wfile = fopen("data/field_vffa_loc_s.dat", "w");
	
	while(!feof(rfile_track)) {
		if(fgets(buf, MAX_CHARINLINE-1, rfile_track) == NULL) break;
		buf[MAX_CHARINLINE-1]='\0';
		if(sscanf(buf, "%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %lf\t %le", &s, &x_track, &y_track, &z_track, &b1, &b2, &b3, &b4, &b5, &b6, &b7, &b8, &b9, &b10) != 14) continue;
		printf("cherche : %le	%le	%le\n", x_track, y_track, z_track);
		do {
			fscanf(rfile_temp, "%le	%le	%le	%le	%le	%le", &x_temp,&y_temp,&z_temp, &bxtemp,&bytemp,&bztemp);
			printf("%le	%le	%le	%le	%le	%le\n", x_temp,y_temp,z_temp,bxtemp,bytemp,bztemp);
		} while(fabs(x_temp-x_track)>TINYLENGTH && fabs(y_temp-y_track)>TINYLENGTH && fabs(z_temp-z_track)>TINYLENGTH);
		printf("OK!\n");
		fprintf(wfile, "%le	%le	%le	%le	%le	%le	%le\n", s, x_track, y_track, z_track, bxtemp, bytemp, bztemp);
	}
	
	fclose(wfile);
	fclose(rfile_temp);
	fclose(rfile_track);
}

/*
extern void field_polygone_shift(char *file_polygone, char *trackout, double angle1, double angle2, double angle_boun)
{
	int i, nblines;
	double s, x, y, z, bx, by, bz, ux, uy, uz, brho, t, r, th, angle, bxp, byp;
	FILE *rfile = NULL;
	FILE *wfile;
	nblines = get_nb_lines_file(trackout);
	rfile = fopen(trackout, "r");
	if(rfile==NULL) errorstop("cannot open trackout");
	wfile = fopen(file_polygone, "w");
	for(i=0;i<nblines;i++) {
		fscanf(wfile, "%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %.10f\t %le", &s, &x, &y, &z, &bx, &by, &bz, &ux, &uy, &uz, &brho, &t, &r, &th);
		if(th<angle_boun) angle = angle1;
		else angle = angle2;
		bxp = by*sin(angle) + bx*cos(angle);
		byp = by*cos(angle) - bx*sin(angle);
		fprintf(wfile, "\n", bxp, byp, bz);
	}
	
	fclose(rfile);
	fclose(wfile);
}
//*/

extern void betafunc_period_vffa(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Particle *reference, struct Lattice *latt)
{
	int i,n;
	double betax, alphax, betaz, alphaz, step_th, qx,qz;
	struct Particle test_part;//, test_part2;
	struct Lattice latt_1cell;
	FILE *wfile;
	//print_cell_para("latt", &(latt->cell[0]));
	
	if(nsteps_cell <= 1) errorstop("!!!ERROR in betafunc_period_vffa: nsteps_cell <= 1!!!");
	//print_cell_para("latt_copy", &(latt_copy.cell[0]));
	gene_latt_1cell(&latt_1cell, &(latt->cell[0])); //copy cell in latt_1cell (no memory allocation)
	
	//find twiss parameters at cup position in the cell
	latt_1cell.cell[0].instrutype = CUP;
	latt_1cell.cell[0].instru.ymax = 0;
	latt_1cell.cell[0].instru.thmax = 0;
	//step_th = latt->cell[0].boun.thmax/(nsteps_cell-1.);
	step_th = comp_step(0, latt->cell[0].boun.thmax, nsteps_cell-1);
	printf("step_th = %le\n", step_th);
	print_cell_para(NULL, "latt_1cell", &(latt_1cell.cell[0]));
	
	
	wfile = fopen(filename, "w");
	for(i = 0; i < nsteps_cell ; i++) {
		//CLRSCR();
		printf("step %i/%i\n", i,nsteps_cell-1);
		//fflush(stdout);
		latt_1cell.cell[0].instrutype = CUP;
		latt_1cell.cell[0].instru.ymax = 0;
		latt_1cell.cell[0].instru.thmax = i*step_th;
		//print_cell_para("latt_1cell", &(latt_1cell.cell[0]));
		test_part = *reference;
		//test_part2 = *reference;
		//print_part_para("\ntest_part", &test_part);
		//printf("\nx_clo = %.10f [m], xprime_clo = %.10e [deg], \nz_clo = %.10f [m], zprime_clo = %.10e [deg]\n", test_part.x, atan(test_part.ux/(test_part.uy))*180./(PI), test_part.z, atan(test_part.uz/(sqrt((test_part.uy)*(test_part.uy)+(test_part.ux)*(test_part.ux))))*180./(PI));
		part_cross_latt(&test_part, &latt_1cell, NULL);
		//print_part_para("test_part apres", &test_part);
		//printf("\nx_clo = %.10f [m], xprime_clo = %.10e [deg], \nz_clo = %.10f [m], zprime_clo = %.10e [deg]\n", test_part.x, atan(test_part.ux/(test_part.uy))*180./(PI), test_part.z, atan(test_part.uz/(sqrt((test_part.uy)*(test_part.uy)+(test_part.ux)*(test_part.ux))))*180./(PI));
		//adjust_particle(&test_part2, &test_part, latt_1cell.cell[0].instru.thmax);
		for(n=0; n<latt_1cell.cell[0].nbcomp; n++) latt_1cell.cell[0].mpara[n][0] -= i*step_th;
		latt_1cell.cell[0].instrutype = NO;
		//print_part_para("test_part2", &test_part2);
		if(put_on_co_nelmin( &latt_1cell, &test_part, 1.e-10, 1.e-3, 1.e-3, 1.e-3, 1.e-3, YES)!= TRUE) errorstop("closed orbit nelmin not found\n");
		if(find_closed_orbite_xxp_zzp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), &(test_part.z), &(test_part.uz), 1.e-10, &latt_1cell, YES) != TRUE) errorstop("closed orbit not found\n");
		//printf("\nx_clo = %.10f [m], xprime_clo = %.10e [deg], \nz_clo = %.10f [m], zprime_clo = %.10e [deg]\n", test_part.x, atan(test_part.ux/(test_part.uy))*180./(PI), test_part.z, atan(test_part.uz/(sqrt((test_part.uy)*(test_part.uy)+(test_part.ux)*(test_part.ux))))*180./(PI));
		
		//part_cross_latt(&test_part2, &latt_copy, "data/betafunc_period_vffa_cross_latt.dat");
		//get_periodic_twiss(&betax, &alphax, &betaz, &alphaz, amp_x, amp_xprime, amp_z, amp_zprime, &test_part, &latt_1cell, 1);
		calc_tune_twiss(&test_part, &qx, &qz, &betax, &alphax, &betaz, &alphaz, amp_x, amp_xprime, amp_z, amp_zprime, &latt_1cell, part_cross_latt, YES, "data/tune_beta_vffa.dat",1);
		fprintf(wfile, "%le  %le  %le  %le  %le\n", test_part.s, betax, alphax, betaz, alphaz);
		for(n=0; n<latt_1cell.cell[0].nbcomp; n++) latt_1cell.cell[0].mpara[n][0] += i*step_th;
		
	}
	fclose(wfile);
	//easyplot("data/track_beta_instru.dat", "($3)", "($2)", "lines lc 7 lw 2", "y [m]", "x [m]", NULL, NULL, "output/traj_beta_xy.eps", "mxtics 5\nset mytics 5");
	//easyplot("data/track_beta_instru.dat", "($1)", "($4)", "lines lc 7 lw 2", "s [m]", "z [m]", NULL, NULL, "output/traj_beta_sz.eps", "mxtics 5\nset mytics 5");
	
}//*/

extern void adjust_particle(struct Particle *part_to_change, struct Particle *part_ori, double angle_rad)
{
	*part_to_change = *part_ori;
	part_to_change->fwk.ae -= angle_rad;
	part_to_change->s = 0;
	part_to_change->x = sqrt(part_ori->x*part_ori->x+part_ori->y*part_ori->y);
	part_to_change->y = 0;
	part_to_change->z = part_ori->z;
	part_to_change->ux = part_ori->ux;//*cos(angle_rad) + part_ori->uy*sin(angle_rad);
	part_to_change->uy = part_ori->uy;//*cos(angle_rad) - part_ori->ux*sin(angle_rad);
	part_to_change->uz = part_ori->uz;
}

extern void compute_map_polar(char *filename, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int),
	double rstart, double rend, int nbr, double rshift, double thstart, double thend, int nbth, double thshift, double zstart, double zend, int nbz, double zshift)
{
	int ir,ith,iz;
	double x,y,z,bx,by,bz,r,th, br,bth, rstep,thstep,zstep;
	
	FILE *wfile;
	wfile = fopen(filename,"w");
	fprintf(wfile, "%i %i %i\n", nbr, nbth, nbz);
	fprintf(wfile,"%i steps in r, from %le to %le (unit: m)\n",nbr-1, rstart+rshift, rend+rshift);
	fprintf(wfile,"%i steps in theta, from %le to %le (unit: deg)\n",nbth-1, thstart+thshift, thend+thshift);
	fprintf(wfile,"%i steps in z, from %le to %le (unit: m)\n",nbz-1, zstart+zshift, zend+zshift);
	fprintf(wfile, "magnetic field unit: T\n");
	fprintf(wfile, "line order: r (horizontal radius), theta (azimuthal angle), z (vertical), B_r, B_theta, B_z\n");
	fprintf(wfile, "loop iteration order: z, theta, r\n\n");

	zstep = comp_step(zstart, zend, nbz);
	rstep = comp_step(rstart, rend, nbr);
	thstep = comp_step(thstart, thend, nbth);
	
	for(ir = 0; ir < nbr; ir++) {
		for(ith = 0; ith < nbth; ith++) {
			for(iz = 0; iz < nbz; iz++) {
				r = rstart + ir*rstep;
				th = thstart + ith*thstep;
				z = zstart + iz*zstep;
				x = r*cos(th);
				y = r*sin(th);
				printf("r=%lf, th=%le, z=%le\n",r,th*180/PI,z);
				if(get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp)==ALIVE) {
					br = bx*cos(th) + by*sin(th);
					bth = -br*sin(th) + by*cos(th);
					r += rshift;
					th += thshift;
					th *= 180./(PI);
					z += zshift;
					fprintf(wfile,"%lf	%lf	%lf	%le	%le	%le\n", r, th, z, br, bth, bz);
				}
				else errorstop("bfield not alive");
			}
		}
	}
	fclose(wfile);
}

extern void compute_fieldheatmap(char *filename_prefix, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int),
							double r0, double rmax, int nbr, double th0, double thmax, int nbth, double z0, double zmax, int nbz, int doyouradial)
{
	char namex[300], namey[300], namez[300];
	int i,j,k;
	double x,y,z,bxtemp,bytemp,bx,by,bz,r,th, rstep,thstep,zstep, w_coor;
	FILE *wfilex, *wfiley, *wfilez;
	
	if(doyouradial==YES) {
		sprintf(namex,"%s_bxshinji.dat", filename_prefix);
		sprintf(namey,"%s_byshinji.dat", filename_prefix);
		//sprintf(namex,"%s_br.dat", filename_prefix);
		//sprintf(namey,"%s_bth.dat", filename_prefix);
	}
	else {
		sprintf(namex,"%s_bx.dat", filename_prefix);
		sprintf(namey,"%s_by.dat", filename_prefix);
	}
	wfilex=fopen(namex, "w");
	wfiley=fopen(namey, "w");
	sprintf(namez,"%s_bz.dat", filename_prefix);
	wfilez=fopen(namez, "w");
	
	if(nbr <= 1 && nbz <= 1) errorstop("map not possible, only 1 dimension is scanned!");
	if(nbr != 1 && nbz != 1 && nbth != 1) errorstop("map not possible, 3 dimensions scanned");
	//if(nbth <= 1) errorstop("map not possible, longitudinal dimension not scanned");
	
	zstep = comp_step(z0, zmax, nbz);
	rstep = comp_step(r0, rmax, nbr);
	thstep = comp_step(th0, thmax, nbth);
	
	if(doyouradial==YES) printf("rstep=%le[m], thstep=%le[deg], zstep=%le[m]\n",rstep,thstep*180./PI,zstep);
	else printf("xstep=%le[m], ystep=%le[m], zstep=%le[m]\n",rstep,thstep,zstep);
	
	for(j=0;j<nbr;j++) {
		r = r0 + j*rstep;
		for(i=0;i<nbz;i++) {
			z = z0 + i*zstep;
			for(k=0;k<nbth;k++) {
				th = k*thstep;
				//th = thmax - k*thstep;
				if(doyouradial==YES) {
					x = r*cos(th);
					y = r*sin(th);
				}
				else {
					x = r;
					y = th;
				}
				if(get_bfield(x, y, z, &bxtemp, &bytemp, &bz, cell, add_contribution_comp)!=ALIVE) {
					//errorstop("field not alive\n");
					bx = 0.;
					by = 0.;
					bz = 0.;
				}
				if(nbz==1) w_coor = x;
				else w_coor = z;
				if(doyouradial==YES) {
					bx = bxtemp*cos(18*PI/180.) + bytemp*sin(18*PI/180.); //br
					by = -bxtemp*sin(18*PI/180.) + bytemp*cos(18*PI/180.); 
					//bx = bxtemp*cos(th) + bytemp*sin(th); //br
					//by = -bxtemp*sin(th) + bytemp*cos(th); //bth
				}
				else {
					bx = bxtemp*cos(18*PI/180.) + bytemp*sin(18*PI/180.); //br
					by = -bxtemp*sin(18*PI/180.) + bytemp*cos(18*PI/180.); 
					//bx = bxtemp;
					//by = bytemp;
				}
				//printf("(%lf, %lf, %lf)\n",x,y,z);
				fprintf(wfilex,"%le\t%le\t%le\n", y, w_coor, bx);
				fprintf(wfiley,"%le\t%le\t%le\n", y, w_coor, by);
				fprintf(wfilez,"%le\t%le\t%le\n", y, w_coor, bz);
			}
			fprintf(wfilex,"\n"); // to build file to plot easyplot3dmap
			fprintf(wfiley,"\n"); // to build file to plot easyplot3dmap
			fprintf(wfilez,"\n"); // to build file to plot easyplot3dmap
		}
	}
	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
}

extern void compute_fieldheatmap_cart_hor_vert(char *filename_prefix, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int),
							double x0, double xmax, int nbx, double y0, double z0, double zmax, int nbz)
{
	char namex[300], namey[300], namez[300];
	int i,j,k;
	double x,y,z,bx,by,bz, xstep,zstep;
	FILE *wfilex, *wfiley, *wfilez;
	
	sprintf(namex,"%s_bx.dat", filename_prefix);
	wfilex=fopen(namex, "w");
	sprintf(namey,"%s_by.dat", filename_prefix);
	wfiley=fopen(namey, "w");
	sprintf(namez,"%s_bz.dat", filename_prefix);
	wfilez=fopen(namez, "w");
	printf("nbx = %i, nbz=%i\n",nbx,nbz);
	if(nbx <= 1 || nbz <= 1) errorstop("map not possible, only 1 dimension is scanned!");
	
	zstep = comp_step(z0, zmax, nbz);
	xstep = comp_step(x0, xmax, nbx);
	
	printf("xstep=%le[m], zstep=%le[m]\n",xstep,zstep);
	
	for(j=0;j<nbx;j++) {
		x = x0 + j*xstep;
		for(i=0;i<nbz;i++) {
			z = z0 + i*zstep;
			if(get_bfield(x, y0, z, &bx, &by, &bz, cell, add_contribution_comp)!=ALIVE) {
				//errorstop("field not alive\n");
				bx = 0.;
				by = 0.;
				bz = 0.;
			}
			fprintf(wfilex,"%le\t%le\t%le\n", x, z, bx);
			fprintf(wfiley,"%le\t%le\t%le\n", x, z, by);
			fprintf(wfilez,"%le\t%le\t%le\n", x, z, bz);
		}
		fprintf(wfilex,"\n"); // to build file to plot easyplot3dmap
		fprintf(wfiley,"\n"); // to build file to plot easyplot3dmap
		fprintf(wfilez,"\n"); // to build file to plot easyplot3dmap
	}
	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
}

extern void scan_vffa(char *filename, double m_step, int nb_m, double b0d_step, int nb_b0d, double angle_step, int nb_angle, double r0d_step, int nb_r0d, struct Lattice *latt, struct Particle *part)
{
	int i, j, i2, i3, i4;
	double b0d, angle, r0d;
	struct Particle test_part, part_m, part_b0, part_r0;
	
	test_part = *part;
	part_m = test_part;
	b0d = latt->cell[0].mpara[1][2];
	angle = latt->cell[0].mpara[0][5];
	r0d = latt->cell[0].mpara[1][1];
	for(i=0;i<nb_m;i++) {
		latt->cell[0].mpara[1][1] = r0d;
		latt->cell[0].mpara[0][5] = angle;
		latt->cell[0].mpara[2][5] = -angle;
		latt->cell[0].mpara[1][2] = b0d;
		test_part = part_m;
		
		for(j=0;j<latt->cell[0].nbcomp;j++) latt->cell[0].mpara[j][4] += m_step; // change m
		put_on_co_nelmin(latt, &(test_part), 1.e-6, 1.e-4, 1.e-4, 1.e-4, 1.e-4, YES);
		if(find_closed_orbite_xxp_zzp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), &(test_part.z), &(test_part.uz), 1.e-8, latt, YES) == TRUE) print_cell_para(filename, "cell", &(latt->cell[0]));
			part_m = test_part;
		part_b0 = part_m;
		for(i2=1;i2<nb_b0d;i2++) {
			latt->cell[0].mpara[1][1] = r0d;
			latt->cell[0].mpara[0][5] = angle;
			latt->cell[0].mpara[2][5] = -angle;
			test_part = part_b0;
			
			latt->cell[0].mpara[1][2] += b0d_step; // change b0d
			put_on_co_nelmin(latt, &(test_part), 1.e-6, 1.e-4, 1.e-4, 1.e-4, 1.e-4, YES);
			if(find_closed_orbite_xxp_zzp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), &(test_part.z), &(test_part.uz), 1.e-8, latt, YES) == TRUE) print_cell_para(filename, "cell", &(latt->cell[0]));
			part_b0 = test_part;
			for(i3=1;i3<nb_angle;i3++) {
				latt->cell[0].mpara[1][1] = r0d;
				test_part = part_r0;
				
				latt->cell[0].mpara[0][5] += angle_step*PI/180.; //change F1 tilt angle
				latt->cell[0].mpara[2][5] -= angle_step*PI/180.; //change F2 tilt angle
				put_on_co_nelmin(latt, &(test_part), 1.e-6, 1.e-4, 1.e-4, 1.e-4, 1.e-4, YES);
				if(find_closed_orbite_xxp_zzp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), &(test_part.z), &(test_part.uz), 1.e-8, latt, YES) == TRUE) print_cell_para(filename, "cell", &(latt->cell[0]));
					part_r0 = test_part;
				for(i4=1;i4<nb_r0d;i4++) {
					latt->cell[0].mpara[1][1] += r0d_step; // change r0d
					put_on_co_nelmin(latt, &(test_part), 1.e-6, 1.e-4, 1.e-4, 1.e-4, 1.e-4, YES);
					if(find_closed_orbite_xxp_zzp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), &(test_part.z), &(test_part.uz), 1.e-8, latt, YES) == TRUE) print_cell_para(filename, "cell", &(latt->cell[0]));
				}
			}
		}
	}
}

// compute acceptance in (u,u') and (v,v') phase spaces with an initial amplitude in u
extern int acceptanceu(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double aumin, double aumax, double austep, double av)
{
	int i;
	double m4d[4][4], decoup_m[4][4], decoup_to_coup[4][4], decoup_to_coup_inv[4][4], u_vect[4], x_vect[4], output[6];
	double php_ref, xprime_ref, zprime_ref, band_bias, amp_u;
	struct Particle test_part;
	FILE *rfile;
	if(outfilename != NULL) rfile = fopen(outfilename,"a");
	
	printf("acceptance in u:");
	php_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	zprime_ref  = atan_ratio(reference->uz, php_ref);
	x_vect[0] = 1.e-5;
	x_vect[1] = 1.e-5;
	x_vect[2] = 1.e-5;
	x_vect[3] = 1.e-5;
	if(get_matrix_firstorder_4d(m4d, &band_bias, reference, x_vect[0], x_vect[1], x_vect[2], x_vect[3], latt, part_cross_latt, YES) != TRUE) return FALSE;
	decouple_matrix_parzen(m4d, decoup_m, decoup_to_coup, output, NULL);
	matrix_inverse_4d(decoup_to_coup, decoup_to_coup_inv);
	
	for(amp_u = aumin; amp_u <= aumax; amp_u += austep) {
		u_vect[0] = amp_u; // u
		u_vect[1] = 0; //uprime
		u_vect[2] = av; //v
		u_vect[3] = 0; //vprime
		mvprod4(x_vect, decoup_to_coup, u_vect);
		test_part = *reference;
		test_part.s = 0;
		test_part.x += x_vect[0]/cos(xprime_ref);
		test_part.z += x_vect[2]/cos(zprime_ref);
		test_part.ux = sin(xprime_ref + x_vect[1])*cos(zprime_ref + x_vect[3]);
		test_part.uy = cos(xprime_ref + x_vect[1])*cos(zprime_ref + x_vect[3]);
		test_part.uz = sin(zprime_ref + x_vect[3]);
		test_part.hat = -2;
		printf("\namp = %lf\n", amp_u);
		
		for(i = 0; i < nbpass; i++) {
			CLRSCR();
			printf("pass number %i\t", i);
			fflush(stdout);
			if (test_part.uy > 0) {
				comp_phase_space_coord(&(x_vect[0]), &(x_vect[1]), &(x_vect[2]), &(x_vect[3]), test_part.x, test_part.z, test_part.ux, test_part.uy, test_part.uz, reference->x, xprime_ref, reference->z, zprime_ref);
				mvprod4(u_vect, decoup_to_coup_inv, x_vect);
			}
			else {
				printf("\tin acceptanceu, test_part.uy <= 0!, particle is going backwards!\n");
				if(outfilename != NULL) fclose(rfile);
				return test_part.status;
			}
			if(outfilename != NULL) fprintf(rfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le\n", x_vect[0], x_vect[1], x_vect[2], x_vect[3], u_vect[0], u_vect[1], u_vect[2], u_vect[3], amp_u);
			part_cross_latt(&(test_part), latt, NULL);
			if(test_part.status != ALIVE) {
				if(outfilename != NULL) fclose(rfile);
				printf("\n");
				return test_part.status;
			}
		}
	}	
	
	if(outfilename != NULL) fclose(rfile);
	printf("\n");
	return test_part.status;
}

// compute acceptance in (u,u') and (v,v') phase spaces with an initial amplitude in u
extern int acceptancev(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double avmin, double avmax, double avstep, double au)
{
	int i;
	double m4d[4][4], decoup_m[4][4], decoup_to_coup[4][4], decoup_to_coup_inv[4][4], u_vect[4], x_vect[4], output[6];
	double php_ref, xprime_ref, zprime_ref, band_bias, amp_v;
	struct Particle test_part;
	FILE *rfile;
	if(outfilename != NULL) rfile = fopen(outfilename,"a");

	printf("acceptance in v:");
	php_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	zprime_ref  = atan_ratio(reference->uz, php_ref);
	x_vect[0] = 1.e-5;
	x_vect[1] = 1.e-5;
	x_vect[2] = 1.e-5;
	x_vect[3] = 1.e-5;
	if(get_matrix_firstorder_4d(m4d, &band_bias, reference, x_vect[0], x_vect[1], x_vect[2], x_vect[3], latt, part_cross_latt, YES) != TRUE) return FALSE;
	decouple_matrix_parzen(m4d, decoup_m, decoup_to_coup, output, NULL);
	matrix_inverse_4d(decoup_to_coup, decoup_to_coup_inv);
	
	for(amp_v = avmin; amp_v <= avmax; amp_v += avstep) {
		u_vect[0] = au; // u
		u_vect[1] = 0; //uprime
		u_vect[2] = amp_v; //v
		u_vect[3] = 0; //vprime
		mvprod4(x_vect, decoup_to_coup, u_vect);
		test_part = *reference;
		test_part.s = 0;
		test_part.x += x_vect[0]/cos(xprime_ref);
		test_part.z += x_vect[2]/cos(zprime_ref);
		test_part.ux = sin(xprime_ref + x_vect[1])*cos(zprime_ref + x_vect[3]);
		test_part.uy = cos(xprime_ref + x_vect[1])*cos(zprime_ref + x_vect[3]);
		test_part.uz = sin(zprime_ref + x_vect[3]);
		test_part.hat = -2;
		printf("\namp = %lf\n", amp_v);
		
		for(i = 0; i < nbpass; i++) {
			CLRSCR();
			printf("pass number %i\t", i);
			fflush(stdout);
			if (test_part.uy > 0) {
				comp_phase_space_coord(&(x_vect[0]), &(x_vect[1]), &(x_vect[2]), &(x_vect[3]), test_part.x, test_part.z, test_part.ux, test_part.uy, test_part.uz, reference->x, xprime_ref, reference->z, zprime_ref);
				mvprod4(u_vect, decoup_to_coup_inv, x_vect);
			}
			else {
				printf("\tin acceptancev, test_part.uy <= 0!, particle is going backwards!\n");
				if(outfilename != NULL) fclose(rfile);
				return test_part.status;
			}
			if(outfilename != NULL) fprintf(rfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le\n", x_vect[0], x_vect[1], x_vect[2], x_vect[3], u_vect[0], u_vect[1], u_vect[2], u_vect[3], amp_v);
			part_cross_latt(&(test_part), latt, NULL);
			if(test_part.status != ALIVE) {
				if(outfilename != NULL) fclose(rfile);
				printf("\n");
				return test_part.status;
			}
		}
	}	
	
	if(outfilename != NULL) fclose(rfile);
	printf("\n");
	return test_part.status;
}

extern int acceptanceu_auto(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double aumin, double *aumax, double austep, double av) 
{
	int i, n;
	double m4d[4][4], decoup_m[4][4], decoup_to_coup[4][4], decoup_to_coup_inv[4][4], u_vect[4], x_vect[4], output[6];
	double au, band_bias, xpart[nbpass], zpart[nbpass], ux[nbpass], uy[nbpass], uz[nbpass], xprime_ref, zprime_ref, php_ref;
	struct Particle test_part;
	FILE *wfile;
	if(outfilename != NULL) wfile = fopen(outfilename, "w");
	
	php_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	zprime_ref  = atan_ratio(reference->uz, php_ref);
	x_vect[0] = 1.e-5;
	x_vect[1] = 1.e-5;
	x_vect[2] = 1.e-5;
	x_vect[3] = 1.e-5;
	if(get_matrix_firstorder_4d(m4d, &band_bias, reference, x_vect[0], x_vect[1], x_vect[2], x_vect[3], latt, part_cross_latt, YES) != TRUE) return FALSE;
	decouple_matrix_parzen(m4d, decoup_m, decoup_to_coup, output, NULL);
	matrix_inverse_4d(decoup_to_coup, decoup_to_coup_inv);
	
	//au = aumin;
	i = 0;
	*aumax = 0.;
	for(;;) {
		au = aumin+i*austep;
		//CLRSCR();
		//printf("amp_u = %lf\t", au);
		//fflush(stdout);
		u_vect[0] = au; // u
		u_vect[1] = 0; //uprime
		u_vect[2] = av; //v
		u_vect[3] = 0; //vprime
		mvprod4(x_vect, decoup_to_coup, u_vect);
		test_part = *reference;
		test_part.s = 0;
		test_part.x += x_vect[0]/cos(xprime_ref);
		test_part.z += x_vect[2]/cos(zprime_ref);
		test_part.ux = sin(xprime_ref + x_vect[1])*cos(zprime_ref + x_vect[3]);
		test_part.uy = cos(xprime_ref + x_vect[1])*cos(zprime_ref + x_vect[3]);
		test_part.uz = sin(zprime_ref + x_vect[3]);
		test_part.hat = -2;
		
		if(track_n_turns_amp(xpart, zpart, ux, uy, uz, &test_part, latt, nbpass, 0, 0) != ALIVE) {
			printf("NO!\n");
			break;
		}
		else {
			printf("OK!\n");
			if(outfilename != NULL) {
				for(n=0;n<nbpass;n++) {
					comp_phase_space_coord(&(x_vect[0]), &(x_vect[1]), &(x_vect[2]), &(x_vect[3]), xpart[n], zpart[n], ux[n], uy[n], uz[n], reference->x, xprime_ref, reference->z, zprime_ref);
					mvprod4(u_vect, decoup_to_coup_inv, x_vect);
					fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le\n", x_vect[0], x_vect[1], x_vect[2], x_vect[3], u_vect[0], u_vect[1], u_vect[2], u_vect[3], au);
				}
			}
		}
		i++;
	}
	if(i==0) *aumax = 0.;
	else *aumax = au-austep;
	if(outfilename != NULL) fclose(wfile);
	printf("final: amp_u = %le\n",*aumax);
	return TRUE;
}

extern int acceptancev_auto(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double avmin, double *avmax, double avstep, double au) 
{
	int i, n;
	double m4d[4][4], decoup_m[4][4], decoup_to_coup[4][4], decoup_to_coup_inv[4][4], u_vect[4], x_vect[4], output[6];
	double av, band_bias, xpart[nbpass], zpart[nbpass], ux[nbpass], uy[nbpass], uz[nbpass], xprime_ref, zprime_ref, php_ref;
	struct Particle test_part;
	FILE *wfile;
	if(outfilename != NULL) wfile = fopen(outfilename, "w");
	
	php_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	zprime_ref  = atan_ratio(reference->uz, php_ref);
	x_vect[0] = 1.e-5;
	x_vect[1] = 1.e-5;
	x_vect[2] = 1.e-5;
	x_vect[3] = 1.e-5;
	if(get_matrix_firstorder_4d(m4d, &band_bias, reference, x_vect[0], x_vect[1], x_vect[2], x_vect[3], latt, part_cross_latt, YES) != TRUE) return FALSE;
	decouple_matrix_parzen(m4d, decoup_m, decoup_to_coup, output, NULL);
	matrix_inverse_4d(decoup_to_coup, decoup_to_coup_inv);
	
	//av = avmin;
	i = 0;
	*avmax = 0.;
	for(;;) {
		av = avmin+i*avstep;
		CLRSCR();
		printf("amp_v = %lf\t", av);
		fflush(stdout);
		u_vect[0] = au; // u
		u_vect[1] = 0; //uprime
		u_vect[2] = av; //v
		u_vect[3] = 0; //vprime
		mvprod4(x_vect, decoup_to_coup, u_vect);
		test_part = *reference;
		test_part.s = 0;
		test_part.x += x_vect[0]/cos(xprime_ref);
		test_part.z += x_vect[2]/cos(zprime_ref);
		test_part.ux = sin(xprime_ref + x_vect[1])*cos(zprime_ref + x_vect[3]);
		test_part.uy = cos(xprime_ref + x_vect[1])*cos(zprime_ref + x_vect[3]);
		test_part.uz = sin(zprime_ref + x_vect[3]);
		test_part.hat = -2;
		
		if(track_n_turns_amp(xpart, zpart, ux, uy, uz, &test_part, latt, nbpass, 0, 0) != ALIVE) {
			printf("NO!\n");
			break;
		}
		else {
			printf("OK!\n");
			if(outfilename != NULL) {
				for(n=0;n<nbpass;n++) {
					comp_phase_space_coord(&(x_vect[0]), &(x_vect[1]), &(x_vect[2]), &(x_vect[3]), xpart[n], zpart[n], ux[n], uy[n], uz[n], reference->x, xprime_ref, reference->z, zprime_ref);
					mvprod4(u_vect, decoup_to_coup_inv, x_vect);
					fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le\n", x_vect[0], x_vect[1], x_vect[2], x_vect[3], u_vect[0], u_vect[1], u_vect[2], u_vect[3], au);
				}
			}
		}
		i++;
	}
	if(i==0) *avmax = 0.;
	else *avmax = av-avstep;
	if(outfilename != NULL) fclose(wfile);
	printf("final: amp_u = %le\n",*avmax);
	return TRUE;
}

extern void vffa_sect_to_rect(char *filename, int n, double theta0_deg, double r0, double b0, double z0, double m, double rho, double theta_en_deg, double theta_ex_deg, double lambda, int interpol_order)
{
	int i;
	double theta1, theta1_deg, rm, theta2, theta2_deg, tilt_angle_deg, ltot, l_1mag, f_pac=0.8;
	
	ltot = fabs(rho*(theta_ex_deg-theta_en_deg))*PI/180.;
	l_1mag = (ltot*f_pac)/n;
	
	for(i=0;i<n;i++) {
		theta1_deg = theta_en_deg+(i+0.5)/n*(theta_ex_deg-theta_en_deg);
		theta1 = theta1_deg*PI/180.;
		rm = sqrt(pow(rho*sin(theta1), 2) + pow(r0+rho*(cos(theta1) - 1.), 2));
		theta2 = asin(rho/rm*sin(theta1));
		theta2_deg = theta2*180./PI;
		tilt_angle_deg = theta1_deg-theta2_deg;
		if(filename!=NULL) {
			FILE *wfile;
			wfile=fopen(filename, "a");
			fprintf(wfile, "%lf	%lf	%lf	%lf	%lf	%lf\n", theta0_deg+theta2_deg, rm, b0, z0, m, tilt_angle_deg);
			fprintf(wfile, "\t%lf	%lf	%i	1\n", -l_1mag/2, lambda, interpol_order);
			fprintf(wfile, "\t%lf	%lf	%i	1\n\n", l_1mag/2, lambda, interpol_order);
			fclose(wfile);
		}
		else {
			printf("%lf	%lf	%lf	%lf	%lf	%lf\n", theta0_deg+theta2_deg, rm, b0, z0, m, tilt_angle_deg);
			printf("\t%lf	%lf	%i	1\n", -l_1mag/2, lambda, interpol_order);
			printf("\t%lf	%lf	%i	1\n\n", l_1mag/2, lambda, interpol_order);
		}
	}
}

extern void set_vffa_rect_convergence_limit(struct Cell *cell, int doyoucompute, double conv_lim, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	
	if(doyoucompute==YES) compute_convergence_limit_vffa(cell, add_contribution_comp);
	else {
		for(i=0;i<cell->nbcomp;i++) cell->mpara[i][6] = conv_lim;
	}
	
}

extern void enge_mult(double b0_mult, double efbe, double efbs, double lambda_add, char *xrange, char *outfile)
{
	char func[500];
	double lambda_mult = lambda_add/2.;
	
	sprintf(func,"(%le)/((1.+exp((x-(%le))/(%le)))*(1.+exp(((%le)-x)/(%le)))) with lines lc 1", b0_mult, efbs, lambda_mult, efbe, lambda_mult);
	plotfunction(func, "x", "f(x)", xrange, NULL, outfile, NULL, NULL);
}

extern void enge_add(double b0_add, double efbe, double efbs, double lambda_add, char *xrange, char *outfile)
{
	char func[500];
	
	sprintf(func,"(%le)/2.*( tanh(((%le)-x)/(%le)) + tanh((x-(%le))/(%le)) ) with lines lc 1", b0_add, efbs, lambda_add, efbe, lambda_add);
	plotfunction(func, "x", "f(x)", xrange, NULL, outfile, NULL, NULL);
}

extern void enge_mult_add(double b0_mult, double b0_add, double efbe, double efbs, double lambda_add, char *xrange, char *outfile)
{
	char func[1000];
	double lambda_mult = lambda_add/2.;
	
	sprintf(func,"(%le)/((1+exp((x-(%le))/(%le)))*(1+exp(((%le)-x)/(%le)))) - ((%le)/2.*(tanh(((%le)-x)/(%le)) + tanh((x-(%le))/(%le)))) with lines lc 1", 
	b0_mult, efbs, lambda_mult, efbe, lambda_mult, b0_add, efbs, lambda_add, efbe, lambda_add);
	plotfunction(func, "x [m]", "B0_m f_m(x) - B0_a f_a(x) [T]", xrange, NULL, outfile, NULL, NULL);
}

extern void adjust_b0_sect_from_rect(double x, double y, double z, double bz_opal, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double precision, int doyouf)
{
	int i,j, n=500;
	double sign,bx,by,bz;
	
	if(doyouf==YES) {
		sign = 1.;
		printf("adjust B0F\n");
	}
	else {
		sign = -1.;
		printf("adjust B0D\n");
	}
	
	for(i=0;i<n;i++) {
		if(get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp)!=ALIVE) errorstop("bfield not alive");
		printf("i=%i, bz = %le, bz_opal = %le, b01=%le, rho1=%le, b02=%le, rho2=%le\n", i, bz, bz_opal, cell->mpara[0][2], cell->mpara[0][5], cell->mpara[cell->nbcomp-1][2], cell->mpara[cell->nbcomp-1][5]);
		if(fabs(bz-bz_opal)<precision) break;
		for(j=0;j<cell->nbcomp;j++) {
			if(sign*cell->mpara[j][5] > 0) cell->mpara[j][2] -= sign*(bz-bz_opal)/2.;
		}
	}
	for(j=0;j<cell->nbcomp;j++) {
		if(sign*cell->mpara[j][5] > 0) {
			printf("B0 = %le\n", cell->mpara[j][2]);
			break;
		}
	}
}

extern void test_add_ff(double efbe, double efbs, double lambda_tan, char *xrange)
{
	char func[1000], func1[1000], func2[1000];
	double lambda_exp = lambda_tan/2.;
	sprintf(func1, "1./2.*(tanh(((%le)-x)/(%le)) + tanh((x-(%le))/(%le)))", efbs, lambda_tan, efbe, lambda_tan);
	sprintf(func2, "1/(1+exp(((%le)-x)/(%le))) + 1/(1+exp((x-(%le))/(%le))) -1", efbe, lambda_exp, efbs, lambda_exp);
	sprintf(func, "%s with lines lc 1", func1);
	plotfunction(func, "x [m]", "f_tan(x)", xrange, NULL, "output/func1_out.eps", NULL, NULL);
	sprintf(func, "%s with lines lc 1", func2);
	plotfunction(func, "x [m]", "f_exp(x)", xrange, NULL, "output/func2_out.eps", NULL, NULL);
	sprintf(func, "(%s)-(%s) with lines lc 1", func2,func1);
	plotfunction(func, "x [m]", "dif", xrange, NULL, "output/func_dif_out.eps", NULL, NULL);
}

extern void vffa_maxwell_test_heatmap(char *txt_prefix, char *txtaspect, double r0, double rmax, int r_step, double th0, double thmax, int th_step, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), int doyoucirc, double order)
{
	char name_div[300], name_curlx[300], name_curly[300],name_curlz[300], xrange[100], yrange[100];
	char eps_div[300], eps_curlx[300], eps_curly[300], eps_curlz[300];
	int i;
	double min,max;
	
	sprintf(name_div,  "%s_div.dat"  , txt_prefix);
	sprintf(name_curlx,"%s_curlx.dat", txt_prefix);
	sprintf(name_curly,"%s_curly.dat", txt_prefix);
	sprintf(name_curlz,"%s_curlz.dat", txt_prefix);
	sprintf(eps_div,  "%s_div.eps"  , txt_prefix);
	sprintf(eps_curlx,"%s_curlx.eps", txt_prefix);
	sprintf(eps_curly,"%s_curly.eps", txt_prefix);
	sprintf(eps_curlz,"%s_curlz.eps", txt_prefix);
	emptyfile(name_div);
	emptyfile(name_curlx);
	emptyfile(name_curly);
	emptyfile(name_curlz);
	
	for(i=0;i<cell->nbcomp;i++) {
		cell->efben[i][2] = order;
		cell->efbex[i][2] = order;
	}
	maxwell_test_map_radius(txt_prefix, r0, rmax, r_step, th0, thmax, th_step, z, dx, dy, dz, cell, add_contribution_comp, doyoucirc);
	find_boundaries_heated_map_file2(name_div, xrange, yrange, &min, &max);
	easyplot3dmap_plusfile(txtaspect, "2", "1", "lines lt 1 lc 0 lw 2", name_div, "x [m]", "y [m]", "div B [T/m]", NULL, NULL, eps_div, NULL, min, max);
	find_boundaries_heated_map_file2(name_curlx, xrange, yrange, &min, &max);
	easyplot3dmap_plusfile(txtaspect, "2", "1", "lines lt 1 lc 0 lw 2", name_curlx, "x [m]", "y [m]", "(curl B)_x [T/m]", NULL, NULL, eps_curlx, NULL, min, max);
	find_boundaries_heated_map_file2(name_curly, xrange, yrange, &min, &max);
	easyplot3dmap_plusfile(txtaspect, "2", "1", "lines lt 1 lc 0 lw 2", name_curly, "x [m]", "y [m]", "(curl B)_y [T/m]", NULL, NULL, eps_curly, NULL, min, max);
	find_boundaries_heated_map_file2(name_curlz, xrange, yrange, &min, &max);
	easyplot3dmap_plusfile(txtaspect, "2", "1", "lines lt 1 lc 0 lw 2", name_curlz, "x [m]", "y [m]", "(curl B)_z [T/m]", NULL, NULL, eps_curlz, NULL, min, max);
}

extern void vffa_compare_scode(char *fileout, char *filescode, double xmin, double xmax, double xstep, double ymin, double ymax, double ystep, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	char namex[300],namey[300],namez[300];
	int nblines, nbx, nby, i,j;
	double x,y,z,bx_scode,by_scode,bz_scode,bx,by,bz,difbx,difby,difbz;
	FILE *rfile=NULL;
	FILE *wfilex,*wfiley,*wfilez;
	
	nblines = get_nb_lines_file(filescode);
	if(fabs(xstep)>TINYLENGTH) nbx = (int) ((xmax - xmin)/xstep) +1;
	else nbx = 1;
	if(fabs(ystep)>TINYLENGTH) nby = (int) ((ymax - ymin)/ystep) +1;
	else nby = 1;
	if(nblines-nbx*nby !=0) errorstop("strange nb of lines");
	
	rfile = fopen(filescode, "r");
	if(rfile==NULL) errorstop("cannot open scode file\n");
	
	sprintf(namex,"%s_x.dat", fileout);
	wfilex = fopen(namex,"w");
	sprintf(namey,"%s_y.dat", fileout);
	wfiley = fopen(namey,"w");
	sprintf(namez,"%s_z.dat", fileout);
	wfilez = fopen(namez,"w");
	
	for(j=0;j<nby;j++) {
		for(i=0;i<nbx;i++) {
			fscanf(rfile, "%le	%le	%le	%le	%le	%le",&x,&z,&y,&bx_scode,&bz_scode,&by_scode);
			//printf("scode: x=%le,y=%le\n",x,y);
			x += cell->mpara[0][1]*cos(cell->mpara[0][0]);
			y += cell->mpara[0][1]*sin(cell->mpara[0][0]);
			//printf("x=%le,y=%le\n",x,y);
			if(get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp)!=ALIVE) errorstop("bfield not alive");
			difbx = (bx_scode + bx);
			difby = (by_scode + by);
			difbz = (bz_scode + bz);
			fprintf(wfilex, "%le	%le	%le\n", y,x,difbx);
			fprintf(wfiley, "%le	%le	%le\n", y,x,difby);
			fprintf(wfilez, "%le	%le	%le\n", y,x,difbz);
		}
		fprintf(wfilex, "\n");
		fprintf(wfiley, "\n");
		fprintf(wfilez, "\n");
	}
	fclose(rfile);
	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
}

extern void compute_fieldheatmap_scode(char *fileout, char *filescode, double xmin, double xmax, double xstep, double ymin, double ymax, double ystep, struct Cell *cell)
{
	char namex[300],namey[300],namez[300];
	int nblines, nbx, nby, i,j;
	double x,y,z,bx,by,bz;
	FILE *rfile=NULL;
	FILE *wfilex,*wfiley,*wfilez;
	
	nblines = get_nb_lines_file(filescode);
	if(fabs(xstep)>TINYLENGTH) nbx = (int) ((xmax - xmin)/xstep) +1;
	else nbx = 1;
	if(fabs(ystep)>TINYLENGTH) nby = (int) ((ymax - ymin)/ystep) +1;
	else nby = 1;
	if(nblines-nbx*nby !=0) errorstop("strange nb of lines");
	
	rfile = fopen(filescode, "r");
	if(rfile==NULL) errorstop("cannot open scode file\n");
	
	sprintf(namex,"%s_bx.dat", fileout);
	wfilex = fopen(namex,"w");
	sprintf(namey,"%s_by.dat", fileout);
	wfiley = fopen(namey,"w");
	sprintf(namez,"%s_bz.dat", fileout);
	wfilez = fopen(namez,"w");
	
	for(j=0;j<nby;j++) {
		for(i=0;i<nbx;i++) {
			fscanf(rfile, "%le	%le	%le	%le	%le	%le",&x,&z,&y,&bx,&bz,&by);
			//printf("scode: x=%le,y=%le\n",x,y);
			x += cell->mpara[0][1]*cos(cell->mpara[0][0]);
			y += cell->mpara[0][1]*sin(cell->mpara[0][0]);
			//printf("x=%le,y=%le\n",x,y);
			fprintf(wfilex, "%le	%le	%le\n", y,x,bx);
			fprintf(wfiley, "%le	%le	%le\n", y,x,by);
			fprintf(wfilez, "%le	%le	%le\n", y,x,bz);
		}
		fprintf(wfilex, "\n");
		fprintf(wfiley, "\n");
		fprintf(wfilez, "\n");
	}
	
	fclose(rfile);
	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
}

extern void wrapper_plot_heatmap_2dfile(char *file, char *file_2d, char *xcolumn, char *ycolumn, char *with, char *xlabel, char *ylabel, char *cblabel, char *psfilename, char *setoption)
{
	char xrange[100],yrange[100];
	double min, max, xmin, xmax, ymin, ymax, xmin_2d, xmax_2d, ymin_2d, ymax_2d;
	find_boundaries_heated_map_file(file, &xmin, &xmax, &ymin, &ymax, &min, &max);
	if(file_2d!=NULL) {
		find_boundaries_2d_file(file_2d, &xmin_2d, &xmax_2d, &ymin_2d, &ymax_2d);
		sprintf(xrange,"[%lf:%lf]", MIN(xmin,ymin_2d), MAX(xmax, ymax_2d));
		sprintf(yrange,"[%lf:%lf]", MIN(ymin,xmin_2d), MAX(ymax, xmax_2d));
		easyplot3dmap_plusfile(file_2d, xcolumn, ycolumn, with, file, xlabel, ylabel, cblabel, xrange, yrange, psfilename, setoption, min, max);
	}
	else {
		sprintf(xrange,"[%lf:%lf]", xmin, xmax);
		sprintf(yrange,"[%lf:%lf]", ymin, ymax);
		easyplot3dmap(file, xlabel, ylabel, cblabel, xrange, yrange, psfilename, setoption, min, max);
	}
}

extern void wrapper_plot_fieldheatmap(char *fileout, char *aspectfile)
{
	char namex[300],namey[300],namez[300], outx[300],outy[300],outz[300];
	sprintf(namex,"%s_bx.dat", fileout);
	sprintf(namey,"%s_by.dat", fileout);
	sprintf(namez,"%s_bz.dat", fileout);
	sprintf(outx,"%s_bx.eps", fileout);
	sprintf(outy,"%s_by.eps", fileout);
	sprintf(outz,"%s_bz.eps", fileout);
	
	wrapper_plot_heatmap_2dfile(namex, aspectfile, "2", "1", "lines lt 1 lc 0 lw 2", "y (long) [m]", "x (hor) [m]", "Bx [T]", outx, NULL);
	wrapper_plot_heatmap_2dfile(namey, aspectfile, "2", "1", "lines lt 1 lc 0 lw 2", "y (long) [m]", "x (hor) [m]", "By [T]", outy, NULL);
	wrapper_plot_heatmap_2dfile(namez, aspectfile, "2", "1", "lines lt 1 lc 0 lw 2", "y (long) [m]", "x (hor) [m]", "Bz [T]", outz, NULL);
}

extern void aspect_rect_vffa3d(char *filename, struct Cell *cell, int doyouerr)
{
	int j, doyouplotcomp=YES;
	double x, y,z, zmin,zmax, x_0, y_0, ffbe, ffbs, lambda_e, angle, xwerr,ywerr,zwerr;//, lambda_s;
	FILE *wfile;
	
	
	wfile = fopen(filename, "a");
	for(j = 0; j < cell->nbcomp; j++) {
		if(cell->boun.thmax!=0) {
			if(cell->mpara[j][0]>0-TINYLENGTH && cell->mpara[j][0]<cell->boun.thmax+TINYLENGTH) doyouplotcomp = YES;
			else doyouplotcomp = NO;
		}
		else {
			if(cell->mpara[j][0] > 0-TINYLENGTH && cell->mpara[j][0] < cell->boun.ymax+TINYLENGTH) doyouplotcomp = YES;
			else doyouplotcomp = NO;
		}
		//printf("comp %i, doyouplot=%i, YES=%i\n", j, doyouplotcomp, YES);
		
		//if(cell->mpara[j][0]>0-TINYLENGTH && cell->mpara[j][0]<cell->boun.thmax+TINYLENGTH) doyouplotcomp = YES;
		//else doyouplotcomp = NO;
		if(doyouplotcomp==YES) {
			if(cell->boun.thmax!=0) {
				x_0 = cell->mpara[j][1]*cos(cell->mpara[j][0]);
				y_0 = cell->mpara[j][1]*sin(cell->mpara[j][0]);
			}
			else {
				x_0 = cell->mpara[j][1];
				y_0 = cell->mpara[j][0];
			}
			zmin = cell->collim.zmin;
			zmax = cell->collim.zmax;

			ffbe = cell->efben[j][0];
			ffbs = cell->efbex[j][0];
			//lambda_e = cell->efben[j][1];
			//lambda_s = cell->efbex[j][1];
			lambda_e = cell->mpara[j][6];
			//angle = cell->mpara[j][0]+cell->mpara[j][5];
			angle = cell->mpara[j][5];
			if(cell->boun.thmax!=0) angle +=cell->mpara[j][0];
			
			x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle); //D
			y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmin;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
		
			x = x_0 - ffbe*sin(angle) + lambda_e*cos(angle); //C
			y = y_0 + ffbe*cos(angle) + lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmin;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbe*sin(angle); //B
			y = y_0 + ffbe*cos(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmin;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbs*sin(angle); //E
			y = y_0 + ffbs*cos(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmin;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbe*sin(angle); //B
			y = y_0 + ffbe*cos(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmin;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbe*sin(angle) - lambda_e*cos(angle); //A
			y = y_0 + ffbe*cos(angle) - lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmin;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
		
			x = x_0 - ffbs*sin(angle) - lambda_e*cos(angle); //F
			y = y_0 + ffbs*cos(angle) - lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmin;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
		
			x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle); //D
			y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmin;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			
			
			x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle); //J
			y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmax;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
		
			x = x_0 - ffbe*sin(angle) + lambda_e*cos(angle); //I
			y = y_0 + ffbe*cos(angle) + lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmax;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbe*sin(angle) + lambda_e*cos(angle); //C
			y = y_0 + ffbe*cos(angle) + lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmin;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbe*sin(angle) + lambda_e*cos(angle); //I
			y = y_0 + ffbe*cos(angle) + lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmax;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbe*sin(angle); //H
			y = y_0 + ffbe*cos(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmax;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbe*sin(angle); //B
			y = y_0 + ffbe*cos(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmin;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbe*sin(angle); //H
			y = y_0 + ffbe*cos(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmax;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbs*sin(angle); //K
			y = y_0 + ffbs*cos(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmax;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbs*sin(angle); //E
			y = y_0 + ffbs*cos(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmin;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbs*sin(angle); //K
			y = y_0 + ffbs*cos(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmax;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbe*sin(angle); //H
			y = y_0 + ffbe*cos(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmax;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbe*sin(angle) - lambda_e*cos(angle); //G
			y = y_0 + ffbe*cos(angle) - lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmax;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbe*sin(angle) - lambda_e*cos(angle); //A
			y = y_0 + ffbe*cos(angle) - lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmin;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbe*sin(angle) - lambda_e*cos(angle); //G
			y = y_0 + ffbe*cos(angle) - lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmax;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbs*sin(angle) - lambda_e*cos(angle); //L
			y = y_0 + ffbs*cos(angle) - lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmax;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbs*sin(angle) - lambda_e*cos(angle); //F
			y = y_0 + ffbs*cos(angle) - lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmin;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbs*sin(angle) - lambda_e*cos(angle); //L
			y = y_0 + ffbs*cos(angle) - lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmax;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle); //J
			y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
			if(doyouerr==YES) {
				apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
				x = xwerr;
				y = ywerr;
				z = zwerr;
			}
			else z = zmax;
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le  %le\n",x,y,z);
			fprintf(wfile, "\n\n");
		}
	}
	fclose(wfile);
}
//*/

extern void write_alierror_para_triplet_vffa(char *filename, struct Lattice *err_latt)
{
	int i,j;
	FILE *wfile;
	
	wfile = fopen(filename, "w");
	for(i=0;i<err_latt->nbcell;i++) {
		for(j=1;j<4;j++) {
			fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le\n", 
			err_latt->cell[i].alierror[j][1],err_latt->cell[i].alierror[j][2],err_latt->cell[i].alierror[j][3],
			err_latt->cell[i].alierror[j][4],err_latt->cell[i].alierror[j][5],err_latt->cell[i].alierror[j][6],
			err_latt->cell[i].alierror[j][7]);
		}
		
	}
	
	fclose(wfile);
}

extern int error_study_triplet_vffa(char *prefix, struct Lattice *latt, struct Particle *ref_part, long seed, double rms_shift_error, double rms_twist_error)
{
	char err_para_name[300], dau_name[300], dav_name[300], res_name[300];
	int doyouerase = NO;
	double qx, qz, beta_x, alpha_x, beta_z, alpha_z, aumax, avmax;
	double prec_clo=1.e-11;
	double da_offset = 0.001;
	struct Lattice latt_err;
	struct Particle test_part;
	FILE *wfile;
	
	sprintf(err_para_name,"%s_shift%lf_twist%lf_seed%ld_errpara.dat", prefix, rms_shift_error, rms_twist_error, seed);
	sprintf(dau_name,"%s_shift%lf_twist%lf_seed%ld_poincarreu.dat", prefix, rms_shift_error, rms_twist_error, seed);
	sprintf(dav_name,"%s_shift%lf_twist%lf_seed%ld_poincarrev.dat", prefix, rms_shift_error, rms_twist_error, seed);
	sprintf(res_name,"%s_shift%lf_twist%lf_results.dat", prefix, rms_shift_error, rms_twist_error);
	if(doyouerase==YES) {
		emptyfile(err_para_name);
		emptyfile(dau_name);
		emptyfile(dav_name);
		emptyfile(res_name);
	}
	
	test_part = *ref_part;
	gene_alignerror_latt(latt, &latt_err);
	comp_error_latt_vffa(&latt_err, seed, rms_shift_error, rms_twist_error);
	print_latt_para(NULL, "latt_alierror", &latt_err);
	write_alierror_para_triplet_vffa(err_para_name, &latt_err);
	if(find_closed_orbite_xxp_zzp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), &(test_part.z), &(test_part.uz), prec_clo, &latt_err, YES) != TRUE) {
		free_latt(&latt_err);
		return FALSE;
	}
	calc_tune_twiss(&test_part, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-5, 1.e-5, 1.e-5, 1.e-5, &latt_err, part_cross_latt, YES, NULL,1);//flag: 0=2d, 1=4d decoupled, anything else=4d undecoupled.
	acceptanceu_auto(dau_name, &test_part, &latt_err, 1000, 0.001, &aumax, 0.002, da_offset); 
	acceptancev_auto(dav_name, &test_part, &latt_err, 1000, 0.001, &avmax, 0.002, da_offset); 
	wfile = fopen(res_name, "a");
	fprintf(wfile, "%ld	%le	%le	%le	%le", seed, test_part.x, test_part.ux, test_part.z, test_part.uz);
	fprintf(wfile, "	%le	%le	%le	%le	%le	%le	%le	%le	%le\n", qx,qz,beta_x, alpha_x, beta_z, alpha_z,aumax,avmax,da_offset);
	fclose(wfile);
	free_latt(&latt_err);
	return TRUE;
}

//chromaticity to addå
extern int m_error_study_triplet_vffa(char *prefix, struct Lattice *latt, struct Particle *ref_part, long seed, double rms_error)
{
	char err_para_name[300], dau_name[300], dav_name[300], res_name[300];
	int i,j, doyouerase = NO;
	double qx, qz, beta_x, alpha_x, beta_z, alpha_z, aumax, avmax;
	double prec_clo=1.e-11;
	double da_offset = 0.001;
	struct Lattice latt_err;
	struct Particle test_part;
	FILE *wfile;
	FILE *wfile_para;
	
	sprintf(err_para_name,"%s_merror%lf_seed%ld_errpara.dat", prefix, rms_error, seed);
	sprintf(dau_name,     "%s_merror%lf_seed%ld_poincarreu.dat", prefix, rms_error, seed);
	sprintf(dav_name,     "%s_merror%lf_seed%ld_poincarrev.dat", prefix, rms_error, seed);
	sprintf(res_name,     "%s_merror%lf_results.dat", prefix, rms_error);
	if(doyouerase==YES) {
		emptyfile(err_para_name);
		emptyfile(dau_name);
		emptyfile(dav_name);
		emptyfile(res_name);
	}
	
	test_part = *ref_part;
	copy_latt_1period(latt, &latt_err);
	comp_m_error_latt_vffa(&latt_err, seed, rms_error);
	print_latt_para(NULL, "latt_merror", &latt_err);
	wfile_para = fopen(err_para_name, "w");
	for(i=0;i<latt_err.nbcell;i++) {
		for(j=1;j<4;j++) fprintf(wfile_para, "%i	%i	%le\n", i,j, latt_err.cell[i].mpara[j][4]);
	}
	fclose(wfile_para);
	if(find_closed_orbite_xxp_zzp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), &(test_part.z), &(test_part.uz), prec_clo, &latt_err, YES) != TRUE) {
		free_latt(&latt_err);
		return FALSE;
	}
	calc_tune_twiss(&test_part, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-5, 1.e-5, 1.e-5, 1.e-5, &latt_err, part_cross_latt, YES, NULL,1);//flag: 0=2d, 1=4d decoupled, anything else=4d undecoupled.
	acceptanceu_auto(dau_name, &test_part, &latt_err, 1000, 0.001, &aumax, 0.002, da_offset); 
	acceptancev_auto(dav_name, &test_part, &latt_err, 1000, 0.001, &avmax, 0.002, da_offset); 
	wfile = fopen(res_name, "a");
	fprintf(wfile, "%ld	%le	%le	%le	%le", seed, test_part.x, test_part.ux, test_part.z, test_part.uz);
	fprintf(wfile, "	%le	%le	%le	%le	%le	%le	%le	%le	%le\n", qx,qz,beta_x, alpha_x, beta_z, alpha_z,aumax,avmax,da_offset);
	fclose(wfile);
	free_latt(&latt_err);
	return TRUE;
}

extern void compare_err_vffa_triplet(char *noerr_file, char *err_file, char *dif_file, double trans_shift, double rot_shift)
{
	long seed;
	int n, nblines;
	
	double x_noerr, ux_noerr, z_noerr, uz_noerr, qx_noerr, qz_noerr, beta_x_noerr, alpha_x_noerr, beta_z_noerr, alpha_z_noerr, aumax_noerr, avmax_noerr, da_offset;
	double x, ux, z, uz, qx, qz, beta_x, alpha_x, beta_z, alpha_z, aumax, avmax, cod;
	FILE *rfile_noerr=NULL;
	FILE *rfile_err=NULL;
	FILE *wfile;
	
	nblines = get_nb_lines_file(err_file);
	
	rfile_noerr = fopen(noerr_file, "r");
	if(rfile_noerr==NULL) errorstop("cannot open no error file");
	rfile_err = fopen(err_file, "r");
	if(rfile_err==NULL) errorstop("cannot open error file");
	wfile = fopen(dif_file, "w");
	fscanf(rfile_noerr, "%ld	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le", &seed, &x_noerr, &ux_noerr, &z_noerr, &uz_noerr, &qx_noerr, &qz_noerr, &beta_x_noerr, &alpha_x_noerr, &beta_z_noerr, &alpha_z_noerr, &aumax_noerr, &avmax_noerr, &da_offset);
	fclose(rfile_noerr);
	
	for(n=0;n<nblines;n++) {
		fscanf(rfile_err, "%ld	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le", &seed, &x, &ux, &z, &uz, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, &aumax, &avmax, &da_offset);
		cod = sqrt(pow(x-x_noerr,2)+pow(z-z_noerr,2));
		fprintf(wfile, "%ld	%le	%le	%le	%le	%le	%le	%le\n", seed, trans_shift, rot_shift, cod, qx-qx_noerr, qz-qz_noerr, aumax-aumax_noerr, avmax-avmax_noerr);
	}
	
	fclose(wfile);
	fclose(rfile_err);
}

extern void gene_matched_beam_emi_4d(char *txtfile, double emit_u, double twiss_betau, double twiss_alphau, double emit_v, double twiss_betav, double twiss_alphav, double r_inv[4][4], int nparts, long seed, int doyouhalo)
{
	long idum;
	int i;
	double uv_vect[4], xz_vect[4];
	double x_ctr, xprime_ctr, z_ctr, zprime_ctr;
	double factor_rms_tot;
	double un, upn, vn, vpn,x,xp,z,zp;
	FILE *wfile;
	FILE *wfile2;
	factor_rms_tot = 6.;
	x_ctr = 0;
	xprime_ctr = 0;
	z_ctr = 0;
	zprime_ctr = 0;
	
	printf("gene_matched_beam_emi_4d: nparts = %i\n", nparts);
	//test errors and warnings
	if(nparts <= 0) errorstop("!!!ERROR in gene_matched_beam_emi_4d:\n nparts <= 0");
	if(twiss_betau <= 0) errorstop("!!!ERROR in gene_matched_beam_emi_4d: twiss_betau <= 0");
	if(twiss_betav <= 0) errorstop("!!!ERROR in gene_matched_beam_emi_4d: twiss_betav <= 0");
	
	//initialise idum (for the random generator)
	if(seed >=0) idum = -seed;
	else idum = seed;
	
	wfile = fopen(txtfile,"w");
	if(doyouhalo==YES) wfile2 = fopen("data/emi/halo_uvspace.dat","w");
	else wfile2 = fopen("data/emi/core_uvspace.dat","w");
	
	//generate a 6D beam, inside ellipses centered on central particle
	for(i = 1; i < nparts+1 ; i++) {
		printf("part number %i\n", i);
		if(doyouhalo==YES) {
			for(;;) {
				un = gasdev(&idum);
				upn = gasdev(&idum);
				//if(fabs(un) < factor_rms_tot*sigmau && fabs(upn) < factor_rms_tot*sigmaup) break;
				if(un*un/16. + upn*upn/16. > 1. && un*un/36. + upn*upn/36. < 1.) break;
				//if(un*un/16. + upn*upn/16. < 1.) break;
			}
			for(;;) {
				vn = gasdev(&idum);
				vpn = gasdev(&idum);
				//if(fabs(vn) < factor_rms_tot*sigmav && fabs(vpn) < factor_rms_tot*sigmavp) break;
				if(vn*vn/16. + vpn*vpn/16. > 1. && vn*vn/36. + vpn*vpn/36. < 1.) break;
				//if(vn*vn/16. + vpn*vpn/16. < 1.) break;
			}
		}
		else {
			for(;;) {
				un = gasdev(&idum);
				upn = gasdev(&idum);
				//if(fabs(un) < factor_rms_tot*sigmau && fabs(upn) < factor_rms_tot*sigmaup) break;
				//if(un*un/16. + upn*upn/16. > 1. && un*un/36. + upn*upn/36. < 1.) break;
				if(un*un/16. + upn*upn/16. < 1.) break;
			}
			for(;;) {
				vn = gasdev(&idum);
				vpn = gasdev(&idum);
				//if(fabs(vn) < factor_rms_tot*sigmav && fabs(vpn) < factor_rms_tot*sigmavp) break;
				//if(vn*vn/16. + vpn*vpn/16. > 1. && vn*vn/36. + vpn*vpn/36. < 1.) break;
				if(vn*vn/16. + vpn*vpn/16. < 1.) break;
			}
		}
		
		uv_vect[0] = sqrt(emit_u*twiss_betau)*un;
		uv_vect[1] = sqrt(emit_u/twiss_betau)*(upn - twiss_alphau*un);
		uv_vect[2] = sqrt(emit_v*twiss_betav)*vn;
		uv_vect[3] = sqrt(emit_v/twiss_betav)*(vpn - twiss_alphav*vn);
		mvprod4(xz_vect, r_inv, uv_vect);
		x = x_ctr + xz_vect[0];
		xp = xprime_ctr + xz_vect[1];
		z = z_ctr + xz_vect[2];
		zp = zprime_ctr + xz_vect[3];
		fprintf(wfile2, "%le	%le	%le	%le\n", uv_vect[0],uv_vect[1],uv_vect[2],uv_vect[3]);
		fprintf(wfile, "%le	%le	%le	%le\n", x,xp,z,zp);
	}// and it is a loop, it does it again nparts times...
	fclose(wfile);
	fclose(wfile2);
}

extern void test_matched_beam_emi(double m4d[4][4], char *textfile, char *textout)
{
	int i,j,nblines;
	double xz_vect[4],x,xp,z,zp;
	FILE *rfile=NULL;
	FILE *wfile;
	
	nblines = get_nb_lines_file(textfile);
	rfile = fopen(textfile, "r");
	if(rfile==NULL) errorstop("cannot open file\n");
	wfile = fopen(textout, "w");
	for(i=0;i<nblines;i++) {
		fscanf(rfile, "%le	%le	%le	%le", &x,&xp,&z,&zp);
		xz_vect[0] = x;
		xz_vect[1] = xp;
		xz_vect[2] = z;
		xz_vect[3] = zp;
		for(j=0;j<10;j++) {
			mvprod4(xz_vect, m4d, xz_vect);
		}
		fprintf(wfile, "%le	%le	%le	%le\n", xz_vect[0],xz_vect[1],xz_vect[2],xz_vect[3]);
	}
	fclose(rfile);
	fclose(wfile);
}

extern void file_mean(char *textout, char *textfile)
{
	int nblines,i,dum;
	double a2,a3,a4,a5,a6,a7,a8,a2t=0,a3t=0,a4t=0,a5t=0,a6t=0,a7t=0,a8t=0;
	FILE *rfile = NULL;
	FILE *wfile;
	nblines = get_nb_lines_file(textfile);
	rfile = fopen(textfile, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	wfile = fopen(textout, "w");

	for(i=0;i<nblines;i++) {
		fscanf(rfile, "%i	%le	%le	%le	%le	%le	%le	%le", &dum,&a2,&a3,&a4,&a5,&a6,&a7,&a8);
		a2t+=a2;
		a3t+=a3;
		a4t+=a4;
		a5t+=a5;
		a6t+=a6;
		a7t+=a7;
		a8t+=a8;
		
	}
	a2t /=nblines;
	a3t /=nblines;
	a4t /=nblines;
	a5t /=nblines;
	a6t /=nblines;
	a7t /=nblines;
	a8t /=nblines;
	
	fprintf(wfile, "0	%le	%le	%le	%le	%le	%le	%le\n", a2t,a3t,a4t,a5t,a6t,a7t,a8t);
	fclose(rfile);
	fclose(wfile);
}

extern void file_histo(char *textout, char *textfile, int columny, double xmin, double xmax, int nbstep)
{
	int nblines,i,j, count[nbstep];
	double a,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,xstep;
	FILE *rfile = NULL;
	FILE *wfile;
	nblines = get_nb_lines_file(textfile);
	rfile = fopen(textfile, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	wfile = fopen(textout, "w");
	for(i=0;i<nbstep;i++) count[i]=0;
	xstep = (xmax-xmin)/(nbstep+1.0);
	printf("xstep=%le\n", xstep);
	for(i=0;i<nblines;i++) {
		fscanf(rfile, "%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le", &a1,&a2,&a3,&a4,&a5,&a6,&a7,&a8,&a9,&a10,&a11,&a12,&a13,&a14);
		if(columny==1) a=a1;
		else if(columny==2) a=a2;
		else if(columny==3) a=a3;
		else if(columny==4) a=a4;
		else if(columny==5) a=a5;
		else if(columny==6) a=a6;
		else if(columny==7) a=a7;
		else if(columny==8) a=a8;
		else if(columny==9) a=a9;
		else if(columny==10) a=a10;
		else if(columny==11) a=a11;
		else if(columny==12) a=a12;
		else if(columny==13) a=a13;
		else if(columny==14) a=a14;
		
		for(j=0;j<nbstep;j++) {
			if(a>xmin+j*xstep && a<xmin+xstep*(j+1)) {
				count[j]++;
				break;
			}
		}
	}
	for(i=0;i<nbstep;i++) fprintf(wfile, "%le	%i\n", xmin+j*xstep, count[i]);
	fclose(rfile);
	fclose(wfile);
}

extern void m_loc_int_by_cell_vffa(char *outfilename, double x, double ymin, double ymax, int nstep_y, double zmin, double zmax, int nstep_z, int nzref, double mref, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double zstep, m[nstep_z+1], z[nstep_z+1], intbdl[nstep_z+1], intbdlabs[nstep_z+1], errorm[nstep_z+1], errorbdl[nstep_z+1], intbdldes[nstep_z+1], btemp1, btemp2;
	FILE *wfile;
	
	zstep = (zmax - zmin)/nstep_z;
	printf("step=%lf\n",zstep);
	for(i = 0; i < nstep_z + 1; i++) {
		z[i] = zmin + i*zstep;
		//printf("%i, z=%lf\n",i,z[i]);
		if(ymax-ymin>TINYLENGTH) int_by_liney_cell(&(intbdl[i]), &(intbdlabs[i]), x, ymin, ymax, nstep_y, z[i], cell, add_contribution_comp);
		else if(get_bfield(x, ymin, z[i], &btemp1, &(intbdl[i]), &btemp2, cell, add_contribution_comp) != ALIVE) errorstop("m_loc_int_by_cell_vffa: get_bfield != ALIVE\n");
	}
	
	for(i = 0; i < nstep_z + 1; i++) {
		m[i] = log(intbdl[i]/intbdl[nzref])/(z[i] - z[nzref]);
		errorm[i] = (m[i] - mref)/mref;
	}
	if(nzref != nstep_z) m[nzref] = (m[nzref - 1] + m[nzref + 1])/2;
	else m[nzref] = m[nzref - 1];
	//plot the desired Bl
	for(i = 0; i < nstep_z + 1; i++) {
		if(i != nzref) {
			intbdldes[i] = intbdl[nzref]*exp(mref*(z[i] - z[nzref]));
			errorbdl[i] = (fabs(intbdl[i]) - fabs(intbdldes[i]))/fabs(intbdldes[i]);
		}
	}
	intbdldes[nzref] = intbdl[nzref];
	errorbdl[nzref] = 0;
	errorm[nzref] = (m[nzref] - mref)/mref;
	
	//data writing
	wfile = fopen(outfilename,"w");
	for(i = 0; i < nstep_z + 1; i++) {
		fprintf(wfile, "%lf  %lf  %lf  %le  %le  %le\n", z[i], intbdl[i], m[i], errorm[i], intbdldes[i], errorbdl[i]);
	}
	fclose(wfile);
}

//vertical FFA m=1/Bz*(dBz/dz)
extern void m_loc_int_bz_cell_vffa(char *outfilename, double x0, double ymin, double ymax, int nstep_y, double zmin, double zmax, int nstep_z, int nzref, double mref, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double zstep, m[nstep_z], z[nstep_z], intbdl[nstep_z], intbdlabs[nstep_z], errorm[nstep_z], errorbdl[nstep_z], intbdldes[nstep_z],btemp1, btemp2, avg_m=0, min_m=1.e10, max_m=-1.e10;
	FILE *wfile;
	
	//calcul de BL
	zstep = (zmax - zmin)/(nstep_z-1);
	for(i = 0; i < nstep_z; i++) {
		z[i] = zmin + i*zstep;
		//printf("%i\n",i);
		if(ymax-ymin>TINYLENGTH) int_bz_liney_cell2(&(intbdl[i]), &(intbdlabs[i]), x0, ymin, ymax, nstep_y, z[i], cell, add_contribution_comp);
		else if(get_bfield(x0, ymin, z[i], &btemp1, &btemp2, &(intbdl[i]), cell, add_contribution_comp) != ALIVE) errorstop("m_loc_int_bz_cell_vffa: get_bfield != ALIVE\n");
	}
	//calcul de m_loc
	for(i = 0; i < nstep_z; i++) {
		if(i != nzref) {
			m[i] = log(intbdl[i]/intbdl[nzref])/(z[i] - z[nzref]);
			errorm[i] = (m[i] - mref)/mref;
		}
	}
	if(nzref != nstep_z-1) m[nzref] = (m[nzref - 1] + m[nzref + 1])/2;
	else m[nzref] = m[nzref - 1];
	//plot the desired Bl
	for(i = 0; i < nstep_z; i++) {
		avg_m+=m[i];
		min_m = MIN(min_m, m[i]);
		max_m = MAX(max_m, m[i]);
		if(i != nzref) {
			intbdldes[i] = intbdl[nzref]*exp(mref*(z[i] - z[nzref]));
			errorbdl[i] = (fabs(intbdl[i]) - fabs(intbdldes[i]))/fabs(intbdldes[i]);
		}
	}
	avg_m/=(double) nstep_z;
	printf("avg vert m: %le, min: %le/100, max: %le/100\n", avg_m, 100*(min_m-avg_m)/avg_m, 100*(max_m-avg_m)/avg_m);
	
	intbdldes[nzref] = intbdl[nzref];
	errorbdl[nzref] = 0;
	errorm[nzref] = (m[nzref] - mref)/mref;
	//data writing
	wfile = fopen(outfilename,"w");
	for(i = 0; i < nstep_z; i++) {
		fprintf(wfile, "%lf  %lf  %lf  %le  %le  %le\n", z[i], intbdl[i], m[i], errorm[i], intbdldes[i], errorbdl[i]);
	}
	fclose(wfile);
}

extern void int_bz_liney_cell2(double *intbdl, double *intbdlabs, double x, double ymin, double ymax, int nstep, double z, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double ystep, y, bprobe[4];
	
	ystep = (ymax - ymin)/nstep; 
	*intbdl = 0;		//int{bz.dl}
	*intbdlabs = 0;		//int{|bz|.dl}
	
	for(i = 0; i < nstep + 1; i++) {
		y = ymin + i*ystep;
		//printf("y=%lf\n",y);
		if(get_bfield(x, y, z, &bprobe[1], &bprobe[2], &bprobe[3], cell, add_contribution_comp) != ALIVE) errorstop("int_bz_liney_cell: get_bfield != ALIVE\n");
		*intbdl += bprobe[3]*ystep;
		*intbdlabs += fabs(bprobe[3])*ystep;
	}
}

extern void int_by_liney_cell(double *intbdl, double *intbdlabs, double x, double ymin, double ymax, int nstep, double z, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double ystep, y, bprobe[4];
	
	ystep = (ymax - ymin)/nstep; 
	*intbdl = 0;		//int{bz.dl}
	*intbdlabs = 0;		//int{|bz|.dl}
	
	for(i = 0; i < nstep + 1; i++) {
		y = ymin + i*ystep;
		if(get_bfield(x, y, z, &bprobe[1], &bprobe[2], &bprobe[3], cell, add_contribution_comp) != ALIVE) errorstop("int_by_liney_cell: get_bfield != ALIVE\n");
		*intbdl += bprobe[2]*ystep;
		*intbdlabs += fabs(bprobe[2])*ystep;
	}
}

//!! change the output for error study21/10/20
extern void file_histo_from_txtfile(char *textout, char *textfile, int nbcolumn, int columny, double xmin, double xmax, int nbstep)
{
	int nblines,i,j;//, count[nbstep];
	double a,ascan[30],xstep;
	FILE *rfile = NULL;
	FILE *wfile = NULL;
	if(nbcolumn>30) errorstop("too many columns");
	nblines = get_nb_lines_file(textfile);
	rfile = fopen(textfile, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	wfile = fopen(textout, "w");
	if(wfile==NULL) errorstop("cannot open wfile");
	
	//for(i=0;i<nbstep;i++) count[i]=0;
	xstep = (xmax-xmin)/(nbstep+1.0);
	printf("xstep=%le\n", xstep);
	for(i=0;i<nblines;i++) {
		for(j=1;j<nbcolumn+1;j++) fscanf(rfile, "%le", &(ascan[j]));
		//for(j=1;j<nbcolumn+1;j++) printf("%le\t", ascan[j]);
		//printf("\n");
		a = ascan[columny];
		fprintf(wfile, "%le\n", a*1.e6*0.080031);
		//fprintf(wfile, "%le\n", a);
		//for(j=0;j<nbstep;j++) {
		//	if(a>xmin+j*xstep && a<xmin+xstep*(j+1)) {
		//		count[j]++;
		//		break;
		//	}
		//}
	}
	//for(i=0;i<nbstep;i++) fprintf(wfile, "%le	%i\n", xmin+i*xstep, count[i]);
	fclose(rfile);
	fclose(wfile);
}

extern void get_range_column(char *txtfile, int nbcolumn, int columnx, double *xmin, double *xmax)
{
	int i,j, nblines;
	double ascan[30];
	FILE *rfile=NULL;
	*xmin = +1.e15;
	*xmax = -1.e15;
	
	if(nbcolumn>30) errorstop("too many columns");
	nblines = get_nb_lines_file(txtfile);
	
	rfile = fopen(txtfile, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(i=0;i<nblines;i++) {
		for(j=1;j<nbcolumn+1;j++) fscanf(rfile, "%le", &(ascan[j]));
		//for(j=1;j<nbcolumn+1;j++) printf("%le\t", ascan[j]);
		//printf("\n");
		*xmin = MIN(ascan[columnx],*xmin);
		*xmax = MAX(ascan[columnx],*xmax);
	}
	fclose(rfile);
	printf("min, max = %le, %le\n", *xmin, *xmax);
}

extern void compute_cod_error_vffa_triplet(char *txtfile, double nsteps_cell, struct Lattice *latt, struct Particle *reference)
{
	int i1, i2;
	double step_th, step_y, x, z, bias, s0,s, th0, th, prec_clo=1.e-10;
	struct Particle part_noerr, part, test_part_noerr, test_part, test_part2_noerr, test_part2;
	struct Lattice tempo_latt, tempo_latt_noerr, latt_1cell, latt_1cell_noerr;
	FILE *wfile;
	
	copy_latt(&tempo_latt, latt);
	copy_error_latt(&tempo_latt, latt);
	copy_latt(&tempo_latt_noerr, latt);
	for(i1 = 0; i1 < tempo_latt_noerr.nbcell; i1++) tempo_latt_noerr.cell[i1].doyou_err = NO;
	part = *reference;
	part.hat = -2; //no acceleration, no output
	part_noerr = *reference;
	part_noerr.hat = -2; //no acceleration, no output
	
	wfile = fopen(txtfile, "w");
	if(find_closed_orbite_xxp_zzp(&part, &(part.x), &(part.ux), &(part.uy), &(part.z), &(part.uz), prec_clo, &tempo_latt, YES) != TRUE) errorstop("closed orbit 1 not found");
	if(find_closed_orbite_xxp_zzp(&part_noerr, &(part_noerr.x), &(part_noerr.ux), &(part_noerr.uy), &(part_noerr.z), &(part_noerr.uz), prec_clo, &tempo_latt_noerr, YES) != TRUE) errorstop("closed orbit 2 not found");
	
	s0 = 0;
	th0 = 0;
	for(i1 = 0; i1 < tempo_latt.nbcell; i1++) {
		test_part = part;
		test_part_noerr = part_noerr;
		tempo_latt.cell[i1].instrutype = CUP;
		tempo_latt.cell[i1].instru.ymax = 0;
		tempo_latt.cell[i1].instru.thmax = 0;
		tempo_latt_noerr.cell[i1].instrutype = CUP;
		tempo_latt_noerr.cell[i1].instru.ymax = 0;
		tempo_latt_noerr.cell[i1].instru.thmax = 0;
		
		step_th = tempo_latt.cell[i1].boun.thmax/(nsteps_cell-1.);
		step_y = tempo_latt.cell[i1].boun.ymax/(nsteps_cell-1.);
		
		part_cross_latt(&test_part, &tempo_latt,NULL); 
		part_cross_latt(&test_part_noerr, &tempo_latt_noerr,NULL); 
		find_x(&x, &bias, &test_part_noerr, &test_part);
		find_z(&z, &bias, &test_part_noerr, &test_part);
		fprintf(wfile, "%le	%le	%le	%le\n", s0, th0, x, z);
		printf("%le	%le	%le	%le\n", s0, th0, x, z);
		
		gene_latt_1cell(&latt_1cell, &tempo_latt.cell[i1]);
		gene_latt_1cell(&latt_1cell_noerr, &tempo_latt_noerr.cell[i1]);
		
		for(i2 = 1; i2 < nsteps_cell ; i2++) {
			test_part2 = test_part;
			test_part2_noerr = test_part_noerr;
			latt_1cell.cell[0].instru.thmax = i2*step_th;
			latt_1cell_noerr.cell[0].instru.thmax = i2*step_th;
			s = (sqrt(pow(test_part.x,2) + pow(test_part.y,2))+latt_1cell.cell[0].deltar)*i2*step_th;
			th = i2*step_th;
			part_cross_latt(&test_part2, &latt_1cell,NULL); 
			part_cross_latt(&test_part2_noerr, &latt_1cell_noerr,NULL); 
			find_x(&x, &bias, &test_part2_noerr, &test_part2);
			find_z(&z, &bias, &test_part2_noerr, &test_part2);
			fprintf(wfile, "%le	%le	%le	%le\n", s0+s, th0+th, x, z);
			printf("%le	%le	%le	%le\n", s0+s, th0+th, x, z);
			
		}
		
		test_part = part;
		test_part_noerr = part_noerr;
		part_cross_latt(&test_part, &tempo_latt,NULL);
		part_cross_latt(&test_part_noerr, &tempo_latt_noerr,NULL); 
		tempo_latt.cell[i1].instrutype = NO;
		tempo_latt_noerr.cell[i1].instrutype = NO;
		s0 += s;
		th0 += tempo_latt.cell[i1].boun.thmax;
		printf("cell number %i/%i (s0 = %.5lf [m])\n", i1+1, latt->nbcell, s0);
	}
	fclose(wfile);
	free(tempo_latt.cell);
	free(tempo_latt_noerr.cell);
	
}

//copy the error parameters from cell 
extern void copy_error_latt(struct Lattice *copy_latt, struct Lattice *ori_latt)
{
	int i, j, npole;
	if(ori_latt->nbcell != copy_latt->nbcell) errorstop("nb of cells do not match");
	for(i=0;i<copy_latt->nbcell;i++) {
		if(ori_latt->cell[i].doyou_err == YES) {
			copy_latt->cell[i].doyou_err = YES;
			if(ori_latt->cell[i].nbcomp != copy_latt->cell[i].nbcomp) errorstop("nb of poles do not match");
			for(npole = 0; npole<ori_latt->cell[i].nbcomp; npole++) {
				for(j=0;j<11;j++) copy_latt->cell[i].alierror[npole][j] = ori_latt->cell[i].alierror[npole][j];
			}
		}
	}
}

extern void da_alignment_shit(char *txtnoerr, char *txtfile, char *txtout, double trans_shift, double rot_shift)
{
	long seed;
	int n, nblines;
	double x_noerr, ux_noerr, z_noerr, uz_noerr, qx_noerr, qz_noerr, beta_x_noerr, alpha_x_noerr, beta_z_noerr, alpha_z_noerr, aumax_noerr, avmax_noerr, da_offset;
	double x, ux, z, uz, qx, qz, beta_x, alpha_x, beta_z, alpha_z, aumax, avmax, cod, eps_umax, eps_vmax;
	FILE *rfile_noerr=NULL;
	FILE *rfile_err=NULL;
	FILE *wfile;
	
	nblines = get_nb_lines_file(txtfile);
	
	rfile_noerr = fopen(txtnoerr, "r");
	if(rfile_noerr==NULL) errorstop("cannot open no error file");
	fscanf(rfile_noerr, "%ld	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le", &seed, &x_noerr, &ux_noerr, &z_noerr, &uz_noerr, &qx_noerr, &qz_noerr, &beta_x_noerr, &alpha_x_noerr, &beta_z_noerr, &alpha_z_noerr, &aumax_noerr, &avmax_noerr, &da_offset);
	fclose(rfile_noerr);
	
	rfile_err = fopen(txtfile, "r");
	if(rfile_err==NULL) errorstop("cannot open error file");
	wfile = fopen(txtout, "w");
	for(n=0;n<nblines;n++) {
		fscanf(rfile_err, "%ld	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le", &seed, &x, &ux, &z, &uz, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, &aumax, &avmax, &da_offset);
		eps_umax = aumax*aumax/beta_x_noerr;
		eps_vmax = avmax*avmax/beta_z_noerr;
		cod = sqrt(pow(x-x_noerr,2)+pow(z-z_noerr,2));
		fprintf(wfile, "%ld	%le	%le	%le	%le	%le	%le	%le\n", seed, trans_shift, rot_shift, cod, qx-qx_noerr, qz-qz_noerr, eps_umax, eps_vmax);
	}
	fclose(rfile_err);
	fclose(wfile);
}

extern void compute_dif_loc_m_vffa(char *txt_prefix, double x, double ymin, double ymax, int nstep_y, double zmin, double zmax, int nstep_z, double m_ref, struct Cell *cell, void (*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	char name_long[300], name_vert[300];
	int i,j;
	double bx1,bx2,by1,by2,bz1,bz2, dif_m_long, dif_m_vert,y,z,dz,ystep,zstep;
	FILE *wfile_long;
	FILE *wfile_vert;
	
	sprintf(name_long,"%sgrad_long.dat", txt_prefix);
	wfile_long = fopen(name_long, "w");
	
	sprintf(name_vert,"%sgrad_vert.dat", txt_prefix);
	wfile_vert = fopen(name_vert, "w");
	
	ystep = comp_step(ymin, ymax, nstep_y);
	zstep = comp_step(zmin, zmax, nstep_z);
	dz = zstep;
	for(i=0;i<nstep_z+1;i++) {
		z = zmin + i*zstep;
		for(j=0;j<nstep_y+1;j++) {
			y = ymin + j*ystep;
			printf("x=%le,y=%le,z=%le\n",x,y,z);
			if(get_bfield(x, y, z, &bx1, &by1, &bz1, cell, add_contribution_comp) != ALIVE) errorstop("compute_vertical_field_gradient: get_bfield != ALIVE\n");
			printf("bx1=%le,by1=%le,bz1=%le\n",bx1,by1,bz1);
			if(get_bfield(x, y, z+dz, &bx2, &by2, &bz2, cell, add_contribution_comp) != ALIVE) errorstop("compute_vertical_field_gradient: get_bfield != ALIVE\n");
			printf("bx=%le,by2=%le,bz2=%le\n",bx2,by2,bz2);
			//dif_m_hor = (bx2-bx1)/dz-m_ref*bx1;
			dif_m_long = ((by2-by1)/dz-m_ref*by1)/exp(m_ref*(z-zmin));
			dif_m_vert = ((bz2-bz1)/dz-m_ref*bz1)/exp(m_ref*(z-zmin));
			fprintf(wfile_long, "%le	%le	%le\n", y, z, dif_m_long);
			fprintf(wfile_vert, "%le	%le	%le\n", y, z, dif_m_vert);
		}
		fprintf(wfile_long, "\n");
		fprintf(wfile_vert, "\n");
	}
	
	fclose(wfile_long);
	fclose(wfile_vert);
}

extern void fets_fieldmap_study(char *latt_file, char *txt_prefix)
{
	char shell_write[500];
	char text_long[300], text_vert[300], out_long[300], out_vert[300], xrange[100], yrange[100];
	char txt_bx[300], txt_by[300],txt_bz[300], out_bx[300], out_by[300],out_bz[300];
	char txt_m_intbz[300], txt_m_bz[300], txt_m_intby[300], txt_m_by[300], out_m_int[300], out_m[300];
	int nby, nbz;
	double x0, ymin, ymax, zmin, zmax, m_ref, min, max;
	struct Lattice latt;
	FILE *shell_read;
	
	sprintf(shell_write,"mkdir %s", txt_prefix);
	shell_read = popen(shell_write,"r");
	pclose(shell_read);
	
	
	sprintf(text_long,"%sgrad_long.dat", txt_prefix);
	sprintf(text_vert,"%sgrad_vert.dat", txt_prefix);
	sprintf(out_long,"%sgrad_long.eps", txt_prefix);
	sprintf(out_vert,"%sgrad_vert.eps", txt_prefix);
	sprintf(out_bx,"%s_bx.eps", txt_prefix);
	sprintf(out_by,"%s_by.eps", txt_prefix);
	sprintf(out_bz,"%s_bz.eps", txt_prefix);
	sprintf(txt_bx,"%s_bx.dat", txt_prefix);
	sprintf(txt_by,"%s_by.dat", txt_prefix);
	sprintf(txt_bz,"%s_bz.dat", txt_prefix);
	
	sprintf(txt_m_intbz,"%sm_intbz.dat", txt_prefix);
	sprintf(txt_m_bz,"%sm_bz.dat", txt_prefix);
	sprintf(txt_m_intby,"%sm_intby.dat", txt_prefix);
	sprintf(txt_m_by,"%sm_by.dat", txt_prefix);
	sprintf(out_m_int,"%sm_loc_int.eps", txt_prefix);
	sprintf(out_m,"%sm_loc.eps", txt_prefix);
	
	
	load_lattice(&latt, latt_file);
	//set_vffa_rect_convergence_limit(&(latt.cell[0]), NO, 0.5, add_field_comp_VFFA_rect_enge_add);
	
	x0 = 0;
	
	ymin = 0;
	ymax = latt.cell[0].boun.ymax;
	nby = 113;
	
	zmin = -0.5;
	zmax = -0.1;
	nbz = 81;
	
	m_ref = 1.25;
	
	m_loc_int_bz_cell_vffa(txt_m_intbz, x0, ymin, ymax, nby, zmin, zmax, nbz, 30, m_ref, &(latt.cell[0]), get_bfield_cartesian_map);
	m_loc_int_bz_cell_vffa(txt_m_bz, x0, ymax/2., ymax/2., 1, zmin, zmax, nbz, 30, m_ref, &(latt.cell[0]), get_bfield_cartesian_map);
	m_loc_int_by_cell_vffa(txt_m_intby, x0, ymin, ymax, nby, zmin, zmax, nbz, 30, m_ref, &(latt.cell[0]), get_bfield_cartesian_map);
	m_loc_int_by_cell_vffa(txt_m_by, x0, 0.625, 0.625, 1, zmin, zmax, nbz, 30, m_ref, &(latt.cell[0]), get_bfield_cartesian_map);
	easyplot2(txt_m_intbz, txt_m_intby, "1", "3", "1", "3", "lines lw 4 lc 1", "lines lt 1 lc 0 lw 4", "vertical", "longitudinal", "z [m]", "local m [/m]", NULL, NULL, out_m_int, "grid\nset mxtics 5\nset mytics 5");
	easyplot2(txt_m_bz, txt_m_by, "1", "3", "1", "3", "lines lw 4 lc 1", "lines lt 1 lc 0 lw 4", "vertical (magnet centre)", "longitudinal (magnet edge)", "z [m]", "local m [/m]", NULL, NULL, out_m, "grid\nset mxtics 5\nset mytics 5");
	
	
	compute_fieldheatmap(txt_prefix, &(latt.cell[0]), get_bfield_cartesian_map, 0, 0, 1, ymin, ymax, nby, zmin, zmax, nbz, NO);
	find_boundaries_heated_map_file2(txt_bx, xrange, yrange, &min, &max);
	easyplot3dmap(txt_bx, "longitudinal [m]", "vertical [m]", "bx [T]", xrange, yrange, out_bx, "size ratio -1\nset mxtics 5\nset mytics 5\nset mxtics 5\nset ytics 0.1", min, max); //\nset cbtics 0.02
	find_boundaries_heated_map_file2(txt_by, xrange, yrange, &min, &max);
	easyplot3dmap(txt_by, "longitudinal [m]", "vertical [m]", "by [T]", xrange, yrange, out_by, "size ratio -1\nset mxtics 5\nset mytics 5\nset mxtics 5\nset ytics 0.1", min, max); //\nset cbtics 0.02
	find_boundaries_heated_map_file2(txt_bz, xrange, yrange, &min, &max);
	easyplot3dmap(txt_bz, "longitudinal [m]", "vertical [m]", "bz [T]", xrange, yrange, out_bz, "size ratio -1\nset mxtics 5\nset mytics 5\nset mxtics 5\nset ytics 0.1", min, max); //\nset cbtics 0.02
	
	
	compute_dif_loc_m_vffa(txt_prefix, x0, ymin, ymax, nby, zmin, zmax, nbz, m_ref, &(latt.cell[0]), get_bfield_cartesian_map);
	find_boundaries_heated_map_file2(text_long, xrange, yrange, &min, &max);
	easyplot3dmap(text_long, "longitudinal [m]", "vertical [m]", "difference local longitudinal gradient [T/m]", xrange, yrange, out_long, "size ratio -1\nset mxtics 5\nset mytics 5\nset mxtics 5\nset ytics 0.1", min, max); //\nset cbtics 0.02
	find_boundaries_heated_map_file2(text_vert, xrange, yrange, &min, &max);
	easyplot3dmap(text_vert, "longitudinal [m]", "vertical [m]", "difference local vertical gradient [T/m]", xrange, yrange, out_vert, "size ratio -1\nset mxtics 5\nset mytics 5\nset mxtics 5\nset ytics 0.1", min, max); //\nset cbtics 0.02
	
	free_latt(&latt);
}

extern void aspect_rect_vffa3d_fig_paper_f(char *filename, struct Cell *cell, int doyouerr)
{
	char namef[300];
	int j, doyouplotcomp=YES;
	double x, y,z, zmin,zmax, x_0, y_0, ffbe, ffbs, lambda_e, angle, xwerr,ywerr,zwerr;//, lambda_s;
	FILE *wfile;
	
	sprintf(namef,"%sf.dat", filename);
	wfile = fopen(namef, "a");
	for(j = 0; j < cell->nbcomp; j++) {
		if(j==1 || j==3) {
			if(cell->mpara[j][0]>0-TINYLENGTH && cell->mpara[j][0]<cell->boun.thmax+TINYLENGTH) doyouplotcomp = YES;
			else doyouplotcomp = NO;
			if(doyouplotcomp==YES) {
				x_0 = cell->mpara[j][1]*cos(cell->mpara[j][0]);
				y_0 = cell->mpara[j][1]*sin(cell->mpara[j][0]);
				zmin = cell->collim.zmin;
				zmax = cell->collim.zmax;

				ffbe = cell->efben[j][0];
				ffbs = cell->efbex[j][0];
				//lambda_e = cell->efben[j][1];
				//lambda_s = cell->efbex[j][1];
				lambda_e = cell->mpara[j][6];
				angle = cell->mpara[j][0]+cell->mpara[j][5];
			
				x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle); //D
				y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
		
				x = x_0 - ffbe*sin(angle) + lambda_e*cos(angle); //C
				y = y_0 + ffbe*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle); //B
				y = y_0 + ffbe*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle); //E
				y = y_0 + ffbs*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle); //B
				y = y_0 + ffbe*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle) - lambda_e*cos(angle); //A
				y = y_0 + ffbe*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
		
				x = x_0 - ffbs*sin(angle) - lambda_e*cos(angle); //F
				y = y_0 + ffbs*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
		
				x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle); //D
				y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			
			
				x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle); //J
				y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
		
				x = x_0 - ffbe*sin(angle) + lambda_e*cos(angle); //I
				y = y_0 + ffbe*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle) + lambda_e*cos(angle); //C
				y = y_0 + ffbe*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle) + lambda_e*cos(angle); //I
				y = y_0 + ffbe*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle); //H
				y = y_0 + ffbe*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle); //B
				y = y_0 + ffbe*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle); //H
				y = y_0 + ffbe*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle); //K
				y = y_0 + ffbs*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle); //E
				y = y_0 + ffbs*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle); //K
				y = y_0 + ffbs*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle); //H
				y = y_0 + ffbe*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle) - lambda_e*cos(angle); //G
				y = y_0 + ffbe*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle) - lambda_e*cos(angle); //A
				y = y_0 + ffbe*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle) - lambda_e*cos(angle); //G
				y = y_0 + ffbe*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle) - lambda_e*cos(angle); //L
				y = y_0 + ffbs*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle) - lambda_e*cos(angle); //F
				y = y_0 + ffbs*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle) - lambda_e*cos(angle); //L
				y = y_0 + ffbs*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle); //J
				y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
				fprintf(wfile, "\n\n");
			}
		}
	}
	fclose(wfile);
}

extern void aspect_rect_vffa3d_fig_paper_d(char *filename, struct Cell *cell, int doyouerr)
{
	char namef[300];
	int j, doyouplotcomp=YES;
	double x, y,z, zmin,zmax, x_0, y_0, ffbe, ffbs, lambda_e, angle, xwerr,ywerr,zwerr;//, lambda_s;
	FILE *wfile;
	
	sprintf(namef,"%sd.dat", filename);
	wfile = fopen(namef, "a");
	for(j = 0; j < cell->nbcomp; j++) {
		if(j==2) {
			if(cell->mpara[j][0]>0-TINYLENGTH && cell->mpara[j][0]<cell->boun.thmax+TINYLENGTH) doyouplotcomp = YES;
			else doyouplotcomp = NO;
			if(doyouplotcomp==YES) {
				x_0 = cell->mpara[j][1]*cos(cell->mpara[j][0]);
				y_0 = cell->mpara[j][1]*sin(cell->mpara[j][0]);
				zmin = cell->collim.zmin;
				zmax = cell->collim.zmax;

				ffbe = cell->efben[j][0];
				ffbs = cell->efbex[j][0];
				//lambda_e = cell->efben[j][1];
				//lambda_s = cell->efbex[j][1];
				lambda_e = cell->mpara[j][6];
				angle = cell->mpara[j][0]+cell->mpara[j][5];
			
				x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle); //D
				y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
		
				x = x_0 - ffbe*sin(angle) + lambda_e*cos(angle); //C
				y = y_0 + ffbe*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle); //B
				y = y_0 + ffbe*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle); //E
				y = y_0 + ffbs*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle); //B
				y = y_0 + ffbe*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle) - lambda_e*cos(angle); //A
				y = y_0 + ffbe*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
		
				x = x_0 - ffbs*sin(angle) - lambda_e*cos(angle); //F
				y = y_0 + ffbs*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
		
				x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle); //D
				y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
			
			
				x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle); //J
				y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
		
				x = x_0 - ffbe*sin(angle) + lambda_e*cos(angle); //I
				y = y_0 + ffbe*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle) + lambda_e*cos(angle); //C
				y = y_0 + ffbe*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle) + lambda_e*cos(angle); //I
				y = y_0 + ffbe*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle); //H
				y = y_0 + ffbe*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle); //B
				y = y_0 + ffbe*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle); //H
				y = y_0 + ffbe*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle); //K
				y = y_0 + ffbs*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle); //E
				y = y_0 + ffbs*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle); //K
				y = y_0 + ffbs*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle); //H
				y = y_0 + ffbe*cos(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle) - lambda_e*cos(angle); //G
				y = y_0 + ffbe*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle) - lambda_e*cos(angle); //A
				y = y_0 + ffbe*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbe*sin(angle) - lambda_e*cos(angle); //G
				y = y_0 + ffbe*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle) - lambda_e*cos(angle); //L
				y = y_0 + ffbs*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle) - lambda_e*cos(angle); //F
				y = y_0 + ffbs*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmin, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmin;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle) - lambda_e*cos(angle); //L
				y = y_0 + ffbs*cos(angle) - lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
			
				x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle); //J
				y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
				if(doyouerr==YES) {
					apply_alierror_comp(cell, j, x, y, zmax, &xwerr, &ywerr, &zwerr);
					x = xwerr;
					y = ywerr;
					z = zwerr;
				}
				else z = zmax;
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le  %le\n",x,y,z);
				fprintf(wfile, "\n\n");
			}
		}
	}
	fclose(wfile);
}

extern void config_maker(char *name, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7)
{
	FILE *wfile=NULL;
	
	wfile = fopen("config.txt", "w");
	
	fprintf(wfile, "Parameters name: %s\n", name);
	fprintf(wfile, "Scaling k or m (m^-1): 1.28\n");
    fprintf(wfile, "\n");
	fprintf(wfile, "Good field height (m): 0.6\n");
	fprintf(wfile, "Vertical end height (m): 0.75\n");
	fprintf(wfile, "Vertical middle join height (m): 0.5\n");
	fprintf(wfile, "Vertical end asymmetry: 0\n");
	fprintf(wfile, "Vertical middle asymmetry: 0\n");
	fprintf(wfile, "Vertical end ramp: 2.5\n");
	fprintf(wfile, "Vertical middle ramp: 1.6\n");
    fprintf(wfile, "\n");
	fprintf(wfile, "Coil meshes: 100\n");
	fprintf(wfile, "FFT meshes (power of two): 1024\n");
	fprintf(wfile, "Use vertical FFT inversion (Y/N): Y\n");
    fprintf(wfile, "\n");
	fprintf(wfile, "Fieldmap transverse cell size (m): 0.005\n");
	fprintf(wfile, "Fieldmap vertical cell size (m): 0.005\n");
	fprintf(wfile, "Fieldmap longitudinal cell size (m): 0.02\n");
	fprintf(wfile, "Fieldmap transverse extension beyond magnet (m): 0.05\n");
	fprintf(wfile, "Fieldmap vertical extension beyond magnet (m): 0.05\n");
	fprintf(wfile, "Fieldmap longitudinal extension beyond magnet (m): 0.625\n");
    fprintf(wfile, "\n");
	fprintf(wfile, "Magnet name: Bf\n");
	fprintf(wfile, "Bmax (T): 1.2\n");
	fprintf(wfile, "Length (m): 1\n");
	fprintf(wfile, "Angle (degrees): 0\n");
	fprintf(wfile, "Winding centre to midplane distance (m): 0.12\n");
	fprintf(wfile, "Winding return half-length (m): 0.04\n");
	fprintf(wfile, "Area per winding (cm^2): 0.25\n");
	fprintf(wfile, "Winding current density (A/mm^2): 250\n");
	fprintf(wfile, "Winding thickness (cm): 2\n");
	fprintf(wfile, "Polynomial adjustment (coeffs x^0,x^1...): %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n\n",x0,x1,x2,x3,x4,x5,x6,x7);
	
	fclose(wfile);
}

extern void write_maplatt_file(char *lattfile, char *title, int periodicity, char *mapfile, double stepsize, double cell_length, double rmin, double rmax, int nbr, double thmin, double thmax, int nbth, double zmin, double zmax, int nbz)
{
	FILE *wfile=NULL;
	
	wfile = fopen(lattfile, "w");
	if(wfile==NULL) errorstop("cannot create wfile\n");
	
	fprintf(wfile, "%s\n\n", title);
	fprintf(wfile, "%i				//lattice periodicity [int]\n", periodicity);
	fprintf(wfile, "1				//number of elements [int]\n");
	fprintf(wfile, "1				//number of cell types [int]\n\n\n");
	//fprintf(wfile, "1 field-cylmap	nosym	cylrad	%s\n", mapfile);
	fprintf(wfile, "1 field-cartmap	nosym	cart %s\n", mapfile);
	fprintf(wfile, "%lf			//tracking step size [m]\n", stepsize);
	fprintf(wfile, "%lf						//cell total length [m]\n", cell_length);
	fprintf(wfile, "1	1		//unitlength/m (1000=[mm], 100=[cm]); B unit (1.e-4 Gauss, 1: Tesla)\n");
	//fprintf(wfile, "312	123			//loop order, line order (123:(x,y,z);231:(y,z,x))\n");
	fprintf(wfile, "132	132			//loop order, line order (123:(x,y,z);231:(y,z,x))\n");
	//fprintf(wfile, "2				// header skip (nb of lines jumped)\n");
	fprintf(wfile, "1				// header skip (nb of lines jumped)\n");
	fprintf(wfile, "0				//delta x (not checked!) (in map unitlength)\n\n");
	fprintf(wfile, "%lf		%lf		%i			// rmin, rmax (in map unitlength), nbsteps_r\n", rmin, rmax, nbr);
	fprintf(wfile, "%lf		%lf		%i			// thmin, thmax (in [deg.]), nbsteps_th\n", thmin, thmax, nbth);
	fprintf(wfile, "%lf		%lf		%i			// zmin, zmax (in map unitlength), nbsteps_z\n\n", zmin, zmax, nbz);
	fclose(wfile);
}

extern void compute_vffa_error_m_int_bz(double *error_min, double *error_max, double x0, double ymin, double ymax, int nstep_y, double zmin, double zmax, int nstep_z, int nzref, double mref, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double zstep, m, z[nstep_z], intbdl[nstep_z], intbdlabs[nstep_z], errorm, btemp1, btemp2;
	
	//calcul de BL
	zstep = (zmax - zmin)/(nstep_z-1);
	for(i = 0; i < nstep_z; i++) {
		z[i] = zmin + i*zstep;
		if(ymax-ymin>TINYLENGTH) int_bz_liney_cell2(&(intbdl[i]), &(intbdlabs[i]), x0, ymin, ymax, nstep_y, z[i], cell, add_contribution_comp);
		else if(get_bfield(x0, ymin, z[i], &btemp1, &btemp2, &(intbdl[i]), cell, add_contribution_comp) != ALIVE) errorstop("m_loc_int_bz_cell_vffa: get_bfield != ALIVE\n");
	}
	
	*error_min = 1.e9;
	*error_max = -1.e9;
	//calcul de m_loc
	for(i = 0; i < nstep_z; i++) {
		if(i != nzref) {
			m = log(intbdl[i]/intbdl[nzref])/(z[i] - z[nzref]);
			errorm = (m - mref)/mref;
			*error_min = MIN(*error_min, errorm);
			*error_max = MAX(*error_max, errorm);
		}
	}
}

extern void compute_vffa_error_m_int_by(double *error_min, double *error_max, double x0, double ymin, double ymax, int nstep_y, double zmin, double zmax, int nstep_z, int nzref, double mref, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double zstep, m, z[nstep_z], intbdl[nstep_z], intbdlabs[nstep_z], errorm, btemp1, btemp2;
	
	//calcul de BL
	zstep = (zmax - zmin)/(nstep_z-1);
	for(i = 0; i < nstep_z; i++) {
		z[i] = zmin + i*zstep;
		if(ymax-ymin>TINYLENGTH) int_by_liney_cell(&(intbdl[i]), &(intbdlabs[i]), x0, ymin, ymax, nstep_y, z[i], cell, add_contribution_comp);
		else if(get_bfield(x0, ymin, z[i], &btemp1, &(intbdl[i]), &btemp2, cell, add_contribution_comp) != ALIVE) errorstop("m_loc_int_bz_cell_vffa: get_bfield != ALIVE\n");
	}
	
	*error_min = 1.e9;
	*error_max = -1.e9;
	//calcul de m_loc
	for(i = 0; i < nstep_z; i++) {
		if(i != nzref) {
			m = log(intbdl[i]/intbdl[nzref])/(z[i] - z[nzref]);
			errorm = (m - mref)/mref;
			*error_min = MIN(*error_min, errorm);
			*error_max = MAX(*error_max, errorm);
		}
	}
}

extern void compute_circmap_from_cartmap_superimpose2(char *textfile, double x0, double y0, double r0, double th0_deg, double xshift, double scale1, double scale2, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, struct Cell *cell)
{
	int i,j,k;
	double step_r, step_th, step_z,z,r,th, theta_min, theta_max, theta_tot,bx,by,bz,br,bth, bx1,by1,bz1,th0, bxs,bys;
	FILE *wfile=NULL;
	
	theta_min = theta_min_deg*PI/180.;
	theta_max = theta_max_deg*PI/180.;
	th0 = th0_deg*PI/180.;
	theta_tot = fabs(theta_min)+ fabs(theta_max);
	step_r =  comp_step(rmin, rmax, rstep);
	step_th = comp_step(theta_min, theta_max, thstep);
	step_z =  comp_step(zmin, zmax, zstep);
	
	//test
	if(test_cell_map(cell)==NO) errorstop("cell not a map");
	//if(cell->map.mapdim[1] < rmax-r0+x0) errorstop("rmax too big");
	
	if(cell->map.mapdim[0] > rmin*cos(th0)-r0+x0 || cell->map.mapdim[0] > rmin*cos(theta_max-th0)-r0+x0) {
		printf("cell->map.mapdim[0]=%le, rmin*cos(theta_min)-r0+x0=%le, rmin*cos(theta_max)-r0+x0=%le\n",cell->map.mapdim[0], rmin*cos(th0)-r0+x0, rmin*cos(theta_max-th0)-r0+x0);
		errorstop("rmin too big2");
	}
	if(rmax*sin(theta_min)+y0 < 0) errorstop("theta_min too small");
	//if(rmax*sin(theta_max)+y0 > cell->map.mapdim[3]) errorstop("theta_max too big");
	printf("zmin = %le, zmax = %le, step_z = %le\n", zmin, zmin+(zstep-1)*step_z, step_z);
	printf("rmin = %le, rmax = %le, step_r = %le\n", rmin, rmin+(rstep-1)*step_r, step_r);
	printf("thmin = %le, thmax = %le, step_th = %le\n", theta_min, theta_min+(thstep-1)*step_th, step_th);
	wfile = fopen(textfile, "w");
	fprintf(wfile, "%i %i %i\n",rstep, thstep, zstep);
	fprintf(wfile, "r(m),theta(deg),z(m),Br(T),Bth,Bz\n");
	for(i=0;i<zstep;i++) {
		z = zmin + i*step_z;
		for(j=0;j<rstep;j++) {
			r = rmin + j*step_r;
			for(k=0;k<thstep;k++) {
				th = theta_min + k*step_th;
				map_neighbours2(&bx, &by, &bz, x0-xshift, y0, r0, th0, theta_tot, r, th, z, scale1, cell, YES);
				map_neighbours2(&bx1, &by1, &bz1, x0+xshift, y0, r0, th0, theta_tot, r, th, z, scale2, cell, YES);
				bx+=bx1;by+=by1;bz+=bz1;
				bxs = bx*cos(th0) - by*sin(th0);
				bys = by*cos(th0) + bx*sin(th0);
				br = bxs*cos(th) + bys*sin(th); 
				bth = bys*cos(th) - bxs*sin(th); 
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r, th*180./PI, z, br, bth, bz);
				//fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r, th*180./PI, z, bxs, bys, bz);
			}
		}
	}
	fclose(wfile);
}

//what happens if x0!=0?
//angle in radiant!
extern void map_neighbours2(double *bx, double *by, double *bz, double x0, double y0, double r0, double th0, double theta_tot, double r, double th, double z, double scale, struct Cell *cell, int kill_notalive)
{
	double x,y,znew,x1,x2,y1,y2,bx0=0,by0=0,bz0=0,bx1=0,by1=0,bz1=0,bx2=0,by2=0,bz2=0;
	
	convert_ff_coord_to_shinji(&x, &y, &znew, r, th, z, r0, th0, y0, 0.);
	//x = r*cos(th-th0)-r0+x0;
	//y = r*sin(th-th0)+y0;
	convert_ff_coord_to_shinji(&x1, &y1, &znew, r, th, z, r0, th0-theta_tot, y0, 0.);
	//x1 = r*cos(th-th0+theta_tot)-r0+x0;
	//y1 = r*sin(th-th0+theta_tot)+y0;
	convert_ff_coord_to_shinji(&x2, &y2, &znew, r, th, z, r0, th0+theta_tot, y0, 0.);
	//x2 = r*cos(th-th0-theta_tot)-r0+x0;
	//y2 = r*sin(th-th0-theta_tot)+y0;
	
	//printf("r=%le, th=%le, (x,y)=%le,%le\n", r, th*180./PI, x,y);
	//printf("r=%le, th=%le, (x,y)=%le,%le, (x1,y1)=%le,%le\n", r, th*180./PI, x,y, x1,y1);
	
	if(get_bfield(x, y, z, &bx0, &by0, &bz0, cell, get_bfield_cartesian_map) != ALIVE) {
		//printf("\nr=%le, th=%le, (x,y)=%le,%le, (x1,y1)=%le,%le, (x2,y2)=%le,%le\n", r, th*180./PI, x,y, x1,y1, x2, y2);
		//printf("b0=(%le,%le,%le)\n",bx0,by0,bz0);
		if(kill_notalive==YES) errorstop("map_neighbours: get_bfield != ALIVE\n");
		else {
			bx0=0;by0=0;bz0=0;
		}
	}
	//printf("\nr=%le, th=%le, (x,y,z)=%le,%le,%le, (x1,y1)=%le,%le, (x2,y2)=%le,%le\n", r, th*180./PI, x,y,z, x1,y1, x2, y2);
	//printf("b0=(%le,%le,%le)\n",bx0,by0,bz0);
	bx0*=scale;by0*=scale;bz0*=scale;
	*bx = bx0; *by = by0; *bz = bz0;
	//printf("b0=(%le,%le,%le)\n",bx0,by0,bz0);
	//printf("r=%le, th=%le, (x,y)=%le,%le, (x1,y1)=%le,%le", r, th*180./PI, x,y, x1,y1);
	if(get_bfield(x1, y1, z, &bx1, &by1, &bz1, cell, get_bfield_cartesian_map) == ALIVE) { //cell before
		//printf("b1=(%le,%le,%le)\n",bx1,by1,bz1);
		bx1*=scale;by1*=scale;bz1*=scale;
		*bx += bx1*cos(theta_tot) + by1*sin(theta_tot);
		*by += by1*cos(theta_tot) - bx1*sin(theta_tot);
		*bz += bz1;
		//printf("   YES");
	}
	//printf("b=(%le,%le,%le)\n", bx,by,bz);
	if(get_bfield(x2, y2, z, &bx2, &by2, &bz2, cell, get_bfield_cartesian_map) == ALIVE) { //cell after
		//printf("b2=(%le,%le,%le)\n",bx2,by2,bz2);
		bx2*=scale;by2*=scale;bz2*=scale;
		*bx += bx2*cos(theta_tot) - by2*sin(theta_tot);
		*by += by2*cos(theta_tot) + bx2*sin(theta_tot);
		*bz += bz2;
	}//*/
	//printf("tot (%le,%le,%le), %le + %le\n", *bx,*by,*bz,  bx2*cos(theta_tot), -by2*sin(theta_tot));
}

// remove the change of coordinates for the central magnet

extern void compute_vffacircmap_doublecoil_singlemagnet2(char *cartmapfile, char *circmapfile, double x0, double y0, double r0, double th0_deg, double xshift, double scale1, double scale2, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep)
{
	struct Lattice latt;
	
	load_lattice(&latt, cartmapfile);
	compute_circmap_from_cartmap_superimpose2(circmapfile, x0, y0, r0, th0_deg, xshift, scale1, scale2, rmin, rmax, rstep, theta_min_deg, theta_max_deg, thstep, zmin, zmax, zstep, &(latt.cell[0]));
	free_latt(&latt);
}

extern void compute_cartmap_superimpose_shift(char *cartmapfile, char *cartmapfile_final, double xshift, double scale1, double scale2, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz)
{
	int i,j,k;
	double xstep, ystep, zstep, x,y,z,bx2,by2,bz2,bx,by,bz;
	struct Lattice latt;
	FILE *wfile;
	
	xstep = comp_step(xmin, xmax, nbx);
	ystep = comp_step(ymin, ymax, nby);
	zstep = comp_step(zmin, zmax, nbz);
	
	load_lattice(&latt, cartmapfile);
	wfile = fopen(cartmapfile_final, "w");
	fprintf(wfile, "%i  %i  %i\n",nbx,nby,nbz);
	for(k=0;k<nbz;k++) {
		z = zmin+k*zstep;
		for(i=0;i<nbx;i++) {
			x = xmin+i*xstep;
			for(j=0;j<nby;j++) {
				y = ymin+j*ystep;
				bx=0;by=0;bz=0;
				if(get_bfield(x-xshift, y, z, &bx, &by, &bz, &(latt.cell[0]), get_bfield_cartesian_map) == ALIVE) {
					bx*=scale1; by*=scale1; bz*=scale1;
				}
				if(get_bfield(x+xshift, y, z, &bx2, &by2, &bz2, &(latt.cell[0]), get_bfield_cartesian_map) == ALIVE) {
					bx+=bx2*scale2; by+=by2*scale2; bz+=bz2*scale2;
				}
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", x,y,z,bx,by,bz);
			}
		}
	}
	fclose(wfile);
	free_latt(&latt);
}

extern void check_shinji_file(char *fileout, char *shinji_file, struct Cell *cell, double zshift, double yshift, double shinji_scale)
{
	int i,nblines;
	double x,px,y,py,z,pz,bx,by,bz,a,b,c,d, bx_ff, by_ff, bz_ff;
	FILE *rfile=NULL;
	FILE *wfile;
	nblines = get_nb_lines_file(shinji_file);
	
	rfile = fopen(shinji_file, "r");
	if(rfile==NULL) errorstop("cannot open shinji's file");
	wfile = fopen(fileout, "w");
	
	newline(rfile);
	for(i=0;i<nblines-1;i++) {
		fscanf(rfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le", &y, &py, &x, &px, &z, &pz, &a, &by, &bx, &bz, &b,&c,&d);
		//printf("%le  %le  %le  %le  %le  %le\n", x,y,z,bx,by,bz);
		bx*=shinji_scale; by*=shinji_scale; bz*=shinji_scale;
		if(get_bfield(x, y-yshift, z-zshift, &bx_ff, &by_ff, &bz_ff, cell, get_bfield_cartesian_map) != ALIVE) errorstop("problem field");
		//printf("z=%le\n",z+zshift);
		//fprintf(wfile, "%le	%le	%le	%le\n", y, bx/bx_ff, by/by_ff, bz/bz_ff);
		//fprintf(wfile, "%le	%le	%le	%le\n", y, (bx_ff-bx)/bx, (by_ff-by)/by, (bz_ff-bz)/bz);
		fprintf(wfile, "%le	%le	%le	%le\n", y, bx_ff, by_ff, bz_ff);
	}
	
	fclose(rfile);
	fclose(wfile);
}

extern void check_shinji_file_circmap(char *fileout, char *shinji_file, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double zshift, double yshift, double shinji_scale, double r0, double th0_deg)
{
	int i,nblines;
	double x,px,y,py,z,pz,bx,by,bz,a,b,c,d, bx_ff, by_ff, bz_ff,r,th,xshinji,yshinji, zshinji, bshinjix,bshinjiy, th0=th0_deg*PI/180.;
	FILE *rfile=NULL;
	FILE *wfile;
	nblines = get_nb_lines_file(shinji_file);
	
	rfile = fopen(shinji_file, "r");
	if(rfile==NULL) errorstop("cannot open shinji's file");
	wfile = fopen(fileout, "w");
	
	newline(rfile);
	for(i=0;i<nblines-1;i++) {
		fscanf(rfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le", &yshinji, &py, &xshinji, &px, &zshinji, &pz, &a, &by, &bx, &bz, &b,&c,&d);
		//printf("%le  %le  %le  %le  %le  %le\n", x,y,z,bx,by,bz);
		bx*=shinji_scale; by*=shinji_scale; bz*=shinji_scale;
		convert_shinji_coord_to_ff(&r, &th, &z, xshinji, yshinji, zshinji, r0, th0, yshift, zshift);
		x = r*cos(th);
		y = r*sin(th);
		printf("r=%le, th=%ledeg, x,y,z=%le,%le,%le, xshinji=%le,yshinji=%le\n",r,th*180./PI, x,y,z,xshinji,yshinji);
		if(get_bfield(x, y, z, &bx_ff, &by_ff, &bz_ff, cell, add_contribution_comp) != ALIVE) errorstop("problem field");
		bshinjix = bx_ff*cos(th0) + by_ff*sin(th0);
		bshinjiy = -bx_ff*sin(th0) + by_ff*cos(th0);
		fprintf(wfile, "%le	%le	%le	%le	%le\n", y,z, bshinjix, bshinjiy, bz_ff);
		//fprintf(wfile, "%le	%le	%le	%le\n", y, (bx_ff-bx)/bx, (by_ff-by)/by, (bz_ff-bz)/bz);
	}
	fclose(rfile);
	fclose(wfile);
}

extern void compute_fieldheatmap_shinji(char *filename_prefix, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int),
							double r0, double rmax, int nbr, double th0, double thmax, int nbth, double z0, double zmax, int nbz, double angle_shinji_deg)
{
	char namex[300], namey[300], namez[300];
	int i,j,k;
	double x,y,z,bxtemp,bytemp,bx,by,bz,r,th, rstep,thstep,zstep, w_coor, angle_shinji=angle_shinji_deg*PI/180.;
	FILE *wfilex, *wfiley, *wfilez;
	
	sprintf(namex,"%s_bxshinji.dat", filename_prefix);
	sprintf(namey,"%s_byshinji.dat", filename_prefix);
	sprintf(namez,"%s_bz.dat", filename_prefix);
	wfilex=fopen(namex, "w");
	wfiley=fopen(namey, "w");
	wfilez=fopen(namez, "w");
	
	if(nbr <= 1 && nbz <= 1) errorstop("map not possible, only 1 dimension is scanned!");
	if(nbr != 1 && nbz != 1) errorstop("map not possible, 3 dimensions scanned");
	if(nbth <= 1) errorstop("map not possible, longitudinal dimension not scanned");
	
	zstep = comp_step(z0, zmax, nbz);
	rstep = comp_step(r0, rmax, nbr);
	thstep = comp_step(th0, thmax, nbth);
	printf("xstep=%le[m], ystep=%le[m], zstep=%le[m]\n",rstep,thstep,zstep);
	for(j=0;j<nbr;j++) {
		r = r0 + j*rstep;
		for(i=0;i<nbz;i++) {
			z = z0 + i*zstep;
			for(k=0;k<nbth;k++) {
				th = k*thstep;
				x = r;
				y = th;
				if(get_bfield(x, y, z, &bxtemp, &bytemp, &bz, cell, add_contribution_comp)!=ALIVE) {
					//xerrorstop("field not alive\n");
					bx = 0.;
					by = 0.;
					bz = 0.;
				}
				if(nbz==1) w_coor = x;
				else w_coor = z;
				bx = bxtemp*cos(angle_shinji) + bytemp*sin(angle_shinji); //x (hor) shinji
				by = -bxtemp*sin(angle_shinji) + bytemp*cos(angle_shinji); //y (long) shinji
				
				fprintf(wfilex,"%le\t%le\t%le\n", y, w_coor, bx);
				fprintf(wfiley,"%le\t%le\t%le\n", y, w_coor, by);
				fprintf(wfilez,"%le\t%le\t%le\n", y, w_coor, bz);
			}
			fprintf(wfilex,"\n"); // to build file to plot easyplot3dmap
			fprintf(wfiley,"\n"); // to build file to plot easyplot3dmap
			fprintf(wfilez,"\n"); // to build file to plot easyplot3dmap
		}
	}
	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
}

extern void compute_m_value_vffa_prototype(char *latt_name, char *out_long, char *out_vert, double x0, double ymin, double ymax_vert, int nby_vert, double ymax_long, int nby_long, double zmin, double zmax, int nstep_z, int nzref, double mref)
{
	struct Lattice latt;
	
	load_lattice(&latt, latt_name);
	printf("latt ymax: %le\n", latt.cell[0].boun.ymax);
	m_loc_int_by_cell_vffa(out_long, x0, ymin, ymax_long, nby_long, zmin, zmax, nstep_z, nzref, mref, &(latt.cell[0]), get_bfield_cartesian_map);
	m_loc_int_bz_cell_vffa(out_vert, x0, ymin, ymax_vert, nby_vert, zmin, zmax, nstep_z, nzref, mref, &(latt.cell[0]), get_bfield_cartesian_map);
	
	free_latt(&latt);
}

extern void compute_m_value_vffa_prototype2(char *latt_name, char *out_long, char *out_vert, double x0, double ymin, double ymax_vert, int nby_vert, double ymax_long, int nby_long, double zmin, double zmax, int nstep_z, int nzref, double mref)
{
	struct Lattice latt;
	
	load_lattice(&latt, latt_name);
	printf("latt ymax: %le\n", latt.cell[0].boun.ymax);
	printf("long\n");
	m_loc_int_by_cell_vffa2(out_long, x0, ymin, ymax_long, nby_long, zmin, zmax, nstep_z, nzref, mref, &(latt.cell[0]), get_bfield_cartesian_map);
	printf("vert\n");
	m_loc_int_bz_cell_vffa(out_vert, x0, ymin, ymax_vert, nby_vert, zmin, zmax, nstep_z, nzref, mref, &(latt.cell[0]), get_bfield_cartesian_map);
	
	free_latt(&latt);
}

extern void m_loc_int_by_cell_vffa2(char *outfilename, double x, double ymin, double ymax, int nstep_y, double zmin, double zmax, int nstep_z, int nzref, double mref, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double zstep, m[nstep_z+1], z[nstep_z+1], intbdl[nstep_z+1], intbdlabs[nstep_z+1], errorm[nstep_z+1], errorbdl[nstep_z+1], intbdldes[nstep_z+1], btemp1, btemp2, avg_m=0, min_m=1.e10, max_m=-1.e10;
	FILE *wfile;
	
	zstep = (zmax - zmin)/nstep_z;
	printf("step=%lf\n", zstep);
	for(i = 0; i < nstep_z + 1; i++) {
		z[i] = zmin + i*zstep;
		//printf("%i,z=%le\n",i,z[i]);
		if(ymax-ymin>TINYLENGTH) int_by_liney_cell(&(intbdl[i]), &(intbdlabs[i]), x, ymin, ymax, nstep_y, z[i], cell, add_contribution_comp);
		else if(get_bfield(x, ymin, z[i], &btemp1, &(intbdl[i]), &btemp2, cell, add_contribution_comp) != ALIVE) errorstop("m_loc_int_by_cell_vffa: get_bfield != ALIVE\n");
	}
	
	for(i = 0; i < nstep_z + 1; i++) {
		m[i] = log(intbdlabs[i]/intbdlabs[nzref])/(z[i] - z[nzref]);
		errorm[i] = (m[i] - mref)/mref;
	}
	if(nzref != nstep_z) m[nzref] = (m[nzref - 1] + m[nzref + 1])/2;
	else m[nzref] = m[nzref - 1];
	//plot the desired Bl
	for(i = 0; i < nstep_z + 1; i++) {
		avg_m+=m[i];
		min_m = MIN(min_m, m[i]);
		max_m = MAX(max_m, m[i]);
		if(i != nzref) {
			intbdlabs[i] = intbdlabs[nzref]*exp(mref*(z[i] - z[nzref]));
			errorbdl[i] = (fabs(intbdlabs[i]) - fabs(intbdlabs[i]))/fabs(intbdlabs[i]);
		}
	}
	avg_m/=(double) (nstep_z+1.0);
	printf("avg long m: %le, min: %le/100, max: %le/100\n", avg_m, 100*(min_m-avg_m)/avg_m, 100*(max_m-avg_m)/avg_m);
	
	intbdldes[nzref] = intbdl[nzref];
	errorbdl[nzref] = 0;
	errorm[nzref] = (m[nzref] - mref)/mref;
	//data writing
	//printf("outfile:%s\n", outfilename);
	wfile = fopen(outfilename,"w");
	fprintf(wfile,"\n");
	for(i = 0; i < nstep_z + 1; i++) {
		fprintf(wfile, "%lf  %lf  %lf  %le  %le  %le\n", z[i], intbdl[i], m[i], errorm[i], intbdldes[i], errorbdl[i]);
	}
	fclose(wfile);
}

extern void translate_shinji_vffa_file(char *shinji_file, char *outfile, double zshift, double yshift, double scale, double r0, double th0_deg, double theta_tot_deg)
{
	int i, flag=0, nblines = get_nb_lines_file(shinji_file);
	double r[nblines],th[nblines],x[nblines],px,y[nblines],py,z[nblines],pz,bx[nblines],by[nblines],bz[nblines],a,b,c,d,xshinji[nblines], yshinji[nblines], zshinji[nblines], theta_tot=theta_tot_deg*PI/180., th0=th0_deg*PI/180.;
	FILE *rfile=NULL;
	FILE *wfile;
	
	rfile = fopen(shinji_file, "r");
	if(rfile==NULL) errorstop("cannot open shinji's file");
	
	for(i=0;i<nblines;i++) {
		fscanf(rfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le", &(yshinji[i]), &py, &(xshinji[i]), &px, &(zshinji[i]), &pz, &a, &by[i], &bx[i], &bz[i], &b,&c,&d);
		convert_shinji_coord_to_ff(&(r[i]), &(th[i]), &(z[i]), xshinji[i], yshinji[i], zshinji[i], r0, th0, yshift, zshift);
		
		//yshinji[i]-=yshift;
		//z[i]=zshinji[i]-zshift;
		//r[i] = sqrt(yshinji[i]*yshinji[i] + (xshinji[i]+r0)*(xshinji[i]+r0));
		//th[i] = atan_ratio(yshinji[i], xshinji[i]+r0) + th0;
		if(th[i]>theta_tot) {
			th[i] -= theta_tot;
			if(flag==0) flag = i;
		}
		x[i] = r[i]*cos(th[i]);
		y[i] = r[i]*sin(th[i]);
		//fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le\n", x[i],y[i],z[i],bx*scale,by*scale,bz*scale,r,th*180./PI);
	}
	fclose(rfile);
	
	wfile = fopen(outfile, "w");
	for(i=flag;i<nblines;i++) fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le\n", x[i],y[i],z[i],bx[i]*scale,by[i]*scale,bz[i]*scale,r[i],th[i]*180./PI);
	for(i=0;i<flag;i++) fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le\n", x[i],y[i],z[i],bx[i]*scale,by[i]*scale,bz[i]*scale,r[i],th[i]*180./PI);
	fclose(wfile);
}

//compute map including surrounding rectangular maps
extern void compute_circmap_from_cartmap2(char *textfile, double x0, double y0, double r0, double th0_deg, double scale, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, struct Cell *cell)
{
	int i,j,k;
	double step_r,step_th,step_z, r,th,z, theta_min,theta_max,theta_tot, bx,by,bz,br,bth, bxs,bys,th0;
	FILE *wfile=NULL;
	
	theta_min = theta_min_deg*PI/180.;
	theta_max = theta_max_deg*PI/180.;
	th0 = th0_deg*PI/180.;
	theta_tot = fabs(theta_min)+ fabs(theta_max);
	step_r =  comp_step(rmin, rmax, rstep);
	step_th = comp_step(theta_min, theta_max, thstep);
	step_z =  comp_step(zmin, zmax, zstep);
	
	//test
	if(test_cell_map(cell)==NO) errorstop("cell not a map");
	//if(cell->map.mapdim[1] < rmax-r0+x0) errorstop("rmax too big");
	
	
	if(cell->map.mapdim[0] > rmin*cos(th0)-r0+x0 || cell->map.mapdim[0] > rmin*cos(theta_max-th0)-r0+x0) {
		printf("cell->map.mapdim[0]=%le, rmin*cos(theta_min)-r0+x0=%le, rmin*cos(theta_max)-r0+x0=%le\n",cell->map.mapdim[0], rmin*cos(th0)-r0+x0, rmin*cos(theta_max-th0)-r0+x0);
		errorstop("rmin too big3");
	}
	printf("zmin = %le, zmax = %le, step_z = %le\n", zmin, zmin+(zstep-1)*step_z, step_z);
	printf("rmin = %le, rmax = %le, step_r = %le\n", rmin, rmin+(rstep-1)*step_r, step_r);
	printf("thmin = %le, thmax = %le, step_th = %le\n", theta_min, theta_min+(thstep-1)*step_th, step_th);
	wfile = fopen(textfile, "w");
	fprintf(wfile, "%i %i %i\n",rstep, thstep, zstep);
	fprintf(wfile, "r(m),theta(deg),z(m),Br(T),Bth,Bz\n");
	for(i=0;i<zstep;i++) {
		z = zmin + i*step_z;
		for(j=0;j<rstep;j++) {
			r = rmin + j*step_r;
			for(k=0;k<thstep;k++) {
				th = theta_min + k*step_th;
				map_neighbours2(&bx, &by, &bz, x0, y0, r0, th0, theta_tot, r, th, z, scale, cell, YES);
				bxs = bx*cos(th0) - by*sin(th0);
				bys = by*cos(th0) + bx*sin(th0);
				br = bxs*cos(th) + bys*sin(th);
				bth = bys*cos(th) - bxs*sin(th);
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r, th*180./PI, z, br, bth, bz);
			}
		}
	}
	fclose(wfile);
}

extern void vffa_map_maker2(char *textfile, double x0, double y0, double r0, double th0_deg, double scale, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, char *lattfile)
{
	struct Lattice latt;
	
	load_lattice(&latt, lattfile);
	printf("lattice loaded 2\n");
	compute_circmap_from_cartmap2(textfile, x0, y0, r0, th0_deg, scale, rmin, rmax, rstep, theta_min_deg, theta_max_deg, thstep, zmin, zmax, zstep, &(latt.cell[0]));
	free_latt(&latt);
}

extern void translate_shinji_vffa_file_wrapper(char *txt_prefix, double zshift, double yshift, double scale, double r0, double th0_deg, double theta_tot_deg)
{
	char shinji_file[500], outfile[500];
	
	sprintf(shinji_file,"%s.dat", txt_prefix);
	sprintf(outfile,"%s_ff_coord.dat", txt_prefix);
	translate_shinji_vffa_file(shinji_file, outfile, zshift, yshift, scale, r0, th0_deg, theta_tot_deg);
	
}

extern void convert_shinji_coord_to_ff(double *r, double *th, double *z, double xshinji, double yshinji, double zshinji, double r0, double th0, double yshift, double zshift)
{
	double y;
	
	y=yshinji-yshift-r0*tan(th0);
	*z=zshinji-zshift;
	*r = sqrt(y*y + (xshinji+r0)*(xshinji+r0));
	*th = atan_ratio(y, xshinji+r0) + th0;
}

//angle [rad]
// of only when angle with radius is 0!
extern void convert_ff_coord_to_shinji(double *xshinji, double *yshinji, double *zshinji, double r, double th, double zff, double r0, double th0, double yshift, double zshift)
{
	*xshinji = r*cos(th-th0) - r0;
	*yshinji = r*sin(th-th0) + yshift;
	*zshinji = zff + zshift;
}

//angle [rad] OK with angle !=0 and horizontal shift
extern void convert_ff_coord_to_shinji2(double *xshinji, double *yshinji, double *zshinji, double x, double y, double zff, double r0, double th0, double angle, double xshift, double yshift, double zshift)
{
	double x0,y0;
	x0 = r0*cos(th0);
	y0 = r0*sin(th0);
	*xshinji = (y-y0)*sin(th0+angle)+(x-x0)*cos(th0+angle) + xshift;
	*yshinji = (y-y0)*cos(th0+angle)-(x-x0)*sin(th0+angle) + yshift;
	*zshinji = zff + zshift;
}

extern void create_cartesian_mesh_tilted(char *textout, double yshift, double zshift, double r0, double th0_deg, double theta_tot_deg, double scale, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz, struct Cell *cell)
{
	int i,j,k;
	double xstep,ystep,zstep,x_tilted=0,y_tilted,z_tilted,bx,by,bz,r,th,x,y,z,bx_tilted,by_tilted,th0=th0_deg*PI/180., theta_tot=theta_tot_deg*PI/180.;
	FILE *wfile;
	
	xstep = comp_step(xmin,xmax,nbx);
	ystep = comp_step(ymin,ymax,nby);
	zstep = comp_step(zmin,zmax,nbz);
	
	wfile=fopen(textout,"w");
	fprintf(wfile, "%i	%i	%i\n", nbx,nby,nbz);
	fprintf(wfile, "x(m),y,z,Bx(T),By,Bz\n");
	for(i=0;i<nbz;i++) {
		z=zmin+i*zstep;
		for(j=0;j<nbx;j++) {
			x=xmin+j*xstep;
			for(k=0;k<nby;k++) {
				y=ymin+k*ystep;
				r = sqrt(x*x+y*y);
				th = atan_ratio(y,x);
				map_neighbours2(&bx_tilted, &by_tilted, &bz, 0., yshift, r0, th0, theta_tot, r, th, z, scale, cell, NO);
				//convert_ff_coord_to_shinji(&x_tilted, &y_tilted, &z_tilted, r, th, z, r0, th0, yshift, zshift);
				//printf("x,y,z=(%lf,%lf,%lf), r,th=(%lf,%lf), tilted=(%lf,%lf,%lf)\n",x,y,z,r,th*180./PI,x_tilted,y_tilted,z_tilted);
				//if(get_bfield(x_tilted, y_tilted, z_tilted, &bx_tilted, &by_tilted, &bz, cell, get_bfield_cartesian_map)!=ALIVE) {
				//	bx_tilted=0;by_tilted=0;bz=0;
					//errorstop("bfield!=ALIVE");
				//}
				//bx_tilted*=scale;by_tilted*=scale;bz*=scale;
				bx = bx_tilted*cos(th0) - by_tilted*sin(th0);
				by = by_tilted*cos(th0) + bx_tilted*sin(th0);
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n",x,y,z,bx,by,bz);
			}
		}
	}
	fclose(wfile);
	
}

extern void create_cartesian_mesh_tilted_wrapper(char *lattfile, char *textout, double yshift, double zshift, double r0, double th0_deg, double theta_tot, double scale, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz)
{
	struct Lattice latt;
	
	load_lattice(&latt, lattfile);
	printf("lattice loaded youhou\n");
	create_cartesian_mesh_tilted(textout, yshift, zshift, r0, th0_deg, theta_tot, scale, xmin, xmax, nbx, ymin, ymax, nby, zmin, zmax, nbz, &(latt.cell[0]));
	free_latt(&latt);
}

extern void create_cartesian_mesh_tilted_superimpose(char *textout, double yshift, double zshift, double r0, double th0_deg, double theta_tot, double xshift, double scale_minus, double scale_plus, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz, struct Cell *cell)
{
	int i,j,k;
	double xstep,ystep,zstep,x_tilted=0,y_tilted,z_tilted,bx,by,bz_1, bz_2,r,th,x,y,z,bx_tilted_1=0,by_tilted_1=0,bx_tilted_2=0,by_tilted_2=0,th0=th0_deg*PI/180.;
	FILE *wfile;
	
	xstep = comp_step(xmin,xmax,nbx);
	ystep = comp_step(ymin,ymax,nby);
	zstep = comp_step(zmin,zmax,nbz);
	
	wfile=fopen(textout,"w");
	fprintf(wfile, "%i	%i	%i\n", nbx,nby,nbz);
	fprintf(wfile, "x(m),y,z,Bx(T),By,Bz\n");
	for(i=0;i<nbz;i++) {
		z=zmin+i*zstep;
		for(j=0;j<nbx;j++) {
			x=xmin+j*xstep;
			for(k=0;k<nby;k++) {
				y=ymin+k*ystep;
				r = sqrt(x*x+y*y);
				th = atan_ratio(y,x);
				map_neighbours2(&bx_tilted_1, &by_tilted_1, &bz_1, -xshift, yshift, r0, th0, theta_tot, r, th, z, scale_minus, cell, NO);
				map_neighbours2(&bx_tilted_2, &by_tilted_2, &bz_2, xshift, yshift, r0, th0, theta_tot, r, th, z, scale_plus, cell, NO);
				bx_tilted_1+=bx_tilted_2;by_tilted_1+=by_tilted_2;bz_1+=bz_2;
				//convert_ff_coord_to_shinji(&x_tilted, &y_tilted, &z_tilted, r, th, z, r0, th0, yshift, zshift);
				//printf("x,y,z=(%lf,%lf,%lf), r,th=(%lf,%lf), tilted=(%lf,%lf,%lf)\n",x,y,z,r,th*180./PI,x_tilted,y_tilted,z_tilted);
				//if(get_bfield(x_tilted, y_tilted, z_tilted, &bx_tilted, &by_tilted, &bz, cell, get_bfield_cartesian_map)!=ALIVE) {
				//	bx_tilted=0;by_tilted=0;bz=0;
					//errorstop("bfield!=ALIVE");
				//}
				//bx_tilted*=scale;by_tilted*=scale;bz*=scale;
				bx = bx_tilted_1*cos(th0) - by_tilted_1*sin(th0);
				by = by_tilted_1*cos(th0) + bx_tilted_1*sin(th0);
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n",x,y,z,bx,by,bz_1);
			}
		}
	}
	fclose(wfile);
	
}

extern void create_cartesian_mesh_tilted_superimpose_wrapper(char *lattfile, char *textout, double yshift, double zshift, double r0, double th0_deg, double theta_tot, double xshift, double scale_minus, double scale_plus, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz)
{
	struct Lattice latt;
	
	load_lattice(&latt, lattfile);
	printf("lattice loaded youhou double coil\n");
	create_cartesian_mesh_tilted_superimpose(textout, yshift, zshift, r0, th0_deg, theta_tot, xshift, scale_minus, scale_plus, xmin, xmax, nbx, ymin, ymax, nby, zmin, zmax, nbz, &(latt.cell[0]));
	free_latt(&latt);
}

extern void compute_field_longit_line(char *textfile, double x0, double z0, double ymin, double ymax, int nby, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double y, bx, by, bz;
	FILE *wfile;
	double ystep = comp_step(ymin, ymax, nby);
	
	
	wfile = fopen(textfile, "w");
	
	for(i=0;i<nby;i++) {
		y=ymin+i*ystep;
		printf("y=%lf\n", y);
		if(get_bfield(x0, y, z0, &bx, &by, &bz, cell, add_contribution_comp) != ALIVE) errorstop("compute_field_longit_line: get_bfield != ALIVE\n");
		fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", x0,y,z0,bx,by,bz);
	}
	fclose(wfile);
}

extern void function_textfile(char *textfile, double a, double b, double c, double m, double ymin, double ymax, int nby)
{
	int i;
	double fit_vert, fit_long, y, ystep = comp_step(ymin, ymax, nby);
	FILE *wfile;
	wfile = fopen(textfile, "w");
	
	for(i=0;i<nby;i++) {
		y=ymin+i*ystep;
		//fit atan
		fit_vert = a/PI*(atan((b-y)/c)+atan((b+y)/c));
		fit_long = a/PI/m/c*(1/(1+((b+y)/c)*((b+y)/c))-1/(1+((b-y)/c)*((b-y)/c)));
		
		//fit tanh c must be divided by 2!
		//fit_vert = a/2.*(tanh((b-y)/c)+tanh((b+y)/c));
		//fit_long = a/2./c/m*(1./(cosh((b+y)/c)*cosh((b+y)/c))-1./(cosh((b-y)/c)*cosh((b-y)/c)));
		fprintf(wfile, "%le	%le	%le\n", y, fit_vert, fit_long);
	}
		
	fclose(wfile);
}

//x0,y0,z0: 0 of the map in latt framework
extern void create_straight_map_fodo(char *map_out, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz, char *latt_d_name, double x0d, double y0d, double z0d, double scaled, char *latt_f_name, double x0f, double y0f, double z0f, double scalef)
{
	int i,j,k;
	double xstep,ystep,zstep,x,y,z,bx,by,bz, bxf,byf,bzf,bxd,byd,bzd;
	struct Lattice latt_f, latt_d;
	FILE *wfile;
	doyou_long_boun_sym = NO;
	xstep = comp_step(xmin,xmax,nbx);
	ystep = comp_step(ymin,ymax,nby);
	zstep = comp_step(zmin,zmax,nbz);
	
	load_lattice(&latt_d, latt_d_name);
	load_lattice(&latt_f, latt_f_name);
	
	wfile=fopen(map_out,"w");
	fprintf(wfile, "%i	%i	%i\n", nbx,nby,nbz);
	fprintf(wfile, "x(m),y,z,Bx(T),By,Bz\n");
	for(i=0;i<nbz;i++) {
		z=zmin+i*zstep;
		for(j=0;j<nbx;j++) {
			x=xmin+j*xstep;
			for(k=0;k<nby;k++) {
				y=ymin+k*ystep;
				bxd=0;byd=0;bzd=0;bxf=0;byf=0;bzf=0;
				map_neighbours_straight(&bxd, &byd, &bzd, x0d, y0d, z0d, ymax-ymin, scaled, x, y, z, &(latt_d.cell[0]), NO);
				map_neighbours_straight(&bxf, &byf, &bzf, x0f, y0f, z0f, ymax-ymin, scalef, x, y, z, &(latt_f.cell[0]), NO);
				bx = bxd+bxf; by = byd+byf; bz = bzd+bzf;
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n",x,y,z,bx,by,bz);
			}
		}
	}
	fclose(wfile);
	free_latt(&latt_d);
	free_latt(&latt_f);
	doyou_long_boun_sym = YES;
}

extern void map_neighbours_straight(double *bx, double *by, double *bz, double x0, double y0, double z0, double ymax, double scale, double x, double y, double z, struct Cell *cell, int kill_notalive)
{
	double bx0=0,by0=0,bz0=0,bx1=0,by1=0,bz1=0,bx2=0,by2=0,bz2=0;
	
	if(get_bfield(x-x0, y-y0, z-z0, &bx0, &by0, &bz0, cell, get_bfield_cartesian_map) != ALIVE) {
		//printf("r=%le, th=%le, (x,y)=%le,%le, (x1,y1)=%le,%le\n", r, th*180./PI, x,y, x1,y1);
		if(kill_notalive==YES) errorstop("map_neighbours: get_bfield != ALIVE\n");
		else {
			bx0=0;by0=0;bz0=0;
		}
	}
	bx0*=scale;by0*=scale;bz0*=scale;
	*bx = bx0; *by = by0; *bz = bz0;
	//printf("b0=(%le,%le,%le)\n",bx0,by0,bz0);
	//printf("r=%le, th=%le, (x,y)=%le,%le, (x1,y1)=%le,%le", r, th*180./PI, x,y, x1,y1);
	if(get_bfield(x-x0, y-y0+ymax, z-z0, &bx1, &by1, &bz1, cell, get_bfield_cartesian_map) == ALIVE) { //cell before
		//printf("b1=(%le,%le,%le)\n",bx1,by1,bz1);
		bx1*=scale;by1*=scale;bz1*=scale;
		*bx += bx1; *by += by1; *bz += bz1;
		//printf("   YES");
	}
	//printf("b=(%le,%le,%le)\n", bx,by,bz);
	if(get_bfield(x-x0, y-y0-ymax, z-z0, &bx2, &by2, &bz2, cell, get_bfield_cartesian_map) == ALIVE) { //cell after
		//printf("b2=(%le,%le,%le)\n",bx2,by2,bz2);
		bx2*=scale;by2*=scale;bz2*=scale;
		*bx += bx2; *by += by2; *bz += bz2;
	}//*/
	//printf("tot (%le,%le,%le), %le + %le\n", *bx,*by,*bz,  bx2*cos(theta_tot), -by2*sin(theta_tot));
}

extern void wrapper_fodo_vffa_from_map(char *foldername, double xmin, double xmax, int nbx, double ymax, int nby, double zmin, double zmax, int nbz, char *lattname_ori, double x0d, double y0d, double z0d, double scaled, double x0f, double y0f, double z0f, double scalef)
{
	char mapfile[500], ori_latt_tot_name[500], lattname[500], heatfile[500], heat_bx[500], heat_by[500], heat_bz[500], heat_bx_eps[500], heat_by_eps[500], heat_bz_eps[500], xrange[100], yrange[100], beamfile[500], trackout[500], traj_xy[500], traj_sz[500], btot[500], tune_scan[500], tune_scan_eps[500];
	char symfile[500];
	int i,j,k;
	double min, max, z_heat=(zmin+zmax)/2.;
	double m4d[4][4], band_bias, t[4][4], r[4][4], output[6], brho0, brhostep, z_0;
	struct Lattice latt;
	struct Beam beam;
	
	sprintf(mapfile,"%s/map_fodo.dat", foldername);
	sprintf(ori_latt_tot_name,"%s/%s", foldername, lattname_ori);
	sprintf(lattname,"%s/fodo_latt.latt", foldername);
	sprintf(heatfile,"%s/fodo_heatmap_horplane", foldername);
	sprintf(heat_bx,"%s/fodo_heatmap_horplane_bx.dat", foldername);
	sprintf(heat_bx_eps,"%s/fodo_heatmap_horplane_bx.eps", foldername);
	sprintf(heat_by,"%s/fodo_heatmap_horplane_by.dat", foldername);
	sprintf(heat_by_eps,"%s/fodo_heatmap_horplane_by.eps", foldername);
	sprintf(heat_bz,"%s/fodo_heatmap_horplane_bz.dat", foldername);
	sprintf(heat_bz_eps,"%s/fodo_heatmap_horplane_bz.eps", foldername);
	sprintf(beamfile,"%s/protontest.beam", foldername);
	sprintf(trackout,"%s/fodo_trackout.dat", foldername);
	sprintf(traj_xy,"%s/fodo_traj_xy.eps", foldername);
	sprintf(traj_sz,"%s/fodo_traj_sz.eps", foldername);
	sprintf(btot,"%s/fodo_btot.eps", foldername);
	sprintf(tune_scan,"%s/tune_scan.dat", foldername);
	sprintf(tune_scan_eps,"%s/tune_scan.eps", foldername);
	sprintf(symfile,"%s/map_symmetry_check.dat", foldername);
	
	create_straight_map_fodo(mapfile, xmin, xmax, nbx, 0, ymax, nby, zmin, zmax, nbz, ori_latt_tot_name, x0d, y0d, z0d, scaled, ori_latt_tot_name, x0f, y0f, z0f, scalef);
	write_maplatt_file(lattname, "straight vffa fodo from map", 1, mapfile, 0.001, ymax, xmin, xmax, nbx, 0, ymax, nby, zmin, zmax, nbz);
	//write_beam_file(beamfile, "proton test", 938.27, 100.e6, 1, 0, z_heat, xp_start_deg, zp_start_deg);
	
	
	load_lattice(&latt,lattname);
	
	compare_symetry_field_straight(symfile, xmin, xmax, nbx, zmin, zmax, nbz, &(latt.cell[0]));
	
	
	compute_fieldheatmap(heatfile, &(latt.cell[0]), get_bfield_cartesian_map, xmin, xmax, nbx, 0, ymax, nby, z_heat, z_heat, 1, NO);
	find_boundaries_heated_map_file2(heat_bx, xrange, yrange, &min, &max);
	easyplot3dmap(heat_bx, "longitudinal [m]", "horizontal [m]", "B horizontal [T]", xrange, yrange, heat_bx_eps, "size ratio -1\nset mxtics 5\nset mytics 5\nset mxtics 5", -2, 2);//-5, 5?
	find_boundaries_heated_map_file2(heat_by, xrange, yrange, &min, &max);
	easyplot3dmap(heat_by, "longitudinal [m]", "horizontal [m]", "B longitudinal [T]", xrange, yrange, heat_by_eps, "size ratio -1\nset mxtics 5\nset mytics 5\nset mxtics 5", min, max);
	find_boundaries_heated_map_file2(heat_bz, xrange, yrange, &min, &max);
	easyplot3dmap(heat_bz, "longitudinal [m]", "horizontal [m]", "B vertical [T]", xrange, yrange, heat_bz_eps, "size ratio -1\nset mxtics 5\nset mytics 5\nset mxtics 5", min, max);
	load_beam(&beam, beamfile, &latt, YES);
	//emptyfile(trackout);
	//part_cross_latt(&(beam.part[0]), &latt, trackout);
	//get_xrange(trackout, xrange);
	////plot_traj(&latt, "inputs/vffa/map/1_25_0/fodo_trackout.dat", "inputs/vffa/map/1_25_0/fodo_aspect.dat", "($3)", "($2)", "($2)", "($1)","lines lc 7 lw 2", "lines lt 1 lc 0 lw 2", "y [m]", "x [m]", NULL, NULL, "inputs/vffa/map/1_25_0/fodo_traj_xy.eps", "size ratio -1\nset mxtics 5\nset mytics 5");
	//easyplot3dmap_plusfile(trackout, "($3)", "($2)", "lines lc 1 lw 2", heat_bx, "y [m]", "x [m]", "B hor [T]", NULL, NULL, traj_xy, "grid\nset size ratio -1\nset mxtics 5\nset mytics 5\nset mxtics 5", -2, 2);
	//easyplot(trackout, "($1)", "($4)", "lines lc 7 lw 2", "s [m]", "z [m]", xrange, NULL, traj_sz, "mxtics 5\nset mytics 5");
	//easyplot3p(trackout, trackout, trackout, "1", "5", "1", "6", "1", "7", "lines lw 3 lc 2", "lines lw 3 lc 3", "lines lw 3 lc 7", 
	//"horizontal", "longitudinal", "vertical", "s [m]", "B [T]", xrange, NULL, btot, "grid\nset mxtics 5\nset mytics 5");
	print_part_para("beam.part", &(beam.part[0]));
	emptyfile(tune_scan);
	brho0=beam.part[0].brho;
	z_0 = beam.part[0].z;
	brhostep = 0.74/30.;
	//beam.part[0].brho=1.189019;
	//beam.part[0].z=-0.123376;
	//if(put_on_co_nelmin(&latt,&(beam.part[0]), 1.e-20, 1.e-3, 1.e-3, 1.e-3, 1.e-3, YES)!=TRUE) errorstop("closed orbit not found");
	
	for(j=0;j<30;j++) {
		beam.part[0].brho=brho0-j*brhostep;
		if(j!=0) beam.part[0].z=log(beam.part[0].brho/brho0)/1.3+z_0;
		printf("j=%i, brho=%lf\t",j, beam.part[0].brho);
		tune_for_tunescan(tune_scan, &(beam.part[0]), &latt, 1.e-9, 1.e-5);
	}
	
	//phase_adv_scan(tune_scan, &latt, &(beam.part[0]), 1.e-8, 20.0e6, 90.e6, 35, 1.e-5);	//tilt test tune scan
	easyplot2(tune_scan, tune_scan, "($1*1.e-6)", "2", "($1*1.e-6)", "3", "lines lt 1 lw 4 lc 0", "lines lt 2 lw 4 lc 0", NULL, NULL, "Momentum [MeV/c]", "{/Symbol n}_u_,_v []", NULL, NULL, tune_scan_eps, NULL);
	//calc_tune_twiss(&(beam.part[0]), &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-4, 1.e-4, 1.e-4, 1.e-4, &latt, part_cross_latt, YES, NULL, 1);//flag: 0=2d, 1=4d decoupled, anything else=4d undecoupled.
	
	free_latt(&latt);
	free_beam(&beam);
}//*/

extern void write_beam_file(char *beamfile, char *title, double e0_mev, double ekin_ev, int nbkeywords, double x0, double z0, double xprime_deg, double zprime_deg)
{
	FILE *wfile;
	
	wfile = fopen(beamfile, "w");
	fprintf(wfile, "%s\n\n", title);
	fprintf(wfile, "%lf	1		//object rest mass [MeV.c^-2], and charge [q.C]\n\n", e0_mev);
	fprintf(wfile, "%i			//number of keywords to be read\n\n", nbkeywords);
	fprintf(wfile, "singlepart %le	// Ekin [eV]\n", ekin_ev);
	fprintf(wfile, "           %lf 0 %lf	// x [m], y [m], z [m] (x:horiz, y:longit, z:verti)\n", x0,z0);
	fprintf(wfile, "           %le %le		// atan_dxdy_deg [deg.], atan_dzdy_deg [deg.]\n\n", xprime_deg, zprime_deg);
	fprintf(wfile, "putonco-xxp 1.e-12		// absolute precision is given in [m]\n");
	//fprintf(wfile, "putonco-xxpzzp 1.e-10		// absolute precision is given in [m]\n");
	//fprintf(wfile, "putonco-nelmin 1.e-20	1.e-3	// absolute precision is given in [m]\n\n");
	fclose(wfile);
}

extern void scattered_point_file_vffa_alignement(char *txtnoerr, char *txtfile, char *txtoutu, char *txtoutv, char *txtoutcod, char *txtouttuneu, char *txtouttunev, double trans_shift, double rot_shift, double scale, double ustep, double vstep, double codstep, double tunestep)
{
	long seed;
	int n, nblines, k,i=0, nb_uarray[1000], j=0, nb_varray[1000], l=0, nb_codarray[1000], m1=0, nb_tuneuarray[1000], m2=0, nb_tunevarray[1000], flag_found;
	double x, ux, z, uz, qx, qz, beta_x, alpha_x, beta_z, alpha_z, aumax, avmax, cod, eps_umax, eps_vmax, uarray[1000], varray[1000], codarray[1000], tuneuarray[1000], tunevarray[1000], tuneu, tunev;
	double x_noerr, ux_noerr, z_noerr, uz_noerr, qx_noerr, qz_noerr, beta_x_noerr, alpha_x_noerr, beta_z_noerr, alpha_z_noerr, aumax_noerr, avmax_noerr, da_offset;
	FILE *rfile_noerr=NULL;
	FILE *rfile_err=NULL;
	FILE *wfileu;
	FILE *wfilev;
	FILE *wfilecod;
	FILE *wfiletuneu;
	FILE *wfiletunev;
	
	for(n=0;n<1000;n++) {
		nb_uarray[n]=0;
		nb_varray[n]=0;
		nb_codarray[n]=0;
		nb_tuneuarray[n]=0;
		uarray[n]=0;
		varray[n]=0;
		codarray[n]=0;
		tuneuarray[n]=0;
	}
	
	rfile_noerr = fopen(txtnoerr, "r");
	if(rfile_noerr==NULL) errorstop("cannot open no error file");
	fscanf(rfile_noerr, "%ld	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le", &seed, &x_noerr, &ux_noerr, &z_noerr, &uz_noerr, &qx_noerr, &qz_noerr, &beta_x_noerr, &alpha_x_noerr, &beta_z_noerr, &alpha_z_noerr, &aumax_noerr, &avmax_noerr, &da_offset);
	fclose(rfile_noerr);
	
	nblines = get_nb_lines_file(txtfile);
		
	rfile_err = fopen(txtfile, "r");
	if(rfile_err==NULL) errorstop("cannot open error file");
	for(n=0;n<nblines;n++) {
		fscanf(rfile_err, "%ld	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le", &seed, &x, &ux, &z, &uz, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, &aumax, &avmax, &da_offset);
		eps_umax = aumax*aumax/beta_x_noerr*scale;
		flag_found=NO;
		for(k=0;k<i;k++) {
			if(uarray[k]-ustep<eps_umax && eps_umax<uarray[k]+ustep) {
				nb_uarray[k] ++;
				flag_found=YES;
				break;
			}
		}
		if(flag_found==NO) {
			uarray[i]=eps_umax;
			nb_uarray[i]=1;
			i++;
		}
		
		eps_vmax = avmax*avmax/beta_z_noerr*scale;
		flag_found=NO;
		for(k=0;k<j;k++) {
			if(varray[k]-vstep<eps_vmax && eps_vmax<varray[k]+vstep) {
				nb_varray[k] ++;
				flag_found=YES;
				break;
			}
		}
		if(flag_found==NO) {
			varray[j]=eps_vmax;
			nb_varray[j]=1;
			j++;
		}
		
		cod = sqrt(pow(x-x_noerr,2)+pow(z-z_noerr,2));
		flag_found=NO;
		for(k=0;k<l;k++) {
			if(codarray[k]-codstep<cod && cod<codarray[k]+codstep) {
				nb_codarray[k] ++;
				flag_found=YES;
				break;
			}
		}
		if(flag_found==NO) {
			codarray[l]=cod;
			nb_codarray[l]=1;
			l++;
		}
		flag_found=NO;
		for(k=0;k<l;k++) {
			if(codarray[k]-codstep<cod && cod<codarray[k]+codstep) {
				nb_codarray[k] ++;
				flag_found=YES;
				break;
			}
		}
		if(flag_found==NO) {
			codarray[l]=cod;
			nb_codarray[l]=1;
			l++;
		}
		
		tuneu = qx-qx_noerr;
		flag_found=NO;
		for(k=0;k<m1;k++) {
			if(tuneuarray[k]-tunestep<tuneu && tuneu<tuneuarray[k]+tunestep) {
				nb_tuneuarray[k] ++;
				flag_found=YES;
				break;
			}
		}
		if(flag_found==NO) {
			tuneuarray[m1]=tuneu;
			nb_tuneuarray[m1]=1;
			m1++;
		}
		tunev = qz-qz_noerr;
		flag_found=NO;
		for(k=0;k<m2;k++) {
			if(tunevarray[k]-tunestep<tunev && tunev<tunevarray[k]+tunestep) {
				nb_tunevarray[k] ++;
				flag_found=YES;
				break;
			}
		}
		if(flag_found==NO) {
			tunevarray[m2]=tunev;
			nb_tunevarray[m2]=1;
			m2++;
		}
	}
	fclose(rfile_err);
	
	wfileu = fopen(txtoutu, "a");
	for(n=0;n<i;n++) {
		fprintf(wfileu, "%le	%le	%le	%i\n", trans_shift, rot_shift, uarray[n], nb_uarray[n]);
	}
	fclose(wfileu);
	
	wfilev = fopen(txtoutv, "a");
	for(n=0;n<j;n++) {
		fprintf(wfilev, "%le	%le	%le	%i\n", trans_shift, rot_shift, varray[n], nb_varray[n]);
	}
	fclose(wfilev);
	
	wfilecod = fopen(txtoutcod, "a");
	for(n=0;n<l;n++) {
		fprintf(wfilecod, "%le	%le	%le	%i\n", trans_shift, rot_shift, codarray[n], nb_codarray[n]);
	}
	fclose(wfilecod);
	wfiletuneu = fopen(txtouttuneu, "a");
	for(n=0;n<m1;n++) {
		fprintf(wfiletuneu, "%le	%le	%le	%i\n", trans_shift, rot_shift, tuneuarray[n], nb_tuneuarray[n]);
	}
	fclose(wfiletuneu);
	wfiletunev = fopen(txtouttunev, "a");
	for(n=0;n<m2;n++) {
		fprintf(wfiletunev, "%le	%le	%le	%i\n", trans_shift, rot_shift, tunevarray[n], nb_tunevarray[n]);
	}
	fclose(wfiletunev);
}

extern void tune_for_tunescan(char *filename, struct Particle *part, struct Lattice *latt, double eps_clo, double amp)
{
	double m4d[4][4], band_bias, t[4][4], r[4][4], output[6], det;
	struct Particle test_part;
	FILE *wfile;
	test_part= *part;
	test_part.hat = -2;
	
	wfile = fopen(filename, "a");
	if(put_on_co_nelmin(latt, &test_part, eps_clo, 1.e-3, 1.e-3, 1.e-3, 1.e-3, YES)!=TRUE) errorstop("closed orbit not found");
	if(find_closed_orbite_xxp_zzp(&(test_part), &(test_part.x), &(test_part.ux), &(test_part.uy), &(test_part.z), &(test_part.uz), eps_clo, latt, YES)!=TRUE) errorstop("closed orbit linear maethod not found\n");
	//print_part_para("test_part", &test_part);
	get_matrix_firstorder_4d(m4d, &band_bias, &(test_part), amp, amp, amp, amp, latt, part_cross_latt, NO);
	det = matrix_det_4d(m4d, 4);
	decouple_matrix_parzen(m4d, t, r, output, NULL);
	fprintf(wfile, "%le	%le	%le	%le	%le\n", test_part.brho*test_part.q*CLIGHT/(UNITCHARGE), output[0], output[1], det, amp);
	fclose(wfile);
}

extern void adjust_triplet_with_two_magnets(char *lattfile_triplet, char *lattfile_twomag)
{
	int i;
	double x,y,z, bx_trip, by_trip, bz_trip, bx_2, by_2, bz_2;
	struct Lattice latt_triplet, latt_twomag;
	
	load_lattice(&latt_triplet,lattfile_triplet);
	set_vffa_rect_convergence_limit(&(latt_triplet.cell[0]), NO, 1.5, add_field_comp_VFFA_rect_atan);
	
	load_lattice(&latt_twomag,lattfile_twomag);
	set_vffa_rect_convergence_limit(&(latt_twomag.cell[0]), NO, 1.5, add_field_comp_VFFA_rect_atan);
	
	x= 0; y=2.5; z=-0.5;
	get_bfield(x,y,z, &bx_trip, &by_trip, &bz_trip, &(latt_triplet.cell[0]),add_field_comp_VFFA_rect_str_atan);
	printf("bz_trip=%le\n", bz_trip);
	for(i=0;i<100;i++) {
		get_bfield(x,y,z, &bx_2, &by_2, &bz_2, &(latt_twomag.cell[0]),add_field_comp_VFFA_rect_str_atan);
		printf("bz_2=%le, b0=%le\n", bz_2, latt_twomag.cell[0].mpara[1][2]);
		if(fabs(bz_2-bz_trip)<1.e-8) {
			printf("found! b0=%le\n", latt_twomag.cell[0].mpara[1][2]);
			break;
		}
		latt_twomag.cell[0].mpara[1][2]+=(bz_trip-bz_2)/bz_trip;
	}
	
	free_latt(&latt_triplet);
	free_latt(&latt_twomag);
}

/*extern void gene_coupled_ellipse_beam(char *textfileu, char *textfile, double emit_u, double twiss_betau, double twiss_alphau, double emit_v, double twiss_betav, double twiss_alphav, double r[4][4], int nparts)
{
	char filename[500], filename2[500], filename3[500], filename4[500], emitu_uuprime[300], emitu_vvprime[300], emitu_xxprime[300], emitu_zzprime[300], emitv_uuprime[300], emitv_vvprime[300], emitv_xxprime[300], emitv_zzprime[300];
	int i;
	double uv_vect[4], xz_vect[4], t, un, upn, vn, vpn,x,xp,z,zp;
	double x_ctr, xprime_ctr, z_ctr, zprime_ctr;
	FILE *wfileu;
	FILE *wfilev;
	
	x_ctr = 0;
	xprime_ctr = 0;
	z_ctr = 0;
	zprime_ctr = 0;
	
	wfileu = fopen(textfileu, "w");
	
	//ellipse in (u, u')
	for(i = 0; i <= nparts; i++) {
		t = i*2.*PI/nparts;
		un = cos(t);
		upn = sin(t);
	
		uv_vect[0] = sqrt(emit_u*twiss_betau)*un;
		uv_vect[1] = sqrt(emit_u/twiss_betau)*(upn - twiss_alphau*un);
		uv_vect[2] = 0;
		uv_vect[3] = 0;
		mvprod4(xz_vect, r, uv_vect);
		x = x_ctr + xz_vect[0];
		xp = xprime_ctr + xz_vect[1];
		z = z_ctr + xz_vect[2];
		zp = zprime_ctr + xz_vect[3];
		fprintf(wfileu, "%le	%le	%le	%le	%le	%le	%le	%le\n", uv_vect[0],uv_vect[1],uv_vect[2],uv_vect[3], x,xp,z,zp);
	}
	fclose(wfileu);
	
	sprintf(filename3,"%s/matched_beam_emitv_xzspace.dat", emplacement);
	sprintf(filename4,"%s/matched_beam_emitv_uvspace.dat", emplacement);
	sprintf(emitv_uuprime,"%s/matched_beam_emitv_uuprime.eps", emplacement);
	sprintf(emitv_vvprime,"%s/matched_beam_emitv_vvprime.eps", emplacement);
	sprintf(emitv_xxprime,"%s/matched_beam_emitv_xxprime.eps", emplacement);
	sprintf(emitv_zzprime,"%s/matched_beam_emitv_zzprime.eps", emplacement);
	wfilev = fopen(textfilev, "w");
	//ellipse in (v, v')
	for(i = 0; i <= nparts; i++) {
		t = i*2.*PI/nparts;
		vn = cos(t);
		vpn = sin(t);
	
		uv_vect[0] = 0;
		uv_vect[1] = 0;
		uv_vect[2] = sqrt(emit_v*twiss_betav)*vn;
		uv_vect[3] = sqrt(emit_v/twiss_betav)*(vpn - twiss_alphav*vn);
		mvprod4(xz_vect, r, uv_vect);
		x = x_ctr + xz_vect[0];
		xp = xprime_ctr + xz_vect[1];
		z = z_ctr + xz_vect[2];
		zp = zprime_ctr + xz_vect[3];
		fprintf(wfilev, "%le	%le	%le	%le	%le	%le	%le	%le\n", uv_vect[0],uv_vect[1],uv_vect[2],uv_vect[3], x,xp,z,zp);
	}
	fclose(wfilev);
	
	printf("\n");
	//easyplot(filename2, "($1*1000)", "($2*1000)", "points pt 7 ps 0.2","u [mm]", "u\251 [mrad]", NULL, NULL, emitu_uuprime, "mxtics 5\nset mytics 5");
	//easyplot(filename2, "($3*1000)", "($4*1000)", "points pt 7 ps 0.2","v [mm]", "v\251 [mrad]", NULL, NULL, emitu_vvprime, "mxtics 5\nset mytics 5");
	//easyplot(filename, "($1*1000)", "($2*1000)", "points pt 7 ps 0.2","x [mm]", "x\251 [mrad]", NULL, NULL, emitu_xxprime, "mxtics 5\nset mytics 5");
	//easyplot(filename, "($3*1000)", "($4*1000)", "points pt 7 ps 0.2","z [mm]", "z\251 [mrad]", NULL, NULL, emitu_zzprime, "mxtics 5\nset mytics 5");
	//easyplot(filename4, "($1*1000)", "($2*1000)", "points pt 7 ps 0.2","u [mm]", "u\251 [mrad]", NULL, NULL, emitv_uuprime, "mxtics 5\nset mytics 5");
	//easyplot(filename4, "($3*1000)", "($4*1000)", "points pt 7 ps 0.2","v [mm]", "v\251 [mrad]", NULL, NULL, emitv_vvprime, "mxtics 5\nset mytics 5");
	//easyplot(filename3, "($1*1000)", "($2*1000)", "points pt 7 ps 0.2","x [mm]", "x\251 [mrad]", NULL, NULL, emitv_xxprime, "mxtics 5\nset mytics 5");
	//easyplot(filename3, "($3*1000)", "($4*1000)", "points pt 7 ps 0.2","z [mm]", "z\251 [mrad]", NULL, NULL, emitv_zzprime, "mxtics 5\nset mytics 5");
}//*/

extern void easyplot_3dtraj(char *trackout1, char *trackout2, char *trackout3, struct Lattice *latt, double xmin, double zmin, char *title1, char *title2, char *title3, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption)
{
	int i, lc_tck1=7, lc_tck2=7, lc_tck3=7;
	FILE *gp;
	
	for(i=0;i<latt->nbcell;i++) aspect_rect_vffa3d("data/aspect3d.dat", &(latt->cell[i]), NO);
	
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	fprintf(gp,"set bmargin 25\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xtics offset -0.5,-1\n");
	fprintf(gp," set ytics offset 2,-0.5\n");
	fprintf(gp," set xlabel '%s' offset 0,-2 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset 6,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set zlabel '%s' offset -6,0 font 'helvetica,25.'\n", zlabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(zrange != NULL) fprintf(gp," set zrange %s \n", zrange);
	if(title1==NULL && title2==NULL && title3==NULL) fprintf(gp,"unset key\n");
	
	//fprintf(gp, "set dgrid3d 50,50\n");
	//fprintf(gp, "set hidden3d\n");
	//fprintf(gp, "set view 60, 60\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set terminal postscript eps enhanced color \n");
	fprintf(gp," set output '%s' \n", psfilename);
	if(title1==NULL) fprintf(gp,"splot '%s' using 2:3:4 with lines lc %i lw 3 notitle, ", trackout1, lc_tck1);
	else fprintf(gp,"splot '%s' using 2:3:4 with lines lc %i lw 3 title '%s', ", trackout1, lc_tck1, title1);
	if(title2==NULL) fprintf(gp,"'%s' using 2:3:4 with lines lc %i lw 3 notitle, ", trackout2, lc_tck2);
	else fprintf(gp,"'%s' using 2:3:4 with lines lc %i lw 3 title '%s', ", trackout2, lc_tck2, title2);
	if(title3==NULL) fprintf(gp,"'%s' using 2:3:4 with lines lc %i lw 3 notitle, ", trackout3, lc_tck3);
	else fprintf(gp,"'%s' using 2:3:4 with lines lc %i lw 3 title '%s', ", trackout3, lc_tck3, title3);
	fprintf(gp,"'data/aspect3d.dat' using 1:2:3 with lines lt 1 lc 0 lw 1.5 notitle, ");
	fprintf(gp,"'data/aspect3d.dat' using (%lf):2:3 with lines lt 1 lc 0 lw 0.5 dt 2 notitle, ", xmin);
	fprintf(gp,"'data/aspect3d.dat' using 1:2:(%lf) with lines lt 1 lc 0 lw 0.5 dt 2 notitle, ", zmin);
	
	fprintf(gp,"'%s' using (%lf):3:4 with lines lc %i lw 1 dt 2  notitle, ", trackout1, xmin, lc_tck1);
	fprintf(gp,"'%s' using 2:3:(%lf) with lines lc %i lw 1 dt 2  notitle, ", trackout1, zmin, lc_tck1);
	
	fprintf(gp,"'%s' using (%lf):3:4 with lines lc %i lw 1 dt 2  notitle, ", trackout2, xmin, lc_tck2);
	fprintf(gp,"'%s' using 2:3:(%lf) with lines lc %i lw 1 dt 2  notitle, ", trackout2, zmin, lc_tck2);
	
	fprintf(gp,"'%s' using (%lf):3:4 with lines lc %i lw 1 dt 2  notitle, ", trackout3, xmin, lc_tck3);
	fprintf(gp,"'%s' using 2:3:(%lf) with lines lc %i lw 1 dt 2  notitle\n", trackout3, zmin, lc_tck3);
	
	//fprintf(gp," set output\n");
	//fprintf(gp," set terminal x11\n");
	//fprintf(gp,"splot '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s\n", 
	//txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2, txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3, txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4, txtfilename5, xcolumn5, ycolumn5, zcolumn5, with5,  txtfilename6, xcolumn6, ycolumn6, zcolumn6, with6);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot_3dtraj_5mom(char *trackout1, char *trackout2, char *trackout3, char *trackout4, char *trackout5, struct Lattice *latt, double xmin, double zmin, 
char *title1, char *title2, char *title3, char *title4, char *title5, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption)
{
	int i, lc_tck1=7, lc_tck2=7, lc_tck3=7, lc_tck4=7, lc_tck5=7;
	FILE *gp;
	
	emptyfile("data/aspect3d.dat");
	for(i=0;i<latt->nbcell;i++) aspect_rect_vffa3d("data/aspect3d.dat", &(latt->cell[i]), NO);
	
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	fprintf(gp,"set bmargin 25\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xtics offset -0.5,-1\n");
	fprintf(gp," set ytics offset 2,-0.5\n");
	fprintf(gp," set xlabel '%s' offset 0,-2 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset 6,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set zlabel '%s' offset -6,0 font 'helvetica,25.'\n", zlabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(zrange != NULL) fprintf(gp," set zrange %s \n", zrange);
	if(title1==NULL && title2==NULL && title3==NULL) fprintf(gp,"unset key\n");
	
	//fprintf(gp, "set dgrid3d 50,50\n");
	//fprintf(gp, "set hidden3d\n");
	//fprintf(gp, "set view 60, 60\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set terminal postscript eps enhanced color \n");
	fprintf(gp," set output '%s' \n", psfilename);
	if(title1==NULL) fprintf(gp,"splot '%s' using 2:3:4 with lines lc %i lw 3 notitle, ", trackout1, lc_tck1);
	else fprintf(gp,"splot '%s' using 2:3:4 with lines lc %i lw 3 title '%s', ", trackout1, lc_tck1, title1);
	if(title2==NULL) fprintf(gp,"'%s' using 2:3:4 with lines lc %i lw 3 notitle, ", trackout2, lc_tck2);
	else fprintf(gp,"'%s' using 2:3:4 with lines lc %i lw 3 title '%s', ", trackout2, lc_tck2, title2);
	if(title3==NULL) fprintf(gp,"'%s' using 2:3:4 with lines lc %i lw 3 notitle, ", trackout3, lc_tck3);
	else fprintf(gp,"'%s' using 2:3:4 with lines lc %i lw 3 title '%s', ", trackout3, lc_tck3, title3);
	fprintf(gp,"'data/aspect3d.dat' using 1:2:3 with lines lt 1 lc 0 lw 1.5 notitle, ");
	fprintf(gp,"'data/aspect3d.dat' using (%lf):2:3 with lines lt 1 lc 0 lw 0.5 dt 2 notitle, ", xmin);
	fprintf(gp,"'data/aspect3d.dat' using 1:2:(%lf) with lines lt 1 lc 0 lw 0.5 dt 2 notitle, ", zmin);
	
	fprintf(gp,"'%s' using (%lf):3:4 with lines lc %i lw 1 dt 2  notitle, ", trackout1, xmin, lc_tck1);
	fprintf(gp,"'%s' using 2:3:(%lf) with lines lc %i lw 1 dt 2  notitle, ", trackout1, zmin, lc_tck1);
	
	fprintf(gp,"'%s' using (%lf):3:4 with lines lc %i lw 1 dt 2  notitle, ", trackout2, xmin, lc_tck2);
	fprintf(gp,"'%s' using 2:3:(%lf) with lines lc %i lw 1 dt 2  notitle, ", trackout2, zmin, lc_tck2);
	
	fprintf(gp,"'%s' using (%lf):3:4 with lines lc %i lw 1 dt 2  notitle, ", trackout3, xmin, lc_tck3);
	fprintf(gp,"'%s' using 2:3:(%lf) with lines lc %i lw 1 dt 2  notitle\n", trackout3, zmin, lc_tck3);
	
	fprintf(gp,"'%s' using (%lf):3:4 with lines lc %i lw 1 dt 2  notitle, ", trackout4, xmin, lc_tck4);
	fprintf(gp,"'%s' using 2:3:(%lf) with lines lc %i lw 1 dt 2  notitle\n", trackout4, zmin, lc_tck4);
	
	fprintf(gp,"'%s' using (%lf):3:4 with lines lc %i lw 1 dt 2  notitle, ", trackout5, xmin, lc_tck5);
	fprintf(gp,"'%s' using 2:3:(%lf) with lines lc %i lw 1 dt 2  notitle\n", trackout5, zmin, lc_tck5);
	
	//fprintf(gp," set output\n");
	//fprintf(gp," set terminal x11\n");
	//fprintf(gp,"splot '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s\n", 
	//txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2, txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3, txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4, txtfilename5, xcolumn5, ycolumn5, zcolumn5, with5,  txtfilename6, xcolumn6, ycolumn6, zcolumn6, with6);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}
//*/
extern void filter_fieldmap_Bf(char *text_fieldmap_Bf, char *outfile, double x0, double z0)
{
	int i, nblines;
	double x,y,z,bx,by,bz;
	FILE *wfile;
	FILE *rfile=NULL;
	
	nblines = get_nb_lines_file(text_fieldmap_Bf);

	rfile = fopen(text_fieldmap_Bf, "r");
	if(rfile==NULL) errorstop("can't open readfile");
	newline(rfile);
	newline(rfile);
	wfile=fopen(outfile, "w");
	
	for(i=0;i<nblines-2;i++) {
		fscanf(rfile, "%le	%le	%le	%le	%le	%le", &x, &z, &y, &bx, &bz, &by);
		if(fabs(x-x0)<TINYLENGTH && fabs(z-z0)<TINYLENGTH) fprintf(wfile, "%le	%le	%le	%.16e	%.16e	%.16e\n", x, z, y, bx, bz, by); 
	}
	fclose(wfile);
	fclose(rfile);
	
}

extern void create_circmap_triplet_from_cartmap(char *textfile, char *latt_f_name, double r0f, double th0f1_deg, double anglef1_deg, double scalef, char *latt_d_name, double r0d, double th0d_deg, double angled_deg, double scaled, double rmin, double rmax, int rstep, double theta_max_deg, int thstep, double zmin, double zmax, int zstep)
{
	int i,j,k;
	double theta_max, th0d,angled,th0f1,anglef1,th0f2,anglef2,x,y,z,r,th,xd,yd,r1,th1,bxd,byd,bzd,bx,by,bz,xf,yf,bxf1,byf1,bzf1,bxf2,byf2,bzf2,br,bth;
	double step_r,step_th,step_z;
	struct Lattice latt_f, latt_d;
	FILE *wfile;
	
	//magnet: r0, th0, angle, x_map=0, y_map=latt.cell[0].boun.ymax /2.;
	theta_max=theta_max_deg*PI/180.;
	th0d = th0d_deg*PI/180.;
	angled = angled_deg*PI/180.;
	th0f1 = th0f1_deg*PI/180.;
	anglef1 = anglef1_deg*PI/180.;
	th0f2 = theta_max - th0f1;
	anglef2 = -anglef1;
	step_r = comp_step(rmin,rmax,rstep);
	step_th = comp_step(0,theta_max,thstep);
	step_z = comp_step(zmin,zmax,zstep);
	
	printf("D MAGNET!\n");
	load_lattice(&latt_d, latt_d_name);
	printf("F MAGNET!\n");
	load_lattice(&latt_f, latt_f_name);
	
	wfile=fopen(textfile,"w");
	fprintf(wfile, "%i	%i	%i\n", rstep,thstep,zstep);
	fprintf(wfile, "r(m),th(rad),z,Bx(T),By,Bz\n");
	
	for(i=0;i<zstep;i++) {
		z=zmin+i*step_z;
		for(j=0;j<rstep;j++) {
			r = rmin + j*step_r;
			for(k=0;k<thstep;k++) {
				th = 0 + k*step_th;
				x=r*cos(th);
				y=r*sin(th);
				//coordinate_global_to_rect_vffa(x_g,y_g,*x_mag,*y_mag,r_mag_cent, th_mag_cent, tilt_mag, x_1,y_1)
				coordinate_global_to_rect_vffa(x, y, &xd, &yd, r0d, th0d, angled, 0, latt_d.cell[0].boun.ymax/2.); //careful not tested with shiftd!=0
				r1 = sqrt(xd*xd+yd*yd);
				th1 = atan_ratio(yd,xd);
				map_neighbours3(&bxd, &byd, &bzd, theta_max, r1, th1, z, scaled, &(latt_d.cell[0]), NO);
				bx = bxd*cos(th0d+angled) - byd*sin(th0d+angled);
				by = bxd*sin(th0d+angled) + byd*cos(th0d+angled);
				bz = bzd;
				coordinate_global_to_rect_vffa(x, y, &xf, &yf, r0f, th0f1, anglef1, 0, latt_f.cell[0].boun.ymax/2.); //careful not tested with shiftf!=0
				r1 = sqrt(xf*xf+yf*yf);
				th1 = atan_ratio(yf,xf);
				map_neighbours3(&bxf1, &byf1, &bzf1, theta_max, r1, th1, z, scalef, &(latt_f.cell[0]), NO);
				bx += bxf1*cos(th0f1+anglef1) - byf1*sin(th0f1+anglef1);
				by += bxf1*sin(th0f1+anglef1) + byf1*cos(th0f1+anglef1);
				bz += bzf1;
				coordinate_global_to_rect_vffa(x, y, &xf, &yf, r0f, th0f2, anglef2, 0, latt_f.cell[0].boun.ymax/2.); //careful not tested with shiftf!=0
				r1 = sqrt(xf*xf+yf*yf);
				th1 = atan_ratio(yf,xf);
				map_neighbours3(&bxf2, &byf2, &bzf2, theta_max, r1, th1, z, scalef, &(latt_f.cell[0]), NO);
				bx += bxf2*cos(th0f2+anglef2) - byf2*sin(th0f2+anglef2);
				by += bxf2*sin(th0f2+anglef2) + byf2*cos(th0f2+anglef2);
				bz += bzf2;
				//fprintf(wfile, "%le	%le	%le	%le	%le	%le\n",r,th*180./PI,z,bx,by,bz);
				br = bx*cos(th) + by*sin(th);
				bth = -bx*sin(th) + by*cos(th);
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n",r,th*180./PI,z,br,bth,bz);
			}
		}
	}
	fclose(wfile);
	free_latt(&latt_d);
	free_latt(&latt_f);
}

/*extern void create_circ_cartmap_triplet_from_cartmap(char *textfile, double cell_length, double theta_max_deg, char *latt_f_name, double xf1, double yf1, double tilt_magf1_deg, double shiftf, double scalef, char *latt_d_name, double xd, double yd, double shiftd, double scaled, double rmin, double rmax, int rstep, double thmin, double thmax, int thstep, double zmin, double zmax, int zstep)
{
	int i,j,k;
	double ac,half_theta_cell,eb,oa,ob,thetad,r0d,phid,thetaf1,r0f,phif1, thetaf2, phif2, bx_tilt,by_tilt, bx,by,bz,bxf1,byf1,bzf1,bxd,byd,bzd,bxf2,byf2,bzf2,theta_max=theta_max_deg*PI/180.,xd,yd,xf,yf,x,y,z,r1,th1,r,th;
	double step_r,step_th,step_z;
	struct Lattice latt_f, latt_d;
	FILE *wfile;
	
	ac = cell_length;
	half_theta_cell = theta_max/2.;
	eb = ac/2.-yf1;
	oa = ac/(2*sin(half_theta_cell));
	ob = ac/(2*tan(half_theta_cell));
	
	thetad = half_theta_cell;
	r0d = ob+shiftd;
	phid = 0.;
	
	thetaf1 = thetad - atan(eb/(ob+shiftf));
	r0f = (ob+shiftf)/(cos(thetad-thetaf1));
	phif1 = tilt_magf1_deg*PI/180. + thetad - thetaf1;
	
	thetaf2 = theta_max - thetaf1;
	phif2 = -phif1;
	
	step_r = comp_step(rmin,rmax,rstep);
	step_th = comp_step(thmin,thmax,thstep);
	step_z = comp_step(zmin,zmax,zstep);
	
	printf("D MAGNET!\n");
	load_lattice(&latt_d, latt_d_name);
	printf("F MAGNET!\n");
	load_lattice(&latt_f, latt_f_name);
	
	wfile=fopen(textfile,"w");
	fprintf(wfile, "%i	%i	%i\n", xstep,ystep,zstep);
	fprintf(wfile, "x(m),y,z,Bx(T),By,Bz\n");
	
	for(i=0;i<zstep;i++) {
		z=zmin+i*step_z;
		for(j=0;j<rstep;j++) {
			x = rmin + j*step_r;
			for(k=0;k<thstep;k++) {
				y = thmin + k*step_th;
				bxf1=0;byf1=0;bzf1=0;bxd=0;byd=0;bzd=0;bxf2=0;byf2=0;bzf2=0;
				//coordinate_global_to_rect_vffa(x, y, &xd, &yd, r0d, thetad, phid, shiftd, latt_d.cell[0].boun.ymax/2.); //careful not tested with shiftd!=0
				//r1 = sqrt(xd*xd+yd*yd);
				//th1 = atan_ratio(yd,xd);
				//map_neighbours3(&bxd, &byd, &bzd, theta_max, r1, th1, z, scaled, &(latt_d.cell[0]), NO);

				coordinate_global_to_rect_vffa(x, y, &xf, &yf, r0f, thetaf1, phif1, shiftf, latt_f.cell[0].boun.ymax/2.); //careful not tested with shiftf!=0
				r1 = sqrt(xf*xf+yf*yf);
				th1 = atan_ratio(yf,xf);
				map_neighbours3(&bxf1, &byf1, &bzf1, theta_max, r1, th1, z, scalef, &(latt_f.cell[0]), NO);

				//coordinate_global_to_rect_vffa(x, y, &xf, &yf, r0f, thetaf2, phif2, shiftf, latt_f.cell[0].boun.ymax/2.); //careful not tested with shiftf!=0
				//r1 = sqrt(xf*xf+yf*yf);
				//th1 = atan_ratio(yf,xf);
				//map_neighbours3(&bxf2, &byf2, &bzf2, theta_max, r1, th1, z, scalef, &(latt_f.cell[0]), NO);
				
				bx_tilt = bxf1+bxd+bxf2; by_tilt = byf1+byd+byf2; bz = bzf1+bzd+bzf2; 
				bx = bx_tilt*cos(half_theta_cell) - by_tilt*sin(half_theta_cell);
				by = by_tilt*cos(half_theta_cell) + bx_tilt*sin(half_theta_cell);
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n",x,y,z,bx,by,bz);
				//fprintf(wfile, "%le	%le	%le	%le	%le	%le\n",r,th*180./PI,z,bx,by,bz);
			}
		}
	}
	fclose(wfile);
	free_latt(&latt_d);
	free_latt(&latt_f);
}
//*/

extern void map_neighbours3(double *bx, double *by, double *bz, double theta_tot, double r, double th, double z, double scale, struct Cell *cell, int kill_notalive)
{
	double x,y,znew,x1,x2,y1,y2,bx0=0,by0=0,bz0=0,bx1=0,by1=0,bz1=0,bx2=0,by2=0,bz2=0;
	
	x = r*cos(th);
	y = r*sin(th);
	x1 = r*cos(th+theta_tot);
	y1 = r*sin(th+theta_tot);
	x2 = r*cos(th-theta_tot);
	y2 = r*sin(th-theta_tot);
	
	//printf("r=%le, th=%le, (x,y)=%le,%le\n", r, th*180./PI, x,y);
	//printf("r=%le, th=%le, (x,y)=%le,%le, (x1,y1)=%le,%le\n", r, th*180./PI, x,y, x1,y1);
	
	if(get_bfield(x, y, z, &bx0, &by0, &bz0, cell, get_bfield_cartesian_map) != ALIVE) {
		//printf("r=%le, th=%le, (x,y)=%le,%le, (x1,y1)=%le,%le\n", r, th*180./PI, x,y, x1,y1);
		if(kill_notalive==YES) errorstop("map_neighbours: get_bfield != ALIVE\n");
		else {
			bx0=0;by0=0;bz0=0;
		}
	}
	bx0*=scale;by0*=scale;bz0*=scale;
	*bx = bx0; *by = by0; *bz = bz0;
	//printf("b0=(%le,%le,%le)\n",bx0,by0,bz0);
	//printf("r=%le, th=%le, (x,y)=%le,%le, (x1,y1)=%le,%le", r, th*180./PI, x,y, x1,y1);
	/*if(get_bfield(x1, y1, z, &bx1, &by1, &bz1, cell, get_bfield_cartesian_map) == ALIVE) { //cell before
		//printf("b1=(%le,%le,%le)\n",bx1,by1,bz1);
		bx1*=scale;by1*=scale;bz1*=scale;
		*bx += bx1*cos(theta_tot) + by1*sin(theta_tot);
		*by += by1*cos(theta_tot) - bx1*sin(theta_tot);
		*bz += bz1;
		//printf("   YES");
	}
	//printf("b=(%le,%le,%le)\n", bx,by,bz);
	if(get_bfield(x2, y2, z, &bx2, &by2, &bz2, cell, get_bfield_cartesian_map) == ALIVE) { //cell after
		//printf("b2=(%le,%le,%le)\n",bx2,by2,bz2);
		bx2*=scale;by2*=scale;bz2*=scale;
		*bx += bx2*cos(theta_tot) - by2*sin(theta_tot);
		*by += by2*cos(theta_tot) + bx2*sin(theta_tot);
		*bz += bz2;
	}//*/
	//printf("tot (%le,%le,%le), %le + %le\n", *bx,*by,*bz,  bx2*cos(theta_tot), -by2*sin(theta_tot));
}

//x1,y1 = centre of the magnet in local (magnet) coordinates in the magnet map coordinates (exemple: in a field map made of 1 magnet, centered on the magnet, y1=ymap_max/2)
//angle in rad
extern void coordinate_global_to_rect_vffa(double x_g, double y_g, double *x_mag, double *y_mag, double r_mag_cent, double th_mag_cent, double tilt_mag, double x_1, double y_1)
{
	double x_0,y_0, angle;
	
	x_0 = r_mag_cent*cos(th_mag_cent);//(=rmagnet*cos(thmagnet))
	y_0 = r_mag_cent*sin(th_mag_cent);
	angle = th_mag_cent+tilt_mag;
	*x_mag = (y_g-y_0)*sin(angle) + (x_g-x_0)*cos(angle)+x_1;
	*y_mag = (y_g-y_0)*cos(angle) - (x_g-x_0)*sin(angle)+y_1;
}

extern void compare_symetry_field_straight(char *textout, double xmin, double xmax, int nbx, double zmin, double zmax, int nbz, struct Cell *cell)
{
	int i,j;
	double bx1,by1,bz1,bx2,by2,bz2,xstep, zstep,x,z;
	FILE *wfile;
	
	xstep = comp_step(xmin,xmax,nbx);
	zstep = comp_step(zmin,zmax,nbz);
	wfile = fopen(textout, "w");

	for(i=0;i<nbx;i++) {
		x = xmin+i*xstep;
		for(j=0;j<nbz;j++) {
			z = zmin+j*zstep;
			if(get_bfield(x, 0, z, &bx1, &by1, &bz1, cell, get_bfield_cartesian_map) != ALIVE) errorstop("compare_symetry_field_straight 1: get_bfield != ALIVE\n");
			if(get_bfield(x, cell->boun.ymax, z, &bx2, &by2, &bz2, cell, get_bfield_cartesian_map) != ALIVE) errorstop("compare_symetry_field_straight 2: get_bfield != ALIVE\n");
			fprintf(wfile, "%le	%le	%le	%le	%le\n", x,z,bx1-bx2,by1-by2,bz1-bz2);
		}
	}
	fclose(wfile);
}

extern int calc_tune_eigenvalue(double *qu, double *qv, double ampx, double ampxp, double ampz, double ampzp, struct Particle *part, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, char *outfile)
{
	int i,j, doyoudecouple,n=4, doyoukeepzeropointfive=YES;
	double e_tot, beta, gamma, det;
	double m4d[4][4], band_bias;
	double complex eigenvaluesW[n], m_1dim[n*n], eigenvectorsVR[n*n];
	FILE *wfile;
	
	if(outfile != NULL || doyouprintf == YES) get_ebg_part(&e_tot, &beta, &gamma, part);
	if(doyouprintf == YES) {
		if(print_color==YES) COLOR("1;35");
		printf("\n--------------- Tunes from eigenvalues of equivalent transfer matrix ---------------");
		if(print_color==YES) COLOR("0");
		printf("\n");
		printf("At E(kin) = %le [eV], x,z_clo = (%le,%le [m], ux,uz_clo = (%le, %le)\n", (gamma - 1)*part->m0*CLIGHT2/(UNITCHARGE), part->x, part->z, part->ux, part->uz);
	}
	
	if(get_matrix_firstorder_4d(m4d, &band_bias, part, ampx, ampxp, ampz, ampzp, latt, *transport_part, doyouprintf) != TRUE) {
	//if(get_matrix_secondorder_4d(m4d, &band_bias, part, ampx, ampxp, ampz, ampzp, latt, *transport_part, doyouprintf, doyousymplect) != TRUE) {
		*qu=0;*qv=0;
		return FALSE;
	}
	det = matrix_det_4d(m4d, 4);
	for(i=0;i<n;i++) {
		for(j=0;j<n;j++) {
			m_1dim[i*n+j] = m4d[i][j];
		}
	}
	doyoudecouple = compute_eigenvector_parzen(eigenvectorsVR, eigenvaluesW, m_1dim);
	if(doyoudecouple == FALSE) {
		printf("warning: degenerated eigenvalues, not stable\n");
		*qu=0;*qv=0;
		return FALSE;
	}
	*qu = atan_ratio(cimag(eigenvaluesW[0]), creal(eigenvaluesW[0]))/(2.*PI);
	*qv = atan_ratio(cimag(eigenvaluesW[2]), creal(eigenvaluesW[2]))/(2.*PI);
	if(doyouprintf == YES) {
		printf("qu: %.15e, 1-qu:%.15e\n", *qu, 1.-fabs(*qu));
		printf("qv: %.15e, 1-qv:%.15e\n", *qv, 1.-fabs(*qv));
	}
	if(doyoukeepzeropointfive==YES) {
		if(*qu>0.5) *qu=1-(*qu);
		if(*qv>0.5) *qv=1-(*qv);
	}
	if(outfile != NULL) {
		wfile = fopen(outfile, "a");
		//fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n", e_tot, *qx/latt->periodicity, *qz/latt->periodicity, *qx, *qz, amp_x, amp_xprime, amp_z, amp_zprime, *betax, *alphax, *betaz, *alphaz, latt->cell[0].mpara[0][6]);
		//fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le\n", e_tot, *qu, *qv, ampx, ampxp, ampz, ampzp);
		//fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le\n", e_tot, *qu, *qv, ampx, ampxp, ampz, ampzp, latt->cell[0].mpara[0][6]); //convergence length study
		//fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le\n", e_tot, *qu, *qv, ampx, ampxp, ampz, ampzp, latt->cell[0].efben[0][6]); //convergence length study
		fprintf(wfile, "%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le\n", e_tot, *qu, *qv, creal(eigenvaluesW[0]), cimag(eigenvaluesW[0]), creal(eigenvaluesW[2]), cimag(eigenvaluesW[2]), ampx, ampxp, ampz, ampzp); //convergence length study
		//fprintf(wfile, "%le  %le  %le  %le  %le\n", -latt->cell[0].mpara[1][2]/latt->cell[0].mpara[2][2], latt->cell[0].efben[0][2], *qu, *qv, det);
		
		//if(latt->cell[0].nbcomp==4) fprintf(wfile, "%le  %le  %le  %le  %le  %le\n", -latt->cell[0].mpara[1][2]/latt->cell[0].mpara[2][2], latt->cell[0].mpara[0][4], *qu, *qv, det, part->brho*part->q*CLIGHT/(UNITCHARGE));
		//else if(latt->cell[0].nbcomp==2) fprintf(wfile, "%le  %le  %le  %le  %le  %le\n", -latt->cell[0].mpara[0][2]/latt->cell[0].mpara[1][2], latt->cell[0].mpara[0][4], *qu, *qv, det, part->brho*part->q*CLIGHT/(UNITCHARGE));
		//else errorstop("nb of components strange\n");
		
		//fprintf(wfile, "%le  %le  %le  %le\n", *qu, *qv, det, part->brho*part->q*CLIGHT/(UNITCHARGE));
		fclose(wfile);
	}
	return TRUE;
}

extern void max_straight_vffa_fdratio_study(char *textout, double fdstart, double fdend, int nbfd, struct Particle *part, struct Lattice *latt)
{
	char trackout[500], tuneout[500];
	int i;
	double qu, qv, fd, fdstep;
	FILE *wfile;
	
	fdstep = comp_step(fdstart, fdend, nbfd);
	sprintf(tuneout, "%s_tunes.dat", textout);
	wfile = fopen(tuneout, "w");
	fprintf(wfile, "-BD0/BF0 	 m [/m] 	 qu 	 qv 	 det(M) 	 Ptot [eV/c]\n");
	fclose(wfile);
	for(i=0;i<nbfd;i++) {
		fd = fdstart+i*fdstep;
		if(latt->cell[0].nbcomp==4) {
			latt->cell[0].mpara[1][2] = -(fd)*latt->cell[0].mpara[2][2];
			latt->cell[0].mpara[3][2] = latt->cell[0].mpara[1][2];
			printf("B0D=%lf, %i\n",latt->cell[0].mpara[1][2], i);
			printf("FD ratio: %lf\n\n", -latt->cell[0].mpara[1][2]/latt->cell[0].mpara[2][2]);
		}
		else if(latt->cell[0].nbcomp==2) {
			printf("2 components\n");
			latt->cell[0].mpara[0][2] = -(fd)*latt->cell[0].mpara[1][2];
			printf("B0D=%lf, F=%lf\n",latt->cell[0].mpara[0][2], latt->cell[0].mpara[1][2]);
			printf("FD ratio: %lf\n\n", -latt->cell[0].mpara[0][2]/latt->cell[0].mpara[1][2]);
		}
		else errorstop("nb of components strange\n");
		
		if(find_closed_orbite_xxp_zzp(part, &(part->x), &(part->ux), &(part->uy), &(part->z), &(part->uz), 1.e-10, latt, YES) != TRUE) errorstop("closed orbit not found");
		sprintf(trackout, "%s_tracking_fd=%.2f.dat", textout, fd);
		wfile = fopen(trackout, "w");
		fprintf(wfile, "s [m]\t hor [m] \t long [m] \t vert [m] \t B_hor [T] \t B_long [T] \t B_vert [T] \t P_hor/Ptot \t P_long/Ptot \t P_vert/Ptot \t Brho [T.m]\n");
		fclose(wfile);
		part_oneturn(part, latt, trackout);
		calc_tune_eigenvalue(&qu, &qv, 1.e-5, 1.e-5, 1.e-5, 1.e-5, part, latt, part_cross_latt, YES, tuneout);
	}
		
}


extern void compute_conv_lim_long_vffa_rect_str(struct Cell *cell, char *textout, int nby, double xmax_start)
{
	int i,j;
	double x,y,ystep, bx2,by2,bz2,bx4,by4,bz4,bx6,by6,bz6,bx8,by8,bz8,bx10,by10,bz10, xmin, xmax,x_x,x_y,x_z;
	double a4,a6,a8,a10;
	FILE *wfile;
	
	wfile = fopen(textout, "w");
	
	ystep = comp_step(0, cell->boun.ymax, nby);
	
	for(i=0;i<nby;i++) {
		y=i*ystep;
		//xmin=0;xmax=xmax_start,x=xmax;
		//for(j=0;j<50;j++) {
		//	cell->efben[0][2] = 2;
		//	//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
		//	get_bfield(x, y, cell->mpara[0][3], &bx2, &by2, &bz2, cell, add_field_comp_VFFA_rect_str_atan);
		//	cell->efben[0][2] = 10;
		//	//max2=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
		//	get_bfield(x, y, cell->mpara[0][3], &bx10, &by10, &bz10, cell, add_field_comp_VFFA_rect_str_atan);
		//	if(fabs(bx10-bx2)>fabs(bx2) || fabs(by10-by2)>fabs(by2) || fabs(bz10-bz2)>fabs(bz2)) xmax=x;
		//	//if(max2 > max1) xmax=xtest;
		//	else xmin=x;
		//	//printf("\txmin=%le,xmax=%le,max1=%le,max2=%le\n",xmin,xmax,max1,max2);
		//	if((xmax-xmin)<1.e-6) {
		//		fprintf(wfile, "%le	%le	%le\n",x,y,xmax-xmin);
		//		break;
		//	}
		//	x = (xmax+xmin)/2.;
		//}
		
		xmin=0;xmax=xmax_start,x=xmax;
		for(j=0;j<50;j++) {
			bx2=0;by2=0;bz2=0;bx10=0;by10=0;bz10=0;
			cell->efben[0][2] = 2;
			//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx2, &by2, &bz2, cell, add_field_comp_VFFA_rect_str_atan);
			cell->efben[0][2] = 4;
			//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx4, &by4, &bz4, cell, add_field_comp_VFFA_rect_str_atan);
			cell->efben[0][2] = 6;
			//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx6, &by6, &bz6, cell, add_field_comp_VFFA_rect_str_atan);
			cell->efben[0][2] = 8;
			//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx8, &by8, &bz8, cell, add_field_comp_VFFA_rect_str_atan);
			cell->efben[0][2] = 10;
			//max2=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx10, &by10, &bz10, cell, add_field_comp_VFFA_rect_str_atan);
			a4=fabs(bx4)-fabs(bx2);
			a6=fabs(bx6)-fabs(a4);
			a8=fabs(bx8)-fabs(a6);
			a10=fabs(bx10)-fabs(a8);
			//if(fabs(a10/a8)<1 && fabs(a8/a6)<1) xmin=x;
			if(fabs(bx10-bx2)>fabs(bx2)) xmax=x;
			//if(max2 > max1) xmax=xtest;
			else xmin=x;
			//if(y>0.874 && y<0.876) {
			//	printf("x,y=%lf,%lf, an+1/an=%le, %le, %le\n",x,y, a6/a4, a8/a6,a10/a8);
			//	printf("a=(%le,%le,%le,%le)\n", a4,a6,a8,a10);
			//	printf("x E [%le, %le]\n",xmin, xmax);
			//}
			//printf("\txmin=%le,xmax=%le,max1=%le,max2=%le\n",xmin,xmax,max1,max2);
			if((xmax-xmin)<1.e-6) {
				x_x = x;
				break;
			}
			x = (xmax+xmin)/2.;
			if(j==49) errorstop("problem!!");
		}
		xmin=0;xmax=xmax_start,x=xmax;
		for(j=0;j<50;j++) {
			bx2=0;by2=0;bz2=0;bx10=0;by10=0;bz10=0;
			cell->efben[0][2] = 2;
			//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx2, &by2, &bz2, cell, add_field_comp_VFFA_rect_str_atan);
			cell->efben[0][2] = 4;
			//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx4, &by4, &bz4, cell, add_field_comp_VFFA_rect_str_atan);
			cell->efben[0][2] = 6;
			//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx6, &by6, &bz6, cell, add_field_comp_VFFA_rect_str_atan);
			cell->efben[0][2] = 8;
			//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx8, &by8, &bz8, cell, add_field_comp_VFFA_rect_str_atan);
			cell->efben[0][2] = 10;
			//max2=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx10, &by10, &bz10, cell, add_field_comp_VFFA_rect_str_atan);
			a4=by4-by2;
			a6=by6-a4;
			a8=by8-a6;
			a10=by10-a8;
			//if(fabs(a10/a8)>1. && fabs(a8/a6)>1.) xmax=x;
			if(fabs(by10-by2)>fabs(by2)) xmax=x;
			//if(max2 > max1) xmax=xtest;
			else xmin=x;
			//printf("\txmin=%le,xmax=%le,max1=%le,max2=%le\n",xmin,xmax,max1,max2);
			if((xmax-xmin)<1.e-6) {
				x_y = x;
				break;
			}
			x = (xmax+xmin)/2.;
			if(j==49) errorstop("problem!!");
		}
		xmin=0;xmax=xmax_start,x=xmax;
		for(j=0;j<50;j++) {
			bx2=0;by2=0;bz2=0;bx10=0;by10=0;bz10=0;
			cell->efben[0][2] = 2;
			//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx2, &by2, &bz2, cell, add_field_comp_VFFA_rect_str_atan);
			cell->efben[0][2] = 4;
			//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx4, &by4, &bz4, cell, add_field_comp_VFFA_rect_str_atan);
			cell->efben[0][2] = 6;
			//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx6, &by6, &bz6, cell, add_field_comp_VFFA_rect_str_atan);
			cell->efben[0][2] = 8;
			//max1=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx8, &by8, &bz8, cell, add_field_comp_VFFA_rect_str_atan);
			cell->efben[0][2] = 10;
			//max2=maxwell_test_max(x, y, temp_cell.mpara[0][3], dx, &temp_cell, add_contribution_comp, doyouglob);
			get_bfield(x, y, cell->mpara[0][3], &bx10, &by10, &bz10, cell, add_field_comp_VFFA_rect_str_atan);
			a4=bz4-bz2;
			a6=bz6-a4;
			a8=bz8-a6;
			a10=bz10-a8;
			//if(fabs(a10/a8)>1. && fabs(a8/a6)>1.) xmax=x;
			if(fabs(bz10-bz2)>fabs(bz2)) xmax=x;
			//if(max2 > max1) xmax=xtest;
			else xmin=x;
			//printf("\txmin=%le,xmax=%le,max1=%le,max2=%le\n",xmin,xmax,max1,max2);
			if((xmax-xmin)<1.e-6) {
				x_z = x;
				break;
			}
			x = (xmax+xmin)/2.;
			if(j==49) errorstop("problem!!");
		}
		fprintf(wfile, "%le	%le	%le	%le\n",y,x_x,x_y,x_z);
	}
	fclose(wfile);
}


extern void compute_str_map(char *filename, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double unitlength, double bscale,
	int loop_order, int line_order, int header_skip, double xstart, double xend, int nbx, double ystart, double yend, int nby, double zstart, double zend, int nbz)
{
	int i,i1,i2,i3,n1,n2,n3,ix,iy,iz;
	double x,y,z,bx,by,bz, xstep,ystep,zstep,coord[4],btemp[4];
	
	FILE *wfile;
	wfile = fopen(filename,"w");
	fprintf(wfile,"%i     %i     %i\n",nbx,nby,nbz);
	for(i=0;i<header_skip-1;i++) fprintf(wfile,"\n"); //jump header_skip -1 lines of header
	
	//if(nbz>1) zstep = (zend-zstart)/(nbz-1);
	//else zstep = 0.;
	zstep = comp_step(zstart, zend, nbz);
	//if(nbr>1) rstep = (rend-rstart)/(nbr-1);
	//else rstep = 0.;
	xstep = comp_step(xstart, xend, nbx);
	//if(nbth>1) thstep = (thend-thstart)/(nbth-1);
	//else thstep = 0.;
	ystep = comp_step(ystart, yend, nby);
	
	fprintf(wfile,"(%i, %i, %i) nodes\n",nbx,nbz,nby);
	fprintf(wfile,"from (%le, %le, %le) to (%le, %le, %le)\n",xstart, zstart, ystart, xend, zend, yend);
	fprintf(wfile,"cartesian grid with step (%le [m], %le [m], %le [m])\n", xstep,zstep,ystep);
	fprintf(wfile, "x (hor) [m], y (vert) [m], z (long) [m], Bx [T], By, Bz\n");
	
	if(loop_order/100%10 == 1) n1=nbx;
	else if(loop_order/100%10 == 2) n1=nby;
	else n1 = nbz;
	if(loop_order/10%10 == 1) n2=nbx;
	else if(loop_order/10%10 == 2) n2=nby;
	else n2 = nbz;
	if(loop_order%10 == 1) n3=nbx;
	else if(loop_order%10 == 2) n3=nby;
	else n3 = nbz;
	
	for(i1 = 0; i1 < n1; i1++) {
		if(loop_order/100%10 == 1) ix=i1;
		else if(loop_order/100%10 == 2) iy=i1;
		else iz = i1;
		for(i2 = 0; i2 < n2; i2++) {
			if(loop_order/10%10 == 1) ix=i2;
			else if(loop_order/10%10 == 2) iy=i2;
			else iz = i2;
			for(i3 = 0; i3 < n3; i3++) {
				if(loop_order%10 == 1) ix=i3;
				else if(loop_order%10 == 2) iy=i3;
				else iz = i3;
				x = xstart + ix*xstep;
				y = ystart + iy*ystep;
				z = zstart + iz*zstep;
				//printf("(x,y,z)=(%lf,%lf,%lf)\n",x,y,z);
				//printf("r=%lf, th=%le, z=%le\n",r,th,z);
				if(get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp)==ALIVE) {
					if(line_order/100%10 == 1) {
						coord[1] = x;
						btemp[1] = bx;
					}
					else if(line_order/100%10 == 2) {
						coord[1] = y;
						btemp[1] = by;
					}
					else {
						coord[1] = z;
						btemp[1] = bz;
					}
					if(line_order/10%10 == 1) {
						coord[2] = x;
						btemp[2] = bx;
					}
					else if(line_order/10%10 == 2) {
						coord[2] = y;
						btemp[2] = by;
					}
					else {
						coord[2] = z;
						btemp[2] = bz;
					}
					if(line_order%10 == 1) {
						coord[3] = x;
						btemp[3] = bx;
					}
					else if(line_order%10 == 2) {
						coord[3] = y;
						btemp[3] = by;
					}
					else {
						coord[3] = z;
						btemp[3] = bz;
					}
					fprintf(wfile,"%lf	%lf	%lf	%lf	%lf	%lf\n", coord[1]*unitlength, coord[2]*unitlength, coord[3]*unitlength, btemp[1]/bscale, btemp[2]/bscale, btemp[3]/bscale);
				}
				else errorstop("bfield not alive");
			}
		}
	}
	fclose(wfile);
}

extern void maxwell_test_across_str_cell(char *textfile, double dx, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double z, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	char name_div[500], name_cx[500],name_cy[500],name_cz[500];
	int i,j;
	double x,y, stepx, stepy,divb,curlbx,curlby,curlbz;
	FILE *wfile_div, *wfile_cx, *wfile_cy, *wfile_cz;
	
	sprintf(name_div,"%s_div.dat", textfile);
	sprintf(name_cx,"%s_curlx.dat", textfile);
	sprintf(name_cy,"%s_curly.dat", textfile);
	sprintf(name_cz,"%s_curlz.dat", textfile);
	
	wfile_div = fopen(name_div, "w");
	wfile_cx  = fopen(name_cx, "w");
	wfile_cy  = fopen(name_cy, "w");
	wfile_cz  = fopen(name_cz, "w");
	stepx = comp_step(xmin, xmax, nbx);
	stepy = comp_step(ymin, ymax, nby);
	
	for(i=0;i<nbx;i++) {
		x = xmin +i*stepx;
		for(j=0;j<nby;j++) {
			y=ymin+j*stepy;
			if(maxwell_test(x, y, z, dx, dx, dx, cell, add_contribution_comp, &divb, &curlbx, &curlby, &curlbz,NO)!=TRUE) { //, &dbzdx, &dbxdz, &dbzdy, &dbydz, &dbxdy, &dbydx);
				divb = 0;
				curlbx = 0;
				curlby = 0;
				curlbz = 0;
			}
			fprintf(wfile_div, "%le	%le	%le\n",y,x,divb);
			fprintf(wfile_cx , "%le	%le	%le\n",y,x,curlbx);
			fprintf(wfile_cy , "%le	%le	%le\n",y,x,curlby);
			fprintf(wfile_cz , "%le	%le	%le\n",y,x,curlbz);
		}
		fprintf(wfile_div, "\n");
		fprintf(wfile_cx , "\n");
		fprintf(wfile_cy , "\n");
		fprintf(wfile_cz , "\n");
	}
	fclose(wfile_div);
	fclose(wfile_cx);
	fclose(wfile_cy);
	fclose(wfile_cz);
}

extern void compute_fodo_circmap(char *textfile, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, 
	char *latt_d_name, double r0d, double th0d_deg, double x0d, double y0d, double z0d, double scaled, char *latt_f_name, double r0f, double th0f_deg, double x0f, double y0f, double z0f, double scalef)
{
	int i,j,k;
	double step_r, step_th, step_z,z,r,th, theta_min, theta_max, theta_tot,bx,by,bz,br,bth, bx1,by1,bz1,th0d=th0d_deg*PI/180.,th0f=th0f_deg*PI/180., bxs,bys;
	struct Lattice latt_f, latt_d;
	FILE *wfile=NULL;
	
	doyou_long_boun_sym = NO;
	
	load_lattice(&latt_d, latt_d_name);
	load_lattice(&latt_f, latt_f_name);
	
	theta_min = theta_min_deg*PI/180.;
	theta_max = theta_max_deg*PI/180.;
	theta_tot = fabs(theta_min)+ fabs(theta_max);
	step_r =  comp_step(rmin, rmax, rstep);
	step_th = comp_step(theta_min, theta_max, thstep);
	step_z =  comp_step(zmin, zmax, zstep);
	
	printf("zmin = %le, zmax = %le, step_z = %le\n", zmin, zmin+(zstep-1)*step_z, step_z);
	printf("rmin = %le, rmax = %le, step_r = %le\n", rmin, rmin+(rstep-1)*step_r, step_r);
	printf("thmin = %le, thmax = %le, step_th = %le\n", theta_min, theta_min+(thstep-1)*step_th, step_th);
	
	wfile = fopen(textfile, "w");
	fprintf(wfile, "%i %i %i\n",rstep, thstep, zstep);
	fprintf(wfile, "r(m),theta(deg),z(m),Br(T),Bth,Bz\n");
	for(i=0;i<zstep;i++) {
		z = zmin + i*step_z;
		for(j=0;j<rstep;j++) {
			r = rmin + j*step_r;
			for(k=0;k<thstep;k++) {
				th = theta_min + k*step_th;
				//if(z>0.025 && z<0.035 && r>7.685 && r<7.695 && th>3.59*PI/180. && th<3.69*PI/180.) {
				
				map_neighbours2(&bx, &by, &bz, x0d, y0d, r0d, th0d, theta_tot, r, th, z-z0d, scaled, &(latt_d.cell[0]), NO);//YES???
				bxs = bx*cos(th0d) - by*sin(th0d);
				bys = by*cos(th0d) + bx*sin(th0d);
				map_neighbours2(&bx1, &by1, &bz1, x0f, y0f, r0f, th0f, theta_tot, r, th, z-z0f, scalef, &(latt_f.cell[0]), NO);//YES???
				bz += bz1;
				bxs += bx1*cos(th0f) - by1*sin(th0f);
				bys += by1*cos(th0f) + bx1*sin(th0f);
				br = bxs*cos(th) + bys*sin(th); 
				bth = bys*cos(th) - bxs*sin(th); 
				//printf("r,th,z=(%le,%le,%le)\n",r,th*PI/180., z);
				//printf("BD:(%le, %le, %le)\n", bx, by, bz);
				//printf("BF:(%le, %le, %le)\n", bx1, by1, bz1);
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r, th*180./PI, z, br, bth, bz);
				
				//}
				//fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r, th*180./PI, z, bxs, bys, bz);
			}
		}
	}
	fclose(wfile);
	free_latt(&latt_d);
	free_latt(&latt_f);
	doyou_long_boun_sym = YES;
}

extern void wrapper_vffa_map_checker(char *pref_out, struct Lattice *latt, struct Particle *part, double emin_ev, double emax_ev, int nbstep)
{
	char name_scan_tune[500], name_track[500], aspect[500], xrange[100],  trajout_hor[500], trajout_vert[500], bout[500], tuneout[500];

	sprintf(name_scan_tune ,"%s_tune_scan.dat",pref_out);
	sprintf(name_track,"%s_trackout.dat",pref_out);
	sprintf(aspect,"%s_aspect.dat",pref_out);
	sprintf(trajout_hor ,"%s_traj_xy.eps",pref_out);
	sprintf(trajout_vert ,"%s_traj_xz.eps",pref_out);
	sprintf(bout ,"%s_btot.eps",pref_out);
	sprintf(tuneout ,"%s_tune.eps",pref_out);
	
	emptyfile(name_track);
	part_cross_latt(part, latt, name_track);
	printf("co\n");
	get_xrange(name_track, xrange);
	plot_traj(latt, name_track, aspect, "($3)", "($2)", "($2)", "($1)","lines lc 7 lw 2", "lines lt 1 lc 0 lw 2", "long [m]", "hor [m]", NULL, NULL, trajout_hor, "size ratio -1\nset mxtics 5\nset mytics 5");
	easyplot(name_track, "($3)", "($4)", "lines lc 7 lw 2", "long [m]", "vert [m]", NULL, NULL, trajout_vert, "grid\n set mxtics 5\nset mytics 5");
	easyplot3p(name_track,name_track,name_track, "1", "15", "1", "16", "1", "7", "lines lw 3 lc 2", "lines lw 3 lc 3", "lines lw 3 lc 7", 
	"horizontal", "longitudinal", "vertical", "s [m]", "B [T]", xrange, NULL, bout, "grid\nset mxtics 5\nset mytics 5");
	
	emptyfile(name_scan_tune);
	phase_adv_scan(name_scan_tune, latt, part, 1.e-10, emin_ev, emax_ev, nbstep, 1.e-5);
	easyplot2(name_scan_tune, name_scan_tune, "($4*1.e-6)", "($1)", "($4*1.e-6)", "2", "lines lt 1 lw 2 lc 0", "lines lt 2 lw 2 lc 0", NULL, NULL, "P [MeV/c]", "{/Symbol n}_u_,_v []", NULL, NULL, tuneout, NULL);
	
}

extern void reduce_alan_map(char *mapout, char *lattout, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
	int i, j, k;
	double newx,newy,newz;
	double x,y,z,bx,by,bz;
	FILE *wfile;
	FILE *rfile=NULL;
	
	//sanity tests
	if(xmin<-1.3-2*TINYLENGTH) errorstop("xmin too small\n");
	if(xmax>1.3+2*TINYLENGTH) errorstop("xmax too big\n");
	if(ymin<-5-2*TINYLENGTH) errorstop("ymin too small\n");
	if(ymax>5+2*TINYLENGTH) errorstop("ymax too big\n");
	if(zmin<0-2*TINYLENGTH) errorstop("zmin too small\n");
	if(zmax>0.6+2*TINYLENGTH) errorstop("zmax too big\n");
	
	newx = (xmax-xmin)/0.01+1.;
	newy = (ymax-ymin)/0.01+1.;
	newz = (zmax-zmin)/0.01+1.;
	
	rfile = fopen("inputs/vffa_map_alan_m1.3/fieldmap-large-m=1.3-g=0.15.dat", "r");
	if(rfile==NULL) errorstop("cannot open file");
	newline(rfile);
	wfile = fopen(mapout, "w");
	fprintf(wfile, "%.0f %.0f %.0f\n",newx,newy,newz);
	
	for(i=0;i<261;i++) {
		for(j=0;j<61;j++) {
			for(k=0;k<1001;k++) {
				fscanf(rfile, "%le %le %le %le %le %le",&x,&z,&y, &bx,&bz,&by);
				//printf("%le %le %le %le %le %le\n",x,z,y, bx,bz,by);
				if(x>xmin && x<xmax && y>ymin && y<ymax && z>zmin && z<zmax) fprintf(wfile, "%le %le %le %le %le %le\n",x,z,y, bx,bz,by);
			}
		}
	}//*/
	fclose(rfile);
	fclose(wfile);
	write_maplatt_file(lattout, "test map", 1, mapout, 0.001, ymax-ymin, xmin, xmax, (int) newx, ymin, ymax, (int) newy, zmin, zmax, (int) newz);
}

extern void reduce_alan_map2(char *mapout, char *lattout, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax)
{
	int i, j, k, count_i=3, count_j=3,count_k=3;
	double newx,newy,newz;
	double x,y,z,bx,by,bz;
	FILE *wfile;
	FILE *rfile=NULL;
	
	//sanity tests
	if(xmin<-1.3-2*TINYLENGTH) errorstop("xmin too small\n");
	if(xmax>1.3+2*TINYLENGTH) errorstop("xmax too big\n");
	if(ymin<-5-2*TINYLENGTH) errorstop("ymin too small\n");
	if(ymax>5+2*TINYLENGTH) errorstop("ymax too big\n");
	if(zmin<0-2*TINYLENGTH) errorstop("zmin too small\n");
	if(zmax>0.6+2*TINYLENGTH) errorstop("zmax too big\n");
	
	newx = (xmax-xmin)/0.04+1.;
	newy = (ymax-ymin)/0.04+1.;
	newz = (zmax-zmin)/0.04+1.;
	
	rfile = fopen("inputs/vffa_map_alan_m1.3/fieldmap-large-m=1.3-g=0.15.dat", "r");
	if(rfile==NULL) errorstop("cannot open file");
	newline(rfile);
	wfile = fopen(mapout, "w");
	fprintf(wfile, "%.0f %.0f %.0f\n",newx,newy,newz);
	
	for(i=0;i<261;i++) {
		count_i++;
		for(j=0;j<61;j++) {
			count_j++;
			for(k=0;k<1001;k++) {
				count_k++;
				fscanf(rfile, "%le %le %le %le %le %le",&x,&z,&y, &bx,&bz,&by);
				//printf("%le %le %le %le %le %le\n",x,z,y, bx,bz,by);
				//printf("count_i=%i, count_j=%i, count_k=%i\n", count_i, count_j, count_k);
				if(x>xmin && x<xmax && y>ymin && y<ymax && z>zmin && z<zmax && count_i==4 && count_j==4 && count_k==4) {
					//printf("\n\nENTER !!\n\n\n");
					fprintf(wfile, "%le %le %le %le %le %le\n",x,z,y, bx,bz,by);
				}
				if(count_k==4) count_k=0;
			}
			if(count_j==4) count_j=0;
		}
		if(count_i==4) count_i=0;
	}//*/
	fclose(rfile);
	fclose(wfile);
	write_maplatt_file(lattout, "test map", 1, mapout, 0.001, ymax-ymin, xmin, xmax, (int) newx, ymin, ymax, (int) newy, zmin, zmax, (int) newz);
}

extern void write_maplatt_cyldeg_file(char *lattfile, char *title, int periodicity, char *mapfile, double stepsize, double cell_length, double rmin, double rmax, int nbr, double thmin, double thmax, int nbth, double zmin, double zmax, int nbz)
{
	FILE *wfile=NULL;
	
	wfile = fopen(lattfile, "w");
	if(wfile==NULL) errorstop("cannot create wfile\n");
	
	fprintf(wfile, "%s\n\n", title);
	fprintf(wfile, "%i				//lattice periodicity [int]\n", periodicity);
	fprintf(wfile, "1				//number of elements [int]\n");
	fprintf(wfile, "1				//number of cell types [int]\n\n\n");
	fprintf(wfile, "1 field-cylmap	nosym	cyldeg	%s\n", mapfile);
	//fprintf(wfile, "1 field-cartmap	nosym	cart %s\n", mapfile);
	fprintf(wfile, "%lf			//tracking step size [m]\n", stepsize);
	fprintf(wfile, "%lf						//cell total length [m]\n", cell_length);
	fprintf(wfile, "1	1		//unitlength/m (1000=[mm], 100=[cm]); B unit (1.e-4 Gauss, 1: Tesla)\n");
	fprintf(wfile, "312	123			//loop order, line order (123:(x,y,z);231:(y,z,x))\n");
	//fprintf(wfile, "132	132			//loop order, line order (123:(x,y,z);231:(y,z,x))\n");
	fprintf(wfile, "2				// header skip (nb of lines jumped)\n");
	//fprintf(wfile, "1				// header skip (nb of lines jumped)\n");
	fprintf(wfile, "0				//delta x (not checked!) (in map unitlength)\n\n");
	fprintf(wfile, "%lf		%lf		%i			// rmin, rmax (in map unitlength), nbsteps_r\n", rmin, rmax, nbr);
	fprintf(wfile, "%lf		%lf		%i			// thmin, thmax (in [deg.]), nbsteps_th\n", thmin, thmax, nbth);
	fprintf(wfile, "%lf		%lf		%i			// zmin, zmax (in map unitlength), nbsteps_z\n\n", zmin, zmax, nbz);
	fclose(wfile);
}

extern void wrapper_kieran_map(char *outfile, double xmax, double ymax)
{
	char mapout[500], lattout[500], fodo_map[500], fodo_latt[500], trackout[500];
	double qx,qz;
	struct Lattice latt;
	struct Beam beam;
	FILE *wfile;
	
	sprintf(mapout, "inputs/vffa_map_alan_m1.3/small_map/x_%lf_y_%lf_1mag.dat",xmax,ymax);
	sprintf(lattout, "inputs/vffa_map_alan_m1.3/small_map/x_%lf_y_%lf_1mag.latt",xmax,ymax);
	sprintf(fodo_map, "inputs/vffa_map_alan_m1.3/small_map/x_%lf_y_%lf_fodo.dat",xmax,ymax);
	sprintf(fodo_latt, "inputs/vffa_map_alan_m1.3/small_map/x_%lf_y_%lf_fodo.latt",xmax,ymax);
	sprintf(trackout, "inputs/vffa_map_alan_m1.3/small_map/x_%lf_y_%lf_trackout.dat",xmax,ymax);
	
	wfile = fopen(outfile, "a");
	printf("\n\n\n\nlong %lf, hor %lf\n\n\n\n", ymax, xmax);
	reduce_alan_map(mapout, lattout, -xmax-TINYLENGTH, xmax+TINYLENGTH, -ymax-TINYLENGTH, ymax+TINYLENGTH, 0-TINYLENGTH, 0.6+TINYLENGTH);
	//compute_fodo_circmap(fodo_map, 7.6, 8.2, 61, 0, 18, 251, 0, 0.6, 61, lattout, 7.8413779601, 4.5, 0, 5, 0, 0.165, "inputs/vffa_map_alan_m1.3/str_1mag_map.latt", 8.0413779601, 13.5, 0, 5, 0, -0.52);
	compute_fodo_circmap(fodo_map, 7.6, 8.2, 61, 0, 18, 251, 0, 0.6, 61, lattout, 7.8413779601, 4.5, 0, ymax, 0, 0.165, lattout, 8.0413779601, 13.5, 0, ymax, 0, -0.52);
	write_maplatt_cyldeg_file(fodo_latt, "test", 1, fodo_map, 0.0001, 18., 7.6, 8.2, 61, 0, 18, 251, 0, 0.6, 61);

	load_lattice(&latt, fodo_latt);
	load_beam(&beam, "inputs/vffa_map_alan_m1.3/proton7mev.beam", &latt, YES);
	write_beam("inputs/vffa_map_alan_m1.3/proton7mev.beam", &(beam.part[0]), 938.27, 1, 7.e6);
	emptyfile(trackout);
	part_cross_latt(&(beam.part[0]), &latt, trackout);
	
	calc_tune_eigenvalue(&qx, &qz, 1.e-5, 1.e-5, 1.e-5, 1.e-5, &(beam.part[0]), &latt, part_cross_latt, YES, NULL);
	fprintf(wfile, "%le	%le	%le	%le\n", xmax, ymax, qx, qz);
	free_latt(&latt);
	free_beam(&beam);
	fclose(wfile);
	remove(mapout);
	remove(lattout);
	remove(fodo_map);
	remove(fodo_latt);
}

extern void write_beam(char *beamfile, struct Particle *part, double m0, int q, double ekin_ev)
{
	FILE *wfile;
	
	wfile = fopen(beamfile, "w");
	fprintf(wfile, "particle %lf MeV\n\n", ekin_ev*1.e6);
	fprintf(wfile, "%lf	%i\n\n", m0, q);
	fprintf(wfile, "2\n\n");
	fprintf(wfile, "singlepart %le\n", ekin_ev);
	fprintf(wfile, "%le 0 %le\n", part->x, part->z);
	fprintf(wfile, "%le %le\n\n", atan_ratio(part->ux,part->uy)*180./PI, atan_ratio(part->uz,part->uy)*180./PI);
	fprintf(wfile, "putonco-xxpzzp 1.e-10\n");
	
	//fprintf(wfile, "putonco-nelmin 1.e-10	1.e-3\n\n");
	
	fclose(wfile);
}

extern void write_beam2(char *beamfile, struct Particle *part)
{
	double gamma, ekin_ev, e_tot,beta;
	FILE *wfile;
	
	get_ebg_part(&e_tot, &beta, &gamma, part);
	ekin_ev = (gamma - 1)*part->m0*CLIGHT2/(UNITCHARGE);
	wfile = fopen(beamfile, "w");
	fprintf(wfile, "particle %lf MeV\n\n", ekin_ev*1.e-6);
	fprintf(wfile, "%lf	%.0f\n\n", part->m0*1.e-6*CLIGHT2/(UNITCHARGE), part->q/UNITCHARGE);
	fprintf(wfile, "2\n\n");
	fprintf(wfile, "singlepart %le\n", ekin_ev);
	fprintf(wfile, "%le 0 %le\n", part->x, part->z);
	fprintf(wfile, "%le %le\n\n", atan_ratio(part->ux,part->uy)*180./PI, atan_ratio(part->uz,part->uy)*180./PI);
	//fprintf(wfile, "putonco-xxpzzp 1.e-10\n");
	fprintf(wfile, "putonco-xxp 1.e-10\n");
	//fprintf(wfile, "putonco-nelmin 1.e-10	1.e-3\n\n");
	
	fclose(wfile);
}

extern void compute_fdf_circmap(char *textfile, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, 
	char *latt_d_name, double r0d, double th0d_deg, double x0d, double y0d, double z0d, double scaled, 
	char *latt_f_name, double r0f, double th0f_deg, double x0f, double y0f, double z0f, double scalef)
{
	int i,j,k;
	double step_r, step_th, step_z,z,r,th, theta_min, theta_max, theta_tot,bx,by,bz,br,bth, bx1,by1,bz1,bx2,by2,bz2,th0d=th0d_deg*PI/180.,th0f=th0f_deg*PI/180., bxs,bys;
	struct Lattice latt_f, latt_d;
	FILE *wfile=NULL;
	
	doyou_long_boun_sym = NO;
	
	load_lattice(&latt_d, latt_d_name);
	load_lattice(&latt_f, latt_f_name);
	
	theta_min = theta_min_deg*PI/180.;
	theta_max = theta_max_deg*PI/180.;
	theta_tot = fabs(theta_min)+ fabs(theta_max);
	step_r =  comp_step(rmin, rmax, rstep);
	step_th = comp_step(theta_min, theta_max, thstep);
	step_z =  comp_step(zmin, zmax, zstep);
	
	printf("zmin = %le, zmax = %le, step_z = %le\n", zmin, zmin+(zstep-1)*step_z, step_z);
	printf("rmin = %le, rmax = %le, step_r = %le\n", rmin, rmin+(rstep-1)*step_r, step_r);
	printf("thmin = %le, thmax = %le, step_th = %le\n", theta_min, theta_min+(thstep-1)*step_th, step_th);
	
	wfile = fopen(textfile, "w");
	fprintf(wfile, "%i %i %i\n",rstep, thstep, zstep);
	fprintf(wfile, "r(m),theta(deg),z(m),Br(T),Bth,Bz\n");
	for(i=0;i<zstep;i++) {
		z = zmin + i*step_z;
		for(j=0;j<rstep;j++) {
			r = rmin + j*step_r;
			for(k=0;k<thstep;k++) {
				th = theta_min + k*step_th;
				//if(z>0.025 && z<0.035 && r>7.685 && r<7.695 && th>3.59*PI/180. && th<3.69*PI/180.) {
				
				map_neighbours2(&bx, &by, &bz, x0d, y0d, r0d, th0d, theta_tot, r, th, z-z0d, scaled, &(latt_d.cell[0]), NO);//YES???
				bxs = bx*cos(th0d) - by*sin(th0d);
				bys = by*cos(th0d) + bx*sin(th0d);
				map_neighbours2(&bx1, &by1, &bz1, x0f, y0f, r0f, th0f, theta_tot, r, th, z-z0f, scalef, &(latt_f.cell[0]), NO);//YES???
				bz += bz1;
				bxs += bx1*cos(th0f) - by1*sin(th0f);
				bys += by1*cos(th0f) + bx1*sin(th0f);
				map_neighbours2(&bx2, &by2, &bz2, x0f, y0f, r0f, theta_tot-th0f, theta_tot, r, th, z-z0f, scalef, &(latt_f.cell[0]), NO);//YES???
				bz += bz2;
				bxs += bx2*cos(theta_tot-th0f) - by2*sin(theta_tot-th0f);
				bys += by2*cos(theta_tot-th0f) + bx2*sin(theta_tot-th0f);
				br = bxs*cos(th) + bys*sin(th); 
				bth = bys*cos(th) - bxs*sin(th); 
				//printf("r,th,z=(%le,%le,%le)\n",r,th*PI/180., z);
				//printf("BD:(%le, %le, %le)\n", bx, by, bz);
				//printf("BF:(%le, %le, %le)\n", bx1, by1, bz1);
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r, th*180./PI, z, br, bth, bz);
				
				//}
				//fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r, th*180./PI, z, bxs, bys, bz);
			}
		}
	}
	fclose(wfile);
	free_latt(&latt_d);
	free_latt(&latt_f);
	doyou_long_boun_sym = YES;
}

extern void compute_fodo_circmap2(char *textfile, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, 
	char *latt_d_name_low, double x0d_low, double y0d_low, char *latt_d_name_high, double x0d_high, double y0d_high, double z0d, double r0d, double th0d_deg, double scaled, 
	char *latt_f_name_low, double x0f_low, double y0f_low, char *latt_f_name_high, double x0f_high, double y0f_high, double z0f, double r0f, double th0f_deg, double scalef)
{
	int i,j,k;
	double step_r, step_th, step_z,z,r,th, theta_min, theta_max, theta_tot,bx,by,bz,br,bth, bx1,by1,bz1,th0d=th0d_deg*PI/180.,th0f=th0f_deg*PI/180., bxs,bys;
	struct Lattice latt_f_high, latt_f_low, latt_d_high, latt_d_low;
	FILE *wfile=NULL;
	
	doyou_long_boun_sym = NO;
	
	load_lattice(&latt_d_high, latt_d_name_high);
	load_lattice(&latt_f_high, latt_f_name_high);
	load_lattice(&latt_d_low, latt_d_name_low);
	load_lattice(&latt_f_low, latt_f_name_low);
	
	theta_min = theta_min_deg*PI/180.;
	theta_max = theta_max_deg*PI/180.;
	theta_tot = fabs(theta_min)+ fabs(theta_max);
	step_r =  comp_step(rmin, rmax, rstep);
	step_th = comp_step(theta_min, theta_max, thstep);
	step_z =  comp_step(zmin, zmax, zstep);
	
	printf("zmin = %le, zmax = %le, step_z = %le\n", zmin, zmin+(zstep-1)*step_z, step_z);
	printf("rmin = %le, rmax = %le, step_r = %le\n", rmin, rmin+(rstep-1)*step_r, step_r);
	printf("thmin = %le, thmax = %le, step_th = %le\n", theta_min, theta_min+(thstep-1)*step_th, step_th);
	
	wfile = fopen(textfile, "w");
	fprintf(wfile, "%i %i %i\n",rstep, thstep, zstep);
	fprintf(wfile, "r(m),theta(deg),z(m),Br(T),Bth,Bz\n");
	for(i=0;i<zstep;i++) {
		z = zmin + i*step_z;
		for(j=0;j<rstep;j++) {
			r = rmin + j*step_r;
			for(k=0;k<thstep;k++) {
				th = theta_min + k*step_th;
				//if(z>0.025 && z<0.035 && r>7.685 && r<7.695 && th>3.59*PI/180. && th<3.69*PI/180.) {
				map_neighbours_low_high(&bx, &by, &bz, x0d_high, y0d_high, &(latt_d_high.cell[0]), x0d_low, y0d_low, &(latt_d_low.cell[0]), r0d, th0d, theta_tot, r, th, z-z0d, scaled, NO);
				bxs = bx*cos(th0d) - by*sin(th0d);
				bys = by*cos(th0d) + bx*sin(th0d);
				map_neighbours_low_high(&bx1, &by1, &bz1, x0f_high, y0f_high, &(latt_f_high.cell[0]), x0f_low, y0f_low, &(latt_f_low.cell[0]), r0f, th0f, theta_tot, r, th, z-z0f, scalef, NO);
				bz += bz1;
				bxs += bx1*cos(th0f) - by1*sin(th0f);
				bys += by1*cos(th0f) + bx1*sin(th0f);
				br = bxs*cos(th) + bys*sin(th); 
				bth = bys*cos(th) - bxs*sin(th); 
				//printf("r,th,z=(%le,%le,%le)\n",r,th*PI/180., z);
				//printf("BD:(%le, %le, %le)\n", bx, by, bz);
				//printf("BF:(%le, %le, %le)\n", bx1, by1, bz1);
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r, th*180./PI, z, br, bth, bz);
				//}
				//fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r, th*180./PI, z, bxs, bys, bz);
			}
		}
	}
	fclose(wfile);
	free_latt(&latt_d_high);
	free_latt(&latt_f_high);
	free_latt(&latt_d_low);
	free_latt(&latt_f_low);
	doyou_long_boun_sym = YES;
}

extern void map_neighbours_low_high(double *bx, double *by, double *bz, double x0_high, double y0_high, struct Cell *cell_high, double x0_low, double y0_low, struct Cell *cell_low, double r0, double th0, double theta_tot, double r, double th, double z, double scale, int kill_notalive)
{
	double bx0=0,by0=0,bz0=0,bx1=0,by1=0,bz1=0,bx2=0,by2=0,bz2=0;
	double x_low, y_low, znew, x1_low, x2_low, y1_low, y2_low, x_high, y_high, x1_high, x2_high, y1_high, y2_high;
	
	convert_ff_coord_to_shinji(&x_high, &y_high, &znew, r, th, z, r0, th0, y0_high, x0_high);
	convert_ff_coord_to_shinji(&x1_high, &y1_high, &znew, r, th, z, r0, th0-theta_tot, y0_high, x0_high);
	convert_ff_coord_to_shinji(&x2_high, &y2_high, &znew, r, th, z, r0, th0+theta_tot, y0_high, x0_high);
	
	convert_ff_coord_to_shinji(&x_low, &y_low, &znew, r, th, z, r0, th0, y0_low, x0_low);
	convert_ff_coord_to_shinji(&x1_low, &y1_low, &znew, r, th, z, r0, th0-theta_tot, y0_low, x0_low);
	convert_ff_coord_to_shinji(&x2_low, &y2_low, &znew, r, th, z, r0, th0+theta_tot, y0_low, x0_low);
	//x = r*cos(th-th0)-r0+x0;
	//y = r*sin(th-th0)+y0;
	//x1 = r*cos(th-th0+theta_tot)-r0+x0;
	//y1 = r*sin(th-th0+theta_tot)+y0;
	//x2 = r*cos(th-th0-theta_tot)-r0+x0;
	//y2 = r*sin(th-th0-theta_tot)+y0;
	
	//printf("r=%le, th=%le, (x,y)=%le,%le\n", r, th*180./PI, x,y);
	//printf("r=%le, th=%le, (x,y)=%le,%le, (x1,y1)=%le,%le\n", r, th*180./PI, x,y, x1,y1);
	
	if(get_bfield(x_high, y_high, z, &bx0, &by0, &bz0, cell_high, get_bfield_cartesian_map) != ALIVE) {
		bx0=0;by0=0;bz0=0;
		if(get_bfield(x_low, y_low, z, &bx0, &by0, &bz0, cell_low, get_bfield_cartesian_map) != ALIVE) {
			if(kill_notalive==YES) errorstop("map_neighbours: get_bfield != ALIVE\n");
			else {
				bx0=0;by0=0;bz0=0;
			}
		}
		//printf("\nr=%le, th=%le, (x,y)=%le,%le, (x1,y1)=%le,%le, (x2,y2)=%le,%le\n", r, th*180./PI, x,y, x1,y1, x2, y2);
		//printf("b0=(%le,%le,%le)\n",bx0,by0,bz0);
	}
	//printf("\nr=%le, th=%le, (x,y,z)=%le,%le,%le, (x1,y1)=%le,%le, (x2,y2)=%le,%le\n", r, th*180./PI, x,y,z, x1,y1, x2, y2);
	//printf("b0=(%le,%le,%le)\n",bx0,by0,bz0);
	bx0*=scale;by0*=scale;bz0*=scale;
	*bx = bx0; *by = by0; *bz = bz0;
	//printf("b0=(%le,%le,%le)\n",bx0,by0,bz0);
	//printf("r=%le, th=%le, (x,y)=%le,%le, (x1,y1)=%le,%le", r, th*180./PI, x,y, x1,y1);
	if(get_bfield(x1_high, y1_high, z, &bx1, &by1, &bz1, cell_high, get_bfield_cartesian_map) != ALIVE) {  //cell before
		bx1=0;by1=0;bz1=0;
		if(get_bfield(x1_low, y1_low, z, &bx1, &by1, &bz1, cell_low, get_bfield_cartesian_map) != ALIVE) {
			bx1=0;by1=0;bz1=0;
		}
		//printf("b1=(%le,%le,%le)\n",bx1,by1,bz1);
	}
	bx1*=scale;by1*=scale;bz1*=scale;
	*bx += bx1*cos(theta_tot) + by1*sin(theta_tot);
	*by += by1*cos(theta_tot) - bx1*sin(theta_tot);
	*bz += bz1;
	//printf("b=(%le,%le,%le)\n", bx,by,bz);
	if(get_bfield(x2_high, y2_high, z, &bx2, &by2, &bz2, cell_high, get_bfield_cartesian_map) == ALIVE) {  //cell after
		bx2=0;by2=0;bz2=0;
		if(get_bfield(x2_low, y2_low, z, &bx2, &by2, &bz2, cell_low, get_bfield_cartesian_map) == ALIVE) {
			bx2=0;by2=0;bz2=0;
		}
		//printf("b2=(%le,%le,%le)\n",bx2,by2,bz2);
		
	}//*/
	bx2*=scale;by2*=scale;bz2*=scale;
	*bx += bx2*cos(theta_tot) - by2*sin(theta_tot);
	*by += by2*cos(theta_tot) + bx2*sin(theta_tot);
	*bz += bz2;
	//printf("tot (%le,%le,%le), %le + %le\n", *bx,*by,*bz,  bx2*cos(theta_tot), -by2*sin(theta_tot));
}

extern void wrapper_kieran_map2(char *outfile, double xmax, double ymax)
{
	char mapout[500], lattout[500], fodo_map[500], fodo_latt[500], trackout[500];
	double qx,qz;
	struct Lattice latt;
	struct Beam beam;
	FILE *wfile;
	
	sprintf(mapout, "inputs/vffa_map_alan_m1.3/doublemap/x_%lf_y_%lf_1mag.dat",xmax,ymax);
	sprintf(lattout, "inputs/vffa_map_alan_m1.3/doublemap/x_%lf_y_%lf_1mag.latt",xmax,ymax);
	sprintf(fodo_map, "inputs/vffa_map_alan_m1.3/doublemap/x_%lf_y_%lf_fodo.dat",xmax,ymax);
	sprintf(fodo_latt, "inputs/vffa_map_alan_m1.3/doublemap/x_%lf_y_%lf_fodo.latt",xmax,ymax);
	sprintf(trackout, "inputs/vffa_map_alan_m1.3/doublemap/x_%lf_y_%lf_trackout.dat",xmax,ymax);
	
	wfile = fopen(outfile, "a");
	printf("\n\n\n\nlong %lf, hor %lf\n\n\n\n", ymax, xmax);
	reduce_alan_map(mapout, lattout, -xmax-TINYLENGTH, xmax+TINYLENGTH, -ymax-TINYLENGTH, ymax+TINYLENGTH, 0-TINYLENGTH, 0.6+TINYLENGTH);
	compute_fodo_circmap2(fodo_map, 7.6, 8.2, 61, 0, 18, 251, 0, 0.6, 61, 
	"inputs/vffa_map_alan_m1.3/map_large_lowres.latt", 0, 4, lattout, 0, ymax, 0, 7.8413779601, 4.5, 0.165, 
	"inputs/vffa_map_alan_m1.3/map_large_lowres.latt", 0, 4, lattout, 0, ymax, 0, 8.0413779601, 13.5, -0.52);
	write_maplatt_cyldeg_file(fodo_latt, "test", 1, fodo_map, 0.0001, 18., 7.6, 8.2, 61, 0, 18, 251, 0, 0.6, 61);

	load_lattice(&latt, fodo_latt);
	load_beam(&beam, "inputs/vffa_map_alan_m1.3/proton7mev.beam", &latt, YES);
	write_beam("inputs/vffa_map_alan_m1.3/proton7mev.beam", &(beam.part[0]), 938.27, 1, 7.e6);
	emptyfile(trackout);
	part_cross_latt(&(beam.part[0]), &latt, trackout);
	
	calc_tune_eigenvalue(&qx, &qz, 1.e-5, 1.e-5, 1.e-5, 1.e-5, &(beam.part[0]), &latt, part_cross_latt, YES, NULL);
	fprintf(wfile, "%le	%le	%le	%le\n", xmax, ymax, qx, qz);
	free_latt(&latt);
	free_beam(&beam);
	fclose(wfile);
	remove(mapout);
	remove(lattout);
	remove(fodo_map);
	remove(fodo_latt);
}

extern void compute_fieldheatmap_btot(char *filename, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int),
							double x0, double xmax, int nbx, double y0, double ymax, int nby, double z0, double zmax, int nbz)
{
	int i,j,k;
	double x,y,z,bxtemp,bytemp,bx,by,bz,xstep,ystep,zstep, btot,w_coor;
	FILE *wfile;
	
	wfile=fopen(filename, "w");
	
		
	zstep = comp_step(z0, zmax, nbz);
	xstep = comp_step(x0, xmax, nbx);
	ystep = comp_step(y0, ymax, nby);
	
	printf("xstep=%le[m], ystep=%le[m], zstep=%le[m]\n",xstep,ystep,zstep);
	
	for(j=0;j<nbx;j++) {
		x = x0 + j*xstep;
		for(i=0;i<nbz;i++) {
			z = z0 + i*zstep;
			for(k=0;k<nby;k++) {
				y = y0 + k*ystep;
				
				if(get_bfield(x, y, z, &bxtemp, &bytemp, &bz, cell, add_contribution_comp)!=ALIVE) {
					//errorstop("field not alive\n");
					bx = 0.;
					by = 0.;
					bz = 0.;
				}
				if(nbz==1) w_coor = x;
				else w_coor = z;
				btot = sqrt(bx*bx+by*by+bz*bz);
				fprintf(wfile,"%le\t%le\t%le\n", y, w_coor, btot);
			}
			fprintf(wfile,"\n"); // to build file to plot easyplot3dmap
		}
	}
	fclose(wfile);
}

/*extern void compute_triplet_circmap2(char *textfile, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, 
	char *latt_d_name_low, double x0d_low, double y0d_low, char *latt_d_name_high, double x0d_high, double y0d_high, double z0d, double r0d, double th0d_deg, double scaled, 
	char *latt_f_name_low, double x0f_low, double y0f_low, char *latt_f_name_high, double x0f_high, double y0f_high, double z0f, double r0f, double th0f_deg, double scalef)
{
	int i,j,k;
	double step_r, step_th, step_z,z,r,th, theta_min, theta_max, theta_tot,bx,by,bz,br,bth, bx1,by1,bz1,th0d=th0d_deg*PI/180.,th0f=th0f_deg*PI/180., bxs,bys;
	struct Lattice latt_f_high, latt_f_low, latt_d_high, latt_d_low;
	FILE *wfile=NULL;
	
	doyou_long_boun_sym = NO;
	
	load_lattice(&latt_d_high, latt_d_name_high);
	load_lattice(&latt_f_high, latt_f_name_high);
	load_lattice(&latt_d_low, latt_d_name_low);
	load_lattice(&latt_f_low, latt_f_name_low);
	
	theta_min = theta_min_deg*PI/180.;
	theta_max = theta_max_deg*PI/180.;
	theta_tot = fabs(theta_min)+ fabs(theta_max);
	step_r =  comp_step(rmin, rmax, rstep);
	step_th = comp_step(theta_min, theta_max, thstep);
	step_z =  comp_step(zmin, zmax, zstep);
	
	printf("zmin = %le, zmax = %le, step_z = %le\n", zmin, zmin+(zstep-1)*step_z, step_z);
	printf("rmin = %le, rmax = %le, step_r = %le\n", rmin, rmin+(rstep-1)*step_r, step_r);
	printf("thmin = %le, thmax = %le, step_th = %le\n", theta_min, theta_min+(thstep-1)*step_th, step_th);
	
	wfile = fopen(textfile, "w");
	fprintf(wfile, "%i %i %i\n",rstep, thstep, zstep);
	fprintf(wfile, "r(m),theta(deg),z(m),Br(T),Bth,Bz\n");
	for(i=0;i<zstep;i++) {
		z = zmin + i*step_z;
		for(j=0;j<rstep;j++) {
			r = rmin + j*step_r;
			for(k=0;k<thstep;k++) {
				th = theta_min + k*step_th;
				//if(z>0.025 && z<0.035 && r>7.685 && r<7.695 && th>3.59*PI/180. && th<3.69*PI/180.) {
				map_neighbours_low_high(&bx, &by, &bz, x0d_high, y0d_high, &(latt_d_high.cell[0]), x0d_low, y0d_low, &(latt_d_low.cell[0]), r0d, th0d, theta_tot, r, th, z-z0d, scaled, NO);
				bxs = bx*cos(th0d) - by*sin(th0d);
				bys = by*cos(th0d) + bx*sin(th0d);
				map_neighbours_low_high(&bx1, &by1, &bz1, x0f_high, y0f_high, &(latt_f_high.cell[0]), x0f_low, y0f_low, &(latt_f_low.cell[0]), r0f, th0f, theta_tot, r, th, z-z0f, scalef, NO);
				bz += bz1;
				bxs += bx1*cos(th0f) - by1*sin(th0f);
				bys += by1*cos(th0f) + bx1*sin(th0f);
				br = bxs*cos(th) + bys*sin(th); 
				bth = bys*cos(th) - bxs*sin(th); 
				//printf("r,th,z=(%le,%le,%le)\n",r,th*PI/180., z);
				//printf("BD:(%le, %le, %le)\n", bx, by, bz);
				//printf("BF:(%le, %le, %le)\n", bx1, by1, bz1);
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r, th*180./PI, z, br, bth, bz);
				//}
				//fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r, th*180./PI, z, bxs, bys, bz);
			}
		}
	}
	fclose(wfile);
	free_latt(&latt_d_high);
	free_latt(&latt_f_high);
	free_latt(&latt_d_low);
	free_latt(&latt_f_low);
	doyou_long_boun_sym = YES;
}//*/

// m = log(b/b0)/(z-z0);
extern double general_m_computation(double b_z, double z, double b_z0, double z0)
{
	double m;
	
	if(z==z0) errorstop("in general_m_computation, z=z0");
	if(fabs(b_z0)>1.e-5) m = log(b_z/b_z0)/(z-z0);
	else m = 0;
	return m;
}

//m = 1/B*dB/dz
extern double local_m_computation(double b_zdz, double b_z, double dz)
{
	double m;
	
	if(dz==0) errorstop("in local_m_computation, dz=0");
	if(fabs(b_z)>1.e-5) m = (b_zdz-b_z)/(dz*b_z);
	else m = 0;
	return m;
}

// compute general m from Btot with z0 = zmin
//it is assumed the map is centered around the magnet, i.e. centre of magnet = (x0=0, y0=ymap_max/2)
extern void compute_m_rect_magnet(char *textout, char *text_lattmap, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz, double halfgap, double long_halflength, double dz)
{
	int i,j,k,flag_out_coil;
	double step_x, step_y, step_z,x,y,z,bx,by,bz,btot,bxdz,bydz,bzdz,btotdz,btotz0[nbx][nby],bxz0[nbx][nby],byz0[nbx][nby],bzz0[nbx][nby],m_gen,m_loc,mx_loc,my_loc,mz_loc, diff;
	//double mx_gen,my_gen,mz_gen;
	
	struct Lattice latt;
	FILE *wfile;
	
	step_x =  comp_step(xmin, xmax, nbx);
	step_y = comp_step(ymin, ymax, nby);
	step_z =  comp_step(zmin, zmax, nbz);
	
	printf("zmin = %le, zmax = %le, step_z = %le\n", zmin, zmin+(nbz-1)*step_z, step_z);
	printf("xmin = %le, xmax = %le, step_x = %le\n", xmin, xmin+(nbx-1)*step_x, step_x);
	printf("ymin = %le, ymax = %le, step_y = %le\n", ymin, ymin+(nby-1)*step_y, step_y);
	
	wfile = fopen(textout, "w");
	load_lattice(&latt, text_lattmap);
	
	if(test_cell_map(&(latt.cell[0]))==NO) errorstop("lattice not a map");
	
	for(i=0;i<nbz;i++) {
		z = zmin + i*step_z;
		for(j=0;j<nbx;j++) {
			x = xmin + j*step_x;
			for(k=0;k<nby;k++) {
				y = ymin + k*step_y;
				if(y>latt.cell[0].boun.ymax/2.-long_halflength && y<latt.cell[0].boun.ymax/2.+long_halflength) { //check if outside of the coils, then m_loc=m_gen=0
					if(x<-halfgap || x>halfgap) flag_out_coil = YES;
					else flag_out_coil = NO;
				}
				else flag_out_coil = NO;
				if(flag_out_coil == YES) {
					m_loc = 0;
					mx_loc = 0;
					my_loc = 0;
					mz_loc = 0;
					diff = 0;
					//m_gen = 0;
					//mx_gen = 0;
					//my_gen = 0;
					//mz_gen = 0;
				}
				else {
					if(get_bfield(x, y, z, &bx, &by, &bz, &(latt.cell[0]), get_bfield_cartesian_map) != ALIVE) errorstop("compute_m_rect_magnet: get_bfield != ALIVE\n");
					btot = sqrt(bx*bx+by*by+bz*bz);
					if(get_bfield(x, y, z+dz, &bxdz, &bydz, &bzdz, &(latt.cell[0]), get_bfield_cartesian_map) != ALIVE) errorstop("compute_m_rect_magnet: get_bfield dz != ALIVE\n");
					btotdz = sqrt(bxdz*bxdz+bydz*bydz+bzdz*bzdz);
					m_loc = local_m_computation(btotdz, btot, dz);
					mx_loc = local_m_computation(bxdz, bx, dz);
					my_loc = local_m_computation(bydz, by, dz);
					mz_loc = local_m_computation(bzdz, bz, dz);
					if(btot>1.e-5) diff = m_loc-(bx*bx*mx_loc+by*by*my_loc+bz*bz*mz_loc)/btot/btot;
					else diff = 0;
					//if(i==0) {
					//	btotz0[j][k] = btot;
					//	bxz0[j][k] = bx;
					//	byz0[j][k] = by;
					//	bzz0[j][k] = bz;
					//	m_gen = 0;
					//	mx_gen = 0;
					//	my_gen = 0;
					//	mz_gen = 0;
					//}
					//else {
					//	m_gen = general_m_computation(btot, z, btotz0[j][k], zmin);
					//	mx_gen = general_m_computation(bx, z, bxz0[j][k], zmin);
					//	my_gen = general_m_computation(by, z, byz0[j][k], zmin);
					//	mz_gen = general_m_computation(bz, z, bzz0[j][k], zmin);
					//}
				}
				//if(i!=0) fprintf(wfile, "%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le\n", x,y,z,m_loc,mx_loc,my_loc,mz_loc,m_gen,mx_gen,my_gen,mz_gen);
				fprintf(wfile, "%le	%le	%le	%le	%le	%le	%le	%le\n", x,y,z,m_loc,mx_loc,my_loc,mz_loc, diff);
			}
			fprintf(wfile,"\n");
		}
	}
	
	free_latt(&latt);
	fclose(wfile);
}

/*extern void compute_full_map_comsol_symmetries(char *text_in, char *text_out, double xmin, double xmax, int nb_x, double ymin, double ymax, int nb_y, double zmin, double zmax, int nb_z)
{
	int i,j,k;
	double x,y,z, xstep,ystep,zstep,a,b,c;
	struct Map map_read;
	FILE *rfile=NULL;
	FILE *wfile;
	xstep = comp_step(xmin, xmax, nb_x);
	ystep = comp_step(ymin, ymax, nb_y);
	zstep = comp_step(zmin, zmax, nb_z);
	
	map_read = allocmap(nb_x, nb_y, nb_z);
	
	
	rfile = fopen(text_in, "r");
	if(rfile==NULL) errorstop("cannot open read file");
	for(i=0;i<9;i++) newline(rfile);
	
	for(i=0;i<nb_y;i++) {
		for(j=0;j<nb_z;j++) {
			for(k=0;k<nb_x;k++) {
				fscanf(rfile, "%le	%le	%le	%le	%le	%le", &x,&z,&y, &a, &b, &c);
				map_read.node[k][i][j].b[0] = a;
				map_read.node[k][i][j].b[1] = c;
				map_read.node[k][i][j].b[2] = b;
				map_read.node[k][i][j].coord[0] = x;
				map_read.node[k][i][j].coord[1] = y;
				map_read.node[k][i][j].coord[2] = z;
				//fscanf(rfile, "%le	%le	%le	%le	%le	%le", &x,&z,&y, &bx[k][i][j], &bz[k][i][j], &by[k][i][j]);
				//printf("%le	%le	%le	%le	%le	%le, %i %i %i\n", x,z,y, a, b, c,i,j,k);
				//printf("%le	%le	%le	%le	%le	%le\n", x,y,z, bx[k][i][j], by[k][i][j], bz[k][i][j]);
			}
		}
	}
	fclose(rfile);

	wfile = fopen(text_out, "w");
	fprintf(wfile,"%i	%i	%i\n",nb_x,nb_y,nb_z);
	fprintf(wfile,"coord:(%lf,%lf, %lf)\n",-map_read.node[nb_x-1-0][0][0].coord[0], map_read.node[nb_x-1-0][0][0].coord[1], map_read.node[nb_x-1-0][0][0].coord[2]);
	//for(j=0;j<nb_z;j++) { //loop order: z,y,x 321, line order xyz 123
		//z = zmin+j*zstep;
	z=zmin;
		for(i=0;i<nb_y;i++) {
			y = ymin +i*ystep; //y from -1 to 0
			for(k=0;k<nb_x;k++) {
				x = -xmax + k*xstep;
				//printf("x=%lf\n",x);
				//fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", x,y,z, -map_read.node[nb_x-1-k][i][j].b[0], map_read.node[nb_x-1-k][i][j].b[1], map_read.node[nb_x-1-k][i][j].b[2]);
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", -map_read.node[nb_x-1-k][i][j].coord[0],map_read.node[nb_x-1-k][i][j].coord[1],map_read.node[nb_x-1-k][i][j].coord[2], -map_read.node[nb_x-1-k][i][j].b[0], map_read.node[nb_x-1-k][i][j].b[1], map_read.node[nb_x-1-k][i][j].b[2]);
			}
			for(k=0;k<nb_x;k++) {
				x = xmin + k*xstep;
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", x,y,z, map_read.node[k][i][j].b[0], map_read.node[k][i][j].b[1], map_read.node[k][i][j].b[2]);
			}
		}
		/*for(i=nb_y-1;i>-1;i--) {
			y = -ymin -i*ystep; // y from +0.01 to 1
			for(k=0;k<nb_x;k++) {
				x = -xmax + k*xstep;
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", x,y,z, -map_read.node[nb_x-k][i][j].b[0], -map_read.node[nb_x-k][i][j].b[1], map_read.node[nb_x-k][i][j].b[2]);
			}
			for(k=0;k<nb_x;k++) {
				x = xmin + k*xstep;
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", x,y,z, map_read.node[k][i][j].b[0], -map_read.node[k][i][j].b[1], map_read.node[k][i][j].b[2]);
			}
		}
		//}
	
	free_map(&(map_read), nb_x, nb_y);
	fclose(wfile);
}//*/

//works only for ymin=0, and y>0
extern void compute_full_map_comsol_symmetries(char *latt_quarter, char *text_out, double xmin, double xmax, int nb_x, double ymin, double ymax, int nb_y, double zmin, double zmax, int nb_z)
{
	int i,j,k;
	double x,y,z, xstep,ystep,zstep,bx,by,bz;
	struct Lattice latt;
	FILE *wfile;
	
	load_lattice(&latt, latt_quarter);
	
	xstep = comp_step(xmin, xmax, nb_x);
	ystep = comp_step(ymin, ymax, nb_y);
	zstep = comp_step(zmin, zmax, nb_z);
	
	wfile = fopen(text_out, "w");
	fprintf(wfile,"%i	%i	%i\n", 2*nb_x-1, 2*nb_y-1, nb_z);
	fprintf(wfile, "loop order z,y,x (321)\n");
	fprintf(wfile, "line order x, y, z [m], Bx, By, Bz [T] (123)\n");
	
	for(i=0;i<nb_z;i++) { //loop order: z,y,x 321, line order xyz 123
		z = zmin+i*zstep; //NO SYMMETRY
		for(j=0;j<nb_y;j++) {
			y = ymax - j*ystep; //y from 2 to 0 because y>0 in the original map
			for(k=0;k<nb_x;k++) {
				x = xmin + k*xstep; //because x<0 in the original map
				//x = xmax -k*xstep; //change to this if x>0 in the original map and change fprintf as well!!
				if(get_bfield(x, y, z, &bx, &by, &bz, &(latt.cell[0]), get_bfield_cartesian_map) != ALIVE) errorstop("compute_full_map_comsol_symmetries: x<0, y<0 get_bfield != ALIVE\n");
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", x, -y, z, bx, -by, bz); //in the original map, x<0, y>0
			}
			for(k=1;k<nb_x;k++) {
				x = xmax - k*xstep;
				if(get_bfield(x, y, z, &bx, &by, &bz, &(latt.cell[0]), get_bfield_cartesian_map) != ALIVE) errorstop("compute_full_map_comsol_symmetries: x>0, y<0 get_bfield != ALIVE\n");
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", -x, -y, z, -bx, -by, bz); //in the original map, x<0, y>0
			}
		}
		for(j=1;j<nb_y;j++) {
			y = ymin + j*ystep; //y from 0 to 2 because y>0 in the original map
			for(k=0;k<nb_x;k++) {
				x = xmin + k*xstep; //because x<0 in the original map
				//x = xmax -k*xstep; //change to this if x>0 in the original map and change fprintf as well!!
				if(get_bfield(x, y, z, &bx, &by, &bz, &(latt.cell[0]), get_bfield_cartesian_map) != ALIVE) errorstop("compute_full_map_comsol_symmetries: x<0, y<0 get_bfield != ALIVE\n");
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", x, y, z, bx, by, bz); //in the original map, x<0, y>0
			}
			for(k=1;k<nb_x;k++) {
				x = xmax - k*xstep;
				if(get_bfield(x, y, z, &bx, &by, &bz, &(latt.cell[0]), get_bfield_cartesian_map) != ALIVE) errorstop("compute_full_map_comsol_symmetries: x>0, y<0 get_bfield != ALIVE\n");
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", -x, y, z, -bx, by, bz); //in the original map, x<0, y>0
			}
		}
	}
	free_latt(&latt);
	fclose(wfile);
}

extern void find_min_max_column_textfile(char *textfile, int nb_col_tot, int col_nb, double *min, double *max, int jump_header)
{
	int i,j,nblines;
	double temp[nb_col_tot];
	FILE *rfile=NULL;
	
	nblines = get_nb_lines_file(textfile);
	
	rfile = fopen(textfile, "r");
	if(rfile==NULL) errorstop("cannot open rfile\n");
	
	*min = 1.e9;
	*max = -1.e9;
	for(i=0;i<jump_header;i++) newline(rfile);
	for(i=0;i<nblines;i++) {
		for(j=0;j<nb_col_tot;j++) {
			fscanf(rfile, "%le", &temp[j]);
			//printf("%le	", temp[j]);
		}
		//printf("\n");
		
		if(sanity_test_number(temp[col_nb-1])==TRUE) {
			if(temp[col_nb-1]<*min) *min = temp[col_nb-1];
			if(temp[col_nb-1]>*max) *max = temp[col_nb-1];
		}
	}
	printf("in %s (%i columns), column %i:\n \t(min, max) = (%le, %le)\n",textfile, nb_col_tot, col_nb, *min,*max);
}

extern void vffa_rect_magnet_create_m_map(char *textout, char *text_lattmap, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz, double halfgap, double long_halflength, double dz)
{
	int i,j,k,flag_out_coil;
	double step_x,step_y,step_z, x,y,z, bx,by,bz,btot,bxdz,bydz,bzdz,btotdz,m_loc;
	//double mx_loc,my_loc,mz_loc;
	
	struct Lattice latt;
	FILE *wfile;
	
	step_x =  comp_step(xmin, xmax, nbx);
	step_y = comp_step(ymin, ymax, nby);
	step_z =  comp_step(zmin, zmax, nbz);
	
	printf("zmin = %le, zmax = %le, step_z = %le\n", zmin, zmin+(nbz-1)*step_z, step_z);
	printf("xmin = %le, xmax = %le, step_x = %le\n", xmin, xmin+(nbx-1)*step_x, step_x);
	printf("ymin = %le, ymax = %le, step_y = %le\n", ymin, ymin+(nby-1)*step_y, step_y);
	
	wfile = fopen(textout, "w");
	fprintf(wfile, "%i	%i	%i\n",nbx,nby,nbz);
	fprintf(wfile, "x [m],y,z,0,0,m_loc [/m]\n");
	load_lattice(&latt, text_lattmap);
	
	if(test_cell_map(&(latt.cell[0]))==NO) errorstop("lattice not a map");
	
	for(i=0;i<nbz;i++) {
		z = zmin + i*step_z;
		for(j=0;j<nbx;j++) {
			x = xmin + j*step_x;
			for(k=0;k<nby;k++) {
				y = ymin + k*step_y;
				if(y>latt.cell[0].boun.ymax/2.-long_halflength && y<latt.cell[0].boun.ymax/2.+long_halflength) { //check if outside of the coils, then m_loc=m_gen=0
					if(x<-halfgap || x>halfgap) flag_out_coil = YES;
					else flag_out_coil = NO;
				}
				else flag_out_coil = NO;
				if(flag_out_coil == YES) m_loc = 0;
				else {
					if(get_bfield(x, y, z, &bx, &by, &bz, &(latt.cell[0]), get_bfield_cartesian_map) != ALIVE) errorstop("compute_m_rect_magnet: get_bfield != ALIVE\n");
					btot = sqrt(bx*bx+by*by+bz*bz);
					if(get_bfield(x, y, z+dz, &bxdz, &bydz, &bzdz, &(latt.cell[0]), get_bfield_cartesian_map) != ALIVE) errorstop("compute_m_rect_magnet: get_bfield dz != ALIVE\n");
					btotdz = sqrt(bxdz*bxdz+bydz*bydz+bzdz*bzdz);
					m_loc = local_m_computation(btotdz, btot, dz);
					//mx_loc = local_m_computation(bxdz, bx, dz);
					//my_loc = local_m_computation(bydz, by, dz);
					//mz_loc = local_m_computation(bzdz, bz, dz);
				}
				fprintf(wfile, "%le	%le	%le	0	0	%le\n", x,y,z,m_loc);
			}
			fprintf(wfile,"\n");
		}
	}
	
	free_latt(&latt);
	fclose(wfile);
}

extern void vffa_str_fodo_create_m_map(char *textout, char *text_lattmap, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz, double halfgap, double long_halflength, double dz)
{
	int i,j,k,flag_out_coil;
	double step_x,step_y,step_z, x,y,z, bx,by,bz,btot,bxdz,bydz,bzdz,btotdz,m_loc;
	//double mx_loc,my_loc,mz_loc;
	
	struct Lattice latt;
	FILE *wfile;
	step_x =  comp_step(xmin, xmax, nbx);
	step_y = comp_step(ymin, ymax, nby);
	step_z =  comp_step(zmin, zmax, nbz);
	
	printf("zmin = %le, zmax = %le, step_z = %le\n", zmin, zmin+(nbz-1)*step_z, step_z);
	printf("xmin = %le, xmax = %le, step_x = %le\n", xmin, xmin+(nbx-1)*step_x, step_x);
	printf("ymin = %le, ymax = %le, step_y = %le\n", ymin, ymin+(nby-1)*step_y, step_y);
	
	printf("Creating file %s...\n", textout);
	wfile = fopen(textout, "w");
	fprintf(wfile, "%i	%i	%i\n",nbx,nby,nbz);
	fprintf(wfile, "x [m],y,z,0,0,m_loc [/m]\n");
	load_lattice(&latt, text_lattmap);
	
	if(test_cell_map(&(latt.cell[0]))==NO) errorstop("lattice not a map");
	
	for(i=0;i<nbz;i++) {
		z = zmin + i*step_z;
		for(j=0;j<nbx;j++) {
			x = xmin + j*step_x;
			for(k=0;k<nby;k++) {
				y = ymin + k*step_y;
				//printf("(%le,%le,%le)\n",x,y,z);
				if(get_bfield(x, y, z, &bx, &by, &bz, &(latt.cell[0]), get_bfield_cartesian_map) != ALIVE) errorstop("vffa_str_fodo_create_m_map: get_bfield != ALIVE\n");
				btot = sqrt(bx*bx+by*by+bz*bz);
				if(get_bfield(x, y, z+dz, &bxdz, &bydz, &bzdz, &(latt.cell[0]), get_bfield_cartesian_map) != ALIVE) errorstop("vffa_str_fodo_create_m_map: get_bfield dz != ALIVE\n");
				btotdz = sqrt(bxdz*bxdz+bydz*bydz+bzdz*bzdz);
				m_loc = local_m_computation(btotdz, btot, dz);
				fprintf(wfile, "%le	%le	%le	0	0	%le\n", x,y,z,m_loc);
			}
			//fprintf(wfile,"\n");
		}
	}
	
	free_latt(&latt);
	fclose(wfile);
}

extern void compute_int_m_track_map(char *outfile, struct Cell *m_cell, struct Lattice *track_latt, struct Particle *part)
{
	int i,nblines;
	double temp[7],s0,s,x,y,z,brho, d1,d2,d3,d4,m, int_m,qu,qv;
	FILE *wfile;
	FILE *rfile=NULL;
	wfile = fopen(outfile, "a");
	
	int_m = 0;
	emptyfile("data/track_for_m_computation.dat");
	part_cross_latt(part, track_latt, "data/track_for_m_computation.dat");
	calc_tune_twiss(part, &qu, &qv, &d1, &d2, &d3, &d4, 1.e-5, 1.e-5, 1.e-5, 1.e-5, track_latt, part_cross_latt, YES, NULL,1);//flag: 0=2d, 1=4d decoupled, anything else=4d undecoupled.
	
	nblines = get_nb_lines_file("data/track_for_m_computation.dat");
	rfile = fopen("data/track_for_m_computation.dat", "r");
	if(rfile==NULL) errorstop("cannot open file");
	fscanf(rfile, "%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le", &s0, &x, &y, &z, &temp[0], &temp[1], &temp[2], &temp[3], &temp[4], &temp[5], &brho);
	newline(rfile);
	printf("%le \t %le \t %le \t %le\n", s0, x, y, z);
	for(i=1;i<nblines-1;i++) {
		fscanf(rfile, "%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le", &s, &x, &y, &z, &temp[0], &temp[1], &temp[2], &temp[3], &temp[4], &temp[5], &brho);
		newline(rfile);
		//printf("%le \t %le \t %le \t %le\n", s, x, y, z);
		
		if(get_bfield(x, y, z, &d1, &d2, &m, m_cell, get_bfield_cartesian_map) != ALIVE) errorstop("compute_int_m_track_map: get_bfield != ALIVE\n");
		int_m += m*(s-s0);
		printf("m=%lf, int_m=%lf, s=%lf, s0=%lf, s-s0=%lf\n",m,int_m,s,s0,s-s0);
		s0 = s;
	}
	fclose(rfile);
	int_m/=s0;
	fprintf(wfile, "%le	%le	%le	%le\n", brho, int_m, qu,qv);
	fclose(wfile);
}

extern int adjust_brho_part_vffa(double z_target, double m_map, double z_step, struct Particle *reference, struct Lattice *latt, int doyouprintf)
{
	int i;
	double brho_step, ini_brho=reference->brho, deltaz, eps_z=1.e-4;
	double e_tot,beta_lor, gamma_lor;
	struct Particle part;
	part= *reference;
	
	for(i=0;i<100;i++) {
		if(doyouprintf==YES) printf("i=%i,z_part=%lf, target %lf\n",i,part.z,z_target);
		if(put_on_co_nelmin(latt, &part, 1.e-20, 1.e-4, 1.e-4, 1.e-4, 1.e-4, doyouprintf) == FALSE) {
			printf("could not find co_nelmin, z=%le\n", part.z);
			reference->s=0;
			return FALSE;
		}
		if(find_closed_orbite_xxp_zzp(&part, &part.x, &part.ux, &part.uy, &part.z, &part.uz, 1.e-10, latt, doyouprintf) == FALSE) {
			printf("could not find co_xxp_zzp, z=%le\n", part.z);
			reference->s=0;
			return FALSE;
		}
		
		deltaz=part.z-z_target;
		if(fabs(deltaz)<eps_z) {
			*reference = part;
			get_ebg_part(&e_tot, &beta_lor, &gamma_lor, reference);
			printf("Momentum adjusted, at z=%lf, Ekin=%le [eV]\n",reference->z,(gamma_lor - 1)*reference->m0*CLIGHT2/(UNITCHARGE));
			reference->s=0;
			return TRUE;
		}
		if(fabs(deltaz)>1.5*z_step) part.brho*=(1-sign(deltaz)*z_step*m_map);
		else part.brho*=(1-deltaz*m_map);
	}
	printf("failed to achieve z_target=%lf, delta z = %le, increase number of tries or decrease precision\n",z_target,deltaz);
	reference->s=0;
	return FALSE;
}

//extern void compute_int_m_track_map_all_mom(char *outfile, struct Cell *m_cell, struct Lattice *track_latt, struct Particle *part, double mom_min, double mom_max)

extern void plot_with_without_iron(char *trackout, char *heatmapfile, char *tuneout, char *epsheatmap, char *eps_long_vert, char *eps_tune, double p0)
{
	char xrange[100], yrange[100], pp0[100];
	double min,max;
	
	find_boundaries_heated_map_file2(heatmapfile, xrange, yrange, &min, &max);
	sprintf(pp0,"($1/(%le))",p0);
	easyplot3dmap_plusfile(trackout, "($3)", "($2)", "lines lc 7 lw 2", heatmapfile, "long [m]", "hor [m]", "B vert [T]", xrange, yrange, epsheatmap, "size ratio -1",  min,max);
	
	easyplot2(tuneout, tuneout, pp0, "2", pp0, "3", "lines lt 1 lw 4 lc 0", "lines lt 2 lw 4 lc 0", NULL, NULL, "P/P0", "{/Symbol n}_u_,_v []", NULL, NULL, eps_tune, NULL);
	easyplot(trackout, "($3)", "($4)", "lines lc 7 lw 2", "long [m]", "vert [m]", xrange, NULL, eps_long_vert, "mxtics 5\nset mytics 5");
	
}

//adjust closed orbit on r_co by changing B0
extern int adjust_b0_ffag_spi(struct Lattice *latt, struct Particle *part, double r_co, double eps_clo, double eps_r0)
{
	int i, j, k;
	double rnew;
	printf("adjust_B0:\n\n");
	for(i = 0; i < 30; i++) {
		printf("i=%i\n",i);
		if(find_closed_orbite_xxp(part, &(part->x), &(part->ux), &(part->uy), eps_clo, latt, NO) == TRUE) {
			rnew = part->x;
			//part->x = rnew;
			//tune_calc_matrix(part, &nux, &nuz, &betax, &alphax, &betaz, &alphaz, 1.e-4, 1.e-5, 1.e-4, 1.e-5, latt, part_cross_latt, NO, NULL);
			if (fabs(rnew - r_co) < eps_r0) {
				//printf("\n\tr0adjust = %.8f [m]\n", latt->cell[0].mpara[0][1]);
				printf("Adjusted B0 at r0=%lf: (%le, %le) [T]\n", latt->cell[0].mpara[0][1], latt->cell[0].mpara[0][2], latt->cell[0].mpara[1][2]);
				return TRUE;
			}
			for(k=0;k<latt->nbcell;k++) {
				for(j = 0; j < latt->cell[k].nbcomp; j++) {
					latt->cell[k].mpara[j][2] *= pow((rnew-r_co + latt->cell[k].mpara[j][1])/latt->cell[k].mpara[j][1], latt->cell[k].mpara[j][3]);
				}
			}
			part->x = r_co;
			//printf("rnew = %lf, r_co = %lf\n", rnew, r_co);
			//for(j = 0; j < latt->cell[0].nbcomp; j++) printf("latt->cell[0].mpara[%i][1] = %lf\n", j, latt->cell[0].mpara[j][1]); 
		}
		else {
			printf("\n \nproblem in adjust_r0, closed orbit not found!\n");
			return FALSE;
		}
	}
	printf("\n \nin adjust_r0 precision not achieved: eps = %le, increase the number of turns or decrease demanded precision\n \n", fabs(rnew - r_co));
	printf("\n\n\t\t\tB0adjust = %lf [m]\n\n\n", latt->cell[0].mpara[0][2]);
	return FALSE;
}

//adjust closed orbit on average r_co by changing B0 (1 cell only!)
extern int adjust_b0_ffag_spi_rav(struct Lattice *latt, struct Particle *part, double r_co_av, double eps_clo, double eps_r0)
{
	int i, j, k, nb_iterations=50;
	double rnew, rnew_av;
	printf("adjust_B0:\n\n");
	for(i = 0; i < nb_iterations; i++) {
		//printf("i=%i\n",i);
		if(find_closed_orbite_xxp(part, &rnew, &(part->ux), &(part->uy), eps_clo, latt, NO) == TRUE) {
			emptyfile("data/r_av_computation.dat");
			part->x = rnew;
			part->s = 0;
			part_cross_latt(part, latt, "data/r_av_computation.dat");
			rnew_av = compute_av_radius_from_trackout_1cell("data/r_av_computation.dat", NO);
			printf("r_co(s0)=%lf, r_co_av = %lf\n", rnew, rnew_av);
			//tune_calc_matrix(part, &nux, &nuz, &betax, &alphax, &betaz, &alphaz, 1.e-4, 1.e-5, 1.e-4, 1.e-5, latt, part_cross_latt, NO, NULL);
			if (fabs(rnew_av - r_co_av) < eps_r0) {
				//printf("\n\tr0adjust = %.8f [m]\n", latt->cell[0].mpara[0][1]);
				printf("iteration %i, adjusted B0 at r0=%lf: (%le, %le) [T]\n", i, latt->cell[0].mpara[0][1], latt->cell[0].mpara[0][2], latt->cell[0].mpara[1][2]);
				part->s = 0;
				return TRUE;
			}
			for(k=0;k<latt->nbcell;k++) {
				for(j = 0; j < latt->cell[k].nbcomp; j++) {
					latt->cell[k].mpara[j][2] *= pow(((rnew_av-r_co_av) + latt->cell[k].mpara[j][1])/latt->cell[k].mpara[j][1], latt->cell[k].mpara[j][3]);
					if(j==0) printf("new B0=%lf T\n", latt->cell[k].mpara[j][2]);
				}
			}
			part->x = rnew - (rnew_av-r_co_av);
			//part->x = r_co_av;
			//printf("rnew = %lf, r_co = %lf\n", rnew, r_co);
			//for(j = 0; j < latt->cell[0].nbcomp; j++) printf("latt->cell[0].mpara[%i][1] = %lf\n", j, latt->cell[0].mpara[j][1]); 
		}
		else {
			printf("\n \niteration %i, problem in adjust_r0, closed orbit not found!\n", i);
			part->s = 0;
			return FALSE;
		}
	}
	printf("\n \nin adjust_r0 precision not achieved: eps = %le, increase the number of iterations (now %i) or decrease demanded precision\n \n", fabs(rnew - r_co_av), nb_iterations);
	printf("\n\n\t\t\tB0adjust = %lf [m]\n\n\n", latt->cell[0].mpara[0][2]);
	part->s = 0;
	return FALSE;
}

extern void write_function_textfile(char *textfile, double xmin, double xmax, int nbx)
{
	int i;
	double y, x, xstep = comp_step(xmin, xmax, nbx);
	double lambda = 0.0395165;
	double c0=0.1455;
	double c1=2.267;
	double c2=-0.6395;
	double c3=1.1558;
	
	FILE *wfile;
	wfile = fopen(textfile, "w");
	
	for(i=0;i<nbx;i++) {
		x=xmin+i*xstep;
		//y = 1/(1+exp(0.1455+2.267*x/0.0395165-0.6395*x/0.0395165*x/0.0395165+1.1558*x/0.0395165*x/0.0395165*x/0.0395165));
		y = 1/(1+exp(c0+c1*x/lambda+c2*x/lambda*x/lambda+c3*x/lambda*x/lambda*x/lambda));
		fprintf(wfile, "%le	%le\n", x, y);
	}
		
	fclose(wfile);
}

//geometric emittance, not normalised emittanc
extern void acceptance_2d_emit_auto(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double emitmin, double *emitmax, double emitstep)
{
	int i, n, doyouprint=YES;
	double qx, qz, beta_x, alpha_x, beta_z, alpha_z, gamma_x, gamma_z, emit;
	double ax, az, xpart[nbpass], zpart[nbpass], ux[nbpass], uy[nbpass], uz[nbpass], x, z, xprime, zprime, php_ref, atanrefx, atanrefz;
	FILE *wfile;
	struct Particle test_part = *reference;
	test_part.hat = -2;
	test_part.status = ALIVE;
	test_part.s = 0;
	
	if(outfilename != NULL) wfile = fopen(outfilename, "w");
	php_ref = sqrt(reference->ux*reference->ux + reference->uy*reference->uy);
	atanrefx = atan_ratio((reference->ux), (reference->uy));
	atanrefz = atan_ratio(reference->uz, php_ref);
	
	calc_tune_twiss(&test_part, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-5, 1.e-5, 1.e-5, 1.e-5, latt, part_cross_latt, YES, NULL,0);//flag: 0=2d, 1=4d decoupled, anything else=4d undecoupled.
	gamma_x = (1+alpha_x*alpha_x)/beta_x;
	gamma_z = (1+alpha_z*alpha_z)/beta_z;
	printf("gamma_x = %lf, gamma_z = %lf\n",gamma_x,gamma_z);
	emit = emitmin;
	i = 0;
	*emitmax = 0.;
	for(;;) {
		emit = emitmin+i*emitstep;
		ax = sqrt(emit/gamma_x);
		az = sqrt(emit/gamma_z);
		//CLRSCR();
		if(doyouprint==YES) printf("unnorm. emit = %le [m], ampx = %lf, ampz = %lf\t", emit, ax, az);
		//fflush(stdout);
		if(track_n_turns_amp(xpart, zpart, ux, uy, uz, &test_part, latt, nbpass, ax, az) != ALIVE) {
			if(doyouprint==YES) printf("NO!\n");
			break;
		}
		else {
			if(doyouprint==YES) printf("OK!\n");
			//*axmax = ax;
			if(outfilename != NULL) {
				for(n=0;n<nbpass;n++) {
					comp_phase_space_coord(&x, &xprime, &z, &zprime, xpart[n], zpart[n], ux[n], uy[n], uz[n], test_part.x, atanrefx, test_part.z, atanrefz);
					fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le\n", x, xprime, z, zprime, emit, ax, az);
				}
			}
		}
		i++;
	}
	if(i==0) *emitmax = 0.;
	else *emitmax = emit-emitstep;
	if(outfilename != NULL) fclose(wfile);
	if(doyouprint==YES) printf("final: emit = %le\n",*emitmax);
}



/*extern void acceptance_2d_turns(char *textout, struct Lattice *latt, struct Particle *part)
{
	int turn;
	double qx, qz, beta_x, alpha_x, beta_z, alpha_z, gamma_x, gamma_z, e_tot, beta_lor, gamma_lor, x_offset, z_offset, emit2d, axmax, emitx, azmax, emitz;
	FILE *wfile;
	
	calc_tune_twiss(part, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-5, 1.e-5, 1.e-5, 1.e-5, latt, part_cross_latt, NO, NULL,0);//flag: 0=2d, 1=4d decoupled, anything else=4d undecoupled.
	gamma_x = (1+alpha_x*alpha_x)/beta_x;
	gamma_z = (1+alpha_z*alpha_z)/beta_z;
	get_ebg_part(&e_tot, &beta_lor, &gamma_lor, part);
	
	x_offset = sqrt(10.e-6/gamma_x/beta_lor/gamma_lor);
	z_offset = sqrt(10.e-6/gamma_z/beta_lor/gamma_lor);
	
	wfile = fopen(textout, "w");
	
	turn = 150;
	emit2d = 0;
	emitx = 0;
	emitz = 0;
	printf("150 cells\n");
	acceptance_2d_emit_auto(NULL, part, latt, turn, 10.e-6, &emit2d, 5.e-5);
	emit2d *= beta_lor*gamma_lor;
	acceptancex_auto(NULL, part, latt, turn, 0.005, &axmax, 0.005, z_offset);
	emitx = axmax*axmax*gamma_x*beta_lor*gamma_lor;
	acceptancez_auto(NULL, part, latt, turn, 0.005, &azmax, 0.002, x_offset); 
	emitz = azmax*azmax*gamma_z*beta_lor*gamma_lor;
	fprintf(wfile, "%i	%le	%le	%le\n", turn, emit2d, emitx, emitz);
	
	turn = 750;
	emit2d = 0;
	emitx = 0;
	emitz = 0;
	printf("750 cells\n");
	acceptance_2d_emit_auto(NULL, part, latt, turn, 10.e-6, &emit2d, 5.e-5);
	emit2d *= beta_lor*gamma_lor;
	acceptancex_auto(NULL, part, latt, turn, 0.005, &axmax, 0.005, z_offset);
	emitx = axmax*axmax*gamma_x*beta_lor*gamma_lor;
	acceptancez_auto(NULL, part, latt, turn, 0.005, &azmax, 0.002, x_offset); 
	emitz = azmax*azmax*gamma_z*beta_lor*gamma_lor;
	fprintf(wfile, "%i	%le	%le	%le\n", turn, emit2d, emitx, emitz);
	
	turn = 1500;
	emit2d = 0;
	emitx = 0;
	emitz = 0;
	printf("1500 cells\n");
	acceptance_2d_emit_auto(NULL, part, latt, turn, 10.e-6, &emit2d, 5.e-5);
	emit2d *= beta_lor*gamma_lor;
	acceptancex_auto(NULL, part, latt, turn, 0.005, &axmax, 0.005, z_offset);
	emitx = axmax*axmax*gamma_x*beta_lor*gamma_lor;
	acceptancez_auto(NULL, part, latt, turn, 0.005, &azmax, 0.002, x_offset); 
	emitz = azmax*azmax*gamma_z*beta_lor*gamma_lor;
	fprintf(wfile, "%i	%le	%le	%le\n", turn, emit2d, emitx, emitz);
	
	turn = 7500;
	emit2d = 0;
	emitx = 0;
	emitz = 0;
	printf("7500 cells\n");
	acceptance_2d_emit_auto(NULL, part, latt, turn, 10.e-6, &emit2d, 5.e-5);
	emit2d *= beta_lor*gamma_lor;
	acceptancex_auto(NULL, part, latt, turn, 0.005, &axmax, 0.005, z_offset);
	emitx = axmax*axmax*gamma_x*beta_lor*gamma_lor;
	acceptancez_auto(NULL, part, latt, turn, 0.005, &azmax, 0.002, x_offset); 
	emitz = azmax*azmax*gamma_z*beta_lor*gamma_lor;
	fprintf(wfile, "%i	%le	%le	%le\n", turn, emit2d, emitx, emitz);
	
	turn = 15000;
	emit2d = 0;
	emitx = 0;
	emitz = 0;
	printf("15000 cells\n");
	acceptance_2d_emit_auto(NULL, part, latt, turn, 10.e-6, &emit2d, 5.e-5);
	emit2d *= beta_lor*gamma_lor;
	acceptancex_auto(NULL, part, latt, turn, 0.005, &axmax, 0.005, z_offset);
	emitx = axmax*axmax*gamma_x*beta_lor*gamma_lor;
	acceptancez_auto(NULL, part, latt, turn, 0.005, &azmax, 0.002, x_offset); 
	emitz = azmax*azmax*gamma_z*beta_lor*gamma_lor;
	fprintf(wfile, "%i	%le	%le	%le\n", turn, emit2d, emitx, emitz);
	
	turn = 45000;
	emit2d = 0;
	emitx = 0;
	emitz = 0;
	printf("45000 cells\n");
	acceptance_2d_emit_auto(NULL, part, latt, turn, 10.e-6, &emit2d, 5.e-5);
	emit2d *= beta_lor*gamma_lor;
	acceptancex_auto(NULL, part, latt, turn, 0.005, &axmax, 0.005, z_offset);
	emitx = axmax*axmax*gamma_x*beta_lor*gamma_lor;
	acceptancez_auto(NULL, part, latt, turn, 0.005, &azmax, 0.002, x_offset); 
	emitz = azmax*azmax*gamma_z*beta_lor*gamma_lor;
	fprintf(wfile, "%i	%le	%le	%le\n", turn, emit2d, emitx, emitz);
	
	turn = 75000;
	emit2d = 0;
	emitx = 0;
	emitz = 0;
	printf("75000 cells\n");
	acceptance_2d_emit_auto(NULL, part, latt, turn, 10.e-6, &emit2d, 5.e-5);
	emit2d *= beta_lor*gamma_lor;
	acceptancex_auto(NULL, part, latt, turn, 0.005, &axmax, 0.005, z_offset);
	emitx = axmax*axmax*gamma_x*beta_lor*gamma_lor;
	acceptancez_auto(NULL, part, latt, turn, 0.005, &azmax, 0.002, x_offset); 
	emitz = azmax*azmax*gamma_z*beta_lor*gamma_lor;
	fprintf(wfile, "%i	%le	%le	%le\n", turn, emit2d, emitx, emitz);
	
	turn = 150000;
	emit2d = 0;
	emitx = 0;
	emitz = 0;
	printf("150000 cells\n");
	acceptance_2d_emit_auto(NULL, part, latt, turn, 10.e-6, &emit2d, 5.e-5);
	emit2d *= beta_lor*gamma_lor;
	acceptancex_auto(NULL, part, latt, turn, 0.005, &axmax, 0.005, z_offset);
	emitx = axmax*axmax*gamma_x*beta_lor*gamma_lor;
	acceptancez_auto(NULL, part, latt, turn, 0.005, &azmax, 0.002, x_offset); 
	emitz = azmax*azmax*gamma_z*beta_lor*gamma_lor;
	fprintf(wfile, "%i	%le	%le	%le\n", turn, emit2d, emitx, emitz);
	
	turn = 300000;
	emit2d = 0;
	emitx = 0;
	emitz = 0;
	printf("300000 cells\n");
	acceptance_2d_emit_auto(NULL, part, latt, turn, 10.e-6, &emit2d, 5.e-5);
	emit2d *= beta_lor*gamma_lor;
	acceptancex_auto(NULL, part, latt, turn, 0.005, &axmax, 0.005, z_offset);
	emitx = axmax*axmax*gamma_x*beta_lor*gamma_lor;
	acceptancez_auto(NULL, part, latt, turn, 0.005, &azmax, 0.002, x_offset); 
	emitz = azmax*azmax*gamma_z*beta_lor*gamma_lor;
	fprintf(wfile, "%i	%le	%le	%le\n", turn, emit2d, emitx, emitz);
	
	fclose(wfile);
}//*/

extern void acceptance_2d_turns(char *textout, struct Lattice *latt, struct Particle *part, int turn)
{
	double qx, qz, beta_x, alpha_x, beta_z, alpha_z, gamma_x, gamma_z, e_tot, beta_lor, gamma_lor, x_offset, z_offset, emit2d, axmax, emitx, azmax, emitz;
	FILE *wfile;
	
	calc_tune_twiss(part, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-5, 1.e-5, 1.e-5, 1.e-5, latt, part_cross_latt, NO, NULL,0);//flag: 0=2d, 1=4d decoupled, anything else=4d undecoupled.
	gamma_x = (1+alpha_x*alpha_x)/beta_x;
	gamma_z = (1+alpha_z*alpha_z)/beta_z;
	get_ebg_part(&e_tot, &beta_lor, &gamma_lor, part);
	
	x_offset = sqrt(10.e-6/gamma_x/beta_lor/gamma_lor);
	z_offset = sqrt(10.e-6/gamma_z/beta_lor/gamma_lor);
	
	wfile = fopen(textout, "a");
	
	emit2d = 0;
	emitx = 0;
	emitz = 0;
	printf("%i lattices\n",turn);
	acceptance_2d_emit_auto(NULL, part, latt, turn, 20.e-6, &emit2d, 1.e-5);
	emit2d *= beta_lor*gamma_lor;
	acceptancex_auto(NULL, part, latt, turn, 0.005, &axmax, 0.001, z_offset);
	emitx = axmax*axmax*gamma_x*beta_lor*gamma_lor;
	acceptancez_auto(NULL, part, latt, turn, 0.005, &azmax, 0.001, x_offset); 
	emitz = azmax*azmax*gamma_z*beta_lor*gamma_lor;
	fprintf(wfile, "%i	%le	%le	%le\n", turn, emit2d, emitx, emitz);
	fclose(wfile);
}

extern void read_change_raccam_fieldmap(char *file_opera, char *fileout)
{
	int i, nblines;
	double x,y,z,bx,by,bz, r, th;
	FILE *rfile=NULL;
	FILE *wfile;
	
	nblines = get_nb_lines_file(file_opera);
	wfile = fopen(fileout, "w");
	rfile = fopen(file_opera, "r");
	if(rfile==NULL) errorstop("cannot open file");
	
	for(i=0;i<nblines;i++) {
		fscanf(rfile, "%le	%le	%le	%le	%le	%le", &x,&y,&z,&bx,&by,&bz);
		printf( "%le	%le	%le	%le	%le	%le\n", x,y,z,bx,by,bz);
		r = sqrt(x*x+y*y);
		th = atan_ratio(y,x);
		printf("%le, %le\n",r,th);
		fprintf(wfile, "%le	%le	%le\n", r,th,bz);
	}
	
	fclose(rfile);
	fclose(wfile);
}

extern void emittance_2d_expansion_order(char *textout, struct Lattice *latt, struct Particle *reference, double order)
{
	int i, j, turn;
	double qx, qz, beta_x, alpha_x, beta_z, alpha_z, gamma_x, gamma_z, e_tot, beta_lor, gamma_lor, x_offset, z_offset, emit2d, axmax, emitx, azmax, emitz;
	double order_ori = latt->cell[0].efben[0][2];
	FILE *wfile;
	struct Particle test_part;
	
	test_part = *reference;
	
	if(find_closed_orbite_xxp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), 1.e-12, latt, YES) == FALSE) errorstop("closed orbit not found\n");
	
	calc_tune_twiss(&test_part, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-5, 1.e-5, 1.e-5, 1.e-5, latt, part_cross_latt, NO, NULL,0);//flag: 0=2d, 1=4d decoupled, anything else=4d undecoupled.
	gamma_x = (1+alpha_x*alpha_x)/beta_x;
	gamma_z = (1+alpha_z*alpha_z)/beta_z;
	get_ebg_part(&e_tot, &beta_lor, &gamma_lor, &test_part);
	
	x_offset = sqrt(10.e-6/gamma_x/beta_lor/gamma_lor);
	z_offset = sqrt(10.e-6/gamma_z/beta_lor/gamma_lor);
	
	wfile = fopen(textout, "a");
	
	//turn = 150000;
	turn = 1000;
	//turn = 10000;
	
	for(j=0;j<latt->nbcell;j++) {
		for(i=0;i<latt->cell[j].nbcomp;i++) latt->cell[j].efben[i][2] = order;
	}
	emit2d = 0;
	emitx = 0;
	emitz = 0;
	printf("order: %.0f\n",order);
	acceptance_2d_emit_auto(NULL, &test_part, latt, turn, 20.e-6, &emit2d, 1.e-5);
	emit2d *= beta_lor*gamma_lor;
	acceptancex_auto(NULL, &test_part, latt, turn, 0.005, &axmax, 0.001, z_offset);
	emitx = axmax*axmax*gamma_x*beta_lor*gamma_lor;
	acceptancez_auto(NULL, &test_part, latt, turn, 0.005, &azmax, 0.001, x_offset); 
	emitz = azmax*azmax*gamma_z*beta_lor*gamma_lor;
	fprintf(wfile, "%lf	%le	%le	%le\n", order, emit2d, emitx, emitz);
	
	for(i=0;i<latt->cell[0].nbcomp;i++) latt->cell[0].efben[i][2] = order_ori;
	fclose(wfile);
}

extern void emittance_2d_stepsize(char *textout, struct Lattice *latt, struct Particle *reference, double stepsize)
{
	int j, turn;
	double qx, qz, beta_x, alpha_x, beta_z, alpha_z, gamma_x, gamma_z, e_tot, beta_lor, gamma_lor, x_offset, z_offset, emit2d, axmax, emitx, azmax, emitz;
	double stepsize_ori = latt->cell[0].stepsize;
	FILE *wfile;
	struct Particle test_part;
	
	test_part = *reference;
	for(j=0;j<latt->nbcell;j++) {
		latt->cell[j].stepsize = stepsize;
	}
	
	if(find_closed_orbite_xxp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), 1.e-12, latt, YES) == FALSE) errorstop("closed orbit not found\n");
	
	calc_tune_twiss(&test_part, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-5, 1.e-5, 1.e-5, 1.e-5, latt, part_cross_latt, YES, NULL,0);//flag: 0=2d, 1=4d decoupled, anything else=4d undecoupled.
	gamma_x = (1+alpha_x*alpha_x)/beta_x;
	gamma_z = (1+alpha_z*alpha_z)/beta_z;
	get_ebg_part(&e_tot, &beta_lor, &gamma_lor, &test_part);
	
	x_offset = sqrt(10.e-6/gamma_x/beta_lor/gamma_lor);
	z_offset = sqrt(10.e-6/gamma_z/beta_lor/gamma_lor);
	
	wfile = fopen(textout, "a");
	
	turn = 1000;
	//turn = 10000;
	
	emit2d = 0;
	emitx = 0;
	emitz = 0;
	printf("step size %lf m\n", stepsize);
	acceptance_2d_emit_auto(NULL, &test_part, latt, turn, 20.e-6, &emit2d, 1.e-5);
	emit2d *= beta_lor*gamma_lor;
	acceptancex_auto(NULL, &test_part, latt, turn, 0.005, &axmax, 0.001, z_offset);
	emitx = axmax*axmax*gamma_x*beta_lor*gamma_lor;
	acceptancez_auto(NULL, &test_part, latt, turn, 0.005, &azmax, 0.001, x_offset); 
	emitz = azmax*azmax*gamma_z*beta_lor*gamma_lor;
	fprintf(wfile, "%lf	%le	%le	%le\n", stepsize, emit2d, emitx, emitz);
	
	for(j=0;j<latt->nbcell;j++) latt->cell[j].stepsize = stepsize_ori;
	fclose(wfile);
}

extern int iterative_calc_tune_twiss(struct Particle *reference, double *qx, double *qz, double *betax, double *alphax, double *betaz, double *alphaz, double amp_x, double amp_z, 
	struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, char *outfile)
{
	int i=0;
	double gammax,gammaz,emitx,emitz, amp_xprime=amp_x, amp_zprime=amp_z, eps=1.e-9;
	
	for(;;) {
		calc_tune_twiss(reference, qx,qz,betax,alphax,betaz,alphaz, amp_x, amp_xprime, amp_z, amp_zprime, latt, transport_part, NO, NULL, 0); //does not work for coupled case!!
		gammax = (1+(*alphax)*(*alphax))/(*betax);
		gammaz = (1+(*alphaz)*(*alphaz))/(*betaz);
		emitx = amp_x*amp_x*gammax;
		emitz = amp_z*amp_z*gammaz;
		printf("emitx %le, amp_xprime %le, emitz %le, amp_zprime %le\n", emitx, amp_xprime, emitz, amp_zprime);
		if(fabs(amp_xprime-sqrt(emitx/(*betax)))<eps && fabs(amp_zprime-sqrt(emitz/(*betaz)))<eps) break;
		amp_xprime=sqrt(emitx/(*betax));
		amp_zprime=sqrt(emitz/(*betaz));
		i++;
		if(i==20) {
			printf("not found!\n");
			return FALSE;
		} 
	}
	printf("found!\n");
	calc_tune_twiss(reference, qx,qz,betax,alphax,betaz,alphaz, amp_x, amp_xprime, amp_z, amp_zprime, latt, transport_part, doyouprintf, outfile, 0);
	return TRUE;
}

extern double compute_av_radius_from_trackout_1cell(char *trackout, int doyouprintf)
{
	int i, nblines;
	double s, x, y, z, bx, by, bz, ux, uy, uz, brho, t, r, rtot=0, s0=0;
	FILE *rfile = NULL;
	
	nblines = get_nb_lines_file(trackout);
	
	rfile = fopen(trackout, "r");
	if(rfile == NULL) errorstop("!!!ERROR in compute_av_radius_from_trackout_1cell: cannot open inputfile!!!");
	
	for(i=0;i<nblines-1;i++) {
		fscanf(rfile, "%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le", &s, &x, &y, &z, &bx, &by, &bz, &ux, &uy, &uz, &brho, &t, &r);
		//printf("r=%le, s=%le\n", r, s);
		rtot += r*(s-s0); 
		newline(rfile);
		s0 = s;
	}
	fclose(rfile);
	rtot = rtot/s;
	//rtot=rtot/nblines;
	if(doyouprintf==YES) printf("average radius %le\n", rtot);
	return rtot;
}

extern void create_textfile_pynaff(char *textfile, int nbpass, double amp_x, double amp_z, struct Particle *reference, struct Lattice *latt)
{
	int i;
	struct Particle test_part;
	FILE *wfile;
	
	test_part = *reference;
	test_part.hat = -2;
	test_part.status = ALIVE;
	test_part.x += amp_x;
	test_part.z += amp_z;
	
	if(textfile!=NULL) wfile = fopen(textfile, "w");
	
	for(i = 0; i < nbpass; i++) {
		part_cross_latt(&test_part, latt,NULL);
		fprintf(wfile, "%le	%le\n", test_part.x, test_part.z);
		fprintf(wfile, "0	0\n");
	}
	fclose(wfile);
}

extern void write_get_bfield(char *textfile, double x, double y, double z, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), int doyoupolar)
{
	double bx, by, bz;
	double br,bth, r, th;
	FILE *wfile;
	
	wfile = fopen(textfile, "a");
	
	get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp);
	printf("at (%lf, %lf, %lf), field (%le, %le, %le)\n",x,y,z,bx,by,bz);
	if(doyoupolar == YES) {
		r= sqrt(x*x + y*y); 
		th = atan_ratio(y, x);
		br = bx*cos(th) + by*sin(th);
		bth = by*cos(th) - bx*sin(th);
		fprintf(wfile, "%le	%le	%le	%le	%le	%le\n",r,th*180./PI,z, br, bth, bz);
	}
	else {
		fprintf(wfile, "%le	%le	%le	%le	%le	%le\n",x,y,z, bx, by, bz);
	}
	fclose(wfile);
}

extern void emittance_2d_brho(char *textout, struct Lattice *latt, struct Particle *reference, double brho)
{
	int i, turn;
	double qx, qz, beta_x, alpha_x, beta_z, alpha_z, gamma_x, gamma_z, e_tot, beta_lor, gamma_lor, x_offset, z_offset, emit2d, axmax, emitx, azmax, emitz, factor, brho_step;
	FILE *wfile;
	struct Particle test_part;
	
	test_part = *reference;
	test_part.hat = -2;
	brho_step = (brho-reference->brho)/10.;
	if(reference->x*brho_step/(reference->brho*(1+latt->cell[0].mpara[0][3]))<0.01) {
		
	}
	
	//printf("brhostep=%lf\n",brho_step);
	for(i=0;i<10;i++) {
		test_part.brho = reference->brho+brho_step*i;
		//printf("brho = %lf\n", test_part.brho);
		//if(i!=0) test_part.x *= pow((reference->brho+brho_step*i)/(reference->brho+brho_step*(i-1)),latt->cell[0].mpara[0][3]);
		if(find_closed_orbite_xxp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), 1.e-12, latt, NO) == FALSE) errorstop("closed orbit not found\n");
	}
		
	//calc_tune_twiss(&test_part, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-5, 1.e-5, 1.e-5, 1.e-5, latt, part_cross_latt, YES, NULL,0);//flag: 0=2d, 1=4d decoupled, anything else=4d undecoupled.
	//gamma_x = (1+alpha_x*alpha_x)/beta_x;
	//gamma_z = (1+alpha_z*alpha_z)/beta_z;
	//get_ebg_part(&e_tot, &beta_lor, &gamma_lor, &test_part);
	
	//x_offset = sqrt(10.e-6/gamma_x/beta_lor/gamma_lor);
	//z_offset = sqrt(10.e-6/gamma_z/beta_lor/gamma_lor);
	
	turn = 1000;
	
	wfile = fopen(textout, "a");
	printf("brho = %lf\n", brho);
	//acceptance_2d_emit_auto(NULL, &test_part, latt, turn, 20.e-6, &emit2d, 1.e-5);
	acceptance_2d_emit_auto(NULL, &test_part, latt, turn, 3.00000e-04, &emit2d, 5.e-6);
	//emit2d *= beta_lor*gamma_lor;
	//acceptancex_auto(NULL, &test_part, latt, turn, 0.005, &axmax, 0.001, z_offset);
	//emitx = axmax*axmax*gamma_x*beta_lor*gamma_lor;
	//acceptancez_auto(NULL, &test_part, latt, turn, 0.005, &azmax, 0.001, x_offset); 
	//emitz = azmax*azmax*gamma_z*beta_lor*gamma_lor;
	//fprintf(wfile, "%lf	%le	%le	%le	%le\n", brho, emit2d, emitx, emitz, factor);
	fprintf(wfile, "%lf	%le	%le\n", brho, emit2d, test_part.x/latt->cell[0].mpara[0][1]);
	
	fclose(wfile);
}

extern void maxwell_heatmap(char *txt_prefix, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double z, double dx, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	char name[300];
	int i,j;
	double xstep, ystep, x, y, divb, curlbx, curlby, curlbz;
	FILE *wfilediv, *wfilex, *wfiley, *wfilez;
	
	sprintf(name,"%s_div.dat", txt_prefix);
	wfilediv = fopen(name, "w");
	sprintf(name,"%s_curlx.dat", txt_prefix);
	wfilex = fopen(name, "w");
	sprintf(name,"%s_curly.dat", txt_prefix);
	wfiley = fopen(name, "w");
	sprintf(name,"%s_curlz.dat", txt_prefix);
	wfilez = fopen(name, "w");
	xstep = comp_step(xmin, xmax, nbx);
	ystep = comp_step(ymin, ymax, nby);
	
	for(i=0;i<nbx;i++) {
		x = xmin+i*xstep;
		for(j=0;j<nby;j++) {
			y = ymin + j*ystep;
			//maxwell_test_txtfile(textfile, x, y, z, dx, dx, dx, cell, add_contribution_comp);
			if(maxwell_test(x, y, z, dx, dx, dx, cell, add_contribution_comp, &divb, &curlbx, &curlby, &curlbz,NO)!=TRUE) {
				divb = 0;
				curlbx = 0;
				curlby = 0;
				curlbz = 0;
			}
			fprintf(wfilediv, "%le	%le	%le\n", y, x, divb);
			fprintf(wfilex, "%le	%le	%le\n", y, x, curlbx);
			fprintf(wfiley, "%le	%le	%le\n", y, x, curlby);
			fprintf(wfilez, "%le	%le	%le\n", y, x, curlbz);
		}
		fprintf(wfilediv, "\n");
		fprintf(wfilex, "\n");
		fprintf(wfiley, "\n");
		fprintf(wfilez, "\n");
	}
	fclose(wfilediv);
	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
}


extern void maxwell_heatmap_wrapper(char *txt_prefix, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double z, double dx, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int)) 
{
	char name[300];
	char nameeps[300];
	char xrange[100], yrange[100];
	double min, max;
	
	maxwell_heatmap(txt_prefix, xmin, xmax, nbx, ymin, ymax, nby, z, dx, cell, add_contribution_comp);
	
	sprintf(name,"%s_div.dat", txt_prefix);
	sprintf(nameeps,"%s_div.eps", txt_prefix);
	find_boundaries_heated_map_file2(name, xrange, yrange, &min, &max);
	easyplot3dmap(name, "y [m]", "x [m]", "div B [T/m]", xrange, yrange, nameeps, "size ratio -1", min, max);
	sprintf(name,"%s_curlx.dat", txt_prefix);
	sprintf(nameeps,"%s_curlx.eps", txt_prefix);
	find_boundaries_heated_map_file2(name, xrange, yrange, &min, &max);
	easyplot3dmap(name, "y [m]", "x [m]", "curl Bx [T/m]", xrange, yrange, nameeps, "size ratio -1", min, max);
	sprintf(name,"%s_curly.dat", txt_prefix);
	sprintf(nameeps,"%s_curly.eps", txt_prefix);
	find_boundaries_heated_map_file2(name, xrange, yrange, &min, &max);
	easyplot3dmap(name, "y [m]", "x [m]", "curl By [T/m]", xrange, yrange, nameeps, "size ratio -1", min, max);
	sprintf(name,"%s_curlz.dat", txt_prefix);
	sprintf(nameeps,"%s_curlz.eps", txt_prefix);
	find_boundaries_heated_map_file2(name, xrange, yrange, &min, &max);
	easyplot3dmap(name, "y [m]", "x [m]", "curl Bz [T/m]", xrange, yrange, nameeps, "size ratio -1", min, max);
}

extern void write_in_txtfile(char *textfile, int nbpoints, double r0, double b0, double rmin1, double rmax1, double k1, double rmin2, double rmax2, double k2)
{
	int i;
	double r1,r2,b1,b2, rstep1, rstep2;
	FILE *wfile;
	
	wfile = fopen(textfile,"w");
	
	rstep1 = comp_step(rmin1, rmax1, nbpoints);
	rstep2 = comp_step(rmin2, rmax2, nbpoints);
	for(i=0;i<nbpoints;i++) {
		r1 = rmin1 + i*rstep1;
		r2 = rmin2 + i*rstep2;
		b1 = b0*pow(r1/r0,k1);
		b2 = b0*pow(r2/r0,k2);
		fprintf(wfile, "%le	%le	%le	%le\n", r1,b1,r2,b2);
	}
	fclose(wfile);
}

extern void write_beta_in_txtfile(char *textfile)
{
	int i,nbpoints=1000;
	double beta, beta1, beta2, betastep, beta0, ff0, k, a, b;
	FILE *wfile;
	
	// f/f0=(beta*c)/(2pi*r)/(beta0*c/(2pi*r0))
	beta1 = 0.1;
	beta2 = 0.4;
	beta0 = 0.151794;
	k = 7.5;
	a = 1-1/(k+1.);
	b = 1/(2*(k+1));
	wfile = fopen(textfile,"w");
	
	betastep = comp_step(beta1, beta2, nbpoints);
	for(i=0;i<nbpoints;i++) {
		beta = beta1 + i*betastep;
		ff0 = pow(beta,a)*pow(1-beta, b)/(exp(log(2)+a*log(beta0)+b*log(1-beta0)));
		fprintf(wfile, "%le	%le\n", beta, ff0);
	}
	fclose(wfile);
}

extern void check_tune_mag_pos(char *textfile, struct Particle *part, struct Lattice *latt, int nbpoints)
{
	int i,j;
	double qx, qz, beta_x, alpha_x, beta_z, alpha_z, step_rad, eps_clo = 1.e-12;
	FILE *wfile;
	
	step_rad = 0.002;
	
	wfile = fopen(textfile, "w");
	for(i=0;i<nbpoints;i++) {
		for(j=0;j<latt->cell[0].nbcomp;j++) latt->cell[0].mpara[j][0] += step_rad;
		if(find_closed_orbite_xxp(part, &(part->x), &(part->ux), &(part->uy), eps_clo, latt, YES) != TRUE) errorstop("closed orbit not found");
		calc_tune_twiss(part, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-5, 1.e-5, 1.e-5, 1.e-5, latt, part_cross_latt, YES, NULL,0);
		fprintf(wfile, "%lf	%le	%le\n", latt->cell[0].mpara[0][0], qx, qz);
	}
	fclose(wfile);
}

extern void compute_tm_collimator_for_emi(char *textfile, struct Lattice *latt, struct Particle *part, double angle_collimator_deg)
{
	char message[300], message2[300];
	int i,j;
	double qx, qz, beta_x, alpha_x, beta_z, alpha_z;
	double m[4][4], m1[4][4], m2[4][4], band_bias;
	double mh[2][2], mv[2][2], m1h[2][2], m1v[2][2];//, m2h[2][2],m2v[2][2];
	double m1_inv[4][4], verif[4][4];
	//struct Framework back_fwk;
	struct Particle test_part;
	
	//back_fwk.xc = 0;
	//back_fwk.yc = 0;
	//back_fwk.ae = -angle_collimator_deg*PI/180.;
	test_part = *part;
	emptyfile(textfile);
	if(find_closed_orbite_xxp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), 1.e-10, latt, YES) != TRUE) errorstop("closed orbit not found");
	
	calc_tune_twiss(&test_part, &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-7, 1.e-7, 1.e-7, 1.e-7, latt, part_cross_latt, YES, NULL,0);
	
	for(i=0;i<4;i++) {
		for(j=0;j<4;j++) {
			m[i][j] = 0;
			m1[i][j] = 0;
			m2[i][j] = 0;
		}
	}
	
	//for(i=0;i<2;i++) {
	//	for(j=0;j<2;j++) {
	//		mh[i][j] = 0;
	//		mv[i][j] = 0;
	//		m1h[i][j] = 0;
	//		m1v[i][j] = 0;
	//		m2h[i][j] = 0;
	//		m2v[i][j] = 0;
	//	}
	//}
	
	
	get_trmatrix_h(mh, &band_bias, &test_part, 1.e-7, 1.e-7, latt, part_cross_latt, NO);
	get_trmatrix_v(mv, &test_part, 1.e-7, 1.e-7, latt, part_cross_latt, NO);
	for(i=0;i<2;i++) {
		for(j=0;j<2;j++) {
			m[i][j] = mh[i][j];
			m[i+2][j+2] = mv[i][j];
		}
	}
	//get_matrix_firstorder_4d(m, &band_bias, &test_part, 1.e-7, 1.e-7, 1.e-7, 1.e-7, latt, part_cross_latt, YES);
	sprintf(message, "full cell FETS-FFA tunes (%lf,%lf)", qx,qz);
	printf_mat_4d_real(message, m, "M",textfile);
	
	latt->cell[0].instrutype = CUP;
	latt->cell[0].instru.ymax = 0;
	latt->cell[0].instru.thmax = angle_collimator_deg*PI/180.;
	
	get_trmatrix_h(m1h, &band_bias, &test_part, 1.e-7, 1.e-7, latt, part_cross_latt, NO);
	get_trmatrix_v(m1v, &test_part, 1.e-7, 1.e-7, latt, part_cross_latt, NO);
	for(i=0;i<2;i++) {
		for(j=0;j<2;j++) {
			m1[i][j] = m1h[i][j];
			m1[i+2][j+2] = m1v[i][j];
		}
	}
	//get_matrix_firstorder_4d(m1, &band_bias, &test_part, 1.e-7, 1.e-7, 1.e-7, 1.e-7, latt, part_cross_latt, YES);
	sprintf(message2, "up to collimator (%lf deg after cell start)", angle_collimator_deg);
	printf_mat_4d_real(message2, m1, "M1",textfile);
	//emptyfile("data/trackout_test.dat");
	//part_cross_latt(&(test_part), latt, "data/trackout_test.dat");
	//test_part.s = 0;
	//changefwk_part(&test_part, &back_fwk);
	//printf("part: (%le, %le, %le)\n", test_part.x,test_part.y, test_part.z);
	latt->cell[0].instrutype = NO;
	latt->cell[0].instru.ymax = 0;
	latt->cell[0].instru.thmax = 0;
	//emptyfile("data/trackout_test.dat");
	//part_cross_latt(&(test_part), latt, "data/trackout_test.dat");
	//plot_traj(latt, "data/trackout_test.dat", "data/aspect.dat", "($3)", "($2)", "($2)", "($1)","lines lc 7 lw 2", "lines lt 1 lc 0 lw 2", "y [m]", "x [m]", NULL, NULL, "output/traj_hor_long_test.eps", "size ratio -1\nset mxtics 5\nset mytics 5");
	
	matrix_inverse_4d(m1, m1_inv);
	printf_mat_4d_real("inverse m1", m1_inv, "M1inv",NULL);
	
	mmprod4(m, m1_inv, m2);
	mmprod4(m2, m1, verif);
	printf_mat_4d_real("verif m2*m1", verif, "V",NULL);
	//get_matrix_firstorder_4d(m2, &band_bias, &test_part, 1.e-7, 1.e-7, 1.e-7, 1.e-7, latt, part_cross_latt, YES);
	printf_mat_4d_real("from the collimator (F magnet entrance)", m2, "M2",textfile);
}

extern void gene_matched_beam_emi_4d_hffa(char *emplacement, double x_co, double xprime_co, double emit_x, double twiss_betax, double twiss_alphax, double emit_z, double twiss_betaz, double twiss_alphaz, int nparts, long seed, int doyouhalo)
{
	char filename[500];//, epsx[500], epsz[500];
	long idum;
	int i;
	double x_ctr, xprime_ctr, z_ctr, zprime_ctr;
	double factor_rms_tot, factor_rms_halo;
	double xn, xpn, zn, zpn,x,xp,z,zp;
	FILE *wfile;
	
	factor_rms_tot = 6.; //maximum sigma of the beam
	factor_rms_halo = 4.; //boundary sigma between core (<factor_rms_halo) and halo (>factor_rms_halo)
	x_ctr = x_co;
	xprime_ctr = xprime_co;
	z_ctr = 0;
	zprime_ctr = 0;
	
	printf("gene_matched_beam_emi_4d for hFFA: nparts = %i\n", nparts);
	//test errors and warnings
	if(nparts <= 0) errorstop("!!!ERROR in gene_matched_beam_emi_4d_hffa:\n nparts <= 0");
	if(twiss_betax <= 0) errorstop("!!!ERROR in gene_matched_beam_emi_4d_hffa: twiss_betax <= 0");
	if(twiss_betaz <= 0) errorstop("!!!ERROR in gene_matched_beam_emi_4d_hffa: twiss_betaz <= 0");
	
	//initialise idum (for the random generator)
	if(seed >=0) idum = -seed;
	else idum = seed;
	
	if(doyouhalo==YES) {
		sprintf(filename,"%s/matched_beam_halo.dat", emplacement);
		//sprintf(epsx,"%s/matched_beam_halo_xxprime.eps", txtfile);
		//sprintf(epsz,"%s/matched_beam_halo_zzprime.eps", txtfile);
	}
	else {
		sprintf(filename,"%s/matched_beam_core.dat", emplacement);
		//sprintf(epsx,"%s/matched_beam_core_xxprime.eps", txtfile);
		//sprintf(epsz,"%s/matched_beam_core_zzprime.eps", txtfile);
	}
	
	wfile = fopen(filename,"w");
	
	//generate a 4D beam, inside ellipses centered on central particle
	for(i = 1; i < nparts+1 ; i++) {
		CLRSCR();
		printf("part number %i", i);
		fflush(stdout);
		if(doyouhalo==YES) { //halo particles between factor_rms_halo sigma and factor_rms_tot sigma
			for(;;) {
				xn = gasdev(&idum); //normalised coordinates
				xpn = gasdev(&idum); //normalised coordinates
				if(xn*xn + xpn*xpn > factor_rms_halo*factor_rms_halo && xn*xn + xpn*xpn < factor_rms_tot*factor_rms_tot) break;
			}
			for(;;) {
				zn = gasdev(&idum); //normalised coordinates
				zpn = gasdev(&idum); //normalised coordinates
				if(zn*zn + zpn*zpn > factor_rms_halo*factor_rms_halo && zn*zn + zpn*zpn < factor_rms_tot*factor_rms_tot) break;
			}
		}
		else { // core particles less than factor_rms_halo sigma
			for(;;) {
				xn = gasdev(&idum); //normalised coordinates
				xpn = gasdev(&idum); //normalised coordinates
				if(xn*xn + xpn*xpn < factor_rms_halo*factor_rms_halo) break;
			}
			for(;;) {
				zn = gasdev(&idum); //normalised coordinates
				zpn = gasdev(&idum); //normalised coordinates
				if(zn*zn + zpn*zpn < factor_rms_halo*factor_rms_halo) break;
			}
		}
		//generation of x, xprime, z, zprime from normalised coordinates
		x = x_ctr + sqrt(emit_x*twiss_betax)*xn;
		xp = xprime_ctr + sqrt(emit_x/twiss_betax)*(xpn - twiss_alphax*xn);
		z = z_ctr + sqrt(emit_z*twiss_betaz)*zn;
		zp = zprime_ctr +sqrt(emit_z/twiss_betaz)*(zpn - twiss_alphaz*zn);
		fprintf(wfile, "%le	%le	%le	%le\n", x,xp,z,zp);
	}// and it is a loop, it does it again nparts times...
	fclose(wfile);
	printf("\n");
	//easyplot(filename, "($1*1000)", "($2*1000)", "points pt 7 ps 0.8","x [mm]", "x\251 [mrad]", NULL, NULL, epsx,"mxtics 5\nset mytics 5");
	//easyplot(filename, "($3*1000)", "($4*1000)", "points pt 7 ps 0.8","z [mm]", "z\251 [mrad]", NULL, NULL, epsz,"mxtics 5\nset mytics 5");
}

extern void fets_ffa_corners_tambouille(char *lattfile, char *beamfile, char *trackout, char *field_file)
{
	double qx, qz, beta_x, alpha_x, beta_z, alpha_z;
	double bfmin, bdmin, bfmax, bdmax, k, r0, b0f, b0d, rmin_f, rmax_f_3mev, rmin_d, rmax_d_3mev, rmax_f_12mev, rmax_d_12mev;
	double th_rmin_f, th_rmax_f, th_rmin_d, th_rmax_d;
	double rav, theta_spi_move;
	double size_beam_3mev = 0.04, size_beam_12mev = 0.03;
	struct Lattice latt;
	struct Beam beam;
	FILE *afile_field;
	
	load_lattice(&latt, lattfile);
	load_beam(&beam, beamfile, &latt, YES);
	k = latt.cell[0].mpara[0][3];
	r0 = latt.cell[0].mpara[0][1];
	b0f = latt.cell[0].mpara[0][2];
	
	calc_tune_twiss(&(beam.part[0]), &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-5, 1.e-5, 1.e-5, 1.e-5, &latt, part_cross_latt, YES, NULL,0);
	emptyfile(trackout);
	part_cross_latt(&(beam.part[0]), &latt, trackout);
	rav = compute_av_radius_from_trackout_1cell(trackout, YES);
	theta_spi_move = tan(latt.cell[0].mpara[0][4])*log(rav/r0)*180./PI;
	printf("move theta_spi: %le\n", theta_spi_move);
	compute_min_max_radius_from_trackout_1cell(trackout, &rmin_f, &th_rmin_f, &rmax_f_3mev, &th_rmax_f, 7.+theta_spi_move, 12.+theta_spi_move, YES);
	compute_min_max_radius_from_trackout_1cell(trackout, &rmin_d, &th_rmin_d, &rmax_d_3mev, &th_rmax_d, 13.75+theta_spi_move, 16.5+theta_spi_move, YES);
	//printf("rmin_f = %le\n", rmin_f);
	//printf("rmax_f = %le\n", rmax_f_3mev);
	rmax_f_12mev = rmax_f_3mev*pow(2.,(1/(k+1.)));
	rmax_d_12mev = rmax_d_3mev*pow(2.,(1/(k+1.)));
	
	rmin_f -= size_beam_3mev;
	rmin_d -= size_beam_3mev;
	rmax_f_12mev += size_beam_12mev;
	rmax_d_12mev += size_beam_12mev;
	if(b0f > 0) {
		b0d = b0f;
		b0f = latt.cell[0].mpara[1][2];
	}
	else b0d = latt.cell[0].mpara[1][2];
	
	bfmax = b0f*pow(rmax_f_12mev/r0, k);
	bdmax = b0d*pow(rmax_d_12mev/r0, k);
	bfmin = b0f*pow(rmin_f/r0, k);
	bdmin = b0d*pow(rmin_d/r0, k);
		
	
	afile_field = fopen(field_file, "a");
	//fprintf(afile_field, "%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le\n", qx, qz, k, rmin_f, bfmin, rmax_f_12mev, bfmax, rmin_d, bdmin, rmax_d_12mev, bdmax);
	//fprintf(afile_field, "%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le\n", qx, qz, k, rmin_f, th_rmin_f, bfmin, rmax_f_12mev, th_rmax_f, bfmax, rmin_d, th_rmin_d, bdmin, rmax_d_12mev, th_rmax_d, bdmax);
	fprintf(afile_field, "%.1f	%.1f	%.1f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f\n", 
	round_nb(qx,1), round_nb(qz,1), round_nb(k,1), round_nb(rmin_f,2), round_nb(bfmin,2), round_nb(rmax_f_12mev,2), round_nb(bfmax,2), round_nb(rmin_d,2), round_nb(bdmin,2), round_nb(rmax_d_12mev,2), round_nb(bdmax,2));
	
	fclose(afile_field);
	free_latt(&latt);
	free_beam(&beam);
}

extern void compute_min_max_radius_from_trackout_1cell(char *trackout, double *rmin, double *th_rmin, double *rmax, double *th_rmax, double thmin_deg, double thmax_deg, int doyouprintf)
{
	int i, nblines;
	double s, x, y, z, bx, by, bz, ux, uy, uz, brho, t, r, th_deg;
	FILE *rfile = NULL;
	*rmin = 1.e9;
	*rmax = 0;
	nblines = get_nb_lines_file(trackout);
	rfile = fopen(trackout, "r");
	if(rfile == NULL) errorstop("!!!ERROR in compute_av_radius_from_trackout_1cell: cannot open inputfile!!!");
	
	for(i=0;i<nblines-1;i++) {
		fscanf(rfile, "%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le", &s, &x, &y, &z, &bx, &by, &bz, &ux, &uy, &uz, &brho, &t, &r, &th_deg);
		//printf("th_deg = %lf, r = %lf\n", th_deg, r);
		if(th_deg>thmin_deg && th_deg<thmax_deg) {
			//printf("YES\n");
			//printf("rmin = %lf, r = %lf\n", *rmin, r);
			if(*rmin>r) {
				*rmin = r;
				*th_rmin = th_deg;
			}
			if(*rmax<r) {
				*rmax = r;
				*th_rmax = th_deg;
			}
		}
		newline(rfile);
		
	}
	fclose(rfile);
	if(doyouprintf==YES) printf("between %lf and %lf, rmin=%lf at th=%lf deg, rmax=%lf at th=%lf deg\n", thmin_deg, thmax_deg, *rmin, *th_rmin, *rmax, *th_rmax);
}

extern void check_rav_fets_ffa(char *lattfile, char *injfile, char *extfile)
{
	double rav_3mev, rav_12mev;
	struct Lattice latt;
	struct Beam beam_3mev, beam_12mev;
	
	load_lattice(&latt, lattfile);
	load_beam(&beam_3mev, injfile, &latt, NO);
	load_beam(&beam_12mev, extfile, &latt, NO);
	
	emptyfile("data/trackout3mev.dat");
	part_cross_latt(&(beam_3mev.part[0]), &latt, "data/trackout3mev.dat");
	compute_av_radius_from_trackout_1cell("data/trackout3mev.dat", YES);
	
	emptyfile("data/trackout12mev.dat");
	part_cross_latt(&(beam_12mev.part[0]), &latt, "data/trackout12mev.dat");
	compute_av_radius_from_trackout_1cell("data/trackout12mev.dat", YES);
	
	
	free_latt(&latt);
	free_beam(&beam_3mev);
	free_beam(&beam_12mev);
}

extern void round_number_fets_ffa_field_file(char *textout, char * textin)
{
	int i, nblines;
	double qx, qz, k, rmin_f, bfmin, rmax_f_12mev, bfmax, rmin_d, bdmin, rmax_d_12mev, bdmax;
	FILE *wfile;
	FILE *rfile = NULL;
	
	errorstop("not working anymore!! rfile changed!!");
	nblines = get_nb_lines_file(textin);
	rfile = fopen(textin, "r");
	if(rfile==NULL) errorstop("cannot open file");
	wfile = fopen(textout, "w");
	for(i=0;i<nblines;i++) {
		fscanf(rfile, "%le	%le	%le	%le	%le	%le	%le	%le	%le	%le	%le", &qx, &qz, &k, &rmin_f, &bfmin, &rmax_f_12mev, &bfmax, &rmin_d, &bdmin, &rmax_d_12mev, &bdmax);
		fprintf(wfile, "%.1f	%.1f	%.1f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f\n", 
		round_nb(qx,1), round_nb(qz,1), round_nb(k,1), round_nb(rmin_f,2), round_nb(bfmin,2), round_nb(rmax_f_12mev,2), round_nb(bfmax,2), round_nb(rmin_d,2), round_nb(bdmin,2), round_nb(rmax_d_12mev,2), round_nb(bdmax,2));
	}
	
	fclose(rfile);
	fclose(wfile);
	
}

extern void plot_latt_fets_ffa(char *latt_name, char *beamname, char *trackout)
{
	struct Lattice latt;
	struct Beam beam;
	
	load_lattice(&latt, latt_name);
	load_beam(&beam, beamname, &latt, NO);
	
	//part_cross_latt(&(beam.part[0]), &latt, trackout);
	part_oneturn(&(beam.part[0]), &latt, trackout);
	
	free_latt(&latt);
	free_beam(&beam);
}

/*extern void write_latex_table_fets_ffa()
{
	
	fprintf(wfile, "\\\begin{table}[!h]\n");
	fprintf(wfile, "\\centering\n");
	fprintf(wfile, " \\caption{Parameters of the magnet with fixed average injection radius.}\n");
	fprintf(wfile, "\\begin{tabular}{|c|c|c|c|c|}\n");
	fprintf(wfile, "\\hline\n");
	fprintf(wfile, "		&	& $Q_x=3.0$ & $Q_x=3.5$ & $Q_x=4.0$ \\\\ \n");
	fprintf(wfile, "\\hline\n");
	fprintf(wfile, " & $k-$value     & %.1f       & %.1f      & %.1f\\\\\n", k[0][0], k[0][1], k[0][2]);
	fprintf(wfile, "  \\cline{2-5}\n");
	fprintf(wfile, "		& F magnet:	&     &  & \\\\\n");
	fprintf(wfile, "		& r$_{\\text{min,max}}$ [m]	 & (4.01, 4.51)  &  (4.01, 4.42) & (4.01, 4.36)\\\\\n");
	fprintf(wfile, "		& $B(\\text{r}_{\\text{min}})$ [T] & -0.46 & -0.47 & -0.47 \\\\\n");
	fprintf(wfile, "$Q_y=4.0$ & $B(\\text{r}_{\\text{max}})$ [T] & -0.94 & -1.02 & -1.09\\\\\n");
	fprintf(wfile, " \\cline{2-5}\n");
	fprintf(wfile, "		& D magnet:	&     &  & \\\\\n");
	fprintf(wfile, "		& r$_{\\text{min,max}}$ [m]	 & (3.95, 4.40)  & (3.94, 4.31)  & (3.94, 4.25) \\\\\n");
	fprintf(wfile, "		& $B(\\text{r}_{\\text{min}})$ [T] & 0.34 & 0.37 & 0.40 \\\\\n");
	fprintf(wfile, "		& $B(\\text{r}_{\\text{max}})$ [T] &0.66  & 0.76 & 0.84 \\\\\n");
	fprintf(wfile, "\\hline\n");
	fprintf(wfile, " & $k-$value     & %.1f       & %.1f      & %.1f\\\\\n", k[1][0], k[1][1], k[1][2]);
	fprintf(wfile, " \\cline{2-5}\n");
	fprintf(wfile, "		& F magnet:	&     &  & \\\\\n");
	fprintf(wfile, "		& r$_{\\text{min,max}}$ [m]	 & (4.01, 4.49)  & (4.01, 4.41) & (4.01, 4.35)\\\\\n");
	fprintf(wfile, "		& $B_0(\\text{r}_{\\text{min}})$ [T] & -0.42 & -0.44 & -0.45 \\\\\n");
	fprintf(wfile, "$Q_y=3.5$ & $B_0(\\text{r}_{\\text{max}})$ [T] & -0.87 & -0.96 & -1.03 \\\\\n");
	fprintf(wfile, " \\cline{2-5}\n");
	fprintf(wfile, "		& D magnet:	&     &  & \\\\\n");
	fprintf(wfile, "		& r$_{\\text{min,max}}$ [m]	 & (3.95, 4.39)  & (3.95, 4.30)  &  (3.95, 4.25) \\\\\n");
	fprintf(wfile, "		& $B_0(\\text{r}_{\\text{min}})$ [T] & 0.27 & 0.31 & 0.34 \\\\\n");
	fprintf(wfile, "		& $B_0(\\text{r}_{\\text{max}})$ [T] & 0.53 & 0.63 & 0.72 \\\\\n");
	fprintf(wfile, "\\hline\n");
	fprintf(wfile, " & $k-$value     & %.1f       & %.1f      & %.1f\\\\\n", k[2][0], k[2][1], k[2][2]);
	fprintf(wfile, "  \\cline{2-5}\n");
	fprintf(wfile, "		& F magnet:	&     &  & \\\\\n");
	fprintf(wfile, "		& r$_{\\text{min,max}}$ [m]	 & (4.01, 4.48)  & (4.01, 4.39) & (4.01, 4.34)\\\\\n");
	fprintf(wfile, "		& $B_0(\\text{r}_{\\text{min}})$ [T] & -0.39 & -0.41 & -0.42 \\\\\n");
	fprintf(wfile, "$Q_y=3.0$ & $B_0(\\text{r}_{\\text{max}})$ [T] & -0.81 & -0.90 & -0.98 \\\\\n");
	fprintf(wfile, " \\cline{2-5}\n");
	fprintf(wfile, "		& D magnet:	&     &  & \\\\\n");
	fprintf(wfile, "		& r$_{\\text{min,max}}$ [m]	 & (3.95, 4.39)  & (3.95, 4.30) & (3.95, 4.24)\\\\\n");
	fprintf(wfile, "		& $B_0(\\text{r}_{\\text{min}})$ [T] & 0.20 & 0.25 & 0.28 \\\\\n");
	fprintf(wfile, "		& $B_0(\\text{r}_{\\text{max}})$ [T] & 0.40 & 0.51 & 0.60 \\\\\n");
	fprintf(wfile, "\\hline\n");
	fprintf(wfile, "\\end{tabular}\n");
	fprintf(wfile, " \\label{general}\n");
}//*/

extern void print_get_bfield_across_cell_radius(char *fileout, double r, int nbsteps, double z, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double x, y, bx, by, bz;
	double br,bth, th, thstep;
	FILE *wfile;
	wfile = fopen(fileout,"a");
	thstep = comp_step(0, cell->boun.thmax, nbsteps);
	
	for(i = 0; i < nbsteps + 1; i++) {
		th = i*thstep;
		x = r*cos(th);
		y = r*sin(th);
		if(get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp)!=ALIVE) errorstop("get_bfield not ALIVE");
		br = bx*cos(th) + by*sin(th);
		bth = by*cos(th) - bx*sin(th);
		fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r, th*180./PI, z, br, bth, bz);
		//fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r, th*180./PI, z, bx, by, bz);
	}
	
}

extern void comp_betafunc_lattices_fets_ffa(char *latt_name, char *beamname, char *betafile)
{
	
	struct Lattice latt;
	struct Beam beam;
	
	load_lattice(&latt, latt_name);
	load_beam(&beam, beamname, &latt, YES);
	compute_periodic_betafunc(betafile, 100, 1.e-4, 1.e-4, 1.e-4, 1.e-4, &(beam.part[0]), &latt, 0);
	
	free_latt(&latt);
	free_beam(&beam);
}


// compute acceptance in (u,u') and (v,v') phase spaces with an initial amplitude in u and v
//return ALIVE if successful
extern int acceptance_uv(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double au, double av, double decoup_to_coup[4][4])
{
	int i;
	double x_vect[4], u_vect[4], decoup_to_coup_inv[4][4];
	double php_ref, xprime_ref, zprime_ref;//, band_bias, amp_v;
	struct Particle test_part;
	FILE *rfile;
	if(outfilename != NULL) rfile = fopen(outfilename,"a");

	printf("acceptance in v:");
	php_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	zprime_ref  = atan_ratio(reference->uz, php_ref);
	
	u_vect[0] = au; // u
	u_vect[1] = 0; //uprime
	u_vect[2] = av; //v
	u_vect[3] = 0; //vprime
	mvprod4(x_vect, decoup_to_coup, u_vect); //xvect physical space vector x,x',z,z'
	test_part = *reference;
	test_part.s = 0;
	test_part.x += x_vect[0]/cos(xprime_ref);
	test_part.z += x_vect[2]/cos(zprime_ref);
	test_part.ux = sin(xprime_ref + x_vect[1])*cos(zprime_ref + x_vect[3]);
	test_part.uy = cos(xprime_ref + x_vect[1])*cos(zprime_ref + x_vect[3]);
	test_part.uz = sin(zprime_ref + x_vect[3]);
	test_part.hat = -2;
	printf("au = %lf, av = %lf\n", au, av);
	
	for(i = 0; i < nbpass; i++) {
		//CLRSCR();
		//printf("pass number %i\t", i);
		//fflush(stdout);
		if (test_part.uy > 0) {
			comp_phase_space_coord(&(x_vect[0]), &(x_vect[1]), &(x_vect[2]), &(x_vect[3]), test_part.x, test_part.z, test_part.ux, test_part.uy, test_part.uz, reference->x, xprime_ref, reference->z, zprime_ref);
			matrix_inverse_4d(decoup_to_coup, decoup_to_coup_inv);
			mvprod4(u_vect, decoup_to_coup_inv, x_vect);
		}
		else {
			printf("\tin acceptancev, test_part.uy <= 0!, particle is going backwards!\n");
			if(outfilename != NULL) fclose(rfile);
			return test_part.status;
		}
		if(outfilename != NULL) fprintf(rfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n", x_vect[0], x_vect[1], x_vect[2], x_vect[3], u_vect[0], u_vect[1], u_vect[2], u_vect[3], au, av);
		part_cross_latt(&(test_part), latt, NULL);
		if(test_part.status != ALIVE) {
			if(outfilename != NULL) fclose(rfile);
			printf("\n");
			return test_part.status;
		}
	}
	
	
	if(outfilename != NULL) fclose(rfile);
	printf("\n");
	return test_part.status;
}

// geometry parameters and name with note on fixfield 03/02/2020, figure 10 (triplet)
extern void change_para_triplet_cell_rect_vffa(struct Cell *cell, double cell_length, double theta_cell_deg, double magnetf_length, double small_drift, double magnetd_length, double shiftf, double shiftd, double tilt_magf_deg, double b0f, double b0d, double m)
{
	int i;
	double r0d, thetad;
	double r0f, thetaf1, thetaf2, phif1, phif2;
	double ac, theta_cell, tilt_f, eb, oa, ob;
	
	ac = cell_length;
	theta_cell = theta_cell_deg*PI/180.;
	tilt_f = tilt_magf_deg*PI/180.;
	eb = magnetf_length/2. + magnetd_length/2. + small_drift;
	oa = ac/(2*sin(theta_cell/2.));
	printf("r0=%lf (=OA in report)\n", oa);
	ob = ac/(2*tan(theta_cell/2.));
	
	thetad = theta_cell/2.;
	r0d = ob+shiftd;
	
	thetaf1 = thetad - atan(eb/(ob+shiftf));
	r0f = (ob+shiftf)/(cos(thetad-thetaf1));
	phif1 = tilt_f + thetad - thetaf1;
	
	thetaf2 = theta_cell - thetaf1;
	phif2 = -phif1;
	printf("phif1 = %lf\n",phif1);
	
	//B0:
	cell->mpara[0][2] = b0f;
	cell->mpara[1][2] = b0d;
	cell->mpara[2][2] = b0f;
	cell->mpara[3][2] = b0f;
	cell->mpara[4][2] = b0d;
	cell->mpara[5][2] = b0f;
	cell->mpara[6][2] = b0f;
	cell->mpara[7][2] = b0d;
	cell->mpara[8][2] = b0f;
	//m:
	for(i=0;i<cell->nbcomp;i++) cell->mpara[i][4] = m;
	//r0:
	cell->mpara[0][1] = r0f;
	cell->mpara[1][1] = r0d;
	cell->mpara[2][1] = r0f;
	cell->mpara[3][1] = r0f;
	cell->mpara[4][1] = r0d;
	cell->mpara[5][1] = r0f;
	cell->mpara[6][1] = r0f;
	cell->mpara[7][1] = r0d;
	cell->mpara[8][1] = r0f;
	//theta0
	cell->mpara[0][0] = thetaf1-theta_cell;
	cell->mpara[1][0] = thetad-theta_cell;
	cell->mpara[2][0] = thetaf2-theta_cell;
	cell->mpara[3][0] = thetaf1;
	cell->mpara[4][0] = thetad;
	cell->mpara[5][0] = thetaf2;
	cell->mpara[6][0] = thetaf1+theta_cell;
	cell->mpara[7][0] = thetad+theta_cell;
	cell->mpara[8][0] = thetaf2+theta_cell;
	//phif:
	cell->mpara[0][5] = phif1;
	cell->mpara[2][5] = phif2;
	cell->mpara[3][5] = phif1;
	cell->mpara[5][5] = phif2;
	cell->mpara[6][5] = phif1;
	cell->mpara[8][5] = phif2;
}


extern void check_fit_iker(char *textfile)
{
	char track[300], fiteps[300], start[100], xrange[100];
	int i,nblines, skiplines=8;
	double s, b, smax=0, bmax = -1.e9, half_length=0.15;
	FILE *rfile=NULL;
	FILE *wfile;
	
	nblines = get_nb_lines_file(textfile);
	printf("opening %s...\n", textfile);
	
	rfile = fopen(textfile, "r");
	for(i=0;i<skiplines;i++) newline(rfile);
	sprintf(track, "%s_clean.dat", textfile);
	wfile = fopen(track, "w");
	
	for(i=0;i<nblines-skiplines-1;i++) {
		fscanf(rfile, "%lf      %le", &s, &b);
		fprintf(wfile, "%lf	%le\n", s, b);
		//printf("%lf	%le\n", s,b);
		if(bmax<b) {
			bmax = b;
			smax = s;
		}
	}
	fclose(rfile);
	fclose(wfile);
	printf("bmax = b(%lf) = %lf\n", smax, bmax);
	
	sprintf(start, "a = %lf\n ben=%lf\nbex=%lf\n c=50", bmax, smax-half_length, smax+half_length);
	sprintf(fiteps, "%s_fit", textfile);
	sprintf(xrange, "[%lf:%lf]", smax-1.,smax+1.);
	gnuplot_fit_general(track, "1", "2", "a/((1+exp(c*(x-bex)))*(1+exp(c*(ben-x))))", "a,ben,bex,c", start, "lines lc 7 lw 2", "lines dt 1 lc 0 lw 2", "s [m]", "B [T]", "Comsol", "fit", xrange, NULL, fiteps, NULL);
	
}


extern void check_fit_iker2(char *textfile)
{
	char track[300], fiteps[300], start[100], equation[300],xrange[100];
	int i,nblines, imax, skiplines=8, ien, iex;
	nblines = get_nb_lines_file(textfile);
	double s[nblines], b[nblines], smax=0, bmax = -1.e9, sen, sex;
	FILE *rfile=NULL;
	FILE *wfile;
	
	printf("opening %s...\n", textfile);
	
	rfile = fopen(textfile, "r");
	for(i=0;i<skiplines;i++) newline(rfile);
	sprintf(track, "%s_clean.dat", textfile);
	wfile = fopen(track, "w");
	
	for(i=0;i<nblines-skiplines-1;i++) {
		fscanf(rfile, "%lf      %le", &(s[i]), &(b[i]));
		fprintf(wfile, "%lf	%le\n", s[i], b[i]);
		//printf("%lf	%le\n", s,b);
		if(bmax<b[i]) {
			bmax = b[i];
			smax = s[i];
			imax = i;
		}
	}
	for(i=imax;i<nblines-skiplines-1;i++) {
		if(b[i]<bmax/2) {
			sex = s[i];
			ien = i;
			break;
		}
	}
	for(i=0;i<imax;i++) {
		if(b[i]>bmax/2) {
			sen = s[i];
			iex = i;
			break;
		}
	}
	fclose(rfile);
	fclose(wfile);
	printf("bmax = b(%lf) = %lf\n", smax, bmax);
	printf("sen = %lf, sex = %lf\n", sen, sex);
	printf("b(sen) = %lf, b(sex) = %lf\n", b[ien], b[iex]);
	
	sprintf(start, "a = %lf\n c=50", bmax);
	sprintf(equation, "a/((1+exp(c*(x-%lf)))*(1+exp(c*(%lf-x))))", sex, sen);
	sprintf(fiteps, "%s_fit2", textfile);
	sprintf(xrange, "[%lf:%lf]", smax-1.,smax+1.);
	gnuplot_fit_general(track, "1", "2", equation, "a,c", start, "lines dt 1 lc 7 lw 4", "lines dt 1 lc 0 lw 2", "s [m]", "B [T]", "Comsol", "fit", xrange, "[0:0.6]", fiteps, NULL);
}

extern void check_fit_iker3(char *textfile)
{
	char track[300], fiteps[300], start[100], equation[300], xrange[100];
	int i,nblines, skiplines=8;
	int ien, iex, imax;
	nblines = get_nb_lines_file(textfile);
	
	double s[nblines], b[nblines], smax=0, bmax = -1.e9, sen, sex;
	//double half_length=0.15;
	FILE *rfile=NULL;
	FILE *wfile;
	
	printf("opening %s...\n", textfile);
	
	rfile = fopen(textfile, "r");
	for(i=0;i<skiplines;i++) newline(rfile);
	sprintf(track, "%s_clean.dat", textfile);
	wfile = fopen(track, "w");
	
	for(i=0;i<nblines-skiplines-1;i++) {
		fscanf(rfile, "%lf      %le", &(s[i]), &(b[i]));
		fprintf(wfile, "%lf	%le\n", s[i], b[i]);
		//printf("%lf	%le\n", s,b);
		if(bmax<b[i]) {
			bmax = b[i];
			smax = s[i];
			imax = i;
		}
	}
	for(i=imax;i<nblines-skiplines-1;i++) {
		if(b[i]<bmax/2) {
			sex = s[i];
			ien = i;
			break;
		}
	}
	for(i=0;i<imax;i++) {
		if(b[i]>bmax/2) {
			sen = s[i];
			iex = i;
			break;
		}
	}
	fclose(rfile);
	fclose(wfile);
	printf("bmax = b(%lf) = %lf\n", smax, bmax);
	// p = c0 + c1*(x/l) + c2*pow(x/l,2) + c3*pow(x/l,3), C0 = 0.1455, C1 = 2.267, C2 = -0.6395, C3 = 1.1558
	//f = B0/((1+exp(P(x-xex)))*(1+exp(P(xen-x)))

	
	//sprintf(start, "a = %lf\n ben=%lf\nbex=%lf\n l=0.12", bmax, smax-half_length, smax+half_length);
	sprintf(start, "a = %lf\n l=0.12", bmax);
	printf("%s\n",start);
	sprintf(equation, "a/((1+exp(0.1455+2.267*((x-(%lf))/l)-0.6395*((x-(%lf))/l)**2+1.1558*((x-(%lf))/l)**3))*(1+exp(0.1455+2.267*(((%lf)-x)/l)-0.6395*(((%lf)-x)/l)**2+1.1558*(((%lf)-x)/l)**3)))", sex, sex, sex, sen, sen, sen);
	sprintf(fiteps, "%s_fit3", textfile);
	sprintf(xrange, "[%lf:%lf]", smax-1.,smax+1.);
	//gnuplot_fit_general(track, "1", "2", "a/((1+exp(0.1455+2.267*((x-bex)/l)-0.6395*((x-bex)/l)**2+1.1558*((x-bex)/l)**3))*(1+exp(0.1455+2.267*((ben-x)/l)-0.6395*((ben-x)/l)**2+1.1558*((ben-x)/l)**3)))", "a,ben,bex,l", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", NULL, NULL, fiteps, NULL);
	gnuplot_fit_general(track, "1", "2", equation, "a,l", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", xrange, NULL, fiteps, NULL);
	
}

extern void check_fit_iker4(char *textfile)
{
	char track[300], fiteps[300], start[100], equation[300],xrange[100];
	int i,nblines, skiplines=8;
	int ien, iex, imax;
	nblines = get_nb_lines_file(textfile);
	
	double s[nblines], b[nblines], smax=0, bmax = -1.e9, sen, sex, l=0.12;
	//double half_length=0.15;
	FILE *rfile=NULL;
	FILE *wfile;
	
	printf("opening %s...\n", textfile);
	
	rfile = fopen(textfile, "r");
	for(i=0;i<skiplines;i++) newline(rfile);
	sprintf(track, "%s_clean.dat", textfile);
	wfile = fopen(track, "w");
	
	for(i=0;i<nblines-skiplines-1;i++) {
		fscanf(rfile, "%lf      %le", &(s[i]), &(b[i]));
		fprintf(wfile, "%lf	%le\n", s[i], b[i]);
		//printf("%lf	%le\n", s,b);
		if(bmax<b[i]) {
			bmax = b[i];
			smax = s[i];
			imax = i;
		}
	}
	for(i=imax;i<nblines-skiplines-1;i++) {
		if(b[i]<bmax/2) {
			sex = s[i];
			ien = i;
			break;
		}
	}
	for(i=0;i<imax;i++) {
		if(b[i]>bmax/2) {
			sen = s[i];
			iex = i;
			break;
		}
	}
	fclose(rfile);
	fclose(wfile);
	printf("bmax = b(%lf) = %lf\n", smax, bmax);
	// p = c0 + c1*(x/l) + c2*pow(x/l,2) + c3*pow(x/l,3), c0 = 0.1455\n c1 = 2.267\n c2 = -0.6395\n c3 = 1.1558
sprintf(xrange, "[%lf:%lf]", smax-1.,smax+1.);
	
	//sprintf(start, "a = %lf\n ben=%lf\nbex=%lf\n l=0.12", bmax, smax-half_length, smax+half_length);
	sprintf(start, "a = %lf\n b = 0.1455\n c = 2.267\n d = -0.6395\n e = 1.1558", bmax);
	sprintf(equation, "a/((1+exp(b+c*((x-%lf)/%lf)+d*((x-%lf)/%lf)**2+e*((x-%lf)/%lf)**3))*(1+exp(b+c*((%lf-x)/%lf)+d*((%lf-x)/%lf)**2+e*((%lf-x)/%lf)**3)))", sex, l, sex, l, sex, l, sen, l, sen, l, sen, l);
	sprintf(fiteps, "%s_fit4", textfile);
	sprintf(xrange, "[%lf:%lf]", smax-1,smax+1);
	//gnuplot_fit_general(track, "1", "2", "a/((1+exp(0.1455+2.267*((x-bex)/l)-0.6395*((x-bex)/l)**2+1.1558*((x-bex)/l)**3))*(1+exp(0.1455+2.267*((ben-x)/l)-0.6395*((ben-x)/l)**2+1.1558*((ben-x)/l)**3)))", "a,ben,bex,l", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", NULL, NULL, fiteps, NULL);
	gnuplot_fit_general(track, "1", "2", equation, "a,b,c,d,e", start, "lines lc 7 lw 2", "lines dt 1 lc 0 lw 2", "s [m]", "B [T]", "Comsol", "fit", xrange, NULL, fiteps, NULL);
	
}

extern void check_fit_iker5(char *textfile, double r0, double opening_angle_deg, double lambda)
{
	char track[300], fiteps[300], start[100], equation[300],xrange[100];
	int i,nblines, skiplines=8;
	int ien, iex, imax;
	nblines = get_nb_lines_file(textfile);
	
	double s[nblines], b[nblines], smax=0, bmax = -1.e9, sen, sex, ben, bex, c0en, c0ex;
	//double half_length=0.15;
	FILE *rfile=NULL;
	FILE *wfile;
	
	printf("opening %s...\n", textfile);
	
	rfile = fopen(textfile, "r");
	for(i=0;i<skiplines;i++) newline(rfile);
	sprintf(track, "%s_clean.dat", textfile);
	wfile = fopen(track, "w");
	
	for(i=0;i<nblines-skiplines-1;i++) {
		fscanf(rfile, "%lf      %le", &(s[i]), &(b[i]));
		fprintf(wfile, "%lf	%le\n", s[i], b[i]);
		//printf("%lf	%le\n", s,b);
		if(bmax<b[i]) {
			bmax = b[i];
			smax = s[i];
			imax = i;
		}
	}
	sex = smax+r0*opening_angle_deg/2.*PI/180.;
	sen = smax-r0*opening_angle_deg/2.*PI/180.;
	for(i=imax;i<nblines-skiplines-1;i++) {
		if(s[i]>sex) {
			bex = b[i];
			c0ex = log(bmax/bex-1);
			iex = i;
			break;
		}
	}
	for(i=0;i<imax;i++) {
		if(s[i]>sen) {
			ben = b[i];
			c0en = log(bmax/ben-1);
			ien = i;
			break;
		}
	}
	fclose(rfile);
	fclose(wfile);
	printf("bmax = b(%lf) = %lf\n", smax, bmax);
	printf("ben = b(%lf) = %lf, c0en=%lf\n", sen, ben, c0en);
	printf("bex = b(%lf) = %lf, c0ex=%lf\n", sex, bex, c0ex);
	// p = c0 + c1*(x/l) + c2*pow(x/l,2) + c3*pow(x/l,3), c0 = 0.1455\n c1 = 2.267\n c2 = -0.6395\n c3 = 1.1558
	sprintf(xrange, "[%lf:%lf]", smax-1.,smax+1.);
	
	//sprintf(start, "a = %lf\n ben=%lf\nbex=%lf\n l=0.12", bmax, smax-half_length, smax+half_length);
	//sprintf(start, "a = %lf\n b = %lf\n c = 2.267\n d = -0.6395\n e = 1.1558", bmax, log(2*bmax/(ben+bex)-1));
	sprintf(start, "a = %lf\n b = -0.27\n c = 2.41\n d = -1\n e = 0.23", bmax);
	//sprintf(start, "a = %lf\n b = %lf\n c = 2.267\n d = -0.6395\n e = 1.1558", bmax, (c0en+c0ex)/2.);
	//sprintf(start, "a = %lf\n b = 0.1455\n c = 2.267\n d = -0.6395\n e = 1.1558", bmax);
	//sprintf(equation, "a/((1+exp(b+c*((x-%lf)/%lf)+d*((x-%lf)/%lf)**2+e*((x-%lf)/%lf)**3))*(1+exp(b+c*((%lf-x)/%lf)+d*((%lf-x)/%lf)**2+e*((%lf-x)/%lf)**3)))", sex, lambda, sex, lambda, sex, lambda, sen, lambda, sen, lambda, sen, lambda);
	sprintf(equation, "a/((1+exp(b+c*((x-%lf)/%lf)+d*((x-%lf)/%lf)**2+e*((x-%lf)/%lf)**3))*(1+exp(b+c*((%lf-x)/%lf)+d*((%lf-x)/%lf)**2+e*((%lf-x)/%lf)**3)))", sex, lambda, sex, lambda, sex, lambda, sen, lambda, sen, lambda, sen, lambda);
	sprintf(fiteps, "%s_fit5", textfile);
	sprintf(xrange, "[%lf:%lf]", smax-0.95,smax+0.95);
	//gnuplot_fit_general(track, "1", "2", "a/((1+exp(0.1455+2.267*((x-bex)/l)-0.6395*((x-bex)/l)**2+1.1558*((x-bex)/l)**3))*(1+exp(0.1455+2.267*((ben-x)/l)-0.6395*((ben-x)/l)**2+1.1558*((ben-x)/l)**3)))", "a,ben,bex,l", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", NULL, NULL, fiteps, NULL);
	gnuplot_fit_general(track, "1", "2", equation, "a,b,c,d,e", start, "lines dt 1 lc 7 lw 4", "lines dt 1 lc 0 lw 2", "s [m]", "B [T]", "Comsol", "fit", xrange, "[0:0.6]", fiteps, NULL);
	
}

extern void check_fit_iker_doublet2(char *textfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda)
{
	char track[300], fiteps[300], start[300], equationf[300], equationd[300],xrange[100], equationtot[500];
	int i,nblines, skiplines=8;
	int ien, iex, ifmax, idmax;
	nblines = get_nb_lines_file(textfile);
	
	double s[nblines], b[nblines], sfmax=0, bfmax = 1.e9, sfen, sfex, bfen, bfex, cf0en, cf0ex, sdmax=0, bdmax = -1.e9, sden, sdex, bden, bdex, cd0en, cd0ex;
	//double half_length=0.15;
	FILE *rfile=NULL;
	FILE *wfile;
	
	printf("opening %s...\n", textfile);
	
	rfile = fopen(textfile, "r");
	for(i=0;i<skiplines;i++) newline(rfile);
	sprintf(track, "%s_clean.dat", textfile);
	wfile = fopen(track, "w");
	
	for(i=0;i<nblines-skiplines-1;i++) {
		fscanf(rfile, "%lf      %le", &(s[i]), &(b[i]));
		fprintf(wfile, "%lf	%le\n", s[i], b[i]);
		//printf("%lf	%le\n", s,b);
		if(bfmax>b[i]) {
			bfmax = b[i];
			sfmax = s[i];
			ifmax = i;
		}
		if(bdmax<b[i]) {
			bdmax = b[i];
			sdmax = s[i];
			idmax = i;
		}
	}
	sfex = sfmax+r0*opening_angle_deg_f/2.*PI/180.;
	sfen = sfmax-r0*opening_angle_deg_f/2.*PI/180.;
	sdex = sdmax+r0*opening_angle_deg_d/2.*PI/180.;
	sden = sdmax-r0*opening_angle_deg_d/2.*PI/180.;
	//for(i=imax;i<nblines-skiplines-1;i++) {
	//	if(s[i]>sfex) {
	//		bex = b[i];
	//		c0ex = log(bmax/bex-1);
	//		iex = i;
	//		break;
	//	}
	//}
	//for(i=0;i<imax;i++) {
	//	if(s[i]>sen) {
	//		ben = b[i];
	//		c0en = log(bmax/ben-1);
	//		ien = i;
	//		break;
	//	}
	//}
	fclose(rfile);
	fclose(wfile);
	printf("bfmax = b(%lf) = %lf\n", sfmax, bfmax);
	printf("bdmax = b(%lf) = %lf\n", sdmax, bdmax);
	//printf("ben = b(%lf) = %lf, c0en=%lf\n", sen, ben, c0en);
	//printf("bex = b(%lf) = %lf, c0ex=%lf\n", sex, bex, c0ex);
	// p = c0 + c1*(x/l) + c2*pow(x/l,2) + c3*pow(x/l,3), c0 = 0.1455\n c1 = 2.267\n c2 = -0.6395\n c3 = 1.1558
	sprintf(xrange, "[%lf:%lf]", sfmax-1.,sdmax+1.);
	
	//sprintf(start, "a = %lf\n ben=%lf\nbex=%lf\n l=0.12", bmax, smax-half_length, smax+half_length);
	//sprintf(start, "a = %lf\n b = %lf\n c = 2.267\n d = -0.6395\n e = 1.1558", bmax, log(2*bmax/(ben+bex)-1));
	//sprintf(start, "a = %lf\n b = -0.27\n c = 2.41\n d = -1\n e = 0.23", bmax);
	sprintf(start, "bf = %lf\n bd = %lf\n c0o = -1.3\n c1o = 3.42\n c2o = -1.2\n c3o = 0.23\n c0i = -1.3\n c1i = 3.42\n c2i = -1.2\n c3i = 0.23", bfmax, bdmax);
	//sprintf(start, "a = %lf\n b = %lf\n c = 2.267\n d = -0.6395\n e = 1.1558", bmax, (c0en+c0ex)/2.);
	//sprintf(start, "a = %lf\n b = 0.1455\n c = 2.267\n d = -0.6395\n e = 1.1558", bmax);
	//sprintf(equation, "a/((1+exp(b+c*((x-%lf)/%lf)+d*((x-%lf)/%lf)**2+e*((x-%lf)/%lf)**3))*(1+exp(b+c*((%lf-x)/%lf)+d*((%lf-x)/%lf)**2+e*((%lf-x)/%lf)**3)))", sex, lambda, sex, lambda, sex, lambda, sen, lambda, sen, lambda, sen, lambda);
	//printf("coucou\n");
	sprintf(equationf, "bf/((1+exp(c0i+c1i*((x-%lf)/%lf)+c2i*((x-%lf)/%lf)**2+c3i*((x-%lf)/%lf)**3))*(1+exp(c0o+c1o*((%lf-x)/%lf)+c2o*((%lf-x)/%lf)**2+c3o*((%lf-x)/%lf)**3)))", sfex, lambda, sfex, lambda, sfex, lambda, sfen, lambda, sfen, lambda, sfen, lambda);
	sprintf(equationd, "bd/((1+exp(c0o+c1o*((x-%lf)/%lf)+c2o*((x-%lf)/%lf)**2+c3o*((x-%lf)/%lf)**3))*(1+exp(c0i+c1i*((%lf-x)/%lf)+c2i*((%lf-x)/%lf)**2+c3i*((%lf-x)/%lf)**3)))", sdex, lambda, sdex, lambda, sdex, lambda, sden, lambda, sden, lambda, sden, lambda);
	sprintf(equationtot, "%s+%s", equationf,equationd);
	sprintf(fiteps, "%s_fit_doublet2", textfile);
	//gnuplot_fit_general(track, "1", "2", "a/((1+exp(0.1455+2.267*((x-bex)/l)-0.6395*((x-bex)/l)**2+1.1558*((x-bex)/l)**3))*(1+exp(0.1455+2.267*((ben-x)/l)-0.6395*((ben-x)/l)**2+1.1558*((ben-x)/l)**3)))", "a,ben,bex,l", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", NULL, NULL, fiteps, NULL);
	gnuplot_fit_general(track, "1", "2", equationtot, "bf,bd,c0o,c1o,c2o,c3o,c0i,c1i,c2i,c3i", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", xrange, NULL, fiteps, NULL);
	
}

extern void check_fit_iker_doublet1(char *textfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda)
{
	char track[300], fiteps[300], start[300], equationf[300], equationd[300],xrange[100], equationtot[500];
	int i,nblines, skiplines=8;
	int ien, iex, ifmax, idmax;
	nblines = get_nb_lines_file(textfile);
	
	double s[nblines], b[nblines], sfmax=0, bfmax = 1.e9, sfen, sfex, bfen, bfex, cf0en, cf0ex, sdmax=0, bdmax = -1.e9, sden, sdex, bden, bdex, cd0en, cd0ex;
	//double half_length=0.15;
	FILE *rfile=NULL;
	FILE *wfile;
	
	printf("opening %s...\n", textfile);
	
	rfile = fopen(textfile, "r");
	for(i=0;i<skiplines;i++) newline(rfile);
	sprintf(track, "%s_clean.dat", textfile);
	wfile = fopen(track, "w");
	
	for(i=0;i<nblines-skiplines-1;i++) {
		fscanf(rfile, "%lf      %le", &(s[i]), &(b[i]));
		fprintf(wfile, "%lf	%le\n", s[i], b[i]);
		//printf("%lf	%le\n", s,b);
		if(bfmax>b[i]) {
			bfmax = b[i];
			sfmax = s[i];
			ifmax = i;
		}
		if(bdmax<b[i]) {
			bdmax = b[i];
			sdmax = s[i];
			idmax = i;
		}
	}
	sfex = sfmax+r0*opening_angle_deg_f/2.*PI/180.;
	sfen = sfmax-r0*opening_angle_deg_f/2.*PI/180.;
	sdex = sdmax+r0*opening_angle_deg_d/2.*PI/180.;
	sden = sdmax-r0*opening_angle_deg_d/2.*PI/180.;
	//for(i=imax;i<nblines-skiplines-1;i++) {
	//	if(s[i]>sfex) {
	//		bex = b[i];
	//		c0ex = log(bmax/bex-1);
	//		iex = i;
	//		break;
	//	}
	//}
	//for(i=0;i<imax;i++) {
	//	if(s[i]>sen) {
	//		ben = b[i];
	//		c0en = log(bmax/ben-1);
	//		ien = i;
	//		break;
	//	}
	//}
	fclose(rfile);
	fclose(wfile);
	printf("bfmax = b(%lf) = %lf\n", sfmax, bfmax);
	printf("bdmax = b(%lf) = %lf\n", sdmax, bdmax);
	//printf("ben = b(%lf) = %lf, c0en=%lf\n", sen, ben, c0en);
	//printf("bex = b(%lf) = %lf, c0ex=%lf\n", sex, bex, c0ex);
	// p = c0 + c1*(x/l) + c2*pow(x/l,2) + c3*pow(x/l,3), c0 = 0.1455\n c1 = 2.267\n c2 = -0.6395\n c3 = 1.1558
	sprintf(xrange, "[%lf:%lf]", sfmax-1.,sdmax+1.);
	
	//sprintf(start, "a = %lf\n ben=%lf\nbex=%lf\n l=0.12", bmax, smax-half_length, smax+half_length);
	//sprintf(start, "a = %lf\n b = %lf\n c = 2.267\n d = -0.6395\n e = 1.1558", bmax, log(2*bmax/(ben+bex)-1));
	//sprintf(start, "a = %lf\n b = -0.27\n c = 2.41\n d = -1\n e = 0.23", bmax);
	sprintf(start, "bf = %lf\n bd = %lf\n c0 = -1.3\n c1 = 3.42\n c2 = -1.2\n c3 = 0.23", bfmax, bdmax);
	//sprintf(start, "a = %lf\n b = %lf\n c = 2.267\n d = -0.6395\n e = 1.1558", bmax, (c0en+c0ex)/2.);
	//sprintf(start, "a = %lf\n b = 0.1455\n c = 2.267\n d = -0.6395\n e = 1.1558", bmax);
	//sprintf(equation, "a/((1+exp(b+c*((x-%lf)/%lf)+d*((x-%lf)/%lf)**2+e*((x-%lf)/%lf)**3))*(1+exp(b+c*((%lf-x)/%lf)+d*((%lf-x)/%lf)**2+e*((%lf-x)/%lf)**3)))", sex, lambda, sex, lambda, sex, lambda, sen, lambda, sen, lambda, sen, lambda);
	//printf("coucou\n");
	sprintf(equationf, "bf/((1+exp(c0+c1*((x-%lf)/%lf)+c2*((x-%lf)/%lf)**2+c3*((x-%lf)/%lf)**3))*(1+exp(c0+c1*((%lf-x)/%lf)+c2*((%lf-x)/%lf)**2+c3*((%lf-x)/%lf)**3)))", sfex, lambda, sfex, lambda, sfex, lambda, sfen, lambda, sfen, lambda, sfen, lambda);
	sprintf(equationd, "bd/((1+exp(c0+c1*((x-%lf)/%lf)+c2*((x-%lf)/%lf)**2+c3*((x-%lf)/%lf)**3))*(1+exp(c0+c1*((%lf-x)/%lf)+c2*((%lf-x)/%lf)**2+c3*((%lf-x)/%lf)**3)))", sdex, lambda, sdex, lambda, sdex, lambda, sden, lambda, sden, lambda, sden, lambda);
	sprintf(equationtot, "%s+%s", equationf,equationd);
	sprintf(fiteps, "%s_fit_doublet1", textfile);
	//gnuplot_fit_general(track, "1", "2", "a/((1+exp(0.1455+2.267*((x-bex)/l)-0.6395*((x-bex)/l)**2+1.1558*((x-bex)/l)**3))*(1+exp(0.1455+2.267*((ben-x)/l)-0.6395*((ben-x)/l)**2+1.1558*((ben-x)/l)**3)))", "a,ben,bex,l", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", NULL, NULL, fiteps, NULL);
	gnuplot_fit_general(track, "1", "2", equationtot, "bf,bd,c0,c1,c2,c3", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", xrange, NULL, fiteps, NULL);
	
}

extern void check_fit_iker_doublet4(char *textfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda)
{
	char track[300], fiteps[300], start[300], equationf[300], equationd[300],xrange[100], equationtot[500];
	int i,nblines, skiplines=8;
	int ien, iex, ifmax, idmax;
	nblines = get_nb_lines_file(textfile);
	
	double s[nblines], b[nblines], sfmax=0, bfmax = 1.e9, sfen, sfex, bfen, bfex, cf0en, cf0ex, sdmax=0, bdmax = -1.e9, sden, sdex, bden, bdex, cd0en, cd0ex;
	//double half_length=0.15;
	FILE *rfile=NULL;
	FILE *wfile;
	
	printf("opening %s...\n", textfile);
	
	rfile = fopen(textfile, "r");
	for(i=0;i<skiplines;i++) newline(rfile);
	sprintf(track, "%s_clean.dat", textfile);
	wfile = fopen(track, "w");
	
	for(i=0;i<nblines-skiplines-1;i++) {
		fscanf(rfile, "%lf      %le", &(s[i]), &(b[i]));
		fprintf(wfile, "%lf	%le\n", s[i], b[i]);
		//printf("%lf	%le\n", s,b);
		if(bfmax>b[i]) {
			bfmax = b[i];
			sfmax = s[i];
			ifmax = i;
		}
		if(bdmax<b[i]) {
			bdmax = b[i];
			sdmax = s[i];
			idmax = i;
		}
	}
	sfex = sfmax+r0*opening_angle_deg_f/2.*PI/180.;
	sfen = sfmax-r0*opening_angle_deg_f/2.*PI/180.;
	sdex = sdmax+r0*opening_angle_deg_d/2.*PI/180.;
	sden = sdmax-r0*opening_angle_deg_d/2.*PI/180.;
	//for(i=imax;i<nblines-skiplines-1;i++) {
	//	if(s[i]>sfex) {
	//		bex = b[i];
	//		c0ex = log(bmax/bex-1);
	//		iex = i;
	//		break;
	//	}
	//}
	//for(i=0;i<imax;i++) {
	//	if(s[i]>sen) {
	//		ben = b[i];
	//		c0en = log(bmax/ben-1);
	//		ien = i;
	//		break;
	//	}
	//}
	fclose(rfile);
	fclose(wfile);
	printf("bfmax = b(%lf) = %lf\n", sfmax, bfmax);
	printf("bdmax = b(%lf) = %lf\n", sdmax, bdmax);
	//printf("ben = b(%lf) = %lf, c0en=%lf\n", sen, ben, c0en);
	//printf("bex = b(%lf) = %lf, c0ex=%lf\n", sex, bex, c0ex);
	// p = c0 + c1*(x/l) + c2*pow(x/l,2) + c3*pow(x/l,3), c0 = 0.1455\n c1 = 2.267\n c2 = -0.6395\n c3 = 1.1558
	sprintf(xrange, "[%lf:%lf]", sfmax-1.,sdmax+1.);
	
	//sprintf(start, "a = %lf\n ben=%lf\nbex=%lf\n l=0.12", bmax, smax-half_length, smax+half_length);
	//sprintf(start, "a = %lf\n b = %lf\n c = 2.267\n d = -0.6395\n e = 1.1558", bmax, log(2*bmax/(ben+bex)-1));
	//sprintf(start, "a = %lf\n b = -0.27\n c = 2.41\n d = -1\n e = 0.23", bmax);
	sprintf(start, "bf = %lf\n bd = %lf\n c01 = -1.3\n c11 = 3.42\n c21 = -1.2\n c31 = 0.23\n c02 = -1.3\n c12 = 3.42\n c22 = -1.2\n c32 = 0.23\n c03 = -1.3\n c13 = 3.42\n c23 = -1.2\n c33 = 0.23\n c04 = -1.3\n c14 = 3.42\n c24 = -1.2\n c34 = 0.23", bfmax, bdmax);
	//sprintf(start, "a = %lf\n b = %lf\n c = 2.267\n d = -0.6395\n e = 1.1558", bmax, (c0en+c0ex)/2.);
	//sprintf(start, "a = %lf\n b = 0.1455\n c = 2.267\n d = -0.6395\n e = 1.1558", bmax);
	//sprintf(equation, "a/((1+exp(b+c*((x-%lf)/%lf)+d*((x-%lf)/%lf)**2+e*((x-%lf)/%lf)**3))*(1+exp(b+c*((%lf-x)/%lf)+d*((%lf-x)/%lf)**2+e*((%lf-x)/%lf)**3)))", sex, lambda, sex, lambda, sex, lambda, sen, lambda, sen, lambda, sen, lambda);
	//printf("coucou\n");
	sprintf(equationf, "bf/((1+exp(c02+c12*((x-%lf)/%lf)+c22*((x-%lf)/%lf)**2+c32*((x-%lf)/%lf)**3))*(1+exp(c01+c11*((%lf-x)/%lf)+c21*((%lf-x)/%lf)**2+c31*((%lf-x)/%lf)**3)))", sfex, lambda, sfex, lambda, sfex, lambda, sfen, lambda, sfen, lambda, sfen, lambda);
	sprintf(equationd, "bd/((1+exp(c04+c14*((x-%lf)/%lf)+c24*((x-%lf)/%lf)**2+c34*((x-%lf)/%lf)**3))*(1+exp(c03+c13*((%lf-x)/%lf)+c23*((%lf-x)/%lf)**2+c33*((%lf-x)/%lf)**3)))", sdex, lambda, sdex, lambda, sdex, lambda, sden, lambda, sden, lambda, sden, lambda);
	sprintf(equationtot, "%s+%s", equationf,equationd);
	sprintf(fiteps, "%s_fit_doublet4", textfile);
	//gnuplot_fit_general(track, "1", "2", "a/((1+exp(0.1455+2.267*((x-bex)/l)-0.6395*((x-bex)/l)**2+1.1558*((x-bex)/l)**3))*(1+exp(0.1455+2.267*((ben-x)/l)-0.6395*((ben-x)/l)**2+1.1558*((ben-x)/l)**3)))", "a,ben,bex,l", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", NULL, NULL, fiteps, NULL);
	gnuplot_fit_general(track, "1", "2", equationtot, "bf,bd,c01,c11,c21,c31,c02,c12,c22,c32,c03,c13,c23,c33,c04,c14,c24,c34", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", xrange, NULL, fiteps, NULL);
	
}
	
//0.5/((1+exp(0.1455+2.267*((2-x)/0.12)-0.6395*((2-x)/0.12)**2+1.1558*((2-x)/0.12)**3))*(1+exp(0.1455+2.267*((x-2.5)/0.12)-0.6395*((x-2.5)/0.12)**2+1.1558*((x-2.5)/0.12)**3)))	


extern void write_field_const_radius_across_cell(char *textfile, double r, double z, int nbsteps, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double x,y, th, thstep;
	
	emptyfile(textfile);
	
	thstep = comp_step(0, cell->boun.thmax, nbsteps);
	for(i=0;i<nbsteps;i++) {
		th = i*thstep;
		x = r*cos(th);
		y = r*sin(th);
		write_get_bfield(textfile, x, y, z, cell, add_contribution_comp, YES);
	}
	
}

extern void shift_bfield_tosca(char *readfile, char *writefile, double *bdmax, double *thdmax, double *bfmax)
{
	int i, nblines;
	nblines = get_nb_lines_file(readfile);
	printf("nblines=%i\n", nblines);
	double r[nblines],th[nblines],z[nblines],br[nblines],bth[nblines],bz[nblines];
	FILE *rfile=NULL;
	FILE *wfile;
	
	rfile = fopen(readfile, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	
	wfile = fopen(writefile, "w");
	
	for(i=0;i<nblines;i++) {
		fscanf(rfile, "%le	%le	%le	%le	%le	%le", &r[i],&th[i],&z[i], &br[i],&bth[i], &bz[i]);
	}
	fclose(rfile);
	
	*bdmax=0;
	*thdmax=0; 
	for(i=0;i<nblines/2;i++) {
		if(bz[i]>0) {
			if(bz[i]>*bdmax) {
				*bdmax = bz[i];
				*thdmax = th[i];
			}
		}
	}
	printf("bdmax = %le = Bz(th=%lf)\n", *bdmax, *thdmax);
	*bfmax = bz[0];
	printf("Bfmax = %le\n", *bfmax);
	for(i=nblines/2;i<nblines;i++) {
		fprintf(wfile, "%le	%le	%le	%lf	%lf	%le\n", r[i],th[i]-th[nblines/2],z[i], br[i],bth[i], bz[i]);
	}
	for(i=1;i<nblines/2+1;i++) {
		fprintf(wfile, "%le	%le	%le	%lf	%lf	%le\n", r[i],th[i]+th[nblines/2],z[i], br[i],bth[i], bz[i]);
	}
	
	fclose(wfile);
}

extern void fit_gnuplot_kurns_mr(char *textfile, char *fiteps, double bd0, double cd1, double cd2, double thcd, double ffbd, double bf0, double cf, double ffbf)
{
	char equation[1000], start[1000], bd1[200], bd2[200], bd3[200], bd4[200], bf[200];
	sprintf(bd1, "bd/((1+exp(cd1*(thcd-ffbd-x)))*(1+exp(cd2*(x-(thcd+ffbd)))))"); // D2 of D1 D2 F D3 D4
	sprintf(bd2, "bd/((1+exp(cd2*((30-thcd)-ffbd-x)))*(1+exp(cd1*(x-((30-thcd)+ffbd)))))"); //D3
	sprintf(bd3, "bd/((1+exp(cd2*((-thcd)-ffbd-x)))*(1+exp(cd1*(x-((-thcd)+ffbd)))))"); //D1
	sprintf(bd4, "bd/((1+exp(cd1*(thcd+30-ffbd-x)))*(1+exp(cd2*(x-(thcd+30+ffbd)))))");
	sprintf(bf,  "bf/((1+exp(cf*(15-ffbf-x)))*(1+exp(cf*(x-15-ffbf))))");
	
	sprintf(equation, "%s + %s + %s + %s + %s", bd1,bd2,bd3,bd4,bf);
	//sprintf(equation, "bd/((1+exp(cd1*(thcd-ffbd-x)))*(1+exp(cd2*(x-thcd+ffbd)))) + bd/((1+exp(cd1*((30-(thcd+ffbd))-x)))*(1+exp(cd2*(x-(30-thcd-ffbd))))) + bf/((1+exp(cf*(15-ffbf-x)))*(1+exp(cf*(x-15-ffbf))))");
	//+ bd/((1+exp(cd2*((-dex)-x)))*(1+exp(cd1*(x+den)))) 
	//+ bd/((1+exp(cd1*(30+thcd-ffbd-x)))*(1+exp(cd2*(30+x-(thcd+ffbd))))) 
	
	sprintf(start, "bd = %le\n cd1 = %le\n cd2 = %le\n thcd = %le\n ffbd = %le\n bf = %le\n cf = %le\n ffbf = %le", bd0,cd1,cd2,thcd,ffbd,bf0,cf,ffbf);
	//sprintf(start, "bd = 2.226237e-01\n cd1 = 3.5\n cd2 = 10\n thcd = 6.19\n ffbd = 1.7\n bf = -4.365629e-01\n cf = 3.5\n ffbf = 5.12");
	gnuplot_fit_general(textfile, "2", "6", equation, "bd,cd1,cd2,thcd,ffbd,bf,cf,ffbf", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", NULL, NULL, fiteps, NULL);
	
}

extern void wrapper_fit_field_toscamap_kurri_mr(char *prefix, double r, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	char textfile[200], textfile_shifted[230], fiteps[230];
	double bdmax, thdmax, bfmax;
	
	sprintf(textfile, "%s_r_%.0fcm.dat", prefix, r*100.);
	sprintf(textfile_shifted, "%s_r_%.0fcm_shifted.dat", prefix, r*100.);
	sprintf(fiteps, "%s_r_%.0fcm_shifted", prefix, r*100.);
	printf("writing file %s...\n", textfile);
	
	write_field_const_radius_across_cell(textfile, r, 0, 3001, cell, add_contribution_comp);
	shift_bfield_tosca(textfile, textfile_shifted, &bdmax, &thdmax, &bfmax);
	fit_gnuplot_kurns_mr(textfile_shifted, fiteps, bdmax, 0.6, 1.4, 6.19, 1.9, bfmax, 2.6, 5.6);
}

extern void fit_iker_clean_file3(char *textfile)
{
	char track50[300], track80[300], track100[300], track120[300];
	int i,nblines, skiplines=8;
	nblines = get_nb_lines_file(textfile);
	
	double s, b50, b80, b100, b120;
	
	FILE *rfile=NULL;
	FILE *wfile50, *wfile80, *wfile100, *wfile120;
	
	printf("opening %s...\n", textfile);
	rfile = fopen(textfile, "r");
	for(i=0;i<skiplines;i++) newline(rfile);
	sprintf(track50, "%s_clean50.dat", textfile);
	sprintf(track80, "%s_clean80.dat", textfile);
	sprintf(track100, "%s_clean100.dat", textfile);
	sprintf(track120, "%s_clean120.dat", textfile);
	wfile50 = fopen(track50, "w");
	wfile80 = fopen(track80, "w");
	wfile100 = fopen(track100, "w");
	wfile120 = fopen(track120, "w");
	
	for(i=0;i<nblines-skiplines-1;i++) {
		fscanf(rfile, "%lf      %le      %le      %le      %le", &(s), &(b50), &(b80), &(b100), &(b120));
		fprintf(wfile50, "%lf	%le\n", s, b50);
		fprintf(wfile80, "%lf	%le\n", s, b80);
		fprintf(wfile100, "%lf	%le\n", s, b100);
		fprintf(wfile120, "%lf	%le\n", s, b120);
	}
	fclose(rfile);
	fclose(wfile50);
	fclose(wfile80);
	fclose(wfile100);
	fclose(wfile120);
}

extern void fit_iker_doublet2_clean(char *cleanfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda)
{
	char fiteps[300], start[300], equationf[300], equationd[300], xrange[100], equationtot[500];
	int i, nblines;
	nblines = get_nb_lines_file(cleanfile);
	
	double s, b, sfmax=0, bfmax = 1.e9, sdmax=0, bdmax = -1.e9, sfen, sfex, sden, sdex;
	FILE *rfile=NULL;
	
	printf("opening %s...\n", cleanfile);
	rfile = fopen(cleanfile, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(i=0;i<nblines-1;i++) {
		fscanf(rfile, "%lf      %le", &(s), &(b));
		//printf("%lf	%le\n", s,b);
		if(bfmax>b) {
			bfmax = b;
			sfmax = s;
		}
		if(bdmax<b) {
			bdmax = b;
			sdmax = s;
		}
	}
	fclose(rfile);
	sfex = sfmax+r0*opening_angle_deg_f/2.*PI/180.;
	sfen = sfmax-r0*opening_angle_deg_f/2.*PI/180.;
	sdex = sdmax+r0*opening_angle_deg_d/2.*PI/180.;
	sden = sdmax-r0*opening_angle_deg_d/2.*PI/180.;
	printf("sfen=%lf, sfex=%lf, sden=%lf, sdex=%lf\n", sfen,sfex,sden,sdex);
	printf("bfmax = b(%lf) = %lf\n", sfmax, bfmax);
	printf("bdmax = b(%lf) = %lf\n", sdmax, bdmax);
	sprintf(xrange, "[%lf:%lf]", sfmax-0.65,sdmax+0.7);
	
	//sprintf(start, "bf = %lf\n bd = %lf\n c0o = -1.7\n c1o = 4.3\n c2o = -2.0\n c3o = 0.5\n c0i = -7.5\n c1i = 11.8\n c2i = -6.5\n c3i = 1.8", bfmax, bdmax);
	sprintf(start, "bf = %lf\n bd = %lf\n c0o = -1.7\n c1o = 4.3\n c2o = -2.0\n c3o = 0.5\n c0i = -1.5\n c1i = 4.8\n c2i = -2.5\n c3i = 0.5", bfmax, bdmax);
	sprintf(equationf, "bf/((1+exp(c0i+c1i*((x-%lf)/%lf)+c2i*((x-%lf)/%lf)**2+c3i*((x-%lf)/%lf)**3))*(1+exp(c0o+c1o*((%lf-x)/%lf)+c2o*((%lf-x)/%lf)**2+c3o*((%lf-x)/%lf)**3)))", sfex, lambda, sfex, lambda, sfex, lambda, sfen, lambda, sfen, lambda, sfen, lambda);
	sprintf(equationd, "bd/((1+exp(c0o+c1o*((x-%lf)/%lf)+c2o*((x-%lf)/%lf)**2+c3o*((x-%lf)/%lf)**3))*(1+exp(c0i+c1i*((%lf-x)/%lf)+c2i*((%lf-x)/%lf)**2+c3i*((%lf-x)/%lf)**3)))", sdex, lambda, sdex, lambda, sdex, lambda, sden, lambda, sden, lambda, sden, lambda);
	sprintf(equationtot, "%s+%s", equationf,equationd);
	sprintf(fiteps, "%s_fit_doublet2", cleanfile);
	gnuplot_fit_general(cleanfile, "1", "2", equationtot, "bf,bd,c0o,c1o,c2o,c3o,c0i,c1i,c2i,c3i", start, "lines lc 7 lw 2", "lines dt 1 lc 0 lw 2", "s [m]", "B [T]", "Comsol", "fit", xrange, NULL, fiteps, NULL);
	//gnuplot_fit_general(cleanfile, "1", "2", equationtot, "bf,bd,c0o,c1o,c2o,c3o,c0i,c1i,c2i,c3i", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", NULL, NULL, fiteps, NULL);
}

extern void fit_iker_doublet1_clean(char *cleanfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda)
{
	char fiteps[300], start[300], equationf[300], equationd[300],xrange[100], equationtot[500];
	int i,nblines;
	nblines = get_nb_lines_file(cleanfile);
	
	double s, b, sfmax=0, bfmax = 1.e9, sdmax=0, bdmax = -1.e9, sfen, sfex, sden, sdex;
	FILE *rfile=NULL;
	
	printf("opening %s...\n", cleanfile);
	rfile = fopen(cleanfile, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(i=0;i<nblines-1;i++) {
		fscanf(rfile, "%lf      %le", &(s), &(b));
		//printf("%lf	%le\n", s,b);
		if(bfmax>b) {
			bfmax = b;
			sfmax = s;
		}
		if(bdmax<b) {
			bdmax = b;
			sdmax = s;
		}
	}
	fclose(rfile);
	sfex = sfmax+r0*opening_angle_deg_f/2.*PI/180.;
	sfen = sfmax-r0*opening_angle_deg_f/2.*PI/180.;
	sdex = sdmax+r0*opening_angle_deg_d/2.*PI/180.;
	sden = sdmax-r0*opening_angle_deg_d/2.*PI/180.;
	
	printf("bfmax = b(%lf) = %lf\n", sfmax, bfmax);
	printf("bdmax = b(%lf) = %lf\n", sdmax, bdmax);
	sprintf(xrange, "[%lf:%lf]", sfmax-1.,sdmax+1.);
	
	sprintf(start, "bf = %lf\n bd = %lf\n c0 = -1.3\n c1 = 3.42\n c2 = -1.2\n c3 = 0.23", bfmax, bdmax);
	sprintf(equationf, "bf/((1+exp(c0+c1*((x-%lf)/%lf)+c2*((x-%lf)/%lf)**2+c3*((x-%lf)/%lf)**3))*(1+exp(c0+c1*((%lf-x)/%lf)+c2*((%lf-x)/%lf)**2+c3*((%lf-x)/%lf)**3)))", sfex, lambda, sfex, lambda, sfex, lambda, sfen, lambda, sfen, lambda, sfen, lambda);
	sprintf(equationd, "bd/((1+exp(c0+c1*((x-%lf)/%lf)+c2*((x-%lf)/%lf)**2+c3*((x-%lf)/%lf)**3))*(1+exp(c0+c1*((%lf-x)/%lf)+c2*((%lf-x)/%lf)**2+c3*((%lf-x)/%lf)**3)))", sdex, lambda, sdex, lambda, sdex, lambda, sden, lambda, sden, lambda, sden, lambda);
	sprintf(equationtot, "%s+%s", equationf,equationd);
	sprintf(fiteps, "%s_fit_doublet1", cleanfile);
	gnuplot_fit_general(cleanfile, "1", "2", equationtot, "bf,bd,c0,c1,c2,c3", start, "lines dt 1 lc 7 lw 4", "lines dt 1 lc 0 lw 2", "s [m]", "B [T]", "Comsol", "fit", xrange, NULL, fiteps, NULL);
}

extern void fit_iker_doublet4_clean(char *cleanfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda)
{
	char fiteps[300], start[300], equationf[300], equationd[300], xrange[100], equationtot[500];
	int i,nblines;
	nblines = get_nb_lines_file(cleanfile);
	
	double s, b, sfmax = 0, bfmax = 1.e9, sdmax = 0, bdmax = -1.e9, sfen, sfex, sden, sdex;
	FILE *rfile=NULL;
	
	printf("opening %s...\n", cleanfile);
	rfile = fopen(cleanfile, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(i=0;i<nblines-1;i++) {
		fscanf(rfile, "%lf      %le", &(s), &(b));
		//printf("%lf	%le\n", s,b);
		if(bfmax>b) {
			bfmax = b;
			sfmax = s;
		}
		if(bdmax<b) {
			bdmax = b;
			sdmax = s;
		}
	}
	fclose(rfile);
	sfex = sfmax+r0*opening_angle_deg_f/2.*PI/180.;
	sfen = sfmax-r0*opening_angle_deg_f/2.*PI/180.;
	sdex = sdmax+r0*opening_angle_deg_d/2.*PI/180.;
	sden = sdmax-r0*opening_angle_deg_d/2.*PI/180.;
	
	printf("bfmax = b(%lf) = %lf\n", sfmax, bfmax);
	printf("bdmax = b(%lf) = %lf\n", sdmax, bdmax);
	sprintf(xrange, "[%lf:%lf]", sfmax-1.,sdmax+1.);
	
	sprintf(start, "bf = %lf\n bd = %lf\n c01 = -1.\n c11 = 3.0\n c21 = -1.\n c31 = 0.2\n c02 = -1.7\n c12 = 3.9\n c22 = -2.0\n c32 = 0.4\n c03 = -1.7\n c13 = 3.9\n c23 = -2.0\n c33 = 0.4\n c04 = -1.\n c14 = 3.\n c24 = -1.\n c34 = 0.2", bfmax, bdmax);
	sprintf(equationf, "bf/((1+exp(c02+c12*((x-%lf)/%lf)+c22*((x-%lf)/%lf)**2+c32*((x-%lf)/%lf)**3))*(1+exp(c01+c11*((%lf-x)/%lf)+c21*((%lf-x)/%lf)**2+c31*((%lf-x)/%lf)**3)))", sfex, lambda, sfex, lambda, sfex, lambda, sfen, lambda, sfen, lambda, sfen, lambda);
	sprintf(equationd, "bd/((1+exp(c04+c14*((x-%lf)/%lf)+c24*((x-%lf)/%lf)**2+c34*((x-%lf)/%lf)**3))*(1+exp(c03+c13*((%lf-x)/%lf)+c23*((%lf-x)/%lf)**2+c33*((%lf-x)/%lf)**3)))", sdex, lambda, sdex, lambda, sdex, lambda, sden, lambda, sden, lambda, sden, lambda);
	sprintf(equationtot, "%s+%s", equationf,equationd);
	sprintf(fiteps, "%s_fit_doublet4", cleanfile);
	gnuplot_fit_general(cleanfile, "1", "2", equationtot, "bf,bd,c01,c11,c21,c31,c02,c12,c22,c32,c03,c13,c23,c33,c04,c14,c24,c34", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", xrange, NULL, fiteps, NULL);
}

/*extern void fit_b_iker_doublet2_clean(char *cleanfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda, double c0[2], double c1[2], double c2[2], double c3[2])
{
	char fiteps[300], start[300], equationf[300], equationd[300], xrange[100], equationtot[500];
	int i, nblines;
	nblines = get_nb_lines_file(cleanfile);
	
	double s, b, sfmax=0, bfmax = 1.e9, sdmax=0, bdmax = -1.e9, sfen, sfex, sden, sdex;
	FILE *rfile=NULL;
	
	printf("opening %s...\n", cleanfile);
	rfile = fopen(cleanfile, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(i=0;i<nblines-1;i++) {
		fscanf(rfile, "%lf      %le", &(s), &(b));
		//printf("%lf	%le\n", s,b);
		if(bfmax>b) {
			bfmax = b;
			sfmax = s;
		}
		if(bdmax<b) {
			bdmax = b;
			sdmax = s;
		}
	}
	fclose(rfile);
	sfex = sfmax+r0*opening_angle_deg_f/2.*PI/180.;
	sfen = sfmax-r0*opening_angle_deg_f/2.*PI/180.;
	sdex = sdmax+r0*opening_angle_deg_d/2.*PI/180.;
	sden = sdmax-r0*opening_angle_deg_d/2.*PI/180.;
	
	printf("bfmax = b(%lf) = %lf\n", sfmax, bfmax);
	printf("bdmax = b(%lf) = %lf\n", sdmax, bdmax);
	sprintf(xrange, "[%lf:%lf]", sfmax-1.,sdmax+1.);
	
	sprintf(start, "bf = %lf\n bd = %lf", bfmax, bdmax);
	sprintf(equationf, "bf/((1+exp(%lf+%lf*((x-%lf)/%lf)+%lf*((x-%lf)/%lf)**2+%lf*((x-%lf)/%lf)**3))*(1+exp(%lf+%lf*((%lf-x)/%lf)+%lf*((%lf-x)/%lf)**2+%lf*((%lf-x)/%lf)**3)))", c0[0], c1[0], sfex, lambda, c2[0], sfex, lambda, c3[0], sfex, lambda, c0[0], c1[0], sfen, lambda, c2[0], sfen, lambda, c3[0], sfen, lambda);
	sprintf(equationd, "bd/((1+exp(%lf+%lf*((x-%lf)/%lf)+%lf*((x-%lf)/%lf)**2+%lf*((x-%lf)/%lf)**3))*(1+exp(%lf+%lf*((%lf-x)/%lf)+%lf*((%lf-x)/%lf)**2+%lf*((%lf-x)/%lf)**3)))", c0[1], c1[1], sdex, lambda, c2[1], sdex, lambda, c3[1], sdex, lambda, c0[1], c1[1], sden, lambda, c2[1], sden, lambda, c3[1], sden, lambda);
	sprintf(equationtot, "%s+%s", equationf,equationd);
	sprintf(fiteps, "%s_b_fit_doublet2", cleanfile);
	gnuplot_fit_general(cleanfile, "1", "2", equationtot, "bf,bd", start, "lines lw 3 lc 1", "lines lw 3 lc 2", "s [m]", "B [T]", "Comsol", "fit", xrange, NULL, fiteps, NULL);
}//*/

extern int set_closed_orbit_xxp(struct Particle *part, double eps_clo, struct Lattice *latt, int j, int doyouprintf)
{
	int i, n;
	double xstart, xend, uxstart, uxend, eps;
	double factor = 1.e-3;
	struct Particle test_part;
	
	
	test_part.hat = -2;	// no acceleration, no output
	eps  = 1.;
	n	 = 0.;
	printf("original B0 = %lf\n", latt->cell[0].mpara[j][2]);
	while (eps > eps_clo) {
		test_part = *part;
		xstart = test_part.x;
		uxstart = test_part.ux;
		test_part.hat = -2;	// no acceleration, no output
		for(i = 0; i < FCLO_CROSS; i++) {
			part_cross_latt(&test_part, latt,NULL);
			xend = test_part.x;
			uxend = test_part.ux;
			if(test_part.status != ALIVE) break;
		}
		eps = sqrt((xend - xstart)*(xend - xstart) + (uxend - uxstart)*(uxend - uxstart));
		latt->cell[0].mpara[j][2] *= (1. + sign(xstart-xend)*eps*factor);
		
		if(doyouprintf == YES) {
			//CLRSCR();
			printf("try # %i, eps = %le, xstart = %lf, xend = %lf, newb0 = %lf\n", n, eps, xstart, xend, latt->cell[0].mpara[j][2]);
			//fflush(stdout);
		}
		if(n > FCLO_TRIES || test_part.status != ALIVE) {
			if(doyouprintf == YES) printf("\n!!! Cannot find closed orbit, eps_clo = %le [m]!!!, part. status = %i\n", eps, test_part.status);
			return FALSE;
		}
		n++;
	}
	
	if(doyouprintf == YES) {
		if(print_color==YES) COLOR("1;34");
		printf("\nClosed orbit found, \neps = %le [m]\n", eps);
		if(print_color==YES) COLOR("0");
	}
	return TRUE;
}

extern void launch_part_for_co(struct Particle *part, double *xbest, double *uxbest, double xmin, double xmax, int nbx, double uxmin_deg, double uxmax_deg, int nbux, struct Lattice *latt)
{
	int i,j;
	double xstep, uxstep, uxmin=uxmin_deg*PI/180., uxmax=uxmax_deg*PI/180.;
	double xstart, uxstart, xend, uxend, eps,epsbest;
	struct Particle test_part;
	
	xstep = comp_step(xmin, xmax, nbx);
	uxstep = comp_step(uxmin, uxmax, nbux);
	test_part.hat = -2;
	
	for(i=0;i<nbx;i++) {
		for(j=0;j<nbux;j++) {
			test_part = *part;
			xstart = xmin + i*xstep;
			uxstart = uxmin + j*uxstep;
			test_part.x = xstart;
			test_part.ux = uxstart;
			printf("i=%i,j=%i, (%lf, %le)\t", i,j, test_part.x, test_part.ux);
			part_cross_latt(&test_part, latt, NULL);
			xend = test_part.x;
			uxend = test_part.ux;
			printf("(%lf, %le) ", xend, uxend);
			if(test_part.status != ALIVE) errorstop("part lost!");
			eps = sqrt((xend - xstart)*(xend - xstart) + (uxend - uxstart)*(uxend - uxstart));
			printf("eps = %le\n", eps);
			if(i==0 && j==0) {
				epsbest = eps;
				*xbest = xstart;
				*uxbest = uxstart;
			}
			else if(eps<epsbest) {
				epsbest = eps;
				*xbest = xstart;
				*uxbest = uxstart;
			}
		}
	}
	printf("xbest = %lf, x'best = %le\n", *xbest, asin(*uxbest)*180./PI);
}

extern void codall_shinji_arrange(char *textfile, char *textout)
{
	int i, nblines;
	
	FILE *rfile=NULL;
	FILE *wfile;
	
	nblines = get_nb_lines_file(textfile);
	double xshinji[nblines], yshinji[nblines], z, dum, bx, by, bz[nblines];
	double xmax=-1.e9, xmin=1e9, ymax=-1.e9, ymin=1e9, xmean, ymean, r, th, newx, newy;
	
	rfile = fopen(textfile, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	
	for(i=0;i<nblines;i++) {
		fscanf(rfile, "%le %le %le %lf %le %le %le", &xshinji[i], &yshinji[i], &z, &dum, &bx, &by, &bz[i]);
		xmax = MAX(xmax,xshinji[i]);
		xmin = MIN(xmin,xshinji[i]);
		ymax = MAX(ymax,yshinji[i]);
		ymin = MIN(ymin,yshinji[i]);
	}
	fclose(rfile);
	xmean = (xmin+xmax)/2.;
	ymean = (ymin+ymax)/2.;
	printf("x: %le to %le, y: %le to %le\n", xmin, xmax, ymin, ymax);
	printf("xmean: %le, ymean: %le\n", xmean, ymean);
	wfile = fopen(textout, "w");
	for(i=0;i<nblines;i++) {
		newx = -xmean + xshinji[i];
		newy = -ymean + yshinji[i];
		r = sqrt(newx*newx + newy*newy);
		th = atan_ratio(newx, newy)*180./PI;
		if(th>90.) break;
		fprintf(wfile, "%le	%le	%le	%le	%le\n", newx, newy, r, th, bz[i]);
	}
	fclose(wfile);
	
	
}


extern void fft_kurns(char *textfile, char *text_fft, double time_start, double amp_fft[], int nbpoints, double sample_step_sec)
{
	
	int i, nblines;//, istart;//, flag=0;
	//double amp_fft[2*nbpoints+1];
	double t_temp, amp_temp;
	FILE *rfile=NULL;
	FILE *wfile;
	nblines = get_nb_lines_file(textfile);
	printf("nblines=%i\n", nblines);
	rfile = fopen(textfile, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(i=0;i<5;i++) newline(rfile);
	
	for(i=0;i<nblines-5;i++) {
		fscanf(rfile, "%le	%le", &t_temp, &amp_temp);
		if(t_temp>time_start) {
			//flag=1;
			//istart=i;
			printf("i=%i\n", i);
			amp_fft[1] = amp_temp;
			amp_fft[2] = 0 ;
			break;
		}
	}
	
	for(i=1;i<nbpoints;i++) {
		fscanf(rfile, "%le	%le", &t_temp, &amp_fft[2*i+1]);
		amp_fft[2*i+2] = 0;
	}
	four1(amp_fft,nbpoints,1);
	fclose(rfile);
	
	wfile = fopen(text_fft, "w");
	for(i=0;i<nbpoints/2.;i++) {
		fprintf(wfile, "%le	%le\n", 2*i/sample_step_sec/nbpoints, amp_fft[2*i+1]*amp_fft[2*i+1]);
	}
	fclose(wfile);
}

extern double compute_max_radius_from_trackout_1cell(char *trackout, int doyouprintf)
{
	int i, nblines;
	double s, x, y, z, bx, by, bz, ux, uy, uz, brho, t, r, rmax=0;
	FILE *rfile = NULL;
	
	nblines = get_nb_lines_file(trackout);
	
	rfile = fopen(trackout, "r");
	if(rfile == NULL) errorstop("!!!ERROR in compute_av_radius_from_trackout_1cell: cannot open inputfile!!!");
	
	for(i=0;i<nblines-1;i++) {
		fscanf(rfile, "%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le", &s, &x, &y, &z, &bx, &by, &bz, &ux, &uy, &uz, &brho, &t, &r);
		//printf("r=%le, s=%le\n", r, s);
		rmax = MAX(r, rmax); 
		newline(rfile);
	}
	fclose(rfile);
	//rtot=rtot/nblines;
	if(doyouprintf==YES) printf("max radius %le\n", rmax);
	return rmax;
}


//adjust closed orbit on max r_co by changing B0 (1 cell only!)
extern int adjust_b0_ffag_spi_rmax(struct Lattice *latt, struct Particle *part, double r_co_max, double eps_clo, double eps_r0)
{
	int i, j, k, nb_iterations=50;
	double rnew, rnew_max;
	printf("adjust_B0:\n\n");
	for(i = 0; i < nb_iterations; i++) {
		//printf("i=%i\n",i);
		if(find_closed_orbite_xxp(part, &rnew, &(part->ux), &(part->uy), eps_clo, latt, NO) == TRUE) {
			emptyfile("data/r_max_computation.dat");
			part->x = rnew;
			part->s = 0;
			part_cross_latt(part, latt, "data/r_max_computation.dat");
			rnew_max = compute_max_radius_from_trackout_1cell("data/r_max_computation.dat", NO);
			printf("r_co(s0)=%lf, r_co_max = %lf\n", rnew, rnew_max);
			//tune_calc_matrix(part, &nux, &nuz, &betax, &alphax, &betaz, &alphaz, 1.e-4, 1.e-5, 1.e-4, 1.e-5, latt, part_cross_latt, NO, NULL);
			if (fabs(rnew_max - r_co_max) < eps_r0) {
				//printf("\n\tr0adjust = %.8f [m]\n", latt->cell[0].mpara[0][1]);
				printf("iteration %i, adjusted B0 at r0=%lf: (%le, %le) [T]\n", i, latt->cell[0].mpara[0][1], latt->cell[0].mpara[0][2], latt->cell[0].mpara[1][2]);
				part->s = 0;
				return TRUE;
			}
			for(k=0;k<latt->nbcell;k++) {
				for(j = 0; j < latt->cell[k].nbcomp; j++) {
					latt->cell[k].mpara[j][2] *= pow(((rnew_max-r_co_max) + latt->cell[k].mpara[j][1])/latt->cell[k].mpara[j][1], latt->cell[k].mpara[j][3]);
					if(j==0) printf("new B0=%lf T\n", latt->cell[k].mpara[j][2]);
				}
			}
			part->x = rnew - (rnew_max-r_co_max);
			//part->x = r_co_av;
			//printf("rnew = %lf, r_co = %lf\n", rnew, r_co);
			//for(j = 0; j < latt->cell[0].nbcomp; j++) printf("latt->cell[0].mpara[%i][1] = %lf\n", j, latt->cell[0].mpara[j][1]); 
		}
		else {
			printf("\n \niteration %i, problem in adjust_r0, closed orbit not found!\n", i);
			part->s = 0;
			return FALSE;
		}
	}
	printf("\n \nin adjust_r0 precision not achieved: eps = %le, increase the number of iterations (now %i) or decrease demanded precision\n \n", fabs(rnew - r_co_max), nb_iterations);
	printf("\n\n\t\t\tB0adjust = %lf [m]\n\n\n", latt->cell[0].mpara[0][2]);
	part->s = 0;
	return FALSE;
}

