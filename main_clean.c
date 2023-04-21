/*
 *  main.c
 *  Fixfield
 *
 *  Created by Jean-Baptiste Lagrange on 12/1/09.
 *  Copyright 2009 Kyoto University. All rights reserved.
 *  all questions jean-baptiste.lagrange@stfc.ac.uk
 */


#include "main.h"

int main(int argc, char *argv[]) {
	struct Lattice latt;
	struct Beam beam;
	
	// ************************************************************************************ //
	//									initialisation										//
	// ************************************************************************************ //
	
	/* global variables initialisation */
	sound_option = NO;
	debug = YES; //full output
	max_angle = 0.0; //used in nuSTORM
	doyou_long_boun_sym = YES; //magnetic symmetry at the end of the cell
	write_enge = NO; // debug option to output field (see get_field.c)
	print_color = YES; // colour in terminal
	
	if(print_color==YES) COLOR("46");
	printf("\nFixField code, Tracking code RK4\n");
	printf("__________    __________  _____ __   _____   \n");
	printf("| ____| \\ \\  / | ____| | |  ___| |   |  _ \\  \n");
	printf("| |_  | |\\ \\/ /| |_  | | | |_  | |   | | \\ \\ \n");
	printf("| __| | | )  ( | __| | | |  _| | |   | |  ) )\n");
	printf("| |   | |/ /\\ \\| |   | | | |___| |___| |_/ / \n");
	printf("|_|   |_/_/  \\_|_|   |_| |_____|_____|____/  \n");
	if(print_color==YES) COLOR("0");
	printf("\n");
	
	/* empty output files */
	FILE *shell_read;
	shell_read = popen("mkdir data","r");
	pclose(shell_read);
	emptyfile("data/pickup.dat");
	
	load_lattice(&latt, "inputs/fets/16cells_fd/extract_opt_newc1_qx3410_qz3390.latt");
	load_beam(&beam, "inputs/fets/16cells_fd/proton3mev.beam", &latt, YES);

	// ************************************************************************************ //
	//									main loop(s)										//
	// ************************************************************************************ //
	
	/////////////// Tracking //////////
	emptyfile("data/trackout.dat");
	part_cross_latt(&(beam.part[0]), &latt, "data/trackout.dat");
	//part_oneturn(&(beam.part[0]), &latt, "data/trackout.dat");
	
	char xrange_track[100];
	get_xrange("data/trackout.dat", xrange_track);
	easyplot("data/trackout.dat", "1", "7", "lines lw 3 lc 7", "s [m]", "B vert [T]", xrange_track, NULL, "output/bz.eps", "grid\nset mxtics 5\nset mytics 5");
	//easyplot3p("data/trackout.dat","data/trackout.dat","data/trackout.dat", "1", "5", "1", "6", "1", "7", "lines lw 3 lc 2", "lines lw 3 lc 3", "lines lw 3 lc 7", "horizontal", "longitudinal", "vertical", "s [m]", "B [T]", xrange_track, NULL, "output/b_tot.eps", "grid\nset mxtics 5\nset mytics 5");
	plot_traj(&latt, "data/trackout.dat", "data/aspect.dat", "($3)", "($2)", "($2)", "($1)","lines lc 7 lw 2", "lines lt 1 lc 0 lw 2", "y [m]", "x [m]", NULL, NULL, "output/traj_hor_long.eps", "size ratio -1\nset mxtics 5\nset mytics 5");
	
	/////////////// Tunes  ////////////
	double qx, qz, beta_x, alpha_x, beta_z, alpha_z;
	calc_tune_twiss(&(beam.part[0]), &qx, &qz, &beta_x, &alpha_x, &beta_z, &alpha_z, 1.e-5, 1.e-5, 1.e-5, 1.e-5, &latt, part_cross_latt, YES, NULL,0);//flag: 0=2d, 1=4d decoupled, anything else=4d undecoupled.
	//double qxmin, qxmax, qzmin, qzmax, qxtotmin, qxtotmax, qztotmin, qztotmax;
	//get_tune_range("data/tunepoints.dat", &qxmin, &qxmax, &qzmin, &qzmax, &qxtotmin, &qxtotmax, &qztotmin, &qztotmax);
	//tune_diag("data/tunepoints.dat", 15, qxmin, qxmax, qzmin, qzmax, 4, 4, YES, NO, "output/tune_diag_lattice.eps");
	//tune_diag("data/tunepoints.dat", 15, qxtotmin, qxtotmax, qztotmin, qztotmax, 5, 5, YES, NO, "output/tune_diag_ring.eps");
	
	/////////////// Beta-functions /////
	emptyfile("data/betafunc.dat");
	//double betax0=32.089, alphax0=0, betaz0=31.697, alphaz0=0;
	//compute_betafunc("data/betafunc.dat", 100, 1.e-4, 1.e-4, 1.e-4, 1.e-4, &(beam.part[0]), &latt, betax0, alphax0, betaz0, alphaz0);
	compute_periodic_betafunc("data/betafunc.dat", 100, 1.e-5, 1.e-5, 1.e-5, 1.e-5, &(beam.part[0]), &latt, 0);
	char xrange[100], betarange[100];
	get_xrange("data/betafunc.dat", xrange);
	get_betarange("data/betafunc.dat", betarange);
	write_efb_txtfile("data/efb_betafunc.dat", &(beam.part[0]), &latt);
	easyplot3p("data/betafunc.dat", "data/betafunc.dat", "data/efb_betafunc.dat", "1", "2", "1", "4", "1", "($2*0.5)", "lines lt 1 lw 4 lc 3","lines lt 2 lw 4 lc 1", "lines lt 1 lc 0 lw 1", "horizontal", "vertical", NULL,"s [m]", "{/Symbol b} [m]", xrange, betarange, "output/betafunc-efb.eps", "grid\nset mxtics 5\nset mytics 5\n");
	
	/////////////// Dispersion  ////////
	//emptyfile("data/dispersion.dat");
	//double disp0=0., dispprime0=0.;
	//compute_dispersion_2d("data/dispersion.dat", 60, 1.e-4, 1.e-4, 1.e-4, &(beam.part[0]), &latt, disp0, dispprime0, YES);
	//compute_periodic_dispersion_2d("data/dispersion.dat", 60, 1.e-4, 1.e-4, 1.e-4, &(beam.part[0]), &latt, NO);
	//char disprange[100];
	//get_disprange("data/dispersion.dat", disprange);
	//easyplot_beta_disp("data/betafunc.dat", "data/dispersion.dat", "data/efb_betafunc.dat", "1", "2", "1", "4", "1", "2", "lines lt 1 lw 4 lc 3","lines lt 2 lw 4 lc 1", "lines lt 5 lw 6 lc rgb 'forest-green'", "s [m]", "{/Symbol b} [m]", " {/Symbol h}[m]", xrange, betarange, "[0:1]", "output/beta_disp.eps", "y2tics 0.25 offset -1.5,0\nset grid", 0.5, NO);
	
	// ************************************************************************************ //
	//										terminate										//
	// ************************************************************************************ //
	free_latt(&latt);
	free_beam(&beam);
	
	if(sound_option==YES) say("computation is done");
	if(print_color==YES) COLOR("0");
	printf("\nDONE\n");
	return 0;
}
