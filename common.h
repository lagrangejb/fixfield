/*
 *  common.h
 *  ringdesign
 *
 *  Copyright 2009 Kyoto University. All rights reserved.
 *
 */
#ifndef COMMON
#define COMMON
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include<complex.h>

/* ******************** strutures ******************** */
struct Framework {
	double xc;						//[m] position of the framwork origin in the global framework
	double yc;						// .
	double ae;						//[rad] angle between the framework and the global framework
};

struct Particle {
	struct Framework fwk;			//local framework in which particle coorinates are given
	int status;						//particle status ("alive" or lost)
	int hat;						//Not accelerated particles (and not kicked in thin-dipole):	-1 write related data in trackout.dat,	-2 don't write related data in trackout.dat, -3 compute max angle
									//Accelerated particles:		1 write related data in trackout.dat,	2 don't write related data in trackout.dat, 3 init t0 of RF cavities
	double m0;						//[kg] tracked particle rest mass
	double q;						//[C] tracked particle charge
	//++//unknowns of the differential equations to integrate://++//
	double x;						//[m] particle coordinates in cartesian coordinates
    double y;						//[m]
	double z;						//[m]
	double ux;						//[\] px/ptot (momentum ratio)
	double uy;						//[\] py/ptot
	double uz;						//[\] pz/ptot
	double brho;					//[T.m] particle rigidity
	double t;						//[s] particle time of flight since the origin of time (generally injection of the "on-phase" particle)
	//++//variable of the differential equation://++//
	double s;						//[m] particle length of flight since the new turn
};

struct Beam {
	int npart;						//number of particles in the bunch (>= 1) 
	struct Particle *part;			//array to store all the bunch's particles - memory has to be allocated for this array
};

struct Collimator {
	double rmin;					//[m]
	double rmax;					//[m]
	double zmin;					//[m]
	double zmax;					//[m]
};

struct Boundary {
	double thmax;					//[rad]
	double ymax;					//[m]
};

struct Mapnode {
	double coord[3];				//coord[0,1,2] =[x,y,z] in cartesian maps, [r,th,z] in cylindrical maps
	double b[3];					//b[0,1,2] =[bx,by,bz] in cartesian maps, [br,bth,bz] in cylindrical maps
};

struct Map {
	int nnodes[3];					//number of nodes in x, y, z /or r, th, z dimensions
	int sym;						//YES if median plane symmetry assumed, NO if not.
	double stepsize[3];				//stepsize[0,1,2] =[step_x,step_y,step_z] in cartesian maps, [step_r,step_th,step_z] in cylindrical maps (all in I.S. unites)
	double mapdim[6];				//mapdim[1...6] = [xmin, xmax, ymin, ymax, zmin, zmax] in cartesian maps, [rmin, rmax, thmin, thmax, zmin, zmax]  in cylindrical maps (all in I.S. unites)
	struct Mapnode ***node;			//memory to be allocated
};

struct Cavity {
	char keyword[KEYWORDL];			//cavity name (not used yet)
	double frf;						//[Hz]
	double v0;						//[V]
	double phi0;					//[rad]
	double t0;						//[s]
};

struct Cell {
	char keyword[KEYWORDL];
	double stepsize;				//[m] stepsize
	struct Framework framework;
	double deltar;					//[m]
	struct Collimator collim;
	struct Boundary boun;
	int nbcomp;						//[/] number of components in the Cell
	
	//instrument (pickup of cup)
	int instrutype;					//=NO for no istrument, =PICKUP for a pickup, =CUP for a Faraday cup.
	struct Boundary instru;
	
	//if field is obtained from analytical field modeld: 
	double **mpara;					//components body parameters, memory to be allocated
	double **efben;					//components entrance effective filed boundary (EFB) parameters, memory to be allocated
	double **efbex;					//components exit effective filed boundary (EFB) parameters, memory to be allocated
	
	//3D translation and rotation of each magnetic pole to simulate alignment errors
	int doyou_err;				//switch, error are included in the tracking only if doyou_errtrack == YES (see get_field.c)
	double **alierror;			//memory to be allocated: 
								//alierror[][0] not used, alierror[][1],[2], and [3] define the pole translation (unit = [m,m,m]); 
								//rotation of the magnetic element of an angle theta = alierror[][4] (unit = [rad], right-hand oriented) 
								//is done around and axe passing by the point ([5],[6],[7]) (unit = [m,m,m]) and parallel to the unit vector ([8],[9],[10]). 
	
	//if field is obtained from interpolation in a field map, memory to be allocated:
	struct Map map;
	
	//if the Cell is an RF cavity
	struct Cavity *cav;				//array of "lattice periodicity" number of cavity. memory to be allocated
};

struct Lattice {
	int nbcell;						//number of Cells in the Lattice
	int periodicity;				//periodicity
	int cellnumber;					//when periodicity > 1, it is necessary to identify the cell number for RF phase computation issues
	struct Cell *cell;				//lattice Cells, memory to be allocated
};

/* ******************** global variables (the fewer the better!) ******************** */
int debug;					//NO: no debug output, YES: maximum output (especially init.c)
double max_angle;			//max horizontal angle during tracking.
int doyou_long_boun_sym;	//YES: field symmetry at longitudinal boundary, NO: field=0 outside of cell
int sound_option;			//sound option YES:ON, NO:OFF 
int write_enge;				//flag to write in spi_enge the fields components (in "data/field_spi_enge.dat")
int print_color;			//flag to use color in terminal output
/* ******************** global functions ******************** */
//m_alloc.c
void alloc_beam(struct Beam *beam, int npart);
void free_beam(struct Beam *beam);
void copy_beam(struct Beam *copy, struct Beam *ini_beam);
struct Cell *alloccell(int nbcell);
double **allocmatrix(int n1, int n2);
struct Map allocmap(int nnode_r, int nnode_th, int nnode_z);
struct Cavity *alloccav(int periodicity);
void free_map(struct Map *map, int nnode_r, int nnode_th);
void free_latt(struct Lattice *cell);
void copy_latt(struct Lattice *copy, struct Lattice *ini_latt);
void copy_latt_1period(struct Lattice *reference_latt, struct Lattice *copy_latt);
void copy_latt_n_first_cells(struct Lattice *copy, struct Lattice *ini_latt, int nb_cell);
void copy_cell_map(struct Cell *copycell, struct Cell *inicell);
void gene_alignerror_latt(struct Lattice *reference_latt, struct Lattice *err_latt);

//frameworks.c
void getready_part(struct Particle *part, struct Cell *cell);
void update_part(int awd_answer, struct Particle *part, struct Cell *cell);
void fwk_pt_loctoglob(double *x, double *y, struct Framework *framework);
void fwk_vect_loctoglob(double *x, double *y, struct Framework *framework);
void fwk_pt_globtoloc(double *x, double *y, struct Framework *framework);
void fwk_vect_globtoloc(double *x, double *y, struct Framework *framework);
void fwk_pt_fwktofwk(double *x, double *y, struct Framework *fwk1, struct Framework *fwk2);
void fwk_vect_fwktofwk(double *x, double *y, struct Framework *fwk1, struct Framework *fwk2);
void changefwk_part(struct Particle *part, struct Framework *newfwk);
void find_fwk_exitface(struct Framework *fwk_exitface, struct Cell *cell);
void find_fwk_instru(struct Framework *fwk_exitface, struct Cell *cell);
void latt_move_exit_to_entrance(struct Lattice *latt);
void shitf_fwk(struct Framework *fwk, double deltax);

//toolbox.c
void errorstop(char *error_text);
void newline(FILE *rfile);
void emptyfile(char *filename);
void print_latt_para(char *filename, char *lattname, struct Lattice *latt);
void print_cell_para(char *filename, char *cellname, struct Cell *cell);
void print_beam_para(char *beamname, struct Beam *beam);
void print_part_para(char *partname, struct Particle *part);
void print_get_bfield(double x, double y, double z, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void printf_mat_1dim_complex(char *message, double complex *m, int nb_dim, char *name_mat, int option, char *txtout);
void printf_mat_4d_real(char *message, double m[4][4], char *name_mat, char *txtout);
double round_nb(double value, int nb_decimal);
int test_power2(unsigned int x);
void get_ebg_part(double *e_tot, double *beta, double *gamma, struct Particle *part);
void twissTMtransform(double m[2][2], double twm[3][3]);
void fix_boundingbox(char *txtfilename, int doyouprint);
void say(char *text_say);
void convert_file(char *ori_file, char *final_file, int resolution, int doyouprint);
double comp_step(double val_start, double val_end, int nvalues);
void write_ellipse(char *wfilename, double emitt_pimmmrad, double beta, double alpha, double xcenter, double ycenter);
void write_ellipse2(char *wfilename, double emitt_pimmmrad, double beta, double alpha, double xcenter, double ycenter);
void write_ellipse_mm_mevc(char *wfilename, double emitt_pimmmrad, double beta, double alpha, double xcenter, double ycenter, double pcenter);
int get_nb_lines_file(char *txtfilename);
void get_smax_file(char *txtfilename, double *smax);
void get_betamax_file(char *txtfilename, double *betamax);
void get_betax_z_max_file(char *txtfilename, double *betax_max, double *betaz_max);
void get_dispmax_file(char *txtfilename, double *dispmax);
void get_dispmin_file(char *txtfilename, double *dispmin);
void get_xrange(char *txtfilename, char *xrange);
void get_betarange(char *txtfilename, char *betarange);
void get_disprange(char *txtfilename, char *disprange);
void get_dispmax_file_4d(char *txtfilename, double *dispxmax, double *dispzmax);
void get_dispmin_file_4d(char *txtfilename, double *dispxmin, double *dispzmin);
void get_disprange_4d(char *txtfilename, char *disprangehor, char *disprangevert);
void get_tune_range(char *txtfilename, double *qx_min, double *qx_max, double *qz_min, double *qz_max, double *qx_tot_min, double *qx_tot_max, double *qz_tot_min, double *qz_tot_max);
void get_space_range(char *aspectfile, double *xmin, double *xmax, double *ymin, double *ymax);
void test_field_maxwell(double x, double y, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void test_field_maxwell2(double x, double y, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
int adjust_r0(struct Lattice *latt, struct Particle *part, double r_co, double *r0, double eps_clo);
void reverse_latt(struct Lattice *new_latt, struct Lattice *latt);
void compute_max_angle_muon_ring(struct Particle *part);
double rec_mul_minus(double a, int i);
double rec_mul_plus(double a, int i);
int fit_ellipse1(int *n, double *x, double *y, double *p);
int fit_ellipse2(int *n, double *x, double *y, double *p, double *e, int *m);
void get_betamin_file(char *txtfilename, double *betamin);
void comp_phase_space_coord(double *x, double *xprime, double *z, double *zprime, double test_x, double test_z, double test_ux, double test_uy, double test_uz, double x_ref, double xprime_ref, double z_ref, double zprime_ref);

//init.c
void gene_singlepart(struct Particle *part, double m0, double q, double x, double y, double z,double px, double py, double pz, double t0, int doyouprintf);
void gene_ellibunch(struct Beam *beam, double dt, double dp_ov_p, double dx, double dxprime, double dz, double dzprime, int nsteps[]);
void gene_rectbunch(struct Beam *beam, double dt, double dp_ov_p, double dx, double dxprime, double dz, double dzprime, int nsteps[]);
void gene_elliran_h(struct Beam *beam, double dt, double dp_ov_p, double emit_x, double twiss_betax, double twiss_alphax, double emit_z, double twiss_betaz, double twiss_alphaz, int nparts, long seed);
void gene_elligasdev_6d(struct Beam *beam, double dt, double dp_ov_p, double emit_x, double twiss_betax, double twiss_alphax, double emit_z, double twiss_betaz, double twiss_alphaz, int nparts, long seed);
void gene_elligasdev_4d(struct Beam *beam, double dt, double dp_ov_p, double emit_x, double twiss_betax, double twiss_alphax, double emit_z, double twiss_betaz, double twiss_alphaz, int nparts, long seed);
void load_beam(struct Beam *beam, char *inputfilename, struct Lattice *latt, int doyouprint);
void load_ao_beam(struct Beam *beam, char *inputfilename, double x_center_pipe);
void load_lattice(struct Lattice *latt, char *inputfilename);
int test_cell_rect_vffa_bend_cell(struct Cell *cell);
int test_cell_vffa(struct Cell *cell);
int test_cell_spiangle(struct Cell *cell);
int test_cell_map(struct Cell *cell);
void gene_alignerror_latt_old(struct Lattice *reference_latt, struct Lattice *err_latt, long *idum, double rms_shift_error, double rms_twist_error, double r_ref_twist);
void add_error_pole(struct Cell *cell, int npole, double xshift, double yshift, double zshift, double angle_rot, double rot_point_x, double rot_point_y, double rot_point_z, double rot_vect_x, double rot_vect_y, double rot_vect_z);
void comp_error_latt_vffa(struct Lattice *err_latt, long seed, double rms_shift_error, double rms_twist_error);
void comp_m_error_latt_vffa(struct Lattice *err_latt, long seed, double rms_error);
void gene_latt_1cell(struct Lattice *new_latt, struct Cell *cell);

//plot.c
void aspectview(char *filename, struct Lattice *latt);
void aspect_fringe_view(char *filename_en, char *filename_ex, struct Lattice *latt);
void easyplot(char *txtfilename, char *xcolumn, char *ycolumn, char *with, char *xlabel,char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void easyplot2(char *txtfilename1, char *txtfilename2, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *with1, char *with2, char *title1, char*title2, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void easyplot3(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *xcolumn1,char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *with1, char *with2,char *with3, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void easyplot3p(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *with1, char *with2, char *with3, char *title1, char *title2, char *title3, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void easyplot4(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *xcolumn4, char *ycolumn4, char *with1, char *with2, char *with3, char *with4, char *title1, char *title2, char *title3, char *title4, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void easyplot5(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *txtfilename5, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *xcolumn4, char *ycolumn4, char *xcolumn5, char *ycolumn5, char *with1, char *with2, char *with3, char *with4, char *with5, char *title1, char *title2, char *title3, char *title4, char *title5, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void easyplot6(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *txtfilename5, char *txtfilename6, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *xcolumn4, char *ycolumn4, char *xcolumn5, char *ycolumn5, char *xcolumn6, char *ycolumn6, char *with1, char *with2, char *with3, char *with4, char *with5, char *with6, char *title1, char *title2, char *title3, char *title4, char *title5, char *title6, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void easyplot8(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *txtfilename5, char *txtfilename6,  char *txtfilename7, char *txtfilename8,char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *xcolumn4, char *ycolumn4, char *xcolumn5, char *ycolumn5, char *xcolumn6, char *ycolumn6, char *xcolumn7, char *ycolumn7, char *xcolumn8, char *ycolumn8, char *with1, char *with2, char *with3, char *with4, char *with5, char *with6,  char *with7, char *with8,char *title1, char *title2, char *title3, char *title4, char *title5, char *title6, char *title7, char *title8, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void easyplot2_left_right(char *txtfilename1, char *txtfilename2, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *with1, char *with2, char *xlabel, char *ylabel, char *y2label, char *xrange, char *yrange, char *y2range, char *psfilename, char *setoption);
void easyplot3_left_right(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *with1, char *with2, char *with3, char *title1, char *title2, char *title3, char *xlabel, char *ylabel, char *y2label, char *xrange, char *yrange, char *y2range, char *psfilename, char *setoption);
void easyplot4_left_right(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *xcolumn4, char *ycolumn4, char *with1, char *with2, char *with3, char *with4, char *title1, char *title2, char *title3, char *title4, char *xlabel, char *ylabel, char *y2label, char *xrange, char *yrange, char *y2range, char *psfilename, char *setoption);
//void easyplot3dmap(char *txtfilename, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_min, double dens_max);
void easyplot3dmap(char *txtfilename, char *xlabel, char *ylabel, char *cblabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_min, double dens_max);
void easyplot3dmap_gen(char *txtfilename, char *xcolumn, char *ycolumn, char *ccolumn, char *xlabel, char *ylabel, char *cblabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_min, double dens_max);
void easyplot3dmap_plusfile(char *txtfile_2d, char *xcolumn_2d, char *ycolumn_2d, char *with_2d, char *txtfilename, char *xlabel, char *ylabel, char *cblabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_min, double dens_max);
void easyplot3dmap_plus2files(char *txtfile_2d, char *xcolumn_2d, char *ycolumn_2d, char *with_2d, char *txtfile_2d2, char *xcolumn_2d2, char *ycolumn_2d2, char *with_2d2, char *txtfilename, char *xlabel, char *ylabel, char *cblabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_min, double dens_max);
void easyplot3dmap_2d(char *filemap, char *file2d, char *para2d, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_min, double dens_max);
//void easyplot3d(char *txtfilename, char *with, char *xlabel, char *ylabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption);
void easyplot3d(char *txtfilename, char *xcolumn, char *ycolumn, char *zcolumn, char *with, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption);
void easyplot2_3d(char *txtfilename1, char *txtfilename2, char *xcolumn1, char *ycolumn1, char *zcolumn1, char *xcolumn2, char *ycolumn2, char *zcolumn2, char *with1, char *with2, char *title1, char *title2, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption);
void easyplot3_3d(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *xcolumn1, char *ycolumn1, char *zcolumn1, char *xcolumn2, char *ycolumn2, char *zcolumn2, char *xcolumn3, char *ycolumn3, char *zcolumn3, char *with1, char *with2, char *with3, char *title1, char *title2, char *title3, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption);
void easyplot4_3d(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *xcolumn1, char *ycolumn1, char *zcolumn1, char *xcolumn2, char *ycolumn2, char *zcolumn2, char *xcolumn3, char *ycolumn3, char *zcolumn3, char *xcolumn4, char *ycolumn4, char *zcolumn4, char *with1, char *with2, char *with3, char *with4, char *title1, char *title2, char *title3, char *title4, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption);
void easyplot5_3d(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *txtfilename5, char *xcolumn1, char *ycolumn1, char *zcolumn1, char *xcolumn2, char *ycolumn2, char *zcolumn2, char *xcolumn3, char *ycolumn3, char *zcolumn3, char *xcolumn4, char *ycolumn4, char *zcolumn4, char *xcolumn5, char *ycolumn5, char *zcolumn5, char *with1, char *with2, char *with3, char *with4, char *with5, char *title1, char *title2, char *title3, char *title4, char *title5, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption);
void easyplot6_3d(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *txtfilename5, char *txtfilename6, char *xcolumn1, char *ycolumn1, char *zcolumn1, char *xcolumn2, char *ycolumn2, char *zcolumn2, char *xcolumn3, char *ycolumn3, char *zcolumn3, char *xcolumn4, char *ycolumn4, char *zcolumn4, char *xcolumn5, char *ycolumn5, char *zcolumn5, char *xcolumn6, char *ycolumn6, char *zcolumn6, char *with1, char *with2, char *with3, char *with4, char *with5, char *with6, char *title1, char *title2, char *title3, char *title4, char *title5, char *title6, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption);
void easyplot7_3d(double x0, double z0, char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *txtfilename5, char *txtfilename6, char *txtfilename7, char *xcolumn1, char *ycolumn1, char *zcolumn1, char *xcolumn2, char *ycolumn2, char *zcolumn2, char *xcolumn3, char *ycolumn3, char *zcolumn3, char *xcolumn4, char *ycolumn4, char *zcolumn4, char *xcolumn5, char *ycolumn5, char *zcolumn5, char *xcolumn6, char *ycolumn6, char *zcolumn6, char *xcolumn7, char *ycolumn7, char *zcolumn7, 
	char *with1, char *with2, char *with3, char *with4, char *with5, char *with6, char *with7, char *title1, char *title2, char *title3, char *title4, char *title5, char *title6, char *title7, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption);
void easyplot_beta_disp(char *txtfilename1, char *txtfilename3, char *txt_filename_efb, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *with1, char *with2, char *with3, char *xlabel, char *ylabel, char *y2label, char *xrange, char *yrange, char *y2range, char *psfilename, char *setoption, double efb_scale, int title);
void write_efb_txtfile(char *filename, struct Particle *reference, struct Lattice *latt);
void plot_twiss_gamma(char *filename, char *with1, char *with2, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void plot_twiss_gamma_w_efb(char *filename, char *file_efb, double efb_scale, char *with1, char *with2, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void plot_one_histogram(char *txtfilename, double xmin, double xmax, double nb_bin, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void plot_two_histogram(char *txtfilename, char *title1, char *txtfilename2, char *title2, double xmin, double xmax, double nb_bin, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void plot_three_histogram(char *txtfilename, char *title1, char *txtfilename2, char *title2, char *txtfilename3, char *title3, double xmin, double xmax, double nb_bin, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void plot_four_histogram(char *txtfilename, char *title1, char *txtfilename2, char *title2, char *txtfilename3, char *title3, char *txtfilename4, char *title4, double xmin, double xmax, double nb_bin, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void plot_five_histogram(char *txtfilename, char *title1, char *txtfilename2, char *title2, char *txtfilename3, char *title3, char *txtfilename4, char *title4, char *txtfilename5, char *title5, double xmin, double xmax, double nb_bin, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void plot_traj(struct Lattice *latt, char *trackout, char *aspectfile, char *xcolumn_out, char *ycolumn_out, char *xcolumn_asp, char *ycolumn_asp, char *with_out, char *with_asp, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void gnuplot_fit(char *txtfilename1, char *xcolumn1, char *ycolumn1, char *a, char *b, char *c, char *with1, 
				 char *with2, char *xlabel, char *ylabel, char *title1, char *title2, char *xrange, char *yrange, char *psfilename, char *setoption);
void gnuplot_fit_general(char *txtfilename, char *xcolumn1, char *ycolumn1, char *fit_function_of_x, char *variables_name, char *init_variables, char *with1, 
 char *with2, char *xlabel, char *ylabel, char *title1, char *title2, char *xrange, char *yrange, char *psfilename, char *setoption);
void tune_diag(char *textfilename, int supersym, double qxmin, double qxmax, double qzmin, double qzmax, int maxorder_s, int maxorder_ns, int doyounorm, int doyouskew, char *output_name);
void gene_elligasdev_h(struct Beam *beam, double dt, double dp_ov_p, double sigmax, double sigmaxp, double twiss_alphax, double twiss_betax, int nparts, long seed);
void plotfunction(char *function, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *title, char *setoption);
void gnuplot_3d_scatered_points(char *txtfilename, char *xcolumn, char *ycolumn, char *cbcolumn, char *xlabel, char *ylabel, char *cblabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_max);

//numerical.c
void rkdrive(struct Particle *part, int *nsteps, double stepsize,void (*derivs)(struct Particle*, double[], double[], struct Cell*, void(*)(double,double,double,double*,double*,double*,struct Cell*,int)),struct Cell* cell, char *txtfile, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void four1(double data[], unsigned long nn, int isign);
int sign(double a);
double atan_ratio(double y, double x);
double bessj0(double x);
double bessj1(double x);
double bessy0(double x);
double bessy1(double x);
double bessel_j(double alpha, double x, int n);
double bessel_y(double alpha, double x, int n);
int overflow_test(double a);
int sanity_test_number(double a);
int factorial(int n);
void mmprod2(double fst[2][2], double sec[2][2], double mul[2][2]);
void mmprod4(double fst[4][4], double sec[4][4], double mul[4][4]);
void mvprod3(double v[], double a[][3], double u[]);
void mvprod4(double v[], double a[][4], double u[]);
void mvprod5(double v[], double a[][5], double u[]);
double matrix_det_4d(double mat[4][4], int n);
double matrix_det_5d(double mat[5][5], int n);
void matrix_transpose_4d(double m[4][4], double tr[4][4], int r);
void matrix_cofactor_4d(double num[4][4], double fac[4][4], int f);
int matrix_inverse_4d(double num[4][4], double inverse[4][4]);
void zgeev(char* jobvl, char* jobvr, int* n, double complex *a, int* lda, double complex *w, double complex *vl, int* ldvl, double complex *vr, int* ldvr, double complex *work, int* lwork, double complex *rwork, int* info);
int zgetrf_(int *m, int *n, double complex *a, int *lda, int *ipiv, int *info);
int zgetri_(int *n, double complex *a, int *lda, int *ipiv, double complex *work, int *lwork, int *info);
void matrix_complex_inverse(double complex *A, double complex *invA, int n);
int invmat(double matrix[5][5], int nfree);
double inimat( double s[5][5], double rl[5], double *x, double *y, int n, double *p, double *e, int ip[5], int nfree);
int inivec( double s[5][5], double s1[5][5], double rl[5], double labda, double *q, double *p, double *e, double *x, double *y, int n, int ip[5], int nfree);
void zge_transpose(double complex *Transposed, double complex *M ,int n);
void matrix_complex_eigensystem(double complex *eigenvectorsVR, double complex *eigenvaluesW, double complex *A, int N);
double comp_print_det_matrix_1dim(char *matrixname, double complex *m_1dim, int n, char *txtout);
int compute_eigenvector_parzen(double complex *eigenvectorsVR, double complex *eigenvaluesW, double complex *m_1dim);
int decouple_matrix_parzen(double m[4][4], double decoup_m[4][4], double decoup_to_coup[4][4], double output[6], char *txtout);
float ran1(long *idum);
float gasdev(long *idum);

//track.c
int find_closed_orbite_x(struct Particle *part, double *x_clo, double eps_clo, struct Lattice *latt, int doyouprintf);
int find_closed_orbite_xz(struct Particle *part, double *x_clo, double *z_clo, double eps_clo, struct Lattice *latt, int doyouprintf);
int find_closed_orbite_xxp(struct Particle *part, double *x_clo, double *ux_clo, double *uy_clo, double eps_clo, struct Lattice *latt, int doyouprintf);
int find_closed_orbite_xxp_zzp(struct Particle *part, double *x_clo, double *ux_clo, double *uy_clo, double *z_clo, double *uz_clo, double eps_clo, struct Lattice *latt, int doyouprintf);
void apply_alierror_comp(struct Cell *cell, int i, double x, double y, double z, double *xwerr, double *ywerr, double *zwerr);
int acceptancex(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double axmin, double axmax, double axstep, double az);
int acceptancez(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double azmin, double azmax, double azstep, double ax);
int acceptancezprime(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double zprimemin, double zprimemax, double zprimestep, double ax);
void acceptancex_auto(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double axmin, double *axmax, double axstep, double az);
void acceptancez_auto(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double azmin, double *azmax, double azstep, double ax);
int track_n_turns_amp(double x[],double z[], double ux[], double uy[], double uz[], struct Particle *reference, struct Lattice *latt, int nbpass, double ax, double az);
void beam_oneturn(struct Beam *bunch, struct Lattice *latt, char *txtfile);
void beam_cross_latt(struct Beam *bunch, struct Lattice *latt, char *txtfile);
void part_oneturn(struct Particle *part, struct Lattice *latt, char *txtfile);
void part_cross_latt(struct Particle *part, struct Lattice *latt, char *txtfile);
int part_cross_cell(struct Particle *part, struct Cell *cell, char *txtfile);
int part_callrk(struct Particle *part, struct Cell *cell, char *txtfile);
void part_cross_thingap(struct Particle *part, struct Cavity *cav);
void part_thindipole(struct Particle *part, struct Cell *cell);
void part_collimator(struct Particle *part, struct Cell *cell);
void part_float_faraday(struct Particle *part, struct Cell *cell);
void derivs(struct Particle *part, double dfds[], double bprobe[], struct Cell *cell,void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));

//get_field.c
int get_bfield(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell,void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void add_field_comp_FFAGradial_he(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGradial_lin(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGradial_2ndorder(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGradial_purescale(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGspiral_he(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGspiral_lin(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGspiral_enge(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGspiral_fullenge(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGstraight_he (double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGstraight_lin(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGstraight_purescale(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGstraighttilt_he (double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGstraighttilt_lin(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_rect_he(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_rect_he_str(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_rect_lin(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_rect_lin_str(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_rect_he(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_rect_enge(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_rect_enge_add(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_rect_str_enge(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_rect_enge_separate_func(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_sect_lin(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_sect_he(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_sect_enge(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_drift(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_realbend(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_sector_purebend(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_sector_purebend_lin(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_quad_he (double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_quad_lin(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_quad_enge(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGtilt_he(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGtilt_lin(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void get_bfield_cylindrical_map(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void get_bfield_cartesian_map(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void get_bfield_cylmap(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void get_bfield_cylmap_nosym(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void get_bfield_cartmap(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void get_bfield_cartmap_nosym(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_FFAGstraight_add(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_rect_bx_add(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_rect_atan(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
void add_field_comp_VFFA_rect_str_atan(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);
//void add_multiple_poles(double x, double y, double z, double *bx, double *by, double *bz, struct Cell *cell, int ic);

//lin_para.c
void tune_calculate_fft(struct Particle *reference, double *nux, double *nuz, int tune_turn, double amp_x, double amp_z, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*,char*), int doyouprintf, char *output_name);
void tune_calc_matrix_no_period(struct Particle *reference, double *qx, double *qz, double amp_x, double amp_xprime, double amp_z, double amp_zprime, double betax0, double alphax0,  double betaz0, double alphaz0, double *betaxfin, double *alphaxfin, double *betazfin, double *alphazfin, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, char *output_name);
void tune_calc_matrix(struct Particle *reference, double *qx, double *qz, double *twiss_bx, double *twiss_ax, double *twiss_bz, double *twiss_az, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*,char*), int doyouprintf, char *output_name);
double calc_phase_advance(double m[2][2], double beta0, double alpha0);
double calc_tune1(double m[2][2]);
double calc_tune2(double m[2][2]);
int calc_phadv_twiss_hor(struct Particle *reference, double *nux, double *qx, double *betaf, double *alphaf, double amp_x, double amp_xprime, double betax0, double alphax0, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf);
int calc_phadv_twiss_vert(struct Particle *reference, double *nuz, double *qz, double *betaf, double *alphaf, double amp_z, double amp_zprime, double betaz0, double alphaz0, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf);
int calc_tune_twiss_4d(struct Particle *reference, double *nux, double *qx, double *nuz, double *qz, double *betax, double *alphax, double *betaz, double *alphaz, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, int together_decoupled);
int calc_phadv_twiss(struct Particle *reference, double *qx, double *qz, double betax0, double alphax0, double betaz0, double alphaz0, double amp_x, double amp_xprime, double amp_z, double amp_zprime, double *betaxf, double *alphaxf, double *betazf, double *alphazf, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, char *outfile);
int calc_tune_twiss_hor(struct Particle *reference, double *nux, double *qx, double *betax, double *alphax, double amp_x, double amp_xprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf);
int calc_tune_twiss_vert(struct Particle *reference, double *nuz, double *qz, double *betaz, double *alphaz, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf);
int calc_tune_twiss_4d(struct Particle *reference, double *nux, double *qx, double *nuz, double *qz, double *betax, double *alphax, double *betaz, double *alphaz, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, int together_decoupled);
int calc_tune_twiss(struct Particle *reference, double *qx, double *qz, double *betax, double *alphax, double *betaz, double *alphaz, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, char *outfile,int together_decoupled);
void phase_adv_scan_woco(char *outfilename, struct Lattice *latt, struct Particle *reference_part, double emin_ev, double emax_ev, int nbstep, double amp);
void phase_adv_scan(char *outfilename, struct Lattice *latt, struct Particle *reference_part, double eps_clo, double emin_ev, double emax_ev, int nbstep, double amp);
double find_s(struct Particle *reference, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*,char*));
void calc_twiss(double m[2][2], double *beta1, double *alpha1, double beta0, double alpha0);
void calc_periodic_twiss(double m[2][2], double nu, double *beta, double *alpha);
void get_periodic_twiss(double *betax0, double *alphax0, double *betaz0, double *alphaz0, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Particle *reference, struct Lattice *latt, int together_decoupled);
void get_twissx_atinstru(double *betax, double *alphax, double beta0, double alpha0, double amp_x, double amp_xprime, struct Particle *reference, struct Lattice *latt);
void get_twissz_atinstru(double *betaz, double *alphaz, double beta0, double alpha0, double amp_x, double amp_xprime, struct Particle *reference, struct Lattice *latt);
void get_twiss_atinstru(double *betax, double *alphax, double *betaz, double *alphaz, double betax1, double alphax1, double betaz1, double alphaz1, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Particle *reference, struct Lattice *latt, int together_decoupled);
void compute_betafunc(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime,struct Particle *reference, struct Lattice *latt, double betax0, double alphax0, double betaz0, double alphaz0, int together_decoupled);
void compute_betafunc_cell(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Particle *reference, struct Cell *cell, double betax0, double alphax0, double betaz0, double alphaz0, int together_decoupled);
void compute_betafunc_latt(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Particle *reference, struct Lattice *latt, double betax0, double alphax0, double betaz0, double alphaz0, int together_decoupled);
void compute_periodic_betafunc(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Particle *reference, struct Lattice *latt, int together_decoupled);
void compute_periodic_dispersion_2d(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_dpovp,struct Particle *reference, struct Lattice *latt, int hor);
void compute_dispersion_2d(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_dpovp,struct Particle *reference, struct Lattice *latt, double disp0, double dispprime0, int hor);
void get_periodic_dispersion_2d(double *disp0, double *dispprime0,  double amp_x, double amp_xprime, double amp_dpovp, struct Particle *reference, struct Lattice *latt, int hor);
void compute_periodic_dispersion_4d(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime, double amp_dpovp, struct Particle *reference, struct Lattice *latt);
void compute_dispersion_4d(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime, double amp_dpovp, struct Particle *reference, struct Lattice *latt, double dispx0, double dispprimex0, double dispz0, double dispprimez0);
void get_periodic_dispersion_4d(double *dispx0, double *dispprimex0, double *dispz0, double *dispprimez0, double amp_x, double amp_xprime, double amp_z, double amp_zprime, double amp_dpovp, struct Particle *reference, struct Lattice *latt);
void find_x(double *x, double *bias, struct Particle *part_ref, struct Particle *part);
void find_z(double *z, double *bias, struct Particle *part_ref, struct Particle *part);
int get_matrix_hor_vert(double mh[2][2], double mv[2][2], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, int together_decoupled);
int get_trmatrix_h(double mh[2][2], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*,char*), int doyouprintf);
int get_trmatrix_v(double mv[2][2], struct Particle *reference, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*,char*), int doyouprintf);
int get_trmatrix_h3x3(double mh[3][3], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, double amp_dpovp, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*,char*), int doyouprintf);
int get_trmatrix_v3x3(double mv[3][3], double *band_bias, struct Particle *reference, double amp_z, double amp_zprime, double amp_dpovp, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*,char*), int doyouprintf);
int get_matrix_firstorder_4d(double m[4][4], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf);
int get_matrix_secondorder_4d(double m[4][4], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf);
int get_matrix_firstorder_5d(double m[5][5], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, double amp_z, double amp_zprime, double amp_dpovp, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf);
int get_matrix_secondorder_5d(double m[5][5], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, double amp_z, double amp_zprime, double amp_dpovp, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf);
int symplectify_2d(double m[2][2], double symplec_m[2][2]);
int symplectify_4d(double m[4][4], double symplec_m[4][4]);

//othercodes_tools.c (other codes tricks)
void plot_magnetic_field_mice(char *name1, char *name2);
void plot_magnetic_field_sam(char *name1, char *name2);
void plot_mice_inspec(char *txtfilename, char *txt_mom1, char *txt_mom2, char *txt_mom3, char *outname, double mom_ref, double mom_spread);
void delta_beta_mice(char *txtfile_center, char *txtfile_off, char *outputfile);
void change_input_zgoubi(char *txtfilename, char *output, int nblines);
void change_input_zgoubi_bfield(char *txtfilename, char *output);
void change_input_zgoubi_acceptance(char *txtfilename, char *output);
void map_opal_circ(struct Cell *cell, char *mapfilename, double rstart, double rend, int nbr, double thstart, double thend, int nbth, double zstart, double zend, int nbz, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void map_opal_cart(struct Cell *cell, char *mapfilename, double xstart, double xend, int nbx, double ystart, double yend, int nby, double zstart, double zend, int nbz, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void adjust_b0_opal(double x, double y, double z, double bz_opal, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double precision, int doyouf);
void opal_magnet_position_angle(char *opal_file, char *opal_gnuplot);
void compare_field_opal_vffa(char *file_dif_prefix, double z_1, struct Cell *cell1, void(*add_contribution_comp1)(double,double,double,double*,double*,double*,struct Cell*,int), double z_2, struct Cell *cell2, void(*add_contribution_comp2)(double,double,double,double*,double*,double*,struct Cell*,int));
void compare_field_map_opal(char *file_opal, char *fileout, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void compare_field_map_opal2(char *file_opal, char *fileout, struct Cell *cell, double z, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double xmin, double xmax, double xstep, double ymin, double ymax, double ystep, double zmin, double zmax, double zstep);
void create_trackout_from_opal_file(char *opal_ori, char *trackoutfile, double rmin, double rmax, double th0_deg, double thmax_deg);
void adjust_map_opal(char *opal_origin_file, char *file_for_fixfield, int nb_skipped_lines, int nbx, int nby, int nbz);
void plot_field_map_opal(char *mapname);
void para_vffa_fodo_scode_to_fixfield(double cell_length, double theta_cell_deg, double magnet_length, double small_drift, double shift, double tilt_d_deg, double tilt_f_deg, double b0d, double b0f, double fringe_add, double m, double z0, int order_expansion);
void para_vffa_triplet_scode_to_fixfield(double cell_length, double theta_cell_deg, double magnet_length, double small_drift, double shift, double tilt_f_deg, double b0d, double b0f, double fringe_add, double m, double z0, int order_expansion);
void para_vffa_triplet_scode_to_fixfield2(double cell_length, double theta_cell_deg, double magnet13_length, double small_drift, double magnet2_length, double shift, double tilt_mag1_deg, double b0f, double b0d, double fringe_add, double m, double z0, int order_expansion);
void para_vffa_triplet_scode_to_fixfield3(double cell_length, double theta_cell_deg, double magnet13_length, double small_drift, double magnet2_length, double shiftf, double shiftd, double tilt_mag1_deg, double b0f, double b0d, double fringe_add, double m, double z0, int order_expansion);
void para_vffa_triplet_scode_to_fixfield_add(double cell_length, double theta_cell_deg, double magnet13_length, double small_drift, double magnet2_length, double shiftf, double shiftd, double tilt_mag1_deg, double b0f, double b0d, double fringe_add, double m, double z0, int order_expansion);
int find_b0_add_to_mult(double b0_add, double efbe, double efbs, double lambda_add, double *b0_mult, double prec);
double find_b0_add_to_mult2(double b0_add, double efbe, double efbs, double lambda_add);
void add_theta_opal_file(struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), char *opal_ori, char *newfile);


//bricabrac.c
void int_bz_liney_cell(double *intbdl, double *intbdlabs, double x, double ymin, double ymax, int nstep, struct Cell *cell,void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void int_bz_radius_cell(double *intbdl, double *intbdlabs, double r, double thmin, double thmax, int nstep, struct Cell *cell,void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void int_bz_max_min_radius_cell(double *intbdl, double *intbdlabs, double *bmax, double *bmin, double r, double thmin, double thmax, int nstep, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void m_loc_bl_cell(char *outfilename, double ymin, double ymax, int nstep_y, double xmin, double xmax, int nstep_x, int nxref, double mref, struct Cell *cell,void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
//void k_loc_bl_cell_direct(char *outfilename, double thmin, double thmax, int nstep_th, double rmin, double rmax, int nstep_r, int nrref, double kref, struct Cell *cell,void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void k_loc_bl_cell(char *outfilename, double thmin, double thmax, int nstep_th, double rmin, double rmax, int nstep_r, int nrref, double kref, struct Cell *cell,void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void k_loc_cell2(char *outfilename, double thmin, double thmax, int nstep_th, double rmin, double rmax, int nstep_r, double deltar, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
double comp_k_value_local(double r, double deltar, double bdeltar, double b);
void m_loc_bl_latt(char *outfilename, int nstep_y, double xmin, double xmax, int nstep_x, int nxref, double mref, struct Lattice *latt,void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void m_loc_y_cell(char *outfilename, double y, double xmin, double xmax, int nstep_x, int nxref, double mref, struct Cell *cell,void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void scan_err_b0(char *outfilename7, char *outfilename3, char *outfile_final, struct Lattice *reference_latt, struct Particle *reference7, struct Particle *reference3, long *idum_err, int nblatt_err);
void gene_b0_error_latt(struct Lattice *err_latt, long *idum, double rms_percentage_error);
void scan_error(char *outfilename7, char *outfilename3, char *outfile_final, struct Lattice *reference_latt, struct Particle *reference7, struct Particle *reference3, long *idum_err, int nblatt_err);
int adjustclosedorbit(struct Lattice *latt, struct Particle *part, double r0, double k, double *b00, double *b01);
void iteration_k2(char *outfilename, struct Lattice *latt, struct Particle *reference_part, double eps_clo, double k2step, double k2max, double xref);
int noangleatend(struct Lattice *latt, struct Particle *reference_part, int print_option);
void scan_alpha_zero(char *outfilename, struct Lattice *latt, struct Particle *reference_part);
void scan_alpha_zero2(char *outfilename, struct Lattice *latt, struct Particle *reference_part);
int adjust_b_maxangle(struct Lattice *latt, struct Particle *reference_part, double auth_angle, double precision, int print_option);
void gene_ellibunch_twiss_x(struct Beam *beam, double emitx, double betax, double alphax, int nbparts);
void gene_ellibunch_twiss_z(struct Beam *beam, double emitz, double betaz, double alphaz, int nbparts);
void gene_4D_beam_gasdev(struct Beam *beam, double emit_x, double twiss_betax, double twiss_alphax, double emit_z, double twiss_betaz, double twiss_alphaz, double factor_rms_tot, int nparts, long seed);
void compute_av_radius(char *trackout, char *outfilename);
double deriv_u_dt(int order, double rlce, double rlcs, double ee, double es);
double deriv_u_dr(int order, double r, double ee, double es, double rlcte, double rlcts);
void derudt(double rlce, double rlcs, double ai1, double bi1, double ci1, double *ai, double *bi, double *ci);
void derudr(int order, double rlcte, double rlcts, double ai1, double bi1, double ci1, double *ai, double *bi, double *ci);
void compute_field_spi_enge(struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double r0, double rmax, int nbr, double th0, double thmax, int nbth, double z0, double zmax, int nbz);
void compute_fieldmap(char *filename, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double r0, double rmax, int nbr, double th0, double thmax, int nbth, double z0, double zmax, int nbz);
void analysis_field_spi_enge(char *rfilename, char *outfilename, int order, int nb_elements);
double maxwell_test_max(double x, double y, double z, double dx, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), int doyouglob);
int maxwell_test(double x, double y, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int),double *divb, double *curlbx, double *curlby, double *curlbz, int doyouglob);
void maxwell_test_txtfile(char *textfile, double x, double y, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void maxwell_test_str_txtfile(char *textfile, double x, double y, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void maxwell_test_dxeffect(char *textfile, double x, double y, double z, double dxmin, double dxmax, int nbstep, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void maxwell_test_across_cell(char *textfile, double r, double z, double dx, int th_step, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
double maxwell_test_across_cell_max(double r, double z, double dx, int th_step, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void maxwell_test_map_radius(char *txt_prefix, double r0, double rmax, int r_step, double th0, double thmax, int th_step, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), int doyoucirc);
double convergence_radius(double dx, int th_step, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void compute_convergence_limit_vffa(struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));//, int doyouglob);
void read_fieldmap_boundaries(char * file, char *fileout);
void remove_useless_lines_eritmap(char * file, char *fileout);
void compare_field(char *file_mod, char *file_map, struct Cell *cell_mod, struct Cell *cell_map, void(*add_contribution_comp_mod)(double,double,double,double*,double*,double*,struct Cell*,int), void(*add_contribution_comp_map)(double,double,double,double*,double*,double*,struct Cell*,int));
void gnuplot_fit_fringe_field(char *txtfilename1, char *xcolumn1, char *ycolumn1, struct Cell *cell, double r, double a, char *with1, char *with2, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption);
void shift_theta_field_map(char *readname, char *writename, int nsteps_r, int nsteps_th, int nsteps_z);
void test_symmetry_map(struct Map *map, int iz);
void check_diff_map_coordinates(struct Map *map_tosca, struct Map *map_mod);
void test_maxwell_cylmap(struct Map *map, int ir, int iz, char *textfile);
void compute_map(char *filename, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double unitlength, double bscale, int loop_order, int line_order, int header_skip, double rstart, double rend, int nbr, double rshift, double thstart, double thend, int nbth, double thshift, double zstart, double zend, int nbz, double zshift);
void create_sym_map_latt(struct Lattice *latt_ori, struct Lattice *latt_sym);
void compare_map_mapfile(struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), char *mapfile, double bscale, double coord_scale, char *outfile);
void filter_fieldmap(char *fieldmap, char *output);
void da_study_fets(struct Lattice *latt, struct Particle *part_ref, char *fileout, double angle);
void scan_k2(char *outfile, struct Lattice *latt, struct Particle *part_ref, double k2_start, double k2_end, int nbstep_k2, double eps_clo, double percentage_mom, double amp);
void compare_diffield_added_separate(char *file_bx, char *file_by, char *file_bz, double r0, double rmax, int r_step, double th0, double thmax, int th_step, double z, struct Lattice *latt1, struct Lattice *latt2, void(*add_contribution_comp1)(double,double,double,double*,double*,double*,struct Cell*,int), 
	void(*add_contribution_comp2)(double,double,double,double*,double*,double*,struct Cell*,int), double *bxmin, double *bxmax, double *bymin, double *bymax, double *bzmin, double *bzmax);
double fn_co(struct Lattice *latt, struct Particle *part, double x[], double betax, double betaz);
void find_co_nelmin(struct Lattice *latt, struct Particle *part, double betax, double betaz, double fn_co(struct Lattice*, struct Particle*, double x[], double, double), int n, double xmin[], double *ynewlo, double reqmin, double step[], int konvge, int kcount, int *icount, int *numres, int *ifault);
int put_on_co_nelmin(struct Lattice *latt, struct Particle *part, double reqmin, double stepx, double stepux, double stepz, double stepuz, int doyouprintf);
void maxwell_test_vffa(char *file_input, struct Lattice *latt, double interpol_order, double r0, double rmax, int r_step, double th0, double thmax, int th_step, double z, double dx, double dy, double dz, int doyoucirc);
void brute_force_search(struct Lattice *latt, struct Particle *reference, double xmin, double xmax, int nbx, double zmin, double zmax, int nbz, double uxmin, double uxmax, int nbux, double uzmin, double uzmax, int nbuz);
void compare_field_point(char *file_dif_prefix, double x_1, double y_1, double z_1, struct Cell *cell1, void(*add_contribution_comp1)(double,double,double,double*,double*,double*,struct Cell*,int), double x_2, double y_2, double z_2, struct Cell *cell2, void(*add_contribution_comp2)(double,double,double,double*,double*,double*,struct Cell*,int), int doyouradial);
void compare_field_cell_heatmap(char *file_dif_prefix, struct Cell *cell1, void(*add_contribution_comp1)(double,double,double,double*,double*,double*,struct Cell*,int), struct Cell *cell2, void(*add_contribution_comp2)(double,double,double,double*,double*,double*,struct Cell*,int), 
  double r0, double rmax, int nbr, double th0, double thmax, int nbth, double z0, double zmax, int nbz, int doyouradial);
void find_boundaries_2d_file(char *file, double *xmin, double *xmax, double *ymin, double *ymax);
void find_boundaries_heated_map_file(char *file, double *xmin, double *xmax, double *ymin, double *ymax, double *val_min, double *val_max);
void find_boundaries_heated_map_file2(char *file, char *xrange, char *yrange, double *val_min, double *val_max);
void compute_fieldmap_special(struct Cell *cell2, char *filename, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int),
							double r0, double rmax, int nbr, double th0, double thmax, int nbth, double z0, double zmax, int nbz);
void file_vffa_btemp(char *trackout);
void reverse_cell(struct Cell *cell_ori, struct Cell *cell_chged);

void betafunc_period_vffa(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Particle *reference, struct Lattice *latt);
void adjust_particle(struct Particle *part_to_change, struct Particle *part_ori, double angle_rad);
void compute_map_polar(char *filename, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), 
	double rstart, double rend, int nbr, double rshift, double thstart, double thend, int nbth, double thshift, double zstart, double zend, int nbz, double zshift);
void compute_fieldheatmap(char *filename_prefix, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double r0, double rmax, int nbr, double th0, double thmax, int nbth, double z0, double zmax, int nbz, int doyouradial);
void compute_fieldheatmap_cart_hor_vert(char *filename_prefix, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double x0, double xmax, int nbx, double y0, double z0, double zmax, int nbz);
void scan_vffa(char *filename, double m_step, int nb_m, double b0d_step, int nb_b0d, double angle_step, int nb_angle, double r0d_step, int nb_r0d, struct Lattice *latt, struct Particle *part);
int acceptanceu(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double aumin, double aumax, double austep, double av);
int acceptancev(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double avmin, double avmax, double avstep, double au);
int acceptanceu_auto(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double aumin, double *aumax, double austep, double av);
int acceptancev_auto(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double avmin, double *avmax, double avstep, double au);
void vffa_sect_to_rect(char *filename, int n, double theta0_deg, double r0, double b0, double z0, double m, double rho, double theta_en_deg, double theta_ex_deg, double lambda, int interpol_order);
void set_vffa_rect_convergence_limit(struct Cell *cell, int doyoucompute, double conv_lim, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void enge_mult(double b0_mult, double efbe, double efbs, double lambda_add, char *xrange, char *outfile);
void enge_add(double b0_add, double efbe, double efbs, double lambda_add, char *xrange, char *outfile);
void enge_mult_add(double b0_mult, double b0_add, double efbe, double efbs, double lambda_add, char *xrange, char *outfile);
void adjust_b0_sect_from_rect(double x, double y, double z, double bz_opal, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double precision, int doyouf);
void test_add_ff(double efbe, double efbs, double lambda_tan, char *xrange);
void vffa_maxwell_test_heatmap(char *txt_prefix, char *txtaspect, double r0, double rmax, int r_step, double th0, double thmax, int th_step, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), int doyoucirc, double order);
void vffa_compare_scode(char *fileout, char *filescode, double xmin, double xmax, double xstep, double ymin, double ymax, double ystep, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void compute_fieldheatmap_scode(char *fileout, char *filescode, double xmin, double xmax, double xstep, double ymin, double ymax, double ystep, struct Cell *cell);
void wrapper_plot_heatmap_2dfile(char *file, char *file_2d, char *xcolumn, char *ycolumn, char *with, char *xlabel, char *ylabel, char *cblabel, char *psfilename, char *setoption);
void wrapper_plot_fieldheatmap(char *fileout, char *aspectfile);
void aspect_rect_vffa3d(char *filename, struct Cell *cell, int doyouerr);
void write_alierror_para_triplet_vffa(char *filename, struct Lattice *err_latt);
int error_study_triplet_vffa(char *prefix, struct Lattice *latt, struct Particle *ref_part, long seed, double rms_shift_error, double rms_twist_error);
void compare_err_vffa_triplet(char *noerr_file, char *err_file, char *dif_file, double trans_shift, double rot_shift);
void gene_matched_beam_emi_4d(char *txtfile, double emit_u, double twiss_betau, double twiss_alphau, double emit_v, double twiss_betav, double twiss_alphav, double r_inv[4][4], int nparts, long seed, int doyouhalo);
void test_matched_beam_emi(double m4d[4][4], char *textfile, char *textout);
void file_histo(char *textout, char *textfile, int columny, double xmin, double xmax, int nbstep);
void file_mean(char *textout, char *textfile);
void m_loc_int_bz_cell_vffa(char *outfilename, double x0, double ymin, double ymax, int nstep_y, double zmin, double zmax, int nstep_z, int nzref, double mref, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void m_loc_int_by_cell_vffa(char *outfilename, double x, double ymin, double ymax, int nstep_y, double zmin, double zmax, int nstep_z, int nzref, double mref, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void int_bz_liney_cell2(double *intbdl, double *intbdlabs, double x, double ymin, double ymax, int nstep, double z, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void int_by_liney_cell(double *intbdl, double *intbdlabs, double x, double ymin, double ymax, int nstep, double z, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void get_range_column(char *txtfile, int nbcolumn, int columnx, double *xmin, double *xmax);
void compute_cod_error_vffa_triplet(char *txtfile, double nsteps_cell, struct Lattice *latt, struct Particle *reference);
void copy_error_latt(struct Lattice *copy_latt, struct Lattice *ori_latt);
void da_alignment_shit(char *txtnoerr, char *txtfile, char *txtout, double trans_shift, double rot_shift);
void compute_dif_loc_m_vffa(char *txt_prefix, double x, double ymin, double ymax, int nstep_y, double zmin, double zmax, int nstep_z, double m_ref, struct Cell *cell, void (*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void fets_fieldmap_study(char *latt_file, char *txt_prefix);
void aspect_rect_vffa3d_fig_paper_f(char *filename, struct Cell *cell, int doyouerr);
void aspect_rect_vffa3d_fig_paper_d(char *filename, struct Cell *cell, int doyouerr);
void config_maker(char *name, double x0, double x1, double x2, double x3, double x4, double x5, double x6, double x7);
void compute_circmap_from_cartmap(char *textfile, double x0, double y0, double r0, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, struct Cell *cell);
void vffa_map_maker(char *textfile, double x0, double y0, double r0, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, char *lattfile);
void write_maplatt_file(char *lattfile, char *title, int periodicity, char *mapfile, double stepsize, double cell_length, double rmin, double rmax, int nbr, double thmin, double thmax, int nbth, double zmin, double zmax, int nbz);
void compute_vffa_error_m_int_bz(double *error_min, double *error_max, double x0, double ymin, double ymax, int nstep_y, double zmin, double zmax, int nstep_z, int nzref, double mref, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void compute_circmap_from_cartmap_superimpose(char *textfile, double x0, double y0, double r0, double xshift, double scale1, double scale2, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, struct Cell *cell);
void map_neighbours(double *bx, double *by, double *bz, double x0, double y0, double r0, double theta_tot, double r, double th, double z, double scale, struct Cell *cell);
void compute_vffacircmap_doublecoil_singlemagnet(char *cartmapfile, char *circmapfile, double x0, double y0, double r0, double xshift, double scale1, double scale2, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep);
void compute_cartmap_superimpose_shift(char *cartmapfile, char *cartmapfile_final, double xshift, double scale1, double scale2, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz);
void check_shinji_file(char *fileout, char *shinji_file, struct Cell *cell, double zshift, double yshift, double shinji_scale);
void compute_fieldheatmap_shinji(char *filename_prefix, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int),double r0, double rmax, int nbr, double th0, double thmax, int nbth, double z0, double zmax, int nbz, double angle_shinji_deg);
void check_shinji_file_circmap(char *fileout, char *shinji_file, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double zshift, double yshift, double shinji_scale, double r0, double th0_deg);
void compute_m_value_vffa_prototype(char *latt_name, char *out_long, char *out_vert, double x0, double ymin, double ymax_vert, int nby_vert, double ymax_long, int nby_long, double zmin, double zmax, int nstep_z, int nzref, double mref);
void compute_circmap_from_cartmap_superimpose2(char *textfile, double x0, double y0, double r0, double th0_deg, double xshift, double scale1, double scale2, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, struct Cell *cell);
void map_neighbours2(double *bx, double *by, double *bz, double x0, double y0, double r0, double th0, double theta_tot, double r, double th, double z, double scale, struct Cell *cell, int kill_notalive);
void map_neighbours3(double *bx, double *by, double *bz, double theta_tot, double r, double th, double z, double scale, struct Cell *cell, int kill_notalive);
void compute_vffacircmap_doublecoil_singlemagnet2(char *cartmapfile, char *circmapfile, double x0, double y0, double r0, double th0_deg, double xshift, double scale1, double scale2, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep);
void translate_shinji_vffa_file(char *shinji_file, char *outfile, double zshift, double yshift, double scale, double r0, double th0_deg, double theta_tot_deg);
void compute_circmap_from_cartmap2(char *textfile, double x0, double y0, double r0, double th0_deg, double scale, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, struct Cell *cell);
void vffa_map_maker2(char *textfile, double x0, double y0, double r0, double th0_deg, double scale, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, char *lattfile);
void translate_shinji_vffa_file_wrapper(char *txt_prefix, double zshift, double yshift, double scale, double r0, double th0_deg, double theta_tot_deg);
void convert_shinji_coord_to_ff(double *r, double *th, double *z, double xshinji, double yshinji, double zshinji, double r0, double th0, double yshift, double zshift);
void convert_ff_coord_to_shinji(double *xshinji, double *yshinji, double *zshinji, double r, double th, double zff, double r0, double th0, double yshift, double zshift);
void convert_ff_coord_to_shinji2(double *xshinji, double *yshinji, double *zshinji, double x, double y, double zff, double r0, double th0, double angle, double xshift, double yshift, double zshift);
void create_cartesian_mesh_tilted(char *textout, double yshift, double zshift, double r0, double th0_deg, double theta_tot_deg, double scale, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz, struct Cell *cell);
void create_cartesian_mesh_tilted_wrapper(char *lattfile, char *textout, double yshift, double zshift, double r0, double th0_deg, double theta_tot, double scale, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz);
void create_cartesian_mesh_tilted_superimpose(char *textout, double yshift, double zshift, double r0, double th0_deg, double theta_tot, double xshift, double scale_minus, double scale_plus, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz, struct Cell *cell);
void create_cartesian_mesh_tilted_superimpose_wrapper(char *lattfile, char *textout, double yshift, double zshift, double r0, double th0_deg, double theta_tot, double xshift, double scale_minus, double scale_plus, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz);
void compute_field_longit_line(char *textfile, double x0, double z0, double ymin, double ymax, int nby, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void function_textfile(char *textfile, double a, double b, double c, double m, double ymin, double ymax, int nby);
void m_loc_int_by_cell_vffa2(char *outfilename, double x, double ymin, double ymax, int nstep_y, double zmin, double zmax, int nstep_z, int nzref, double mref, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void compute_m_value_vffa_prototype2(char *latt_name, char *out_long, char *out_vert, double x0, double ymin, double ymax_vert, int nby_vert, double ymax_long, int nby_long, double zmin, double zmax, int nstep_z, int nzref, double mref);
void create_straight_map_fodo(char *map_out, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz, char *latt_d_name, double x0d, double y0d, double z0d, double scaled, char *latt_f_name, double x0f, double y0f, double z0f, double scalef);
void map_neighbours_straight(double *bx, double *by, double *bz, double x0, double y0, double z0, double ymax, double scale, double x, double y, double z, struct Cell *cell, int kill_notalive);
void write_beam_file(char *beamfile, char *title, double e0_mev, double ekin_ev, int nbkeywords, double x0, double z0, double xprime_deg, double zprime_deg);
void scattered_point_file_vffa_alignement(char *txtnoerr, char *txtfile, char *txtoutu, char *txtoutv, char *txtoutcod, char *txtouttuneu, char *txtouttunev, double trans_shift, double rot_shift, double scale, double ustep, double vstep, double codstep, double tunestep);
void wrapper_fodo_vffa_from_map(char *foldername, double xmin, double xmax, int nbx, double ymax, int nby, double zmin, double zmax, int nbz, char *lattname_ori, double x0d, double y0d, double z0d, double scaled, double x0f, double y0f, double z0f, double scalef);
void tune_for_tunescan(char *filename, struct Particle *part, struct Lattice *latt, double eps_clo, double amp);
void adjust_triplet_with_two_magnets(char *lattfile_triplet, char *lattfile_twomag);
void easyplot_3dtraj(char *trackout1, char *trackout2, char *trackout3, struct Lattice *latt, double xmin, double zmin, char *title1, char *title2, char *title3, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption);
void easyplot_3dtraj_5mom(char *trackout1, char *trackout2, char *trackout3, char *trackout4, char *trackout5, struct Lattice *latt, double xmin, double zmin, 
char *title1, char *title2, char *title3, char *title4, char *title5, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption);
void filter_fieldmap_Bf(char *text_fieldmap_Bf, char *outfile, double x0, double z0);
void coordinate_global_to_rect_vffa(double x_g, double y_g, double *x_mag, double *y_mag, double r_mag_cent, double th_mag_cent, double tilt_mag, double x_1, double y_1);
void compare_symetry_field_straight(char *textout, double xmin, double xmax, int nbx, double zmin, double zmax, int nbz, struct Cell *cell);
int calc_tune_eigenvalue(double *qu, double *qv, double ampx, double ampxp, double ampz, double ampzp, struct Particle *part, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, char *outfile);
void max_straight_vffa_fdratio_study(char *textout, double fdstart, double fdend, int nbfd, struct Particle *part, struct Lattice *latt);
void compute_conv_lim_long_vffa_rect_str(struct Cell *cell, char *textout, int nby, double xmax_start);
void compute_str_map(char *filename, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double unitlength, double bscale,
	int loop_order, int line_order, int header_skip, double xstart, double xend, int nbx, double ystart, double yend, int nby, double zstart, double zend, int nbz);
void maxwell_test_across_str_cell(char *textfile, double dx, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double z, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void compute_fodo_circmap(char *textfile, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, 
	char *latt_d_name, double r0d, double th0d_deg, double x0d, double y0d, double z0d, double scaled, char *latt_f_name, double r0f, double th0f_deg, double x0f, double y0f, double z0f, double scalef);
void wrapper_vffa_map_checker(char *pref_out, struct Lattice *latt, struct Particle *part, double emin_ev, double emax_ev, int nbstep);
void reduce_alan_map(char *mapout, char *lattout, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
void reduce_alan_map2(char *mapout, char *lattout, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
void write_maplatt_cyldeg_file(char *lattfile, char *title, int periodicity, char *mapfile, double stepsize, double cell_length, double rmin, double rmax, int nbr, double thmin, double thmax, int nbth, double zmin, double zmax, int nbz);
void wrapper_kieran_map(char *outfile, double xmax, double ymax);
void write_beam(char *beamfile, struct Particle *part, double m0, int q, double ekin_ev);
void write_beam2(char *beamfile, struct Particle *part);
void compute_fdf_circmap(char *textfile, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, char *latt_d_name, double r0d, double th0d_deg, double x0d, double y0d, double z0d, double scaled, char *latt_f_name, double r0f, double th0f_deg, double x0f, double y0f, double z0f, double scalef);
void compute_fodo_circmap2(char *textfile, double rmin, double rmax, int rstep, double theta_min_deg, double theta_max_deg, int thstep, double zmin, double zmax, int zstep, 
	char *latt_d_name_low, double x0d_low, double y0d_low, char *latt_d_name_high, double x0d_high, double y0d_high, double z0d, double r0d, double th0d_deg, double scaled, 
	char *latt_f_name_low, double x0f_low, double y0f_low, char *latt_f_name_high, double x0f_high, double y0f_high, double z0f, double r0f, double th0f_deg, double scalef);
void map_neighbours_low_high(double *bx, double *by, double *bz, double x0_high, double y0_high, struct Cell *cell_high, double x0_low, double y0_low, struct Cell *cell_low, double r0, double th0, double theta_tot, double r, double th, double z, double scale, int kill_notalive);
void wrapper_kieran_map2(char *outfile, double xmax, double ymax);
void create_circmap_triplet_from_cartmap(char *textfile, char *latt_f_name, double r0f, double th0f1_deg, double anglef1_deg, double scalef, char *latt_d_name, double r0d, double th0d_deg, double angled_deg, double scaled, double rmin, double rmax, int rstep, double theta_max_deg, int thstep, double zmin, double zmax, int zstep);
void compute_fieldheatmap_btot(char *filename, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int),
							double x0, double xmax, int nbx, double y0, double ymax, int nby, double z0, double zmax, int nbz);
double general_m_computation(double b_z, double z, double b_z0, double z0);
double local_m_computation(double b_zdz, double b_z, double dz);
void compute_m_rect_magnet(char *textout, char *text_lattmap, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz, double halfgap, double long_halflength, double dz);
//void compute_full_map_comsol_symmetries(char *text_in, char *text_out, double xmin, double xmax, int nb_x, double ymin, double ymax, int nb_y, double zmin, double zmax, int nb_z);
void compute_full_map_comsol_symmetries(char *latt_quarter, char *text_out, double xmin, double xmax, int nb_x, double ymin, double ymax, int nb_y, double zmin, double zmax, int nb_z);
void find_min_max_column_textfile(char *textfile, int nb_col_tot, int col_nb, double *min, double *max, int jump_header);
void vffa_rect_magnet_create_m_map(char *textout, char *text_lattmap, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz, double halfgap, double long_halflength, double dz);
void vffa_str_fodo_create_m_map(char *textout, char *text_lattmap, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double zmin, double zmax, int nbz, double halfgap, double long_halflength, double dz);
void compute_int_m_track_map(char *outfile, struct Cell *m_cell, struct Lattice *track_latt, struct Particle *part);
int adjust_brho_part_vffa(double z_target, double m_map, double z_step, struct Particle *reference, struct Lattice *latt, int doyouprintf);
void plot_with_without_iron(char *trackout, char *heatmapfile, char *tuneout, char *epsheatmap, char *eps_long_vert, char *eps_tune, double p0);
int adjust_b0_ffag_spi(struct Lattice *latt, struct Particle *part, double r_co, double eps_clo, double eps_r0);
int adjust_b0_ffag_spi_rav(struct Lattice *latt, struct Particle *part, double r_co_av, double eps_clo, double eps_r0);
void write_function_textfile(char *textfile, double xmin, double xmax, int nbx);
void acceptance_2d_emit_auto(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double emitmin, double *emitmax, double emitstep);
void acceptance_2d_turns(char *textout, struct Lattice *latt, struct Particle *part, int turns);
void read_change_raccam_fieldmap(char *file_opera, char *fileout);
void emittance_2d_expansion_order(char *textout, struct Lattice *latt, struct Particle *reference, double order);
void emittance_2d_stepsize(char *textout, struct Lattice *latt, struct Particle *reference, double stepsize);
int iterative_calc_tune_twiss(struct Particle *reference, double *qx, double *qz, double *betax, double *alphax, double *betaz, double *alphaz, double amp_x, double amp_z, 
	struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, char *outfile);
double compute_av_radius_from_trackout_1cell(char *trackout, int doyouprintf);
void create_textfile_pynaff(char *textfile, int nbpass, double amp_x, double amp_z, struct Particle *reference, struct Lattice *latt);
void write_get_bfield(char *textfile, double x, double y, double z, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), int doyoupolar);
void emittance_2d_brho(char *textout, struct Lattice *latt, struct Particle *reference, double brho);
void maxwell_heatmap(char *txt_prefix, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double z, double dx, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void maxwell_heatmap_wrapper(char *txt_prefix, double xmin, double xmax, int nbx, double ymin, double ymax, int nby, double z, double dx, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void write_in_txtfile(char *textfile, int nbpoints, double r0, double b0, double rmin1, double rmax1, double k1, double rmin2, double rmax2, double k2);
void write_beta_in_txtfile(char *textfile);
void check_tune_mag_pos(char *textfile, struct Particle *part, struct Lattice *latt, int nbpoints);
void compute_tm_collimator_for_emi(char *textfile, struct Lattice *latt, struct Particle *part, double angle_collimator_deg);
void gene_matched_beam_emi_4d_hffa(char *emplacement, double x_co, double xprime_co, double emit_x, double twiss_betax, double twiss_alphax, double emit_z, double twiss_betaz, double twiss_alphaz, int nparts, long seed, int doyouhalo);
void fets_ffa_corners_tambouille(char *lattfile, char *beamfile, char *trackout, char *field_file);
void compute_min_max_radius_from_trackout_1cell(char *trackout, double *rmin, double *th_rmin, double *rmax, double *th_rmax, double thmin_deg, double thmax_deg, int doyouprintf);
void check_rav_fets_ffa(char *lattfile, char *injfile, char *extfile);
void round_number_fets_ffa_field_file(char *textout, char * textin);
void plot_latt_fets_ffa(char *latt_name, char *beamname, char *trackout);
void print_get_bfield_across_cell_radius(char *fileout, double r, int nbsteps, double z, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void comp_betafunc_lattices_fets_ffa(char *latt_name, char *beamname, char *betafile);
int acceptance_uv(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double au, double av, double decoup_to_coup[4][4]);
void change_para_triplet_cell_rect_vffa(struct Cell *cell, double cell_length, double theta_cell_deg, double magnetf_length, double small_drift, double magnetd_length, double shiftf, double shiftd, double tilt_magf_deg, double b0f, double b0d, double m);
void check_fit_iker(char *textfile);
void check_fit_iker2(char *textfile);
void check_fit_iker3(char *textfile);
void check_fit_iker4(char *textfile);
void check_fit_iker5(char *textfile, double r0, double opening_angle_deg, double lambda);
void write_field_const_radius_across_cell(char *textfile, double r, double z, int nbsteps, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void shift_bfield_tosca(char *readfile, char *writefile, double *bdmax, double *thdmax, double *bfmax);
void fit_gnuplot_kurns_mr(char *textfile, char *fiteps, double bd0, double cd1, double cd2, double thcd, double ffbd, double bf0, double cf, double ffbf);
void wrapper_fit_field_toscamap_kurri_mr(char *prefix, double r, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int));
void check_fit_iker_doublet2(char *textfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda);
void check_fit_iker_doublet1(char *textfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda);
void check_fit_iker_doublet4(char *textfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda);
void fit_iker_clean_file3(char *textfile);
void fit_iker_doublet2_clean(char *cleanfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda);
void fit_iker_doublet1_clean(char *cleanfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda);
void fit_iker_doublet4_clean(char *cleanfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda);
void fit_b_iker_doublet2_clean(char *cleanfile, double r0, double opening_angle_deg_f, double opening_angle_deg_d, double lambda, double c0[2], double c1[2], double c2[2], double c3[2]);
int set_closed_orbit_xxp(struct Particle *part, double eps_clo, struct Lattice *latt, int i, int doyouprintf);
void launch_part_for_co(struct Particle *part, double *xbest, double *uxbest, double xmin, double xmax, int nbx, double uxmin_deg, double uxmax_deg, int nbux, struct Lattice *latt);
void codall_shinji_arrange(char *textfile, char *textout);
void fft_kurns(char *textfile, char *text_fft, double time_start, double amp_fft[], int nbpoints, double sample_step_sec);
double compute_max_radius_from_trackout_1cell(char *trackout, int doyouprintf);
int adjust_b0_ffag_spi_rmax(struct Lattice *latt, struct Particle *part, double r_co_max, double eps_clo, double eps_r0);

/* ****************** common MACROs ********************* */
#define MAX(a,b)    (( a > b ) ? a : b )
#define MIN(a,b)    (( a < b ) ? a : b )
#define SIGN(a,b)	((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SWAP(a,b)	tempr=(a);(a)=(b);(b)=tempr 

#endif
