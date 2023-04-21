/*
 *  toolbox.c
 *  ringdesign
 *
 *  Copyright 2009 Kyoto University. All rights reserved.
 *
 */

#include "toolbox.h"

extern void errorstop(char *error_text)
{
	char shell_write[100];
	FILE *shell_read;
	if(print_color==YES) COLOR("1;31");
	printf("\nerrorstop procedure...\n");
	printf("%s\n", error_text);
	printf("...now exiting to system...\n");
	if(print_color==YES) COLOR("0");
	if(sound_option==YES) {
		sprintf(shell_write,"say unexpected error");
		shell_read = popen(shell_write,"r");
		pclose(shell_read);
	}
	exit(1);
}

//when read a file, read until the end of the line
extern void newline(FILE *rfile)
{
	char temp_str[MAX_CHARINLINE+2];
	
	fgets(temp_str, MAX_CHARINLINE+2, rfile);
	if(strlen(temp_str) > MAX_CHARINLINE) errorstop("!!!ERROR in newline (toolbox.c): maximum number of char per line in input file exceeded!!!");
}

//create an empty data file
extern void emptyfile(char *filename)
{
	FILE *fill;
	
	fill = fopen(filename, "w");
	fclose(fill);
}

//print lattice parameters
extern void print_latt_para(char *filename, char *lattname, struct Lattice *latt)
{
	char cellname[KEYWORDL];
	int i;
	
	if(filename!=NULL) {
		FILE *wfile;
		
		printf("writing lattice in %s\n",filename);
		wfile = fopen(filename, "w");
		fprintf(wfile, "%s\n\n",lattname);
		fprintf(wfile, "%i				//lattice periodicity [int]\n", latt->periodicity);
		fprintf(wfile, "%i				//number of elements [int]\n", latt->nbcell);
		fprintf(wfile, "%i				//number of cell types [int]\n\n\n", latt->nbcell);
		fclose(wfile);
		for(i=0; i<latt->nbcell; i++) {
			sprintf(cellname, "cell[%i]", i);
			print_cell_para(filename, cellname, &(latt->cell[i]));
		}
	}
	else {
		printf("%s: \n", lattname);
		printf("nb of cells: %i\n", latt->nbcell);
		printf("superperiod: %i\n", latt->periodicity);
		for(i=0; i<latt->nbcell; i++) {
			sprintf(cellname, "cell[%i]", i);
			print_cell_para(filename, cellname, &(latt->cell[i]));
		}
		printf("\n");
	}
}

//print cell parameters does not work with ffag-tiltpol
extern void print_cell_para(char *filename, char *cellname, struct Cell *cell)
{
	int j, n;
	
	if(filename!=NULL) {
		FILE *wfile;
		wfile = fopen(filename, "a");
		fprintf(wfile, "1 %s %i	//%s\n", cell->keyword, cell->nbcomp, cellname);
		fprintf(wfile, "%lf			//step size [m]\n", cell->stepsize);
		if(test_cell_vffa(cell)==TRUE) fprintf(wfile, "%lf  %lf  %lf  %lf	//collimators rmin [m], rmax [m], zmin [m], zmax [m]\n", cell->collim.rmin, cell->collim.rmax, cell->collim.zmin, cell->collim.zmax);
		else fprintf(wfile, "%lf  %lf  %lf	//collimators rmin [m], rmax [m], zmax [m]\n", cell->collim.rmin, cell->collim.rmax, cell->collim.zmax);
		if(cell->boun.thmax != 0) fprintf(wfile, "%lf			//element total opening angle [deg]\n", cell->boun.thmax*180./PI);
		else if(cell->boun.ymax != 0) fprintf(wfile, "%lf			//element total length [m]\n", cell->boun.ymax);
		fprintf(wfile, "%lf				//deltar\n\n", cell->deltar);
		if(test_cell_map(cell)==NO) {
			//printf("not a map!\n");
			for(n = 0; n < cell->nbcomp; n++) {
				if(cell->boun.thmax != 0) fprintf(wfile, "%lf	", cell->mpara[n][0]*180./PI);
				else if(cell->boun.ymax != 0) fprintf(wfile, "%lf	", cell->mpara[n][0]);
				for(j = 1; j < SIZE_MPARA; j++) {
					if(	(test_cell_spiangle(cell)==6 && j==5) ||
						(test_cell_spiangle(cell)==1 && j==4) ||
						(test_cell_spiangle(cell)==3 && j==4) ||
						(test_cell_spiangle(cell)==4 && j==4)) {
						fprintf(wfile, "%lf	", cell->mpara[n][j]*180./PI);
					}
					else fprintf(wfile, "%lf	", cell->mpara[n][j]);
				}
				fprintf(wfile,"\n\t");
				if(cell->boun.thmax != 0 && test_cell_rect_vffa_bend_cell(cell)!=TRUE) fprintf(wfile,"%lf	", cell->efben[n][0]*180./PI);
				else if(cell->boun.ymax != 0 || test_cell_rect_vffa_bend_cell(cell)==TRUE) fprintf(wfile,"%lf	", cell->efben[n][0]);
				for(j = 1; j < SIZE_EFB; j++) fprintf(wfile,"%lf	", cell->efben[n][j]);
				fprintf(wfile,"\n\t");
				if(cell->boun.thmax != 0 && test_cell_rect_vffa_bend_cell(cell)!=TRUE) fprintf(wfile,"%lf	", cell->efbex[n][0]*180./PI);
				else if(cell->boun.ymax != 0 || test_cell_rect_vffa_bend_cell(cell)==TRUE) fprintf(wfile,"%lf	", cell->efbex[n][0]);
				for(j = 1; j < SIZE_EFB; j++) fprintf(wfile,"%lf	", cell->efbex[n][j]);
				fprintf(wfile,"\n\n");
			}
		}
		else if(test_cell_map(cell)==YES) {
			fprintf(wfile, "map symmetry: %i, (%i=YES, %i=NO)\n", cell->map.sym, YES, NO);
			fprintf(wfile, "map dimensions min:(%le, %le, %le), max: (%le, %le, %le)\n", cell->map.mapdim[0], cell->map.mapdim[2], cell->map.mapdim[4], cell->map.mapdim[1], cell->map.mapdim[3], cell->map.mapdim[5]);
			fprintf(wfile, "map step size (%le, %le, %le)\n", cell->map.stepsize[0], cell->map.stepsize[1], cell->map.stepsize[2]);
			fprintf(wfile, "nb of nodes: (%i, %i, %i)\n", cell->map.nnodes[0], cell->map.nnodes[1], cell->map.nnodes[2]);
								}
		else errorstop("strange, map or not?\n");
		fclose(wfile);
	}
	else {
		printf("%s: \n", cellname);
		printf("keyword: %s, nb of components: %i\n", cell->keyword, cell->nbcomp);
		printf("stepsize = %le\n", cell->stepsize);
		printf("framework: xc = %lf[m] yc = %lf[m] ae = %lf[deg]\n", cell->framework.xc, cell->framework.yc, cell->framework.ae*180./(PI));
		printf("deltar = %lf\n", cell->deltar);
		printf("boundary (thmax, ymax): (%le[deg], %le[m])\n", cell->boun.thmax*180./PI, cell->boun.ymax);
		if(cell->doyou_err==NO) printf("alignment error: NO\n");
		else if(cell->doyou_err==YES) printf("alignment error: YES\n");
		else errorstop("alignment error flag has a strange value!\n");
		if(test_cell_map(cell)==NO) {
			for(n = 0; n < cell->nbcomp; n++) {
				printf("\t");
				for(j = 0; j < SIZE_MPARA; j++) printf("%lf ", cell->mpara[n][j]);
				printf("\n\t\t");
				for(j = 0; j < SIZE_EFB; j++) printf("%lf ", cell->efben[n][j]);
				printf("\n\t\t");
				for(j = 0; j < SIZE_EFB; j++) printf("%lf ", cell->efbex[n][j]);
				printf("\n\t");
				if(cell->doyou_err==YES) {
					printf("ali_err: ");
					for(j = 1; j < 4; j++) printf("%lf ", cell->alierror[n][j]);
					printf("\n\t");//%lf\n\t\t", cell->alierror[n][4]);
					for(j = 4; j < 11; j++) printf("%lf ", cell->alierror[n][j]);
					printf("\n");
				}
				printf("\n");
			}
		}
		else if(test_cell_map(cell)==YES) {
			printf("map symmetry: %i, (%i=YES, %i=NO)\n", cell->map.sym, YES, NO);
			printf("map dimensions min:(%le, %le, %le), max: (%le, %le, %le)\n", cell->map.mapdim[0], cell->map.mapdim[2], cell->map.mapdim[4], cell->map.mapdim[1], cell->map.mapdim[3], cell->map.mapdim[5]);
			printf("map step size (%le, %le, %le)\n", cell->map.stepsize[0], cell->map.stepsize[1], cell->map.stepsize[2]);
			printf("nb of nodes: (%i, %i, %i)\n", cell->map.nnodes[0], cell->map.nnodes[1], cell->map.nnodes[2]);
		}
		else errorstop("strange, map or not?\n");
	}
}

//print particle parameters
extern void print_beam_para(char *beamname, struct Beam *beam)
{
	char partname[15];
	int i;
	
	printf("%s: \n", beamname);
	for(i=0;i<beam->npart;i++) {
		sprintf(partname,"part[%i]", i);
		print_part_para(partname, &(beam->part[i]));
	}
}

extern void print_part_para(char *partname, struct Particle *part)
{
	printf("%s: \n", partname);
	printf("framework: xc = %lf[m] yc = %lf[m] ae = %lf[deg]\n", part->fwk.xc, part->fwk.yc, part->fwk.ae*180./(PI));
	printf("(status,m0,q) = (%i, %le, %le)\n", part->status, part->m0, part->q);
	printf("(brho,s,t) = (%le, %le, %le)\n", part->brho, part->s, part->t);
	printf("(x,y,z) = (%.8e, %.8e, %.8e)\n", part->x, part->y, part->z);
	printf("(ux,uy,uz) = (%.8e, %.8e, %.8e)\n", part->ux, part->uy, part->uz);
	printf("x' = %.8e[deg], z' = %.8e[deg]\n\n", atan(part->ux/(part->uy))*180./(PI), atan(part->uz/(sqrt((part->uy)*(part->uy)+(part->ux)*(part->ux))))*180./(PI));
}

extern void print_get_bfield(double x, double y, double z, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	double bx, by, bz;
	double br,bth, r, th;
	
	get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp);
	//printf("at (%lf, %lf, %lf), field (%le, %le, %le)\n",x,y,z,bx,by,bz);
	if(cell->boun.thmax != 0) {
		r= sqrt(x*x + y*y); 
		th = atan_ratio(y, x);
		br = bx*cos(th) + by*sin(th);
		bth = by*cos(th) - bx*sin(th);
		printf("r=%lf [m], th=%lf [deg], z=%lf [m], (Br,Bth, Bz)=(%le, %le, %le)\n",r,th*180./PI, z, br, bth, bz);
	}
}

extern void printf_mat_1dim_complex(char *message, double complex *m, int nb_dim, char *name_mat, int option, char *txtout)
{
	int i,j;
	FILE *wfile;
	
	if(txtout !=NULL) {
		wfile = fopen(txtout, "a");
		fprintf(wfile, "%s\n", message);
	}
	
	printf("%s\n", message);
	for(i=0;i<nb_dim;i++) {
		if(option==1) {
			for(j=0;j<nb_dim;j++) {
				printf("%s[%i]=(%le, %le)\t",name_mat, i*nb_dim+j, creal(m[i*nb_dim+j]), cimag(m[i*nb_dim+j]));
				if(txtout != NULL) fprintf(wfile, "%s[%i]=(%le, %le)\t",name_mat, i*nb_dim+j, creal(m[i*nb_dim+j]), cimag(m[i*nb_dim+j]));
			}
			printf("\\\\\n");
			if(txtout != NULL) fprintf(wfile, "\\\\\n");
		}
		else if(option==2) {
			for(j=0;j<nb_dim;j++) {
				printf("%lf & ",creal(m[i*nb_dim+j]));
				if(txtout != NULL) fprintf(wfile, "%lf & ",creal(m[i*nb_dim+j]));
			}
			printf("\\\\\n");
			if(txtout != NULL) fprintf(wfile, "\\\\\n");
		}
		else if(option==3) {
			for(j=0;j<nb_dim;j++) {
				printf("%lf\t ",creal(m[i*nb_dim+j]));
				if(txtout != NULL) fprintf(wfile, "%lf\t ",creal(m[i*nb_dim+j]));
			}
			printf("\n");
			if(txtout != NULL) fprintf(wfile, "\n");
		}
	}
	//printf("\n");
	if(txtout !=NULL) {
		fprintf(wfile, "\n");
		fclose(wfile);
	}

}

extern void printf_mat_4d_real(char *message, double m[4][4], char *name_mat, char *txtout)
{
	int i,j,n=4;
	FILE *wfile;
	
	if(txtout !=NULL) {
		wfile = fopen(txtout, "a");
		fprintf(wfile, "%s\n", message);
		fprintf(wfile, "%s\n", name_mat);
	}
	
	printf("%s\n", message);
	printf("%s\n", name_mat);
	for(i=0;i<n;i++) {
		for(j=0;j<n;j++) {
			printf("%le 	",m[i][j]);
			if(txtout != NULL) fprintf(wfile, "%le  ",m[i][j]);
		}
		printf("\n");
		if(txtout != NULL) fprintf(wfile, "\n");
	}
	//printf("\n");
	if(txtout !=NULL) {
		fprintf(wfile, "\n");
		fclose(wfile);
	}
}

extern double round_nb(double value, int nb_decimal)
{
	double coef = pow(10., nb_decimal);
	return (round(value*coef))/coef;
}

// test if number is 2^n
extern int test_power2(unsigned int x)
{
 while (((x % 2) == 0) && x > 1) x /= 2; // While x is even and > 1
 if(x == 1) return TRUE;
 else return FALSE;
}

//get total energy, beta, gamma of a given Particle
extern void get_ebg_part(double *e_tot, double *beta, double *gamma, struct Particle *part)
{
	*e_tot	= sqrt(pow(part->m0*CLIGHT2, 2) + (part->brho*part->q*part->brho*part->q)*CLIGHT2);
	*gamma	= *e_tot/(part->m0*CLIGHT2);
	*beta	= sqrt((*gamma)*(*gamma) - 1)/(*gamma);
}

//generate the Twiss transfer matrix from a usual 2x2 transfer matrix
extern void twissTMtransform(double m[2][2], double twm[3][3]) {
	twm[0][0] = m[0][0]*m[0][0];
	twm[0][1] = -2.*m[0][0]*m[0][1];
	twm[0][2] = m[0][1]*m[0][1];
	twm[1][0] = -m[0][0]*m[1][0];
	twm[1][1] = m[0][0]*m[1][1] + m[0][1]*m[1][0];
	twm[1][2] = -m[0][1]*m[1][1];
	twm[2][0] = m[1][0]*m[1][0];
	twm[2][1] = -2.*m[1][0]*m[1][1];
	twm[2][2] = m[1][1]*m[1][1];
}

//warning: doyouprint must be YES otherwise doesn't work...
extern void fix_boundingbox(char *txtfilename, int doyouprint) {
	char shell_write[500], buf[MAX_CHARINLINE];
	FILE *shell_read;

	sprintf(shell_write,"./fixbb.sh %s", txtfilename);
	shell_read = popen(shell_write,"r");
	if(doyouprint==YES) {
		if(fgets(buf, MAX_CHARINLINE-1, shell_read) == NULL) printf("no output in fix_boundingbox");
		printf("%s",buf);
	}
	pclose(shell_read);
}

extern void say(char *text_say) {
	char shell_write[500];
	FILE *shell_read;
	sprintf(shell_write,"say %s", text_say);
	shell_read = popen(shell_write,"r");
	pclose(shell_read);
}

extern void convert_file(char *ori_file, char *final_file, int resolution, int doyouprint)
{
	char shell_write[500], buf[MAX_CHARINLINE];
	FILE *shell_read;
	
	if(resolution != 0) sprintf(shell_write,"convert -density %i %s %s", resolution, ori_file, final_file);
	else sprintf(shell_write,"convert %s %s", ori_file, final_file);
	shell_read = popen(shell_write,"r");
	if(doyouprint==YES) {
		if(fgets(buf, MAX_CHARINLINE-1, shell_read) == NULL) printf("no output in convert_file");
		printf("%s",buf);
	}
	pclose(shell_read);
	
}

extern double comp_step(double val_start, double val_end, int nvalues)
{
	double step;
	
	if(nvalues>1) {
		step=(val_end-val_start)/(nvalues-1);
		if(debug == YES) printf("start=%lf, end=%lf, n=%i, step=%le)\n",val_start, val_end, nvalues, step);
		return step;
	}
	else {
		if(debug == YES) printf("n=1, step=0\n");
		return 0;
	}
}

extern void write_ellipse(char *wfilename, double emitt_pimmmrad, double beta, double alpha, double xcenter, double ycenter) {
	int npoints = 100, i;
	double x, y, xn, yn, t;//, gamma_cs, angle;
	FILE *wfile;
	if(beta==0) errorstop("beta=0, cannot write ellipse");
	//gamma_cs = (1+alpha*alpha)/beta;
	//angle = 0.5*atan(alpha/(gamma_cs-beta));
	wfile = fopen(wfilename, "w");
	
	for(i = 0; i <= npoints; i++) {
		t = i*2.*PI/npoints;
		//xn = sqrt(emitt_pimmmrad*1.e-6)*cos(t)*sqrt(beta);
		//yn = sqrt(emitt_pimmmrad*1.e-6)*sin(t)/sqrt(beta);
		//x = xn*cos(angle)-yn*sin(angle);
		//y = yn*cos(angle) + xn*sin(angle);
		xn = sqrt(emitt_pimmmrad*1.e-6)*cos(t);
		yn = sqrt(emitt_pimmmrad*1.e-6)*sin(t);
		x = sqrt(beta)*xn;
		y = (yn - alpha*xn)/sqrt(beta);
		fprintf(wfile, "%le  %le\n", x + xcenter, y + ycenter);
	}
	
	fprintf(wfile, "\n%le %le\n", xcenter, ycenter);
	
	fclose(wfile);
	
}

extern void write_ellipse2(char *wfilename, double emitt_pimmmrad, double beta, double alpha, double xcenter, double ycenter) {
	int npoints = 100, i;
	double x, y, xmin, xstep, xmax, gamma_cs;
	FILE *wfile;
	if(beta==0) errorstop("beta=0, cannot write ellipse");
	gamma_cs = (1+alpha*alpha)/beta;
	//angle = 0.5*asin(alpha/(gamma_cs-beta));
	wfile = fopen(wfilename, "w");
	xmin = -sqrt(emitt_pimmmrad*1.e-6/gamma_cs);
	xmax = sqrt(emitt_pimmmrad*1.e-6/gamma_cs);
	xstep = (xmax-xmin)/npoints;
	
	for(i = 0; i <= npoints; i++) {
		x = xmin+i*xstep;
		y = (-alpha*x+sqrt(beta*emitt_pimmmrad*1.e-6-x*x))/beta;
		fprintf(wfile, "%le  %le\n", x + xcenter, y + ycenter);
	}
	for(i=0;i<=npoints;i++) {
		x = xmax-i*xstep;
		y = (-alpha*x-sqrt(beta*emitt_pimmmrad*1.e-6-x*x))/beta;
		fprintf(wfile, "%le  %le\n", x + xcenter, y + ycenter);
	}
	
	fprintf(wfile, "\n%le %le\n", xcenter, ycenter);
	
	fclose(wfile);
	
}

extern void write_ellipse_mm_mevc(char *wfilename, double emitt_pimmmrad, double beta, double alpha, double xcenter, double ycenter, double pcenter) {
	int npoints = 1000, i;
	double x, y, xmin, xstep, xmax, gamma_cs;
	FILE *wfile;
	if(beta==0) errorstop("beta=0, cannot write ellipse");
	gamma_cs = (1+alpha*alpha)/beta;
	//angle = 0.5*asin(alpha/(gamma_cs-beta));
	wfile = fopen(wfilename, "w");
	xmin = -sqrt(emitt_pimmmrad*1.e-6/gamma_cs);
	xmax = sqrt(emitt_pimmmrad*1.e-6/gamma_cs);
	xstep = (xmax-xmin)/npoints;
	
	for(i = 0; i <= npoints; i++) {
		x = xmin+i*xstep;
		if(beta*emitt_pimmmrad*1.e-6-x*x<0.) {
			printf("sqrt neg (%le), force to 0\n", beta*emitt_pimmmrad*1.e-6-x*x);
			y = (-alpha*x)/beta*pcenter;
		}
		else y = (-alpha*x+sqrt(beta*emitt_pimmmrad*1.e-6-x*x))/beta*pcenter;
		fprintf(wfile, "%le  %le\n", x*1000. + xcenter, y + ycenter);
		//printf("%le, %le, %le\n",x, y, sqrt(beta*emitt_pimmmrad*1.e-6-x*x));
	}
	for(i=0;i<=npoints;i++) {
		x = xmax-i*xstep;
		if(beta*emitt_pimmmrad*1.e-6-x*x<0.) {
			printf("sqrt neg (%le), force to 0\n", beta*emitt_pimmmrad*1.e-6-x*x);
			y = (-alpha*x)/beta*pcenter;
		}
		y = (-alpha*x-sqrt(beta*emitt_pimmmrad*1.e-6-x*x))/beta*pcenter;
		fprintf(wfile, "%le  %le\n", x*1000. + xcenter, y + ycenter);
	}
	x = xmin;
	if(beta*emitt_pimmmrad*1.e-6-x*x<0.) {
		printf("sqrt neg (%le), force to 0\n", beta*emitt_pimmmrad*1.e-6-x*x);
		y = (-alpha*x)/beta*pcenter;
	}
	else y = (-alpha*x+sqrt(beta*emitt_pimmmrad*1.e-6-x*x))/beta*pcenter;
	fprintf(wfile, "%le  %le\n", x*1000. + xcenter, y + ycenter);
	
	fprintf(wfile, "\n%le %le\n", xcenter, ycenter);
	
	fclose(wfile);
}

extern int get_nb_lines_file(char *txtfilename) {
	char name[500], shell_write[500];
	double nblines_double;
	FILE *shell_read;

	sprintf(shell_write,"wc -l %s", txtfilename);
	shell_read = popen(shell_write,"r");
	fscanf(shell_read, "%lf %s", &nblines_double, name);
	pclose(shell_read);
	return (int) nblines_double;
}

extern void get_smax_file(char *txtfilename, double *smax) {
	int n, nblines;
	FILE *rfile = NULL;
	
	nblines = get_nb_lines_file(txtfilename);
	rfile = fopen(txtfilename, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(n=0;n<nblines-2;n++) newline(rfile);	
	fscanf(rfile, "%le", smax);
	fclose(rfile);
}

extern void get_betamax_file(char *txtfilename, double *betamax) {
	int n, nblines;
	double s, bx, ax, bz, az, bmax=0.;
	FILE *rfile = NULL;
	
	nblines = get_nb_lines_file(txtfilename);
	rfile = fopen(txtfilename, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(n=0;n<nblines-1;n++) {
		fscanf(rfile, "%le	%le	%le	%le	%le", &s, &bx, &ax, &bz, &az);
		if(bmax < bx || bmax < bz) bmax = MAX(bx, bz);
	}
	*betamax = bmax;
	//printf("betamax = %lf\n", bmax);
	fclose(rfile);
}

extern void get_betax_z_max_file(char *txtfilename, double *betax_max, double *betaz_max) {
	int n, nblines;
	double s, bx, ax, bz, az, bxmax=0.,bzmax=0.;
	FILE *rfile = NULL;
	
	nblines = get_nb_lines_file(txtfilename);
	rfile = fopen(txtfilename, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(n=0;n<nblines-1;n++) {
		fscanf(rfile, "%le	%le	%le	%le	%le", &s, &bx, &ax, &bz, &az);
		if(bxmax < bx) bxmax = bx;
		if(bzmax < bz) bzmax = bz;
		
	}
	*betax_max = bxmax;
	*betaz_max = bzmax;
	//printf("betaxmax = %lf, betazmax=%lf\n", bxmax,bzmax);
	fclose(rfile);
}

extern void get_xrange(char *txtfilename, char *xrange) {
	double smax;
	get_smax_file(txtfilename, &smax);
	//printf("smax = %lf\n", smax);
	sprintf(xrange,"[0:%lf]", smax);
}

extern void get_betarange(char *txtfilename, char *betarange) {
	double bmin, bmax, brange;
	get_betamax_file(txtfilename, &bmax);
	brange = floor(bmax)+1;
	//sprintf(betarange,"[0:%lf]", brange);
	
	get_betamin_file(txtfilename, &bmin);
	bmin = floor(bmin);
	sprintf(betarange,"[%lf:%lf]", bmin, brange);
}

extern void get_dispmax_file(char *txtfilename, double *dispmax) {
	int n, nblines;
	double s, d=-1.0e18, dp;
	FILE *rfile = NULL;
	
	nblines = get_nb_lines_file(txtfilename);
	rfile = fopen(txtfilename, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(n=0;n<nblines-1;n++) {
		fscanf(rfile, "%le	%le	%le", &s, &d, &dp);
		//printf("%le\n", d);
		if(*dispmax < d) *(dispmax) = d;
	}
	//printf("in get_dispmax_file, dmax=%lf", *dispmax);
	fclose(rfile);
}

extern void get_dispmin_file(char *txtfilename, double *dispmin) {
	int n, nblines;
	double s, d=1.0e18, dp;
	FILE *rfile = NULL;
	
	nblines = get_nb_lines_file(txtfilename);
	rfile = fopen(txtfilename, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(n=0;n<nblines-2;n++) {
		fscanf(rfile, "%le	%le	%le", &s, &d, &dp);
		//printf("%le\n", d);
		if(*dispmin > d) *(dispmin) = d;
	}
	//printf("in get_dispmax_file, dmax=%lf", *dispmax);
	fclose(rfile);
}

extern void get_disprange(char *txtfilename, char *disprange) {
	double dmax, dmin, drange;
	get_dispmax_file(txtfilename, &dmax);
	get_dispmin_file(txtfilename, &dmin);
	printf("dmin = %lf, dmax = %lf\n",dmin, dmax);
	drange = dmax*1.02;
	sprintf(disprange,"[%lf:%lf]", dmin, drange);
}

extern void get_dispmax_file_4d(char *txtfilename, double *dispxmax, double *dispzmax) {
	int n, nblines;
	double s, dx=-1.0e18, dpx, dz=-1.0e18, dpz;
	FILE *rfile = NULL;
	
	nblines = get_nb_lines_file(txtfilename);
	rfile = fopen(txtfilename, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(n=0;n<nblines-1;n++) {
		fscanf(rfile, "%le	%le	%le	%le	%le", &s, &dx, &dpx, &dz, &dpz);
		//printf("dx=%le, dz=%le\n", dx, dz);
		if(*dispxmax < dx) *(dispxmax) = dx;
		if(*dispzmax < dz) *(dispzmax) = dz;
	}
	//printf("in get_dispmax_file, dmax=%lf", *dispmax);
	fclose(rfile);
}

extern void get_dispmin_file_4d(char *txtfilename, double *dispxmin, double *dispzmin) {
	int n, nblines;
	double s, dx=1.0e18, dpx, dz=1.0e18, dpz;
	FILE *rfile = NULL;
	
	nblines = get_nb_lines_file(txtfilename);
	rfile = fopen(txtfilename, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(n=0;n<nblines-2;n++) {
		fscanf(rfile, "%le	%le	%le	%le	%le", &s, &dx, &dpx, &dz, &dpz);
		//printf("%le\n", d);
		if(*dispxmin > dx) *(dispxmin) = dx;
		if(*dispzmin > dz) *(dispzmin) = dz;
	}
	//printf("in get_dispmax_file, dmax=%lf", *dispmax);
	fclose(rfile);
}

extern void get_disprange_4d(char *txtfilename, char *disprangehor, char *disprangevert) {
	double dxmax, dxmin, dzmax, dzmin, dxrange, dzrange;
	get_dispmax_file_4d(txtfilename, &dxmax, &dzmax);
	get_dispmin_file_4d(txtfilename, &dxmin, &dzmin);
	printf("dxmin = %lf, dxmax = %lf, dzmin = %lf, dzmax = %lf\n",dxmin, dxmax, dzmin, dzmax);
	dxrange = dxmax*1.02;
	dzrange = dzmax*1.02;
	sprintf(disprangehor,"[%lf:%lf]", dxmin, dxrange);
	sprintf(disprangevert,"[%lf:%lf]", dzmin, dzrange);
}

extern void get_tune_range(char *txtfilename, double *qx_min, double *qx_max, double *qz_min, double *qz_max, double *qx_tot_min, double *qx_tot_max, double *qz_tot_min, double *qz_tot_max) {
	char buf[MAX_CHARINLINE];
	double e_tot, qx=0.0, qz=0.0, qx_tot=0.0, qz_tot=0.0, temp;
	FILE *rfile = NULL;
	
	rfile = fopen(txtfilename, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	if(fgets(buf, MAX_CHARINLINE-1, rfile) == NULL) errorstop("strange in get_tune_range");
	if(sscanf(buf, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le", &e_tot, qx_min, qz_min, qx_tot_min, qz_tot_min, &temp, &temp, &temp, &temp, &temp, &temp, &temp, &temp) != 13) errorstop("strange2 in get_tune_range");
	printf("qx=%lf, qz=%lf\n",*qx_min, *qz_min);
	if(isnan(qx)==TRUE || isnan(qz)==TRUE || isnan(qx_tot)==TRUE || isnan(qz_tot)==TRUE) printf("ole!\n");
	*qx_max = *qx_min;
	*qx_tot_max = *qx_tot_min;
	*qz_max = *qz_min;
	*qz_tot_max = *qz_tot_min;
	while(!feof(rfile)) {
	    if(fgets(buf, MAX_CHARINLINE-1, rfile) == NULL) break;
	    buf[MAX_CHARINLINE-1]='\0';

	    if(sscanf(buf, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le", &e_tot, &qx, &qz, &qx_tot, &qz_tot, &temp, &temp, &temp, &temp, &temp, &temp, &temp, &temp) != 13) errorstop("strange...");
		if(isnan(qx)==TRUE || isnan(qz)==TRUE || isnan(qx_tot)==TRUE || isnan(qz_tot)==TRUE) printf("ole!\n");
		*qx_min = MIN(qx,*qx_min);
		*qx_max = MAX(qx,*qx_max);
		*qz_min = MIN(qz,*qz_min);
		*qz_max = MAX(qz,*qz_max);
		*qx_tot_min = MIN(qx_tot,*qx_tot_min);
		*qx_tot_max = MAX(qx_tot,*qx_tot_max);
		*qz_tot_min = MIN(qz_tot,*qz_tot_min);
		*qz_tot_max = MAX(qz_tot,*qz_tot_max);
	}
	fclose(rfile);
	*qx_min = floor(*qx_min);
	*qx_max = floor(*qx_max)+1;
	*qz_min = floor(*qz_min);
	*qz_max = floor(*qz_max)+1;
	*qx_tot_min = floor(*qx_tot_min);
	*qx_tot_max = floor(*qx_tot_max)+1;
	*qz_tot_min = floor(*qz_tot_min);
	*qz_tot_max = floor(*qz_tot_max)+1;
}

extern void get_space_range(char *aspectfile, double *xmin, double *xmax, double *ymin, double *ymax) {
	char buf[MAX_CHARINLINE];
	double x,y;
	FILE *rfile = NULL;
	
	rfile = fopen(aspectfile, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	
	while(!feof(rfile)) {
	    
	    if(fgets(buf, MAX_CHARINLINE-1, rfile) == NULL) break;
	    buf[MAX_CHARINLINE-1]='\0';

	    if(sscanf(buf, "%le %le", &x, &y) != 2) continue;
		*xmin = MIN(x, *xmin);
		*xmax = MAX(x, *xmax);
		*ymin = MIN(y, *ymin);
		*ymax = MAX(y, *ymax);
	}
	fclose(rfile);
}

extern void test_field_maxwell(double x, double y, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	double bx, by, bz, bxpdx, bxpdy, bxpdz, bypdy, bypdx, bypdz, bzpdz, bzpdx, bzpdy;
	 
	get_bfield(x, y, z, &bx, &by, &bz, cell,add_contribution_comp);
	get_bfield(x+dx, y, z, &bxpdx, &bypdx, &bzpdx, cell,add_contribution_comp);
	get_bfield(x, y+dy, z, &bxpdy, &bypdy, &bzpdy, cell,add_contribution_comp);
	get_bfield(x, y, z+dz, &bxpdz, &bypdz, &bzpdz, cell,add_contribution_comp);
	
	printf("TEST IF FIELD SATISFIES DivB = 0\n");
	printf("divB = %le\n", (bxpdx - bx)/dx+(bypdy - by)/dy+(bzpdz - bz)/dz);
	printf("TEST IF FIELD SATISFIES RotB = 0\n");
	printf("dBzdx - dBxdz = %le, while dBzdx = %le\n", (bzpdx - bz)/dx - (bxpdz - bx)/dz, (bzpdx - bz)/dx);
	printf("dBzdy - dBydz = %le, while dBzdy = %le\n", (bzpdy - bz)/dy - (bypdz - by)/dz, (bzpdy - bz)/dy);
	printf("dBxdy - dBydx = %le, while dBxdy = %le\n", (bxpdy - bx)/dy - (bypdx - by)/dx, (bxpdy - bx)/dy);
}

extern void test_field_maxwell2(double x, double y, double z, double dx, double dy, double dz, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	double bx, by, bz, bxpdx, bxpdy, bxpdz, bypdy, bypdx, bypdz, bzpdz, bzpdx, bzpdy, bxmdx, bxmdy, bxmdz, bymdy, bymdx, bymdz, bzmdz, bzmdx, bzmdy;
	 
	get_bfield(x, y, z, &bx, &by, &bz, cell,add_contribution_comp);
	get_bfield(x+dx, y, z, &bxpdx, &bypdx, &bzpdx, cell,add_contribution_comp);
	get_bfield(x, y+dy, z, &bxpdy, &bypdy, &bzpdy, cell,add_contribution_comp);
	get_bfield(x, y, z+dz, &bxpdz, &bypdz, &bzpdz, cell,add_contribution_comp);
	get_bfield(x-dx, y, z, &bxmdx, &bymdx, &bzmdx, cell,add_contribution_comp);
	get_bfield(x, y-dy, z, &bxmdy, &bymdy, &bzmdy, cell,add_contribution_comp);
	get_bfield(x, y, z-dz, &bxmdz, &bymdz, &bzmdz, cell,add_contribution_comp);
	
	printf("TEST IF FIELD SATISFIES DivB = 0\n");
	printf("divB = %le\n", (bxpdx - bx)/dx+(bypdy - by)/dy+(bzpdz - bz)/dz);
	printf("TEST 2 IF FIELD SATISFIES DivB = 0\n");
	printf("divB = %le\n", (bxpdx - bxmdx)/(2*dx) + (bypdy - bymdy)/(2*dy) + (bzpdz - bzmdz)/(2*dz));
	printf("TEST IF FIELD SATISFIES RotB = 0\n");
	printf("dBzdx - dBxdz = %le, while dBzdx = %le\n", (bzpdx - bz)/dx - (bxpdz - bx)/dz, (bzpdx - bz)/dx);
	printf("dBzdy - dBydz = %le, while dBzdy = %le\n", (bzpdy - bz)/dy - (bypdz - by)/dz, (bzpdy - bz)/dy);
	printf("dBxdy - dBydx = %le, while dBxdy = %le\n", (bxpdy - bx)/dy - (bypdx - by)/dx, (bxpdy - bx)/dy);
	printf("TEST 2 IF FIELD SATISFIES RotB = 0\n");
	printf("dBzdx - dBxdz = %le, while dBzdx = %le\n", (bzpdx - bzmdx)/(2*dx) - (bxpdz - bxmdz)/(2*dz), (bzpdx - bzmdx)/(2*dx));
	printf("dBzdy - dBydz = %le, while dBzdy = %le\n", (bzpdy - bzmdy)/(2*dy) - (bypdz - bymdz)/(2*dz), (bzpdy - bzmdy)/(2*dy));
	printf("dBxdy - dBydx = %le, while dBxdy = %le\n", (bxpdy - bxmdy)/(2*dy) - (bypdx - bymdx)/(2*dx), (bxpdy - bxmdy)/(2*dy));
}

// adjust the fields of the element to keep the same FD ratio, with the particle on the wanted closed orbit (r0)
//warning!! this box only works with a 1-element lattice!
// works also with straights : r0 = x0 and k = novero.
extern int adjust_r0(struct Lattice *latt, struct Particle *part, double r_co, double *r0, double eps_clo)
{
	int i, j;
	double rnew;
	printf("adjust_r0:\n\n");
	for(i = 0; i < 10; i++) {
		if(find_closed_orbite_xxp(part, &rnew, &(part->ux), &(part->uy), eps_clo, latt, YES) == TRUE) {
			//part->x = rnew;
			part->x = r_co;
			//tune_calc_matrix(part, &nux, &nuz, &betax, &alphax, &betaz, &alphaz, 1.e-4, 1.e-5, 1.e-4, 1.e-5, latt, part_cross_latt, NO, NULL);
			if (fabs(rnew - r_co) < eps_clo*10.) {
				printf("\n\tr0adjust = %.8f [m]\n", latt->cell[0].mpara[0][1]);
				*r0 = latt->cell[0].mpara[0][1];
				return TRUE;
			}
			
			for(j = 0; j < latt->cell[0].nbcomp; j++) {
				latt->cell[0].mpara[j][1] -= (rnew-r_co);
			}
			//printf("rnew = %lf, r_co = %lf\n", rnew, r_co);
			//for(j = 0; j < latt->cell[0].nbcomp; j++) printf("latt->cell[0].mpara[%i][1] = %lf\n", j, latt->cell[0].mpara[j][1]); 
		}
		else {
			printf("\n \nproblem in adjust_r0, closed orbit not found!\n");
			return FALSE;
		}
	}
	printf("\n \nin adjust_r0 precision not achieved: eps = %le, increase the number of turns or decrease demanded precision\n \n", fabs(rnew - r_co));
	printf("\n\n\t\t\tr0adjust = %lf [m]\n\n\n", latt->cell[0].mpara[0][1]);
	return FALSE;
}

// reverse the lattice
// ! works with one cell lattices !
extern void reverse_latt(struct Lattice *new_latt, struct Lattice *latt)
{
	int n,l;
	
	new_latt->periodicity = 1;
	new_latt->nbcell = 1;
	new_latt->cell = alloccell(new_latt->nbcell);
	new_latt->cell[0].deltar = 0; //default value
	new_latt->cell[0].nbcomp = 1; //default value
	new_latt->cell[0].instrutype = NO; //default value
	//set cell.framework
	new_latt->cell[0].framework = latt->cell[0].framework;
	
	strcpy((new_latt->cell[0].keyword), (latt->cell[0].keyword));
	//printf("  --> Cell[%i] type: \"%s\"\n",nb_refcell+m, latt->cell[nb_refcell+m].keyword);
	if(strcmp(latt->cell[0].keyword, "ffag-r-he") == 0 ||
	   strcmp(latt->cell[0].keyword, "ffag-r-lin") == 0 || 
	   strcmp(latt->cell[0].keyword, "ffag-r-enge") == 0 ||
	   strcmp(latt->cell[0].keyword, "ffag-s-lin") == 0 || 
	   strcmp(latt->cell[0].keyword, "ffag-s-he") == 0 ||
	   strcmp(latt->cell[0].keyword, "ffag-s-enge") == 0 ||
	   strcmp(latt->cell[0].keyword, "ffag-sdl-lin") == 0) {
		new_latt->cell[0].nbcomp = latt->cell[0].nbcomp;
		//allocate memory
		new_latt->cell[0].mpara = allocmatrix(new_latt->cell[0].nbcomp, SIZE_MPARA);
		new_latt->cell[0].efben = allocmatrix(new_latt->cell[0].nbcomp, SIZE_EFB);
		new_latt->cell[0].efbex = allocmatrix(new_latt->cell[0].nbcomp, SIZE_EFB);
		
		//allocate memory for errors
		new_latt->cell[0].alierror = allocmatrix(new_latt->cell[0].nbcomp, SIZE_ALIERR);
		
		new_latt->cell[0].stepsize = latt->cell[0].stepsize;
		new_latt->cell[0].deltar = latt->cell[0].deltar;
		new_latt->cell[0].collim.rmin = latt->cell[0].collim.rmin;
		new_latt->cell[0].collim.rmax = latt->cell[0].collim.rmax;
		new_latt->cell[0].collim.zmin = latt->cell[0].collim.zmin;
		new_latt->cell[0].collim.zmax = latt->cell[0].collim.zmax;
		new_latt->cell[0].boun.thmax = latt->cell[0].boun.thmax;
		new_latt->cell[0].boun.ymax = latt->cell[0].boun.ymax;
		for(n = 0; n < latt->cell[0].nbcomp; n++) {
			new_latt->cell[0].mpara[n][0] = latt->cell[0].boun.ymax - latt->cell[0].mpara[latt->cell[0].nbcomp-1-n][0];
			for(l=1;l<SIZE_MPARA;l++) new_latt->cell[0].mpara[n][l] = latt->cell[0].mpara[latt->cell[0].nbcomp-1-n][l];
			for(l=0;l<SIZE_EFB;l++) new_latt->cell[0].efben[n][l] = latt->cell[0].efben[latt->cell[0].nbcomp-1-n][l];
			for(l=0;l<SIZE_EFB;l++) new_latt->cell[0].efbex[n][l] = latt->cell[0].efbex[latt->cell[0].nbcomp-1-n][l];
		}
	}
}

extern void compute_max_angle_muon_ring(struct Particle *part)
{
	double x, y;
	//printf("max_angle=%lf\n", max_angle);
	if(part->status == ALIVE) {
		x = part->x;
		y = part->y;
		fwk_pt_loctoglob(&x, &y, &(part->fwk));
		if(max_angle < fabs(part->ux)) max_angle = fabs(part->ux);
	}
}

extern double rec_mul_minus(double a, int i)
{
	return (i==0) ? a : rec_mul_minus(a,i-1)*(a-i);
}

extern double rec_mul_plus(double a, int i)
{
	return (i==0) ? a : rec_mul_plus(a,i-1)*(a+i);
}

//            ellipse1.dc2, Purpose:      Routine which gives an initial estimate for an ellipse to N points in the array X,Y.
//Author:       K.G. Begeman
//Updates:      Jan 23, 1991: KGB Document created.
//ELLIPSE1   Returns 0 when successfull, and 1 on error 
// N          Number of points in X and Y.
// X          Array with X-coordinates.
// Y          Array with Y-coordinates.
// P          Estimated ellipse parameters: P(1) = radius, P(2) = inclination (0..90 degrees), P(3) = X0 (centre), 
//P(4) = Y0 (centre), P(5) = position angle of major axis (0..180 degrees) w.r.t. Y axis.
extern int fit_ellipse1(int *n, double *x, double *y, double *p)
{
	double	ellips[5];	// ellipse parameters 
	double	matrix[5][5];// the matrix 
	double	vector[5];	// the vector 
	double	xc, yc;		// centre of ellipse 
	int		i, j, k;
	
	//Due to stability problems with ellipses where the origin is far away from the centre of the ellipse, we first have to find an approximate
	for (xc = yc = 0.0, k = 0; k < (*n); k++) {
		xc += x[k];
		yc += y[k];
	}
	xc /= (double) (*n);	// mean x position
	yc /= (double) (*n);	// mean y position
	
	//Next: initialize the matrix to be inverted by invmat and the vector
	for (i = 0; i < 5; i++) {
		for (j = 0; j < 5; j++) matrix[i][j] = 0.0;
		vector[i] = 0.0;
	}
	
	//Then: fill matrix
	for (k = 0; k < (*n); k++) {
		double	xs = x[k] - xc;	// x coord. w.r.t. centre 
		double	ys = y[k] - yc;	// y coord. w.r.t. centre 
		double	xx;
		double	xy;
		double	yy;
		
		xx = xs * xs;	// x * x 
		xy = xs * ys;	// cross product 
		yy = ys * ys;	// y * y 
		matrix[0][0] += xx * xx;
		matrix[0][1] += 2.0 * xx * xy;
		matrix[0][2] += xy * xy;
		matrix[0][3] += xx * xs;
		matrix[0][4] += xy * xs;
		matrix[1][1] += 4.0 * xy * xy;
		matrix[1][2] += 2.0 * xy * yy;
		matrix[1][3] += 2.0 * xy * xs;
		matrix[1][4] += 2.0 * xy * ys;
		matrix[2][2] += yy * yy;
		matrix[2][3] += xy * ys;
		matrix[2][4] += yy * ys;
		matrix[3][3] += xx;
		matrix[3][4] += xy;
		matrix[4][4] += yy;
		vector[0] += xx;
		vector[1] += 2.0 * xy;
		vector[2] += yy;
		vector[3] += xs;
		vector[4] += ys;
	}
	//make the matrix symmetric since this saves a lot of typing
	for(i = 0; i < 5; i++) for (j = 0; j < i; j++) matrix[i][j] = matrix[j][i];
	
	//invert the matrix
	invmat( matrix, 5 );
	//if(invmat( matrix, 5 )) return( 1 ); // cannot invert matrix 
	
	//solve MATRIX * ELLIPS = VECTOR
	for (i = 0; i < 5; i++) {
		ellips[i] = 0.0;	// reset this ellipse parm.
		for (j = 0; j < 5; j++) ellips[i] += matrix[i][j] * vector[j];
	}
	//NOTE: The ellipse equation taken is AA.x^2 + 2.BB.x.y + CC.y^2 + DD.x + EE.y = 1,
	//Where AA..EE were now solved for from which now the ellipse-parameters are derived
	{
		double	aa = ellips[0];
		double	bb = ellips[1];
		double	cc = ellips[2];
		double	dd = ellips[3];
		double	ee = ellips[4];
		double	pa, pp;
		double	cospa, sinpa, sinpp;
		double	s1, s2, s3, y1, y2, y3;
		double	x0, y0;
		double	ab, al, r;
		
		pp = atan( 2.0 * bb / ( aa - cc ) );	// estimate of position angle 
		pa = 0.5 * pp;				// p.a. of an (UNDETERMINED) axis 
		cospa = cos( pa );			// cosine 
		sinpp = sin( pp );			// sine of double angle 
		sinpa = sin( pa );			// sine 
		al = 2.0 * bb / sinpp / ( aa + cc );	// auxiliary 
		r = sqrt( ( 1.0 + al ) / ( 1.0 - al ) );	// axial ratio (OR ITS RECIPROCAL)
		s1 = bb * ee - cc * dd;			// three other auxiliaries 
		s2 = bb * dd - aa * ee;
		s3 = aa * cc - bb * bb;
		x0 = s1 / s3;				// X-centre of ellipse
		y0 = s2 / s3;				// Y-centre of ellipse
		y1 = sinpa * sinpa + r * r * cospa * cospa;
		y2 = x0 * y0 * ( r * r - 1.0 ) * sinpp;
		y3 = y0 * y0 * ( cospa * cospa + r * r * sinpa * sinpa );
		ab = sqrt( y1 * ( x0 * x0 + 1.0 / aa ) + y2 + y3 ); //length of (yet undetermined) axis (A or B)
		pa = pa * 45.0 / atan( 1.0 );		// convert to degrees 
		
		//determination which axis is the long axis
		if (r < 1.0) {	// the other is 
			ab /= r;	// this one is the major axis 
			pa -= 90.0;	// so change position angle 
		} 
		else r = 1.0 / r; // the right one is, so change axial ratio 
		if (pa < 0.0) pa += 180.0;	// position angle in range 0.....180
		p[0] = ab;				// radius 
		p[1] = acos( r ) * 45.0 / atan( 1.0 );	// inclination assuming projected circle
		p[2] = x0 + xc;	// new x-position of centre
		p[3] = y0 + yc;	// new y-position of centre
		p[4] = pa;		// position angle major axis 
	}
	return 0;
}

//            ellipse2.dc2, Purpose:      Fits an ellipse to a set of X and Y positions. 
//Author:       K.G. Begeman
//Updates:      Jan 25, 1991: KGB Document created.
//ELLIPSE2     Returns number of iterations on success, else, -1: No free parameters, -2: Exceede iteration limit (50), -3: Diagonal of matrix has zeroes, -4: Matrix could not be inverted
// N            Number of positions.
// X            Array with X coordinates (in units).
// Y            Array with Y coordinates (in units).
// P            On input contains the initial estimates of the ellipse parameters. On output the fitted parameters. P contains:
//              P(1) = radius (in units), P(2) = inclination (degrees), P(3) = X centre of ellipse (units), P(4) = Y centre of ellipse (units), P(5) = Position angle (degrees)
// E            Contains the errors in the fitted parameters.
// M            Mask for free (1) or fixed (0) parameters.
extern int fit_ellipse2(int *n, double *x, double *y, double *p, double *e, int *m)
{
	double chi;		// red. chi-squared 
	double labda;	// mixing parameter
	double q;		// red. chi-squared 
	double rl[5];	// vector 
	double s[5][5], s1[5][5];	// matrices
	double mix_fac = 10.0; //mixing * factor
	double mix_start = 0.01; //start mixing parameter
	double lab_max = 1.0E+10; // max. mixing parameter
	double lab_min = 1.0E-10; // min. mixing parameter
	double tol = 0.00001; //tolerance
	int h = 0;		// iteration counter
	int max_iter = 50; //max. number of iterations
	int i;			// counter
	int ip[5];		// permutation array
	int nfree = 0;	// number of free parameters
	int r = 0;		// return value 
	
	labda = mix_start * mix_fac;	// start value
	for (i = 0; i < 5; i++) if (m[i]) ip[nfree++] = i; // fit this parameter
	if (nfree == 0 || nfree >= (*n)) errorstop("number of free parameters is wrong!");
	do {
		if (++h > max_iter) r = -2; break; // too many iterations
		chi = inimat( s, rl, x, y, (*n), p, e, ip, nfree );
		if (labda > lab_min) labda /= mix_fac;	// new labda
		r = inivec( s, s1, rl, labda, &q, p, e, x, y, (*n), ip, nfree );
		//if (r) break;				// error from inivec 
		while (q >= chi) {			// interpolation loop 
			if (labda > lab_max) break;	// leave loop 
			labda *= mix_fac;				// new labda 
			r = inivec( s, s1, rl, labda, &q, p, e, x, y, (*n), ip, nfree );
			if (r) break;				// error from inivec
		}
		if (labda <= lab_max) for (i = 0; i < 5; i++) p[i] = e[i];
		if (fabs( chi - q ) / q < tol || labda > lab_max) {
			labda = 0.0;				// for Taylor solution 
			chi = inimat( s, rl, x, y, (*n), p, e, ip, nfree );
			r = inivec( s, s1, rl, labda, &q, p, e, x, y, (*n), ip, nfree );
			if (!r) {
				r = h;				// number of iterations
				for (i = 0; i < 5; i++) {
					p[i] = e[i];			// the parameters
					e[i] = 0.0;			// reset errors
				}
				q = sqrt( q / (double) ( (*n) - nfree ) );
				for (i = 0; i < nfree; i++) e[ip[i]] = q * sqrt( s1[i][i] / s[i][i] );
			}
		}
	} while (!r);	// until error or finished 
	return r;
}

/*extern void get_tune_range(char *txtfilename, double *qx_min, double *qx_max, double *qz_min, double *qz_max, double *qx_tot_min, double *qx_tot_max, double *qz_tot_min, double *qz_tot_max) 
{
	int n, nblines;
	double e_tot, qx, qz, qx_tot, qz_tot, temp, qx_min_temp, qz_min_temp, qx_max_temp, qz_max_temp, qx_tot_min_temp, qz_tot_min_temp, qx_tot_max_temp, qz_tot_max_temp;
	FILE *rfile = NULL;
	
	nblines = get_nb_lines_file(txtfilename);
	rfile = fopen(txtfilename, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	fscanf(rfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le", 
	&e_tot, &qx_min_temp, &qz_min_temp, &qx_tot_min_temp, &qz_tot_min_temp, &temp, &temp, &temp, &temp, &temp, &temp, &temp, &temp, &temp, &temp, &temp, &temp);
	newline(rfile);
	qx_max_temp = qx_min_temp;
	qz_max_temp = qz_min_temp;
	qx_tot_max_temp = qx_tot_min_temp;
	qz_tot_max_temp = qz_tot_min_temp;
	for(n=0;n<nblines-2;n++) {
		fscanf(rfile, "%le  %le  %le  %le  %le", &e_tot, &qx, &qz, &qx_tot, &qz_tot);
		printf("%le  %le  %le  %le  %le\n", e_tot, qx, qz, qx_tot, qz_tot);
		
		fscanf(rfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le", 
		&e_tot, &qx, &qz, &qx_tot, &qz_tot, &temp, &temp, &temp, &temp, &temp, &temp, &temp, &temp, &temp, &temp, &temp, &temp);
		//fscanf(rfile, "%le  %le  %le  %le  %le  %le  %le", &ekin_ev, &qx, &qz, &ptot, &ptot2, &qx2, &qz2);
		if(qx<qx_min_temp) qx_min_temp = qx;
		else if(qx>qx_max_temp) qx_max_temp = qx;
		if(qz<qz_min_temp) qz_min_temp = qz;
		else if(qz>qz_max_temp) qz_min_temp = qz;
		if(qx_tot<qx_tot_min_temp) qx_tot_min_temp = qx_tot;
		else if(qx_tot>qx_tot_max_temp) qx_tot_max_temp = qx_tot;
		if(qz_tot<qz_tot_min_temp) qz_tot_min_temp = qz_tot;
		else if(qz_tot>qz_tot_max_temp) qz_tot_min_temp = qz_tot;
	}
	fclose(rfile);
	*qx_min = floor(qx_min_temp);
	*qx_max = floor(qx_max_temp)+1;
	*qz_min = floor(qz_min_temp);
	*qz_max = floor(qz_max_temp)+1;
	*qx_tot_min = floor(qx_tot_min_temp);
	*qx_tot_max = floor(qx_tot_max_temp)+1;
	*qz_tot_min = floor(qz_tot_min_temp);
	*qz_tot_max = floor(qz_tot_max_temp)+1;
}//*/

extern void get_betamin_file(char *txtfilename, double *betamin) {
	int n, nblines;
	double s, bx, ax, bz, az, bmin=1.0e18;
	FILE *rfile = NULL;
	
	nblines = get_nb_lines_file(txtfilename);
	rfile = fopen(txtfilename, "r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(n=0;n<nblines-1;n++) {
		fscanf(rfile, "%le	%le	%le	%le	%le", &s, &bx, &ax, &bz, &az);
		if(bmin > bx || bmin > bz) bmin = MIN(bx, bz);
	}
	*betamin = bmin;
	//printf("betamax = %lf\n", bmax);
	fclose(rfile);
}

// compute x, x', z, z' w/ respect to the closed orbit from the particle parameters 
extern void comp_phase_space_coord(double *x, double *xprime, double *z, double *zprime, double test_x, double test_z, double test_ux, double test_uy, double test_uz, double x_ref, double xprime_ref, double z_ref, double zprime_ref)
{
	double php_test;
	 
	*xprime = atan_ratio(test_ux, test_uy) - xprime_ref;
	//*x = xpart[n] - reference->x;
	*x = (test_x - x_ref)*(cos(xprime_ref) + sin(xprime_ref)*tan(*xprime));
	php_test = sqrt(test_ux*test_ux + test_uy*test_uy);
	*zprime = atan_ratio(test_uz, php_test) - zprime_ref;
	//z = zpart[n] - reference->z;
	*z = (test_z - z_ref)*(cos(zprime_ref) + sin(zprime_ref)*tan(*zprime));
}

