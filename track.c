/*
 *  track.c
 *  ringdesign
 *
 *  Copyright 2009 Kyoto University. All rights reserved.
 *
 */

#include "track.h"

extern int find_closed_orbite_x(struct Particle *part, double *x_clo, double eps_clo, struct Lattice *latt, int doyouprintf)
{
	int i, n;
	double xmin, xmax, eps;
	struct Particle test_part;
	test_part = *part;
	test_part.hat = -2;	// no acceleration, no output
	eps  = 1.;
	n	 = 0.;
	
	while (eps > eps_clo) {
		xmin = xmax = test_part.x;
		test_part = *part;
		test_part.hat = -2;	// no acceleration, no output
		test_part.x = xmin;
		
		for(i = 0; i < FCLO_CROSS; i++) {
			part_cross_latt(&test_part, latt,NULL);
			xmax = MAX(test_part.x,xmax);
			xmin = MIN(test_part.x,xmin);
			if(test_part.status != ALIVE) break;
		}
		
		test_part.x = (xmin+xmax)/2.;
		eps = xmax - xmin;
		if(doyouprintf == YES) {
			CLRSCR();
			//printf("%c[2K\r", 33);
			printf("try # %i, eps = %le, xmin = %lf, xmax = %lf", n, eps, xmin, xmax);
			fflush(stdout);
		}
		if(n > FCLO_TRIES || test_part.status != ALIVE) {
			if(doyouprintf == YES) printf("\n!!! Cannot find closed orbit, eps_clo = %le [m]!!!, part. status = %i\n", eps, test_part.status);
			return FALSE;
		}
		n++;
	}
	*x_clo = test_part.x;
	if(doyouprintf == YES) {
		if(print_color==YES) COLOR("1;34");
		printf("\nClosed orbit found: x_clo = %.8f [m], xprime_clo forced to be 0, \neps = %le [m]\n", *x_clo, eps);
		if(print_color==YES) COLOR("0");
	}
	return TRUE;
}

extern int find_closed_orbite_xz(struct Particle *part, double *x_clo, double *z_clo, double eps_clo, struct Lattice *latt, int doyouprintf)
{
	int i, n;
	double xmin, xmax, zmin, zmax, eps;
	struct Particle test_part;
	test_part = *part;
	test_part.hat = -2;	// no acceleration, no output
	eps  = 1.;
	n	 = 0.;
	
	while (eps > eps_clo) {
		xmin = xmax = test_part.x;
		zmin = zmax = test_part.z;
		test_part = *part;
		test_part.hat = -2;	// no acceleration, no output
		test_part.x = xmin;
		test_part.z = zmin;
		
		for(i = 0; i < FCLO_CROSS; i++) {
			part_cross_latt(&test_part, latt,NULL);
			xmax = MAX(test_part.x,xmax);
			xmin = MIN(test_part.x,xmin);
			zmax = MAX(test_part.z,zmax);
			zmin = MIN(test_part.z,zmin);
			if(test_part.status != ALIVE) break;
		}
		
		test_part.x = (xmin+xmax)/2.; 
		test_part.z = (zmin+zmax)/2.;
		eps = sqrt((xmax - xmin)*(xmax - xmin) + (zmax - zmin)*(zmax - zmin));
		//eps = zmax-zmin;
		if(doyouprintf == YES) {
			//CLRSCR();
			//printf("%c[2K\r", 33);
			printf("try # %i, eps = %le, x:min,max: (%lf, %lf), z:min,max: (%lf, %lf)", n, eps, xmin, xmax, zmin, zmax);
			printf("\n");
			//fflush(stdout);
		}
		if(n > FCLO_TRIES || test_part.status != ALIVE) {
			if(doyouprintf == YES) printf("\n!!! Cannot find closed orbit, eps_clo = %le [m]!!!, part. status = %i\n", eps, test_part.status);
			return FALSE;
		}
		n++;
	}
	*x_clo = test_part.x;
	*z_clo = test_part.z;
	if(doyouprintf == YES) {
		if(print_color==YES) COLOR("1;34");
		printf("\nClosed orbit found: x_clo = %.8f [m], z_clo = %.8f [m], xprime_clo forced to be 0, \neps = %le [m]\n", *x_clo, *z_clo, eps);
		if(print_color==YES) COLOR("0");
	}
	return TRUE;
}

extern int find_closed_orbite_xxp(struct Particle *part, double *x_clo, double *ux_clo, double *uy_clo, double eps_clo, struct Lattice *latt, int doyouprintf)
{
	int i, n;
	double xmin, xmax, uxmin, uxmax, utot_xy, eps;
	struct Particle test_part;
	
	
	test_part = *part;
	utot_xy = sqrt(part->ux*part->ux + part->uy*part->uy);
	test_part.hat = -2;	// no acceleration, no output
	eps  = 1.;
	n	 = 0.;
	while (eps > eps_clo) {
		xmin = xmax = test_part.x;
		uxmin = uxmax = test_part.ux;
		test_part = *part;
		test_part.hat = -2;	// no acceleration, no output
		test_part.x = xmin;
		test_part.ux = uxmin;
		test_part.uy = sqrt(utot_xy*utot_xy - uxmin*uxmin);
		for(i = 0; i < FCLO_CROSS; i++) {
			part_cross_latt(&test_part, latt,NULL);
			xmax = MAX(test_part.x, xmax);
			xmin = MIN(test_part.x,xmin);
			uxmax = MAX(test_part.ux, uxmax);
			uxmin = MIN(test_part.ux, uxmin);
			if(test_part.status != ALIVE) break;
		}
		
		test_part.x = (xmin+xmax)/2.;
		eps = sqrt((xmax - xmin)*(xmax - xmin) + (uxmax - uxmin)*(uxmax - uxmin));
		//eps = xmax - xmin;
		test_part.ux = (uxmin + uxmax)/2.;
		test_part.uy = sqrt(utot_xy*utot_xy - uxmin*uxmin);
		if(doyouprintf == YES) {
			CLRSCR();
			printf("try # %i, eps = %le, x = %lf, %lf, x' = %le, %le", n, eps, xmin, xmax, asin(uxmin)*180./PI, asin(uxmax)*180./PI);
			fflush(stdout);
		}
		if(n > FCLO_TRIES || test_part.status != ALIVE) {
			if(doyouprintf == YES) printf("\n!!! Cannot find closed orbit, eps_clo = %le [m]!!!, part. status = %i\n", eps, test_part.status);
			return FALSE;
		}
		n++;
	}
	*x_clo = test_part.x;
	*ux_clo = test_part.ux;
	*uy_clo = test_part.uy;
	if(doyouprintf == YES) {
		if(print_color==YES) COLOR("1;34");
		printf("\nClosed orbit found: x_clo = %.8f [m], xprime_clo = %le [deg], \neps = %le [m]\n", *x_clo, atan(*ux_clo/(*uy_clo))*180./(PI), eps);
		if(print_color==YES) COLOR("0");
	}
	return TRUE;
}

extern int find_closed_orbite_xxp_zzp(struct Particle *part, double *x_clo, double *ux_clo, double *uy_clo, double *z_clo, double *uz_clo, double eps_clo, struct Lattice *latt, int doyouprintf)
{
	int i, n;
	double xmin, xmax, uxmin, uxmax, eps, zmin, zmax, uzmin, uzmax;
	struct Particle test_part;
	
	test_part = *part;
	//utot_xy = sqrt(part->ux*part->ux + part->uy*part->uy + part->uz*part->uz);
	test_part.hat = -2;	// no acceleration, no output
	eps  = 1.;
	n	 = 0.;
	
	while (eps > eps_clo) {
		xmin = xmax = test_part.x;
		uxmin = uxmax = test_part.ux;
		zmin = zmax = test_part.z;
		uzmin = uzmax = test_part.uz;
		test_part = *part;
		test_part.hat = -2;	// no acceleration, no output
		test_part.x = xmin;
		test_part.ux = uxmin;
		test_part.z = zmin;
		test_part.uz = uzmin;
		//test_part.uy = sqrt(utot_xy*utot_xy - uxmin*uxmin- uzmin*uzmin);
		test_part.uy = sqrt(1.0 - uxmin*uxmin- uzmin*uzmin);
		
		for(i = 0; i < FCLO_CROSS; i++) {
			part_cross_latt(&test_part,latt,NULL);
			xmax = MAX(test_part.x, xmax);
			xmin = MIN(test_part.x,xmin);
			uxmax = MAX(test_part.ux, uxmax);
			uxmin = MIN(test_part.ux, uxmin);
			zmax = MAX(test_part.z, zmax);
			zmin = MIN(test_part.z,zmin);
			uzmax = MAX(test_part.uz, uzmax);
			uzmin = MIN(test_part.uz, uzmin);
			if(test_part.status != ALIVE) break;
		}
		
		test_part.x = (xmin+xmax)/2.;
		//eps = xmax - xmin + zmax - zmin;
		eps = sqrt((xmax - xmin)*(xmax - xmin) + (zmax - zmin)*(zmax - zmin) + (uxmax - uxmin)*(uxmax - uxmin) + (uzmax - uzmin)*(uzmax - uzmin));
		//eps = fabs(xmax - xmin) + fabs(zmax - zmin) + (fabs(uxmax - uxmin) + fabs(uzmax - uzmin))/(test_part.uy);
		test_part.ux = (uxmin + uxmax)/2.;
		test_part.z = (zmin+zmax)/2.;
		test_part.uz = (uzmin+uzmax)/2.;
		test_part.uy = sqrt(1.0 - uxmin*uxmin-uzmin*uzmin);
		if(doyouprintf == YES) {
			CLRSCR();
			printf("try # %i, eps = %le, xmin = %lf, xmax = %lf, zmin = %lf, zmax = %lf", n, eps, xmin, xmax, zmin, zmax);
			fflush(stdout);
		}
		
		if(n > FCLO_TRIES || test_part.status != ALIVE) {
			if(doyouprintf == YES) printf("\n!!! Cannot find closed orbite, eps_clo = %le [m]!!!, part. status = %i\n", eps, test_part.status);
			return FALSE;
		}
		n++;
	}
	*x_clo = test_part.x;
	*ux_clo = test_part.ux;
	*uy_clo = test_part.uy;
	*z_clo = test_part.z;
	*uz_clo = test_part.uz;
	if(doyouprintf == YES) {
		if(print_color==YES) COLOR("1;34");
		printf("\nClosed orbit found: x_clo = %.10f [m], xprime_clo = %.10e [deg], \nz_clo = %.10f [m], zprime_clo = %.10e [deg], eps = %le [m]\n", *x_clo, atan(*ux_clo/(*uy_clo))*180./(PI), *z_clo, atan(*uz_clo/(sqrt((*uy_clo)*(*uy_clo)+(*ux_clo)*(*ux_clo))))*180./(PI), eps);
		print_part_para("test_part to find c.o.", &test_part);
		if(print_color==YES) COLOR("0");
	}
	
	return TRUE;
}

extern int find_co_ellipse(struct Particle *part, double *x_clo, double *z_clo, double *ux_clo, double *uy_clo, double *uz_clo, double eps_clo, int nb_points, struct Lattice *latt, int doyouprintf)
{
	int i, n, m[5], a[5];
	double utot_xy, eps, x[nb_points], ux[nb_points], z[nb_points], uz[nb_points], p[5], e[5], q[5], f[5];
	struct Particle test_part;
	
	for(i=0;i<5;i++) {
		p[i] = 0.;
		e[i] = 0.;
		m[i] = 1;
		q[i] = 0.;
		f[i] = 0.;
		a[i] = 1;
	}
	
	test_part = *part;
	utot_xy = sqrt(part->ux*part->ux + part->uy*part->uy);
	test_part.hat = -2;	// no acceleration, no output
	p[2] = test_part.x;
	p[3] = test_part.ux;
	q[2] = test_part.z;
	q[3] = test_part.uz;
	eps  = 1.;
	n	 = 0.;
	
	while (eps > eps_clo) {
		test_part = *part;
		test_part.hat = -2;	// no acceleration, no output
		test_part.x = p[2];
		test_part.ux = p[3];
		test_part.z = q[2];
		test_part.uz = q[3];
		if(fabs(test_part.z)<TINYLENGTH && fabs(test_part.uz)<TINYDIMLESS) test_part.uy = sqrt(utot_xy - p[3]*p[3]);
		else test_part.uy = sqrt(1. - p[3]*p[3] - q[3]*q[3]);
		for(i = 0; i < nb_points; i++) {
			part_cross_latt(&test_part, latt, NULL);
			x[i] = test_part.x;
			ux[i] = test_part.ux;
			z[i] = test_part.z;
			uz[i] = test_part.uz;
			if(test_part.status != ALIVE) break;
		}
		fit_ellipse1(&nb_points, x, ux, p);
		fit_ellipse2(&nb_points, x, ux, p, e, m);
		if(fabs(test_part.z)<TINYLENGTH && fabs(test_part.uz)<TINYDIMLESS) eps = fabs(test_part.x - p[2]);
		else {
			fit_ellipse1(&nb_points, z, uz, q);
			fit_ellipse2(&nb_points, z, uz, q, f, a);
			eps = fabs(test_part.x - p[2]) + fabs(test_part.z - q[2]);
		} 
		
		if(doyouprintf == YES) {
			CLRSCR();
			printf("try # %i, eps = %le, x = %lf, ux = %lf, z = %lf, uz = %lf", n, eps, p[2], p[3], q[2], q[3]);
			fflush(stdout);
		}
		if(n > FCLO_TRIES || test_part.status != ALIVE) {
			if(doyouprintf == YES) printf("\n!!! Cannot find closed orbit, eps_clo = %le [m]!!!, part. status = %i\n", eps, test_part.status);
			return FALSE;
		}
		n++;
	}
	*x_clo = test_part.x;
	*ux_clo = test_part.ux;
	*uy_clo = test_part.uy;
	if(fabs(test_part.z)<TINYLENGTH && fabs(test_part.uz)<TINYDIMLESS) {
		*z_clo = 0.;
		*uz_clo = 0.;
	}
	else {
		*z_clo = test_part.z;
		*uz_clo = test_part.uz;
	}
	
	if(doyouprintf == YES) {
		if(print_color==YES) COLOR("1;34");
		printf("\nClosed orbit found: x_clo = %.8f [m], xprime_clo = %le [deg], z_clo = %.8f [m], zprime_clo = %le [deg], \neps = %le [m]\n", *x_clo, atan(*ux_clo/(*uy_clo))*180./(PI), *z_clo, atan(*uz_clo/(*uy_clo))*180./(PI), eps);
		if(print_color==YES) COLOR("0");
	}
	return TRUE;
}

extern void apply_alierror_comp(struct Cell *cell, int i, double x, double y, double z, double *xwerr, double *ywerr, double *zwerr)
{
	double err_angle, ax,ay,az,ux,uy,uz,co,si,xtemp,ytemp,ztemp;
	
	err_angle = cell->alierror[i][4];
	ax = cell->alierror[i][5];
	ay = cell->alierror[i][6];
	az = cell->alierror[i][7];
	ux = cell->alierror[i][8];
	uy = cell->alierror[i][9];
	uz = cell->alierror[i][10];
	//debug test
	if(fabs(ux*ux + uy*uy + uz*uz - 1) > 1.e-12) errorstop("ERROR in apply_alierror. vector to define rotation axe is not a unit vector. This is forbidden!!");
	//rotation
	co = cos(-err_angle); //"-" signe because to simulate a rotation of the magnetic element of an angle +err_angle, you need to ask the field in the non-rotated magnetic element after a rotation of angle -err_angle.
	si = sin(-err_angle); //"-" signe because to simulate a rotation of the magnetic element of an angle +err_angle, you need to ask the field in the non-rotated magnetic element after a rotation of angle -err_angle.
	xtemp = x - ax;
	ytemp = y - ay;
	ztemp = z - az;
	*xwerr = xtemp*(ux*ux+(1-ux*ux)*co)		+ ytemp*(ux*uy*(1-co) - uz*si)	+ ztemp*(ux*uz*(1-co) + uy*si);
	*ywerr = xtemp*(ux*uy*(1-co) + uz*si)	+ ytemp*(uy*uy + (1-uy*uy)*co)	+ ztemp*(uy*uz*(1-co) - ux*si);
	*zwerr = xtemp*(ux*uz*(1-co) - uy*si)	+ ytemp*(uy*uz*(1-co) + ux*si)	+ ztemp*(uz*uz + (1-uz*uz)*co);
	*xwerr += ax;
	*ywerr += ay;
	*zwerr += az;
	//translation
	*xwerr -= cell->alierror[i][1]; //"-" signe: same reason than for the rotation angle signe, see few lines above.
	*ywerr -= cell->alierror[i][2];
	*zwerr -= cell->alierror[i][3];
}

extern int acceptancex(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double axmin, double axmax, double axstep, double az)
{
	int i;
	double x, xprime, z, zprime, ampx, php_ref, atanrefx, atanrefz;
	struct Particle test_part;
	FILE *rfile;
	if(outfilename != NULL) rfile = fopen(outfilename,"a");
	
	php_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	atanrefx = atan_ratio((reference->ux), (reference->uy));
	atanrefz = atan_ratio(reference->uz, php_ref);
	printf("acceptancex:");
	
	for (ampx = axmin; ampx <= axmax; ampx += axstep) {
		test_part = *reference;
		test_part.s = 0;
		test_part.z += az;
		test_part.x += ampx;
		test_part.hat = -2;
		printf("\namp = %lf\n", ampx);
		
		for (i = 0; i < nbpass; i++) {
			CLRSCR();
			printf("pass number %i\t", i);
			fflush(stdout);
			if (test_part.uy > 0) {
				comp_phase_space_coord(&x, &xprime, &z, &zprime,  test_part.x, test_part.z, test_part.ux, test_part.uy, test_part.uz, reference->x, atanrefx, reference->z, atanrefz);
			}
			else {
				printf("\tin acceptancex, test_part.uy <= 0!, particle is going backwards!\n");
				if(outfilename != NULL) fclose(rfile);
				return test_part.status;
			}
			if(outfilename != NULL) fprintf(rfile, "%le  %le  %le  %le  %le\n", x, xprime, z, zprime, ampx);
			part_cross_latt(&(test_part), latt,NULL);
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

extern int acceptancez(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double azmin, double azmax, double azstep, double ax)
{
	int i;
	double x, xprime, z, zprime, ampz, php_ref, atanrefx, atanrefz;
	struct Particle test_part;
	FILE *rfile;
	if(outfilename != NULL) rfile = fopen(outfilename,"a");
	php_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	atanrefx = atan_ratio((reference->ux), (reference->uy));
	atanrefz = atan_ratio(reference->uz, php_ref);
	printf("acceptancez:");
	for (ampz = azmin; ampz <= azmax; ampz += azstep) {
		test_part = *reference;
		test_part.s = 0;
		test_part.x += ax;
		test_part.z += ampz;
		test_part.hat = -2;
		printf("\namp = %lf\n", ampz);
		
		for (i = 0; i < nbpass; i++) {
			CLRSCR();
			printf("pass %i\t", i);
			fflush(stdout);
			if (test_part.uy > 0) {
				comp_phase_space_coord(&x, &xprime, &z, &zprime,  test_part.x, test_part.z, test_part.ux, test_part.uy, test_part.uz, reference->x, atanrefx, reference->z, atanrefz);
			}
			else {
				printf("\tin acceptancez, test_part.uy <= 0!, particle is going backwards!\n");
				if(outfilename != NULL) fclose(rfile);
				return test_part.status;
			}
			if(outfilename != NULL) fprintf(rfile, "%le  %le  %le  %le  %le\n", x, xprime, z, zprime, ampz);
			part_cross_latt(&(test_part), latt,NULL);
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

extern int acceptancezprime(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double zprimemin, double zprimemax, double zprimestep, double ax)
{
	int i;
	double asquare, x, xprime, z, zprime, ampz, php_ref, atanrefx, atanrefz;
	struct Particle test_part;
	FILE *rfile;
	if(outfilename != NULL) rfile = fopen(outfilename,"a");
	php_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	atanrefx = atan_ratio((reference->ux), (reference->uy));
	atanrefz = atan_ratio(reference->uz, php_ref);
	asquare = reference->ux*reference->ux/(reference->uy*reference->uy);
	
	for (ampz = zprimemin; ampz <= zprimemax; ampz += zprimestep) {
		test_part = *reference;
		test_part.s = 0;
		test_part.x += ax;
		test_part.hat = -2;
		test_part.uy = sign(reference->uy)*(1/(sqrt((asquare+1)*(tan(ampz)*tan(ampz)+1))));
		test_part.ux = sign(reference->ux)*sqrt(asquare)*test_part.uy;
		test_part.uz = test_part.uy*sqrt(1+asquare)*tan(ampz);
		
		for (i = 0; i < nbpass; i++) {
			printf("pass number %i\n", i);
			if (test_part.uy > 0) {
				comp_phase_space_coord(&x, &xprime, &z, &zprime,  test_part.x, test_part.z, test_part.ux, test_part.uy, test_part.uz, reference->x, atanrefx, reference->z, atanrefz);
			}
			else {
				printf("test_part.uy <= 0!, particle is going backwards!\n");
				if(outfilename != NULL) fclose(rfile);
				return FALSE;
			}
			if(outfilename != NULL) fprintf(rfile, "%le  %le  %le  %le  %le\n", x, xprime, z, zprime, ampz);
			part_cross_latt(&(test_part), latt,NULL);
			if(test_part.status != ALIVE) {
				if(outfilename != NULL) fclose(rfile);
				return test_part.status;
			}
		}
	}	
	
	if(outfilename != NULL) fclose(rfile);
	return test_part.status;
}

extern void acceptancex_auto(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double axmin, double *axmax, double axstep, double az) 
{
	int i, n, doyouprint=YES;
	double ax, xpart[nbpass], zpart[nbpass], ux[nbpass], uy[nbpass], uz[nbpass], x, z, xprime, zprime, php_ref, atanrefx, atanrefz;
	FILE *wfile;
	if(outfilename != NULL) wfile = fopen(outfilename, "w");
	
	php_ref = sqrt(reference->ux*reference->ux + reference->uy*reference->uy);
	atanrefx = atan_ratio((reference->ux), (reference->uy));
	atanrefz = atan_ratio(reference->uz, php_ref);
	
	ax = axmin;
	i = 0;
	*axmax = 0.;
	for(;;) {
		ax = axmin+i*axstep;
		CLRSCR();
		if(doyouprint==YES) printf("ax = %lf\t", ax);
		fflush(stdout);
		if(track_n_turns_amp(xpart, zpart, ux, uy, uz, reference, latt, nbpass, ax, az) != ALIVE) {
			if(doyouprint==YES) printf("NO!\n");
			break;
		}
		else {
			if(doyouprint==YES) printf("OK!\n");
			//*axmax = ax;
			if(outfilename != NULL) {
				for(n=0;n<nbpass;n++) {
					comp_phase_space_coord(&x, &xprime, &z, &zprime, xpart[n], zpart[n], ux[n], uy[n], uz[n], reference->x, atanrefx, reference->z, atanrefz);
					fprintf(wfile, "%le  %le  %le  %le  %le\n", x, xprime, z, zprime, ax);
				}
			}
		}
		i++;
	}
	if(i==0) *axmax = 0.;
	else *axmax = ax-axstep;
	if(outfilename != NULL) fclose(wfile);
	if(doyouprint==YES) printf("final: ax = %le\n",*axmax);
}

extern int track_n_turns_amp(double x[],double z[], double ux[], double uy[], double uz[], struct Particle *reference, struct Lattice *latt, int nbpass, double ax, double az) 
{
	int i;
	struct Particle test_part;
	test_part = *reference;
	test_part.s = 0;
	test_part.z += az;
	test_part.x += ax;
	for(i=0;i<nbpass;i++) {
		if(test_part.status != ALIVE) return test_part.status;
		else {
			x[i] = test_part.x;
			z[i] = test_part.z;
			ux[i] = test_part.ux;
			uy[i] = test_part.uy;
			uz[i] = test_part.uz;
		}
		part_cross_latt(&(test_part), latt,NULL);
	}
	return test_part.status;
}

extern void acceptancez_auto(char *outfilename, struct Particle *reference, struct Lattice *latt, int nbpass, double azmin, double *azmax, double azstep, double ax) 
{
	int i, n, doyouprint=YES;
	double az, xpart[nbpass], zpart[nbpass], ux[nbpass], uy[nbpass], uz[nbpass], x, z, xprime, zprime, php_ref, atanrefx, atanrefz;
	FILE *wfile;
	if(outfilename != NULL) wfile = fopen(outfilename, "w");
	
	php_ref = sqrt(reference->ux*reference->ux + reference->uy*reference->uy);
	atanrefx = atan_ratio((reference->ux), (reference->uy));
	atanrefz = atan_ratio(reference->uz, php_ref);
	
	az = azmin;
	i = 0;
	*azmax = 0.;
	for(;;) {
		az = azmin+i*azstep;
		CLRSCR();
		if(doyouprint==YES) printf("az = %lf\t", az);
		fflush(stdout);
		if(track_n_turns_amp(xpart, zpart, ux, uy, uz, reference, latt, nbpass, ax, az) != ALIVE) {
			if(doyouprint==YES) printf("NO!\n");
			break;
		}
		else {
			if(doyouprint==YES) printf("OK!\n");
			//*azmax = az;
			if(outfilename != NULL) {
				for(n=0;n<nbpass;n++) {
					comp_phase_space_coord(&x, &xprime, &z, &zprime, xpart[n], zpart[n], ux[n], uy[n], uz[n], reference->x, atanrefx, reference->z, atanrefz);
					fprintf(wfile, "%le  %le  %le  %le  %le\n", x, xprime, z, zprime, az);
				}
			}
		}
		i++;
	}
	if(i==0) *azmax = 0.;
	else *azmax = az-azstep;
	if(outfilename != NULL) fclose(wfile);
	if(doyouprint==YES) printf("final: az = %le\n",*azmax);
}

extern void beam_oneturn(struct Beam *bunch, struct Lattice *latt, char *txtfile)
{
	int i;
	
	for(i = 0; i < bunch->npart ; i++) {
		part_oneturn(&(bunch->part[i]), latt, txtfile);
	}
}

extern void beam_cross_latt(struct Beam *bunch, struct Lattice *latt, char *txtfile)
{
	int i;
	
	for(i = 0; i < bunch->npart ; i++) {
		part_cross_latt(&(bunch->part[i]), latt, txtfile);
	}
}

extern void part_oneturn(struct Particle *part, struct Lattice *latt, char *txtfile)
{
	int i;
	struct Lattice tempo_latt;

	part->s = 0; //part->s reset to 0 at the beginnig of each new Lattice
	
	//copy a Lattice, with cell[] stored at a diffent memory address, 
	//!!!and DO NOT FORGET to free the allocated memory before exiting part_oneturn!!!
	copy_latt(&tempo_latt, latt);
	if(latt->periodicity == 0) errorstop("!!!ERROR in part_oneturn: cell->periodicity = 0!!!");
	
	if(part->hat == 3) part->t = 0; //force t = 0 for rfscheme.cavity[].dt0 initialization turn
	
	for(i = 0; i < tempo_latt.periodicity; i++) {
		tempo_latt.cellnumber = i;
		//printf("cell %i\n",i);
		part_cross_latt(part, &tempo_latt,txtfile);		
		
		//move framework
		latt_move_exit_to_entrance(&tempo_latt);
	}
	
	//free allocated memory
	free(tempo_latt.cell);
}

extern void part_cross_latt(struct Particle *part, struct Lattice *latt, char *txtfile)
{
	int i;

	if(latt->nbcell <= 0) errorstop("!!!ERROR in part_cross_latt: latt->nbcell <= 0!!!");
	
	for(i = 0; i < latt->nbcell; i++) {
		if(strcmp(latt->cell[i].keyword, "rf-thingap") == 0) {	//rf cavity
			if(part->hat > 0) {
				if (latt->cellnumber < 0) {
						printf("latt->cellnumber = %i\n", latt->cellnumber);
						errorstop("!!!ERROR in part_cross_latt: cannot call part_cross_thingap() since latt->cellnumber < 0!!!\n");
					}
				part_cross_thingap(part, &(latt->cell[i].cav[latt->cellnumber]));
			}
		}
		else if(strcmp(latt->cell[i].keyword, "thin-hdipole") == 0) {
			if(part->hat > 0) part_thindipole(part,  &(latt->cell[i])); //localized dipole kick
		}
		else if(strcmp(latt->cell[i].keyword, "collimator") == 0) {
			if(part->status == ALIVE) part_collimator(part,  &(latt->cell[i])); //localized collimator
		}
		else if(strcmp(latt->cell[i].keyword, "float-faraday") == 0) {
			if(part->status == ALIVE) part_float_faraday(part,  &(latt->cell[i])); //localized collimator
		}
		else {	//in any other case call part_cross_cell
			//printf("\n\n\ncell %i nb of poles %i\n\n\n", i, latt->cell[i].nbcomp);
			if(part_cross_cell(part, &(latt->cell[i]),txtfile) == CUP) break;	//and if part is stop on a cup, break.
		}
	}
	latt->cellnumber = -1; //set cellnumber to a < 0 value, to be sure that a value will be given explicitely later. (other wise errorstop procedure is called, see above)
}

extern int part_cross_cell(struct Particle *part, struct Cell *cell, char *txtfile)
{
	int awd_answer;

	//debug tests
	if(part->status != ALIVE) return NO;
	if(cell->stepsize <= 0) errorstop("!!!ERROR in part_cross_cell: stepsize <= 0!!!");

	//change framework and do some checks
	getready_part(part, cell);
	
	//track untill reach a pickup, a cup or the cell boundary 
	awd_answer = part_callrk(part,cell,txtfile);
	
	if(awd_answer == PICKUP) { //reached a pickup, do not return yet, but swich off the pickup and call part_callrk again.
		if(part->status == ALIVE) {
			cell->instrutype = NO;
			awd_answer = part_callrk(part,cell,txtfile); //call again part_callrk
			cell->instrutype = PICKUP; //restore initial status
		}
	}
	if(awd_answer == YES) { //reached boundary, return
		if(part->status == ALIVE) update_part(awd_answer, part, cell);
		return awd_answer;
	}
	if(awd_answer == CUP) { //reached a cup, return
		//print_part_para("part before update", part);
		if(part->status == ALIVE) update_part(awd_answer, part, cell);
		//print_part_para("part after update", part);
		return awd_answer;
	}

	//if(part->status == ALIVE) errorstop("!!!ERROR in part_cross_cell: End of function reached but no boundary reached, while part->status is still =  ALIVE!!!");
	if(part->status == ALIVE) printf("!!!WARNING in part_cross_cell: End of function reached but no boundary reached, while part->status is still =  ALIVE!!!");
	return awd_answer;
}

extern int part_callrk(struct Particle *part, struct Cell *cell, char *txtfile)
{
	int awd_answer = NO;

	//track until a pickup of the cell boundary, choosing between the different methods to get the B field
	if (strcmp(cell->keyword, "drift") == 0)			rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_drift);
	else if (strcmp(cell->keyword, "ffag-r-lin") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGradial_lin);
	//else if (strcmp(cell->keyword, "ffag-r-enge") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGradial_2ndorder);
	else if (strcmp(cell->keyword, "ffag-r-enge") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGradial_purescale);
	else if (strcmp(cell->keyword, "ffag-r-he") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGradial_he);
	else if (strcmp(cell->keyword, "ffag-s-he") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGstraight_he);
	else if (strcmp(cell->keyword, "ffag-s-lin") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGstraight_lin);
	else if (strcmp(cell->keyword, "ffag-s-enge") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGstraight_purescale);
	else if (strcmp(cell->keyword, "ffag-s-enge-add") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGstraight_add);
	else if (strcmp(cell->keyword, "field-cylmap") == 0) rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, get_bfield_cylindrical_map);
	else if (strcmp(cell->keyword, "field-cartmap") == 0) rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, get_bfield_cartesian_map);
	else if (strcmp(cell->keyword, "ffag-tilt-he") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGtilt_he);
	else if (strcmp(cell->keyword, "ffag-tilt-lin") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGtilt_lin);
	//else if (strcmp(cell->keyword, "ffag-tilt-enge") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGtilt_enge);
	else if (strcmp(cell->keyword, "ffag-spi-he") == 0)		rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGspiral_he);
	else if (strcmp(cell->keyword, "ffag-spi-lin") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGspiral_lin);
	else if (strcmp(cell->keyword, "ffag-spi-enge") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGspiral_enge);
	else if (strcmp(cell->keyword, "ffag-spi-fullenge") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGspiral_fullenge);
	else if (strcmp(cell->keyword, "ffag-sti-he") == 0)		rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGstraighttilt_he);
	else if (strcmp(cell->keyword, "ffag-sti-lin") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGstraighttilt_lin);
	//else if (strcmp(cell->keyword, "ffag-sti-enge") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_FFAGstraighttilt_purescale);
	else if (strcmp(cell->keyword, "realbend") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_realbend);
	else if (strcmp(cell->keyword, "quad-he") == 0)		rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_quad_he);
	else if (strcmp(cell->keyword, "quad-lin") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_quad_lin);
	else if (strcmp(cell->keyword, "quad-enge") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_quad_enge);
	//else if (strcmp(cell->keyword, "sext-he") == 0)		errorstop("sext-he not implemented yet\n");
	//else if (strcmp(cell->keyword, "sext-lin") == 0)	errorstop("sext-lin not implemented yet\n");
	//else if (strcmp(cell->keyword, "sext-enge") == 0)	errorstop("sext-enge not implemented yet\n");
	//else if (strcmp(cell->keyword, "oct-he") == 0)		errorstop("sext-he not implemented yet\n");
	//else if (strcmp(cell->keyword, "oct-lin") == 0)	errorstop("sext-lin not implemented yet\n");
	//else if (strcmp(cell->keyword, "oct-enge") == 0)	errorstop("sext-enge not implemented yet\n");
	//else if (strcmp(cell->keyword, "purebend-he") == 0)		errorstop("purebend-he not implemented yet\n");
	else if (strcmp(cell->keyword, "purebend-lin") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_sector_purebend_lin);
	else if (strcmp(cell->keyword, "purebend-enge") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_sector_purebend);

	else if (strcmp(cell->keyword, "vffa-rect-he") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_rect_he);
	else if (strcmp(cell->keyword, "vffa-rect-lin") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_rect_lin);
	else if (strcmp(cell->keyword, "vffa-rect-enge") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_rect_enge);
	else if (strcmp(cell->keyword, "vffa-rect-atan") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_rect_atan);
	else if (strcmp(cell->keyword, "vffa-rect-enge-add") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_rect_enge_add);
	else if (strcmp(cell->keyword, "vffa-rect-enge-separ") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_rect_enge_separate_func);
	else if (strcmp(cell->keyword, "vffa-rect-enge-bx") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_rect_bx_add);
	
	else if (strcmp(cell->keyword, "vffa-rect-str-he") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_rect_he_str);
	else if (strcmp(cell->keyword, "vffa-rect-str-lin") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_rect_lin_str);
	else if (strcmp(cell->keyword, "vffa-rect-str-enge") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_rect_str_enge);
	else if (strcmp(cell->keyword, "vffa-rect-atan-str") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_rect_str_atan);
	
	//else if (strcmp(cell->keyword, "vffa-sect-he") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_sect_he);
	else if (strcmp(cell->keyword, "vffa-sect-lin") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_sect_lin);
	else if (strcmp(cell->keyword, "vffa-sect-enge") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_field_comp_VFFA_sect_enge);
	//else if (strcmp(cell->keyword, "multiple-poles") == 0)	rkdrive(part, &awd_answer, cell->stepsize,derivs,cell,txtfile, add_multiple_poles);
	
	else{
		printf("\n!!!ERROR 2 in part_callrk: unexpected keyword : \"%s\"!!!\n", cell->keyword);
		errorstop("unexpected Cell's keywork");
	}
	return awd_answer;
}

extern void part_cross_thingap(struct Particle *part, struct Cavity *cav)
{
	double e_tot_part, wt, v_rf, brho_before;//, u_hp_before;
	
	wt = TWOPI*(part->t - cav->t0)*cav->frf;
	wt	= wt - floor(wt/(TWOPI))*TWOPI;
	//printf("wt=%lf\n", wt);
	//if(part->hat == 3) {	//if part->hat == 3, it means that we are initializing cav->t0... 
	//	cav->t0 = part->t;
	//	wt= 0;
	//}
	
	v_rf= cav->v0*sin(wt + cav->phi0);
	//printf("wt+phi0=%lf [deg]\n", (wt+cav->phi0)*180./(PI));
	//printf("vrf=%le [V], phi = %lf\n", v_rf, (wt+cav->phi0)*180./PI);
	
	e_tot_part	= sqrt(pow(part->m0*CLIGHT2,2) + (part->brho*part->q*part->brho*part->q)*CLIGHT2) + v_rf*part->q;
	printf("etotpart=%le, delta e = %le\n", e_tot_part, v_rf*part->q);
	
	if(e_tot_part - part->m0*CLIGHT2 <= 0) { //kill particle if stoped or strats going bakward
		printf("\n!Warning in part_cross_thingap: particle is going backward, particle marked as lost (LOST_BACK)!\n");
		part->status = LOST_BACK;
		return;
	}
	
	//u_hp_before = sqrt(part->ux*part->ux + part->uy*part->uy);
	brho_before = part->brho;
	part->brho = sqrt(e_tot_part*e_tot_part/(CLIGHT2) - part->m0*CLIGHT2*part->m0)/part->q;
	
	//assuming that Px and Pz are unchanged:
	
	//part->uy = brho_before/part->brho*sqrt(pow(part->brho/brho_before,2) - part->ux*part->ux - part->uz*part->uz);
	part->ux = part->ux*brho_before/part->brho;
	part->uz = part->uz*brho_before/part->brho;
	part->uy = sqrt(1. - part->ux*part->ux - part->uz*part->uz);
	
	if(isnan(part->uy)) {
		printf("\n!Warning in part_cross_thingap: error has been done computing rf kick, particle marked as (LOST_ERR_RF),\nv_rf = %le, part->ux =%le, brho_before = %le, part->brho =%le\n", v_rf, part->ux, brho_before, part->brho);
		part->status = LOST_ERR_RF;
	}
	
}

extern void part_thindipole(struct Particle *part, struct Cell *cell)
{
	double xprime, uhp;
	
	uhp = sqrt(part->ux*part->ux + part->uy*part->uy);
	xprime = atan(part->ux/part->uy);
	xprime += cell->mpara[0][0]/part->brho; // kick of an angle BL*q/p
	part->ux = uhp*sin(xprime);
	part->uy = uhp*cos(xprime);
}

extern void part_collimator(struct Particle *part, struct Cell *cell)
{
	double r;
	
	if(cell->boun.thmax != 0) r = sqrt(part->x*part->x + part->y*part->y);
	else r = part->x;
	
	
	if(r < cell->collim.rmin) {
			printf("part_collimator: particle hits Horizontal aperture,  r = %lf < rmin = %lf [m]\n", r, cell->collim.rmin);
			part->status = LOST_H_MIN;
			return;
	}
	if(r > cell->collim.rmax) {
			printf("part_collimator: particle hits Horizontal aperture,  r = %lf > rmax = %lf [m]\n", r, cell->collim.rmax);
			part->status = LOST_H_MAX;
			return;
	}
	if(part->z < cell->collim.zmin) {
			printf("part_collimator: particle hits Vertical aperture, z = %lf < zmin = %lf [m]\n", part->z, cell->collim.zmin);
			part->status = LOST_V_MIN;
			return;
	}
	if(part->z > cell->collim.zmax) {
			printf("part_collimator: particle hits Vertical aperture, z = %lf > zmax = %lf [m]\n", part->z, cell->collim.zmax);
			part->status = LOST_V_MAX;
			return;
	}
}

// careful with vffa?
extern void part_float_faraday(struct Particle *part, struct Cell *cell)
{
	double r;
	
	if(cell->boun.thmax != 0) r = sqrt(part->x*part->x + part->y*part->y);
	else r = part->x;
	
	
	if(r > cell->collim.rmin && r < cell->collim.rmax) {
		printf("part_float_faraday: particle hits faraday cup, rmin = %lf < r = %lf < rmax = %lf [m]\n", cell->collim.rmin, r, cell->collim.rmax);
		part->status = LOST_FARADAY;
		return;
	}
}

extern void derivs(struct Particle *part, double dfds[], double bprobe[], struct Cell *cell,
	void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	double gamma;
	
	//Note: bprobe[0] = bx, bprobe[1] = by, and bprobe[2] = bz;
	part->status = get_bfield(part->x,part->y, part->z, &bprobe[0], &bprobe[1], &bprobe[2], cell, add_contribution_comp);
	//if(part->status!=ALIVE) printf("part status: %i\n", part->status);
	//printf("part_status=%i\n",part->status);
	//printf("part:x=%lf,y=%lf,z=%lf,ux=%lf,uy=%lf,uz=%lf\n",part->x,part->y, part->z, part->ux,part->uy, part->uz);
	//printf("b:%lf,%lf,%lf\n", bprobe[0],bprobe[1],bprobe[2]);
	
	dfds[0] = part->ux; //dx/ds
	dfds[1] = part->uy; //dy/ds
	dfds[2] = part->uz; //dz/ds
	dfds[3] = (part->uy*bprobe[2] - part->uz*bprobe[1])/part->brho; //dux/ds
	dfds[4] = (part->uz*bprobe[0] - part->ux*bprobe[2])/part->brho; //duy/ds
	dfds[5] = (part->ux*bprobe[1] - part->uy*bprobe[0])/part->brho; //duz/ds
	gamma = sqrt(1. + pow(part->q*part->brho/(part->m0*CLIGHT),2)); //reltivistic gamma
	dfds[6] = gamma/(CLIGHT*sqrt(gamma*gamma - 1.)); //dt/ds (= 1/(beta*CLIGHT))
	//printf("%lf,%lf,%lf,%lf,%lf,%lf,%lf\n", dfds[0], dfds[1],dfds[2],dfds[3],dfds[4],dfds[5],dfds[6]);
}

