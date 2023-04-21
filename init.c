/*
 *  init.c
 *  ringdesign
 *
 *  Copyright 2009 Kyoto University. All rights reserved.
 *
 */

#include "init.h"

// ************************************************************************************ //
//									Beams and particles									//
// ************************************************************************************ //
extern void gene_singlepart(struct Particle *part, double m0, double q, double x, double y, double z,
	double px, double py, double pz, double t0, int doyouprintf)
{
	double p = sqrt(px*px + py*py + pz*pz);
	
	if(doyouprintf == YES) printf("   p_part = %le [USI] = %le [eV/c], Brho = %.7f[T.m]\n", p, p*CLIGHT/(UNITCHARGE), p/q);
	if(p == 0) errorstop("!!!ERROR in gene_singlepart: total momentum of generated particle must not be 0");
	
	part->status	= ALIVE;
	//part->hat	= 3; //maximum angle search.
	//part->hat	= -2; //no acceleration, no output.
	//part->hat	= -1; //no acceleration, output.
	part->hat	= 1; //acceleration, output.
	part->m0	= m0;
	part->q		= q;
	part->x		= x;
	part->y		= y;
	part->z		= z;
	part->ux	= px/p;
	part->uy	= py/p;
	part->uz	= pz/p;
	part->brho	= p/(part->q);
	part->t		= t0;
	part->s		= 0;
	
	part->fwk.xc = 0;	// particle coordinates are
	part->fwk.yc = 0;	// initialy expressed in
	part->fwk.ae = 0;	// the global framework.
}

extern void gene_ellibunch(struct Beam *beam, double dt, double dp_ov_p, double dx, double dxprime, double dz, double dzprime, int nsteps[])
//here: d"A" = half variation of "A"
//it is assumed here that beam->part[0] has already been initialized, and this will be taken as the center of the bunch of particle that is going to be loaded in beam.
{
	int i1, i2, i3, i4, i5, i6, n;
	double t_ctr, p_ctr, x_ctr, xprime_ctr, z_ctr, zprime_ctr,
		t, p, dpp, x, xprime, y, z, zprime, px, py, pz,
		tl, tx, tz, al, ax, az, bl, bx, bz;
	struct Particle central_part;
	
	//test errors and warnings
	if(nsteps[1] <= 0) errorstop("!!!ERROR in gene_ellibunch:\n nsteps[1] <= 0");
	if(nsteps[2] <= 0) errorstop("!!!ERROR in gene_ellibunch:\n nsteps[2] <= 0");
	if(nsteps[3] <= 0) errorstop("!!!ERROR in gene_ellibunch:\n nsteps[3] <= 0");
	if(nsteps[4] <= 0) errorstop("!!!ERROR in gene_ellibunch:\n nsteps[4] <= 0");
	if(nsteps[5] <= 0) errorstop("!!!ERROR in gene_ellibunch:\n nsteps[5] <= 0");
	if(nsteps[6] <= 0) errorstop("!!!ERROR in gene_ellibunch:\n nsteps[6] <= 0");
	if(beam->part[0].uy == 0) errorstop("!!!ERROR in gene_ellibunch: beam->part[0].uy must not be 0!!!\n");
	
	//save beam->part[0] in central part and free the memory already alocated for beam
	central_part = beam->part[0];
	free_beam(beam);
	alloc_beam(beam, nsteps[1]*nsteps[2]*nsteps[3]*nsteps[4]*nsteps[5]*nsteps[6]+1);
	beam->part[0] = central_part;
	
	//compute central particle clue parameters (total momentum, x, xprime, z, zprime)
	t_ctr = beam->part[0].t;
	p_ctr = fabs(beam->part[0].brho*beam->part[0].q);
	x_ctr = beam->part[0].x;
	xprime_ctr = atan(beam->part[0].ux/beam->part[0].uy);
	z_ctr = beam->part[0].z;
	zprime_ctr = atan(beam->part[0].uz/sqrt(pow(beam->part[0].ux, 2) + pow(beam->part[0].uy, 2)));
	
	//generate a 6D beam, along ellipses centered on central particle
	for(i1 = 0; i1 < nsteps[1] ; i1++) {
		for(i2 = 0; i2 < nsteps[2] ; i2++) {
			for(i3 = 0; i3 < nsteps[3] ; i3++) {
				for(i4 = 0; i4 < nsteps[4] ; i4++) {
					for(i5 = 0; i5 < nsteps[5]; i5++) {
						for(i6 = 0; i6 < nsteps[6]; i6++) {
							
							tl = TWOPI*i1/(nsteps[1]);
							al = (i2 + 1.)/(nsteps[2])*dt;
							bl = (i2 + 1.)/(nsteps[2])*dp_ov_p;
							
							tx = TWOPI*i3/(nsteps[3]);
							ax = (i4 + 1.)/(nsteps[4])*dx;
							bx = (i4 + 1.)/(nsteps[4])*dxprime;
					
							tz = TWOPI*i5/(nsteps[5]);
							az = (i6 + 1.)/(nsteps[6])*dz;
							bz = (i6 + 1.)/(nsteps[6])*dzprime;
							
							if (nsteps[1] > 1) {
								t = t_ctr + al*cos(tl);
								dpp = bl*sin(tl);
							}
							else {
								t = t_ctr;
								dpp = 0.;
							}
							p = p_ctr*(1. + dpp);
							
							if (nsteps[3] > 1) {
								x = x_ctr + ax*cos(tx);
								xprime	= xprime_ctr + bx*sin(tx);
							}
							else {
								x = x_ctr;
								xprime = xprime_ctr;
							}
							
							if (nsteps[5] > 1) {
								z		= z_ctr + az*cos(tz);
								zprime	= zprime_ctr + bz*sin(tz);
							}
							else {
								z = z_ctr;
								zprime = zprime_ctr;
							}
							
							y = 0;
							px = p*sin(xprime)*cos(zprime);
							py = p*cos(xprime)*cos(zprime);
							pz = p*sin(zprime);
							
							n = 1 + i1*nsteps[6]*nsteps[5]*nsteps[4]*nsteps[3]*nsteps[2] + i2*nsteps[6]*nsteps[5]*nsteps[4]*nsteps[3] + i3*nsteps[6]*nsteps[5]*nsteps[4] + i4*nsteps[6]*nsteps[5] + i5*nsteps[6] + i6;
							
							gene_singlepart(&(beam->part[n]),beam->part[0].m0,beam->part[0].q, x, y, z, px, py, pz, t, NO);
						}
					}
				}
			}
		}
	}
}

extern void gene_rectbunch(struct Beam *beam, double dt, double dp_ov_p, double dx, double dxprime, double dz, double dzprime, int nsteps[])
//here: d"A" = half variation of "A"
//it is assumed here that beam->part[0] has already been initialized, and this will be taken as the center of the bunch of particle that is going to be loaded in beam.
{
	int i1, i2, i3, i4, i5, i6, n, flag;
	double t_ctr, p_ctr, x_ctr, xprime_ctr, z_ctr, zprime_ctr,
	t, p, dpp, x, xprime, y, z, zprime, px, py, pz;
	struct Particle central_part;
	
	//test errors and warnings
	if(nsteps[1] <= 0) errorstop("!!!ERROR in gene_rectbunch:\n nsteps[1] <= 0");
	if(nsteps[2] <= 0) errorstop("!!!ERROR in gene_rectbunch:\n nsteps[2] <= 0");
	if(nsteps[3] <= 0) errorstop("!!!ERROR in gene_rectbunch:\n nsteps[3] <= 0");
	if(nsteps[4] <= 0) errorstop("!!!ERROR in gene_rectbunch:\n nsteps[4] <= 0");
	if(nsteps[5] <= 0) errorstop("!!!ERROR in gene_rectbunch:\n nsteps[5] <= 0");
	if(nsteps[6] <= 0) errorstop("!!!ERROR in gene_rectbunch:\n nsteps[6] <= 0");
	if(beam->part[0].uy == 0) errorstop("!!!ERROR in gene_rectbunch: beam->part[0].uy must not be 0!!!\n");
	
	//save beam->part[0] in central part and free the memory already alocated for beam
	central_part = beam->part[0];
	free_beam(beam);
	alloc_beam(beam, nsteps[1]*nsteps[2]*nsteps[3]*nsteps[4]*nsteps[5]*nsteps[6]/2 + 2);
	beam->part[0] = central_part;
	
	//compute central particle clue parameters (total momentum, x, xprime, z, zprime)
	t_ctr = beam->part[0].t;
	p_ctr = fabs(beam->part[0].brho*beam->part[0].q);
	x_ctr = beam->part[0].x;
	xprime_ctr = atan(beam->part[0].ux/beam->part[0].uy);
	z_ctr = beam->part[0].z;
	zprime_ctr = atan(beam->part[0].uz/sqrt(pow(beam->part[0].ux, 2) + pow(beam->part[0].uy, 2)));
	
	//generate a 6D beam, along grids centered on central particle
	n = 1;
	for(i1 = 0; i1 < nsteps[1] ; i1++) {
		for(i2 = 0; i2 < nsteps[2] ; i2++) {
			for(i3 = 0; i3 < nsteps[3] ; i3++) {
				for(i4 = 0; i4 < nsteps[4] ; i4++) {
					for(i5 = 0; i5 < nsteps[5]; i5++) {
						for(i6 = 0; i6 < nsteps[6]; i6++) {
							
							if(nsteps[1] == 1) t = t_ctr;
							else t = t_ctr + dt*(2.*i1/(nsteps[1]-1) - 1);
							
							if(nsteps[2] == 1) dpp = 0;
							else dpp = dp_ov_p*(2.*i2/(nsteps[2]-1) - 1);
							
							if(nsteps[3] == 1) x = x_ctr;
							else x = x_ctr + dx*(2.*i3/(nsteps[3]-1) - 1);
							
							if(nsteps[4] == 1) xprime = xprime_ctr;
							else xprime = xprime_ctr + dxprime*(2.*i4/(nsteps[4]-1) - 1);
							//else xprime = xprime_ctr + dxprime*(1-(x - x_ctr)/(dx-x_ctr));
							
							if(nsteps[5] == 1) z = z_ctr;
							else z = z_ctr + dz*(2.*i5/(nsteps[5]-1) - 1);
							
							if(nsteps[6] == 1) zprime = zprime_ctr;
							else zprime = zprime_ctr + dzprime*(2.*i6/(nsteps[6]-1) - 1);
							
							p = p_ctr*(1. + dpp);
							y = 0;
							px = p*sin(xprime)*cos(zprime);
							py = p*cos(xprime)*cos(zprime);
							pz = p*sin(zprime);
							//printf("x-x_ctr = %lf, xprime-xprime_ctr = %lf\n", x-x_ctr, xprime-xprime_ctr);
							//printf("(xprime-xprime_ctr) - dxprime*(1 + (x-x_ctr)/dx) = %le, (xprime-xprime_ctr) - dxprime*(-1 + (x-x_ctr)/dx) = %le, dxprime*(1 + (x-x_ctr)/dx) = %le\n", (xprime-xprime_ctr) - dxprime*(1 + (x-x_ctr)/dx), (xprime-xprime_ctr) - dxprime*(-1 + (x-x_ctr)/dx), dxprime*(1 + (x-x_ctr)/dx));
							//printf("(xprime-xprime_ctr) - dxprime*(1 + (x-x_ctr)/dx) = %le, (xprime-xprime_ctr) - dxprime*(-1 - (x-x_ctr)/dx) = %le\n", (xprime-xprime_ctr) - dxprime*(1 + (x-x_ctr)/dx), (xprime-xprime_ctr) - dxprime*(-1 - (x-x_ctr)/dx));
							if((x - x_ctr) > 0) {
								if( ( (xprime-xprime_ctr) - (dxprime*(1 - (x-x_ctr)/dx)) <= TINYLENGTH ) && (( (xprime-xprime_ctr) - (dxprime*(-1 + (x-x_ctr)/dx))) >= -1*TINYLENGTH ) ) flag = 1;
								else flag = 0;
							}
							else if((x - x_ctr) < 0) {
								if( ( (xprime-xprime_ctr) - (dxprime*(1 + (x-x_ctr)/dx)) <= TINYLENGTH ) && ( (xprime-xprime_ctr) - (dxprime*(-1 - (x-x_ctr)/dx)) >= -1*TINYLENGTH ) ) flag = 1;
								else flag = 0;
							}
							if((x - x_ctr) == 0) {
								flag = 1;
							}
							//n = 1 + i1*nsteps[6]*nsteps[5]*nsteps[4]*nsteps[3]*nsteps[2] + i2*nsteps[6]*nsteps[5]*nsteps[4]*nsteps[3] + i3*nsteps[6]*nsteps[5]*nsteps[4] + i4*nsteps[6]*nsteps[5] + i5*nsteps[6] + i6;
							
							//printf("flag = %i\n\n", flag);
							if(flag == 1) {
								gene_singlepart(&(beam->part[n]),beam->part[0].m0,beam->part[0].q, x, y, z, px, py, pz, t, NO);
								n++;
							}
						}
					}
				}
			}
		}
	}
}

extern void gene_elliran_h(struct Beam *beam, double dt, double dp_ov_p, double emit_x, double twiss_betax, double twiss_alphax, double emit_z, double twiss_betaz, double twiss_alphaz, int nparts, long seed)
//here: d"A" = half variation of "A"
//it is assumed here that beam->part[0] has already been initialized, and this will be taken as the center of the bunch of particle that is going to be loaded in beam.
//emit_x and _z are geometrical emittances given in pi.m.rad
{
	int i;
	long idum;
	double t_ctr, p_ctr, x_ctr, xprime_ctr, z_ctr, zprime_ctr,
		t=0., p, dpp, x, xprime, y, z, zprime, px, py, pz,
		emit_max_x, emit_max_z, xn, xpn, zn, zpn;
	struct Particle central_part;
	
	if(debug == YES) printf("gene_elliran_h: nparts = %i\n", nparts);
	
	//test errors and warnings
	if(nparts <= 0) errorstop("!!!ERROR in gene_elliran_h:\n nparts <= 0");
	if(beam->part[0].uy == 0) errorstop("!!!ERROR in gene_elliran_h: beam->part[0].uy = 0!!!\n");
	if(twiss_betax <= 0) errorstop("!!!ERROR in gene_elliran_h: twiss_betax <= 0");
	if(twiss_betaz <= 0) errorstop("!!!ERROR in gene_elliran_h: twiss_betaz <= 0");
	
	//compute corresponding emit_max
	//emit_max_x = dx*dx/twiss_betax;
	//emit_max_z = dz*dz/twiss_betaz;
	emit_max_x = emit_x;
	emit_max_z = emit_z;
	
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
	
	//generate a 6D beam, inside ellipses centered on central particle
	for(i = 1; i < nparts+1 ; i++) {
		
		//elliptical distribution in longitudinal plane
		dpp = 2*(ran1(&idum)-0.5);
		//for(;;) {
		//	t = 2*(ran1(&idum)-0.5);
		//	dpp = 2*(ran1(&idum)-0.5);
		//	if((t*t + dpp*dpp) <= 1) break; //if the particle randomly generated is not inside the 2D longitudinal ellipse, ask for another one
		//}
		t *= dt;
		dpp *= dp_ov_p;
		
		//4D elliptical distribution in transverse planes. 
		//!!!We use NORMALIZED transverse coordinates (xn = x/sqrt(twiss_beta*emit_max_x), xpn = ((twiss_alpha*x + twiss_beta*xp)/sqrt(twiss_beta*emit_max_x), same thing for zn and zpn)
		for(;;) {
			xn = 2.*(ran1(&idum)-0.5);
			xpn = 2.*(ran1(&idum)-0.5);
			zn = 2.*(ran1(&idum)-0.5);
			zpn = 2.*(ran1(&idum)-0.5);
			if(sqrt(xn*xn + xpn*xpn + zn*zn + zpn*zpn) <= 1) break; //if the particle randomly generated is not inside 4D transverse ellipse, ask for another one
		}
		
		//transfor normalized transverse coordinates (xn,xpn,zn,zpn) into real transverse coordinates (x,xp,z,zp)
		x = sqrt(twiss_betax*emit_max_x)*xn;
		xprime = (sqrt(twiss_betax*emit_max_x)*xpn - twiss_alphax*x)/twiss_betax;
		z = sqrt(twiss_betaz*emit_max_z)*zn;
		zprime = (sqrt(twiss_betaz*emit_max_z)*zpn - twiss_alphaz*z)/twiss_betaz;
		
		t += t_ctr;
		p = p_ctr*(1. + dpp);
		
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

extern void gene_elligasdev_6d(struct Beam *beam, double dt, double dp_ov_p, double emit_x, double twiss_betax, double twiss_alphax, double emit_z, double twiss_betaz, double twiss_alphaz, int nparts, long seed)
//here: d"A" = half variation of "A"
//it is assumed here that beam->part[0] has already been initialized, and this will be taken as the center of the bunch of particle that is going to be loaded in beam.
//emit_x and _z are geometrical emittances given in pi.m.rad
//emittance is rms one
{
	int i;
	long idum;
	double t_ctr, p_ctr, x_ctr, xprime_ctr, z_ctr, zprime_ctr,
	t, p, dpp, x, xprime, y, z, zprime, px, py, pz,
	xn, xpn, zn, zpn, sigmax, sigmaxp, sigmaz, sigmazp, factor_rms_tot;
	struct Particle central_part;
	
	if(debug == YES) printf("gene_elligasdev_6d: nparts = %i\n", nparts);
	
	//test errors and warnings
	if(nparts <= 0) errorstop("!!!ERROR in gene_elligasdev_6d:\n nparts <= 0");
	if(beam->part[0].uy == 0) errorstop("!!!ERROR in gene_elligasdev_6d: beam->part[0].uy = 0!!!\n");
	if(twiss_betax <= 0) errorstop("!!!ERROR in gene_elligasdev_6d: twiss_betax <= 0");
	if(twiss_betaz <= 0) errorstop("!!!ERROR in gene_elligasdev_6d: twiss_betaz <= 0");
	if(fabs(twiss_alphax) > TINYDIMLESS || fabs(twiss_alphaz) > TINYDIMLESS) errorstop("!!!ERROR in gene_elligasdev_6d: twiss_alpha != 0");
	if(fabs(dt) > TINYDIMLESS) errorstop("!!!ERROR in gene_elligasdev_6d: dt != 0, change the code");
	
	factor_rms_tot = 6.;
	
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
	
	sigmax = sqrt(emit_x*twiss_betax);
	sigmaxp = sqrt(emit_x*(1+twiss_alphax*twiss_alphax)/twiss_betax);
	sigmaz = sqrt(emit_z*twiss_betaz);
	sigmazp = sqrt(emit_z*(1+twiss_alphaz*twiss_alphaz)/twiss_betaz);
	
	//printf("sigma: x=%le, xp=%le, \nz=%le, zp=%le\n", sigmax, sigmaxp, sigmaz, sigmazp);
	
	//initialise idum (for the random generator)
	if(seed >=0) idum = -seed;
	else idum = seed;
	
	//generate a 6D beam, inside ellipses centered on central particle
	for(i = 1; i < nparts+1 ; i++) {
		//printf("part number %i\n", i);
		//elliptical distribution in longitudinal plane
		//for(;;) {
		//	t = 2*(ran1(&idum)-0.5);
		//	dpp = 2*(ran1(&idum)-0.5);
		//	if((t*t + dpp*dpp) <= 1) break; //if the particle randomly generated is not inside the 2D longitudinal ellipse, ask for another one
		//}
		//GAUSSIAN DISTRIBUTION IN LONGITUDINAL !!!!!!!!!!!!!!!!!!!! 
		for(;;) {
			t = gasdev(&idum)*dt;
			dpp = gasdev(&idum)*dp_ov_p;
			//printf("t = %le, dpp = %le, dpovp = %le\n", t, dpp, dp_ov_p);
			if(fabs(dpp) < factor_rms_tot*dp_ov_p) break;
			//if(fabs(t) < factor_rms_tot*dt && fabs(dpp) < factor_rms_tot*dp_ov_p) break;
		}
		//4D elliptical distribution in transverse planes. 
		//!!!We use NORMALIZED transverse coordinates (xn = x/sqrt(twiss_beta*emit_max_x), xpn = ((twiss_alpha*x + twiss_beta*xp)/sqrt(twiss_beta*emit_max_x), same thing for zn and zpn)
		/*for(;;) {
		 xn = 2.*(ran1(&idum)-0.5);
		 xpn = 2.*(ran1(&idum)-0.5);
		 zn = 2.*(ran1(&idum)-0.5);
		 zpn = 2.*(ran1(&idum)-0.5);
		 if(sqrt(xn*xn + xpn*xpn + zn*zn + zpn*zpn) <= 1) break; //if the particle randomly generated is not inside 4D transverse ellipse, ask for another one
		 }*/
		for(;;) {
			xn = gasdev(&idum)*sigmax;
			xpn = gasdev(&idum)*sigmaxp;
			if(fabs(xn) < factor_rms_tot*sigmax && fabs(xpn) < factor_rms_tot*sigmaxp) break;
		}
		for(;;) {
			zn = gasdev(&idum)*sigmaz;
			zpn = gasdev(&idum)*sigmazp;
			if(fabs(zn) < factor_rms_tot*sigmaz && fabs(zpn) < factor_rms_tot*sigmazp) break;
		}
		
		//transform normalized transverse coordinates (xn,xpn,zn,zpn) into real transverse coordinates (x,xp,z,zp)
		//x = sqrt(twiss_betax*emit_max_x)*xn;
		//xprime = (sqrt(twiss_betax*emit_max_x)*xpn - twiss_alphax*x)/twiss_betax;
		//z = sqrt(twiss_betaz*emit_max_z)*zn;
		//zprime = (sqrt(twiss_betaz*emit_max_z)*zpn - twiss_alphaz*z)/twiss_betaz;
		
		t += t_ctr;
		p = p_ctr*(1. + dpp);
		
		x = x_ctr + xn;
		xprime	= xprime_ctr + xpn;
		
		z = z_ctr + zn;
		zprime	= zprime_ctr + zpn;
		
	//	xprime -=twiss_alphax/twiss_betax*xn; //tilt the ellipse in the horizontal plane.
		
	//	zprime -=twiss_alphaz/twiss_betaz*zn; //tilt the ellipse in the horizontal plane.
		
		y = 0;
		px = p*sin(xprime)*cos(zprime);
		py = p*cos(xprime)*cos(zprime);
		pz = p*sin(zprime);
		
		//generate a particle at (x,y,z) with momentum (px,py,pz) and time t.
		gene_singlepart(&(beam->part[i]),beam->part[0].m0,beam->part[0].q, x, y, z, px, py, pz, t, NO);
	}// and it is a loop, it does it again nparts times...
}

extern void gene_elligasdev_4d(struct Beam *beam, double dt, double dp_ov_p, double emit_x, double twiss_betax, double twiss_alphax, double emit_z, double twiss_betaz, double twiss_alphaz, int nparts, long seed)
//here: d"A" = half variation of "A"
//it is assumed here that beam->part[0] has already been initialized, and this will be taken as the center of the bunch of particle that is going to be loaded in beam.
//emit_x and _z are geometrical emittances given in pi.m.rad
//emittance is rms one
{
	int i;
	long idum;
	double t_ctr, p_ctr, x_ctr, xprime_ctr, z_ctr, zprime_ctr,
	t=0., p, dpp, x, xprime, y, z, zprime, px, py, pz,
	xn, xpn, zn, zpn, sigmax, sigmaxp, sigmaz, sigmazp, factor_rms_tot;
	struct Particle central_part;
	
	if(debug == YES) printf("gene_elligasdev_4d: nparts = %i\n", nparts);
	
	//test errors and warnings
	if(nparts <= 0) errorstop("!!!ERROR in gene_elligasdev_4d:\n nparts <= 0");
	if(beam->part[0].uy == 0) errorstop("!!!ERROR in gene_elligasdev_4d: beam->part[0].uy = 0!!!\n");
	if(twiss_betax <= 0) errorstop("!!!ERROR in gene_elligasdev_4d: twiss_betax <= 0");
	if(twiss_betaz <= 0) errorstop("!!!ERROR in gene_elligasdev_4d: twiss_betaz <= 0");
	if(fabs(twiss_alphax) > TINYDIMLESS || fabs(twiss_alphaz) > TINYDIMLESS) errorstop("!!!ERROR in gene_elligasdev_4d: twiss_alpha != 0");
	if(fabs(dt) > TINYDIMLESS) errorstop("!!!ERROR in gene_elligasdev_4d: dt != 0, change the code");
	
	factor_rms_tot = 3.;
	
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
	
	sigmax = sqrt(emit_x*twiss_betax);
	sigmaxp = sqrt(emit_x*(1+twiss_alphax*twiss_alphax)/twiss_betax);
	sigmaz = sqrt(emit_z*twiss_betaz);
	sigmazp = sqrt(emit_z*(1+twiss_alphaz*twiss_alphaz)/twiss_betaz);
	
	//printf("sigma: x=%le, xp=%le, \nz=%le, zp=%le\n", sigmax, sigmaxp, sigmaz, sigmazp);
	
	//initialise idum (for the random generator)
	if(seed >=0) idum = -seed;
	else idum = seed;
	
	//generate a 6D beam, inside ellipses centered on central particle
	for(i = 1; i < nparts+1 ; i++) {
		//printf("part number %i\n", i);
		//elliptical distribution in longitudinal plane
		//for(;;) {
		//	t = 2*(ran1(&idum)-0.5);
		//	dpp = 2*(ran1(&idum)-0.5);
		//	if((t*t + dpp*dpp) <= 1) break; //if the particle randomly generated is not inside the 2D longitudinal ellipse, ask for another one
		//}
		
		//UNIFORM DISTRIBUTION IN LONGITUDINAL !!!!!!!!!!!!!!!!!!!!!!
		dpp = 2*(ran1(&idum)-0.5);
		t *= dt;
		dpp *= dp_ov_p;
		//printf("dpp = %lf\n", dpp);
		
		//4D elliptical distribution in transverse planes. 
		//!!!We use NORMALIZED transverse coordinates (xn = x/sqrt(twiss_beta*emit_max_x), xpn = ((twiss_alpha*x + twiss_beta*xp)/sqrt(twiss_beta*emit_max_x), same thing for zn and zpn)
		/*for(;;) {
		 xn = 2.*(ran1(&idum)-0.5);
		 xpn = 2.*(ran1(&idum)-0.5);
		 zn = 2.*(ran1(&idum)-0.5);
		 zpn = 2.*(ran1(&idum)-0.5);
		 if(sqrt(xn*xn + xpn*xpn + zn*zn + zpn*zpn) <= 1) break; //if the particle randomly generated is not inside 4D transverse ellipse, ask for another one
		 }*/
		for(;;) {
			xn = gasdev(&idum)*sigmax;
			xpn = gasdev(&idum)*sigmaxp;
			if(fabs(xn) < factor_rms_tot*sigmax && fabs(xpn) < factor_rms_tot*sigmaxp) break;
		}
		for(;;) {
			zn = gasdev(&idum)*sigmaz;
			zpn = gasdev(&idum)*sigmazp;
			if(fabs(zn) < factor_rms_tot*sigmaz && fabs(zpn) < factor_rms_tot*sigmazp) break;
		}
		
		//transform normalized transverse coordinates (xn,xpn,zn,zpn) into real transverse coordinates (x,xp,z,zp)
		//x = sqrt(twiss_betax*emit_max_x)*xn;
		//xprime = (sqrt(twiss_betax*emit_max_x)*xpn - twiss_alphax*x)/twiss_betax;
		//z = sqrt(twiss_betaz*emit_max_z)*zn;
		//zprime = (sqrt(twiss_betaz*emit_max_z)*zpn - twiss_alphaz*z)/twiss_betaz;
		
		t += t_ctr;
		p = p_ctr*(1. + dpp);
		
		x = x_ctr + xn;
		xprime	= xprime_ctr + xpn;
		
		z = z_ctr + zn;
		zprime	= zprime_ctr + zpn;
		
	//	xprime -=twiss_alphax/twiss_betax*xn; //tilt the ellipse in the horizontal plane.
		
	//	zprime -=twiss_alphaz/twiss_betaz*zn; //tilt the ellipse in the horizontal plane.
		
		y = 0;
		px = p*sin(xprime)*cos(zprime);
		py = p*cos(xprime)*cos(zprime);
		pz = p*sin(zprime);
		
		//generate a particle at (x,y,z) with momentum (px,py,pz) and time t.
		gene_singlepart(&(beam->part[i]),beam->part[0].m0,beam->part[0].q, x, y, z, px, py, pz, t, NO);
	}// and it is a loop, it does it again nparts times...
}

extern void load_beam(struct Beam *beam, char *inputfilename, struct Lattice *latt, int doyouprint)
{
	int i,j, maxop, accept_flag, nsteps[7], tot_nb_parts=0;
	double m0, q, ekin_ev, x, y, z, px, py, pz, t, ptot, atan_dxdy_deg, atan_dzdy_deg, 
		dt, dp_ov_p, dx, dxprime, dz, dzprime, eps_clo, stepx;
	char temp_str[MAX_CHARINLINE+2], keyword[MAX_CHARINLINE+2];
	struct Beam beam_temp;
	FILE *rfile = NULL;
	
	rfile = fopen(inputfilename, "r");
	if(rfile == NULL) errorstop ("!!!ERROR in load_beam: cannot open inputfile!!!");

	//read input file
	printf("\nBEAM INPUT: ");
	if(print_color==YES) COLOR("32");
	printf("%s ", inputfilename);
	if(print_color==YES) COLOR("0");
	printf("...\n\n");
	fgets(temp_str, MAX_CHARINLINE+2, rfile);
	newline(rfile);
	if(strlen(temp_str) > MAX_CHARINLINE) errorstop("!!!ERROR in load_beam: maximum number of char per line in input file exceeded!!!");
	if(debug == YES) printf("... header: \t\t%s", temp_str);
	
	//read MassCharge
	fscanf(rfile, "%lf  %lf", &m0, &q );
	newline(rfile);
	if(debug == YES) printf("... MassCharge: \t%le [MeV*c^-2]\t%le [q.C]\n", m0, q);
	m0 = m0*1.e6*UNITCHARGE/(CLIGHT2); //convert to [kg]
	q = q*UNITCHARGE; //convert to [C]
	
	//number of keywords to be read
	fscanf(rfile, "%i", &maxop);
	newline(rfile);
	if(debug == YES) printf("... number of keywords to be read: %i\n", maxop);
	
	//read up to maxop keywords
	for(i = 1; i <= maxop; i++) {
		
		
		//read the keyword (and check that we have not reached the end of the input file)
		if(fscanf(rfile, "%s", keyword) == EOF) {
			printf("\n!!!ERROR in load_beam: reaching end of file %s\n", inputfilename);
			errorstop("End of file reached unexpectedly, you may try to read more lines than the number of lines in your input file. Please check your input file.");
		}
		
		//generate singlepart
		if(strcmp(keyword, "singlepart") == 0) {
			tot_nb_parts++;
			//allocate memory in beam for a single particle
			if(tot_nb_parts==1) alloc_beam(beam, 1);
			else {
				copy_beam(&beam_temp, beam);
				free_beam(beam);
				alloc_beam(beam, tot_nb_parts);
				for(j=0;j<tot_nb_parts-1;j++) beam->part[j] = beam_temp.part[j];
				free_beam(&beam_temp);
			}
			
			if(debug == YES) printf("  --> single particle generation routine: \"%s\"\n", keyword);
			fscanf(rfile, "%le", &ekin_ev);
			newline(rfile);
			fscanf(rfile, "%lf %lf %lf", &x, &y, &z);
			newline(rfile);
			fscanf(rfile, "%lf %lf", &atan_dxdy_deg, &atan_dzdy_deg); //atan_dxds_deg: angle between the particle traj. and y axis [deg], atan_dzdy_deg: angle between the particle traj. and the horizontal plane [deg]
			newline(rfile);
			ptot = sqrt(pow(ekin_ev*UNITCHARGE/(CLIGHT) + m0*CLIGHT,2) - pow(m0*CLIGHT,2));
			px = ptot*sin(atan_dxdy_deg*PI/180.)*cos(atan_dzdy_deg*PI/180.);
			py = ptot*cos(atan_dxdy_deg*PI/180.)*cos(atan_dzdy_deg*PI/180.);
			pz = ptot*sin(atan_dzdy_deg*PI/180.);
			t = 0;
			if(debug == YES) printf("in load_beam, generate singlepart, py = %le\n", py);
			gene_singlepart(&(beam->part[tot_nb_parts-1]), m0, q, x, y, z, px, py, pz, t, doyouprint);
			//beam->part[0].hat = 2; //acceleration, no output
			beam->part[tot_nb_parts-1].hat = 1; //acceleration and output
		}
		
		//put central particle on the closed orbit (x, xprime = 0) (and exit process if closed orbit in not found)
		else if(strcmp(keyword, "putonco-x") == 0) {
			fscanf(rfile, "%le", &eps_clo);
			newline(rfile);
			if(doyouprint==YES) printf("\nNOW LOOKING FOR CLOSED ORBIT ...\n");
			accept_flag = find_closed_orbite_x(&(beam->part[0]), &(beam->part[0].x), eps_clo, latt, doyouprint);
			if(accept_flag != TRUE) errorstop("!!!ERROR in load_beam, closed orbit not found");
			if(doyouprint==YES) printf("\n");
		}
		
		//put central particle on the closed orbit  (x and xprime) (and exit process if closed orbit in not found)
		else if(strcmp(keyword, "putonco-xxp") == 0) {
			fscanf(rfile, "%le", &eps_clo);
			newline(rfile);
			if(doyouprint==YES) printf("\nNOW LOOKING FOR CLOSED ORBIT ...\n");
			accept_flag = find_closed_orbite_xxp(&(beam->part[0]), &(beam->part[0].x), &(beam->part[0].ux), &(beam->part[0].uy), eps_clo, latt, doyouprint);
			if(accept_flag != TRUE) errorstop("!!!ERROR in load_beam, closed orbite not found!!!");
			if(doyouprint==YES) printf("\n");
		}
	
		//put central particle on the closed orbit  (x and xprime, z and zprime) (and exit process if closed orbit in not found)
		else if(strcmp(keyword, "putonco-xxpzzp") == 0) {
			fscanf(rfile, "%le", &eps_clo);
			newline(rfile);
			if(doyouprint==YES) printf("\nNOW LOOKING FOR CLOSED ORBIT ...\n");
			accept_flag = find_closed_orbite_xxp_zzp(&(beam->part[0]), &(beam->part[0].x), &(beam->part[0].ux), &(beam->part[0].uy), &(beam->part[0].z), &(beam->part[0].uz), eps_clo, latt, doyouprint);
			if(accept_flag != TRUE) errorstop("!!!ERROR in load_beam, closed orbite not found!!!");
			if(doyouprint==YES) printf("\n");
		}
		
		else if(strcmp(keyword, "putonco-nelmin") == 0) {
			fscanf(rfile, "%le %le", &eps_clo, &stepx);
			newline(rfile);
			if(doyouprint==YES) printf("\nNOW LOOKING FOR CLOSED ORBIT ...\n");
			accept_flag = put_on_co_nelmin(latt, &(beam->part[0]), eps_clo, stepx, stepx, stepx, stepx, doyouprint);
			if(accept_flag != TRUE) errorstop("!!!ERROR in load_beam, closed orbite not found!!!");
			if(doyouprint==YES) printf("\n");
		}
		
		//generate ellibunch
		else if(strcmp(keyword, "ellibunch") == 0) {
			if(debug == YES) printf("  --> bunch generation routine: \"%s\"\n", keyword);
			fscanf(rfile, "%i %i %i %i %i %i", &(nsteps[1]), &(nsteps[2]), &(nsteps[3]), &(nsteps[4]), &(nsteps[5]), &(nsteps[6]));
			if(debug == YES) printf("\tbunch number of particles in the bunch: [%i*%i]*[%i*%i]*[%i*%i]\n", nsteps[1], nsteps[2], nsteps[3], nsteps[4], nsteps[5], nsteps[6]);
			newline(rfile);
			fscanf(rfile, "%lf %lf", &dt, &dp_ov_p);
			if(debug == YES) printf("\t1/2.Delta(t) = %lf [s]\t1/2.Delta(dp/p) = %lf [percent]\n", dt, dp_ov_p*100);
			newline(rfile);
			fscanf(rfile, "%lf %lf", &dx, &dxprime);
			if(debug == YES) printf("\t1/2.Delta(x) = %lf [mm]\t1/2.Delta(xprime) = %lf [mrad]\n", dx*1000, dxprime*1000);
			newline(rfile);
			fscanf(rfile, "%lf %lf", &dz, &dzprime);
			if(debug == YES) printf("\t1/2.Delta(z) = %lf [mm]\t1/2.Delta(zprime) = %lf [mrad]\n", dz*1000, dzprime*1000);
			newline(rfile);			
			//fill the beam with nsteps[1]*nsteps[2]*..*nsteps[6] + 1 particles. Memory allocation is done inside gene_ellibunch
			gene_ellibunch(beam, dt, dp_ov_p, dx, dxprime, dz, dzprime, nsteps);
		}
		
		//generate rectangular bunch
		else if(strcmp(keyword, "rectbunch") == 0) {
			if(debug == YES) printf("  --> bunch generation routine: \"%s\"\n", keyword);
			fscanf(rfile, "%i %i %i %i %i %i", &(nsteps[1]), &(nsteps[2]), &(nsteps[3]), &(nsteps[4]), &(nsteps[5]), &(nsteps[6]));
			if(debug == YES) printf("\tbunch number of particles in the bunch: [%i*%i]*[%i*%i]*[%i*%i]\n", nsteps[1], nsteps[2], nsteps[3], nsteps[4], nsteps[5], nsteps[6]);
			newline(rfile);
			fscanf(rfile, "%lf %lf", &dt, &dp_ov_p);
			if(debug == YES) printf("\t1/2.Delta(t) = %lf [s]\t1/2.Delta(dp/p) = %lf [percent]\n", dt, dp_ov_p*100);
			newline(rfile);
			fscanf(rfile, "%lf %lf", &dx, &dxprime);
			if(debug == YES) printf("\t1/2.Delta(x) = %lf [mm]\t1/2.Delta(xprime) = %lf [mrad]\n", dx*1000, dxprime*1000);
			newline(rfile);
			fscanf(rfile, "%lf %lf", &dz, &dzprime);
			if(debug == YES) printf("\t1/2.Delta(z) = %lf [mm]\t1/2.Delta(zprime) = %lf [mrad]\n", dz*1000, dzprime*1000);
			newline(rfile);
			//fill the beam with nsteps[1]*nsteps[2]*..*nsteps[6] + 1 particles. Memory allocation is done inside gene_rectbunch		
			gene_rectbunch(beam, dt, dp_ov_p, dx, dxprime, dz, dzprime, nsteps);
		}
		
		else {
			printf("unrecognized keyword = %s while reading %s\n", keyword, inputfilename);
			errorstop("!!!ERROR in load_beam: unreconized keyword");
		}
	}
		
	fclose(rfile);
}

extern void load_ao_beam(struct Beam *beam, char *inputfilename, double x_center_pipe)
{
	char name[500], shell_write[100];
	int i, n, id, evnum, trkid, parent, weight;
	double m0, q, x, y, z, px, py, pz, t;
	FILE *shell_read;
	FILE *rfile = NULL;
	
	printf("\nREADING INPUT AO BEAM: ");
	if(print_color==YES) COLOR("32");
	printf("%s \n", inputfilename);
	if(print_color==YES) COLOR("0");
	
	sprintf(shell_write,"wc -l %s", inputfilename);
	
	shell_read = popen(shell_write,"r");
	fscanf(shell_read, "%i %s", &n, name);
	pclose(shell_read);
	printf("nb_lines = %i in %s\n", n, name);
	
	rfile = fopen(inputfilename,"r");
	if(rfile == NULL) errorstop("cannot open file in load_ao_beam\n");
	//printf("n = %i\n\n\n", n);
	alloc_beam(beam, n);
	
	m0 = MUONM;
	q = +1*UNITCHARGE;
	
	for(i = 0; i < n;i++) {
		fscanf(rfile, "%le %le %le %le %le %le %le %i %i %i %i %i\n", &x, &z, &y, &px, &pz, &py, &t, &id, &evnum, &trkid, &parent, &weight);
		//printf("%le %le %le %le %le %le %le %i %i %i %i %i\n", x, z, y, px, pz, py, t, id, evnum, trkid, parent, weight);
		if(id == -13) {
			x *= 1e-3; //in [m]
			x += x_center_pipe; // + center of the pipe
			y *= 1e-3; //in [m]
			z *= 1e-3; //in [m]
			px *= 1.e6*UNITCHARGE/CLIGHT; //in SI
			py *= 1.e6*UNITCHARGE/CLIGHT; //in SI
			pz *= 1.e6*UNITCHARGE/CLIGHT; //in SI
			t *= 1.e-9; //in [s]
			
			gene_singlepart(&(beam->part[i]), m0, q, x, y, z, px, py, pz, t, NO);
		}
		else errorstop("not a muon\n!");
	}
	fclose(rfile);
	
}

// ************************************************************************************ //
//									lattices, cells										//
// ************************************************************************************ //
extern void load_lattice(struct Lattice *latt, char *inputfilename)
{
	int i, nb_idcells, cell_nb, nb_typecell;
	char temp_str[MAX_CHARINLINE+2];
	struct Framework fwk;
	FILE *rfile = NULL;
	if(debug == YES) printf("enter load_lattice\n");
	//printf("lattice file: %s\n", inputfilename);
	fwk.xc = 0;
	fwk.yc = 0;
	fwk.ae = 0;
	
	//read input file
	printf("\nLATTICE INPUT: ");
	if(print_color==YES) COLOR("32");
	printf("%s ", inputfilename);
	if(print_color==YES) COLOR("0");
	printf("...\n");
	if(debug == YES) printf("\n");
	rfile = fopen(inputfilename, "r");
	if(rfile == NULL) errorstop("!!!ERROR in load_lattice: cannot open inputfile!!!");
	
	fgets(temp_str, MAX_CHARINLINE+2, rfile);
	newline(rfile);
	if(strlen(temp_str) > MAX_CHARINLINE) errorstop("!!!ERROR maximum number of char per line in input file exceeded!!!");
	if(debug == YES) printf("... header: \t\t%s", temp_str);
	
	//lattice periodicity
	fscanf(rfile, "%i", &(latt->periodicity));
	newline(rfile);
	if(debug == YES) printf("... periodicity: \t%i\n", latt->periodicity);
	
	//number of Cells in the Lattice
	fscanf(rfile, "%i", &(latt->nbcell));
	newline(rfile);
	if(debug == YES) printf("... nbcell: \t\t%i\n", latt->nbcell);
	
	//allocate memory for latt->cell[]
	latt->cell = alloccell(latt->nbcell);
	
	//number of Cells in the Lattice
	fscanf(rfile, "%i", &(nb_typecell));
	newline(rfile);
	if(debug == YES) printf("... nb of cell type: \t\t%i\n", nb_typecell);
	
	//load every single Cell
	cell_nb = 0;
	for(i = 0; i < nb_typecell; i++) {
		if(debug == YES) printf("... load Cell type number %i:\n", i);
		if(debug == YES) printf("... load Cell[%i]:\n", cell_nb);
		//initialize deltar, nbcomp, instrutype, and doyou_err
		latt->cell[cell_nb].deltar = 0; //default value
		latt->cell[cell_nb].nbcomp = 1; //default value
		latt->cell[cell_nb].instrutype = NO; //default value
		latt->cell[cell_nb].doyou_err = NO; //default value
		
		//set cell.framework
		latt->cell[cell_nb].framework = fwk;
		if(debug == YES) printf("  --> framework:\n\txc = %lf[m] yc = %lf[m] ae = %lf[deg]\n", latt->cell[i].framework.xc, latt->cell[i].framework.yc, latt->cell[i].framework.ae*180./(PI));
		
		//read the keyword (and check that we have not reached the end of the input file)
		fscanf(rfile, "%i", &nb_idcells);
		if(fscanf(rfile, "%s", latt->cell[cell_nb].keyword) == EOF) {
			printf("\n!!!ERROR in load_lattice:  reaching end of file %s\n", inputfilename);
			errorstop("End of file reached unexpectedly, you may try to read more lines than the number of lines in your input file. Please check your input file.");
		}
		if(debug == YES) printf("  --> Cell type: \"%s\"", latt->cell[cell_nb].keyword);
		
		//if keyword == drift, generate drift space
		if(strcmp(latt->cell[cell_nb].keyword, "drift") == 0 || strcmp(latt->cell[cell_nb].keyword, "drift-cyl") == 0) genecell_drift(rfile, &fwk, cell_nb, latt);
		//if keyword == ffag-r-*, generate radial ffag latt with Hard-Edges, linear or Enge fringe field fall-off
		else if(strcmp(latt->cell[cell_nb].keyword, "ffag-r-he") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "ffag-r-lin") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "ffag-r-enge") == 0) {
			genecell_ffagr(rfile, &fwk, cell_nb, latt);
		}
		//if keyword == ffag-s-*, generate ffag-straight latt (keyword ffag-straight)
		else if(strcmp(latt->cell[cell_nb].keyword, "ffag-s-lin") == 0 || 
				strcmp(latt->cell[cell_nb].keyword, "ffag-s-he") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "ffag-s-enge") == 0) {
			genecell_ffagstr(rfile, &fwk, cell_nb, latt);
		}
		//if keyword == field-cylmap or field-cartmap, load field map from data file, with point taken along a cylindrical or cartesian mesh but expressed in cartesian coordinates
		else if(strcmp(latt->cell[cell_nb].keyword, "field-cylmap") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "field-cartmap") == 0) {
			load_map(rfile, &fwk, cell_nb, latt);
		}
		//if keyword == field-cylmap-old or field-cartmap-old, load field map from data file, with point taken along a cylindrical or cartesian mesh but expressed in cartesian coordinates
		else if(strcmp(latt->cell[cell_nb].keyword, "field-cylmap-old") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "field-cartmap-old") == 0) {
			load_map_legacy(rfile, &fwk, cell_nb, latt);
		}
		//if keyword == ffag-sti-*, generate ffag-straight tilt latt (keyword ffag-sti)
		else if(strcmp(latt->cell[cell_nb].keyword, "ffag-sti-lin") == 0 || 
				strcmp(latt->cell[cell_nb].keyword, "ffag-sti-he") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "ffag-sti-enge") == 0) {
			genecell_ffagstrtilt(rfile, &fwk, cell_nb, latt);
		}
		//if keyword == ffag-spi-*, generate ffag-spiral latt (keyword ffag-spi)
		else if(strcmp(latt->cell[cell_nb].keyword, "ffag-spi-lin") == 0 || 
				strcmp(latt->cell[cell_nb].keyword, "ffag-spi-he") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "ffag-spi-enge") == 0 ||
				strcmp(latt->cell[i].keyword, "ffag-spi-fullenge")  == 0) {
			genecell_ffagspiral(rfile, &fwk, cell_nb, latt);
		}
		//if keyword == ffag-tilt-*, generate ffag-tilt latt (keyword ffag-tilt)
		else if(strcmp(latt->cell[cell_nb].keyword, "ffag-tilt-lin") == 0 || 
				strcmp(latt->cell[cell_nb].keyword, "ffag-tilt-he") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "ffag-tilt-enge") == 0) {
			genecell_ffagtilt(rfile, &fwk, cell_nb, latt);
		}
		//if keyword == rf-thingap, generate a RF cavity gap
		else if(strcmp(latt->cell[cell_nb].keyword, "rf-thingap") == 0) {
			genecell_cavity(rfile, cell_nb, latt->periodicity, latt);
		}
		//if keyword == collimator, generate a localized collimator
		else if(strcmp(latt->cell[cell_nb].keyword, "collimator") == 0) {
			genecell_collimator(rfile, cell_nb, latt);
		}
		//if keyword == realbend, generate a bending magnet
		else if(strcmp(latt->cell[cell_nb].keyword, "realbend") == 0) {
			genecell_realdipole(rfile, &fwk, cell_nb, latt);
		}
		//if keyword == purebend, generate a bending magnet
		else if(strcmp(latt->cell[cell_nb].keyword, "purebend-he") == 0 ||
		strcmp(latt->cell[cell_nb].keyword, "purebend-lin") == 0 ||
		strcmp(latt->cell[cell_nb].keyword, "purebend-enge") == 0) {
			genecell_puredipole(rfile, &fwk, cell_nb, latt);
		}
		//if keyword == quad, generate a quadrupole magnet
		else if(strcmp(latt->cell[cell_nb].keyword, "quad-lin") == 0 || 
				strcmp(latt->cell[cell_nb].keyword, "quad-he") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "quad-enge") == 0) {
			genecell_quad(rfile, &fwk, cell_nb, latt);
		}
		//if keyword == sext, generate a sextupole magnet
		else if(strcmp(latt->cell[cell_nb].keyword, "sext-lin") == 0 || 
				strcmp(latt->cell[cell_nb].keyword, "sext-he") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "sext-enge") == 0) {
			genecell_sext(rfile, &fwk, cell_nb, latt);
		}
		//if keyword == oct, generate a sextupole magnet
		else if(strcmp(latt->cell[cell_nb].keyword, "oct-lin") == 0 || 
				strcmp(latt->cell[cell_nb].keyword, "oct-he") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "oct-enge") == 0) {
			genecell_oct(rfile, &fwk, cell_nb, latt);
		}
		//if keyword == thindipole, generate a bending magnet
		else if(strcmp(latt->cell[cell_nb].keyword, "thin-hdipole") == 0) {
			genecell_thindipole(rfile, &fwk, cell_nb, latt);
		}
		//if keyword == float-faraday, generate a floating collimator
		else if(strcmp(latt->cell[cell_nb].keyword, "float-faraday") == 0) {
			genecell_float_faradaycup(rfile, cell_nb, latt);
		}
		//if keyword == vffa-rect, generate a quadrupole magnet
		else if(strcmp(latt->cell[cell_nb].keyword, "vffa-rect-lin") == 0 || 
				strcmp(latt->cell[cell_nb].keyword, "vffa-rect-he") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "vffa-rect-enge-bx") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "vffa-rect-enge-separ") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "vffa-rect-enge") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "vffa-rect-enge-add") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "vffa-rect-atan") == 0) {
			genecell_vffa_rect(rfile, &fwk, cell_nb, latt);
		}
		else if(strcmp(latt->cell[cell_nb].keyword, "vffa-rect-str-he") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "vffa-rect-str-lin") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "vffa-rect-str-enge") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "vffa-rect-atan-str") == 0) {
			genecell_vffa_rect_str(rfile, &fwk, cell_nb, latt);
		}
		else if(strcmp(latt->cell[cell_nb].keyword, "vffa-sect-lin") == 0 || 
				strcmp(latt->cell[cell_nb].keyword, "vffa-sect-he") == 0 ||
				strcmp(latt->cell[cell_nb].keyword, "vffa-sect-enge") == 0) {
			genecell_vffa_sect(rfile, &fwk, cell_nb, latt);
		}
		else {
			printf("\nunrecognized keyword = %s   (while reading %s)\n", latt->cell[cell_nb].keyword, inputfilename);
			errorstop("!!!ERROR in load_lattice: unrecognized keyword");
		}
		
		if(nb_idcells>1) {
			if(debug == YES) printf("  --> Cell type: \"%s\"", latt->cell[cell_nb].keyword);
			duplicate_cell_fmod(cell_nb, &fwk, nb_idcells, latt);
			if(debug == YES) printf("\n");
		}
		cell_nb+=nb_idcells;
	}
	if(cell_nb != latt->nbcell) {
		printf("cell_nb = %i, latt->nbcell = %i\n", cell_nb, latt->nbcell);
		errorstop("problem in load_lattice, nbcells do not match !");
	}
	fclose(rfile);
}

static void genecell_collimator(FILE *rfile, int i, struct Lattice *latt)
{
	newline(rfile);
	if(debug == YES) printf("\n Localized collimator:\n");
	
	fscanf(rfile, "%lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmax));
	latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	if(debug == YES) printf("\n");
}

static void genecell_float_faradaycup(FILE *rfile, int i, struct Lattice *latt)
{
	newline(rfile);
	if(debug == YES) printf("\n Localized faraday cup:\n");
	
	fscanf(rfile, "%lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmax));
	latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	if(debug == YES) printf("WARNING ! These collimators are reversed (kill the beam between xmin and xmax)!!\n");

	if(debug == YES) printf("\n");
}

static void genecell_drift(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	double amax;
	
	newline(rfile);
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("\n  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmax));
	latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	fscanf(rfile, "%lf", &amax);
	newline(rfile);
	if(strcmp(latt->cell[i].keyword, "drift") == 0) {
		latt->cell[i].boun.ymax = amax;
		latt->cell[i].boun.thmax = 0; //disable on th
	}
	else if(strcmp(latt->cell[i].keyword, "drift-cyl") == 0) {
		latt->cell[i].boun.thmax = amax*PI/180.;
		latt->cell[i].boun.ymax = 0; //disable on y
	}
	else errorstop("in genecell_drift, unexpected keyword\n");
	//fscanf(rfile, "%lf", &(latt->cell[i].boun.ymax));
	//newline(rfile);
	//latt->cell[i].boun.thmax = 0; //disable test on theta
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf  (if 0.0: disabled) \t Max_length = %lf [m]\n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));	
	if(debug == YES) printf("\n");
}

static void genecell_thindipole(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	int j;
	//allocate memory for latt->cell[i].mpara
	latt->cell[i].mpara = allocmatrix(1, SIZE_MPARA);
	
	fscanf(rfile, "%lf", &(latt->cell[i].mpara[0][0])); //Int(Bdl)
	newline(rfile);
	for(j=1;j<SIZE_MPARA;j++) latt->cell[i].mpara[0][j]=0;
	if(debug == YES) printf("\n  --> BL = %lf [T.m]\n", latt->cell[i].mpara[0][0]);
	
	latt->cell[i].boun.ymax = 0; //disable test on y
	latt->cell[i].boun.thmax = 0; //disable test on theta
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));	
	if(debug == YES) printf("\n");
}

static void genecell_realdipole(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	double at_deg;
	
	latt->cell[i].nbcomp = 1;
	//allocate memory for latt->cell[i].mpara, latt->cell[i].efben and latt->cell[i].efbex
	latt->cell[i].mpara = allocmatrix(1, SIZE_MPARA);
	latt->cell[i].efben = allocmatrix(1, SIZE_EFB);
	latt->cell[i].efbex = allocmatrix(1, SIZE_EFB);
	//allocate memory
	latt->cell[i].alierror = allocmatrix(latt->cell[i].nbcomp, SIZE_ALIERR);
	
	newline(rfile);
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmax));
	latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	fscanf(rfile, "%lf", &at_deg);
	newline(rfile);
	latt->cell[i].boun.thmax = at_deg*PI/180.;
	latt->cell[i].boun.ymax = 0; //disable on y
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf [rad] \t Max_length = %lf (if 0.0: disabled) \n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].deltar));
	newline(rfile);
	shitf_fwk(&(latt->cell[i].framework), -1*latt->cell[i].deltar);
	//printf("genecell_ffagr: latt->cell[i].framework.xc = %le, latt->cell[i].framework.yc = %le, latt->cell[i].framework.ae = %le\n", latt->cell[i].framework.xc, latt->cell[i].framework.yc, latt->cell[i].framework.ae*180./(PI));
	if(debug == YES) printf("  --> Delta r (origin shift):%le [m] \n", latt->cell[i].deltar);
	
	fscanf(rfile, "%lf", &(latt->cell[i].mpara[0][0]));	//[m] rcentre
	fscanf(rfile, "%lf", &(latt->cell[i].mpara[0][1]));	//[m] FFB = x + rcentre
	fscanf(rfile, "%lf", &(latt->cell[i].mpara[0][2]));	//[m] ffe
	fscanf(rfile, "%lf", &(latt->cell[i].mpara[0][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
	fscanf(rfile, "%lf", &(latt->cell[i].mpara[0][4]));	//[T] B0
	newline(rfile);
	latt->cell[i].mpara[0][5] = 0;	//[] not used here
	latt->cell[i].mpara[0][6] = 0;	//[] not used here
	
	
	fscanf(rfile, "%lf", &(latt->cell[i].efben[0][0])); //[m] FFBen
	fscanf(rfile, "%lf", &(latt->cell[i].efben[0][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
	fscanf(rfile, "%lf", &(latt->cell[i].efben[0][2]));	//[/] kappa in lambda = ffe*(r0/r)^kappa
	fscanf(rfile, "%lf", &(latt->cell[i].efben[0][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
	newline(rfile);
	
	fscanf(rfile, "%lf", &(latt->cell[i].efbex[0][0]));	//[m] FFBex
	fscanf(rfile, "%lf", &(latt->cell[i].efbex[0][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
	fscanf(rfile, "%lf", &(latt->cell[i].efbex[0][2]));	//[/] kappa in lambda = ffe*(r0/r)^kappa
	fscanf(rfile, "%lf", &(latt->cell[i].efbex[0][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
	newline(rfile);
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));	
	if(debug == YES) {
		char cellname[50];
		sprintf(cellname,"cell #%i", i);
		print_cell_para(NULL, cellname, &(latt->cell[i]));
		printf("\n");
	}
}

static void genecell_puredipole(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	int n;
	double ac_deg, at_deg, wm_deg, thfringe_deg;
			
	fscanf(rfile, "%i", &(latt->cell[i].nbcomp));
	newline(rfile);
	if(debug == YES) printf(", %i poles\n", latt->cell[i].nbcomp);
	
	//allocate memory for latt->cell[i].mpara, latt->cell[i].efben and latt->cell[i].efbex
	latt->cell[i].mpara = allocmatrix(latt->cell[i].nbcomp, SIZE_MPARA);
	latt->cell[i].efben = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	latt->cell[i].efbex = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	
	//allocate memory
	latt->cell[i].alierror = allocmatrix(latt->cell[i].nbcomp, SIZE_ALIERR);
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmax));
	latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	fscanf(rfile, "%lf", &at_deg);
	newline(rfile);
	latt->cell[i].boun.thmax = at_deg*PI/180.;
	latt->cell[i].boun.ymax = 0; //disable on y
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf [rad] \t Max_length = %lf (if 0.0: disabled) \n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].deltar));
	newline(rfile);
	shitf_fwk(&(latt->cell[i].framework), -1*latt->cell[i].deltar);
	//printf("genecell_ffagr: latt->cell[i].framework.xc = %le, latt->cell[i].framework.yc = %le, latt->cell[i].framework.ae = %le\n", latt->cell[i].framework.xc, latt->cell[i].framework.yc, latt->cell[i].framework.ae*180./(PI));
	if(debug == YES) printf("  --> Delta r (origin shift):%le [m] \n", latt->cell[i].deltar);
	
	//load nbcomp FFAG poles
	if(debug == YES) printf("  --> subCells parameters (units [m], [rad], [T]):\n");
	for(n = 0; n < latt->cell[i].nbcomp; n++) {
		fscanf(rfile, "%lf", &ac_deg);								//[deg] magnet center position
		latt->cell[i].mpara[n][0] = ac_deg*PI/180.;
		
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][1]));	//[m] r0 in 
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][2]));	//[T] B0
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][3]));	//[/] rc
		latt->cell[i].mpara[n][4] = 0;	//[] not used here
		latt->cell[i].mpara[n][5] = 0;	//[] not used here
		latt->cell[i].mpara[n][6] = 0;	//[] not used here
		
		newline(rfile);
		//Read entrance EFB parameters
		
		fscanf(rfile, "%lf", &wm_deg); //[deg] magnet half oppening angle (position of the effective field boundary (EFB))
		latt->cell[i].efben[n][0] = wm_deg*PI/180.;
		
		if(strcmp(latt->cell[i].keyword, "purebend-he") == 0) {			// hard edge field
			latt->cell[i].efben[n][1] = 0; // not used
			latt->cell[i].efben[n][2] = 0; // not used
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "purebend-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &thfringe_deg);
			latt->cell[i].efben[n][1] = thfringe_deg*PI/180.; //[rad] fringe field (half-)opening angle
			latt->cell[i].efben[n][2] = 0; // not used
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "purebend-enge") == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][2]));	//[/] kappa in lambda = ffe*(r0/r)^kappa
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		newline(rfile);
		//Read exit EFB parameters
		
		fscanf(rfile, "%lf", &wm_deg); //[deg] magnet half oppening angle (position of the effective field boundary (EFB))
		latt->cell[i].efbex[n][0] = wm_deg*PI/180.;
		
		if(strcmp(latt->cell[i].keyword, "purebend-he") == 0) {			// hard edge field
			latt->cell[i].efbex[n][1] = 0; // not used
			latt->cell[i].efbex[n][2] = 0; // not used
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "purebend-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &thfringe_deg);
			latt->cell[i].efbex[n][1] = thfringe_deg*PI/180.; //[rad] fringe field (half-)opening angle
			latt->cell[i].efbex[n][2] = 0; // not used
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "purebend-enge") == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][2]));	//[/] kappa in lambda = ffe*(r0/r)^kappa
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		newline(rfile);
	}
	
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == YES) {
		char cellname[50];
		sprintf(cellname,"cell #%i", i);
		print_cell_para(NULL, cellname, &(latt->cell[i]));
		printf("\n");
	}
}

static void genecell_ffagr(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	int n;
	double ac_deg, at_deg, wm_deg, thfringe_deg;
			
	fscanf(rfile, "%i", &(latt->cell[i].nbcomp));
	newline(rfile);
	if(debug == YES) printf(", %i poles\n", latt->cell[i].nbcomp);
	
	//allocate memory for latt->cell[i].mpara, latt->cell[i].efben and latt->cell[i].efbex
	latt->cell[i].mpara = allocmatrix(latt->cell[i].nbcomp, SIZE_MPARA);
	latt->cell[i].efben = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	latt->cell[i].efbex = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	
	//allocate memory
	latt->cell[i].alierror = allocmatrix(latt->cell[i].nbcomp, SIZE_ALIERR);
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmax));
	latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	fscanf(rfile, "%lf", &at_deg);
	newline(rfile);
	latt->cell[i].boun.thmax = at_deg*PI/180.;
	latt->cell[i].boun.ymax = 0; //disable on y
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf [rad] \t Max_length = %lf (if 0.0: disabled) \n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].deltar));
	newline(rfile);
	shitf_fwk(&(latt->cell[i].framework), -1*latt->cell[i].deltar);
	//printf("genecell_ffagr: latt->cell[i].framework.xc = %le, latt->cell[i].framework.yc = %le, latt->cell[i].framework.ae = %le\n", latt->cell[i].framework.xc, latt->cell[i].framework.yc, latt->cell[i].framework.ae*180./(PI));
	if(debug == YES) printf("  --> Delta r (origin shift):%le [m] \n", latt->cell[i].deltar);
	
	//load nbcomp FFAG poles
	if(debug == YES) printf("  --> subCells parameters (units [m], [rad], [T]):\n");
	for(n = 0; n < latt->cell[i].nbcomp; n++) {
		fscanf(rfile, "%lf", &ac_deg);								//[deg] magnet center position
		latt->cell[i].mpara[n][0] = ac_deg*PI/180.;
		
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][1]));	//[m] r0 in Bz = B0*(r/r0)^k and in lambda = ffe*(r0/r)^kappa
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][2]));	//[T] B0 in Bz = B0*(r/r0)^k
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][3]));	//[/] k in  Bz = B0*(r/r0)^k (geometrical field index)
		latt->cell[i].mpara[n][4] = 0;	//[] not used here
		latt->cell[i].mpara[n][5] = 0;	//[] not used here
		latt->cell[i].mpara[n][6] = 0;	//[] not used here
		
		newline(rfile);
		//Read entrance EFB parameters
		
		fscanf(rfile, "%lf", &wm_deg); //[deg] magnet half oppening angle (position of the effective field boundary (EFB))
		latt->cell[i].efben[n][0] = wm_deg*PI/180.;
		
		if(strcmp(latt->cell[i].keyword, "ffag-r-he") == 0) {			// hard edge field
			latt->cell[i].efben[n][1] = 0; // not used
			latt->cell[i].efben[n][2] = 0; // not used
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-r-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &thfringe_deg);
			latt->cell[i].efben[n][1] = thfringe_deg*PI/180.; //[rad] fringe field (half-)opening angle
			latt->cell[i].efben[n][2] = 0; // not used
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-r-enge") == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][2]));	//[/] kappa in lambda = ffe*(r0/r)^kappa
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		
		newline(rfile);
		//Read exit EFB parameters
		
		fscanf(rfile, "%lf", &wm_deg); //[deg] magnet half oppening angle (position of the effective field boundary (EFB))
		latt->cell[i].efbex[n][0] = wm_deg*PI/180.;
		
		if(strcmp(latt->cell[i].keyword, "ffag-r-he") == 0) {			// hard edge field
			latt->cell[i].efbex[n][1] = 0; // not used
			latt->cell[i].efbex[n][2] = 0; // not used
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-r-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &thfringe_deg);
			latt->cell[i].efbex[n][1] = thfringe_deg*PI/180.; //[rad] fringe field (half-)opening angle
			latt->cell[i].efbex[n][2] = 0; // not used
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-r-enge") == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][2]));	//[/] kappa in lambda = ffe*(r0/r)^kappa
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		
		newline(rfile);
	}
	
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == YES) {
		char cellname[50];
		sprintf(cellname,"cell #%i", i);
		print_cell_para(NULL, cellname, &(latt->cell[i]));
		printf("\n");
	}
}

static void genecell_ffagstr(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	int n;
	
	fscanf(rfile, "%i", &(latt->cell[i].nbcomp));
	newline(rfile);
	if(debug == YES) printf(", %i poles\n", latt->cell[i].nbcomp);
	
	//allocate memory for latt->cell[i].mpara, latt->cell[i].efben and latt->cell[i].efbex
	latt->cell[i].mpara = allocmatrix(latt->cell[i].nbcomp, SIZE_MPARA);
	latt->cell[i].efben = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	latt->cell[i].efbex = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	
	//allocate memory
	latt->cell[i].alierror = allocmatrix(latt->cell[i].nbcomp, SIZE_ALIERR);
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmax));
	latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	latt->cell[i].boun.thmax = 0; //disable on y
	fscanf(rfile, "%lf", &(latt->cell[i].boun.ymax));
	newline(rfile);
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf (if 0.0: disabled) \t Max_length = %lf [m]\n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].deltar));
	newline(rfile);
	shitf_fwk(&(latt->cell[i].framework), -1*latt->cell[i].deltar);
	if(debug == YES) printf("  --> Delta r (origin shift):%le [m] \n", latt->cell[i].deltar);
	
	//load nbcomp FFAG poles
	if(debug == YES) printf("  --> subCells parameters (units [m], [rad], [T]):\n");
	for(n = 0; n < latt->cell[i].nbcomp; n++) {
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][0]));	//[m] magnet center position
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][1]));	//[m] x0 in Bz = B0*exp(kovero*(x-x0))
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][2]));	//[T] B0 in Bz = B0*exp(kovero*(x-x0))
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][3]));	//[1/m] kovero in  Bz = B0*exp(kovero*(x-x0))
		latt->cell[i].mpara[n][4] = 0;	//[] not used here
		latt->cell[i].mpara[n][5] = 0;	//[] not used here
		latt->cell[i].mpara[n][6] = 0;	//[] not used here
		
		newline(rfile);
		//Read entrance EFB parameters
		if(strcmp(latt->cell[i].keyword, "ffag-s-he") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
			latt->cell[i].efben[n][1] = 0;						//[m] not used here
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-s-lin") == 0) {	
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda (=contsant)
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-s-enge") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda
			latt->cell[i].efben[n][2] = 0;						//[/] not used
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)			
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-s-enge-add") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][2]));	//[/] order of extrapolation off-midplane
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)			
		}
			
		newline(rfile);
		//Read exit EFB parameters
		if(strcmp(latt->cell[i].keyword, "ffag-s-he") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
			latt->cell[i].efben[n][1] = 0;						//[m] not used here
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-s-lin") == 0) {	
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda (=contsant)
			latt->cell[i].efbex[n][2] = 0;						//[/] not used here
			latt->cell[i].efbex[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-s-enge") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda
			latt->cell[i].efben[n][2] = 0;						//[/] not used
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-s-enge-add") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][2]));	//[/] order of extrapolation off-midplane
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)			
		}
		newline(rfile);
	}
			
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == YES) {
		char cellname[50];
		sprintf(cellname,"cell #%i", i);
		print_cell_para(NULL, cellname, &(latt->cell[i]));
		printf("\n");
	}
}

static void genecell_ffagspiral(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	int n;
	double ac_deg, at_deg, wm_deg, thfringe_deg, zeta_deg;
	//double rmin, rmax, zmin, zmax;
	
	fscanf(rfile, "%i", &(latt->cell[i].nbcomp));
	newline(rfile);
	if(debug == YES) printf(", %i poles\n", latt->cell[i].nbcomp);
	
	//allocate memory for latt->cell[i].mpara, latt->cell[i].efben and latt->cell[i].efbex
	latt->cell[i].mpara = allocmatrix(latt->cell[i].nbcomp, SIZE_MPARA);
	latt->cell[i].efben = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	latt->cell[i].efbex = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	
	//allocate memory
	latt->cell[i].alierror = allocmatrix(latt->cell[i].nbcomp, SIZE_ALIERR);
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmax));
	latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	//fscanf(rfile, "%lf %lf %lf", &(rmin), &(rmax), &(zmax));
	//zmin = -zmax;
	//printf("coll: (%lf,%lf,%lf,%lf)\n",rmin, rmax, zmin, zmax);
	//latt->cell[i].collim.rmin = rmin;
	//latt->cell[i].collim.rmax = rmax;
	//latt->cell[i].collim.zmin = zmin;
	//latt->cell[i].collim.zmax = zmax;
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	fscanf(rfile, "%lf", &at_deg);
	newline(rfile);
	latt->cell[i].boun.thmax = at_deg*PI/180.;
	latt->cell[i].boun.ymax = 0; //disable on y
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf [rad] \t Max_length = %lf (if 0.0: disabled) \n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].deltar));
	newline(rfile);
	shitf_fwk(&(latt->cell[i].framework), -1*latt->cell[i].deltar);
	//printf("genecell_ffagr: latt->cell[i].framework.xc = %le, latt->cell[i].framework.yc = %le, latt->cell[i].framework.ae = %le\n", latt->cell[i].framework.xc, latt->cell[i].framework.yc, latt->cell[i].framework.ae*180./(PI));
	if(debug == YES) printf("  --> Delta r (origin shift):%le [m] \n", latt->cell[i].deltar);
		
	//load nbcomp FFAG poles
	if(debug == YES) printf("  --> subCells parameters (units [m], [rad], [T]):\n");
	for(n = 0; n < latt->cell[i].nbcomp; n++) {
		fscanf(rfile, "%lf", &ac_deg);						//[deg] magnet center position
		latt->cell[i].mpara[n][0] = ac_deg*PI/180.;
		
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][1]));	//[m] r0 in Bz = B0*(r/r0)^k and in lambda = ffe*(r0/r)^kappa
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][2]));	//[T] B0 in Bz = B0*(r/r0)^k
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][3]));	//[/] k in  Bz = B0*(r/r0)^k (geometrical field index)
		fscanf(rfile, "%lf", &zeta_deg);
		latt->cell[i].mpara[n][4] = zeta_deg*PI/180.;		// spiral/tilt angle zeta0 [rad]
		newline(rfile);
		latt->cell[i].mpara[n][5] = 0;	//[] not used here
		latt->cell[i].mpara[n][6] = 0;	//[] not used here
		
		//Read entrance EFB parameters
		fscanf(rfile, "%lf", &wm_deg); //[deg] magnet half oppening angle (position of the effective field boundary (EFB))
		latt->cell[i].efben[n][0] = wm_deg*PI/180.;
		
		if(strcmp(latt->cell[i].keyword, "ffag-spi-he") == 0) {			// hard edge field
			latt->cell[i].efben[n][1] = 0; // not used
			latt->cell[i].efben[n][2] = 0; // not used
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-spi-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &thfringe_deg);
			latt->cell[i].efben[n][1] = thfringe_deg*PI/180.; //[rad] fringe field (half-)opening angle
			latt->cell[i].efben[n][2] = 0; // not used
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-spi-enge")  == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][2]));	//[/]  Order of interpolation off the mid-plane
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-spi-fullenge")  == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][2]));	//[/]  Order of interpolation off the mid-plane
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C0 coefficient (related to the fringe field shape)
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][4]));	//[/] Enge C1 coefficient (related to the fringe field shape)
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][5]));	//[/] Enge C2 coefficient (related to the fringe field shape)
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][6]));	//[/] Enge C3 coefficient (related to the fringe field shape)
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][7]));	//[/] Enge C4 coefficient (related to the fringe field shape)
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][8]));	//[/] Enge C5 coefficient (related to the fringe field shape)
		}
		newline(rfile);
		//Read exit EFB parameters
		
		fscanf(rfile, "%lf", &wm_deg); //[deg] magnet half oppening angle (position of the effective field boundary (EFB))
		latt->cell[i].efbex[n][0] = wm_deg*PI/180.;
		
		if(strcmp(latt->cell[i].keyword, "ffag-spi-he") == 0) {			// hard edge field
			latt->cell[i].efbex[n][1] = 0; // not used
			latt->cell[i].efbex[n][2] = 0; // not used
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-spi-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &thfringe_deg);
			latt->cell[i].efbex[n][1] = thfringe_deg*PI/180.; //[rad] fringe field (half-)opening angle
			latt->cell[i].efbex[n][2] = 0; // not used
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-spi-enge")  == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][2]));	//[/]  Order of interpolation off the mid-plane
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-spi-fullenge")  == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][2]));	//[/]  Order of interpolation off the mid-plane
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C0 coefficient (related to the fringe field shape)
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][4]));	//[/] Enge C1 coefficient (related to the fringe field shape)
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][5]));	//[/] Enge C2 coefficient (related to the fringe field shape)
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][6]));	//[/] Enge C3 coefficient (related to the fringe field shape)
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][7]));	//[/] Enge C4 coefficient (related to the fringe field shape)
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][8]));	//[/] Enge C5 coefficient (related to the fringe field shape)
		}
		newline(rfile);
	}
	
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == YES) {
		char cellname[50];
		sprintf(cellname,"cell #%i", i);
		print_cell_para(NULL, cellname, &(latt->cell[i]));
		printf("\n");
	}
}

static void genecell_cavity(FILE *rfile, int i, int periodicity, struct Lattice *latt)
{
	int j;
	double phi0_deg;
	
	//allocate memory
	latt->cell[i].cav = alloccav(periodicity);
	
	//fscanf(rfile, "%s", temp_str);
	//printf("  --> cavity keywork: %s\n", temp_str);
	
	for(j = 0; j < periodicity; j++) {
		//strcpy(latt->cell[i].cav[j].keyword, temp_str);
		newline(rfile);
		fscanf(rfile, "%lf", &(latt->cell[i].cav[j].frf));
		newline(rfile);
		fscanf(rfile, "%lf", &(latt->cell[i].cav[j].v0));
		newline(rfile);
		fscanf(rfile, "%lf", &(phi0_deg));
		latt->cell[i].cav[j].phi0 = phi0_deg*(PI)/180.;
		newline(rfile);
		fscanf(rfile, "%lf", &(latt->cell[i].cav[j].t0));
		newline(rfile);
		//latt->cell[i].cav[j].frf	= 74.991737e06;
		//latt->cell[i].cav[j].v0		= -0.68e06;
		//latt->cell[i].cav[j].phi0	= 30*PI/180.;
		//latt->cell[i].cav[j].t0		= 0;
	}
}

static void genecell_ffagtilt(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	int n;
	double ac_deg, at_deg, wm_deg, thfringe_deg, zeta_deg;
	
	fscanf(rfile, "%i", &(latt->cell[i].nbcomp));
	newline(rfile);
	printf(", %i poles\n", latt->cell[i].nbcomp);
	
	//allocate memory for latt->cell[i].mpara, latt->cell[i].efben and latt->cell[i].efbex
	latt->cell[i].mpara = allocmatrix(latt->cell[i].nbcomp, SIZE_MPARA);
	latt->cell[i].efben = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	latt->cell[i].efbex = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	
	//allocate memory
	latt->cell[i].alierror = allocmatrix(latt->cell[i].nbcomp, SIZE_ALIERR);
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	printf("  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmax));
	latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
		
	fscanf(rfile, "%lf", &at_deg);
	newline(rfile);
	latt->cell[i].boun.thmax = at_deg*PI/180.;
	latt->cell[i].boun.ymax = 0; //disable on y
	printf("  --> Cell boundary:\n\tMax_angle = %lf [rad] \t Max_length = %lf (if 0.0: disabled) \n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].deltar));
	newline(rfile);
	shitf_fwk(&(latt->cell[i].framework), -1*latt->cell[i].deltar);
	//printf("genecell_ffagr: latt->cell[i].framework.xc = %le, latt->cell[i].framework.yc = %le, latt->cell[i].framework.ae = %le\n", latt->cell[i].framework.xc, latt->cell[i].framework.yc, latt->cell[i].framework.ae*180./(PI));
	printf("  --> Delta r (origin shift):%le [m] \n", latt->cell[i].deltar);
	
	//load nbcomp FFAG poles
	printf("  --> subCells parameters (units [m], [rad], [T]):\n");
	for(n = 0; n < latt->cell[i].nbcomp; n++) {
		fscanf(rfile, "%lf", &ac_deg);						//[deg] magnet center position
		latt->cell[i].mpara[n][0] = ac_deg*PI/180.;
		
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][1]));	//[m] r0 in Bz = B0*(r/r0)^k and in lambda = ffe*(r0/r)^kappa
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][2]));	//[T] B0 in Bz = B0*(r/r0)^k
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][3]));	//[/] k1 in  Bz = B0*(r/r0)^k1 (geometrical field index)
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][4]));	//[/] k2 in  Bz = B0*exp^(k2*(zeta0 - arcsin())) (geometrical field index)
		fscanf(rfile, "%lf", &zeta_deg);
		latt->cell[i].mpara[n][5] = zeta_deg*PI/180.;		// spiral/tilt angle zeta0 [rad]
		newline(rfile);
		latt->cell[i].mpara[n][6] = 0;	//[] not used here
		
		//Read entrance EFB parameters
		fscanf(rfile, "%lf", &wm_deg); //[deg] magnet half oppening angle (position of the effective field boundary (EFB))
		latt->cell[i].efben[n][0] = wm_deg*PI/180.;
		
		if(strcmp(latt->cell[i].keyword, "ffag-tilt-he") == 0) {			// hard edge field
			latt->cell[i].efben[n][1] = 0; // not used
			latt->cell[i].efben[n][2] = 0; // not used
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-tilt-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &thfringe_deg);
			latt->cell[i].efben[n][1] = thfringe_deg*PI/180.; //[rad] fringe field (half-)opening angle
			latt->cell[i].efben[n][2] = 0; // not used
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-tilt-enge") == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][2]));	//[/] kappa in lambda = ffe*(r0/r)^kappa
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		
		newline(rfile);
		//Read exit EFB parameters
		
		fscanf(rfile, "%lf", &wm_deg); //[deg] magnet half oppening angle (position of the effective field boundary (EFB))
		latt->cell[i].efbex[n][0] = wm_deg*PI/180.;
		
		if(strcmp(latt->cell[i].keyword, "ffag-tilt-he") == 0) {			// hard edge field
			latt->cell[i].efbex[n][1] = 0; // not used
			latt->cell[i].efbex[n][2] = 0; // not used
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-tilt-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &thfringe_deg);
			latt->cell[i].efbex[n][1] = thfringe_deg*PI/180.; //[rad] fringe field (half-)opening angle
			latt->cell[i].efbex[n][2] = 0; // not used
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-tilt-enge") == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][2]));	//[/] kappa in lambda = ffe*(r0/r)^kappa
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		newline(rfile);
	}
	
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == YES) {
		char cellname[50];
		sprintf(cellname,"cell #%i", i);
		print_cell_para(NULL, cellname, &(latt->cell[i]));
		printf("\n");
	}
}

static void genecell_ffagstrtilt(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	int n;
	double zeta_deg;
	
	fscanf(rfile, "%i", &(latt->cell[i].nbcomp));
	newline(rfile);
	if(debug == YES) printf(", %i poles\n", latt->cell[i].nbcomp);
	
	//allocate memory for latt->cell[i].mpara, latt->cell[i].efben and latt->cell[i].efbex
	latt->cell[i].mpara = allocmatrix(latt->cell[i].nbcomp, SIZE_MPARA);
	latt->cell[i].efben = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	latt->cell[i].efbex = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	
	//allocate memory
	latt->cell[i].alierror = allocmatrix(latt->cell[i].nbcomp, SIZE_ALIERR);
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmax));
	latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	latt->cell[i].boun.thmax = 0; //disable on th
	fscanf(rfile, "%lf", &(latt->cell[i].boun.ymax));
	newline(rfile);
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf (if 0.0: disabled) \t Max_length = %lf [m]\n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].deltar));
	newline(rfile);
	shitf_fwk(&(latt->cell[i].framework), -1*latt->cell[i].deltar);
	if(debug == YES) printf("  --> Delta r (origin shift):%le [m] \n", latt->cell[i].deltar);
	
	//load nbcomp FFAG poles
	if(debug == YES) printf("  --> subCells parameters (units [m], [rad], [T]):\n");
	for(n = 0; n < latt->cell[i].nbcomp; n++) {
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][0]));	//[m] magnet center position at x0
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][1]));	//[m] x0 in Bz = B0*exp(kovero*(x-x0))
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][2]));	//[T] B0 in Bz = B0*exp(kovero*(x-x0))
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][3]));	//[1/m] kovero in  Bz = B0*exp(kovero*(x-x0))
		fscanf(rfile, "%lf", &zeta_deg); // [deg] zeta angle of the tilted magnet  
		latt->cell[i].mpara[n][4] = zeta_deg*PI/180.;
		newline(rfile);
		latt->cell[i].mpara[n][5] = 0;	//[] not used here
		latt->cell[i].mpara[n][6] = 0;	//[] not used here
		
		//Read entrance EFB parameters
		fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
		
		if(strcmp(latt->cell[i].keyword, "ffag-sti-he") == 0) {
			latt->cell[i].efben[n][1] = 0;						//[m] not used here
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-sti-lin") == 0) {	
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda (=contsant)
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-sti-enge") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda
			latt->cell[i].efben[n][2] = 0;						//[/] not used
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)			
		}
		
		newline(rfile);
		//Read exit EFB parameters
		fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
		
		if(strcmp(latt->cell[i].keyword, "ffag-sti-he") == 0) {
			latt->cell[i].efben[n][1] = 0;						//[m] not used here
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-sti-lin") == 0) {	
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda (=contsant)
			latt->cell[i].efbex[n][2] = 0;						//[/] not used here
			latt->cell[i].efbex[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "ffag-sti-enge") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda
			latt->cell[i].efben[n][2] = 0;						//[/] not used
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		newline(rfile);
	}
	
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == YES) {
		char cellname[50];
		sprintf(cellname,"cell #%i", i);
		print_cell_para(NULL, cellname, &(latt->cell[i]));
		printf("\n");
	}
}

static void genecell_quad(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	int n;
			
	fscanf(rfile, "%i", &(latt->cell[i].nbcomp));
	newline(rfile);
	if(debug == YES) printf(", %i quads\n", latt->cell[i].nbcomp);
	
	//allocate memory for latt->cell[i].mpara, latt->cell[i].efben and latt->cell[i].efbex
	latt->cell[i].mpara = allocmatrix(latt->cell[i].nbcomp, SIZE_MPARA);
	latt->cell[i].efben = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	latt->cell[i].efbex = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	
	//allocate memory
	latt->cell[i].alierror = allocmatrix(latt->cell[i].nbcomp, SIZE_ALIERR);

	
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmax));
	latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	latt->cell[i].boun.thmax = 0.; //disable
	fscanf(rfile, "%lf", &(latt->cell[i].boun.ymax));
	newline(rfile);
	
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf [rad] \t Max_length = %lf (if 0.0: disabled) \n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].deltar));
	newline(rfile);
	shitf_fwk(&(latt->cell[i].framework), -1*latt->cell[i].deltar);
	//printf("genecell_ffagr: latt->cell[i].framework.xc = %le, latt->cell[i].framework.yc = %le, latt->cell[i].framework.ae = %le\n", latt->cell[i].framework.xc, latt->cell[i].framework.yc, latt->cell[i].framework.ae*180./(PI));
	if(debug == YES) printf("  --> Delta r (origin shift):%le [m] \n", latt->cell[i].deltar);
	
	//load nb quads
	if(debug == YES) printf("  --> subCells parameters (units [m], [rad], [T]):\n");
	for(n = 0; n < latt->cell[i].nbcomp; n++) {
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][0]));	//[m] magnet center position
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][1]));	//[m] x0 of B(x0)=0
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][2]));	//[T/m] gradient g in Bz=-g*(x-x0) and Bx=-g*z
		latt->cell[i].mpara[n][3] = 0;	//[] not used here
		latt->cell[i].mpara[n][4] = 0;	//[] not used here
		latt->cell[i].mpara[n][5] = 0;	//[] not used here
		latt->cell[i].mpara[n][6] = 0;	//[] not used here
		
		newline(rfile);
		//Read entrance EFB parameters
		if(strcmp(latt->cell[i].keyword, "quad-he") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
			latt->cell[i].efben[n][1] = 0;						//[m] not used here
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "quad-lin") == 0) {	
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda (=constant)
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "quad-enge") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda
			latt->cell[i].efben[n][2] = 0;						//[/] not used
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)			
		}
			
		newline(rfile);
		//Read exit EFB parameters
		if(strcmp(latt->cell[i].keyword, "quad-he") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
			latt->cell[i].efben[n][1] = 0;						//[m] not used here
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "quad-lin") == 0) {	
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda (=contsant)
			latt->cell[i].efbex[n][2] = 0;						//[/] not used here
			latt->cell[i].efbex[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "quad-enge") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda
			latt->cell[i].efben[n][2] = 0;						//[/] not used
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		newline(rfile);
	}
	
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == YES) {
		char cellname[50];
		sprintf(cellname,"cell #%i", i);
		print_cell_para(NULL, cellname, &(latt->cell[i]));
		printf("\n");
	}
}

static void genecell_sext(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	int n;
	
	fscanf(rfile, "%i", &(latt->cell[i].nbcomp));
	newline(rfile);
	if(debug == YES) printf(", %i sextupoles\n", latt->cell[i].nbcomp);
	
	//allocate memory for latt->cell[i].mpara, latt->cell[i].efben and latt->cell[i].efbex
	latt->cell[i].mpara = allocmatrix(latt->cell[i].nbcomp, SIZE_MPARA);
	latt->cell[i].efben = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	latt->cell[i].efbex = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	
	//allocate memory
	latt->cell[i].alierror = allocmatrix(latt->cell[i].nbcomp, SIZE_ALIERR);

	
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmax));
	latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	latt->cell[i].boun.thmax = 0.; //disable
	fscanf(rfile, "%lf", &(latt->cell[i].boun.ymax));
	newline(rfile);
	
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf [rad] \t Max_length = %lf (if 0.0: disabled) \n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].deltar));
	newline(rfile);
	shitf_fwk(&(latt->cell[i].framework), -1*latt->cell[i].deltar);
	if(debug == YES) printf("  --> Delta r (origin shift):%le [m] \n", latt->cell[i].deltar);
	
	//load nb sexts
	if(debug == YES) printf("  --> subCells parameters (units [m], [rad], [T]):\n");
	for(n = 0; n < latt->cell[i].nbcomp; n++) {
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][0]));	//[m] magnet center position
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][1]));	//[m] x0 of B(x0)=0
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][2]));	//[T/m^2] gradient g in Bz=g/2*(x^2-z^2) and Bx=g*xz
		latt->cell[i].mpara[n][3] = 0;	//[] not used here
		latt->cell[i].mpara[n][4] = 0;	//[] not used here
		latt->cell[i].mpara[n][5] = 0;	//[] not used here
		latt->cell[i].mpara[n][6] = 0;	//[] not used here
		
		newline(rfile);
		//Read entrance EFB parameters
		if(strcmp(latt->cell[i].keyword, "sext-he") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
			latt->cell[i].efben[n][1] = 0;						//[m] not used here
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "sext-lin") == 0) {	
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda (=constant)
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "sext-enge") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda
			latt->cell[i].efben[n][2] = 0;						//[/] not used
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)			
		}
			
		newline(rfile);
		//Read exit EFB parameters
		if(strcmp(latt->cell[i].keyword, "sext-he") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
			latt->cell[i].efben[n][1] = 0;						//[m] not used here
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "sext-lin") == 0) {	
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda (=contsant)
			latt->cell[i].efbex[n][2] = 0;						//[/] not used here
			latt->cell[i].efbex[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "sext-enge") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda
			latt->cell[i].efben[n][2] = 0;						//[/] not used
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		newline(rfile);
	}
	
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == YES) {
		char cellname[50];
		sprintf(cellname,"cell #%i", i);
		print_cell_para(NULL, cellname, &(latt->cell[i]));
		printf("\n");
	}
}

static void genecell_oct(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	int n;
			
	fscanf(rfile, "%i", &(latt->cell[i].nbcomp));
	newline(rfile);
	if(debug == YES) printf(", %i octupoles\n", latt->cell[i].nbcomp);
	
	//allocate memory for latt->cell[i].mpara, latt->cell[i].efben and latt->cell[i].efbex
	latt->cell[i].mpara = allocmatrix(latt->cell[i].nbcomp, SIZE_MPARA);
	latt->cell[i].efben = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	latt->cell[i].efbex = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	
	//allocate memory
	latt->cell[i].alierror = allocmatrix(latt->cell[i].nbcomp, SIZE_ALIERR);

	
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmax));
	latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	latt->cell[i].boun.thmax = 0.; //disable
	fscanf(rfile, "%lf", &(latt->cell[i].boun.ymax));
	newline(rfile);
	
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf [rad] \t Max_length = %lf (if 0.0: disabled) \n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].deltar));
	newline(rfile);
	shitf_fwk(&(latt->cell[i].framework), -1*latt->cell[i].deltar);
	if(debug == YES) printf("  --> Delta r (origin shift):%le [m] \n", latt->cell[i].deltar);
	
	//load nb sexts
	if(debug == YES) printf("  --> subCells parameters (units [m], [rad], [T]):\n");
	for(n = 0; n < latt->cell[i].nbcomp; n++) {
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][0]));	//[m] magnet center position
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][1]));	//[m] x0 of B(x0)=0
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][2]));	//[T/m^3] gradient g in (zgoubi: BX =0, BY =g*(3Y^2−Z^2)*Z, BZ =g*(Y^2−3*Z^2)*Y.)
		latt->cell[i].mpara[n][3] = 0;	//[] not used here
		latt->cell[i].mpara[n][4] = 0;	//[] not used here
		latt->cell[i].mpara[n][5] = 0;	//[] not used here
		latt->cell[i].mpara[n][6] = 0;	//[] not used here
		
		newline(rfile);
		//Read entrance EFB parameters
		if(strcmp(latt->cell[i].keyword, "sext-he") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
			latt->cell[i].efben[n][1] = 0;						//[m] not used here
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "sext-lin") == 0) {	
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda (=constant)
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "sext-enge") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda
			latt->cell[i].efben[n][2] = 0;						//[/] not used
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)			
		}
			
		newline(rfile);
		//Read exit EFB parameters
		if(strcmp(latt->cell[i].keyword, "sext-he") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
			latt->cell[i].efben[n][1] = 0;						//[m] not used here
			latt->cell[i].efben[n][2] = 0;						//[/] not used here
			latt->cell[i].efben[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "sext-lin") == 0) {	
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda (=contsant)
			latt->cell[i].efbex[n][2] = 0;						//[/] not used here
			latt->cell[i].efbex[n][3] = 0;						//[/] not used here
		}
		else if(strcmp(latt->cell[i].keyword, "sext-enge") == 0) {
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0]));	//[m] magnet half length
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda
			latt->cell[i].efben[n][2] = 0;						//[/] not used
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		newline(rfile);
	}
	
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == YES) {
		char cellname[50];
		sprintf(cellname,"cell #%i", i);
		print_cell_para(NULL, cellname, &(latt->cell[i]));
		printf("\n");
	}
}

static void genecell_vffa_rect(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	int n;
	double ac_deg, at_deg, tilt_angle_deg;
	
	fscanf(rfile, "%i", &(latt->cell[i].nbcomp));
	newline(rfile);
	if(debug == YES) printf(", %i poles\n", latt->cell[i].nbcomp);
	
	//allocate memory for latt->cell[i].mpara, latt->cell[i].efben and latt->cell[i].efbex
	latt->cell[i].mpara = allocmatrix(latt->cell[i].nbcomp, SIZE_MPARA);
	latt->cell[i].efben = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	latt->cell[i].efbex = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	
	//allocate memory
	latt->cell[i].alierror = allocmatrix(latt->cell[i].nbcomp, SIZE_ALIERR);
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmin), &(latt->cell[i].collim.zmax));
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	fscanf(rfile, "%lf", &at_deg);
	newline(rfile);
	latt->cell[i].boun.thmax = at_deg*PI/180.;
	latt->cell[i].boun.ymax = 0; //disable on y
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf [rad] \t Max_length = %lf (if 0.0: disabled) \n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].deltar));
	newline(rfile);
	shitf_fwk(&(latt->cell[i].framework), -1*latt->cell[i].deltar);
	if(debug == YES) printf("  --> Delta r (origin shift):%le [m] \n", latt->cell[i].deltar);
	
	//load nbcomp VFFA poles
	if(debug == YES) printf("  --> subCells parameters (units [m], [rad], [T]):\n");
	for(n = 0; n < latt->cell[i].nbcomp; n++) {
		fscanf(rfile, "%lf", &ac_deg);								//[deg] magnet center position
		latt->cell[i].mpara[n][0] = ac_deg*PI/180.;
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][1]));	//[m] r0 [m] magnet center radius (also used in lambda = ffe*(r0/r)^kappa)
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][2]));	//[T] B0 in Bz = B0*exp(m(z-z0))
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][3]));	//[/m] z0 in  Bz = B0*exp(m(z-z0)) 
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][4]));	//[/m] m in  Bz = B0*exp(m(z-z0)) (normalized field gradient)
		fscanf(rfile, "%lf", &tilt_angle_deg);				//[deg] tilt angle
		latt->cell[i].mpara[n][5] = tilt_angle_deg*PI/180.;
		latt->cell[i].mpara[n][6] = 0;	//[] convergence limit, computed later
		newline(rfile);
		//Read entrance EFB parameters
		fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0])); //[m] magnet half oppening length (position of the effective field boundary (EFB))
		if(strcmp(latt->cell[i].keyword, "vffa-rect-he") == 0) {			// hard edge field
			latt->cell[i].efben[n][1] = 0; // not used
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][2]));	//[/]  Order of interpolation off the mid-plane
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "vffa-rect-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1])); //[m] fringe field (half-)opening length
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][2]));
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "vffa-rect-enge") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-enge-bx") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-enge-add") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-atan") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-enge-separ") == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][2]));	//[/]  Order of interpolation off the mid-plane
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		newline(rfile);
		//Read exit EFB parameters
		fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0])); //[m] magnet half oppening length (position of the effective field boundary (EFB))
		if(strcmp(latt->cell[i].keyword, "vffa-rect-he") == 0) {			// hard edge field
			latt->cell[i].efbex[n][1] = 0; // not used
			latt->cell[i].efbex[n][2] = 0; // not used
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "vffa-rect-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1])); //[m] fringe field (half-)opening length
			latt->cell[i].efbex[n][2] = 0; // not used			
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "vffa-rect-enge") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-enge-bx") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-enge-add") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-atan") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-enge-separ") == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][2]));	//[/] Order of interpolation off the mid-plane
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		newline(rfile);
	}
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == YES) {
		char cellname[50];
		sprintf(cellname,"cell #%i", i);
		print_cell_para(NULL, cellname, &(latt->cell[i]));
		printf("\n");
	}
}

static void genecell_vffa_rect_str(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	int n;
	double tilt_angle_deg;
	
	fscanf(rfile, "%i", &(latt->cell[i].nbcomp));
	newline(rfile);
	if(debug == YES) printf(", %i poles\n", latt->cell[i].nbcomp);
	
	//allocate memory for latt->cell[i].mpara, latt->cell[i].efben and latt->cell[i].efbex
	latt->cell[i].mpara = allocmatrix(latt->cell[i].nbcomp, SIZE_MPARA);
	latt->cell[i].efben = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	latt->cell[i].efbex = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	
	//allocate memory
	latt->cell[i].alierror = allocmatrix(latt->cell[i].nbcomp, SIZE_ALIERR);
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmin), &(latt->cell[i].collim.zmax));
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\txmin = %lf [m] xmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].boun.ymax));
	newline(rfile);
	latt->cell[i].boun.thmax = 0;
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf [rad] \t Max_length = %lf (if 0.0: disabled) \n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].deltar));
	newline(rfile);
	shitf_fwk(&(latt->cell[i].framework), -1*latt->cell[i].deltar);
	if(debug == YES) printf("  --> Delta r (origin shift):%le [m] \n", latt->cell[i].deltar);
	
	//load nbcomp VFFA poles
	if(debug == YES) printf("  --> subCells parameters (units [m], [rad], [T]):\n");
	for(n = 0; n < latt->cell[i].nbcomp; n++) {
		fscanf(rfile, "%lf", &latt->cell[i].mpara[n][0]);	//[m] magnet center position
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][1]));	//[m] x0 [m] magnet center radius (also used in lambda = ffe*(r0/r)^kappa)
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][2]));	//[T] B0 in Bz = B0*exp(m(z-z0))
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][3]));	//[m] z0 in  Bz = B0*exp(m(z-z0)) 
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][4]));	//[/m] m in  Bz = B0*exp(m(z-z0)) (normalized field gradient)
		fscanf(rfile, "%lf", &tilt_angle_deg);				//[deg] tilt angle
		latt->cell[i].mpara[n][5] = tilt_angle_deg*PI/180.;
		latt->cell[i].mpara[n][6] = 0;	//[] convergence limit, computed later
		newline(rfile);
		//Read entrance EFB parameters
		fscanf(rfile, "%lf", &(latt->cell[i].efben[n][0])); //[m] magnet half oppening length (position of the effective field boundary (EFB))
		if(strcmp(latt->cell[i].keyword, "vffa-rect-str-he") == 0) {			// hard edge field
			latt->cell[i].efben[n][1] = 0; // not used
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][2]));	//[/]  Order of interpolation off the mid-plane
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "vffa-rect-str-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1])); //[m] fringe field (half-)opening length
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][2])); //[/]  Order of interpolation off the mid-plane
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "vffa-rect-str-enge") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-atan-str") == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][2]));	//[/]  Order of interpolation off the mid-plane
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		newline(rfile);
		//Read exit EFB parameters
		fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][0])); //[m] magnet half oppening length (position of the effective field boundary (EFB))
		if(strcmp(latt->cell[i].keyword, "vffa-rect-str-he") == 0) {			// hard edge field
			latt->cell[i].efbex[n][1] = 0; // not used
			latt->cell[i].efbex[n][2] = 0; // not used
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "vffa-rect-str-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1])); //[m] fringe field (half-)opening length
			latt->cell[i].efbex[n][2] = 0; // not used			
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "vffa-rect-str-enge") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-atan-str") == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][2]));	//[/] Order of interpolation off the mid-plane
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		newline(rfile);
	}
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == YES) {
		char cellname[50];
		sprintf(cellname,"cell #%i", i);
		print_cell_para(NULL, cellname, &(latt->cell[i]));
		printf("\n");
	}
}

static void genecell_vffa_sect(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	int n;
	double ac_deg, at_deg, wm_deg, thfringe_deg;
			
	fscanf(rfile, "%i", &(latt->cell[i].nbcomp));
	newline(rfile);
	if(debug == YES) printf(", %i poles\n", latt->cell[i].nbcomp);
	
	//allocate memory for latt->cell[i].mpara, latt->cell[i].efben and latt->cell[i].efbex
	latt->cell[i].mpara = allocmatrix(latt->cell[i].nbcomp, SIZE_MPARA);
	latt->cell[i].efben = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	latt->cell[i].efbex = allocmatrix(latt->cell[i].nbcomp, SIZE_EFB);
	
	//allocate memory
	latt->cell[i].alierror = allocmatrix(latt->cell[i].nbcomp, SIZE_ALIERR);
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("  --> stepsize: %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf %lf %lf %lf", &(latt->cell[i].collim.rmin), &(latt->cell[i].collim.rmax), &(latt->cell[i].collim.zmin), &(latt->cell[i].collim.zmax));
	newline(rfile);
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	fscanf(rfile, "%lf", &at_deg);
	newline(rfile);
	latt->cell[i].boun.thmax = at_deg*PI/180.;
	latt->cell[i].boun.ymax = 0; //disable on y
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf [rad] \t Max_length = %lf (if 0.0: disabled) \n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	fscanf(rfile, "%lf", &(latt->cell[i].deltar));
	newline(rfile);
	shitf_fwk(&(latt->cell[i].framework), -1*latt->cell[i].deltar);
	//printf("genecell_ffagr: latt->cell[i].framework.xc = %le, latt->cell[i].framework.yc = %le, latt->cell[i].framework.ae = %le\n", latt->cell[i].framework.xc, latt->cell[i].framework.yc, latt->cell[i].framework.ae*180./(PI));
	if(debug == YES) printf("  --> Delta r (origin shift):%le [m] \n", latt->cell[i].deltar);
	
	//load nbcomp VFFA poles
	if(debug == YES) printf("  --> subCells parameters (units [m], [rad], [T]):\n");
	for(n = 0; n < latt->cell[i].nbcomp; n++) {
		fscanf(rfile, "%lf", &ac_deg);								//[deg] magnet center position
		latt->cell[i].mpara[n][0] = ac_deg*PI/180.;
		
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][1]));	//[m] r0 of the centre of the sector magnet ()
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][2]));	//[T] B0 in Bz = B0*exp(m(z-z0))
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][3]));	//[/m] z0 in  Bz = B0*exp(m(z-z0)) 
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][4]));	//[/m] m in  Bz = B0*exp(m(z-z0)) (normalized field gradient)
		fscanf(rfile, "%lf", &(latt->cell[i].mpara[n][5]));	//[m] rho_0 in Bz = B0*exp(m(z-z0)) and in lambda = ffe*(rho_0/r)^kappa
		latt->cell[i].mpara[n][6] = 0;	//[] not used here
		
		newline(rfile);
		//Read entrance EFB parameters
		
		fscanf(rfile, "%lf", &wm_deg); //[deg] magnet half oppening angle (position of the effective field boundary (EFB))
		latt->cell[i].efben[n][0] = wm_deg*PI/180.;
		
		if(strcmp(latt->cell[i].keyword, "vffa-sect-he") == 0) {			// hard edge field
			latt->cell[i].efben[n][1] = 0; // not used
			latt->cell[i].efben[n][2] = 0; // not used
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "vffa-sect-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &thfringe_deg);
			latt->cell[i].efben[n][1] = thfringe_deg*PI/180.; //[rad] fringe field (half-)opening angle
			latt->cell[i].efben[n][2] = 0; // not used
			latt->cell[i].efben[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "vffa-sect-enge") == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][2]));	//[/] kappa in lambda = ffe*(r0/r)^kappa
			fscanf(rfile, "%lf", &(latt->cell[i].efben[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		
		newline(rfile);
		//Read exit EFB parameters
		
		fscanf(rfile, "%lf", &wm_deg); //[deg] magnet half oppening angle (position of the effective field boundary (EFB))
		latt->cell[i].efbex[n][0] = wm_deg*PI/180.;
		
		if(strcmp(latt->cell[i].keyword, "vffa-sect-he") == 0) {			// hard edge field
			latt->cell[i].efbex[n][1] = 0; // not used
			latt->cell[i].efbex[n][2] = 0; // not used
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "vffa-sect-lin") == 0) {		// soft edge with linear fringe field falloff
			fscanf(rfile, "%lf", &thfringe_deg);
			latt->cell[i].efbex[n][1] = thfringe_deg*PI/180.; //[rad] fringe field (half-)opening angle
			latt->cell[i].efbex[n][2] = 0; // not used
			latt->cell[i].efbex[n][3] = 0; // not used
		}
		else if(strcmp(latt->cell[i].keyword, "vffa-sect-enge") == 0) {	// soft edge with Enge fringe falloff
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][1]));	//[m] ffe in lambda = ffe*(r0/r)^kappa (lambda = fringe field extend)
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][2]));	//[/] kappa in lambda = ffe*(r0/r)^kappa
			fscanf(rfile, "%lf", &(latt->cell[i].efbex[n][3]));	//[/] Enge C1 coefficient (related to the fringe field shape)
		}
		
		newline(rfile);
	}
	
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == YES) {
		char cellname[50];
		sprintf(cellname,"cell #%i", i);
		print_cell_para(NULL, cellname, &(latt->cell[i]));
		printf("\n");
	}
}

extern int test_cell_rect_vffa_bend_cell(struct Cell *cell) {
	if( strcmp(cell->keyword, "vffa-rect-lin") == 0 || 
		strcmp(cell->keyword, "vffa-rect-he") == 0 ||
		strcmp(cell->keyword, "vffa-rect-enge-bx") == 0 ||
		strcmp(cell->keyword, "vffa-rect-enge-separ") == 0 ||
		strcmp(cell->keyword, "vffa-rect-enge") == 0 ||
		strcmp(cell->keyword, "vffa-rect-enge-add") == 0 ||
		strcmp(cell->keyword, "vffa-rect-atan") == 0) return YES;
	else return NO;
}

extern int test_cell_vffa(struct Cell *cell) {
	if( strcmp(cell->keyword, "vffa-rect-lin") == 0 || 
		strcmp(cell->keyword, "vffa-rect-he") == 0 ||
		strcmp(cell->keyword, "vffa-rect-enge-bx") == 0 ||
		strcmp(cell->keyword, "vffa-rect-enge-separ") == 0 ||
		strcmp(cell->keyword, "vffa-rect-enge") == 0 ||
		strcmp(cell->keyword, "vffa-rect-enge-add") == 0 ||
		strcmp(cell->keyword, "vffa-rect-atan") == 0 ||
		strcmp(cell->keyword, "vffa-rect-str-he") == 0 ||
		strcmp(cell->keyword, "vffa-rect-str-lin") == 0 ||
		strcmp(cell->keyword, "vffa-rect-str-enge") == 0 ||
		strcmp(cell->keyword, "vffa-rect-atan-str") == 0 ||
		strcmp(cell->keyword, "vffa-sect-lin") == 0 || 
		strcmp(cell->keyword, "vffa-sect-he") == 0 ||
		strcmp(cell->keyword, "vffa-sect-enge") == 0) return YES;
	else return NO;
}

//careful no YES and NO used in this function!!
extern int test_cell_spiangle(struct Cell *cell) {
	if( 	strcmp(cell->keyword, "ffag-spi-lin") == 0 || 
			strcmp(cell->keyword, "ffag-spi-he") == 0 ||
			strcmp(cell->keyword, "ffag-spi-enge") == 0 ||
			strcmp(cell->keyword, "ffag-spi-fullenge") == 0) return 1; 
	else if(strcmp(cell->keyword, "ffag-sti-lin") == 0 || 
			strcmp(cell->keyword, "ffag-sti-he") == 0 ||
			strcmp(cell->keyword, "ffag-sti-enge") == 0) return 3;
	else if(strcmp(cell->keyword, "ffag-tilt-lin") == 0 ||
			strcmp(cell->keyword, "ffag-tilt-he") == 0 ||
			strcmp(cell->keyword, "ffag-tilt-enge") == 0) return 4;
	else if(strcmp(cell->keyword, "ffag-tiltpol-lin") == 0 ||
			strcmp(cell->keyword, "ffag-tiltpol-he") == 0 ||
			strcmp(cell->keyword, "ffag-tiltpol-enge") == 0) return 5;
	else if(test_cell_rect_vffa_bend_cell(cell)==TRUE) return 6;
	else return 0;
}


// ************************************************************************************ //
//										maps											//
// ************************************************************************************ //
//load map from mapfile to struct Map *map
// map point taken along a cartesian or cylindrical mesh and expressed in
// - cartesian coordinates (cart), or 
// - cylindrical with angle in radiant (cylrad), or 
// - cylindrical with angle in degrees (cyldeg)
//things have changed write up needed!!
static void load_map(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	char sym_word[10], bcyl_word[10];
	int sym_mid_plane, r_inc,th_inc,z_inc, loop_order,line_order, jump_line, i1,i2,i3,n1,n2,n3,node1,node2,node3, ir,ith,iz, j, nsteps_r,nsteps_th,nsteps_z, nsteps_r_read,nsteps_th_read,nsteps_z_read;
	double unitlength, scale, deltax, rstart,rend,thstart,thend,zstart,zend,step_r,step_th,step_z, r,z,th, coord[4], b[4], x_test, y_test, thmin, angle_tot;
	char mapfilename[FILENAME_L];
	FILE *mapfile = NULL;
	
	if(debug==YES) {
		if(strcmp(latt->cell[i].keyword, "field-cylmap") == 0) printf("  --> map with cylindrical mesh");
		else if(strcmp(latt->cell[i].keyword, "field-cartmap") == 0) printf("  --> map with cartesian mesh");
		else errorstop("strange keyword in load_map.");
	}
	
	fscanf(rfile, "%s", (sym_word));
	if(strcmp(sym_word, "sym") == 0) {
		sym_mid_plane = YES;
		if(debug == YES) printf("  --> symmetric map w/ respect to horizontal mid-plane\n");
	}
	else if(strcmp(sym_word, "nosym") == 0) {
		sym_mid_plane = NO;
		if(debug == YES) printf("\n  --> no symmetry assumed in the map\n");
	}
	else errorstop("keyword neither sym nor nosym, check your input file.");
	
	// 
	fscanf(rfile, "%s", (bcyl_word));
	if(strcmp(bcyl_word, "cart") == 0) {
		if(strcmp(latt->cell[i].keyword, "field-cylmap") == 0) errorstop("Straight cell in cylindrical mesh, check everything is OK in init!!");
		else {if(debug == YES) printf("\n  --> Straight cell expressed in cartesian coordinates\n");}
	}
	else if(strcmp(bcyl_word, "cylrad") == 0) {
		if(strcmp(latt->cell[i].keyword, "field-cartmap") == 0) errorstop("circular cell in cartesian mesh with radiant angle, check everything is OK in init!!");
		else {
			if(debug == YES) {
				printf("\n  --> circular cell in cylindrical mesh, angle in rad\n");
			}
		}
	}
	else if(strcmp(bcyl_word, "cyldeg") == 0) {
		if(debug == YES) {
			if(strcmp(latt->cell[i].keyword, "field-cartmap") == 0) printf("circular cell in cartesian mesh with deg angle, check everything is OK in init!!");
			else printf("\n  --> circular cell in cylindrical mesh, angle in deg\n");
		}
	}
	
	else {errorstop("keyword neither cylrad, cyldeg nor cart, check your input file.");}
	
	fscanf(rfile, "%s", (mapfilename));
	newline(rfile);
	if(debug == YES) printf("  --> file: %s\n", mapfilename);
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == YES) printf("  --> stepsize = %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf", &(angle_tot));
	newline(rfile);
	if(strcmp(bcyl_word, "cyldeg") == 0) angle_tot*=PI/180.; //circ map angle in rad
	if(debug == YES) {
		if(strcmp(bcyl_word, "cyldeg") == 0 || strcmp(bcyl_word, "cylrad") == 0) printf("  --> circular cell total angle %lf [deg]\n", angle_tot*180./PI);
		else printf("  --> straight cell total length %lf [m]\n", angle_tot);
	}
	
	fscanf(rfile, "%lf	%le", &(unitlength),&scale);
	newline(rfile);
	if(debug == YES) {
		printf("  --> unitlength: %lf (1 if m are used in the map file, 100 if cm, 1000 if mm)\n", unitlength);
		printf("  --> scaling factor : B = %lf*B_read_in_the_mapfile\n", scale);
	}
	
	fscanf(rfile, "%i	%i", &loop_order, &line_order);
	newline(rfile);
	if(debug == YES) printf("  --> Loop order: %i, and Line order: %i (1: r or x, 2: th or y, 3: z)\n", loop_order, line_order);
	if(loop_order!=123 && loop_order!=132 && loop_order!=213 && loop_order!=231 && loop_order!=312 && loop_order!=321) errorstop("loop_order has a strange value\n");
	if(line_order!=123 && line_order!=132 && line_order!=213 && line_order!=231 && line_order!=312 && line_order!=321) errorstop("loop_order has a strange value\n");
	
	fscanf(rfile, "%i", &(jump_line));
	newline(rfile);
	if(debug == YES) printf("  --> jumping %i lines...\n", jump_line);
	
	fscanf(rfile, "%lf", &(deltax));
	newline(rfile);
	latt->cell[i].deltar = deltax/(unitlength); // stored value in [m]
	shitf_fwk(&(latt->cell[i].framework), -1.0*latt->cell[i].deltar);
	if(debug == YES) printf("  --> Delta x (origin shift): %le [m] \n\tframework: xc = %lf[m] yc = %lf[m] ae = %lf[deg]\n", latt->cell[i].deltar, latt->cell[i].framework.xc, latt->cell[i].framework.yc, latt->cell[i].framework.ae*180./(PI));
	
	fscanf(rfile, "%lf %lf %i", &(rstart), &(rend), &(nsteps_r)); //r (cylindrical case) or x (cartesian case)
	newline(rfile);
	rstart /= (unitlength); //stored value in [m]
	rend /= (unitlength); //stored value in [m]
	if(debug == YES) printf("\n r:  %lf \tto %lf [m],\tnsteps_r=%i\n", rstart, rend, nsteps_r);
	
	fscanf(rfile, "%lf %lf %i", &(thstart), &(thend), &(nsteps_th)); //th (cylindrical case) or y (cartesian case)
	newline(rfile);
	if(strcmp(latt->cell[i].keyword, "field-cylmap") == 0 && strcmp(bcyl_word, "cyldeg") == 0) { // if cylindrical mesh
		thstart *= PI/180.; //stored value in [rad]
		thend *= PI/180.; //stored value in [rad]
		if(debug == YES) printf(" th:  %lf \tto %lf [deg],\tnsteps_th=%i\n", thstart*180./(PI), thend*180./(PI), nsteps_th);
	}
	else if(strcmp(latt->cell[i].keyword, "field-cartmap") == 0) { //if cartesian mesh 
		thstart /= (unitlength); //stored value in [m]
		thend /= (unitlength); //stored value in [m]
		if(debug == YES) printf(" y:  %lf \tto %lf [m],\tn_step_y=%i\n", thstart, thend, nsteps_th);
	}
	else errorstop("strange keyword\n");
	
	fscanf(rfile, "%lf %lf %i", &(zstart), &(zend), &(nsteps_z));
	newline(rfile);
	zstart /= (unitlength); //stored value in [m]
	zend /= (unitlength); //stored value in [m]
	if(debug == YES) printf(" z:  %lf \t\tto %lf [m],\tnsteps_z=%i\n", zstart, zend, nsteps_z);
	
	//if(MIN(zstart,zend) > 0) errorstop("!!!ERROR in load_map: zmin > 0\n");
	//if(nsteps_r <= 1 || nsteps_th <= 1 || nsteps_z <= 1) errorstop("!!!ERROR in load_cylmap_tosca: number of steps <= 1 !!!");
	step_r = comp_step(rstart, rend, nsteps_r);
	step_th = comp_step(thstart, thend, nsteps_th);
	step_z = comp_step(zstart, zend, nsteps_z);
	
	if(rend>rstart) r_inc = YES;
	else r_inc = NO;
	if(thend>thstart) th_inc = YES; 
	else th_inc = NO;
	if(zend>zstart) z_inc = YES;
	else z_inc = NO;
	
	//Cell collimators set at the limits in r and z of the map
	if(strcmp(latt->cell[i].keyword, "field-cartmap") == 0 && strcmp(bcyl_word, "cyldeg") == 0) {
		latt->cell[i].collim.rmin = MIN(rstart,rend)/cos(angle_tot);
		latt->cell[i].collim.rmax = MAX(rstart,rend)/sin(angle_tot);
		if(latt->cell[i].collim.rmin<0 || latt->cell[i].collim.rmax<0) errorstop("circular map in cartesian mesh with negative radius!");
	}
	else {
		latt->cell[i].collim.rmin = MIN(rstart,rend)-TINYLENGTH;
		latt->cell[i].collim.rmax = MAX(rstart,rend)+TINYLENGTH;
	}
	latt->cell[i].collim.zmin = MIN(zstart,zend)-TINYLENGTH;
	latt->cell[i].collim.zmax = MAX(zstart,zend)+TINYLENGTH;
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	
	//Cell boundaries set at the limits in th of the map
	if(strcmp(bcyl_word, "cyldeg") == 0 || strcmp(bcyl_word, "cylrad") == 0) { // if circular map
		latt->cell[i].boun.thmax = angle_tot;//fabs(thend-thstart); //thmin always=0, thmax = fabs(thend-thstart).
		latt->cell[i].boun.ymax = 0; //disable boundary on y
		if(strcmp(latt->cell[i].keyword, "field-cartmap") == 0) printf("\n\nWARNING!!!!!!\nCircular map in cartesian mesh\n\n");
		else {
			if(fabs(angle_tot-fabs(thend-thstart))>TINYDIMLESS) errorstop("cell total angle different than map total angle");
			thmin = MIN(thend,thstart);
		}
	}
	else { //straight map
		latt->cell[i].boun.thmax = 0; //disable boundary on angle
		latt->cell[i].boun.ymax = angle_tot;//fabs(thend-thstart); //here th = y
		if(fabs(angle_tot-fabs(thend-thstart))>TINYDIMLESS) errorstop("cell total length different than map total length");
	}
	if(debug == YES) printf("  --> Cell boundary:\n\tMax_angle = %lf [deg] \t Max_length = %lf [m] (if 0.0: disabled) \n", latt->cell[i].boun.thmax*180./PI, latt->cell[i].boun.ymax);
	
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//++++++++++++++++++++++++++++ read map file +++++++++++++++++++++++++++++++//
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	
	//open the map file and load the map
	mapfile = fopen(mapfilename,"r");
	if(mapfile == NULL) errorstop("\n!!!ERROR connot open mapfile!!!");
	if(debug == YES)  printf("\n... now loading field map ...\n");
	
	//read Tosca file type header
	fscanf(mapfile,"%i %i %i",&nsteps_r_read, &nsteps_th_read, &nsteps_z_read);
	if(debug == YES) printf("mapfile header...%i	%i	%i\n",nsteps_r_read, nsteps_th_read, nsteps_z_read);
	for(j = 0; j < jump_line; j++) newline(mapfile); //and then skip jump_line lines of header
	if(nsteps_r_read != nsteps_r) errorstop("!!!ERROR in load_map: the field map you are loading does not contain the required number of steps in r");
	if(nsteps_th_read != nsteps_th) errorstop("!!!ERROR in load_map: the field map you are loading does not contain the required number of steps in th");
	if(nsteps_z_read != nsteps_z) errorstop("!!!ERROR in load_map: the field map you are loading does not contain the required number of steps in z");
	
	//allocate memory for the map
	latt->cell[i].map = allocmap(nsteps_r, nsteps_th, nsteps_z);
	//flag sym_mid_plane put in the map structure here:
	latt->cell[i].map.sym = sym_mid_plane;
	
	if(debug == YES) {
		if(strcmp(latt->cell[i].keyword, "field-cylmap") == 0) printf("  --> step-size in r = %.3e [mm], th = %.3e [deg], z = %.3e [mm]\n", step_r*1000, step_th*180./(PI), step_z*1000);
		else if(strcmp(latt->cell[i].keyword, "field-cartmap") == 0) printf("  --> step-size in x = %.3e [mm], y = %.3e [mm], z = %.3e [mm]\n", step_r*1000, step_th*1000, step_z*1000);
	}
	latt->cell[i].map.nnodes[0] = nsteps_r;
	latt->cell[i].map.nnodes[1] = nsteps_th;
	latt->cell[i].map.nnodes[2] = nsteps_z;
	
	latt->cell[i].map.stepsize[0]	= fabs(step_r);
	latt->cell[i].map.stepsize[1]	= fabs(step_th);
	latt->cell[i].map.stepsize[2]	= fabs(step_z);
	
	latt->cell[i].map.mapdim[0]	= MIN(rstart,rend);
	latt->cell[i].map.mapdim[1]	= MAX(rstart,rend);
	latt->cell[i].map.mapdim[2]	= 0.;
	latt->cell[i].map.mapdim[3]	= fabs(thend-thstart);
	latt->cell[i].map.mapdim[4]	= MIN(zstart,zend);
	latt->cell[i].map.mapdim[5]	= MAX(zstart,zend);
	
	if(loop_order/100%10 == 1) n1=nsteps_r;
	else if(loop_order/100%10 == 2) n1=nsteps_th;
	else n1 = nsteps_z;
	if(loop_order/10%10 == 1) n2=nsteps_r;
	else if(loop_order/10%10 == 2) n2=nsteps_th;
	else n2 = nsteps_z;
	if(loop_order%10 == 1) n3=nsteps_r;
	else if(loop_order%10 == 2) n3=nsteps_th;
	else n3 = nsteps_z;
	if(debug==YES) {
		printf("1 = x or r, 2 = y or theta, 3 = z\n");
		printf("loop for order 1: %i\n", loop_order/100%10);
		printf("loop for order 2: %i\n", loop_order/10%10);
		printf("loop for order 3: %i\n", loop_order%10);
	}
	
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
				r = rstart + ir*step_r;
				th = thstart + ith*step_th;
				z = zstart + iz*step_z;
				fscanf(mapfile,"%le %le %le %le %le %le", &(coord[(line_order/100)%10]), &(coord[(line_order/10)%10]), &(coord[line_order%10]), &(b[(line_order/100)%10]), &(b[(line_order/10)%10]), &(b[line_order%10]));
				printf("%le %le %le %le %le %le\n",(coord[(line_order/100)%10]), (coord[(line_order/10)%10]), (coord[line_order%10]), (b[(line_order/100)%10]), (b[(line_order/10)%10]), (b[line_order%10]));
				//printf("r=%lf, th=%lf,z=%lf\n",r,th,z);
				
				////////      check map node coordinates      //////////
				if(strcmp(latt->cell[i].keyword, "field-cylmap") == 0) { // if cylindrical mesh map
					//if(strcmp(bcyl_word, "cart") == 0) {
					//	x_test = unitlength*r*cos(th);
					//	y_test = unitlength*r*sin(th);
					//}
					//else {
					x_test = unitlength*r; //cylindrical coordinates in 
					if(strcmp(bcyl_word, "cylrad") == 0) y_test = th;
					else y_test = th*180./(PI);
					//}
					if(fabs(coord[1] - x_test) > EPS_MAPPOSI*unitlength) { 
						printf("dif=%le, EPS_MAPPOSI=%le\n", coord[1] - x_test, EPS_MAPPOSI*unitlength);
						printf("r=%lf, th=%lf,z=%lf, xtest=%le\n",r,th,z, x_test);
						printf ("\n!!!ERROR in load_map, node X (or r) coordinates does not match (%le vs %le)\n", coord[1], x_test);
						errorstop("Please check your field map file for error, or increase the maximum allowed discrepancy EPS_MAPPOSI (in constant.h)");
					}
					if(fabs(coord[2] - y_test) > EPS_MAPPOSI*unitlength) { 
						printf("dif=%le, EPS_MAPPOSI=%le\n", coord[2] - y_test, EPS_MAPPOSI*unitlength);
						printf("r=%lf, th=%lf,z=%lf, ytest=%le\n",r,th,z, y_test);
						printf ("\n!!!ERROR in load_map, node Y (or theta) coordinates does not match (%le vs %le)\n", coord[2], y_test);
						errorstop("Please check your field map file for error, or increase the maximum allowed discrepancy EPS_MAPPOSI (in constant.h)");
					}
				}
				else { //if cartesian mesh
					if(strcmp(bcyl_word, "cart") == 0) {
						if(fabs(coord[1] - unitlength*r) > EPS_MAPPOSI*unitlength) {
							printf ("\n!!!ERROR in load_map, node X coordinates does not match (%le vs %le)\n", coord[1], unitlength*r);
							printf("x=%lf, y=%lf,z=%lf\n",r,th,z);
							errorstop("Please check your field map file for error, or increase the maximum allowed discrepancy EPS_MAPPOSI (in constant.h)");
						}
					}
					else if(strcmp(bcyl_word, "cylrad") == 0) {
						if(fabs(coord[2] - th) > EPS_MAPPOSI) {
							printf ("\n!!!ERROR in load_map, node Y coordinates does not match (%le vs %le)\n", coord[2], unitlength*th);
							printf("r=%lf, th=%lf,z=%lf\n",r,th,z);
							errorstop("Please check your field map file for error, or increase the maximum allowed discrepancy EPS_MAPPOSI (in constant.h)");
						}
					}
					else if(strcmp(bcyl_word, "cyldeg") == 0) {
						if(fabs(coord[2] - th) > EPS_MAPPOSI) {
							printf ("\n!!!ERROR in load_map, node Y coordinates does not match (%le vs %le)\n", coord[2], unitlength*th);
							printf("r=%lf, th=%lf,z=%lf\n",r,th,z);
							errorstop("Please check your field map file for error, or increase the maximum allowed discrepancy EPS_MAPPOSI (in constant.h)");
						}
					}
					else errorstop("what's goingg on?\n");
				}
				if(fabs(coord[3] - unitlength*z) > EPS_MAPPOSI*unitlength) {
					printf ("\n!!!ERROR in load_map, node Z coordinates does not match (%le vs %le)\n", coord[3], unitlength*z);
					printf("x or r=%lf, y or th=%lf,z=%lf\n",r,th,z);
					errorstop("Please check your field map file for error, or increase the maximum allowed discrepancy EPS_MAPPOSI (in constant.h)");
				}
				//take care of increase or decrease values in the mapfile
				if(r_inc==YES) node1 = ir; //start to fill the nodes from the beginning
				else node1 = nsteps_r-1-ir; //start to fill the nodes from the end
				if(th_inc==YES) node2 = ith;
				else node2 = nsteps_th-1-ith;
				if(z_inc==YES) node3 = iz;
				else node3 = nsteps_z-1-iz;
				
				latt->cell[i].map.node[node1][node2][node3].coord[0] = r; //r (cylindrical case) or x (cartesian case)
				latt->cell[i].map.node[node1][node2][node3].coord[1] = th-thstart; //th (cylindrical case) or y (cartesian case)
				latt->cell[i].map.node[node1][node2][node3].coord[2] = z;
				if(strcmp(latt->cell[i].keyword, "field-cartmap") == 0) {
					latt->cell[i].map.node[node1][node2][node3].b[0]	= scale*b[1];	//bx
					latt->cell[i].map.node[node1][node2][node3].b[1]	= scale*b[2];	//by
				}
				else {
					//latt->cell[i].map.node[node1][node2][node3].b[0]	= scale*b[1];	//bx
					//latt->cell[i].map.node[node1][node2][node3].b[1]	= scale*b[2];	//by
					latt->cell[i].map.node[node1][node2][node3].b[0]	= scale*(b[1]*cos(th-thmin) - b[2]*sin(th-thmin));	//bx=br*cos(th)+bth*sin(th)
					latt->cell[i].map.node[node1][node2][node3].b[1]	= scale*(b[2]*cos(th-thmin) + b[1]*sin(th-thmin));	//by=-br*sin(th)+bth*cos(th)
				}
				latt->cell[i].map.node[node1][node2][node3].b[2]	= scale*b[3];	//bz
				//printf("%i,%i,%i: %lf, %lf, %lf, %le, %le, %le\n", node1,node2,node3, latt->cell[i].map.node[node1][node2][node3].coord[0], latt->cell[i].map.node[node1][node2][node3].coord[1], latt->cell[i].map.node[node1][node2][node3].coord[2],
				//latt->cell[i].map.node[node1][node2][node3].b[0], latt->cell[i].map.node[node1][node2][node3].b[1],latt->cell[i].map.node[node1][node2][node3].b[2]);
			}
		}
	}
	fclose(mapfile);
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == YES) printf("\n");
}

//load map from mapfile to struct Map *map
//used for old input, in particular for tosca maps with cylindrical mesh expressed in cartesian coordinates
// map point taken along a cartesian or cylindrical mesh and expressed in cartesian coordinates
static void load_map_legacy(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt)
{
	char sym_word[10];
	int sym_mid_plane, r_inc,th_inc,z_inc, loop_order,line_order, jump_line, i1,i2,i3,n1,n2,n3,node1,node2,node3, ir,ith,iz, j, nsteps_r,nsteps_th,nsteps_z, nsteps_r_read,nsteps_th_read,nsteps_z_read;
	double unitlength, scale, deltax, rstart,rend,thstart,thend,zstart,zend,step_r,step_th,step_z, r,z,th, coord[4], b[4];
	char mapfilename[FILENAME_L];
	FILE *mapfile = NULL;
	
	if(strcmp(latt->cell[i].keyword, "field-cylmap-old") == 0) sprintf(latt->cell[i].keyword, "field-cylmap");
	else if(strcmp(latt->cell[i].keyword, "field-cartmap-old") == 0) sprintf(latt->cell[i].keyword, "field-cartmap");
	fscanf(rfile, "%s", (sym_word));
	if(strcmp(sym_word, "sym") == 0) {
		sym_mid_plane = YES;
		if(debug == YES) printf("\n  --> symmetric map w/ respect to mid-plane\n");
	}
	else if(strcmp(sym_word, "nosym") == 0) {
		latt->cell[i].map.sym = NO;
		if(debug == YES) printf("\n  --> no symmetry assumed in the map\n");
	}
	else errorstop("keyword neither sym nor nosym, check your input file.");
	
	fscanf(rfile, "%s", (mapfilename));
	newline(rfile);
	if(debug == 1) printf("  --> file: %s\n", mapfilename);
	
	fscanf(rfile, "%lf", &(latt->cell[i].stepsize));
	newline(rfile);
	if(debug == 1) printf("  --> stepsize = %lf [m]\n", latt->cell[i].stepsize);
	
	fscanf(rfile, "%lf	%le", &(unitlength),&scale);
	newline(rfile);
	if(debug == 1) {
		printf("  --> unitlength: %lf (1 if m are used in the map file, 100 if cm, 1000 if mm)\n", unitlength);
		printf("  --> scaling factor : B = %lf*B_read_in_the_mapfile\n", scale);
	}
	
	fscanf(rfile, "%i	%i", &loop_order, &line_order);
	newline(rfile);
	if(debug == 1) printf("  --> Loop order: %i, and Line order: %i (1: r or x, 2: th or y, 3: z)\n", loop_order, line_order);
	if(loop_order!=123 && loop_order!=132 && loop_order!=213 && loop_order!=231 && loop_order!=312 && loop_order!=321) errorstop("loop_order has a strange value\n");
	if(line_order!=123 && line_order!=132 && line_order!=213 && line_order!=231 && line_order!=312 && line_order!=321) errorstop("loop_order has a strange value\n");
	
	fscanf(rfile, "%i", &(jump_line));
	newline(rfile);
	if(debug == 1) printf("  --> jumping %i lines...\n", jump_line);
	
	fscanf(rfile, "%lf", &(deltax));
	newline(rfile);
	latt->cell[i].deltar = deltax/(unitlength); // stored value in [m]
	shitf_fwk(&(latt->cell[i].framework), -1.0*latt->cell[i].deltar);
	if(debug == 1) printf("  --> Delta x (origin shift): %le [m] \n\tframework: xc = %lf[m] yc = %lf[m] ae = %lf[deg]\n", latt->cell[i].deltar, latt->cell[i].framework.xc, latt->cell[i].framework.yc, latt->cell[i].framework.ae*180./(PI));
	
	fscanf(rfile, "%lf %lf %i", &(rstart), &(rend), &(nsteps_r)); //r (cylindrical case) or x (cartesian case)
	newline(rfile);
	rstart /= (unitlength); //stored value in [m]
	rend /= (unitlength); //stored value in [m]
	if(debug == 1) printf("\n r:  %lf \tto %lf [m],\tnsteps_r=%i\n", rstart, rend, nsteps_r);
	
	fscanf(rfile, "%lf %lf %i", &(thstart), &(thend), &(nsteps_th)); //th (cylindrical case) or y (cartesian case)
	newline(rfile);
	if(strcmp(latt->cell[i].keyword, "field-cylmap") == 0) { // if cylindrical map
		thstart *= PI/180.; //stored value in [rad]
		thend *= PI/180.; //stored value in [rad]
		if(debug == 1) printf(" th:  %lf \tto %lf [deg],\tnsteps_th=%i\n", thstart*180./(PI), thend*180./(PI), nsteps_th);
	}
	else if(strcmp(latt->cell[i].keyword, "field-cartmap") == 0) { //if cartesian map 
		thstart /= (unitlength); //stored value in [m]
		thend /= (unitlength); //stored value in [m]
		if(debug == 1) printf(" y:  %lf \tto %lf [m],\tn_step_y=%i\n", thstart, thend, nsteps_th);
	}
	else errorstop("strange keyword\n");
	
	fscanf(rfile, "%lf %lf %i", &(zstart), &(zend), &(nsteps_z));
	newline(rfile);
	zstart /= (unitlength); //stored value in [m]
	zend /= (unitlength); //stored value in [m]
	if(debug == 1) printf(" z:  %lf \t\tto %lf [m],\tnsteps_z=%i\n", zstart, zend, nsteps_z);
	
	if(MIN(zstart,zend) > 0) errorstop("!!!ERROR in load_map: zmin > 0\n");
	if(nsteps_r <= 1 || nsteps_th <= 1 || nsteps_z <= 1) errorstop("!!!ERROR in load_cylmap_tosca: number of steps <= 1 !!!");
	step_r = (rend - rstart)/(nsteps_r - 1.);
	step_th = (thend - thstart)/(nsteps_th - 1.);
	step_z = (zend - zstart)/(nsteps_z - 1.);
	if(rend>rstart) r_inc = YES;
	else r_inc = NO;
	if(thend>thstart) th_inc = YES; 
	else th_inc = NO;
	if(zend>zstart) z_inc = YES;
	else z_inc = NO;
	
	//Cell collimators set at the limits in r and z of the map
	latt->cell[i].collim.rmin = MIN(rstart,rend)-TINYLENGTH;
	latt->cell[i].collim.rmax = MAX(rstart,rend)+TINYLENGTH;
	latt->cell[i].collim.zmax = MAX(zstart,zend)+TINYLENGTH;
	if(sym_mid_plane == YES) latt->cell[i].collim.zmin = -latt->cell[i].collim.zmax;
	else latt->cell[i].collim.zmin = MIN(zstart,zend)-TINYLENGTH;
	if(debug == YES) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmin = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmin, latt->cell[i].collim.zmax);
	//if(debug == 1) printf("  --> collimators:\n\trmin = %lf [m] rmax = %lf [m] zmax = %lf [m]\n", latt->cell[i].collim.rmin, latt->cell[i].collim.rmax, latt->cell[i].collim.zmax);
	
	//Cell boundaries set at the limits in th of the map
	if(strcmp(latt->cell[i].keyword, "field-cylmap") == 0) { // if cylindrical map
		latt->cell[i].boun.thmax = fabs(thend-thstart); //thmin always=0, thmax = fabs(thend-thstart).
		latt->cell[i].boun.ymax = 0; //disable boundary on y
	}
	else {
		latt->cell[i].boun.thmax = 0; //disable boundary on angle
		latt->cell[i].boun.ymax = fabs(thend-thstart); //here th = y
	}
	if(debug == 1) printf("  --> Cell boundary:\n\tMax_angle = %lf [rad] \t Max_length = %lf (if 0.0: disabled) \n", latt->cell[i].boun.thmax, latt->cell[i].boun.ymax);
	
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	//++++++++++++++++++++++++++++ read map file +++++++++++++++++++++++++++++++//
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++//
	
	//open the map file and load the map
	mapfile = fopen(mapfilename,"r");
	if(mapfile == NULL) errorstop("\n!!!ERROR connot open mapfile!!!");
	if(debug == 1)  printf("\n... now loading field map ...\n");
	
	//read Tosca file type header
	fscanf(mapfile,"%i %i %i",&nsteps_r_read, &nsteps_th_read, &nsteps_z_read);
	if(debug == 1) printf("mapfile header...%i	%i	%i\n",nsteps_r_read, nsteps_th_read, nsteps_z_read);
	for(j = 0; j < jump_line; j++) newline(mapfile); //and then skip jump_line lines of header
	if(nsteps_r_read != nsteps_r) errorstop("!!!ERROR in load_map: the field map you are loading does not contain the required number of steps in r");
	if(nsteps_th_read != nsteps_th) errorstop("!!!ERROR in load_map: the field map you are loading does not contain the required number of steps in th");
	if(nsteps_z_read != nsteps_z) errorstop("!!!ERROR in load_map: the field map you are loading does not contain the required number of steps in z");
	
	//allocate memory for the map
	latt->cell[i].map = allocmap(nsteps_r, nsteps_th, nsteps_z);
	//flag sym_mid_plane put in the map structure here:
	latt->cell[i].map.sym = sym_mid_plane;
	
	if(debug == 1) {
		if(strcmp(latt->cell[i].keyword, "field-cylmap") == 0) printf("  --> step-size in r = %.3e [mm], th = %.3e [deg], z = %.3e [mm]\n", step_r*1000, step_th*180./(PI), step_z*1000);
		else if(strcmp(latt->cell[i].keyword, "field-cartmap") == 0) printf("  --> step-size in x = %.3e [mm], y = %.3e [mm], z = %.3e [mm]\n", step_r*1000, step_th*1000, step_z*1000);
	}
	latt->cell[i].map.nnodes[0] = nsteps_r;
	latt->cell[i].map.nnodes[1] = nsteps_th;
	latt->cell[i].map.nnodes[2] = nsteps_z;
	
	latt->cell[i].map.stepsize[0]	= fabs(step_r);
	latt->cell[i].map.stepsize[1]	= fabs(step_th);
	latt->cell[i].map.stepsize[2]	= fabs(step_z);
	
	latt->cell[i].map.mapdim[0]	= MIN(rstart,rend);
	latt->cell[i].map.mapdim[1]	= MAX(rstart,rend);
	latt->cell[i].map.mapdim[2]	= 0.;
	latt->cell[i].map.mapdim[3]	= fabs(thend-thstart);
	latt->cell[i].map.mapdim[4]	= MIN(zstart,zend);
	latt->cell[i].map.mapdim[5]	= MAX(zstart,zend);
	
	if(loop_order/100%10 == 1) n1=nsteps_r;
	else if(loop_order/100%10 == 2) n1=nsteps_th;
	else n1 = nsteps_z;
	if(loop_order/10%10 == 1) n2=nsteps_r;
	else if(loop_order/10%10 == 2) n2=nsteps_th;
	else n2 = nsteps_z;
	if(loop_order%10 == 1) n3=nsteps_r;
	else if(loop_order%10 == 2) n3=nsteps_th;
	else n3 = nsteps_z;
	
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
				r = rstart + ir*step_r;
				th = thstart + ith*step_th;
				z = zstart + iz*step_z;
				fscanf(mapfile,"%le %le %le %le %le %le", &(coord[(line_order/100)%10]), &(coord[(line_order/10)%10]), &(coord[line_order%10]), &(b[(line_order/100)%10]), &(b[(line_order/10)%10]), &(b[line_order%10]));
				//printf("%le %le %le %le %le %le\n",(coord[(line_order/100)%10]), (coord[(line_order/10)%10]), (coord[line_order%10]), (b[(line_order/100)%10]), (b[(line_order/10)%10]), (b[line_order%10]));
				//printf("r=%lf, th=%lf,z=%lf\n",r,th,z);
				
				////////      check map node coordinates      //////////
				if(strcmp(latt->cell[i].keyword, "field-cylmap") == 0) { // if cylindrical mesh map
					if(fabs(coord[1] - unitlength*r*cos(th)) > EPS_MAPPOSI*unitlength) { 
						printf("dif=%le, EPS_MAPPOSI=%le\n", coord[1] - unitlength*r*cos(th), EPS_MAPPOSI*unitlength);
						printf ("\n!!!ERROR in load_map, node X coordinates does not match (%le vs %le)\n", coord[1], unitlength*r*cos(th));
						errorstop("Please check your field map file for error, or increase the maximum allowed discrepancy EPS_MAPPOSI (in constant.h)");
					}
					if(fabs(coord[2] - unitlength*r*sin(th)) > EPS_MAPPOSI*unitlength) { 
						printf("dif=%le, EPS_MAPPOSI=%le\n", coord[2] - unitlength*r*sin(th), EPS_MAPPOSI*unitlength);
						printf ("\n!!!ERROR in load_map, node Y coordinates does not match (%le vs %le)\n", coord[2], unitlength*r*sin(th));
						errorstop("Please check your field map file for error, or increase the maximum allowed discrepancy EPS_MAPPOSI (in constant.h)");
					}
				}
				else { //if cartesian mesh
					if(fabs(coord[1] - unitlength*r) > EPS_MAPPOSI*unitlength) {
						printf ("\n!!!ERROR in load_map, node X coordinates does not match (%le vs %le)\n", coord[1], unitlength*r);
						errorstop("Please check your field map file for error, or increase the maximum allowed discrepancy EPS_MAPPOSI (in constant.h)");
					}
					if(fabs(coord[2] - unitlength*th) > EPS_MAPPOSI*unitlength) {
						printf ("\n!!!ERROR in load_map, node Y coordinates does not match (%le vs %le)\n", coord[2], unitlength*th);
						errorstop("Please check your field map file for error, or increase the maximum allowed discrepancy EPS_MAPPOSI (in constant.h)");
					}
				}
				if(fabs(coord[3] - unitlength*z) > EPS_MAPPOSI*unitlength) {
					printf ("\n!!!ERROR in load_map, node Z coordinates does not match (%le vs %le)\n", coord[3], unitlength*z);
					errorstop("Please check your field map file for error, or increase the maximum allowed discrepancy EPS_MAPPOSI (in constant.h)");
				}
				//take care of increase or decrease values in the mapfile
				if(r_inc==YES) node1 = ir; //start to fill the nodes from the beginning
				else node1 = nsteps_r-1-ir; //start to fill the nodes from the end
				if(th_inc==YES) node2 = ith;
				else node2 = nsteps_th-1-ith;
				if(z_inc==YES) node3 = iz;
				else node3 = nsteps_z-1-iz;
				
				latt->cell[i].map.node[node1][node2][node3].coord[0] = r; //r (cylindrical case) or x (cartesian case)
				latt->cell[i].map.node[node1][node2][node3].coord[1] = th;//thstart+th; //th (cylindrical case) or y (cartesian case)
				latt->cell[i].map.node[node1][node2][node3].coord[2] = z;
				latt->cell[i].map.node[node1][node2][node3].b[0]	= scale*b[1];	//bx
				latt->cell[i].map.node[node1][node2][node3].b[1]	= scale*b[2];	//by
				latt->cell[i].map.node[node1][node2][node3].b[2]	= scale*b[3];	//bz
			}
		}
	}
	fclose(mapfile);
	//modify fwk to be the framework attached to the exit face of this Cell
	find_fwk_exitface(fwk, &(latt->cell[i]));
	if(debug == 1) printf("\n");
}

//test if the cell is a field map (return YES) or not (return NO)
extern int test_cell_map(struct Cell *cell) {
	if(	strcmp(cell->keyword, "field-cylmap") == 0 ||
		strcmp(cell->keyword, "field-cartmap") == 0 ||
		strcmp(cell->keyword, "cylmap-tosca") == 0 ||
		strcmp(cell->keyword, "cartmap-tosca") == 0 ||
		strcmp(cell->keyword, "cartmap-meas") == 0 ||
		strcmp(cell->keyword, "cylmap-fieldtest") == 0 ||
		strcmp(cell->keyword, "cylmap-fieldtest2") == 0) return YES;
	else return NO;
}

// ************************************************************************************ //
//							    manip of lattices										//
// ************************************************************************************ //
// bug to fix if superperiod>1 in reference_latt, framework not updated!
//obsolete, do not use!!
extern void gene_alignerror_latt_old(struct Lattice *reference_latt, struct Lattice *err_latt, long *idum, double rms_shift_error, double rms_twist_error, double r_ref_twist)
//rms_shift_error in [m], rms_twist_error in [rad]
{
	int i, j, nbcell, n1, n2;
	double dir_x, dir_y, dir_z, dir_norm, shift, ax, ay, az, dirtw_x, dirtw_y, dirtw_z, dirtw_norm, angle;
	
	//initialize err_latt
	nbcell = reference_latt->nbcell*reference_latt->periodicity;
	err_latt->periodicity = 1;
	err_latt->nbcell = nbcell;
	err_latt->cellnumber = 0;
	
	//allocate memory for err_latt->cell
	err_latt->cell = alloccell(nbcell);
	for(j = 0; j < reference_latt->periodicity; j++) {
		for(i = 0; i < reference_latt->nbcell; i++) {
			
			err_latt->cell[i+reference_latt->nbcell*j] = reference_latt->cell[i];
			//Warning! check if fwk is passed in the line above!
			if(test_cell_map(&(reference_latt->cell[i]))==YES) errorstop("Error in gene_alignerror_latt: memory allocation in the case of a tracking in a field map with error is not yet implemented");
			//if(strcmp(reference_latt->cell[i].keyword, "cylmap-tosca") == 0) errorstop("Error in gene_alignerror_latt: memory allocation in the case of a tracking in a field map with error is not yet implemented");
			
			else if(strcmp(reference_latt->cell[i].keyword, "rf-thingap") == 0) {
					//allocate memory for the err_latt->cell[i+reference_latt->nbcell*j].cav
					err_latt->cell[i+reference_latt->nbcell*j].cav = alloccav(1);
					//initialize err_latt->cell[i+reference_latt->nbcell*j].cav
					err_latt->cell[i+reference_latt->nbcell*j].cav[0] = reference_latt->cell[i].cav[j];
			}
			
			else if(strcmp(reference_latt->cell[i].keyword, "drift") == 0 || strcmp(reference_latt->cell[i].keyword, "collimator") == 0); //do nothing
			
			else {	
					//allocate memory and copy mpara, efben and efbex
					err_latt->cell[i+reference_latt->nbcell*j].mpara	= allocmatrix(reference_latt->cell[i].nbcomp, SIZE_MPARA);
					for(n1 = 0; n1 < reference_latt->cell[i].nbcomp; n1++) {
						for(n2 = 0; n2 < SIZE_MPARA; n2++) {
							err_latt->cell[i+reference_latt->nbcell*j].mpara[n1][n2] = reference_latt->cell[i].mpara[n1][n2];
						}
					}
					err_latt->cell[i+reference_latt->nbcell*j].efben	= allocmatrix(reference_latt->cell[i].nbcomp, SIZE_EFB);
					for(n1 = 0; n1 < reference_latt->cell[i].nbcomp; n1++) {
						for(n2 = 0; n2 < SIZE_EFB; n2++) {
							err_latt->cell[i+reference_latt->nbcell*j].efben[n1][n2] = reference_latt->cell[i].efben[n1][n2];
						}
					}
					err_latt->cell[i+reference_latt->nbcell*j].efbex	= allocmatrix(reference_latt->cell[i].nbcomp, SIZE_EFB);
					for(n1 = 0; n1 < reference_latt->cell[i].nbcomp; n1++) {
						for(n2 = 0; n2 < SIZE_EFB; n2++) {
							err_latt->cell[i+reference_latt->nbcell*j].efbex[n1][n2] = reference_latt->cell[i].efbex[n1][n2];
						}
					}
					
				//allocate memory for alierror (initialization is done several lines after...)
				err_latt->cell[i+reference_latt->nbcell*j].alierror	= allocmatrix(reference_latt->cell[i].nbcomp, SIZE_ALIERR);
			}
		}
	}
	
	//initialization of err_latt->cell[].alierror with random distribution of the errors
	for(i = 0; i < err_latt->nbcell; i++) {
		if(strcmp(err_latt->cell[i].keyword, "cylmap-tosca") != 0 && strcmp(err_latt->cell[i].keyword, "rf-thingap") != 0  && strcmp(err_latt->cell[i].keyword, "drift") != 0  && strcmp(err_latt->cell[i].keyword, "collimator") != 0) {
			
			//determination of the translation parameters (in the cell's local cartesian framework)
			for(;;) {
				dir_x = ran1(idum);
				dir_y = ran1(idum);
				dir_z = ran1(idum);
				dir_norm = sqrt(dir_x*dir_x + dir_y*dir_y + dir_z*dir_z);
				if(dir_norm <= 1 && dir_norm > 0) break; //if the vector(dir_x, dir_y, dir_z), that will give the direction of the alignment error, is inside the unit sphere (nul vector excluded), break. Else, do it again to find another one. 
			}
			shift = gasdev(idum)*rms_shift_error;
			
			//determine the rotation parameters (in the cell's local cartesian framework)
			if(err_latt->cell[i].boun.ymax != 0) {
				ax = r_ref_twist;
				ay = err_latt->cell[i].boun.ymax/2.;
				az = 0;
			}
			else if(err_latt->cell[i].boun.thmax != 0) {
				ax = r_ref_twist*cos(err_latt->cell[i].boun.thmax/2.);
				ay = r_ref_twist*sin(err_latt->cell[i].boun.thmax/2.);
				az = 0;
			}
			else errorstop("ERROR: unexpected case in gene_alignerror_latt. Give a look to init.c and check");
			
			for(;;) {
				dirtw_x = ran1(idum);
				dirtw_y = ran1(idum);
				dirtw_z = ran1(idum);
				dirtw_norm = sqrt(dirtw_x*dirtw_x + dirtw_y*dirtw_y + dirtw_z*dirtw_z);
				if(dirtw_norm <= 1 && dirtw_norm > 0) break; //if the vector(dirtw_x, dirtw_y, dirtw_z), 
				//that will give the direction of rotation axe, is inside the unit sphere (nul vector excluded), break. Else, do it again to find another one. 
			}
			angle = gasdev(idum)*rms_twist_error;
			
			
			for(n1 = 0; n1 < err_latt->cell[i].nbcomp; n1++) {
				
				err_latt->cell[i].doyou_err = YES;
				
				//translation vector ([m,m,m])
				err_latt->cell[i].alierror[n1][1] = dir_x*shift/dir_norm;
				err_latt->cell[i].alierror[n1][2] = dir_y*shift/dir_norm;
				err_latt->cell[i].alierror[n1][3] = dir_z*shift/dir_norm;
				
				//rotation angle [rad]
				err_latt->cell[i].alierror[n1][4] = angle;
				
				//the rotation axe passes by this point...
				err_latt->cell[i].alierror[n1][5] = ax;
				err_latt->cell[i].alierror[n1][6] = ay;
				err_latt->cell[i].alierror[n1][7] = az;
				//... and is parallel to this unit vector (must be a unit vector!)
				err_latt->cell[i].alierror[n1][8] = dirtw_x/dirtw_norm;
				err_latt->cell[i].alierror[n1][9] = dirtw_y/dirtw_norm;
				err_latt->cell[i].alierror[n1][10] = dirtw_z/dirtw_norm;
				
			}
		}
	}
}

extern void add_error_pole(struct Cell *cell, int npole, double xshift, double yshift, double zshift, 
		double angle_rot, double rot_point_x, double rot_point_y, double rot_point_z, double rot_vect_x, double rot_vect_y, double rot_vect_z)
{
	//translation vector ([m,m,m])
	cell->alierror[npole][1] = xshift;
	cell->alierror[npole][2] = yshift;
	cell->alierror[npole][3] = zshift;
	//rotation angle [rad]
	cell->alierror[npole][4] = angle_rot;
	//the rotation axe passes by this point...
	cell->alierror[npole][5] = rot_point_x;
	cell->alierror[npole][6] = rot_point_y;
	cell->alierror[npole][7] = rot_point_z;
	//... and is parallel to this unit vector (must be a unit vector!)
	cell->alierror[npole][8] =  rot_vect_x;
	cell->alierror[npole][9] =  rot_vect_y;
	cell->alierror[npole][10] = rot_vect_z;
}


extern void comp_error_latt_vffa(struct Lattice *err_latt, long seed, double rms_shift_error, double rms_twist_error)
//rms_shift_error in [m], rms_twist_error in [rad]
{
	int i, n1,n2, doyourandom;
	long idum;
	double dir_x, dir_y, dir_z, dir_norm, shift, shift_x, shift_y, shift_z, ax, ay, az, dirtw_x, dirtw_y, dirtw_z, dirtw_norm, angle, glob_ax, glob_ay, prev_ax[100], prev_ay[100], prev_ax0[100], prev_ay0[100];
	double latt_angle_tot=0, factor_rms_tot=5;
	
	for(i = 0; i < err_latt->nbcell; i++) latt_angle_tot += err_latt->cell[i].boun.thmax;
	printf("error lattice total opening angle: %lf deg\n", latt_angle_tot*180./PI);
	
	//initialise idum (for the random generator)
	if(seed >=0) idum = -seed;
	else idum = seed;
	
	for(i = 0; i < err_latt->nbcell; i++) {
		err_latt->cell[i].doyou_err = YES;
		if(err_latt->cell[i].nbcomp>100) printf("warning in comp_error_latt_vffa, maximum number of components is 100!");
		for(n1 = 0; n1 < err_latt->cell[i].nbcomp; n1++) {
			//rotation around magnet centre
			ax = err_latt->cell[i].mpara[n1][1]*cos(err_latt->cell[i].mpara[n1][0]);
			ay = err_latt->cell[i].mpara[n1][1]*sin(err_latt->cell[i].mpara[n1][0]);
			az = err_latt->cell[i].mpara[n1][3];
			doyourandom = YES;
			if(i!=0) { //find if the pole has been given an error in previous cell, and if so copy it in the current cell
				glob_ax = ax;
				glob_ay = ay;
				fwk_pt_loctoglob(&glob_ax, &glob_ay, &(err_latt->cell[i].framework));
				if(i==err_latt->nbcell-1 && fabs(latt_angle_tot-TWOPI)<TINYDIMLESS) {
					for(n2 = 0; n2 < err_latt->cell[0].nbcomp; n2++) {
						if(fabs(prev_ax0[n2]-glob_ax)<TINYLENGTH && fabs(prev_ay0[n2]-glob_ay)<TINYLENGTH) {
							printf("cell %i, comp %i same than in cell 0 comp %i (last cell case)\n", i, n1, n2);
							doyourandom = NO;
							shift_x = err_latt->cell[0].alierror[n2][1];
							shift_y = err_latt->cell[0].alierror[n2][2];
							fwk_pt_fwktofwk(&shift_x, &shift_y, &(err_latt->cell[0].framework), &(err_latt->cell[i].framework));
							shift_z = err_latt->cell[0].alierror[n2][3];
							angle = err_latt->cell[0].alierror[n2][4];
							dirtw_x = err_latt->cell[0].alierror[n2][8];
							dirtw_y = err_latt->cell[0].alierror[n2][9];
							fwk_vect_fwktofwk(&dirtw_x, &dirtw_y, &(err_latt->cell[0].framework), &(err_latt->cell[i].framework));
							dirtw_z = err_latt->cell[0].alierror[n2][10];
							add_error_pole(&(err_latt->cell[i]), n1, shift_x, shift_y, shift_z, angle, ax, ay, az, dirtw_x, dirtw_y, dirtw_z);
							break;
						}
					}
				}
				for(n2 = 0; n2 < err_latt->cell[i-1].nbcomp; n2++) {
					if(fabs(prev_ax[n2]-glob_ax)<TINYLENGTH && fabs(prev_ay[n2]-glob_ay)<TINYLENGTH) {
						printf("cell %i, comp %i same than in cell %i comp %i (prev cell case)\n", i, n1, i-1, n2);
						doyourandom = NO;
						shift_x = err_latt->cell[i-1].alierror[n2][1];
						shift_y = err_latt->cell[i-1].alierror[n2][2];
						fwk_pt_fwktofwk(&shift_x, &shift_y, &(err_latt->cell[i-1].framework), &(err_latt->cell[i].framework));
						shift_z = err_latt->cell[i-1].alierror[n2][3];
						angle = err_latt->cell[i-1].alierror[n2][4];
						dirtw_x = err_latt->cell[i-1].alierror[n2][8];
						dirtw_y = err_latt->cell[i-1].alierror[n2][9];
						fwk_vect_fwktofwk(&dirtw_x, &dirtw_y, &(err_latt->cell[i-1].framework), &(err_latt->cell[i].framework));
						dirtw_z = err_latt->cell[i-1].alierror[n2][10];
						add_error_pole(&(err_latt->cell[i]), n1, shift_x, shift_y, shift_z, angle, ax, ay, az, dirtw_x, dirtw_y, dirtw_z);
						break;
					}
				}
			}
			if(doyourandom==YES) {
				//determination of the translation parameters (in the cell's local cartesian framework)
				for(;;) {
					dir_x = ran1(&idum);
					dir_y = ran1(&idum);
					dir_z = ran1(&idum);
					dir_norm = sqrt(dir_x*dir_x + dir_y*dir_y + dir_z*dir_z);
					if(dir_norm <= 1 && dir_norm > 0) break; //if the vector(dir_x, dir_y, dir_z), that will give the direction of the alignment error, is inside the unit sphere (nul vector excluded), break. Else, do it again to find another one. 
				}
				if(fabs(rms_shift_error)<TINYLENGTH) shift = 0;
				else {
					for(;;) {
						shift = gasdev(&idum)*rms_shift_error;
						if(fabs(shift) < factor_rms_tot*rms_shift_error) break;
					}
				}
				
				shift_x = dir_x*shift/dir_norm; 
				shift_y = dir_y*shift/dir_norm; 
				shift_z = dir_z*shift/dir_norm;
				
				//determine the rotation vector (in the cell's local cartesian framework)
				for(;;) {
					dirtw_x = ran1(&idum);
					dirtw_y = ran1(&idum);
					dirtw_z = ran1(&idum);
					dirtw_norm = sqrt(dirtw_x*dirtw_x + dirtw_y*dirtw_y + dirtw_z*dirtw_z);
					if(dirtw_norm <= 1 && dirtw_norm > 0) break; //if the vector(dirtw_x, dirtw_y, dirtw_z), that will give the direction of rotation axe, is inside the unit sphere (nul vector excluded), break. Else, do it again to find another one. 
				}
				if(fabs(rms_twist_error)<TINYDIMLESS) angle = 0;
				else {
					for(;;) {
						angle = gasdev(&idum)*rms_twist_error;
						if(fabs(angle) < factor_rms_tot*rms_twist_error) break;
					}
				}
				add_error_pole(&(err_latt->cell[i]), n1, shift_x, shift_y, shift_z, angle, ax, ay, az, dirtw_x/dirtw_norm, dirtw_y/dirtw_norm, dirtw_z/dirtw_norm);
				printf("cell %i, comp %i alierror added\n", i, n1);
			}
		}
		for(n1 = 0; n1 < err_latt->cell[i].nbcomp; n1++) {
			prev_ax[n1] = err_latt->cell[i].mpara[n1][1]*cos(err_latt->cell[i].mpara[n1][0]);
			prev_ay[n1] = err_latt->cell[i].mpara[n1][1]*sin(err_latt->cell[i].mpara[n1][0]);
			fwk_pt_loctoglob(&(prev_ax[n1]), &(prev_ay[n1]), &(err_latt->cell[i].framework));
			if(i==0 && fabs(latt_angle_tot-TWOPI)<TINYDIMLESS) {
				prev_ax0[n1] = prev_ax[n1];
				prev_ay0[n1] = prev_ay[n1];
			}
		}
	}
}

extern void comp_m_error_latt_vffa(struct Lattice *err_latt, long seed, double rms_error)
//rms_error in (delta m)/m
{
	int i, n1,n2, doyourandom;
	long idum;
	double dmovm, ax, ay, glob_ax, glob_ay, prev_ax[100], prev_ay[100], prev_ax0[100], prev_ay0[100];
	double latt_angle_tot=0, factor_rms_tot=5;
	
	for(i = 0; i < err_latt->nbcell; i++) latt_angle_tot += err_latt->cell[i].boun.thmax;
	printf("error lattice total opening angle: %lf deg\n", latt_angle_tot*180./PI);
	
	//initialise idum (for the random generator)
	if(seed >=0) idum = -seed;
	else idum = seed;
	
	for(i = 0; i < err_latt->nbcell; i++) {
		if(err_latt->cell[i].nbcomp>100) printf("warning in comp_m_error_latt_vffa, maximum number of components is 100!");
		for(n1 = 0; n1 < err_latt->cell[i].nbcomp; n1++) {
			//rotation around magnet centre
			ax = err_latt->cell[i].mpara[n1][1]*cos(err_latt->cell[i].mpara[n1][0]);
			ay = err_latt->cell[i].mpara[n1][1]*sin(err_latt->cell[i].mpara[n1][0]);
			doyourandom = YES;
			if(i!=0) { //find if the pole has been given an error in previous cell, and if so copy it in the current cell
				glob_ax = ax;
				glob_ay = ay;
				fwk_pt_loctoglob(&glob_ax, &glob_ay, &(err_latt->cell[i].framework));
				if(i==err_latt->nbcell-1 && fabs(latt_angle_tot-TWOPI)<TINYDIMLESS) {
					for(n2 = 0; n2 < err_latt->cell[0].nbcomp; n2++) {
						if(fabs(prev_ax0[n2]-glob_ax)<TINYLENGTH && fabs(prev_ay0[n2]-glob_ay)<TINYLENGTH) {
							printf("cell %i, comp %i same than in cell 0 comp %i (last cell case)\n", i, n1, n2);
							doyourandom = NO;
							err_latt->cell[i].mpara[n1][4] = err_latt->cell[0].mpara[n2][4];
							break;
						}
					}
				}
				for(n2 = 0; n2 < err_latt->cell[i-1].nbcomp; n2++) {
					if(fabs(prev_ax[n2]-glob_ax)<TINYLENGTH && fabs(prev_ay[n2]-glob_ay)<TINYLENGTH) {
						printf("cell %i, comp %i same than in cell %i comp %i (prev cell case)\n", i, n1, i-1, n2);
						doyourandom = NO;
						err_latt->cell[i].mpara[n1][4] = err_latt->cell[i-1].mpara[n2][4];
						break;
					}
				}
			}
			if(doyourandom==YES) {
				for(;;) {
					dmovm = gasdev(&idum)*rms_error;
					if(fabs(dmovm) < factor_rms_tot*rms_error) break;
				}
				err_latt->cell[i].mpara[n1][4] *= (1+dmovm);
				printf("cell %i, comp %i error in m added\n", i, n1);
			}
		}
		for(n1 = 0; n1 < err_latt->cell[i].nbcomp; n1++) {
			prev_ax[n1] = err_latt->cell[i].mpara[n1][1]*cos(err_latt->cell[i].mpara[n1][0]);
			prev_ay[n1] = err_latt->cell[i].mpara[n1][1]*sin(err_latt->cell[i].mpara[n1][0]);
			fwk_pt_loctoglob(&(prev_ax[n1]), &(prev_ay[n1]), &(err_latt->cell[i].framework));
			if(i==0 && fabs(latt_angle_tot-TWOPI)<TINYDIMLESS) {
				prev_ax0[n1] = prev_ax[n1];
				prev_ay0[n1] = prev_ay[n1];
			}
		}
	}
}

//duplicate cells nb_idcells times from the last charged cell
static void duplicate_cell_fmod(int nb_refcell, struct Framework *fwk, int nb_idcells, struct Lattice *latt)
{
	int m, n, l;
	if(debug == YES) printf(" duplicating cell[%i] %i times\n", nb_refcell, nb_idcells-1);
	for(m = 1; m < nb_idcells; m++) {
		//initialize deltar, nbcomp, and instrutype
		latt->cell[nb_refcell+m].deltar = 0; //default value
		latt->cell[nb_refcell+m].nbcomp = 1; //default value
		latt->cell[nb_refcell+m].instrutype = NO; //default value
		latt->cell[nb_refcell+m].doyou_err = NO; //default value
		
		//set cell.framework
		latt->cell[nb_refcell+m].framework = *fwk;
		//printf("  --> Cell[%i] type: \"%s\"\n",nb_refcell, latt->cell[nb_refcell].keyword);
		strcpy((latt->cell[nb_refcell+m].keyword), (latt->cell[nb_refcell].keyword));
		//printf("  --> Cell[%i] type: \"%s\"\n",nb_refcell+m, latt->cell[nb_refcell+m].keyword);
		if(strcmp(latt->cell[nb_refcell].keyword, "ffag-r-he") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "ffag-r-lin") == 0 || 
		   strcmp(latt->cell[nb_refcell].keyword, "ffag-r-enge") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "ffag-spi-he") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "ffag-spi-lin") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "ffag-spi-enge") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "ffag-spi-fullenge") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "ffag-s-lin") == 0 || 
		   strcmp(latt->cell[nb_refcell].keyword, "ffag-s-he") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "ffag-s-enge") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "ffag-sdl-lin") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "quad-he") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "quad-enge") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "quad-lin") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "vffa-rect-he") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "vffa-rect-enge") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "vffa-rect-enge-add") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "vffa-rect-str-he") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "vffa-rect-str-lin") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "vffa-rect-str-enge") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "vffa-rect-lin") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "vffa-sect-he") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "vffa-sect-enge") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "vffa-sect-lin") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "vffa-rect-atan") == 0 ||
		   strcmp(latt->cell[nb_refcell].keyword, "vffa-rect-atan-str") == 0) {
			latt->cell[nb_refcell+m].nbcomp = latt->cell[nb_refcell].nbcomp;
			//allocate memory
			latt->cell[nb_refcell+m].mpara = allocmatrix(latt->cell[nb_refcell+m].nbcomp, SIZE_MPARA);
			latt->cell[nb_refcell+m].efben = allocmatrix(latt->cell[nb_refcell+m].nbcomp, SIZE_EFB);
			latt->cell[nb_refcell+m].efbex = allocmatrix(latt->cell[nb_refcell+m].nbcomp, SIZE_EFB);
			//allocate memory for errors
			latt->cell[nb_refcell+m].alierror = allocmatrix(latt->cell[nb_refcell+m].nbcomp, SIZE_ALIERR);
			
			
			latt->cell[nb_refcell+m].stepsize = latt->cell[nb_refcell].stepsize;
			latt->cell[nb_refcell+m].deltar = latt->cell[nb_refcell].deltar;
			latt->cell[nb_refcell+m].collim.rmin = latt->cell[nb_refcell].collim.rmin;
			latt->cell[nb_refcell+m].collim.rmax = latt->cell[nb_refcell].collim.rmax;
			latt->cell[nb_refcell+m].collim.zmin = latt->cell[nb_refcell].collim.zmin;
			latt->cell[nb_refcell+m].collim.zmax = latt->cell[nb_refcell].collim.zmax;
			latt->cell[nb_refcell+m].boun.thmax = latt->cell[nb_refcell].boun.thmax;
			latt->cell[nb_refcell+m].boun.ymax = latt->cell[nb_refcell].boun.ymax;
			for(n = 0; n < latt->cell[nb_refcell].nbcomp; n++) {
				for(l=0;l<SIZE_MPARA;l++) latt->cell[nb_refcell+m].mpara[n][l] = latt->cell[nb_refcell].mpara[n][l];
				for(l=0;l<SIZE_EFB;l++) latt->cell[nb_refcell+m].efben[n][l] = latt->cell[nb_refcell].efben[n][l];
				for(l=0;l<SIZE_EFB;l++) latt->cell[nb_refcell+m].efbex[n][l] = latt->cell[nb_refcell].efbex[n][l];
			}
			//modify fwk to be the framework attached to the exit face of this Cell
			find_fwk_exitface(fwk, &(latt->cell[nb_refcell+m]));
			if(debug == YES) printf(" cell[%i] identical as cell[%i] is created. \n", nb_refcell+m, nb_refcell);
		}
		else {
			printf("\nnon implemented keyword for duplication = %s\n", latt->cell[nb_refcell].keyword);
			errorstop("error in duplicate_cell_fmod, try to duplicate a non implemented keyword cell");
		}
	}
}

//generate a lattice that contains only the cell *Cell (!! no new memory allocation for new_latt->cell, if you modify it, take care!!)
extern void gene_latt_1cell(struct Lattice *new_latt, struct Cell *cell) {
	new_latt->nbcell = 1;
	new_latt->periodicity = 1;
	new_latt->cellnumber = 1;
	new_latt->cell = cell;
}

