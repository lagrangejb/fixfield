/*
 *  lin_para.c
 *  ringdesign
 *
 *  Copyright 2010 Kyoto University. All rights reserved.
 *
 */

#include "lin_para.h"

// ************************************************************************************ //
//							phase advances/tunes										//
// ************************************************************************************ //

//tune calculation using Fast Fourier Transform
extern void tune_calculate_fft(struct Particle *reference, double *nux, double *nuz, int tune_turn, double amp_x, double amp_z, 
	struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, char *output_name)
//IMPORTANT!: it is assumed that Particle *reference is on the closed orbit
{
	int i,ixmax, izmax;
	double x_posi[2*(tune_turn+1)]; 
	double x_power=0, x_power_new, x_power_2=0, x_power_new_2, 
		z_posi[2*(tune_turn+1)], z_power=0, z_power_new, z_power_2=0, z_power_new_2,
		e_tot, beta, gamma, xmin=0, xmax=0, ref, nux_2, nuz_2;
	struct Particle test_part;

	printf("number of passes: %i\n", tune_turn);
	if(test_power2(tune_turn)==FALSE) errorstop("in tune_calculate_fft, number of turns is not a power of 2!");

	FILE *wfile;
	if(output_name != NULL) wfile = fopen(output_name, "a");
	
	if(reference->z != 0 || reference->uz != 0 ) printf("!!!ERROR in tune_calculate_fft: unexpected parameter(s)!!!\n z = %le[m], uz = %le[/]\n", reference->z, reference->uz);
	test_part = *reference;
	test_part.hat = -2;
	test_part.status = ALIVE;
	test_part.x += amp_x;
	test_part.z += amp_z;
	
	// tune calculation 
	for(i = 0; i < tune_turn; i++){
		//printf("pass %i \n", i);
		//CLRSCR();
		//printf("pass %i\t", i);
		//fflush(stdout);
		(*transport_part)(&test_part, latt,NULL);
		x_posi[2*i+1] = test_part.x;
		z_posi[2*i+1] = test_part.z;
		x_posi[2*i+2] = 0;
		z_posi[2*i+2] = 0;
		if(i==0) xmin = test_part.x;
		else if(xmin>test_part.x) xmin = test_part.x;
		if(xmax<test_part.x) xmax = test_part.x;
		// output data
		//if(output_debug == YES) fprintf( envfp,"%d	%le	%le	%le	%le\n", i+1, x_posi[2*i+1], test_part.ux, z_posi[2*i+1], test_part.uz);
	}
	ref = (xmin+xmax)/2.;
	//printf("xmin = %lf, xmax = %lf, ref = %lf, co = %lf\n",xmin, xmax, ref,reference->x);
	for(i = 0; i < tune_turn; i++) x_posi[2*i+1] -= ref;
	// call fft 
	four1(x_posi,tune_turn,1);
	four1(z_posi,tune_turn,1);

	//primary peaks
	for(i = 0; i < tune_turn/2; i++) {
		x_power_new = x_posi[2*i+1]*x_posi[2*i+1];
		if(x_power_new > x_power) {
			x_power = x_power_new;
			*nux = ((double) i)/tune_turn;
			ixmax = i;
		}
		z_power_new = z_posi[2*i+1]*z_posi[2*i+1];
		if(z_power_new > z_power){
			z_power = z_power_new;
			*nuz = ((double) i)/tune_turn;
			izmax = i;
		}
		//if(output_debug == YES) fprintf( envfp2,"%d	%le	%le\n", i, x_posi[2*i+1], z_posi[2*i+1]);
	}
	//if(output_debug == YES) fclose(envfp);
	//if(output_debug == YES) fclose(envfp2);
	
	//secondary peaks
	for(i = 1; i < tune_turn/2; i++) {
		if(i>ixmax-4 && i<ixmax+4) continue;
		else {
			x_power_new_2 = x_posi[2*i+1]*x_posi[2*i+1];
			if(x_power_new_2 > x_power_2) {
				x_power_2 = x_power_new_2;
				nux_2 = ((double) i)/tune_turn;
			}
		}
	}
	for(i = 1; i < tune_turn/2; i++) {
		if(i>izmax-4 && i<izmax+4) continue;
		else {
			z_power_new_2 = z_posi[2*i+1]*z_posi[2*i+1];
			if(z_power_new_2 > z_power_2) {
				z_power_2 = z_power_new_2;
				nuz_2 = ((double) i)/tune_turn;
			}
		}
	}
	if(izmax > 0) {
		if((ixmax-1)%(izmax)==0 || (ixmax)%(izmax)==0 || (ixmax+1)%(izmax)==0) {
	//if(fabs(remainder(*nux, *nuz)) < 1.e-2) {
			printf("nux is taken from secondary peak\n");
		//SWAP(*nux,nux_2);
			*nux = nux_2;
	//}
		}	
	}
	if(ixmax > 0) {
		if((izmax-1)%(ixmax)==0 || (izmax)%(ixmax)==0 || (izmax+1)%(ixmax)==0) {
	//if(fabs(remainder(*nuz, *nux)) < 1.e-2) {
			printf("nuz is taken from secondary peak\n");
		//SWAP(*nuz,nuz_2);
			*nuz = nuz_2;
		}
	}
	if(ixmax == 0) *nux = nux_2;
	if(izmax == 0) *nuz = nuz_2;
	
	//}
	/*
	if(z_power_2 > z_power/8. && x_power_2 > x_power/8.) printf("important secondary peaks in both planes\n");
	else if(x_power_2 > x_power/8.) {
		//printf("%lf, %lf\n", fmod(*nux, *nuz), fmod(*nuz, *nux));
		if(fabs(remainder(*nux, *nuz)) < 1.e-2 || fabs(remainder(*nuz, *nux)) < 1.e-2) {
			printf("nux is taken from secondary peak\n");
			//SWAP(*nux,nux_2);
			*nux = nux_2;
		}
	}
	else if(z_power_2 > z_power/8.) {
		if(fabs(remainder(*nux, *nuz)) < 1.e-2 || fabs(remainder(*nuz, *nux)) < 1.e-2) {
			printf("nuz is taken from secondary peak\n");
			//SWAP(*nuz,nuz_2);
			*nuz = nuz_2;
		}
	}//*/
	
	
	if(output_name != NULL || doyouprintf == YES) get_ebg_part(&e_tot, &beta, &gamma, reference);
	//if(output_name != NULL) fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le\n", e_tot, *nux, *nuz, 6.+*nux*latt->periodicity, 4-(*nuz*latt->periodicity), amp_x, amp_z);
	//if(output_name != NULL) fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le\n", reference->brho, *nux, *nuz, *nux*latt->periodicity, *nuz*latt->periodicity, amp_x, amp_z);
	// $6: hffa 3mev normalised emittance ampx^2*gammax*gamma_lorentz*beta_lorentz
	if(output_name != NULL) fprintf(wfile, "%le  %le  %le  %le  %le  %le\n", reference->brho, *nux, *nuz, *nux*latt->periodicity, *nuz*latt->periodicity, amp_x*amp_x);
	//fprintf(wfile, "%le  %le\t%le  %le\n", *nux, *nuz, amp_x, amp_z);
	if(doyouprintf == YES) { 
		//get_ebg_part(&e_tot, &beta, &gamma, reference);
		if(print_color==YES) COLOR("1;35");
		printf("\n-------------------- Tune calculated using FFT (%i passes) -------------------", tune_turn);
		if(print_color==YES) COLOR("0");
		printf("\n");
		printf("@ E(kin) = %le [eV], x_clo = %le [m]\n", (gamma - 1)*reference->m0*CLIGHT2/(UNITCHARGE), reference->x);
		printf("(nux, nuz) : (%lf, %lf)\t", *nux, *nuz);
		printf("(1 - nux,z)   : (%lf, %lf)\n", (1 - *nux), (1 - *nuz));
		printf("(nux_2, nuz_2) : (%lf, %lf)\t", nux_2, nuz_2);
		printf("(1 - nux,z_2)   : (%lf, %lf)\n", (1 - nux_2), (1 - nuz_2));
		printf("max amps (x1,2)   : (%lf, %lf)\n", x_power, x_power_2);
		printf("max amps (z1,2)   : (%lf, %lf)\n", z_power, z_power_2);
		if(print_color==YES) COLOR("1;35");
		printf("--------------------------------------------------------------------------------");
		if(print_color==YES) COLOR("0");
		printf("\n");
	}
	if(output_name != NULL) fclose(wfile);
}

//legacy, use calc_phadv_twiss
extern void tune_calc_matrix_no_period(struct Particle *reference, double *qx, double *qz, double amp_x, double amp_xprime, double amp_z, double amp_zprime, 
	double betax0, double alphax0,  double betaz0, double alphaz0, double *betaxfin, double *alphaxfin, double *betazfin, double *alphazfin, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, char *output_name)
{
	double mh[2][2], mv[2][2], twmh[3][3], v_h[3], twmv[3][3], v_v[3], band_bias, e_tot, beta, gamma,
	//betaxfin, alphaxfin, betazfin, alphazfin,
	halfintx, halfintz, nux, nuz;//intx, intz;
	FILE *wfile;
	
	printf("beam line phase advance calculation:\n");
	if(output_name != NULL) wfile = fopen(output_name, "a");
	
	get_trmatrix_h(mh, &band_bias, reference, amp_x, amp_xprime, latt, *transport_part, NO);
	twissTMtransform(mh, twmh);
	v_h[0] = betax0;
	v_h[1] = alphax0;
	v_h[2] = (1. + alphax0*alphax0)/betax0;
	mvprod3(v_h, twmh, v_h);
	*betaxfin = v_h[0];
	*alphaxfin = v_h[1];
	//printf("betax0 = %lf, alphax0 = %lf, betaxfin = %lf, alphaxfin = %lf\n", betax0, alphax0, *betaxfin, *alphaxfin);
	
	get_trmatrix_v(mv, reference, amp_z, amp_zprime, latt, *transport_part, NO);
	twissTMtransform(mv, twmv);
	v_v[0] = betaz0;
	v_v[1] = alphaz0;
	v_v[2] = (1. + alphaz0*alphaz0)/betaz0;
	mvprod3(v_v, twmv, v_v);
	*betazfin = v_v[0];
	*alphazfin = v_v[1];
	//a revoir !!!!!!!!!!!!!
	nux = atan( ((alphax0 - (*alphaxfin))*mh[0][0] - (*betaxfin)*mh[1][0])/((1 + alphax0*(*alphaxfin))*mh[0][0] + alphax0*(*betaxfin)*mh[1][0]) )/(TWOPI);
	if(nux < 0) nux = 0.5+(nux);
	//nux = asin(mh[1][2]/(sqrt(betax0*betaxfin)))/(TWOPI);

	nuz = atan( ((alphaz0 - (*alphazfin))*mv[0][0] - (*betazfin)*mv[1][0])/((1 + alphaz0*(*alphazfin))*mv[0][0] + alphaz0*(*betazfin)*mv[1][0]) )/(TWOPI);
	if(nuz < 0) nuz = 0.5+nuz;
	//nuz = asin(mv[1][2]/(sqrt(betaz0*betazfin)))/(TWOPI);
	
	//intx = 0;
	//intz = 0;
	//get integer part...
	halfintx = halfintegertunex(reference, latt, amp_xprime);
	halfintz = halfintegertunez(reference, latt, amp_zprime);
	//intx = floor(halfintx);
	//intz = floor(halfintz);
	//printf("intx = %i\n", intx);
	
	//if(mh[1][2] < 0) *qx = intx + (1 - nux);
	//if(mh[1][2] >= 0) *qx = intx + nux;
	//if(mv[1][2] < 0) *qz = intz + (1 - nuz);
	//if(mv[1][2] >= 0) *qz = intz + nuz;
	*qx = halfintx + nux;
	*qz = halfintz +nuz;
	
	*qx = latt->periodicity*(*qx);
	*qz = latt->periodicity*(*qz);
	if(output_name != NULL || doyouprintf == YES) get_ebg_part(&e_tot, &beta, &gamma, reference);
	//if(output_name != NULL) fprintf(wfile, "%le  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n", reference->ux, latt->cell[0].mpara[0][2], latt->cell[0].mpara[1][2], latt->cell[0].mpara[2][2], latt->cell[0].mpara[3][2], latt->cell[0].mpara[0][3], betax0, alphax0, betaz0, alphaz0, *betaxfin, *alphaxfin, *betazfin, *alphazfin, *qx, *qz);
	if(output_name != NULL) fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n", e_tot, *qx/latt->periodicity, *qz/latt->periodicity, *qx, *qz, amp_x, amp_xprime, amp_z, amp_zprime, betax0, alphax0, betaz0, alphaz0, *betaxfin, *alphaxfin, *betazfin, *alphazfin);
	
	if(doyouprintf == YES) {
		if(print_color==YES) COLOR("1;35");
		printf("\n------------------- Get equivalent transfer matrix no period-------------------");
		if(print_color==YES) COLOR("0");
		printf("\n");
		printf("At E(kin) = %le [eV], x_clo = %le [m], x'_clo = %le\n", (gamma - 1)*reference->m0*CLIGHT2/(UNITCHARGE), reference->x, reference->ux);
		printf("MH: %.6f\t%.6f\tMV: %.6f\t%.6f\tdet(MH)=%.6f\n", mh[0][0], mh[0][1], mv[0][0], mv[0][1], mh[0][0]*mh[1][1] - mh[1][0]*mh[0][1]);
		printf("   %.6f\t%.6f\t    %.6f\t%.6f\tdet(MV)=%.6f\n", mh[1][0], mh[1][1], mv[1][0], mv[1][1], mv[0][0]*mv[1][1] - mv[1][0]*mv[0][1]);
		//printf("MH: %lf  %lf\t\tMV: %lf  %lf\tdet(MH) = %lf\n", mh[1][1], mh[1][2], mv[1][1], mv[1][2], mh[1][1]*mh[2][2] - mh[2][1]*mh[1][2]);
		//printf("    %lf  %lf \t   %lf  %lf\tdet(MV) = %lf\n", mh[2][1], mh[2][2], mv[2][1], mv[2][2], mv[1][1]*mv[2][2] - mv[2][1]*mv[1][2]);
		printf("(nux, nuz) : (%lf, %lf)\n", nux, nuz);
		printf("(psix, psiz) : (%lf, %lf) [deg]\n", nux*360., nuz*360.);
		//printf("(1 - nux,z)   : (%lf, %lf)\n", (1 - nux), (1 - nuz));
		printf("(qx, qz)   : (%lf, %lf)  (w/o supersymmetry)\n", *qx/latt->periodicity, *qz/latt->periodicity);
		printf("(qx, qz)   : (%lf, %lf) (w/ supersymmetry)\n", *qx, *qz);
		
		printf("Begin:(beta_x[m], alpha_x)=(%.3f, %.3f),(beta_z[m], alpha_z)=(%.3f, %.3f)\n", betax0, alphax0, betaz0, alphaz0);
		if(nux == nux) printf("End:(beta_x[m], alpha_x)=(%.3f, %.3f), ", *betaxfin, *alphaxfin);
		else printf("!Warning in tune_calc_matrix: cannot compute horizontal twiss parameters because nux = %lf\n", nux);
		if(nuz == nuz) printf("(beta_z[m], alpha_z)=(%.3f, %.3f)\n", *betazfin, *alphazfin);
		else printf("\n!Warning in tune_calc_matrix: cannot compute vertical twiss parameters because nuz = %lf\n", nuz);
		if(print_color==YES) COLOR("1;35");
		printf("-------------------------------------------------------------------------------");
		if(print_color==YES) COLOR("0");
		printf("\n\n");
	}
	if(output_name != NULL) fclose(wfile);
}

//legacy, use calc_tune_twiss, tune calculation from the equivalent transfer matrix
extern void tune_calc_matrix(struct Particle *reference, double *qx, double *qz, double *twiss_bx, double *twiss_ax, double *twiss_bz, double *twiss_az, double amp_x, double amp_xprime, double amp_z, double amp_zprime, 
	struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, char *output_name)
{
	double mh[2][2], mv[2][2], band_bias, e_tot, beta, gamma,
		halfintx, halfintz, intx, intz, nux, nuz;
	//double twiss_gx, twiss_gz;
	FILE *wfile;
	
	printf("periodic tune calculation:\n");
	if(output_name != NULL) wfile = fopen(output_name, "a");
	get_trmatrix_h(mh, &band_bias, reference, amp_x, amp_xprime, latt, *transport_part, NO);
	get_trmatrix_v(mv, reference, amp_z, amp_zprime, latt, *transport_part, NO);
//	if(fabs(mh[1][1] + mh[2][2])<2.+TINYDIMLESS) {
		nux = acos((mh[0][0] + mh[1][1])/2.)/(TWOPI);
		halfintx = halfintegertunex(reference, latt, amp_xprime);
		intx = floor(halfintx);
		if(mh[0][1] >= 0) *qx = intx + nux;
		else if(mh[0][1] < 0) *qx = intx + (1 - nux);
		*qx = latt->periodicity*(*qx);
		*twiss_bx = fabs(mh[0][1]/sin(nux*TWOPI));
		*twiss_ax = (mh[0][0]-mh[1][1])/(2.*sin(*qx*TWOPI));
//	}
//	else *qx = 0;
//	if(fabs(mv[1][1] + mv[2][2])<2.+TINYDIMLESS) {
		nuz = acos((mv[0][0] + mv[1][1])/2.)/(TWOPI);
		halfintz = halfintegertunez(reference, latt, amp_zprime);
		intz = floor(halfintz);
		if(mv[0][1] >= 0) *qz = intz + nuz;
		else if(mv[0][1] < 0) *qz = intz + (1 - nuz);
		*qz = latt->periodicity*(*qz);
		*twiss_bz = fabs(mv[0][1]/sin(nuz*TWOPI));
		*twiss_az = (mv[0][0]-mv[1][1])/(2.*sin(*qz*TWOPI));
//	}
//	else *qz = 0;
	
	if(output_name != NULL || doyouprintf == YES) get_ebg_part(&e_tot, &beta, &gamma, reference);

	if(output_name != NULL) fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n", (gamma - 1)*reference->m0*CLIGHT2/(UNITCHARGE), *qx/latt->periodicity, *qz/latt->periodicity, *qx, *qz, amp_x, amp_xprime, amp_z, amp_zprime, *twiss_bx, *twiss_ax, *twiss_bz, *twiss_az);
	if(doyouprintf == YES) {
		//printf("gamma=%lf, beta=%lf\n", gamma, beta);
		if(print_color==YES) COLOR("1;35");
		printf("\n------------------------ Get equivalent transfer matrix -----------------------");
		if(print_color==YES) COLOR("0");
		printf("\n");
		printf("At E(kin) = %.4e [eV], x_clo = %.4e [m], x'_clo = %.4e [rad]\n", (gamma - 1)*reference->m0*CLIGHT2/(UNITCHARGE), reference->x, reference->ux);
		printf("MH: %.6f\t%.6f\tMV: %.6f\t%.6f\tdet(MH)=%.6f\n", mh[0][0], mh[0][1], mv[0][0], mv[0][1], mh[0][0]*mh[1][1] - mh[1][0]*mh[0][1]);
		printf("   %.6f\t%.6f\t    %.6f\t%.6f\tdet(MV)=%.6f\n", mh[1][0], mh[1][1], mv[1][0], mv[1][1], mv[0][0]*mv[1][1] - mv[1][0]*mv[0][1]);
		printf("(mh11+mh22, mv11+mv22) : (%.7f, %.7f) \n", mh[0][0] + mh[1][1], mv[0][0] + mv[1][1]);
		printf("(nux, nuz) : (%lf, %lf)\n", nux, nuz);
		//printf("(1 - nux,z)   : (%lf, %lf)\n", (1 - nux), (1 - nuz));
		printf("(psix, psiz) : (%lf, %lf) [deg]\n", nux*360., nuz*360.);
		printf("(qx, qz)   : (%lf, %lf)\n", *qx, *qz);
		//printf("muonring-> (qx/m, qz/m)   : (%lf, %lf)\n", *qx/(latt->cell[0].boun.ymax), *qz/(latt->cell[0].boun.ymax));
		
		
		if(nux == nux && sin(nux*TWOPI) != 0) printf("(beta_x[m], alpha_x) = (%.3f, %.3f), ", *twiss_bx, *twiss_ax);
		else printf("!Warning in tune_calc_matrix: cannot compute horizontal twiss parameters because nux = %lf\n", nux);
		
		if(nuz == nuz && sin(nuz*TWOPI) != 0) printf("(beta_z[m], alpha_z) = (%.3f, %.3f)\n", *twiss_bz, *twiss_az);
		else printf("\n!Warning in tune_calc_matrix: cannot compute vertical twiss parameters because nuz = %lf\n", nuz);
	
		if(print_color==YES) COLOR("1;35");
		printf("-------------------------------------------------------------------------------");
		if(print_color==YES) COLOR("0");
		printf("\n");
	}
	
	if(output_name != NULL) fclose(wfile);
}

extern double calc_phase_advance(double m[2][2], double beta0, double alpha0)
{
	return atan_ratio(m[0][1],(m[0][0]*beta0 - m[0][1]*alpha0))/(TWOPI);
}

//periodic tune computation, method 1:
extern double calc_tune1(double m[2][2]) 
{
	double nu;
	nu = acos((m[0][0] + m[1][1])/2.)/(TWOPI);
	if(m[0][1] >= 0) return nu;
	else if(m[0][1] < 0) return (1.0 - nu);
	else return 0;
}

//periodic tune computation, method 2 (better accuracy around 0.5):
extern double calc_tune2(double m[2][2])
{
	double sinnu,cosnu;
	cosnu = (m[0][0] + m[1][1])/2.;
	sinnu = sign(m[0][1])*sqrt(-m[0][1]*m[1][0]-(m[0][0]-m[1][1])*(m[0][0]-m[1][1])/4.);
	printf("sinnu = %le, cosnu = %le\n", sinnu, cosnu);
	return atan_ratio(sinnu,cosnu)/(TWOPI);
}

//return 0.5 if the value of horizontal phase advance exceed 180 deg. (return 1 if exceed 360 deg, and so on...). 
//Warning: does not work if the phase advance of a single cell is larger than 180 deg. If so you must cut your Cell(s) into pieces.
static double halfintegertunex(struct Particle *reference, struct Lattice *latt, double ampxprime)
{
	int i, sgn;
	double halfinteger, dif, xprime_ref, uhp_ref, bias;
	struct Particle part;
	struct Particle partref;
	
	//printf("halfintegertunex\n");
	halfinteger = 0;
	part.s = 0; //part->s reset to 0 at the beginnig of each new lattice
	
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	uhp_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	
	part	= *reference;
	partref = *reference;
	part.hat = -2;		// no acceleration
	partref.hat = -2;	// no acceleration
	part.ux = uhp_ref*sin(xprime_ref + ampxprime);
	part.uy = uhp_ref*cos(xprime_ref + ampxprime);
	
	dif = 0;
	
	if(latt->nbcell <= 0) errorstop("!!!ERROR in halfintegertunex: latt->nbcell <= 0!!!\n");
	for(i = 0; i < latt->nbcell; i++) {
		sgn = sign(dif);
		part_cross_cell(&part, &(latt->cell[i]),NULL);
		part_cross_cell(&partref, &(latt->cell[i]),NULL);
		
		find_x(&dif, &bias, &partref, &part);
		//printf("dif = %le\n", dif);
		if(sign(dif) != sgn && sgn != 0) halfinteger += 0.5;
		//printf("integer = %le\n", halfinteger);
	}
	latt->cellnumber = -1; //set cellnumber to a < 0 value, to be sure that a value will be given explicitely later. (other wise errorstop procedure is called, see above)
	return halfinteger;
}

//return 0.5 if the value of vertical phase advance exceed 180 deg. (return 1 if exceed 360 deg, and so on...). 
//Warning: does not work if the phase advance of a single cell is larger than 180 deg. If so you must cut your Cell(s) into pieces.
static double halfintegertunez(struct Particle *reference, struct Lattice *latt, double ampzprime)
{
	int i, sgn;
	double halfinteger, xprime_ref;
	struct Particle part;
	struct Particle partref;
	
	//printf("halfintegertunez\n");
	halfinteger = 0;
	part.s = 0; //part->s reset to 0 at the beginnig of each new Lattice
	
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	//output_option = 1;
	part	= *reference;
	partref = *reference;
	part.hat = -2;
	partref.hat = -2;
	part.ux = sin(xprime_ref)*cos(ampzprime);
	part.uy = cos(xprime_ref)*cos(ampzprime);
	part.uz = sin(ampzprime);
	part.z = 0;
	//printf("ampzprime = %le, part.z = %le, part.ux = %le, part.uy = %le, part.uz = %le \n", ampzprime, part.z, part.ux, part.uy, part.uz);
	
	if(latt->nbcell <= 0) errorstop("!!!ERROR in halfintegertunez: latt->nbcell <= 0!!!");
	for(i = 0; i < latt->nbcell; i++) {
		sgn = sign(part.z);
		part_cross_cell(&part, &(latt->cell[i]),NULL);
		part_cross_cell(&partref, &(latt->cell[i]),NULL);
		//printf("part.z = %le\n", part.z);
		if (i != 0){
			if(sign(part.z) != sgn && sgn != 0) halfinteger += 0.5;
		}
		//printf("integer = %le\n", halfinteger);
	}
	latt->cellnumber = -1; //set cellnumber to a < 0 value, to be sure that a value will be given explicitely later. (other wise errorstop procedure is called, see above)
	//output_option = 0;
	return halfinteger;
}

extern int calc_phadv_twiss_hor(struct Particle *reference, double *nux, double *qx, double *betaf, double *alphaf, double amp_x, double amp_xprime, double betax0, double alphax0, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf)
{
	int intx;
	double mh[2][2], band_bias,halfintx;
	
	if(get_trmatrix_h(mh, &band_bias, reference, amp_x, amp_xprime, latt, *transport_part, NO) != TRUE) return FALSE;
	halfintx = halfintegertunex(reference, latt, amp_xprime);
	intx = floor(halfintx);
	*nux = calc_phase_advance(mh, betax0, alphax0);
	*qx = latt->periodicity*(*nux+intx);
	calc_twiss(mh, betaf, alphaf, betax0, alphax0);
	if(doyouprintf == YES) {
		printf("MH: %.6f\t%.6f\tdet(MH)=%.6f\n", mh[0][0], mh[0][1], mh[0][0]*mh[1][1] - mh[1][0]*mh[0][1]);
		printf("   %.6f\t%.6f\tnux=%lf, qx=%lf\n", mh[1][0], mh[1][1], *nux, *qx);
	}
	if(*nux!=*nux) return FALSE;
	return TRUE;
}

extern int calc_phadv_twiss_vert(struct Particle *reference, double *nuz, double *qz, double *betaf, double *alphaf, double amp_z, double amp_zprime, double betaz0, double alphaz0, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf)
{
	int intz;
	double mv[2][2],halfintz;
	
	if(get_trmatrix_v(mv, reference, amp_z, amp_zprime, latt, *transport_part, NO) != TRUE) return FALSE;
	halfintz = halfintegertunez(reference, latt, amp_zprime);
	intz = floor(halfintz);
	*nuz = calc_phase_advance(mv, betaz0, alphaz0);
	*qz = latt->periodicity*(*nuz+intz);
	calc_twiss(mv, betaf, alphaf, betaz0, alphaz0);
	if(doyouprintf == YES) {
		printf("MV: %.6f\t%.6f\tdet(MV)=%.6f\n", mv[0][0], mv[0][1], mv[0][0]*mv[1][1] - mv[1][0]*mv[0][1]);
		printf("   %.6f\t%.6f\tnuz=%lf, qz=%lf\n", mv[1][0], mv[1][1], *nuz, *qz);
	}
	if(*nuz!=*nuz) return FALSE;
	return TRUE;
}

extern int calc_phadv_twiss(struct Particle *reference, double *qx, double *qz, double betax0, double alphax0, double betaz0, double alphaz0, double amp_x, double amp_xprime, double amp_z, double amp_zprime, double *betaxf, double *alphaxf, double *betazf, double *alphazf,
	struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, char *outfile)
{
	double e_tot, beta, gamma, nux, nuz;
	FILE *wfile;
	
	if(outfile != NULL || doyouprintf == YES) get_ebg_part(&e_tot, &beta, &gamma, reference);
	if(doyouprintf == YES) {
		if(print_color==YES) COLOR("1;35");
		printf("\n---------- Phase advances and Twiss from equivalent transfer matrix -----------");
		if(print_color==YES) COLOR("0");
		printf("\n");
		printf("At E(kin) = %le [eV], x_clo = %le [m], x'_clo = %le\n", (gamma - 1)*reference->m0*CLIGHT2/(UNITCHARGE), reference->x, reference->ux);
	}
	if(calc_phadv_twiss_hor(reference, &nux, qx, betaxf, alphaxf, amp_x, amp_xprime, betax0, alphax0, latt, *transport_part, doyouprintf) != TRUE) return FALSE;
	if(calc_phadv_twiss_vert(reference, &nuz, qz, betazf, alphazf, amp_z, amp_zprime, betaz0, alphaz0, latt, *transport_part, doyouprintf) != TRUE) return FALSE;
	if(outfile != NULL) {
		wfile = fopen(outfile, "a");
		fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n", e_tot, *qx/latt->periodicity, *qz/latt->periodicity, *qx, *qz, amp_x, amp_xprime, amp_z, amp_zprime, betax0, alphax0, betaz0, alphaz0, *betaxf, *alphaxf, *betazf, *alphazf);
		fclose(wfile);
	}
	if(doyouprintf == YES) {
		printf("(nux, nuz) : (%lf, %lf)\n", nux, nuz);
		printf("(psix, psiz) : (%lf, %lf) [deg]\n", nux*360., nuz*360.);
		printf("(qx, qz)   : (%lf, %lf) (w/ supersymmetry)\n", *qx, *qz);
		
		printf("Begin:(beta_x[m], alpha_x)=(%.3f, %.3f),(beta_z[m], alpha_z)=(%.3f, %.3f)\n", betax0, alphax0, betaz0, alphaz0);
		if(*qx!=0) printf("End:(beta_x[m], alpha_x)=(%.3f, %.3f), ", *betaxf, *alphaxf);
		else printf("Not stable in horizontal\n");
		if(*qz!=0) printf("(beta_z[m], alpha_z)=(%.3f, %.3f)\n", *betazf, *alphazf);
		else printf("Not stable in vertical\n");
		if(print_color==YES) COLOR("1;35");
		printf("-------------------------------------------------------------------------------");
		if(print_color==YES) COLOR("0");
		printf("\n\n");
	}
	return TRUE;
}

extern int calc_tune_twiss_hor(struct Particle *reference, double *nux, double *qx, double *betax, double *alphax, double amp_x, double amp_xprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf)
{
	int intx;
	double mh[2][2], band_bias,halfintx;
	
	if(get_trmatrix_h(mh, &band_bias, reference, amp_x, amp_xprime, latt, *transport_part, NO) != TRUE) return FALSE;
	halfintx = halfintegertunex(reference, latt, amp_xprime);
	intx = floor(halfintx);
	//*nux = calc_tune1(mh);
	*nux = calc_tune2(mh);
	*qx = latt->periodicity*(*nux+intx);
	*betax = fabs(mh[0][1]/sin((*nux)*TWOPI));
	*alphax = (mh[0][0]-mh[1][1])/(2.*sin((*qx)*TWOPI));
	if(doyouprintf == YES) {
		printf("MH: %.6f\t%.6f\tdet(MH)=%.6f\n", mh[0][0], mh[0][1], mh[0][0]*mh[1][1] - mh[1][0]*mh[0][1]);
		printf("   %.6f\t%.6f\tnux=%lf, qx=%lf\n", mh[1][0], mh[1][1], *nux, *qx);
	}
	if(*nux!=*nux) return FALSE;
	return TRUE;
}

extern int calc_tune_twiss_vert(struct Particle *reference, double *nuz, double *qz, double *betaz, double *alphaz, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf)
{
	int intz;
	double mv[2][2],halfintz;
	
	if(get_trmatrix_v(mv, reference, amp_z, amp_zprime, latt, *transport_part, NO) != TRUE) return FALSE;
	halfintz = halfintegertunez(reference, latt, amp_zprime);
	intz = floor(halfintz);
	//*nuz = calc_tune1(mv);
	*nuz = calc_tune2(mv);
	*qz = latt->periodicity*(*nuz+intz);
	*betaz = fabs(mv[0][1]/sin((*nuz)*TWOPI));
	*alphaz = (mv[0][0]-mv[1][1])/(2.*sin((*qz)*TWOPI));
	if(doyouprintf == YES) {
		printf("MV: %.6f\t%.6f\tdet(MV)=%.6f\n", mv[0][0], mv[0][1], mv[0][0]*mv[1][1] - mv[1][0]*mv[0][1]);
		printf("   %.6f\t%.6f\tnuz=%lf, qz=%lf\n", mv[1][0], mv[1][1], *nuz, *qz);
	}
	if(*nuz!=*nuz) return FALSE;
	return TRUE;
}

extern int calc_tune_twiss_4d(struct Particle *reference, double *nux, double *qx, double *nuz, double *qz, double *betax, double *alphax, double *betaz, double *alphaz, 
	double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, int together_decoupled)
{
	int intx, intz, print_matrix=NO;
	double mh[2][2], mv[2][2], band_bias, halfintx, halfintz;
	*nux = 0;
	*nuz = 0;
	
	if(get_matrix_hor_vert(mh, mv, &band_bias, reference, amp_x, amp_xprime, amp_z, amp_zprime, latt, transport_part, print_matrix, together_decoupled) != TRUE) {
		//if(doyouprintf == YES) printf("get_matrix_hor_vert return FALSE\n");
		return FALSE;
	}
	
	if(together_decoupled != 1) {
		halfintx = halfintegertunex(reference, latt, amp_xprime);
		intx = floor(halfintx);
		halfintz = halfintegertunez(reference, latt, amp_zprime);
		intz = floor(halfintz);
	}
	else {
		if(doyouprintf == YES) printf("Warning! uncoupled spaces case: integrales of tune not computed!\n");
		intx = 0;
		intz = 0;
	}
	
	*nux = calc_tune1(mh);
	*nuz = calc_tune1(mv);
	//*nux = calc_tune2(mh);
	//*nuz = calc_tune2(mv);
	*qx = latt->periodicity*(*nux+intx);
	*qz = latt->periodicity*(*nuz+intz);
	calc_periodic_twiss(mh, *nux, betax, alphax);
	//*betax = fabs(mh[0][1]/sin((*nux)*TWOPI));
	//*alphax = (mh[0][0]-mh[1][1])/(2.*sin((*nux)*TWOPI));
	calc_periodic_twiss(mv, *nuz, betaz, alphaz);
	//*betaz = fabs(mv[0][1]/sin((*nuz)*TWOPI));
	//*alphaz = (mv[0][0]-mv[1][1])/(2.*sin((*nuz)*TWOPI));
	if(doyouprintf == YES) {
		printf("MH: %.6f\t%.6f\tdet(MH)=%.6f\n", mh[0][0], mh[0][1], mh[0][0]*mh[1][1] - mh[1][0]*mh[0][1]);
		printf("   %.6f\t%.6f\tnux=%lf, qx=%lf\n", mh[1][0], mh[1][1], *nux, *qx);
		printf("MV: %.6f\t%.6f\tdet(MV)=%.6f\n", mv[0][0], mv[0][1], mv[0][0]*mv[1][1] - mv[1][0]*mv[0][1]);
		printf("   %.6f\t%.6f\tnuz=%lf, qz=%lf\n", mv[1][0], mv[1][1], *nuz, *qz);
	}
	if(*nux!=*nux || *nuz!=*nuz) return FALSE;
	return TRUE;
}

extern int calc_tune_twiss(struct Particle *reference, double *qx, double *qz, double *betax, double *alphax, double *betaz, double *alphaz, double amp_x, double amp_xprime, double amp_z, double amp_zprime,
	struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, char *outfile, int together_decoupled)
{
	double e_tot, beta, gamma, nux, nuz;
	FILE *wfile;
	
	if(outfile != NULL || doyouprintf == YES) get_ebg_part(&e_tot, &beta, &gamma, reference);
	if(doyouprintf == YES) {
		if(print_color==YES) COLOR("1;35");
		printf("\n--------------- Tunes and Twiss from equivalent transfer matrix ---------------");
		if(print_color==YES) COLOR("0");
		printf("\n");
		printf("At E(kin) = %le [eV], x_clo = %le [m], x'_clo = %le\n", (gamma - 1)*reference->m0*CLIGHT2/(UNITCHARGE), reference->x, reference->ux);
	}
	
	if(calc_tune_twiss_4d(reference, &nux, qx, &nuz, qz, betax, alphax, betaz, alphaz, amp_x, amp_xprime, amp_z, amp_zprime, latt, *transport_part, doyouprintf, together_decoupled) != TRUE) return FALSE;

	if(outfile != NULL) {
		wfile = fopen(outfile, "a");
		//fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n", e_tot, *qx/latt->periodicity, *qz/latt->periodicity, *qx, *qz, amp_x, amp_xprime, amp_z, amp_zprime, *betax, *alphax, *betaz, *alphaz, latt->cell[0].mpara[0][6]);
		//fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n", e_tot, *qx/latt->periodicity, *qz/latt->periodicity, *qx, *qz, amp_x, amp_xprime, amp_z, amp_zprime, *betax, *alphax, *betaz, *alphaz);
		fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n", reference->brho, *qx/latt->periodicity, *qz/latt->periodicity, *qx, *qz, amp_x, amp_xprime, amp_z, amp_zprime, *betax, *alphax, *betaz, *alphaz);
		fclose(wfile);
	}
	if(doyouprintf == YES) {
		printf("(nux, nuz) : (%lf, %lf)\n", nux, nuz);
		printf("(1-nux, 1-nuz) : (%lf, %lf)\n", 1.-nux, 1.-nuz);
		//printf("(psix, psiz) : (%lf, %lf) [deg]\n", nux*360., nuz*360.);
		printf("(qx, qz)   : (%lf, %lf) (w/ supersymmetry)\n", *qx, *qz);
		
		if(*qx!=0) printf("Periodic:(beta_x[m], alpha_x)=(%.3f, %.3f), ", *betax, *alphax);
		else printf("Not stable in horizontal\n");
		if(*qz!=0) printf("(beta_z[m], alpha_z)=(%.3f, %.3f)\n", *betaz, *alphaz);
		else printf("Not stable in vertical\n");
		if(print_color==YES) COLOR("1;35");
		printf("-------------------------------------------------------------------------------\n");
		if(print_color==YES) COLOR("0");
		printf("\n");
	}
	return TRUE;
}

// scan phase advance without finding closed orbit !!!!!! LEGACY, TO MODIFY !!!!!!
extern void phase_adv_scan_woco(char *outfilename, struct Lattice *latt, struct Particle *reference_part, double emin_ev, double emax_ev, int nbstep, double amp)
{
	int i;
	double ptot, qx, qz, betax, alphax, betaz, alphaz, estep_ev, ekin_ev;
	struct Particle test_part;
	test_part = *reference_part;
	test_part.hat = -2;
		
	estep_ev = (emax_ev - emin_ev)/nbstep;
	
	for(i = 0; i < nbstep +1; i++) {
		ekin_ev = emax_ev - i*estep_ev;
		
		
		//FODO ONLY!!!!!!
		test_part = *reference_part;
		test_part.hat = -2;
		
		//printf("ekin = %lf\n", ekin_ev);
		ptot = sqrt(pow(ekin_ev*UNITCHARGE/(CLIGHT) + test_part.m0*CLIGHT,2) - pow(test_part.m0*CLIGHT,2));
		printf("\n\nptot = %le [eV/c]\n", ptot*CLIGHT/(UNITCHARGE));
		test_part.brho	= ptot/(test_part.q);
		tune_calc_matrix(&test_part, &qx, &qz, &betax, &alphax, &betaz, &alphaz, amp, amp, amp, amp, latt, part_cross_latt, YES, outfilename);
	}
}

// scan phase advance with finding closed orbit xxp !!!!!! LEGACY, TO MODIFY !!!!!!
extern void phase_adv_scan(char *outfilename, struct Lattice *latt, struct Particle *reference_part, double eps_clo, double emin_ev, double emax_ev, int nbstep, double amp)
{
	int i;//,j;
	double ptot, qx, qz, betax, alphax, betaz, alphaz, estep_ev, ekin_ev;//, m4d[4][4], band_bias, t[4][4], r[4][4], output[6];
	struct Particle test_part;
	//FILE *wfile;
	test_part = *reference_part;
	test_part.hat = -2;
		
	estep_ev = (emax_ev - emin_ev)/nbstep;
	//wfile = fopen(outfilename, "w");
	for(i = 0; i < nbstep +1; i++) {
		ekin_ev = emax_ev - i*estep_ev;
		ptot = sqrt(pow(ekin_ev*UNITCHARGE/(CLIGHT) + test_part.m0*CLIGHT,2) - pow(test_part.m0*CLIGHT,2));
		printf("\n\nptot = %le [eV/c]\n",ptot*CLIGHT/(UNITCHARGE));
		test_part.brho	= ptot/(test_part.q);
		//if(put_on_co_nelmin(latt, &test_part, eps_clo, amp, amp, amp, amp, YES)!=TRUE) errorstop("!!!ERROR in phase_adv_scan, closed orbite not found");
		if(find_closed_orbite_xxp_zzp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), &(test_part.z), &(test_part.uz), eps_clo, latt, YES) != TRUE) errorstop("closed orbit not found\n");
		//calc_tune_eigenvalue(&qx, &qz, amp, amp, amp, amp, &test_part, latt, part_cross_latt, YES, outfilename);
		
	//	for(j=0;j<6;j++) output[j]=0;
	//	get_matrix_firstorder_4d(m4d, &band_bias, &test_part, amp, amp, amp, amp, latt, part_cross_latt, YES);
	//	decouple_matrix_parzen(m4d, t, r, output, NULL);
		//fprintf(wfile, "%le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le  %le\n", ekin_ev, output[0], output[1], output[0], output[1], amp, amp, amp, amp, output[2], output[3], output[4], output[5]);
		calc_tune_twiss(&test_part, &qx, &qz, &betax, &alphax, &betaz, &alphaz, amp, amp, amp, amp, latt, part_cross_latt, YES, outfilename, 1);//flag: 0=2d, 1=4d decoupled, anything else=4d undecoupled.
		//if(find_closed_orbite_xxp(&test_part, &(test_part.x), &(test_part.ux), &(test_part.uy), eps_clo, latt, YES) != TRUE) errorstop("!!!ERROR in phase_adv_scan, closed orbite not found");
		//tune_calc_matrix(&test_part, &qx, &qz, &betax, &alphax, &betaz, &alphaz, amp, amp, amp, amp, latt, part_cross_latt, YES, outfilename);
		//tune_calc_matrix_no_period(&test_part, &qx, &qz, amp, amp, amp, amp, 1.53/2., 0, 1.53/2., 0, latt, part_cross_latt, YES, outfilename);
	}
	//fclose(wfile);
}

// ************************************************************************************ //
//									Twiss - Betafunctions								//
// ************************************************************************************ //

extern double find_s(struct Particle *reference, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*,char*))
{
	struct Particle part;
	
	part = *reference;
	part.hat = -2;
	(*transport_part)(&part, latt, NULL);
	return part.s;
}

extern void calc_twiss(double m[2][2], double *beta1, double *alpha1, double beta0, double alpha0)
{
	double twm[3][3], v[3];
	
	twissTMtransform(m, twm);
	v[0] = beta0;
	v[1] = alpha0;
	v[2] = (1. + alpha0*alpha0)/beta0;
	mvprod3(v, twm, v);
	*beta1 = v[0];
	*alpha1 = v[1];
}

extern void calc_periodic_twiss(double m[2][2], double nu, double *beta, double *alpha)
{
	if(sin(nu*TWOPI) != 0 && nu==nu) {
		*beta = fabs(m[0][1]/sin(nu*TWOPI));
		*alpha = (m[0][0]-m[1][1])/(2.*sin(nu*TWOPI));
	}
	else {
		*beta = 0; *alpha = 0;
		printf("!Warning, periodic Twiss not computed, nu = %le\n", nu);
	}
}

//get the periodic values of the twiss parameters at the entrance of the Lattice
extern void get_periodic_twiss(double *betax0, double *alphax0, double *betaz0, double *alphaz0, 
	double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Particle *reference, struct Lattice *latt, int together_decoupled)
{
	int print_out=NO;
	double mh[2][2], mv[2][2], band_bias, nux, nuz;
	
	get_matrix_hor_vert(mh, mv, &band_bias, reference, amp_x, amp_xprime, amp_z, amp_zprime, latt, part_cross_latt, print_out, together_decoupled);
	if(print_out==YES) printf("in get_periodic_twiss, band_bias = %le\n", band_bias);
	nux = calc_tune1(mh);
	nuz = calc_tune1(mv);
	calc_periodic_twiss(mh, nux, betax0, alphax0);
	calc_periodic_twiss(mv, nuz, betaz0, alphaz0);
}

//get the value of horizontal twiss parameters at an arbitrary point
extern void get_twissx_atinstru(double *betax, double *alphax, double beta1, double alpha1, double amp_x, double amp_xprime, struct Particle *reference, struct Lattice *latt)
{
	double m[2][2], band_bias;
	//double dmu;
	
	if(beta1 == 0) errorstop("!!!ERROR in get_twissx_atinstru: beta1 = 0");
	
	get_trmatrix_h(m, &band_bias, reference, amp_x, amp_xprime, latt, part_cross_latt, NO);
	calc_twiss(m, betax, alphax, beta1, alpha1);
	
	//dmu = atan_ratio(m[1][2], beta1*m[1][1] - m[1][2]*alpha1);
	//*betax = pow(m[1][2]/sin(dmu),2)/beta1;
	//*alphax = (cos(dmu) - m[2][2]*sqrt(*betax/beta1))/(sin(dmu));
	
	//twissTMtransform(m, twm);
	//v_h[0] = beta1;
	//v_h[1] = alpha1;
	//v_h[2] = (1. + alpha1*alpha1)/beta1;
	//mvprod3(v_h, twm, v_h);
	//*betax = v_h[0];
	//*alphax = v_h[1];
}

//get the value of vertical twiss parameters at an arbitrary point
extern void get_twissz_atinstru(double *betaz, double *alphaz, double beta1, double alpha1, double amp_z, double amp_zprime, struct Particle *reference, struct Lattice *latt)
{
	double m[2][2];
	//double dmu;
	
	if(beta1 == 0) errorstop("!!!ERROR in get_twissx_atinstru: beta1 = 0");
	
	get_trmatrix_v(m, reference, amp_z, amp_zprime, latt, part_cross_latt, NO);
	calc_twiss(m, betaz, alphaz, beta1, alpha1);
	
	//dmu = atan_ratio(1., beta1*m[1][1]/m[1][2] - alpha1);
	//*betaz = pow(m[1][2]/sin(dmu),2)/(1.*beta1);
	//*alphaz = (cos(dmu) - m[2][2]*sqrt(*betaz/beta1))/(sin(dmu));
	
	//twissTMtransform(m, twm);
	//v_v[0] = beta1;
	//v_v[1] = alpha1;
	//v_v[2] = (1. + alpha1*alpha1)/beta1;
	//mvprod3(v_v, twm, v_v);
	//*betaz = v_v[0];
	//*alphaz = v_v[1];
}

// flag together_decoupled =0: decoupled, !=0 4d, =1 4d and decoupling from Parzen
extern void get_twiss_atinstru(double *betax, double *alphax, double *betaz, double *alphaz, double betax1, double alphax1, double betaz1, double alphaz1, 
	double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Particle *reference, struct Lattice *latt, int together_decoupled)
{
	int print_out = NO;
	double mh[2][2], mv[2][2], band_bias;
	
	get_matrix_hor_vert(mh, mv, &band_bias, reference, amp_x, amp_xprime, amp_z, amp_zprime, latt, part_cross_latt, print_out, together_decoupled);
	if(print_out==YES) printf("in get_twiss_atinstru, band_bias = %le\n", band_bias);
	calc_twiss(mh, betax, alphax, betax1, alphax1);
	calc_twiss(mv, betaz, alphaz, betaz1, alphaz1);
}

//legacy - use compute_betafunc_latt
extern void compute_betafunc(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime,
	struct Particle *reference, struct Lattice *latt, double betax0, double alphax0, double betaz0, double alphaz0, int together_decoupled)
{
	int i1, i2;
	double betax1, alphax1, betaz1, alphaz1, betax, alphax, betaz, alphaz, s0, s, step_th, step_y;
	//double betxmin=betax0, betxmax=betax0, betzmin=betaz0, betzmax=betaz0;
	struct Particle test_part;
	struct Lattice tempo_latt, latt_1cell;
	FILE *wfile;
	
	//debug test
	if(latt->nbcell <= 0) errorstop("!!!ERROR in compute_betafunc: latt->nbcell <= 0!!!");
	if(nsteps_cell <= 1) errorstop("!!!ERROR in compute_betafunc: nsteps_cell <= 1!!!");
	
	//copy a Lattice, with cell[] stored at a diffent memory address, 
	//!!!and DO NOT FORGET to free the allocated memory befor exiting compute_betafunc!!!
	copy_latt(&tempo_latt, latt);
	
	//init.
//	emptyfile("efb_betafunc.dat");
	if(filename !=NULL) wfile = fopen(filename, "a");
	test_part = *reference;
	test_part.hat = -2;//no acceleration, no output;
	if(debug==YES) printf("compute_betafunc: \n");
	//write first point in wfile
	//if(filename !=NULL) fprintf(wfile, "%le  %le  %le  %le  %le\n", 0.0, betax0, alphax0, betaz0, alphaz0);
	
	//printf("betax0 = %lf, betaz0 = %lf, alphax0 = %lf, alphaz0 = %lf\n", betax0, betaz0, alphax0, alphaz0);
	
	//find twiss parameters at different position in the lattice
	s0 = 0;
	for(i1 = 0; i1 < tempo_latt.nbcell; i1++) {
		
		tempo_latt.cell[i1].instrutype = CUP;
		tempo_latt.cell[i1].instru.ymax = 0;
		tempo_latt.cell[i1].instru.thmax = 0;
		
		step_th = tempo_latt.cell[i1].boun.thmax/(nsteps_cell-1.);
		step_y = tempo_latt.cell[i1].boun.ymax/(nsteps_cell-1.);
		//printf("step_th = %le, step_y = %le\n", step_th, step_y);
		//printf("get_twiss_atinstru:\n");
		//get_twissx_atinstru(&betax,&alphax,betax0,alphax0,amp_x,amp_xprime,reference,&tempo_latt);
		//get_twissz_atinstru(&betaz,&alphaz,betaz0,alphaz0,amp_z,amp_zprime,reference,&tempo_latt);
		get_twiss_atinstru(&betax, &alphax, &betaz, &alphaz, betax0, alphax0, betaz0, alphaz0, amp_x, amp_xprime, amp_z, amp_zprime, reference, &tempo_latt, together_decoupled);
		
		//if(betax<betxmin) betxmin=betax;
		//if(betax>betxmax) betxmax=betax;
		//if(betaz<betzmin) betzmin=betaz;
		//if(betaz>betzmax) betzmax=betaz;
		
		if(filename !=NULL) fprintf(wfile, "%le  %le  %le  %le  %le\n", s0, betax, alphax, betaz, alphaz);
		betax1 = betax;
		alphax1 = alphax;
		betaz1 = betaz;
		alphaz1 = alphaz;
		
		gene_latt_1cell(&latt_1cell, &tempo_latt.cell[i1]);
		
		for(i2 = 0; i2 < nsteps_cell ; i2++) {
			if(latt_1cell.cell[0].boun.thmax != 0) {
				latt_1cell.cell[0].instru.thmax = i2*step_th;
				s = (sqrt(pow(test_part.x,2) + pow(test_part.y,2))+latt_1cell.cell[0].deltar)*i2*step_th;
			}
			else if(latt_1cell.cell[0].boun.ymax != 0) {
				latt_1cell.cell[0].instru.ymax = i2*step_y;
				s = i2*step_y;
			}
			//printf("s=%le\n", s);
			
			//printf("get_twiss_atinstru\n");
			//get_twissx_atinstru(&betax,&alphax,betax1,alphax1,amp_x,amp_xprime,&test_part, &latt_1cell);
			//get_twissz_atinstru(&betaz,&alphaz,betaz1,alphaz1,amp_z,amp_zprime,&test_part, &latt_1cell);
			get_twiss_atinstru(&betax, &alphax, &betaz, &alphaz, betax1, alphax1, betaz1, alphaz1, amp_x, amp_xprime, amp_z, amp_zprime, &test_part, &latt_1cell, together_decoupled);
			
			//if(betax<betxmin) betxmin=betax;
			//if(betax>betxmax) betxmax=betax;
			//if(betaz<betzmin) betzmin=betaz;
			//if(betaz>betzmax) betzmax=betaz;
			
			if(filename !=NULL) fprintf(wfile, "%le  %le  %le  %le  %le\n", s0+s, betax, alphax, betaz, alphaz);
		}
		
		test_part = *reference;
		test_part.hat = -2; //no acceleration, no output;
		part_cross_latt(&test_part, &tempo_latt,NULL); //(only to calculate s from sqrt(pow(test_part.x,2) + pow(test_part.y,2)))
		tempo_latt.cell[i1].instrutype = NO;
		
		//write_efbposition("data/efb_betafunc.dat", s0, test_part.x, &latt_1cell);
		
		s0 += s;
		CLRSCR();
		if(debug==YES) printf("cell number %i/%i (s0 = %.5lf [m])", i1+1, latt->nbcell, s0);
		fflush(stdout);
	}
	if(filename !=NULL) fprintf(wfile,"\n");
	if(filename !=NULL) fclose(wfile);
	if(debug==YES) printf("\n");
	
	//free allocated memory
	free(tempo_latt.cell);
	//printf("result: %le	%le	%le	%le\n", betxmin, betxmax, betzmin, betzmax);
}

//compute beta functions over a cell, from equivalent transfer matrices
extern void compute_betafunc_cell(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime,
	struct Particle *reference, struct Cell *cell, double betax0, double alphax0, double betaz0, double alphaz0, int together_decoupled)
{
	int i;
	double betax, alphax, betaz, alphaz, step_th, step_y, deltar,s;
	struct Particle test_part;
	struct Lattice latt_1cell;
	FILE *wfile;
	
	if(nsteps_cell <= 1) errorstop("!!!ERROR in compute_betafunc: nsteps_cell <= 1!!!");
	
	gene_latt_1cell(&latt_1cell, cell); //copy cell in latt_1cell (no memory allocation)
	test_part = *reference;
	test_part.hat = -2;//no acceleration, no output;
	
	//find twiss parameters at cup position in the cell
	latt_1cell.cell[0].instrutype = CUP;
	latt_1cell.cell[0].instru.ymax = 0;
	latt_1cell.cell[0].instru.thmax = 0;
	deltar = latt_1cell.cell[0].deltar; //store value, set to zero in the function and put it back at the end
	latt_1cell.cell[0].deltar = 0;
	step_th = cell->boun.thmax/(nsteps_cell-1.);
	step_y = cell->boun.ymax/(nsteps_cell-1.);
	
	wfile = fopen(filename, "a");
	fprintf(wfile, "%le  %le  %le  %le  %le\n", test_part.s, betax0, alphax0, betaz0, alphaz0);
	for(i = 1; i < nsteps_cell ; i++) {
		if(latt_1cell.cell[0].boun.thmax != 0) {
			latt_1cell.cell[0].instru.thmax = i*step_th;
			s = test_part.s + (sqrt(pow(test_part.x,2) + pow(test_part.y,2))+latt_1cell.cell[0].deltar)*i*step_th;
		}
		else if(latt_1cell.cell[0].boun.ymax != 0) {
			latt_1cell.cell[0].instru.ymax = i*step_y;
			s = test_part.s + i*step_y;
		}
		//s = find_s(&test_part, &latt_1cell, part_cross_latt);
		//get_twissx_atinstru(&betax, &alphax, betax0, alphax0, amp_x, amp_xprime, &test_part, &latt_1cell);
		//get_twissz_atinstru(&betaz, &alphaz, betaz0, alphaz0, amp_z, amp_zprime, &test_part, &latt_1cell);
		get_twiss_atinstru(&betax, &alphax, &betaz, &alphaz, betax0, alphax0, betaz0, alphaz0, amp_x, amp_xprime, amp_z, amp_zprime, &test_part, &latt_1cell, together_decoupled);
		fprintf(wfile, "%le  %le  %le  %le  %le\n", s, betax, alphax, betaz, alphaz);
	}
	fclose(wfile);
	cell->instrutype = NO;
	cell->deltar = deltar;
}

//compute beta functions over the whole Lattice latt, from equivalent transfer matrices
extern void compute_betafunc_latt(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime,
	struct Particle *reference, struct Lattice *latt, double betax0, double alphax0, double betaz0, double alphaz0, int together_decoupled)
{
	int i;
	double betax1, alphax1, betaz1, alphaz1;
	struct Particle test_part;
	struct Lattice tempo_latt;
	
	if(latt->nbcell <= 0) errorstop("!!!ERROR in compute_betafunc: latt->nbcell <= 0!!!");
	if(nsteps_cell <= 1) errorstop("!!!ERROR in compute_betafunc: nsteps_cell <= 1!!!");
	
	//copy a Lattice, with cell[] stored at a diffent memory address, 
	//!!!and DO NOT FORGET to free the allocated memory befor exiting compute_betafunc!!!
	copy_latt(&tempo_latt, latt);
	
	emptyfile(filename);
	printf("compute_betafunc: \n");
	//printf("betax0 = %lf, betaz0 = %lf, alphax0 = %lf, alphaz0 = %lf\n", betax0, betaz0, alphax0, alphaz0);
	
	//find twiss parameters at different position in the lattice
	//s0 = 0.;
	for(i = 0; i < tempo_latt.nbcell; i++) {
		tempo_latt.cell[i].instrutype = CUP;
		tempo_latt.cell[i].instru.ymax = 0;
		tempo_latt.cell[i].instru.thmax = 0;
		//get_twissx_atinstru(&betax1,&alphax1,betax0,alphax0,amp_x,amp_xprime,reference,&tempo_latt);
		//get_twissz_atinstru(&betaz1,&alphaz1,betaz0,alphaz0,amp_z,amp_zprime,reference,&tempo_latt);
		get_twiss_atinstru(&betax1, &alphax1, &betaz1, &alphaz1, betax0, alphax0, betaz0, alphaz0, amp_x, amp_xprime, amp_z, amp_zprime, reference, &tempo_latt, together_decoupled);
		test_part = *reference;
		test_part.hat = -2; //no acceleration, no output;
		part_cross_latt(&test_part, &tempo_latt, NULL); //to prepare the particle for the next cell (deltar,...)
		//compute_betafunc_cell(filename, nsteps_cell, amp_x, amp_xprime, amp_z, amp_zprime, &test_part, &(tempo_latt.cell[i]), s0, &s, betax1, alphax1, betaz1, alphaz1);
		compute_betafunc_cell(filename, nsteps_cell, amp_x, amp_xprime, amp_z, amp_zprime, &test_part, &(tempo_latt.cell[i]), betax1, alphax1, betaz1, alphaz1, together_decoupled);
		tempo_latt.cell[i].instrutype = NO;
		//s0 += s;
		//CLRSCR();
		//printf("cell number %i/%i (s0 = %.5lf [m])", i+1, latt->nbcell, s0);
		printf("cell number %i/%i (s0 = %.5lf [m])", i+1, latt->nbcell, test_part.s);
		printf("\n");
		//fflush(stdout);
	}
	//free allocated memory
	free(tempo_latt.cell);
}

//compute periodic beta functions over the whole Lattice latt, from equivalent transfer matrices
extern void compute_periodic_betafunc(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime, 
	struct Particle *reference, struct Lattice *latt, int together_decoupled)
{
	double betax0, alphax0, betaz0, alphaz0;
	
	//debug test
	if(latt->nbcell <= 0) errorstop("!!!ERROR in compute_periodic_betafunc: latt->nbcell <= 0!!!");
	if(nsteps_cell <= 1) errorstop("!!!ERROR in compute_periodic_betafunc: nsteps_cell <= 1!!!");
	
	//get periodic beta and alpha
	get_periodic_twiss(&betax0,&alphax0,&betaz0,&alphaz0,amp_x,amp_xprime,amp_z,amp_zprime,reference,latt, together_decoupled);
	
	//get betafunctions starting from periadic twiss parameters
	compute_betafunc(filename, nsteps_cell, amp_x, amp_xprime, amp_z, amp_zprime, reference, latt, betax0, alphax0, betaz0, alphaz0, together_decoupled);
	//compute_betafunc_latt(filename, nsteps_cell, amp_x, amp_xprime, amp_z, amp_zprime, reference, latt, betax0, alphax0, betaz0, alphaz0,together_decoupled);
}

/*
static void write_efbposition(char *filename, double s0, double x_clo, struct Lattice *latt)
{
	int i, j;
	double efb_ymax, x, y;
	FILE *wfile;
	
	efb_ymax = 1;
	
		
	wfile = fopen(filename, "a");
		
	for(i = 0; i < latt->nbcell; i++) {
		for(j = 0; j < latt->cell[i].nbcomp; j++) {
			if(strcmp(latt->cell[i].keyword, "ffag-r-he") == 0 || 
				strcmp(latt->cell[i].keyword, "ffag-r-lin") == 0 || 
				strcmp(latt->cell[i].keyword, "ffag-r-enge") == 0) {
				x = s0 + x_clo*(latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0]);
				y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				
			//	x = s0 + x_clo*(latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0]);
				y = efb_ymax;
				fprintf(wfile, "%le  %le\n",x,y);
				
				x = s0 + x_clo*(latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0]);
			//	y = efb_ymax;
				fprintf(wfile, "%le  %le\n",x,y);
				
			//	x = s0 + x_clo*(latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0]);
				y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				
				x = s0 + x_clo*(latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0]);
			//	y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				
				fprintf(wfile, "\n");
			}
			else if(strcmp(latt->cell[i].keyword, "ffag-s-he") == 0 || 
					strcmp(latt->cell[i].keyword, "ffag-s-lin") == 0 || 
					strcmp(latt->cell[i].keyword, "ffag-s-enge") == 0 ||
					strcmp(latt->cell[i].keyword, "ffag-sdl-lin") == 0) {
				x = s0 + (latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0]);
				y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				
			//	x = s0 + (latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0]);
				y = efb_ymax;
				fprintf(wfile, "%le  %le\n",x,y);
				
				x = s0 + (latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0]);
			//	y = efb_ymax;
				fprintf(wfile, "%le  %le\n",x,y);
				
			//	x = s0 + (latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0]);
				y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				
				x = s0 + (latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0]);
			//	y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				
				fprintf(wfile, "\n");
			}
		}
	}
		
	fclose(wfile);
}//*/

// ************************************************************************************ //
//									dispersion											//
// ************************************************************************************ //

extern void compute_periodic_dispersion_2d(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_dpovp, struct Particle *reference, struct Lattice *latt, int hor)
{
	double disp0, dispprime0;
	
	//debug test
	if(latt->nbcell <= 0) errorstop("!!!ERROR in compute_periodic_dispersion_2d: latt->nbcell <= 0!!!");
	if(nsteps_cell <= 1) errorstop("!!!ERROR in compute_periodic_dispersion_2d: nsteps_cell <= 1!!!");
	
	//get periodic dispersion value and its derivee
	get_periodic_dispersion_2d(&disp0, &dispprime0, amp_x, amp_xprime, amp_dpovp, reference, latt, hor);	
	
	//get betafunctions starting from periadic twiss parameters
	compute_dispersion_2d(filename, nsteps_cell, amp_x, amp_xprime, amp_dpovp, reference, latt, disp0, dispprime0, hor);
}

extern void compute_dispersion_2d(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_dpovp, struct Particle *reference, struct Lattice *latt, double disp0, double dispprime0, int hor)
{
	int i1, i2;
	double disp1, dispprime1, s0, s, step_th, step_y, band_bias;
	double mh[3][3], vect_disp[3];
	struct Particle test_part;
	struct Lattice tempo_latt, latt_1cell;
	FILE *wfile;
	
	//emptyfile("data/efb_dispersion.dat");
	//debug test
	if(latt->nbcell <= 0) errorstop("!!!ERROR in compute_dispersion_2d: latt->nbcell <= 0!!!");
	if(nsteps_cell <= 1) errorstop("!!!ERROR in compute_dispersion_2d: nsteps_cell <= 1!!!");
	
	//copy a Lattice, with cell[] stored at a diffent memory address, 
	//!!!and DO NOT FORGET to free the allocated memory before exiting part_oneturn!!!
	copy_latt(&tempo_latt, latt);
	
	//init.
	wfile = fopen(filename, "w");
	test_part = *reference;
	test_part.hat = -2;//no acceleration, no output;
	printf("compute_dispersion: \n");
	//write first point in wfile
	fprintf(wfile, "%lf  %lf  %lf\n", 0.0, disp0, dispprime0);
	
	//find dispersion at different positions in the lattice
	s0 = 0;
	for(i1 = 0; i1 < tempo_latt.nbcell; i1++) {
		tempo_latt.cell[i1].instrutype = CUP;
		tempo_latt.cell[i1].instru.ymax = 0;
		tempo_latt.cell[i1].instru.thmax = 0;
		
		step_th = tempo_latt.cell[i1].boun.thmax/(nsteps_cell-1.);
		step_y = tempo_latt.cell[i1].boun.ymax/(nsteps_cell-1.);
		
		if(hor==YES) get_trmatrix_h3x3(mh, &band_bias, reference, amp_x, amp_xprime, amp_dpovp,  &tempo_latt, part_cross_latt, NO);
		else if(hor==NO) get_trmatrix_v3x3(mh, &band_bias, reference, amp_x, amp_xprime, amp_dpovp,  &tempo_latt, part_cross_latt, NO);

		vect_disp[0] = disp0;
		vect_disp[1] = dispprime0;
		vect_disp[2] = 1;
		
		mvprod3(vect_disp, mh, vect_disp);
		
		fprintf(wfile, "%le  %le  %le\n", s0, vect_disp[0], vect_disp[1]);
		disp1 = vect_disp[0];
		dispprime1 = vect_disp[1];
		
		gene_latt_1cell(&latt_1cell, &tempo_latt.cell[i1]);
		
		
		for(i2 = 0; i2 < nsteps_cell ; i2++) {
			
			if(latt_1cell.cell[0].boun.thmax != 0) {
				latt_1cell.cell[0].instru.thmax = i2*step_th;
				s = (sqrt(pow(test_part.x,2) + pow(test_part.y,2))+latt_1cell.cell[0].deltar)*i2*step_th;
			}
			else if(latt_1cell.cell[0].boun.ymax != 0) {
				latt_1cell.cell[0].instru.ymax = i2*step_y;
				s = i2*step_y;
			}
			
			//if(hor==YES) get_trmatrix_h3x3(mh, &band_bias, reference, amp_x, amp_xprime, amp_dpovp, &tempo_latt, part_cross_latt, NO);
			//else if(hor==NO) get_trmatrix_v3x3(mh, &band_bias, reference, amp_x, amp_xprime, amp_dpovp, &tempo_latt, part_cross_latt, NO);
			if(hor==YES) get_trmatrix_h3x3(mh, &band_bias, &test_part, amp_x, amp_xprime, amp_dpovp, &latt_1cell, part_cross_latt, NO);
			else if(hor==NO) get_trmatrix_v3x3(mh, &band_bias, &test_part, amp_x, amp_xprime, amp_dpovp, &latt_1cell, part_cross_latt, NO);

			vect_disp[0] = disp1;
			vect_disp[1] = dispprime1;
			vect_disp[2] = 1;
			mvprod3(vect_disp, mh, vect_disp);
			
			fprintf(wfile, "%le  %le  %le\n", s0+s, vect_disp[0], vect_disp[1]);
		}
		
		test_part = *reference;
		test_part.hat = -2; //no acceleration, no output;
		part_cross_latt(&test_part, &tempo_latt,NULL); //(only to calculate s from sqrt(pow(test_part.x,2) + pow(test_part.y,2)))
		tempo_latt.cell[i1].instrutype = NO;
		
		
		//write_efbposition("data/efb_dispersion.dat", s0, test_part.x, &latt_1cell);
		s0 += s;
		CLRSCR();
		printf("cell number %i/%i (s0 = %.5lf [m])", i1+1, latt->nbcell, s0);
		fflush(stdout);
	}
	fprintf(wfile,"\n");
	fclose(wfile);
	printf("\n");
	//free allocated memory
	free(tempo_latt.cell);
}

extern void get_periodic_dispersion_2d(double *disp0, double *dispprime0, double amp_x, double amp_xprime, double amp_dpovp, struct Particle *reference, struct Lattice *latt, int hor)
{
	double mh[3][3], band_bias, det;
	
	if(hor==YES) get_trmatrix_h3x3(mh, &band_bias, reference, amp_x, amp_xprime, amp_dpovp, latt, part_cross_latt, NO);
	else if(hor==NO) get_trmatrix_v3x3(mh, &band_bias, reference, amp_x, amp_xprime, amp_dpovp, latt, part_cross_latt, NO);
	det = mh[0][0]*mh[1][1]-mh[0][1]*mh[1][0];
	*disp0 = ((1 - mh[1][1])*mh[0][2] + mh[0][1]*mh[1][2])/(1 + det - mh[0][0] - mh[1][1]);
	*dispprime0 = ((1 - mh[0][0])*mh[1][2] + mh[1][0]*mh[0][2])/(1 + det - mh[0][0] - mh[1][1]);
}

extern void compute_periodic_dispersion_4d(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime, double amp_dpovp, struct Particle *reference, struct Lattice *latt)
{
	double dispx0, dispprimex0, dispz0, dispprimez0;
	
	//debug test
	if(latt->nbcell <= 0) errorstop("!!!ERROR in compute_periodic_dispersion_4d: latt->nbcell <= 0!!!");
	if(nsteps_cell <= 1) errorstop("!!!ERROR in compute_periodic_dispersion_4d: nsteps_cell <= 1!!!");
	
	get_periodic_dispersion_4d(&dispx0, &dispprimex0, &dispz0, &dispprimez0, amp_x, amp_xprime, amp_z, amp_zprime, amp_dpovp, reference, latt);
	compute_dispersion_4d(filename, nsteps_cell, amp_x, amp_xprime, amp_z, amp_zprime, amp_dpovp, reference, latt, dispx0, dispprimex0, dispz0, dispprimez0);
}

extern void compute_dispersion_4d(char *filename, int nsteps_cell, double amp_x, double amp_xprime, double amp_z, double amp_zprime, double amp_dpovp, struct Particle *reference, struct Lattice *latt, double dispx0, double dispprimex0, double dispz0, double dispprimez0)
{
	int i1, i2;
	double dispx1, dispprimex1, dispz1, dispprimez1, s0, s, step_th, step_y, band_bias;
	double m5d[5][5], vect_disp[5];
	struct Particle test_part;
	struct Lattice tempo_latt, latt_1cell;
	FILE *wfile;
	
	//emptyfile("data/efb_dispersion.dat");
	//debug test
	if(latt->nbcell <= 0) errorstop("!!!ERROR in compute_dispersion_4d: latt->nbcell <= 0!!!");
	if(nsteps_cell <= 1) errorstop("!!!ERROR in compute_dispersion_4d: nsteps_cell <= 1!!!");
	
	//copy a Lattice, with cell[] stored at a diffent memory address, 
	//!!!and DO NOT FORGET to free the allocated memory before exiting part_oneturn!!!
	copy_latt(&tempo_latt, latt);
	
	//init.
	wfile = fopen(filename, "w");
	test_part = *reference;
	test_part.hat = -2;//no acceleration, no output;
	printf("compute_dispersion: \n");
	//write first point in wfile
	fprintf(wfile, "%lf  %lf  %lf  %lf  %lf\n", 0.0, dispx0, dispprimex0, dispz0, dispprimez0);
	
	//find dispersion at different positions in the lattice
	s0 = 0;
	for(i1 = 0; i1 < tempo_latt.nbcell; i1++) {
		tempo_latt.cell[i1].instrutype = CUP;
		tempo_latt.cell[i1].instru.ymax = 0;
		tempo_latt.cell[i1].instru.thmax = 0;
		
		step_th = tempo_latt.cell[i1].boun.thmax/(nsteps_cell-1.);
		step_y = tempo_latt.cell[i1].boun.ymax/(nsteps_cell-1.);
		
		get_matrix_firstorder_5d(m5d, &band_bias, reference, amp_x, amp_xprime, amp_z, amp_zprime, amp_dpovp, &tempo_latt, part_cross_latt, NO);
		
		vect_disp[0] = dispx0;
		vect_disp[1] = dispprimex0;
		vect_disp[2] = dispz0;
		vect_disp[3] = dispprimez0;
		vect_disp[4] = 1;
		
		mvprod5(vect_disp, m5d, vect_disp);
		
		fprintf(wfile, "%le  %le  %le  %le  %le\n", s0, vect_disp[0], vect_disp[1], vect_disp[2], vect_disp[3]);
		dispx1 = vect_disp[0];
		dispprimex1 = vect_disp[1];
		dispz1 = vect_disp[2];
		dispprimez1 = vect_disp[3];
		
		gene_latt_1cell(&latt_1cell, &tempo_latt.cell[i1]);
		
		
		for(i2 = 0; i2 < nsteps_cell ; i2++) {
			
			if(latt_1cell.cell[0].boun.thmax != 0) {
				latt_1cell.cell[0].instru.thmax = i2*step_th;
				s = (sqrt(pow(test_part.x,2) + pow(test_part.y,2))+latt_1cell.cell[0].deltar)*i2*step_th;
			}
			else if(latt_1cell.cell[0].boun.ymax != 0) {
				latt_1cell.cell[0].instru.ymax = i2*step_y;
				s = i2*step_y;
			}
			
			get_matrix_firstorder_5d(m5d, &band_bias, &test_part, amp_x, amp_xprime, amp_z, amp_zprime, amp_dpovp, &latt_1cell, part_cross_latt, NO);
			
			vect_disp[0] = dispx1;
			vect_disp[1] = dispprimex1;
			vect_disp[2] = dispz1;
			vect_disp[3] = dispprimez1;
			vect_disp[4] = 1;
			mvprod5(vect_disp, m5d, vect_disp);
			
			fprintf(wfile, "%le  %le  %le  %le  %le\n", s0+s, vect_disp[0], vect_disp[1], vect_disp[2], vect_disp[3]);
		}
		
		test_part = *reference;
		test_part.hat = -2; //no acceleration, no output;
		part_cross_latt(&test_part, &tempo_latt,NULL); //(only to calculate s from sqrt(pow(test_part.x,2) + pow(test_part.y,2)))
		tempo_latt.cell[i1].instrutype = NO;
		
		
		//write_efbposition("data/efb_dispersion.dat", s0, test_part.x, &latt_1cell);
		s0 += s;
		CLRSCR();
		printf("cell number %i/%i (s0 = %.5lf [m])", i1+1, latt->nbcell, s0);
		fflush(stdout);
	}
	fprintf(wfile,"\n");
	fclose(wfile);
	printf("\n");
	//free allocated memory
	free(tempo_latt.cell);
}

extern void get_periodic_dispersion_4d(double *dispx0, double *dispprimex0, double *dispz0, double *dispprimez0, double amp_x, double amp_xprime, double amp_z, double amp_zprime, double amp_dpovp, struct Particle *reference, struct Lattice *latt)
{
	double m5d[5][5], band_bias;//, test_vect[5];
	//printf("in get_periodic_dispersion:\n");
	
	get_matrix_firstorder_5d(m5d, &band_bias, reference, amp_x, amp_xprime, amp_z, amp_zprime, amp_dpovp, latt, part_cross_latt, NO);
	
	*dispx0 = -((m5d[0][3]*m5d[0][3]*((m5d[0][3]*(-1 + m5d[1][1] + m5d[1][2]*m5d[2][1] + m5d[2][2] - m5d[1][1]*m5d[2][2]) - 
	        m5d[0][2]*(m5d[1][3]*m5d[2][1] + m5d[2][2] - m5d[1][1]*m5d[2][2]) + m5d[0][1]*(m5d[1][3]*(-1 + m5d[2][2]) - m5d[1][2]*m5d[2][2]))*
	      (m5d[3][4]*(m5d[0][3]*m5d[1][2] - m5d[0][2]*m5d[1][3]) - m5d[1][4]*(m5d[0][2] + m5d[0][3]*m5d[3][2] - m5d[0][2]*m5d[3][3]) + 
	        m5d[0][4]*(m5d[1][2] + m5d[1][3]*m5d[3][2] - m5d[1][2]*m5d[3][3])) - 
	     (m5d[2][4]*m5d[0][3]*m5d[1][2] - m5d[0][4]*m5d[1][3] - m5d[2][4]*m5d[0][2]*m5d[1][3] + m5d[0][4]*m5d[1][3]*m5d[2][2] - m5d[0][4]*m5d[1][2]*m5d[2][2] + 
	        m5d[1][4]*(m5d[0][3] - m5d[0][3]*m5d[2][2] + m5d[0][2]*m5d[2][2]))*
	      (m5d[0][3]*(m5d[1][2]*m5d[3][1] + m5d[3][2] - m5d[1][1]*m5d[3][2]) - m5d[0][2]*(-1 + m5d[1][1] + m5d[1][3]*m5d[3][1] + m5d[3][3] - m5d[1][1]*m5d[3][3]) + 
	        m5d[0][1]*(m5d[1][2] + m5d[1][3]*m5d[3][2] - m5d[1][2]*m5d[3][3]))))/
	 (m5d[0][3]*(m5d[0][3]*(-1 + m5d[1][1] + m5d[1][2]*m5d[2][1] + m5d[2][2] - m5d[1][1]*m5d[2][2]) - m5d[0][2]*(m5d[1][3]*m5d[2][1] + m5d[2][2] - m5d[1][1]*m5d[2][2]) + 
	      m5d[0][1]*(m5d[1][3]*(-1 + m5d[2][2]) - m5d[1][2]*m5d[2][2]))*
	    ((m5d[0][3]*m5d[1][2] - m5d[0][2]*m5d[1][3])*(-1 + m5d[0][0] + m5d[0][3]*m5d[3][0] + m5d[3][3] - m5d[0][0]*m5d[3][3]) - 
	      (m5d[0][3]*m5d[1][0] + m5d[1][3] - m5d[0][0]*m5d[1][3])*(m5d[0][2] + m5d[0][3]*m5d[3][2] - m5d[0][2]*m5d[3][3])) - 
	   m5d[0][3]*((m5d[0][3]*m5d[1][2] - m5d[0][2]*m5d[1][3])*(m5d[0][3]*m5d[2][0] + m5d[2][2] - m5d[0][0]*m5d[2][2]) - 
	      (m5d[0][3]*m5d[1][0] + m5d[1][3] - m5d[0][0]*m5d[1][3])*(m5d[0][3]*(-1 + m5d[2][2]) - m5d[0][2]*m5d[2][2]))*
	    (m5d[0][3]*(m5d[1][2]*m5d[3][1] + m5d[3][2] - m5d[1][1]*m5d[3][2]) - m5d[0][2]*(-1 + m5d[1][1] + m5d[1][3]*m5d[3][1] + m5d[3][3] - m5d[1][1]*m5d[3][3]) + 
	      m5d[0][1]*(m5d[1][2] + m5d[1][3]*m5d[3][2] - m5d[1][2]*m5d[3][3]))));
	
	*dispprimex0 = -((m5d[2][4]*m5d[0][2]*m5d[1][0] + m5d[3][4]*m5d[0][3]*m5d[1][0] + m5d[2][4]*m5d[1][2] - m5d[2][4]*m5d[0][0]*m5d[1][2] + m5d[3][4]*m5d[1][3] - m5d[3][4]*m5d[0][0]*m5d[1][3] + 
	   m5d[3][4]*m5d[0][3]*m5d[1][2]*m5d[2][0] - m5d[3][4]*m5d[0][2]*m5d[1][3]*m5d[2][0] - m5d[3][4]*m5d[0][3]*m5d[1][0]*m5d[2][2] - m5d[3][4]*m5d[1][3]*m5d[2][2] + 
	   m5d[3][4]*m5d[0][0]*m5d[1][3]*m5d[2][2] + m5d[3][4]*m5d[0][2]*m5d[1][0]*m5d[2][2] + m5d[3][4]*m5d[1][2]*m5d[2][2] - m5d[3][4]*m5d[0][0]*m5d[1][2]*m5d[2][2] - 
	   m5d[2][4]*m5d[0][3]*m5d[1][2]*m5d[3][0] + m5d[2][4]*m5d[0][2]*m5d[1][3]*m5d[3][0] + m5d[2][4]*m5d[0][3]*m5d[1][0]*m5d[3][2] + m5d[2][4]*m5d[1][3]*m5d[3][2] - m5d[2][4]*m5d[0][0]*m5d[1][3]*m5d[3][2] - 
	   m5d[2][4]*m5d[0][2]*m5d[1][0]*m5d[3][3] - m5d[2][4]*m5d[1][2]*m5d[3][3] + m5d[2][4]*m5d[0][0]*m5d[1][2]*m5d[3][3] - 
	   m5d[1][4]*(-1 + m5d[2][2] + m5d[0][3]*m5d[3][0] - m5d[0][3]*m5d[2][2]*m5d[3][0] + m5d[0][3]*m5d[2][0]*m5d[3][2] + m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[2][2]*m5d[3][3] + 
	      m5d[0][2]*(m5d[2][0] + m5d[2][2]*m5d[3][0] - m5d[2][0]*m5d[3][3]) - m5d[0][0]*(-1 + m5d[2][2] + m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[2][2]*m5d[3][3])) + 
	   m5d[0][4]*(m5d[1][3]*(m5d[3][0] - m5d[2][2]*m5d[3][0] + m5d[2][0]*m5d[3][2]) + m5d[1][2]*(m5d[2][0] + m5d[2][2]*m5d[3][0] - m5d[2][0]*m5d[3][3]) - 
	      m5d[1][0]*(-1 + m5d[2][2] + m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[2][2]*m5d[3][3])))/
	 (-1 + m5d[1][1] + m5d[0][2]*m5d[2][0] - m5d[0][2]*m5d[1][1]*m5d[2][0] + m5d[0][2]*m5d[1][0]*m5d[2][1] + m5d[1][2]*m5d[2][1] + m5d[2][2] - m5d[1][1]*m5d[2][2] + 
	   m5d[0][3]*m5d[3][0] - m5d[0][3]*m5d[1][1]*m5d[3][0] - m5d[0][3]*m5d[1][2]*m5d[2][1]*m5d[3][0] + m5d[0][2]*m5d[1][3]*m5d[2][1]*m5d[3][0] - m5d[0][3]*m5d[2][2]*m5d[3][0] + 
	   m5d[0][3]*m5d[1][1]*m5d[2][2]*m5d[3][0] + m5d[0][2]*m5d[2][2]*m5d[3][0] - m5d[0][2]*m5d[1][1]*m5d[2][2]*m5d[3][0] + m5d[0][3]*m5d[1][0]*m5d[3][1] + m5d[1][3]*m5d[3][1] + 
	   m5d[0][3]*m5d[1][2]*m5d[2][0]*m5d[3][1] - m5d[0][2]*m5d[1][3]*m5d[2][0]*m5d[3][1] - m5d[0][3]*m5d[1][0]*m5d[2][2]*m5d[3][1] - m5d[1][3]*m5d[2][2]*m5d[3][1] + 
	   m5d[0][2]*m5d[1][0]*m5d[2][2]*m5d[3][1] + m5d[1][2]*m5d[2][2]*m5d[3][1] + m5d[0][3]*m5d[2][0]*m5d[3][2] - m5d[0][3]*m5d[1][1]*m5d[2][0]*m5d[3][2] + m5d[0][3]*m5d[1][0]*m5d[2][1]*m5d[3][2] + 
	   m5d[1][3]*m5d[2][1]*m5d[3][2] + m5d[2][2]*m5d[3][2] - m5d[1][1]*m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[1][1]*m5d[3][3] - m5d[0][2]*m5d[2][0]*m5d[3][3] + 
	   m5d[0][2]*m5d[1][1]*m5d[2][0]*m5d[3][3] - m5d[0][2]*m5d[1][0]*m5d[2][1]*m5d[3][3] - m5d[1][2]*m5d[2][1]*m5d[3][3] - m5d[2][2]*m5d[3][3] + m5d[1][1]*m5d[2][2]*m5d[3][3] + 
	   m5d[0][1]*(m5d[1][3]*(m5d[3][0] - m5d[2][2]*m5d[3][0] + m5d[2][0]*m5d[3][2]) + m5d[1][2]*(m5d[2][0] + m5d[2][2]*m5d[3][0] - m5d[2][0]*m5d[3][3]) - 
	      m5d[1][0]*(-1 + m5d[2][2] + m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[2][2]*m5d[3][3])) - 
	   m5d[0][0]*(-1 + m5d[2][2] + m5d[1][3]*m5d[3][1] - m5d[1][3]*m5d[2][2]*m5d[3][1] + m5d[1][3]*m5d[2][1]*m5d[3][2] + m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[2][2]*m5d[3][3] + 
	      m5d[1][2]*(m5d[2][1] + m5d[2][2]*m5d[3][1] - m5d[2][1]*m5d[3][3]) - m5d[1][1]*(-1 + m5d[2][2] + m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[2][2]*m5d[3][3]))));
	
	*dispz0 = -((m5d[1][4]*m5d[0][1]*m5d[2][0] + m5d[3][4]*m5d[0][3]*m5d[2][0] - m5d[3][4]*m5d[0][3]*m5d[1][1]*m5d[2][0] + m5d[3][4]*m5d[0][1]*m5d[1][3]*m5d[2][0] + m5d[1][4]*m5d[2][1] - 
	   m5d[1][4]*m5d[0][0]*m5d[2][1] + m5d[3][4]*m5d[0][3]*m5d[1][0]*m5d[2][1] + m5d[3][4]*m5d[1][3]*m5d[2][1] - m5d[3][4]*m5d[0][0]*m5d[1][3]*m5d[2][1] + m5d[3][4]*m5d[2][2] - 
	   m5d[3][4]*m5d[0][0]*m5d[2][2] - m5d[3][4]*m5d[0][1]*m5d[1][0]*m5d[2][2] - m5d[3][4]*m5d[1][1]*m5d[2][2] + m5d[3][4]*m5d[0][0]*m5d[1][1]*m5d[2][2] - m5d[1][4]*m5d[0][3]*m5d[2][1]*m5d[3][0] + 
	   m5d[1][4]*m5d[0][1]*m5d[2][2]*m5d[3][0] + m5d[1][4]*m5d[0][3]*m5d[2][0]*m5d[3][1] + m5d[1][4]*m5d[2][2]*m5d[3][1] - m5d[1][4]*m5d[0][0]*m5d[2][2]*m5d[3][1] - 
	   m5d[1][4]*m5d[0][1]*m5d[2][0]*m5d[3][3] - m5d[1][4]*m5d[2][1]*m5d[3][3] + m5d[1][4]*m5d[0][0]*m5d[2][1]*m5d[3][3] - 
	   m5d[2][4]*(-1 + m5d[1][1] + m5d[0][3]*m5d[3][0] - m5d[0][3]*m5d[1][1]*m5d[3][0] + m5d[0][3]*m5d[1][0]*m5d[3][1] + m5d[1][3]*m5d[3][1] + m5d[3][3] - m5d[1][1]*m5d[3][3] + 
	      m5d[0][1]*(m5d[1][0] + m5d[1][3]*m5d[3][0] - m5d[1][0]*m5d[3][3]) - m5d[0][0]*(-1 + m5d[1][1] + m5d[1][3]*m5d[3][1] + m5d[3][3] - m5d[1][1]*m5d[3][3])) + 
	   m5d[0][4]*((m5d[1][3]*m5d[2][1] + m5d[2][2] - m5d[1][1]*m5d[2][2])*m5d[3][0] - m5d[2][0]*(-1 + m5d[1][1] + m5d[1][3]*m5d[3][1] + m5d[3][3] - m5d[1][1]*m5d[3][3]) + 
	      m5d[1][0]*(m5d[2][1] + m5d[2][2]*m5d[3][1] - m5d[2][1]*m5d[3][3])))/
	 (-1 + m5d[1][1] + m5d[0][2]*m5d[2][0] - m5d[0][2]*m5d[1][1]*m5d[2][0] + m5d[0][2]*m5d[1][0]*m5d[2][1] + m5d[1][2]*m5d[2][1] + m5d[2][2] - m5d[1][1]*m5d[2][2] + 
	   m5d[0][3]*m5d[3][0] - m5d[0][3]*m5d[1][1]*m5d[3][0] - m5d[0][3]*m5d[1][2]*m5d[2][1]*m5d[3][0] + m5d[0][2]*m5d[1][3]*m5d[2][1]*m5d[3][0] - m5d[0][3]*m5d[2][2]*m5d[3][0] + 
	   m5d[0][3]*m5d[1][1]*m5d[2][2]*m5d[3][0] + m5d[0][2]*m5d[2][2]*m5d[3][0] - m5d[0][2]*m5d[1][1]*m5d[2][2]*m5d[3][0] + m5d[0][3]*m5d[1][0]*m5d[3][1] + m5d[1][3]*m5d[3][1] + 
	   m5d[0][3]*m5d[1][2]*m5d[2][0]*m5d[3][1] - m5d[0][2]*m5d[1][3]*m5d[2][0]*m5d[3][1] - m5d[0][3]*m5d[1][0]*m5d[2][2]*m5d[3][1] - m5d[1][3]*m5d[2][2]*m5d[3][1] + 
	   m5d[0][2]*m5d[1][0]*m5d[2][2]*m5d[3][1] + m5d[1][2]*m5d[2][2]*m5d[3][1] + m5d[0][3]*m5d[2][0]*m5d[3][2] - m5d[0][3]*m5d[1][1]*m5d[2][0]*m5d[3][2] + m5d[0][3]*m5d[1][0]*m5d[2][1]*m5d[3][2] + 
	   m5d[1][3]*m5d[2][1]*m5d[3][2] + m5d[2][2]*m5d[3][2] - m5d[1][1]*m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[1][1]*m5d[3][3] - m5d[0][2]*m5d[2][0]*m5d[3][3] + 
	   m5d[0][2]*m5d[1][1]*m5d[2][0]*m5d[3][3] - m5d[0][2]*m5d[1][0]*m5d[2][1]*m5d[3][3] - m5d[1][2]*m5d[2][1]*m5d[3][3] - m5d[2][2]*m5d[3][3] + m5d[1][1]*m5d[2][2]*m5d[3][3] + 
	   m5d[0][1]*(m5d[1][3]*(m5d[3][0] - m5d[2][2]*m5d[3][0] + m5d[2][0]*m5d[3][2]) + m5d[1][2]*(m5d[2][0] + m5d[2][2]*m5d[3][0] - m5d[2][0]*m5d[3][3]) - 
	      m5d[1][0]*(-1 + m5d[2][2] + m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[2][2]*m5d[3][3])) - 
	   m5d[0][0]*(-1 + m5d[2][2] + m5d[1][3]*m5d[3][1] - m5d[1][3]*m5d[2][2]*m5d[3][1] + m5d[1][3]*m5d[2][1]*m5d[3][2] + m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[2][2]*m5d[3][3] + 
	      m5d[1][2]*(m5d[2][1] + m5d[2][2]*m5d[3][1] - m5d[2][1]*m5d[3][3]) - m5d[1][1]*(-1 + m5d[2][2] + m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[2][2]*m5d[3][3]))));
	
	*dispprimez0 = -((-(m5d[3][4]*(-1 + m5d[1][1] + m5d[0][2]*m5d[2][0] - m5d[0][2]*m5d[1][1]*m5d[2][0] + m5d[0][2]*m5d[1][0]*m5d[2][1] + m5d[1][2]*m5d[2][1] + m5d[2][2] - m5d[1][1]*m5d[2][2] + 
	        m5d[0][1]*(m5d[1][0] + m5d[1][2]*m5d[2][0] - m5d[1][0]*m5d[2][2]) - m5d[0][0]*(-1 + m5d[1][1] + m5d[1][2]*m5d[2][1] + m5d[2][2] - m5d[1][1]*m5d[2][2]))) + 
	   m5d[1][4]*m5d[0][1]*m5d[3][0] + m5d[2][4]*m5d[0][2]*m5d[3][0] - m5d[2][4]*m5d[0][2]*m5d[1][1]*m5d[3][0] + m5d[2][4]*m5d[0][1]*m5d[1][2]*m5d[3][0] + m5d[1][4]*m5d[0][2]*m5d[2][1]*m5d[3][0] - 
	   m5d[1][4]*m5d[0][1]*m5d[2][2]*m5d[3][0] + m5d[1][4]*m5d[3][1] - m5d[1][4]*m5d[0][0]*m5d[3][1] + m5d[2][4]*m5d[0][2]*m5d[1][0]*m5d[3][1] + m5d[2][4]*m5d[1][2]*m5d[3][1] - 
	   m5d[2][4]*m5d[0][0]*m5d[1][2]*m5d[3][1] - m5d[1][4]*m5d[0][2]*m5d[2][0]*m5d[3][1] - m5d[1][4]*m5d[2][2]*m5d[3][1] + m5d[1][4]*m5d[0][0]*m5d[2][2]*m5d[3][1] + m5d[2][4]*m5d[3][2] - 
	   m5d[2][4]*m5d[0][0]*m5d[3][2] - m5d[2][4]*m5d[0][1]*m5d[1][0]*m5d[3][2] - m5d[2][4]*m5d[1][1]*m5d[3][2] + m5d[2][4]*m5d[0][0]*m5d[1][1]*m5d[3][2] + m5d[1][4]*m5d[0][1]*m5d[2][0]*m5d[3][2] + 
	   m5d[1][4]*m5d[2][1]*m5d[3][2] - m5d[1][4]*m5d[0][0]*m5d[2][1]*m5d[3][2] + 
	   m5d[0][4]*((1 - m5d[1][1] - m5d[1][2]*m5d[2][1] - m5d[2][2] + m5d[1][1]*m5d[2][2])*m5d[3][0] + m5d[2][0]*(m5d[1][2]*m5d[3][1] + m5d[3][2] - m5d[1][1]*m5d[3][2]) + 
	      m5d[1][0]*(m5d[3][1] - m5d[2][2]*m5d[3][1] + m5d[2][1]*m5d[3][2])))/
	 (-1 + m5d[1][1] + m5d[0][2]*m5d[2][0] - m5d[0][2]*m5d[1][1]*m5d[2][0] + m5d[0][2]*m5d[1][0]*m5d[2][1] + m5d[1][2]*m5d[2][1] + m5d[2][2] - m5d[1][1]*m5d[2][2] + 
	   m5d[0][3]*m5d[3][0] - m5d[0][3]*m5d[1][1]*m5d[3][0] - m5d[0][3]*m5d[1][2]*m5d[2][1]*m5d[3][0] + m5d[0][2]*m5d[1][3]*m5d[2][1]*m5d[3][0] - m5d[0][3]*m5d[2][2]*m5d[3][0] + 
	   m5d[0][3]*m5d[1][1]*m5d[2][2]*m5d[3][0] + m5d[0][2]*m5d[2][2]*m5d[3][0] - m5d[0][2]*m5d[1][1]*m5d[2][2]*m5d[3][0] + m5d[0][3]*m5d[1][0]*m5d[3][1] + m5d[1][3]*m5d[3][1] + 
	   m5d[0][3]*m5d[1][2]*m5d[2][0]*m5d[3][1] - m5d[0][2]*m5d[1][3]*m5d[2][0]*m5d[3][1] - m5d[0][3]*m5d[1][0]*m5d[2][2]*m5d[3][1] - m5d[1][3]*m5d[2][2]*m5d[3][1] + 
	   m5d[0][2]*m5d[1][0]*m5d[2][2]*m5d[3][1] + m5d[1][2]*m5d[2][2]*m5d[3][1] + m5d[0][3]*m5d[2][0]*m5d[3][2] - m5d[0][3]*m5d[1][1]*m5d[2][0]*m5d[3][2] + m5d[0][3]*m5d[1][0]*m5d[2][1]*m5d[3][2] + 
	   m5d[1][3]*m5d[2][1]*m5d[3][2] + m5d[2][2]*m5d[3][2] - m5d[1][1]*m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[1][1]*m5d[3][3] - m5d[0][2]*m5d[2][0]*m5d[3][3] + 
	   m5d[0][2]*m5d[1][1]*m5d[2][0]*m5d[3][3] - m5d[0][2]*m5d[1][0]*m5d[2][1]*m5d[3][3] - m5d[1][2]*m5d[2][1]*m5d[3][3] - m5d[2][2]*m5d[3][3] + m5d[1][1]*m5d[2][2]*m5d[3][3] + 
	   m5d[0][1]*(m5d[1][3]*(m5d[3][0] - m5d[2][2]*m5d[3][0] + m5d[2][0]*m5d[3][2]) + m5d[1][2]*(m5d[2][0] + m5d[2][2]*m5d[3][0] - m5d[2][0]*m5d[3][3]) - 
	      m5d[1][0]*(-1 + m5d[2][2] + m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[2][2]*m5d[3][3])) - 
	   m5d[0][0]*(-1 + m5d[2][2] + m5d[1][3]*m5d[3][1] - m5d[1][3]*m5d[2][2]*m5d[3][1] + m5d[1][3]*m5d[2][1]*m5d[3][2] + m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[2][2]*m5d[3][3] + 
	      m5d[1][2]*(m5d[2][1] + m5d[2][2]*m5d[3][1] - m5d[2][1]*m5d[3][3]) - m5d[1][1]*(-1 + m5d[2][2] + m5d[2][2]*m5d[3][2] + m5d[3][3] - m5d[2][2]*m5d[3][3]))));

	// test
	//test_vect[0] = *dispx0;
	//test_vect[1] = *dispprimex0;
	//test_vect[2] = *dispz0;
	//test_vect[3] = *dispprimez0;
	//test_vect[4] = 1;
	//mvprod5(test_vect, m5d, test_vect);
	//printf("get_periodic_dispersion_4d: dif = %le, %le, %le, %le\n", test_vect[0]-*dispx0, test_vect[1]-*dispprimex0, test_vect[2]-*dispz0, test_vect[3]-*dispprimez0);
}
//*/

// ************************************************************************************ //
//								linear matrix functions									//
// ************************************************************************************ //

//find actual distance in horizontal from the closed orbit for acurate transfer matrix computation
extern void find_x(double *x, double *bias, struct Particle *part_ref, struct Particle *part)
{
	double xprime, xprime_ref, dr;
	
	if(part->y > 2*TINYLENGTH) printf("!Warning1 in find_x!, y=%le\n", part->y);
	if(part_ref->y > 2*TINYLENGTH) printf("!Warning2 in find_x!, yref=%le\n", part_ref->y);
	
	xprime = atan_ratio(part->ux, part->uy);
	xprime_ref = atan_ratio(part_ref->ux, part_ref->uy);
	dr = part->x - part_ref->x;
	
	*bias = dr*sin(xprime_ref);//<>=0
	*x = dr*cos(xprime_ref) + *bias*tan(xprime - xprime_ref);//<>=0
	//*x = dr;
}

//find actual distance in vertical from the closed orbit for acurate transfer matrix computation
extern void find_z(double *z, double *bias, struct Particle *part_ref, struct Particle *part)
{
	double uhp, uhp_ref, zprime, zprime_ref, dr;
	
	if(part->y > 2*TINYLENGTH) printf("!Warning1 in find_z!, y=%le\n", part->y);
	if(part_ref->y > 2*TINYLENGTH) printf("!Warning2 in find_z!, yref=%le\n", part_ref->y);
	
	uhp = sqrt(part->ux*part->ux + part->uy*part->uy);
	uhp_ref = sqrt(part_ref->ux*part_ref->ux + part_ref->uy*part_ref->uy);
	zprime = atan_ratio(part->uz, uhp);
	zprime_ref = atan_ratio(part_ref->uz, uhp_ref);
	dr = part->z - part_ref->z;
	
	*bias = dr*sin(zprime_ref);//<>=0
	*z = dr*cos(zprime_ref) + *bias*tan(zprime - zprime_ref);//<>=0
	//*z = dr;
}

//general function to output horizontal & vertical (or uncoupled planes) transfer matrix
extern int get_matrix_hor_vert(double mh[2][2], double mv[2][2], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf, int together_decoupled)
{
	int i,j,n=4;
	double m4d[4][4], decoup_m4d[4][4], output[6];
	double decoup_to_coup[4][4];
	
	if(together_decoupled == 0) {
		if(get_trmatrix_h(mh, band_bias, reference, amp_x, amp_xprime, latt, *transport_part, doyouprintf) != TRUE) return FALSE;
		if(get_trmatrix_v(mv, reference, amp_z, amp_zprime, latt, *transport_part, doyouprintf) != TRUE) return FALSE;
	}
	else {
		if(get_matrix_firstorder_4d(m4d, band_bias, reference, amp_x, amp_xprime, amp_z, amp_zprime, latt, *transport_part, doyouprintf) != TRUE) return FALSE;
		//if(get_matrix_secondorder_4d(m4d, band_bias, reference, amp_x, amp_xprime, amp_z, amp_zprime, latt, *transport_part, doyouprintf) != TRUE) return FALSE;
		//for(i=0;i<n;i++) {
		//	for(j=0;j<n;j++) printf("%le\t", m4d[i][j]);
		//	printf("\n");
		//}
		//printf("\n");
		if(together_decoupled == 1) {
			if(decouple_matrix_parzen(m4d, decoup_m4d, decoup_to_coup, output, NULL) != TRUE) return FALSE;
			for(i=0;i<n;i++) for(j=0;j<n;j++) m4d[i][j] = decoup_m4d[i][j];
		}
		for(i=0;i<2;i++) {
			for(j=0;j<2;j++) {
				mh[i][j] = m4d[i][j];
				mv[i][j] = m4d[i+2][j+2];
			}
		}
	}
	return TRUE;
}

//get equivalent horizontal transfer matrix
extern int get_trmatrix_h(double mh[2][2], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf)
//output the linear transfer matrix from the computation of the transfer matrix up to the second order.
//return FALSE if a problem had happened during the computation of matrix coefficients, else return TRUE.
//IMPORTANT!: transfer matrix calculated around the orbit of *reference (if you want *reference to be on the closed orbit, you have to make sure of that before)
//IMPORTANT2!:  if your xprime_ref ~= 0 at both entrance and exit faces of the latt, you are not concerned by this warning. 
//			    Else, notice that this procedure works even when xprime_ref !=0, BUT it is then assumed 
//				that the particle goes straight on a distance <= band_bias around the latt boundary.
//				If fields are actually non-zero between +/- band_bias around the latt boundary, the result will be biased (more or less depending on "how much" the particle does not go straight over this distance).
//				Check the value of band_bias for safety. Must be "small" if you have a "significant" field around the boundary...
//reference is assumed to stay in the horizontal mid-plane. 
{
	double xprime_ref, uhp_ref,
		x_pdx, x_mdx, x_pdxp, x_mdxp, bias,
		xprime_pdx, xprime_mdx, xprime_pdxp, xprime_mdxp;
	struct Particle part_ref, part_pdx, part_mdx, part_pdxp, part_mdxp;
//	FILE *wfile;
//	wfile = fopen("data/transfermatrix.dat", "w");
	*band_bias = 0;
	mh[0][0] = mh[0][1] = mh[1][0] = mh[1][1] = 0;
	
	if(amp_x == 0 || amp_xprime == 0) errorstop("!!!ERROR: you call the procedure ''get_trmatrix_h'' with an invalid argument: amp_x(or)xprime = 0");
	if(fabs(reference->z) > 1.e-5 || fabs(reference->uz) > 1.e-5) {
		printf("!WARNING in get_trmatrix_h: reference particle z or uz are not 0. This case is not implemented yet. HORIZONTAL transfer matrix won't be calculated\n");
		return FALSE;
	}
	
	
	//reference particle initial angle in the Cell local (entrance face) framework
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	if(cos(xprime_ref) == 0) {
		printf("!Warning in get_trmatrix_h: cos(xprime_ref) == 0, no transfer matrix will be computed\n");
		return FALSE;
	}
	
	uhp_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	
	
	//++//Get 2X2 horizontal transfer matrix
	//init
	part_ref = *reference;
	part_ref.hat = -2; //no acceleration, no output in trackout.dat
	
	part_pdx = part_ref;
	part_pdx.x += amp_x/cos(xprime_ref);
	
	part_mdx = part_ref;
	part_mdx.x -= amp_x/cos(xprime_ref);
	
	part_pdxp= part_ref;
	part_pdxp.ux = uhp_ref*sin(xprime_ref + amp_xprime);
	part_pdxp.uy = uhp_ref*cos(xprime_ref + amp_xprime);
	
	part_mdxp = part_ref;
	part_mdxp.ux = uhp_ref*sin(xprime_ref - amp_xprime);
	part_mdxp.uy = uhp_ref*cos(xprime_ref - amp_xprime);
	
	*band_bias = MAX(*band_bias, amp_x/cos(xprime_ref));
	
	//transport the 5 particles
	(*transport_part)(&part_ref, latt,NULL);
	if(part_ref.status != ALIVE) {
		printf("part_ref not alive\n");
		return FALSE;
	}
	(*transport_part)(&part_pdx, latt,NULL);
	if(part_pdx.status != ALIVE) {
		printf("part_pdx not alive\n");
		return FALSE;	
	}
	(*transport_part)(&part_mdx, latt,NULL);
	if(part_mdx.status != ALIVE) {
		printf("part_mdx not alive\n");
		return FALSE;	
	}
	(*transport_part)(&part_pdxp, latt,NULL);
	if(part_pdxp.status != ALIVE) {
		printf("part_pdxp not alive\n");
		return FALSE;	
	}
	(*transport_part)(&part_mdxp, latt,NULL);
	if(part_mdxp.status != ALIVE) {
		printf("part_mdxp not alive\n");
		return FALSE;	
	}
	
	//compute matrix coef.
	find_x(&x_pdx, &bias, &part_ref, &part_pdx);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdx, &bias, &part_ref, &part_mdx);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdxp, &bias, &part_ref, &part_pdxp);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdxp, &bias, &part_ref, &part_mdxp);
	*band_bias = MAX(*band_bias, bias);
	
	xprime_ref = atan_ratio(part_ref.ux, part_ref.uy);
	xprime_pdx = atan_ratio(part_pdx.ux, part_pdx.uy) - xprime_ref;
	xprime_mdx = atan_ratio(part_mdx.ux, part_mdx.uy) - xprime_ref;
	xprime_pdxp = atan_ratio(part_pdxp.ux, part_pdxp.uy) - xprime_ref;
	xprime_mdxp = atan_ratio(part_mdxp.ux, part_mdxp.uy) - xprime_ref;
	
//	fprintf(wfile, "%lf	%lf\n", 0., xprime_ref);
//	fprintf(wfile, "%lf	%lf\n", x_pdx, xprime_pdx);
//	fprintf(wfile, "%lf	%lf\n", x_mdx, xprime_mdx);
//	fprintf(wfile, "%lf	%lf\n", x_pdxp, xprime_pdxp);
//	fprintf(wfile, "%lf	%lf\n", x_mdxp, xprime_mdxp);
//	fclose(wfile);
	mh[0][0] = (x_pdx - x_mdx)/(2.*amp_x); //second order canceled out
	mh[0][1] = (x_pdxp - x_mdxp)/(2.*amp_xprime); //second order canceled out
	mh[1][0] = (xprime_pdx - xprime_mdx)/(2.*amp_x); //second order canceled out
	mh[1][1] = (xprime_pdxp - xprime_mdxp)/(2.*amp_xprime); //second order canceled out
	
	if(doyouprintf == YES) {
		printf("\n------------------------ Get equivalent HORIZONTAL transfer matrix ---------------------------\n");
		printf("mh[0][0] = %lf, mh[0][1] = %lf \nmh[1][0] = %lf, mh[1][1] = %lf \t", mh[0][0], mh[0][1], mh[1][0], mh[1][1]);
		printf("Determinant = %lf\n", mh[0][0]*mh[1][1] - mh[0][1]*mh[1][0]);
	}
	
	return TRUE;
}

//get equivalent vertical transfer matrix
extern int get_trmatrix_v(double mv[2][2], struct Particle *reference, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*,char*), int doyouprintf)
//output the linear transfer matrix from the computation of the transfer matrix up to the second order.
//return FALSE if a problem had happened during the computation of matrix coefficients, else return TRUE.
//reference is assumed to stay in the horizontal mid-plane. 
{
	double xprime_ref, uhp_ref,
		z_pdz, z_mdz, z_pdzp, z_mdzp,
		zprime_pdz, zprime_mdz, zprime_pdzp, zprime_mdzp;
	//double not_symplec[2][2];
	struct Particle part_ref, part_pdz, part_mdz, part_pdzp, part_mdzp;
	
	mv[0][0] = mv[0][1] = mv[1][0] = mv[1][1] = 0;
	
	if(amp_z == 0 || amp_zprime == 0) errorstop("!!!ERROR: you call the procedure ''get_trmatrix_v'' with an invalid argument: amp_z(or)zprime = 0");
	if(fabs(reference->z) > 1.e-5 || fabs(reference->uz) > 1.e-5) {
		printf("!WARNING in get_trmatrix_v: reference particle z or uz are not 0. This case is not implemented yet. VERTICAL transfer matrix won't be calculated\n");
		return FALSE;
	}
	
	//reference particle initial angle in the Cell local (entrance face) framework
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	if(cos(xprime_ref) == 0) {
		printf("!Warning in get_trmatrix_v: cos(xprime_ref) == 0, no transfer matrix will be computed\n");
		return FALSE;
	}
	//reference particle total momentum
	uhp_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	
	//++//Get 2X2 vertical transfer matrix
	//init
	part_ref = *reference;
	part_ref.hat = -2; //no acceleration, no output in trackout.dat
	
	part_pdz = part_ref;
	part_pdz.z += amp_z;
	
	part_mdz = part_ref;
	part_mdz.z -= amp_z;
	
	part_pdzp= part_ref;
	part_pdzp.ux = sin(xprime_ref)*cos(amp_zprime);
	part_pdzp.uy = cos(xprime_ref)*cos(amp_zprime);
	part_pdzp.uz = sin(amp_zprime);
	
	part_mdzp = part_ref;
	part_mdzp.ux = sin(xprime_ref)*cos(-amp_zprime);
	part_mdzp.uy = cos(xprime_ref)*cos(-amp_zprime);
	part_mdzp.uz = sin(-amp_zprime);
	
	//transport the 5 particles
	(*transport_part)(&part_ref, latt,NULL);
	if(part_ref.status != ALIVE) {
		printf("part_ref not alive\n");
		return FALSE;	
	}
	(*transport_part)(&part_pdz, latt,NULL);
	if(part_pdz.status != ALIVE) {
		printf("part_pdz not alive\n");
		return FALSE;	
	}	
	(*transport_part)(&part_mdz, latt,NULL);
	if(part_mdz.status != ALIVE) {
		printf("part_mdz not alive\n");
		return FALSE;	
	}
	(*transport_part)(&part_pdzp, latt,NULL);
	if(part_pdzp.status != ALIVE) {
		printf("part_pdzp not alive\n");
		return FALSE;	
	}
	(*transport_part)(&part_mdzp, latt,NULL);
	if(part_mdzp.status != ALIVE) {
		printf("part_mdzp not alive\n");
		return FALSE;	
	}
	
	
	//compute matrix coef.
	if(fabs(reference->z) > 1.e-3 || fabs(reference->uz) > 1.e-3) {
		printf("!WARNING in get_trmatrix_v: after transport, reference particle z or uz are not 0. This case is not implemented yet. VERTICAL transfer matrix won't be calculated\n");
		return FALSE;
	}
	
	z_pdz = part_pdz.z;
	z_mdz = part_mdz.z;
	z_pdzp = part_pdzp.z;
	z_mdzp = part_mdzp.z;
	
	zprime_pdz = atan_ratio(part_pdz.uz, uhp_ref);
	zprime_mdz = atan_ratio(part_mdz.uz, uhp_ref);
	zprime_pdzp = atan_ratio(part_pdzp.uz, uhp_ref);
	zprime_mdzp = atan_ratio(part_mdzp.uz, uhp_ref);
	
	mv[0][0] = (z_pdz - z_mdz)/(2.*amp_z); //second order canceled out
	mv[0][1] = (z_pdzp - z_mdzp)/(2.*amp_zprime); //second order canceled out
	mv[1][0] = (zprime_pdz - zprime_mdz)/(2.*amp_z); //second order canceled out
	mv[1][1] = (zprime_pdzp - zprime_mdzp)/(2.*amp_zprime); //second order canceled out
	
	if(doyouprintf == YES) {
		printf("\n------------------------ Get equivalent VERTICAL transfer matrix ---------------------------\n");
		printf("mv[0][0] = %lf, mv[0][1] = %lf \nmv[1][0] = %lf, mv[1][1] = %lf \t", mv[0][0], mv[0][1], mv[1][0], mv[1][1]);
		printf("Determinant = %lf\n", mv[0][0]*mv[1][1] - mv[0][1]*mv[1][0]);
	}
	
	return TRUE;
}

extern int get_trmatrix_h3x3(double mh[3][3], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, double amp_dpovp, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*,char*), int doyouprintf)
//output the linear (first order) horizontal transfer matrix from the computation of the transfer matrix up to the second order.
//return FALSE if a problem had happened during the computation of matrix coefficients, else return TRUE.
//IMPORTANT!: transfer matrix calculated around the orbit of *reference (if you want *reference to be on the closed orbit, you have to make sure of that before)
//IMPORTANT2!: if your xprime_ref ~= 0 at both entrance and exit faces of the latt, you are not concerned by this warning. 
//			   Else, notice that this procedure works even when xprime_ref !=0, BUT it is then assumed that the particle goes straight on a distance <= band_bias around the latt boundary.
//			   If fields are actually non-zero between +/- band_bias around the latt boundary, the result will be biased (more or less depending on "how much" the particle does not go straight over this distance).
//			   Check the value of band_bias for safety. Must be "small" if you have a "significant" field around the boundary...
//reference is assumed to stay in the horizontal mid-plane. 
{
	int i,j,n=3;
	double xprime_ref, uhp_ref,
	x_pdx, x_mdx, x_pdxp, x_mdxp, bias,
	xprime_pdx, xprime_mdx, xprime_pdxp, xprime_mdxp,
	x_pdpovp, x_mdpovp, xprime_pdpovp, xprime_mdpovp;
	struct Particle part_ref, part_pdx, part_mdx, part_pdxp, part_mdxp, part_pdpovp, part_mdpovp;
	
	*band_bias = 0;
	for(i=0;i<n;i++) for(j=0;j<n;j++) mh[i][j]=0;
	
	if(amp_x == 0 || amp_xprime == 0 || amp_dpovp == 0) errorstop("!!!ERROR: you call the procedure ''get_trmatrix_h3x3'' with an invalid argument: amp_x(or)xprime(or)dpovp = 0");
	if(fabs(reference->z) > 0 || fabs(reference->uz) > 0) {
		printf("!WARNING in get_trmatrix_h: reference particle z or uz are not 0. This case is not implemented yet. HORIZONTAL transfer matrix won't be calculated\n");
		return FALSE;
	}
	
	
	//reference particle initial angle in the Cell local (entrance face) framework
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	if(cos(xprime_ref) == 0) {
		printf("!Warning in get_trmatrix_h: cos(xprime_ref) == 0, no transfer matrix will be computed\n");
		return FALSE;
	}
	
	uhp_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	
	
	//++//Get 3X3 horizontal transfer matrix
	//init
	part_ref = *reference;
	part_ref.hat = -2; //no acceleration, no output in trackout.dat
	
	part_pdx = part_ref;
	part_pdx.x += amp_x/cos(xprime_ref);
	
	part_mdx = part_ref;
	part_mdx.x -= amp_x/cos(xprime_ref);
	
	part_pdxp= part_ref;
	part_pdxp.ux = uhp_ref*sin(xprime_ref + amp_xprime);
	part_pdxp.uy = uhp_ref*cos(xprime_ref + amp_xprime);
	
	part_mdxp = part_ref;
	part_mdxp.ux = uhp_ref*sin(xprime_ref - amp_xprime);
	part_mdxp.uy = uhp_ref*cos(xprime_ref - amp_xprime);
	
	part_pdpovp = part_ref;
	part_pdpovp.brho *= (1 + amp_dpovp);
	
	part_mdpovp = part_ref;
	part_mdpovp.brho *= (1 - amp_dpovp);
	
	*band_bias = MAX(*band_bias, amp_x/cos(xprime_ref));
	
	//transport the 7 particles
	(*transport_part)(&part_ref, latt,NULL);
	if(part_ref.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdx, latt,NULL);
	if(part_pdx.status != ALIVE) return FALSE;	
	(*transport_part)(&part_mdx, latt,NULL);
	if(part_mdx.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdxp, latt,NULL);
	if(part_pdxp.status != ALIVE) return FALSE;
	(*transport_part)(&part_mdxp, latt,NULL);
	if(part_mdxp.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdpovp, latt,NULL);
	if(part_pdpovp.status != ALIVE) return FALSE;
	(*transport_part)(&part_mdpovp, latt,NULL);
	if(part_mdpovp.status != ALIVE) return FALSE;
	
	//compute matrix coef.
	find_x(&x_pdx, &bias, &part_ref, &part_pdx);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdx, &bias, &part_ref, &part_mdx);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdxp, &bias, &part_ref, &part_pdxp);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdxp, &bias, &part_ref, &part_mdxp);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdpovp, &bias, &part_ref, &part_pdpovp);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdpovp, &bias, &part_ref, &part_mdpovp);
	*band_bias = MAX(*band_bias, bias);
	
	xprime_ref = atan_ratio(part_ref.ux, part_ref.uy);
	xprime_pdx = atan_ratio(part_pdx.ux, part_pdx.uy) - xprime_ref;
	xprime_mdx = atan_ratio(part_mdx.ux, part_mdx.uy) - xprime_ref;
	xprime_pdxp = atan_ratio(part_pdxp.ux, part_pdxp.uy) - xprime_ref;
	xprime_mdxp = atan_ratio(part_mdxp.ux, part_mdxp.uy) - xprime_ref;
	xprime_pdpovp = atan_ratio(part_pdpovp.ux, part_pdpovp.uy) - xprime_ref;
	xprime_mdpovp = atan_ratio(part_mdpovp.ux, part_mdpovp.uy) - xprime_ref;
	
	mh[0][0] = (x_pdx - x_mdx)/(2.*amp_x); //second order cancelled out
	mh[0][1] = (x_pdxp - x_mdxp)/(2.*amp_xprime); //second order cancelled out
	mh[0][2] = (x_pdpovp - x_mdpovp)/(2.*amp_dpovp); //second order cancelled out
	mh[1][0] = (xprime_pdx - xprime_mdx)/(2.*amp_x); //second order cancelled out
	mh[1][1] = (xprime_pdxp - xprime_mdxp)/(2.*amp_xprime); //second order cancelled out
	mh[1][2] = (xprime_pdpovp - xprime_mdpovp)/(2.*amp_dpovp); //second order cancelled out
	mh[2][0] = 0;
	mh[2][1] = 0;
	mh[2][2] = 1.;
	
	if(doyouprintf == YES) {
		double m_det[4][4], determinant=0.0;
		for(i=0;i<n;i++) for (j=0;j<n;j++) m_det[i][j] = mh[i][j];
		determinant = matrix_det_4d(m_det, n);
		printf("\n------------------------ Get equivalent HORIZONTAL transfer matrix ---------------------------\n");
		printf("mh[0][0] = %lf, mh[0][1] = %lf, mh[0][2] = %lf\nmh[1][0] = %lf, mh[1][1] = %lf, mh[1][2] = %lf\nmh[2][0] = %lf, mh[2][1] = %lf, mh[2][2] = %lf\t", mh[0][0], mh[0][1], mh[0][2], mh[1][0], mh[1][1], mh[1][2], mh[2][0], mh[2][1], mh[2][2]);
		printf("Determinant = %lf\n", determinant);
	}
	
	return TRUE;
}

extern int get_trmatrix_v3x3(double mv[3][3], double *band_bias, struct Particle *reference, double amp_z, double amp_zprime, double amp_dpovp, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*,char*), int doyouprintf)
//output the linear (first order) vertical transfer matrix from the computation of the transfer matrix up to the second order.
//return FALSE if a problem had happened during the computation of matrix coefficients, else return TRUE.
//IMPORTANT!: transfer matrix calculated around the orbit of *reference (if you want *reference to be on the closed orbit, you have to make sure of that before)
//IMPORTANT2!: if your xprime_ref ~= 0 at both entrance and exit faces of the latt, you are not concerned by this warning. 
//			   Else, notice that this procedure works even when xprime_ref !=0, BUT it is then assumed that the particle goes straight on a distance <= band_bias around the latt boundary.
//			   If fields are actually non-zero between +/- band_bias around the latt boundary, the result will be biased (more or less depending on "how much" the particle does not go straight over this distance).
//			   Check the value of band_bias for safety. Must be "small" if you have a "significant" field around the boundary...
{
	int i,j,n=3;
	double xprime_ref, zprime_ref, uhp_ref, z_pdz, z_mdz, z_pdzp, z_mdzp, bias, uhp, zprime_pdz, zprime_mdz, zprime_pdzp, zprime_mdzp,
	z_pdpovp, z_mdpovp, zprime_pdpovp, zprime_mdpovp;
	struct Particle part_ref, part_pdz, part_mdz, part_pdzp, part_mdzp, part_pdpovp, part_mdpovp;
	
	*band_bias = 0;
	for(i=0;i<n;i++) for(j=0;j<n;j++) mv[i][j]=0;
	
	if(amp_z == 0 || amp_zprime == 0 || amp_dpovp == 0) errorstop("!!!ERROR: you call the procedure ''get_trmatrix_v3x3'' with an invalid argument: amp_z(or)zprime(or)dpovp = 0");
	
	
	//reference particle initial angle in the Cell local (entrance face) framework
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	uhp_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	zprime_ref = atan_ratio(reference->uz, uhp_ref);
	if(cos(xprime_ref) == 0) {
		printf("!Warning in get_trmatrix_v3x3: cos(xprime_ref) == 0, no transfer matrix will be computed\n");
		return FALSE;
	}
	
	
	//++//Get 3X3 horizontal transfer matrix
	//init
	part_ref = *reference;
	part_ref.hat = -2; //no acceleration, no output in trackout.dat
	
	part_pdz = part_ref;
	part_pdz.z += amp_z/cos(zprime_ref);
	
	part_mdz = part_ref;
	part_mdz.z -= amp_z/cos(zprime_ref);
	
	part_pdzp= part_ref;
	part_pdzp.ux = sin(xprime_ref)*cos(zprime_ref + amp_zprime);
	part_pdzp.uy = cos(xprime_ref)*cos(zprime_ref + amp_zprime);
	part_pdzp.uz = sin(zprime_ref + amp_zprime);
	
	part_mdzp = part_ref;
	part_mdzp.ux = sin(xprime_ref)*cos(zprime_ref - amp_zprime);
	part_mdzp.uy = cos(xprime_ref)*cos(zprime_ref - amp_zprime);
	part_mdzp.uz = sin(zprime_ref - amp_zprime);
	
	part_pdpovp = part_ref;
	part_pdpovp.brho *= (1 + amp_dpovp);
	
	part_mdpovp = part_ref;
	part_mdpovp.brho *= (1 - amp_dpovp);
	
	*band_bias = MAX(*band_bias, amp_z/cos(zprime_ref));
	
	//transport the 7 particles
	(*transport_part)(&part_ref, latt,NULL);
	if(part_ref.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdz, latt,NULL);
	if(part_pdz.status != ALIVE) return FALSE;	
	(*transport_part)(&part_mdz, latt,NULL);
	if(part_mdz.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdzp, latt,NULL);
	if(part_pdzp.status != ALIVE) return FALSE;
	(*transport_part)(&part_mdzp, latt,NULL);
	if(part_mdzp.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdpovp, latt,NULL);
	if(part_pdpovp.status != ALIVE) return FALSE;
	(*transport_part)(&part_mdpovp, latt,NULL);
	if(part_mdpovp.status != ALIVE) return FALSE;
	
	//compute matrix coef.
	find_z(&z_pdz, &bias, &part_ref, &part_pdz);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_mdz, &bias, &part_ref, &part_mdz);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdzp, &bias, &part_ref, &part_pdzp);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_mdzp, &bias, &part_ref, &part_mdzp);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdpovp, &bias, &part_ref, &part_pdpovp);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_mdpovp, &bias, &part_ref, &part_mdpovp);
	*band_bias = MAX(*band_bias, bias);
	
	uhp_ref = sqrt(part_ref.ux*part_ref.ux + part_ref.uy*part_ref.uy);
	zprime_ref  = atan_ratio(part_ref.uz, uhp_ref);
	uhp = sqrt(part_pdz.ux*part_pdz.ux + part_pdz.uy*part_pdz.uy);
	zprime_pdz = atan_ratio(part_pdz.uz, uhp) - zprime_ref;
	uhp = sqrt(part_mdz.ux*part_mdz.ux + part_mdz.uy*part_mdz.uy);
	zprime_mdz = atan_ratio(part_mdz.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdzp.ux*part_pdzp.ux + part_pdzp.uy*part_pdzp.uy);
	zprime_pdzp = atan_ratio(part_pdzp.uz, uhp) - zprime_ref;
	uhp = sqrt(part_mdzp.ux*part_mdzp.ux + part_mdzp.uy*part_mdzp.uy);
	zprime_mdzp = atan_ratio(part_mdzp.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdpovp.ux*part_pdpovp.ux + part_pdpovp.uy*part_pdpovp.uy);
	zprime_pdpovp = atan_ratio(part_pdpovp.uz, uhp) - zprime_ref; 
	uhp = sqrt(part_mdpovp.ux*part_mdpovp.ux + part_mdpovp.uy*part_mdpovp.uy);
	zprime_mdpovp = atan_ratio(part_mdpovp.uz, uhp) - zprime_ref; 
	
	mv[0][0] = (z_pdz - z_mdz)/(2.*amp_z); //second order cancelled out
	mv[0][1] = (z_pdzp - z_mdzp)/(2.*amp_zprime); //second order cancelled out
	mv[0][2] = (z_pdpovp - z_mdpovp)/(2.*amp_dpovp); //second order cancelled out
	mv[1][0] = (zprime_pdz - zprime_mdz)/(2.*amp_z); //second order cancelled out
	mv[1][1] = (zprime_pdzp - zprime_mdzp)/(2.*amp_zprime); //second order cancelled out
	mv[1][2] = (zprime_pdpovp - zprime_mdpovp)/(2.*amp_dpovp); //second order cancelled out
	mv[2][0] = 0;
	mv[2][1] = 0;
	mv[2][2] = 1.;
	
	if(doyouprintf == YES) {
		double m_det[4][4], determinant=0.;
		for(i=0;i<n;i++) for (j=0;j<n;j++) m_det[i][j] = mv[i][j];
		determinant = matrix_det_4d(m_det, n);
		printf("\n------------------------ Get equivalent VERTICAL transfer matrix ---------------------------\n");
		printf("mv[0][0] = %lf, mv[0][1] = %lf, mv[0][2] = %lf\nmv[1][0] = %lf, mv[1][1] = %lf, mv[1][2] = %lf\nmv[2][0] = %lf, mv[2][1] = %lf, mv[2][2] = %lf\t", mv[0][0], mv[0][1], mv[0][2], mv[1][0], mv[1][1], mv[1][2], mv[2][0], mv[2][1], mv[2][2]);
		printf("Determinant = %lf\n", determinant);
	}
	return TRUE;
}

extern int get_matrix_firstorder_4d(double m[4][4], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf)
//output the linear 4D transfer matrix from the computation of the transfer matrix to the first order.
//return FALSE if a problem had happened during the computation of matrix coefficients, else return TRUE.
//IMPORTANT!: transfer matrix calculated around the orbit of *reference (if you want *reference to be on the closed orbit, you have to make sure of that before)
//IMPORTANT2!:  if your xprime_ref ~= 0 at both entrance and exit faces of the latt, you are not concerned by this warning. 
//			    Else, notice that this procedure works even when xprime_ref !=0, BUT it is then assumed 
//				that the particle goes straight on a distance <= band_bias around the latt boundary.
//				If fields are actually non-zero between +/- band_bias around the latt boundary, the result will be biased (more or less depending on "how much" the particle does not go straight over this distance).
//				Check the value of band_bias for safety. Must be "small" if you have a "significant" field around the boundary...
{
	int i,j,n=4;
	double xprime_ref, zprime_ref, uhp, uhp_ref, bias;
	double x_pdx, x_pdux, x_pdz, x_pduz, xprime_pdx, xprime_pdux, xprime_pdz, xprime_pduz;
	double z_pdx, z_pdux, z_pdz, z_pduz, zprime_pdx, zprime_pdux, zprime_pdz, zprime_pduz;
	//double not_symplec[4][4];
	struct Particle part_ref, part_pdx, part_pdux, part_pdz, part_pduz;
	
	for(i=0;i<n;i++) for(j=0;j<n;j++) m[i][j] = 0.;
	if(amp_x == 0 || amp_xprime == 0 || amp_z == 0 || amp_zprime == 0) errorstop("!!!ERROR: you call the procedure ''get_trmatrix_firstorder'' with an invalid argument: amp_x/z(or)ux/uz = 0");
	
	*band_bias = 0;
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	if(cos(xprime_ref) == 0) {
		printf("!Warning in get_matrix_firstorder_4d: cos(xprime_ref) == 0, no transfer matrix will be computed\n");
		return FALSE;
	}
	uhp_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	zprime_ref = atan_ratio(reference->uz, uhp_ref);
	if(cos(zprime_ref) == 0) {
		printf("!Warning in get_matrix_firstorder_4d: cos(zprime_ref) == 0, no transfer matrix will be computed\n");
		return FALSE;
	}
	//init
	part_ref = *reference;
	part_ref.hat = -2; //no acceleration, no output in trackout.dat
	
	part_pdx = part_ref;
	part_pdx.x += amp_x/cos(xprime_ref);
	
	part_pdux = part_ref;
	part_pdux.ux = sin(xprime_ref + amp_xprime)*cos(zprime_ref); 
	part_pdux.uy = cos(xprime_ref + amp_xprime)*cos(zprime_ref); 
	part_pdux.uz = sin(zprime_ref);
	
	part_pdz = part_ref;
	part_pdz.z += amp_z/cos(zprime_ref);
	
	part_pduz = part_ref;
	part_pduz.ux = sin(xprime_ref)*cos(zprime_ref+amp_zprime);
	part_pduz.uy = cos(xprime_ref)*cos(zprime_ref+amp_zprime);
	part_pduz.uz = sin(zprime_ref+amp_zprime);
	
	*band_bias = MAX(amp_z/cos(zprime_ref), amp_x/cos(xprime_ref));
	
	//transport the 5 particles
	//printf("part_ref\n");
	(*transport_part)(&part_ref, latt,NULL);
	if(part_ref.status != ALIVE) return FALSE;
	//printf("part_pdx\n");
	(*transport_part)(&part_pdx, latt,NULL);
	if(part_pdx.status != ALIVE) return FALSE;
	//printf("part_pdux\n");
	(*transport_part)(&part_pdux, latt,NULL);
	if(part_pdux.status != ALIVE) return FALSE;
	//printf("part_pdz\n");
	(*transport_part)(&part_pdz, latt,NULL);
	if(part_pdz.status != ALIVE) return FALSE;	
	//printf("part_pduz\n");
	(*transport_part)(&part_pduz, latt,NULL);
	if(part_pduz.status != ALIVE) return FALSE;
	
	find_x(&x_pdx, &bias, &part_ref, &part_pdx);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdux, &bias, &part_ref, &part_pdux);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdz, &bias, &part_ref, &part_pdz);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pduz, &bias, &part_ref, &part_pduz);
	*band_bias = MAX(*band_bias, bias);
	
	xprime_ref  = atan_ratio(part_ref.ux, part_ref.uy);
	xprime_pdx  = atan_ratio(part_pdx.ux, part_pdx.uy) - xprime_ref;
	xprime_pdux = atan_ratio(part_pdux.ux, part_pdux.uy) - xprime_ref;
	xprime_pdz  = atan_ratio(part_pdz.ux, part_pdz.uy) - xprime_ref;
	xprime_pduz = atan_ratio(part_pduz.ux, part_pduz.uy) - xprime_ref;
	
	find_z(&z_pdx, &bias, &part_ref, &part_pdx);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdux, &bias, &part_ref, &part_pdux);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdz, &bias, &part_ref, &part_pdz);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pduz, &bias, &part_ref, &part_pduz);
	*band_bias = MAX(*band_bias, bias);
	
	uhp_ref = sqrt(part_ref.ux*part_ref.ux + part_ref.uy*part_ref.uy);
	zprime_ref  = atan_ratio(part_ref.uz, uhp_ref);
	uhp = sqrt(part_pdx.ux*part_pdx.ux + part_pdx.uy*part_pdx.uy);
	zprime_pdx  = atan_ratio(part_pdx.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdux.ux*part_pdux.ux + part_pdux.uy*part_pdux.uy);
	zprime_pdux = atan_ratio(part_pdux.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdz.ux*part_pdz.ux + part_pdz.uy*part_pdz.uy);
	zprime_pdz  = atan_ratio(part_pdz.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pduz.ux*part_pduz.ux + part_pduz.uy*part_pduz.uy);
	zprime_pduz = atan_ratio(part_pduz.uz, uhp) - zprime_ref;
	
	
	m[0][0] = x_pdx/amp_x;     //(part_pdx.x  - part_ref.x)/(amp_x);
	m[1][0] = xprime_pdx/amp_x;//(part_pdx.ux - part_ref.ux)/(amp_x);
	m[2][0] = z_pdx/amp_x;     //(part_pdx.z  - part_ref.z)/(amp_x);
	m[3][0] = zprime_pdx/amp_x;//(part_pdx.uz - part_ref.uz)/(amp_x);
	
	m[0][1] = x_pdux/amp_xprime;     //(part_pdux.x  - part_ref.x)/(amp_xprime);
	m[1][1] = xprime_pdux/amp_xprime;//(part_pdux.ux - part_ref.ux)/(amp_xprime);
	m[2][1] = z_pdux/amp_xprime;     //(part_pdux.z  - part_ref.z)/(amp_xprime);
	m[3][1] = zprime_pdux/amp_xprime;//(part_pdux.uz - part_ref.uz)/(amp_xprime);
	
	m[0][2] = x_pdz/amp_z;     //(part_pdz.x  - part_ref.x)/(amp_z);
	m[1][2] = xprime_pdz/amp_z;//(part_pdz.ux - part_ref.ux)/(amp_z);
	m[2][2] = z_pdz/amp_z;     //(part_pdz.z  - part_ref.z)/(amp_z);
	m[3][2] = zprime_pdz/amp_z;//(part_pdz.uz - part_ref.uz)/(amp_z);
	
	m[0][3] = x_pduz/amp_zprime;     //(part_pduz.x  - part_ref.x)/(amp_zprime);
	m[1][3] = xprime_pduz/amp_zprime;//(part_pduz.ux - part_ref.ux)/(amp_zprime);
	m[2][3] = z_pduz/amp_zprime;     //(part_pduz.z  - part_ref.z)/(amp_zprime);
	m[3][3] = zprime_pduz/amp_zprime;//(part_pduz.uz - part_ref.uz)/(amp_zprime);
	
	if(doyouprintf == YES) {
		double det;
		det = matrix_det_4d(m, n);
		printf("\n--------------------------------- Get matrix first order -------------------------------------\n");
		for(i=0;i<n;i++) {
			for(j=0;j<n;j++) printf("%lf\t", m[i][j]);
			printf("\n");
		}
		printf("bias = %le\n", *band_bias);
		printf("determinant of m = %lf\n", det);
	}
	return TRUE;
}

extern int get_matrix_secondorder_4d(double m[4][4], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, double amp_z, double amp_zprime, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf)
//output the linear transfer matrix from the computation of the transfer matrix up to the second order.
//return FALSE if a problem had happened during the computation of matrix coefficients, else return TRUE.
//IMPORTANT!: transfer matrix calculated around the orbit of *reference (if you want *reference to be on the closed orbit, you have to make sure of that before)
//IMPORTANT2!:	if your xprime_ref ~= 0 at both entrance and exit faces of the latt, you are not concerned by this warning. 
//				Else, notice that this procedure works even when xprime_ref !=0, BUT it is then assumed 
//				that the particle goes straight on a distance <= band_bias around the latt boundary.
//				If fields are actually non-zero between +/- band_bias around the latt boundary, the result will be biased (more or less depending on "how much" the particle does not go straight over this distance).
//				Check the value of band_bias for safety. Must be "small" if you have a "significant" field around the boundary...
{
	int i,j,n=4;
	double xprime_ref, zprime_ref, uhp_ref, bias, x_pdx, xprime_pdx, z_pdx, zprime_pdx, x_mdx, xprime_mdx, z_mdx, zprime_mdx, x_pdxp, xprime_pdxp, z_pdxp, zprime_pdxp, x_mdxp, xprime_mdxp, z_mdxp, zprime_mdxp,
		x_pdz, xprime_pdz, z_pdz, zprime_pdz, x_mdz, xprime_mdz, z_mdz, zprime_mdz, x_pdzp, xprime_pdzp, z_pdzp, zprime_pdzp, x_mdzp, xprime_mdzp, z_mdzp, zprime_mdzp, uhp;
	//double not_symplec[4][4];
	struct Particle part_ref, part_pdx, part_mdx, part_pdxp, part_mdxp, part_pdz, part_mdz, part_pdzp, part_mdzp;
//	FILE *wfile;
//	wfile = fopen("data/transfermatrix.dat", "w");
	*band_bias = 0;
	
	for(i=0;i<n;i++) for(j=0;j<4;j++) m[i][j] = 0.;
	if(amp_x == 0 || amp_xprime == 0 || amp_z == 0 || amp_zprime == 0) errorstop("!!!ERROR: you call the procedure ''get_trmatrix'' with an invalid argument: amp_x/z(or)xprime/zprime = 0");
	
	//reference particle initial angle in the Cell local (entrance face) framework
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	if(cos(xprime_ref) == 0) {
		printf("!Warning in get_trmatrix: cos(xprime_ref) == 0, no transfer matrix will be computed\n");
		return FALSE;
	}
	uhp_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	zprime_ref = atan_ratio(reference->uz, uhp_ref);
	if(cos(zprime_ref) == 0) {
		printf("!Warning in get_trmatrix: cos(zprime_ref) == 0, no transfer matrix will be computed\n");
		return FALSE;
	}
		
	//init
	part_ref = *reference;
	part_ref.hat = -2; //no acceleration, no output in trackout.dat
	
	part_pdx = part_ref;
	part_pdx.x += amp_x/cos(xprime_ref);
	
	part_mdx = part_ref;
	part_mdx.x -= amp_x/cos(xprime_ref);
	
	part_pdxp= part_ref;
	part_pdxp.ux = sin(xprime_ref + amp_xprime)*cos(zprime_ref); 
	part_pdxp.uy = cos(xprime_ref + amp_xprime)*cos(zprime_ref); 
	part_pdxp.uz = sin(zprime_ref);
	
	part_mdxp = part_ref;
	part_mdxp.ux = sin(xprime_ref - amp_xprime)*cos(zprime_ref);
	part_mdxp.uy = cos(xprime_ref - amp_xprime)*cos(zprime_ref);
	part_mdxp.uz = sin(zprime_ref);
	
	part_pdz = part_ref;
	part_pdz.z += amp_z/cos(zprime_ref);
	
	part_mdz = part_ref;
	part_mdz.z -= amp_z/cos(zprime_ref);
	
	part_pdzp= part_ref;
	part_pdzp.ux = sin(xprime_ref)*cos(zprime_ref + amp_zprime);
	part_pdzp.uy = cos(xprime_ref)*cos(zprime_ref + amp_zprime);
	part_pdzp.uz = sin(zprime_ref + amp_zprime);
	
	part_mdzp = part_ref;
	part_mdzp.ux = sin(xprime_ref)*cos(zprime_ref - amp_zprime);
	part_mdzp.uy = cos(xprime_ref)*cos(zprime_ref - amp_zprime);
	part_mdzp.uz = sin(zprime_ref - amp_zprime);
	
	*band_bias = MAX(amp_z/cos(zprime_ref), amp_x/cos(xprime_ref));
	
	//transport the 5 particles
	(*transport_part)(&part_ref, latt,NULL);
	if(part_ref.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdx, latt,NULL);
	if(part_pdx.status != ALIVE) return FALSE;	
	(*transport_part)(&part_mdx, latt,NULL);
	if(part_mdx.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdxp, latt,NULL);
	if(part_pdxp.status != ALIVE) return FALSE;
	(*transport_part)(&part_mdxp, latt,NULL);
	if(part_mdxp.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdz, latt,NULL);
	if(part_pdz.status != ALIVE) return FALSE;	
	(*transport_part)(&part_mdz, latt,NULL);
	if(part_mdz.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdzp, latt,NULL);
	if(part_pdzp.status != ALIVE) return FALSE;
	(*transport_part)(&part_mdzp, latt,NULL);
	if(part_mdzp.status != ALIVE) return FALSE;
	
	//compute matrix coef.
	find_x(&x_pdx, &bias, &part_ref, &part_pdx);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdx, &bias, &part_ref, &part_mdx);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdxp, &bias, &part_ref, &part_pdxp);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdxp, &bias, &part_ref, &part_mdxp);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdz, &bias, &part_ref, &part_pdz);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdz, &bias, &part_ref, &part_mdz);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdzp, &bias, &part_ref, &part_pdzp);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdzp, &bias, &part_ref, &part_mdzp);
	*band_bias = MAX(*band_bias, bias);
	
	xprime_ref = atan_ratio(part_ref.ux, part_ref.uy);
	xprime_pdx = atan_ratio(part_pdx.ux, part_pdx.uy) - xprime_ref;
	xprime_mdx = atan_ratio(part_mdx.ux, part_mdx.uy) - xprime_ref;
	xprime_pdxp = atan_ratio(part_pdxp.ux, part_pdxp.uy) - xprime_ref;
	xprime_mdxp = atan_ratio(part_mdxp.ux, part_mdxp.uy) - xprime_ref;
	xprime_pdz = atan_ratio(part_pdz.ux, part_pdz.uy) - xprime_ref;
	xprime_mdz = atan_ratio(part_mdz.ux, part_mdz.uy) - xprime_ref;
	xprime_pdzp = atan_ratio(part_pdzp.ux, part_pdzp.uy) - xprime_ref;
	xprime_mdzp = atan_ratio(part_mdzp.ux, part_mdzp.uy) - xprime_ref;
	
	find_z(&z_pdx, &bias, &part_ref, &part_pdx);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_mdx, &bias, &part_ref, &part_mdx);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdxp, &bias, &part_ref, &part_pdxp);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_mdxp, &bias, &part_ref, &part_mdxp);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdz, &bias, &part_ref, &part_pdz);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_mdz, &bias, &part_ref, &part_mdz);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdzp, &bias, &part_ref, &part_pdzp);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_mdzp, &bias, &part_ref, &part_mdzp);
	*band_bias = MAX(*band_bias, bias);
	
	uhp_ref = sqrt(part_ref.ux*part_ref.ux + part_ref.uy*part_ref.uy);
	zprime_ref  = atan_ratio(part_ref.uz, uhp_ref);
	uhp = sqrt(part_pdx.ux*part_pdx.ux + part_pdx.uy*part_pdx.uy);
	zprime_pdx = atan_ratio(part_pdx.uz, uhp) - zprime_ref;
	uhp = sqrt(part_mdx.ux*part_mdx.ux + part_mdx.uy*part_mdx.uy);
	zprime_mdx = atan_ratio(part_mdx.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdxp.ux*part_pdxp.ux + part_pdxp.uy*part_pdxp.uy);
	zprime_pdxp = atan_ratio(part_pdxp.uz, uhp) - zprime_ref;
	uhp = sqrt(part_mdxp.ux*part_mdxp.ux + part_mdxp.uy*part_mdxp.uy);
	zprime_mdxp = atan_ratio(part_mdxp.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdz.ux*part_pdz.ux + part_pdz.uy*part_pdz.uy);
	zprime_pdz = atan_ratio(part_pdz.uz, uhp) - zprime_ref;
	uhp = sqrt(part_mdz.ux*part_mdz.ux + part_mdz.uy*part_mdz.uy);
	zprime_mdz = atan_ratio(part_mdz.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdzp.ux*part_pdzp.ux + part_pdzp.uy*part_pdzp.uy);
	zprime_pdzp = atan_ratio(part_pdzp.uz, uhp) - zprime_ref;
	uhp = sqrt(part_mdzp.ux*part_mdzp.ux + part_mdzp.uy*part_mdzp.uy);
	zprime_mdzp = atan_ratio(part_mdzp.uz, uhp) - zprime_ref;
	
	//second order canceled out
	m[0][0] = (x_pdx - x_mdx)/(2.*amp_x); 
	m[1][0] = (xprime_pdx - xprime_mdx)/(2.*amp_x); 
	m[2][0] = (z_pdx - z_mdx)/(2.*amp_x); 
	m[3][0] = (zprime_pdx - zprime_mdx)/(2.*amp_x); 
	
	m[0][1] = (x_pdxp - x_mdxp)/(2.*amp_xprime); 
	m[1][1] = (xprime_pdxp - xprime_mdxp)/(2.*amp_xprime);
	m[2][1] = (z_pdxp - z_mdxp)/(2.*amp_xprime); 
	m[3][1] = (zprime_pdxp - zprime_mdxp)/(2.*amp_xprime);
	
	m[0][2] = (x_pdz - x_mdz)/(2.*amp_z); 
	m[1][2] = (xprime_pdz - xprime_mdz)/(2.*amp_z); 
	m[2][2] = (z_pdz - z_mdz)/(2.*amp_z); 
	m[3][2] = (zprime_pdz - zprime_mdz)/(2.*amp_z); 
	
	m[0][3] = (x_pdzp - x_mdzp)/(2.*amp_zprime); 
	m[1][3] = (xprime_pdzp - xprime_mdzp)/(2.*amp_zprime);
	m[2][3] = (z_pdzp - z_mdzp)/(2.*amp_zprime); 
	m[3][3] = (zprime_pdzp - zprime_mdzp)/(2.*amp_zprime);
	
	//for(i=0;i<4;i++) for(j=0;j<4;j++) det_mat[i][j] = m[i][j];
	
	if(doyouprintf == YES) {
		double det;
		det = matrix_det_4d(m, n);
		for(i=0;i<n;i++) {
			for(j=0;j<n;j++) printf("%lf\t", m[i][j]);
			printf("\n");
		}
		printf("bias = %le\n", *band_bias);
		printf("determinant of m = %lf\n", det);
	}
	
	return TRUE;
}

extern int get_matrix_firstorder_5d(double m[5][5], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, double amp_z, double amp_zprime, double amp_dpovp, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf)
//output the linear 5D matrix from the computation of the transfer matrix to the first order.
//return FALSE if a problem had happened during the computation of matrix coefficients, else return TRUE.
//IMPORTANT!: transfer matrix calculated around the orbit of *reference (if you want *reference to be on the closed orbit, you have to make sure of that before)
//IMPORTANT2!:  if your xprime_ref ~= 0 at both entrance and exit faces of the latt, you are not concerned by this warning. 
//			    Else, notice that this procedure works even when xprime_ref !=0, BUT it is then assumed 
//				that the particle goes straight on a distance <= band_bias around the latt boundary.
//				If fields are actually non-zero between +/- band_bias around the latt boundary, the result will be biased (more or less depending on "how much" the particle does not go straight over this distance).
//				Check the value of band_bias for safety. Must be "small" if you have a "significant" field around the boundary...
{
	int i,j,n=5;
	double xprime_ref, zprime_ref, uhp, uhp_ref, bias;
	double x_pdx, x_pdux, x_pdz, x_pduz, xprime_pdx, xprime_pdux, xprime_pdz, xprime_pduz;
	double z_pdx, z_pdux, z_pdz, z_pduz, zprime_pdx, zprime_pdux, zprime_pdz, zprime_pduz;
	double x_pdp, xprime_pdp, z_pdp, zprime_pdp;
	struct Particle part_ref, part_pdx, part_pdux, part_pdz, part_pduz, part_pdp;
	
	for(i=1;i<n;i++) for(j=1;j<n;j++) m[i][j] = 0.;
	if(amp_x == 0 || amp_xprime == 0 || amp_z == 0 || amp_zprime == 0) errorstop("!!!ERROR: you call the procedure ''get_trmatrix_firstorder_5d'' with an invalid argument: amp_x/z(or)ux/uz = 0");
	
	*band_bias = 0;
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	uhp_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	zprime_ref = atan_ratio(reference->uz, uhp_ref);
	//init
	part_ref = *reference;
	//part_ref.hat = -2; //no acceleration, no output in trackout.dat
	part_pdx = part_ref;
	part_pdx.x += amp_x/cos(xprime_ref);
	
	part_pdux = part_ref;
	part_pdux.ux = sin(xprime_ref + amp_xprime)*cos(zprime_ref); 
	part_pdux.uy = cos(xprime_ref + amp_xprime)*cos(zprime_ref); 
	part_pdux.uz = sin(zprime_ref);
	
	part_pdz = part_ref;
	part_pdz.z += amp_z/cos(zprime_ref);
	
	part_pduz = part_ref;
	part_pduz.ux = sin(xprime_ref)*cos(zprime_ref+amp_zprime);
	part_pduz.uy = cos(xprime_ref)*cos(zprime_ref+amp_zprime);
	part_pduz.uz = sin(zprime_ref+amp_zprime);
	
	part_pdp = part_ref;
	part_pdp.brho *= (1 + amp_dpovp);
	
	*band_bias = MAX(amp_z/cos(zprime_ref), amp_x/cos(xprime_ref));
	
	//transport the 5 particles
	(*transport_part)(&part_ref, latt,NULL);
	if(part_ref.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdx, latt,NULL);
	if(part_pdx.status != ALIVE) return FALSE;	
	(*transport_part)(&part_pdux, latt,NULL);
	if(part_pdux.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdz, latt,NULL);
	if(part_pdz.status != ALIVE) return FALSE;	
	(*transport_part)(&part_pduz, latt,NULL);
	if(part_pduz.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdp, latt,NULL);
	if(part_pdp.status != ALIVE) return FALSE;
	
	find_x(&x_pdx, &bias, &part_ref, &part_pdx);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdux, &bias, &part_ref, &part_pdux);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdz, &bias, &part_ref, &part_pdz);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pduz, &bias, &part_ref, &part_pduz);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdp, &bias, &part_ref, &part_pdp); 
	*band_bias = MAX(*band_bias, bias);
	
	xprime_ref  = atan_ratio(part_ref.ux, part_ref.uy);
	xprime_pdx  = atan_ratio(part_pdx.ux, part_pdx.uy) - xprime_ref;
	xprime_pdux = atan_ratio(part_pdux.ux, part_pdux.uy) - xprime_ref;
	xprime_pdz  = atan_ratio(part_pdz.ux, part_pdz.uy) - xprime_ref;
	xprime_pduz = atan_ratio(part_pduz.ux, part_pduz.uy) - xprime_ref;
	xprime_pdp = atan_ratio(part_pdp.ux, part_pdp.uy) - xprime_ref; 
	
	find_z(&z_pdx, &bias, &part_ref, &part_pdx);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdux, &bias, &part_ref, &part_pdux);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdz, &bias, &part_ref, &part_pdz);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pduz, &bias, &part_ref, &part_pduz);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdp, &bias, &part_ref, &part_pdp); 
	*band_bias = MAX(*band_bias, bias);
	
	uhp_ref = sqrt(part_ref.ux*part_ref.ux + part_ref.uy*part_ref.uy);
	zprime_ref  = atan_ratio(part_ref.uz, uhp_ref);
	uhp = sqrt(part_pdx.ux*part_pdx.ux + part_pdx.uy*part_pdx.uy);
	zprime_pdx  = atan_ratio(part_pdx.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdux.ux*part_pdux.ux + part_pdux.uy*part_pdux.uy);
	zprime_pdux = atan_ratio(part_pdux.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdz.ux*part_pdz.ux + part_pdz.uy*part_pdz.uy);
	zprime_pdz  = atan_ratio(part_pdz.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pduz.ux*part_pduz.ux + part_pduz.uy*part_pduz.uy);
	zprime_pduz = atan_ratio(part_pduz.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdp.ux*part_pdp.ux + part_pdp.uy*part_pdp.uy);
	zprime_pdp = atan_ratio(part_pdp.uz, uhp) - zprime_ref; 
	
	
	m[0][0] = x_pdx/amp_x;
	m[1][0] = xprime_pdx/amp_x;
	m[2][0] = z_pdx/amp_x;
	m[3][0] = zprime_pdx/amp_x;
	m[4][0] = 0;
	
	m[0][1] = x_pdux/amp_xprime;
	m[1][1] = xprime_pdux/amp_xprime;
	m[2][1] = z_pdux/amp_xprime;
	m[3][1] = zprime_pdux/amp_xprime;
	m[4][1] = 0;
	
	m[0][2] = x_pdz/amp_z;
	m[1][2] = xprime_pdz/amp_z;
	m[2][2] = z_pdz/amp_z;
	m[3][2] = zprime_pdz/amp_z;
	m[4][2] = 0;
	
	m[0][3] = x_pduz/amp_zprime;
	m[1][3] = xprime_pduz/amp_zprime;
	m[2][3] = z_pduz/amp_zprime;
	m[3][3] = zprime_pduz/amp_zprime;
	m[4][3] = 0;
	
	m[0][4] = x_pdp/amp_dpovp;
	m[1][4] = xprime_pdp/amp_dpovp;
	m[2][4] = z_pdp/amp_dpovp;
	m[3][4] = zprime_pdp/amp_dpovp;
	m[4][4] = 1;
	
	if(doyouprintf == YES) {
		double det;
		det = matrix_det_5d(m, n);
		printf("\n--------------------------------- Get matrix first order -------------------------------------\n");
		for(i=0;i<n;i++) {
			for(j=0;j<n;j++) printf("%lf\t", m[i][j]);
			printf("\n");
		}
		printf("bias = %le\n", *band_bias);
		printf("determinant of m = %lf\n", det);
	}
	return TRUE;
}

extern int get_matrix_secondorder_5d(double m[5][5], double *band_bias, struct Particle *reference, double amp_x, double amp_xprime, double amp_z, double amp_zprime, double amp_dpovp, struct Lattice *latt, void(*transport_part)(struct Particle*, struct Lattice*, char*), int doyouprintf)
//output the linear 5D transfer matrix from the computation of the transfer matrix up to the second order.
//return FALSE if a problem had happened during the computation of matrix coefficients, else return TRUE.
//IMPORTANT!: transfer matrix calculated around the orbit of *reference (if you want *reference to be on the closed orbit, you have to make sure of that before)
//IMPORTANT2!:	if your xprime_ref ~= 0 at both entrance and exit faces of the latt, you are not concerned by this warning. 
//				Else, notice that this procedure works even when xprime_ref !=0, BUT it is then assumed 
//				that the particle goes straight on a distance <= band_bias around the latt boundary.
//				If fields are actually non-zero between +/- band_bias around the latt boundary, the result will be biased (more or less depending on "how much" the particle does not go straight over this distance).
//				Check the value of band_bias for safety. Must be "small" if you have a "significant" field around the boundary...
{
	int i,j,n=5;
	double xprime_ref, zprime_ref, uhp_ref, bias, x_pdx, xprime_pdx, z_pdx, zprime_pdx, x_mdx, xprime_mdx, z_mdx, zprime_mdx, x_pdxp, xprime_pdxp, z_pdxp, zprime_pdxp, x_mdxp, xprime_mdxp, z_mdxp, zprime_mdxp,
		x_pdz, xprime_pdz, z_pdz, zprime_pdz, x_mdz, xprime_mdz, z_mdz, zprime_mdz, x_pdzp, xprime_pdzp, z_pdzp, zprime_pdzp, x_mdzp, xprime_mdzp, z_mdzp, zprime_mdzp, uhp;
	double x_pdp, xprime_pdp, z_pdp, zprime_pdp, x_mdp, xprime_mdp, z_mdp, zprime_mdp;
	struct Particle part_ref, part_pdx, part_mdx, part_pdxp, part_mdxp, part_pdz, part_mdz, part_pdzp, part_mdzp, part_pdp, part_mdp;
//	FILE *wfile;
//	wfile = fopen("data/transfermatrix.dat", "w");
	*band_bias = 0;
	
	for(i=0;i<n;i++) for(j=0;j<4;j++) m[i][j] = 0.;
	if(amp_x == 0 || amp_xprime == 0 || amp_z == 0 || amp_zprime == 0) errorstop("!!!ERROR: you call the procedure ''get_trmatrix'' with an invalid argument: amp_x/z(or)xprime/zprime = 0");
	
	//reference particle initial angle in the Cell local (entrance face) framework
	xprime_ref = atan_ratio(reference->ux, reference->uy);
	if(cos(xprime_ref) == 0) {
		printf("!Warning in get_trmatrix: cos(xprime_ref) == 0, no transfer matrix will be computed\n");
		return FALSE;
	}
	uhp_ref = sqrt(pow(reference->ux, 2) + pow(reference->uy, 2));
	zprime_ref = atan_ratio(reference->uz, uhp_ref);
	if(cos(zprime_ref) == 0) {
		printf("!Warning in get_trmatrix: cos(zprime_ref) == 0, no transfer matrix will be computed\n");
		return FALSE;
	}
		
	//init
	part_ref = *reference;
	part_ref.hat = -2; //no acceleration, no output in trackout.dat
	
	part_pdx = part_ref;
	part_pdx.x += amp_x/cos(xprime_ref);
	
	part_mdx = part_ref;
	part_mdx.x -= amp_x/cos(xprime_ref);
	
	part_pdxp= part_ref;
	part_pdxp.ux = sin(xprime_ref + amp_xprime)*cos(zprime_ref); 
	part_pdxp.uy = cos(xprime_ref + amp_xprime)*cos(zprime_ref); 
	part_pdxp.uz = sin(zprime_ref);
	
	part_mdxp = part_ref;
	part_mdxp.ux = sin(xprime_ref - amp_xprime)*cos(zprime_ref);
	part_mdxp.uy = cos(xprime_ref - amp_xprime)*cos(zprime_ref);
	part_mdxp.uz = sin(zprime_ref);
	
	part_pdz = part_ref;
	part_pdz.z += amp_z/cos(zprime_ref);
	
	part_mdz = part_ref;
	part_mdz.z -= amp_z/cos(zprime_ref);
	
	part_pdzp= part_ref;
	part_pdzp.ux = sin(xprime_ref)*cos(zprime_ref + amp_zprime);
	part_pdzp.uy = cos(xprime_ref)*cos(zprime_ref + amp_zprime);
	part_pdzp.uz = sin(zprime_ref + amp_zprime);
	
	part_mdzp = part_ref;
	part_mdzp.ux = sin(xprime_ref)*cos(zprime_ref - amp_zprime);
	part_mdzp.uy = cos(xprime_ref)*cos(zprime_ref - amp_zprime);
	part_mdzp.uz = sin(zprime_ref - amp_zprime);
	
	part_pdp = part_ref;
	part_pdp.brho *= (1 + amp_dpovp);
	
	part_mdp = part_ref;
	part_mdp.brho *= (1 - amp_dpovp);
	
	*band_bias = MAX(amp_z/cos(zprime_ref), amp_x/cos(xprime_ref));
	
	//transport the 5 particles
	(*transport_part)(&part_ref, latt,NULL);
	if(part_ref.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdx, latt,NULL);
	if(part_pdx.status != ALIVE) return FALSE;	
	(*transport_part)(&part_mdx, latt,NULL);
	if(part_mdx.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdxp, latt,NULL);
	if(part_pdxp.status != ALIVE) return FALSE;
	(*transport_part)(&part_mdxp, latt,NULL);
	if(part_mdxp.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdz, latt,NULL);
	if(part_pdz.status != ALIVE) return FALSE;	
	(*transport_part)(&part_mdz, latt,NULL);
	if(part_mdz.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdzp, latt,NULL);
	if(part_pdzp.status != ALIVE) return FALSE;
	(*transport_part)(&part_mdzp, latt,NULL);
	if(part_mdzp.status != ALIVE) return FALSE;
	(*transport_part)(&part_pdp, latt,NULL);
	if(part_pdp.status != ALIVE) return FALSE;
	(*transport_part)(&part_mdp, latt,NULL);
	if(part_mdp.status != ALIVE) return FALSE;
	
	//compute matrix coef.
	find_x(&x_pdx, &bias, &part_ref, &part_pdx);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdx, &bias, &part_ref, &part_mdx);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdxp, &bias, &part_ref, &part_pdxp);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdxp, &bias, &part_ref, &part_mdxp);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdz, &bias, &part_ref, &part_pdz);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdz, &bias, &part_ref, &part_mdz);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdzp, &bias, &part_ref, &part_pdzp);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdzp, &bias, &part_ref, &part_mdzp);
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_pdp, &bias, &part_ref, &part_pdp); 
	*band_bias = MAX(*band_bias, bias);
	find_x(&x_mdp, &bias, &part_ref, &part_mdp); 
	*band_bias = MAX(*band_bias, bias);
	
	xprime_ref = atan_ratio(part_ref.ux, part_ref.uy);
	xprime_pdx = atan_ratio(part_pdx.ux, part_pdx.uy) - xprime_ref;
	xprime_mdx = atan_ratio(part_mdx.ux, part_mdx.uy) - xprime_ref;
	xprime_pdxp = atan_ratio(part_pdxp.ux, part_pdxp.uy) - xprime_ref;
	xprime_mdxp = atan_ratio(part_mdxp.ux, part_mdxp.uy) - xprime_ref;
	xprime_pdz = atan_ratio(part_pdz.ux, part_pdz.uy) - xprime_ref;
	xprime_mdz = atan_ratio(part_mdz.ux, part_mdz.uy) - xprime_ref;
	xprime_pdzp = atan_ratio(part_pdzp.ux, part_pdzp.uy) - xprime_ref;
	xprime_mdzp = atan_ratio(part_mdzp.ux, part_mdzp.uy) - xprime_ref;
	xprime_pdp = atan_ratio(part_pdp.ux, part_pdp.uy) - xprime_ref; 
	xprime_mdp = atan_ratio(part_mdp.ux, part_mdp.uy) - xprime_ref; 
	
	
	find_z(&z_pdx, &bias, &part_ref, &part_pdx);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_mdx, &bias, &part_ref, &part_mdx);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdxp, &bias, &part_ref, &part_pdxp);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_mdxp, &bias, &part_ref, &part_mdxp);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdz, &bias, &part_ref, &part_pdz);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_mdz, &bias, &part_ref, &part_mdz);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdzp, &bias, &part_ref, &part_pdzp);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_mdzp, &bias, &part_ref, &part_mdzp);
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_pdp, &bias, &part_ref, &part_pdp); 
	*band_bias = MAX(*band_bias, bias);
	find_z(&z_mdp, &bias, &part_ref, &part_mdp); 
	*band_bias = MAX(*band_bias, bias);
	
	uhp_ref = sqrt(part_ref.ux*part_ref.ux + part_ref.uy*part_ref.uy);
	zprime_ref  = atan_ratio(part_ref.uz, uhp_ref);
	uhp = sqrt(part_pdx.ux*part_pdx.ux + part_pdx.uy*part_pdx.uy);
	zprime_pdx = atan_ratio(part_pdx.uz, uhp) - zprime_ref;
	uhp = sqrt(part_mdx.ux*part_mdx.ux + part_mdx.uy*part_mdx.uy);
	zprime_mdx = atan_ratio(part_mdx.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdxp.ux*part_pdxp.ux + part_pdxp.uy*part_pdxp.uy);
	zprime_pdxp = atan_ratio(part_pdxp.uz, uhp) - zprime_ref;
	uhp = sqrt(part_mdxp.ux*part_mdxp.ux + part_mdxp.uy*part_mdxp.uy);
	zprime_mdxp = atan_ratio(part_mdxp.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdz.ux*part_pdz.ux + part_pdz.uy*part_pdz.uy);
	zprime_pdz = atan_ratio(part_pdz.uz, uhp) - zprime_ref;
	uhp = sqrt(part_mdz.ux*part_mdz.ux + part_mdz.uy*part_mdz.uy);
	zprime_mdz = atan_ratio(part_mdz.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdzp.ux*part_pdzp.ux + part_pdzp.uy*part_pdzp.uy);
	zprime_pdzp = atan_ratio(part_pdzp.uz, uhp) - zprime_ref;
	uhp = sqrt(part_mdzp.ux*part_mdzp.ux + part_mdzp.uy*part_mdzp.uy);
	zprime_mdzp = atan_ratio(part_mdzp.uz, uhp) - zprime_ref;
	uhp = sqrt(part_pdp.ux*part_pdp.ux + part_pdp.uy*part_pdp.uy);
	zprime_pdp = atan_ratio(part_pdp.uz, uhp) - zprime_ref; 
	uhp = sqrt(part_mdp.ux*part_mdp.ux + part_mdp.uy*part_mdp.uy);
	zprime_mdp = atan_ratio(part_mdp.uz, uhp) - zprime_ref; 
	
	//second order canceled out
	m[0][0] = (x_pdx - x_mdx)/(2.*amp_x); 
	m[1][0] = (xprime_pdx - xprime_mdx)/(2.*amp_x); 
	m[2][0] = (z_pdx - z_mdx)/(2.*amp_x); 
	m[3][0] = (zprime_pdx - zprime_mdx)/(2.*amp_x); 
	m[4][0] = 0;
	
	m[0][1] = (x_pdxp - x_mdxp)/(2.*amp_xprime); 
	m[1][1] = (xprime_pdxp - xprime_mdxp)/(2.*amp_xprime);
	m[2][1] = (z_pdxp - z_mdxp)/(2.*amp_xprime); 
	m[3][1] = (zprime_pdxp - zprime_mdxp)/(2.*amp_xprime);
	m[4][1] = 0;
	
	m[0][2] = (x_pdz - x_mdz)/(2.*amp_z); 
	m[1][2] = (xprime_pdz - xprime_mdz)/(2.*amp_z); 
	m[2][2] = (z_pdz - z_mdz)/(2.*amp_z); 
	m[3][2] = (zprime_pdz - zprime_mdz)/(2.*amp_z); 
	m[4][2] = 0;
	
	m[0][3] = (x_pdzp - x_mdzp)/(2.*amp_zprime); 
	m[1][3] = (xprime_pdzp - xprime_mdzp)/(2.*amp_zprime);
	m[2][3] = (z_pdzp - z_mdzp)/(2.*amp_zprime); 
	m[3][3] = (zprime_pdzp - zprime_mdzp)/(2.*amp_zprime);
	m[4][3] = 0;
	
	m[0][4] = (x_pdp - x_mdp)/(2*amp_dpovp);
	m[1][4] = (xprime_pdp - xprime_mdp)/(2*amp_dpovp);
	m[2][4] = (z_pdp - z_mdp)/(2*amp_dpovp);
	m[3][4] = (zprime_pdp - zprime_mdp)/(2*amp_dpovp);
	m[4][4] = 1;
	
	if(doyouprintf == YES) {
		double det;
		det = matrix_det_5d(m, n);
		for(i=0;i<n;i++) {
			for(j=0;j<n;j++) printf("%lf\t", m[i][j]);
			printf("\n");
		}
		printf("bias = %le\n", *band_bias);
		printf("determinant of m = %lf\n", det);
	}
	
	return TRUE;
}

extern int symplectify_2d(double m[2][2], double symplec_m[2][2])
{
	int i,j,n=2; //careful, n must be even since I divide by 2 later....
	double s[n][n], v[n][n], imm[n][n],ipm[n][n],inv_ipm[n][n], det_ipm; //imm:I-M, ipm=I+M, inv_ipm=(I+M)^-1
	double w[n][n],ipsw[n][n],imsw[n][n],inv_imsw[n][n], det_imsw;//ipsw=I+SW, imsw=I-SW;
	
	for(i=0;i<n;i++) {
		for(j=0;j<n;j++) {
			s[i][j]=0;
			v[i][j]=0;
			imm[i][j]=-m[i][j];
			ipm[i][j]=m[i][j];
			if(i==j) {
				ipm[i][j]+=1;
				imm[i][j]+=1;
			}
		}
	}
	for(i=0;i<n/2;i++) {
		s[2*i][2*i+1]=-1;
		s[2*i+1][2*i]=1;
	}
	
	//inverse of ipm:
	//matrix_inverse_4d(ipm, inv_ipm,n);
	//det_ipm = ipm[0][0]*ipm[1][1] - ipm[1][0]*ipm[0][1];
	det_ipm = matrix_det_4d(ipm, n);
	if(det_ipm==0) return FALSE;
	inv_ipm[0][0]=ipm[1][1]/det_ipm;
	inv_ipm[0][1]=-ipm[0][1]/det_ipm;
	inv_ipm[1][0]=-ipm[1][0]/det_ipm;
	inv_ipm[1][1]=ipm[0][0]/det_ipm;
	
	//computation of V
	mmprod2(s, imm, v);
	mmprod2(v,inv_ipm,v);
	
	for(i=0;i<n;i++) {
		for(j=0;j<n;j++) {
			w[i][j]=0.5*(v[i][j]+v[j][i]);
		}
	}
	
	mmprod2(s,w,ipsw);
	for(i=0;i<n;i++) {
		for(j=0;j<n;j++) {
			imsw[i][j] = -ipsw[i][j];
		}
		ipsw[i][i]+=1;
		imsw[i][i]+=1;
	}
	
	//inverse of imsv:
	//matrix_inverse_4d(imsv, inv_imsv,n);
	det_imsw = matrix_det_4d(imsw, n);
	if(det_imsw==0) return FALSE;
	inv_imsw[0][0]= imsw[1][1]/det_imsw;
	inv_imsw[0][1]=-imsw[0][1]/det_imsw;
	inv_imsw[1][0]=-imsw[1][0]/det_imsw;
	inv_imsw[1][1]= imsw[0][0]/det_imsw;
	
	mmprod2(ipsw,inv_imsw,symplec_m);
	
	return TRUE;
}

//symplectification algotithm from McKay et al (EPAC2006, WEPCH152)
extern int symplectify_4d(double m[4][4], double symplec_m[4][4])
{
	int i,j,n=4; //careful, n must be even since I divide by 2 later....
	double s[n][n], v[n][n], imm[n][n],ipm[n][n],inv_ipm[n][n]; //imm:I-M, ipm=I+M, inv_ipm=(I+M)^-1
	double w[n][n],ipsw[n][n],imsw[n][n],inv_imsw[n][n]; //ipsw=I+SW, imsw=I-SW;
	//double ipsv[n][n],imsv[n][n],inv_imsv[n][n];//ipsv=I+SV, imsv=I-SV;
	double test[n][n];
	
	for(i=0;i<n;i++) {
		for(j=0;j<n;j++) {
			s[i][j]=0;
			ipm[i][j]=m[i][j];
			imm[i][j]=-m[i][j];
			if(i==j) {
				ipm[i][j]+=1;
				imm[i][j]+=1;
			}
		}
	}
	for(i=0;i<n/2;i++) {
		s[2*i][2*i+1]=-1;
		s[2*i+1][2*i]=1;
	}
	//printf_mat_4d_real("s:", s, "s", NULL);
	printf_mat_4d_real("M:", m, "M", NULL);
	printf("det M=%le\n", matrix_det_4d(m,4));
	//printf_mat_4d_real("I-M:", imm, "I-M", NULL);
	//printf("det I-M=%le\n", matrix_det_4d(imm,4));
	//printf_mat_4d_real("I+M:", ipm, "I+M", NULL);
	//printf("det I+M=%le\n", matrix_det_4d(ipm,4));
	
	if(matrix_inverse_4d(ipm, inv_ipm)!=TRUE) return FALSE; //inverse of ipm
	//printf_mat_4d_real("(I+M)^-1:", inv_ipm, "(I+M)^-1", NULL);
	//printf("det (I+M)^-1=%le\n", matrix_det_4d(inv_ipm,4));
	//mmprod4(ipm, inv_ipm, test);
	//printf_mat_4d_real("(I+M)*(I+M)^-1:", test, "(I+M)^-1", NULL);

	//computation of V
	mmprod4(s, imm, v);
	//printf_mat_4d_real("S(I-M):", v, "V", NULL);
	//printf("det S(I-M)=%le\n", matrix_det_4d(v,4));
	mmprod4(v,inv_ipm,v);
	//printf_mat_4d_real("V:", v, "V", NULL);
	//printf("det V=%le\n", matrix_det_4d(v,4));
	
	//mmprod4(imm, inv_ipm, v);
	//mmprod4(s, v, v);
	
	//computation of W
	//matrix_transpose_4d(v, v_tr, n);
	for(i=0;i<n;i++) {
		for(j=0;j<n;j++) {
			w[i][j]=0.5*(v[i][j]+v[j][i]);
		}
	}
	
	mmprod4(s,w,ipsw);
	//printf_mat_4d_real("SW:", ipsw, "SW", NULL);
	//printf("det SW=%le\n", matrix_det_4d(ipsw,4));
	for(i=0;i<n;i++) {
		for(j=0;j<n;j++) {
			imsw[i][j] = -ipsw[i][j];
		}
		ipsw[i][i]+=1;
		imsw[i][i]+=1;
	}
	//printf_mat_4d_real("I+SW:", ipsw, "I+SW", NULL);
	//printf("det I+SW=%le\n", matrix_det_4d(ipsw,4));
	//printf_mat_4d_real("I-SW:", imsw, "I-SW", NULL);
	//printf("det I-SW=%le\n", matrix_det_4d(imsw,4));
	if(matrix_inverse_4d(imsw, inv_imsw)!=TRUE) return FALSE; //inverse of imsw
	//printf_mat_4d_real("(I-SW)^-1:", inv_imsw, "(I-SW)^-1", NULL);
	//printf("det (I-SW)^-1=%le\n", matrix_det_4d(inv_imsw,4));
	mmprod4(imsw, inv_imsw, test);
	//printf_mat_4d_real("(I-SW)*(I-SW)^-1:", test, "(I+M)^-1", NULL);
	
	mmprod4(ipsw,inv_imsw,symplec_m);
	//printf_mat_4d_real("SYMPLEC:", symplec_m, "symp", NULL);
	//printf("det SYMP=%le\n", matrix_det_4d(symplec_m,4));
	
	
	return TRUE;
}
