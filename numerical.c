/*
 *  numerical.c
 *  ringdesign
 *
 *  Routines below are taken from "Numerical Recipes in C"
 *  wed site www.fizyka.umk.pl/nrbook/bookcpdf.html
 *  most of them have however been modified by T.P.
 * zgeev zgetrf, zgetri

 */

#include "numerical.h"

// ************************************************************************************ //
//									RK drivers 											//
// ************************************************************************************ //
extern void rkdrive(struct Particle *part, int *awd_answer, double stepsize,
	void (*derivs)(struct Particle*, double[], double[], struct Cell*, void(*)(double,double,double,double*,double*,double*,struct Cell*,int)),
	struct Cell* cell, char *txtfile, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i, plot3d=NO;
	double dfds[NVAR], bprobe[3];
	FILE *wtrackout, *wpickup;
	//printf("x=%lf,y=%lf\n",part->x,part->y);		

	if(part->status == ALIVE) {
		if(stepsize < TINYLENGTH) nrerror("Stepsize too small (< TINYLENGTH) in routine rkdrive"); //check that the stepsize is not infinitly small
		if(NVAR != 7) nrerror("!!!ERROR in rkdrive: NVAR != 7 while rkdrive is only implemented for one type of physics (tracking only in B field), this type of physics requires NVAR = 7!!!");
		//gestion de l'exception instru.thmax&ymax = 0
		if(cell->instrutype != NO && cell->instru.thmax == 0 && cell->instru.ymax == 0) {
			*awd_answer = cell->instrutype;
			if(cell->instrutype == PICKUP) {
				wpickup = fopen("data/pickup.dat","a");
				write_trackout(wpickup, part, bprobe);
				fclose(wpickup);
			}
			return;
		}
		if(txtfile != NULL) wtrackout = fopen(txtfile,"a");
		
		//loop until the particle reach the Cell boundary (taking at most MAXSTP steps)
		for(i = 0;i < MAXSTP;i++) {
			(*derivs)(part, dfds, bprobe, cell, add_contribution_comp);	//initialization of dfds and bprobe

			if(part->status != ALIVE) {
				//printf("part not alive, break\n");
				break; //if particle has hit a collimator, break
			}
			if(txtfile != NULL) write_trackout(wtrackout, part, bprobe); //if doyoutrackout == YES, write coordinates in trackout.dat
			//MUON RING
			if(part->hat == 3) compute_max_angle_muon_ring(part);
			rk4(part, dfds, stepsize, derivs, cell, add_contribution_comp);
			if(part->status != ALIVE) break; //if particle has hit a collimator, break
		
			*awd_answer  = arewedone(part, cell);
			
			if(*awd_answer == CUP || *awd_answer == PICKUP) { //Did we reach a pickup or a cup?
				//printf("cup or pickup?\n");
				moveback_totheboun(part,&(cell->instru),derivs,cell, add_contribution_comp); //moveback last point on the exact Cell edge
				(*derivs)(part, dfds, bprobe, cell, add_contribution_comp); //get bprobe for the next output in trackout.dat
				break;
			}
			else if(*awd_answer == YES) { //Did we reach the cell boundary?
				//printf("boundary reached\n");
				if(i == 0) errorstop("!!!ERROR in rkdrive: Cell boundary reached after only on step!!!\nyou may have done a mistake in the geometrical parameter(s) of at least one lattice Cell");
				moveback_totheboun(part,&(cell->boun),derivs,cell, add_contribution_comp); //moveback last point on the exact Cell edge
				(*derivs)(part, dfds, bprobe, cell, add_contribution_comp); //get bprobe for the next output in trackout.dat
				break;
			}
		}
		
		//terminate
		if(txtfile != NULL) {
			write_trackout(wtrackout, part, bprobe);
			if(plot3d==NO) fprintf(wtrackout,"\n"); //jump one line in trackout.dat
			fclose(wtrackout);
		}
		if(*awd_answer != YES && (part->hat*part->hat) == 1) {
			wpickup = fopen("data/pickup.dat","a");
			write_trackout(wpickup, part, bprobe);
			fclose(wpickup);
		}
		
		if(i+1 >= MAXSTP) nrerror("too many steps in routine rkdrive,\nplease increase MAXSTP or increase the step size");
	}
}

// ************************************************************************************ //
//									simple RK4 routine									//
// ************************************************************************************ //
static int rk4(struct Particle *part, double dfds[], double h,
	void derivs(struct Particle*, double[], double[], struct Cell*, void(*)(double,double,double,double*,double*,double*,struct Cell*,int)), struct Cell* cell,
	void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int i;
	double hh, h6, dfm[NVAR], dft[NVAR], fake_bprobe[3];
	struct Particle ft;
		
	//initialization:
	hh		= h*0.5;
	h6		= h/6.0;
	
	//Runge-Kutta-4 integation step:
	//first step:
	ft = *part;
	ft.x	= part->x	+ hh*dfds[0];
	ft.y	= part->y	+ hh*dfds[1];
	ft.z	= part->z	+ hh*dfds[2];
	ft.ux	= part->ux	+ hh*dfds[3];
	ft.uy	= part->uy	+ hh*dfds[4];
	ft.uz	= part->uz	+ hh*dfds[5];
	ft.t	= part->t	+ hh*dfds[6];
	//second step
	//if(write_enge==YES) printf("sec step\n");
	derivs(&ft, dft, fake_bprobe, cell, add_contribution_comp);
	if(ft.status != ALIVE) {	//if ft has hit a collimator, *part is marked as lost, then return.
		part->status = ft.status;
		return 1;
	}
	ft.x	= part->x	+ hh*dft[0];
	ft.y	= part->y	+ hh*dft[1];
	ft.z	= part->z	+ hh*dft[2];
	ft.ux	= part->ux	+ hh*dft[3];
	ft.uy	= part->uy	+ hh*dft[4];
	ft.uz	= part->uz	+ hh*dft[5];
	ft.t	= part->t	+ hh*dft[6];
	//third step
	//if(write_enge==YES) printf("third step\n");
	derivs(&ft, dfm, fake_bprobe, cell, add_contribution_comp);
	if(ft.status != ALIVE) {	//if ft has hit a collimator, *part is marked as lost, then return.
		part->status = ft.status;
		return 1;
	}
	ft.x	= part->x	+ h*dfm[0];
	ft.y	= part->y	+ h*dfm[1];
	ft.z	= part->z	+ h*dfm[2];
	ft.ux	= part->ux	+ h*dfm[3];
	ft.uy	= part->uy	+ h*dfm[4];
	ft.uz	= part->uz	+ h*dfm[5];
	ft.t	= part->t	+ h*dfm[6];
	for(i=0;i<=NVAR-1;i++) dfm[i] += dft[i];
	//fourth step
	//if(write_enge==YES) printf("4th step\n");
	derivs(&ft, dft, fake_bprobe, cell, add_contribution_comp);
	if(ft.status != ALIVE) {	//if ft has hit a collimator, *part is marked as lost, then return.
		part->status = ft.status;
		return 1;
	}
	part->x		= part->x+h6*(dfds[0]+dft[0]+2.0*dfm[0]);//accumulate increments with proper weights.
	part->y		= part->y+h6*(dfds[1]+dft[1]+2.0*dfm[1]);
	part->z		= part->z+h6*(dfds[2]+dft[2]+2.0*dfm[2]);
	part->ux	= part->ux+h6*(dfds[3]+dft[3]+2.0*dfm[3]);
	part->uy	= part->uy+h6*(dfds[4]+dft[4]+2.0*dfm[4]);
	part->uz	= part->uz+h6*(dfds[5]+dft[5]+2.0*dfm[5]);
	part->t		= part->t+h6*(dfds[6]+dft[6]+2.0*dfm[6]);
	
	part->s += h;	//increment part->s of the step-size h
	return 0;
}

// ************************************************************************************ //
//							"arewedone?" routines										//
// ************************************************************************************ //
//test if the particle reached the Cell boundary
static int arewedone(struct Particle *part, struct Cell *cell)
{
	double th;
	
	th = atan_ratio(part->y,part->x);
	if(th<0 && cell->boun.thmax > 1.5*PI) th += 2.*PI;
	
	
	//did you reach a pickup or a cup?
	if(cell->instrutype != NO) {
		if(cell->instru.thmax > 0 && th >= cell->instru.thmax) return cell->instrutype;
		if(cell->instru.ymax > 0 && part->y >= cell->instru.ymax) return cell->instrutype;	
	}
	
	//did you reach cell boundary?
	if(cell->boun.thmax > 0) if(th >= cell->boun.thmax) return YES; //radial Cell, test if th >= th boundary
	if(cell->boun.ymax > 0 && part->y >= cell->boun.ymax) return YES; //straight Cell, test if y >= Cell length
	
	return NO;
}

static void moveback_totheboun(struct Particle *part, struct Boundary *theboun, void (*derivs)(struct Particle*, double[], double[], struct Cell*,
	void(*)(double,double,double,double*,double*,double*,struct Cell*,int)), struct Cell *cell,
	void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int j, maxtry;
	double th, r, d, dfds[NVAR], fake_bprobe[NVAR];
	
	maxtry = 100; //arbitrary chosen, but you must not need more than 10 steps to get close enough to the boundary
	
	for(j = 1; j <= maxtry; j++) {
	
		if(theboun->thmax > 0) {
			th = atan_ratio(part->y,part->x);
			if(th<0 && cell->boun.thmax > 1.5*PI) th += 2.*PI;
			r = sqrt(part->x*part->x + part->y*part->y);
			d = r*(theboun->thmax - th);
		}
		else if(theboun->ymax > 0) {
			d = theboun->ymax - part->y;
		}
		else nrerror("!!! ERROR in moveback_lastpoint: unexpected case: both boundary conditions are = 0!!!");
		
		if (fabs(d) < TINYLENGTH) return; //if d is small enougth, we reached the boundary!! retrun
		
		//else do one more rk4 step with step size = d (can be ether > 0 or < 0)
		(*derivs)(part, dfds, fake_bprobe, cell, add_contribution_comp);
		rk4(part, dfds, d, derivs, cell, add_contribution_comp);
	}
	
	if(j >= maxtry) { //failed, display error message for debug
		printf("!!! ERROR in moveback_totheboun: cannot set particle on the the boundary (either Cell boundary or instru (cup of pickup))\n");
		printf("after %i steps in moveback_lastpoint, the particle is still further than TINYLENGTH=%le [m] from the boundary (d = %le [m])\n", maxtry, TINYLENGTH, d);
		printf("possible solution: improve the computation accuracy of increase TINYLENGTH (see constant.h)\n");
		nrerror("ERROR in moveback_totheboun, see above for details");
	}
}

// ************************************************************************************ //
//							write in trackout.dat routine								//
// ************************************************************************************ //
static void write_trackout(FILE *wfile, struct Particle *part, double bprobe[])
{
	//int mod;
	double x, y, r, th, br, bth, bshinjix, bshinjiy;
	double th_glob;
	//double phi, phimod, frf;
	
		if(part->status == ALIVE) {
			x = part->x;
			y = part->y;
			//frf = 74.991737e06;
			//phi = part->t*360.*frf;
			//mod = (int)(phi/360.);
			//phimod = phi - mod*360.;
			fwk_pt_loctoglob(&x, &y, &(part->fwk));
			r = sqrt(x*x+y*y);
			th=atan_ratio(y,x);
			br = bprobe[0]*cos(th) + bprobe[1]*sin(th);
			bth = -bprobe[0]*sin(th) + bprobe[1]*cos(th);
			bshinjix = bprobe[0]*cos(4.5*PI/180.) + bprobe[1]*sin(4.5*PI/180.);
			bshinjiy = -bprobe[0]*sin(4.5*PI/180.) + bprobe[1]*cos(4.5*PI/180.);
			th_glob = th - tan(30.*PI/180.)*log(r/4.);
			//bshinjix = 0;
			//bshinjiy = 0;
			//MUON RING :
			//if(max_angle < fabs(part->ux)) max_angle = fabs(part->ux);
			//fprintf(wfile, "%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le\n", part->s, x, y, part->z, bprobe[0], bprobe[1], bprobe[2], part->ux, part->uy, part->uz, sqrt(part->brho*part->brho*89880 + 0.261121)-0.511, part->t, phimod);
			//fprintf(wfile, "%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %.10f\t %le\t %le\t %le\n", part->s, x, y, part->z, bprobe[0], bprobe[1], bprobe[2], part->ux, part->uy, part->uz, part->brho, part->t, r, th*180./PI, br, bth);
			//fprintf(wfile, "%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %.10f\t %le\t %le\t %le\t %le\t %le\n", part->s, x, y, part->z, bprobe[0], bprobe[1], bprobe[2], part->ux, part->uy, part->uz, part->brho, part->t, r, th*180./PI, th_glob*180./PI, bshinjiy, br, bth);
			fprintf(wfile, "%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %.10f\t %le\t %le\n", part->s, x, y, part->z, bprobe[0], bprobe[1], bprobe[2], part->ux, part->uy, part->uz, r, th*180./PI, th_glob*180./PI, br, bth);
			//fprintf(wfile, "%le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le \t %le\n", part->s, x, y, part->z, bprobe[0], bprobe[1], bprobe[2], part->ux, part->uy, part->uz, part->brho);
			//fprintf(wfile, "%i  %.3e \t %.3e \t %.3e \t %.3e \t %.3e \t %.3e \t %.3e\n", ao_partnb, (x-36.15)*1000., 1000.*part->z, y*1000., part->brho*299.792458*part->ux, part->brho*299.792458*part->uz, part->brho*299.792458*part->uy, 1.e9*part->t);
		}
}

// ************************************************************************************ //
//									Fast Fourier Transform								//
// ************************************************************************************ //
// start data from i=1, not normalised.
extern void four1(double data[], unsigned long nn, int isign)
{
    unsigned long n,mmax,m,j,istep,i; 
    double wtemp,wr,wpr,wpi,wi,theta; 
    double tempr,tempi; 
    n=nn << 1; 
    j=1; 
    for (i=1;i<n;i+=2) { 
        if (j> i) { 
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]); 
        } 
        m=nn; 
        while (m >= 2 && j> m) { 
            j-= m; 
            m>>= 1; 
        } 
        j+= m; 
    } 
    mmax=2; 
    while (n> mmax) { 
        istep=mmax << 1; 
        theta=isign*(6.28318530717959/mmax); 
        wtemp=sin(0.5*theta); 
        wpr =-2.0*wtemp*wtemp;
        
        
        wpi=sin(theta); 
        wr=1.0; 
        wi=0.0; 
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) { 
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1]; 
                tempi=wr*data[j+1]+wi*data[j]; 
                data[j]=data[i]-tempr; 
                data[j+1]=data[i+1]-tempi; 
                data[i] +=tempr; 
                data[i+1] += tempi; 
            } 
            wr=(wtemp=wr)*wpr-wi*wpi+wr; 
            wi=wi*wpr+wtemp*wpi+wi; 
        } 
        mmax=istep; 
    } 
}

// ************************************************************************************ //
//									Various math functions								//
// ************************************************************************************ //

//return -1 if a < 0, 1 if > 0, 0 if a = 0;
extern int sign(double a)
{
	if(a < 0) return -1;
	else if(a > 0) return 1;
	else if (a == 0) return 0;
	else {
		errorstop("unexpeced value in sign function");
		return 999; //never happen!
	}
}

//get value of atan(y/x), and return PI/2 if x = 0
extern double atan_ratio(double y, double x)
{
	//printf("in atan_ratio: x = %le y = %le\n", x, y);
	if(x < 0) return PI+atan(y/x);
	if(x > 0) return atan(y/x);
	else if(x == 0 && y > 0) return PI/2.;
	else if(x == 0 && y < 0) return 3*PI/2.;
	else if(x == 0 && y == 0) {
		printf("!WARNING in atan_ratio: x and y are null, 0 is returned!\n");
		return 0;
	}
	else {
		printf("\nIn atan_ratio, ERROR: x = %le y = %le\n", x, y);
		errorstop("ERROR in atan_ratio: unexpected parameter value(s)");
		return 999;
	}
}

// Evaluate Bessel function of first kind and order 0 at input x (from Kapteyn Institute Groningen)
extern double bessj0(double x)
{
	double ax,z;
	double xx,y,ans,ans1,ans2;
	
	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	}
	else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3+y*(-0.6911147651e-5+y*(0.7621095161e-6-y*0.934935152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}

//Evaluate Bessel function of first kind and order 1 at input x (from Kapteyn Institute Groningen)
extern double bessj1( double x )
{
	double ax,z;
	double xx,y,ans,ans1,ans2;
	
	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74+y*(99447.43394+y*(376.9991397+y*1.0))));
		ans=ans1/ans2;
	}
	else {
		z=8.0/ax;
		y=z*z;
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}

// Evaluate Bessel function of second kind and order 0 at input x (from Kapteyn Institute Groningen)
extern double bessy0(double x)
{
	double z;
	double xx,y,ans,ans1,ans2;
	
	if (x < 8.0) {
		y=x*x;
		ans1 = -2957821389.0+y*(7062834065.0+y*(-512359803.6+y*(10879881.29+y*(-86327.92757+y*228.4622733))));
		ans2=40076544269.0+y*(745249964.8+y*(7189466.438+y*(47447.26470+y*(226.1030244+y*1.0))));
		ans=(ans1/ans2)+0.636619772*bessj0(x)*log(x);
	}
	else {
		z=8.0/x;
		y=z*z;
		xx=x-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3+y*(-0.6911147651e-5+y*(0.7621095161e-6+y*(-0.934945152e-7))));
		ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
	}
	return ans;
}

// Evaluate Bessel function of second kind and order 1 at input x (from Kapteyn Institute Groningen)
extern double bessy1(double x)
{
	double z;
	double xx,y,ans,ans1,ans2;
	
	if (x < 8.0) {
		y=x*x;
		ans1=x*(-0.4900604943e13+y*(0.1275274390e13+y*(-0.5153438139e11+y*(0.7349264551e9+y*(-0.4237922726e7+y*0.8511937935e4)))));
		ans2=0.2499580570e14+y*(0.4244419664e12+y*(0.3733650367e10+y*(0.2245904002e8+y*(0.1020426050e6+y*(0.3549632885e3+y)))));
		ans=(ans1/ans2)+0.636619772*(bessj1(x)*log(x)-1.0/x);
	} else {
		z=8.0/x;
		y=z*z;
		xx=x-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
		ans=sqrt(0.636619772/x)*(sin(xx)*ans1+z*cos(xx)*ans2);
	}
	return ans;
}

// will not work for negative x and a non integer negative alpha!
extern double bessel_j(double alpha, double x, int n)
{
	int i,k;
	double alp, fact, gam, res=0;
	if(x<0 && alpha<0 && fabs(round(alpha)-alpha)>TINYDIMLESS) errorstop(" bessel_j does not work if x<0 and alpha<0, not integer");
	if(alpha<0 && fabs(round(alpha)-alpha)<TINYDIMLESS) alp = -alpha;
	else alp=alpha;
	for(k=0;k<n;k++) {
		fact=1;
		for(i=1;i<k+1;i++) fact*=i;
		gam = tgamma(alp+k+1);
		if(gam!=gam || fact!=fact) {
			printf("Bessel 1st kind, order %i out of range\n", k);
			break;
		}
		if(x<0 && alpha<0) res += pow(x/2.,(int) (alp)+2*k)*pow(-1,k)/(fact*gam);
		else res += pow(x/2.,alp+2.*k)*pow(-1,k)/(fact*gam);
	}
	return res;
}

extern double bessel_y(double alpha, double x, int n)
{
	int i,j, alp;
	double sum, fact,res;
	
	if(fabs(round(alpha)-alpha)>TINYDIMLESS) return (bessel_j(alpha,x,n)*cos(alpha*PI)-bessel_j(-alpha,x,n))/(sin(alpha*PI));
	else {
		alp= (int) round(fabs(alpha));
		res = 2*bessel_j(alp,x,n)*(log(x/2.)+EULERGAMMA);
		for(i=0;i<alp;i++) res -= (alp-i-1)*pow(x/2.,2*i-alp);
		for(i=0;i<n;i++) {
			sum=0;
			for(j=1;j<i+1;j++) sum+=2/j;
			for(j=i+1;j<alp+i+1;j++) sum+=1/j;
			fact=1;
			for(j=1;j<i+1;j++) fact*=j;
			for(j=1;j<alp+i+1;j++) fact*=j;
			res -= pow(-1,i)*pow(x/2.,2*i+alp)*sum/fact;
		}
		res /=PI;
		if(alpha<0) res *=pow(-1,-alp);
		return res;
	}
}

extern int overflow_test(double a)
{
	if(a + 1 != a) return TRUE;
	else return FALSE;
}

extern int sanity_test_number(double a)
{
	if(overflow_test(a)==FALSE || isfinite(a)==0) return FALSE;
	else return TRUE;
}


extern int factorial(int n)
{
	int i, fact=1;
	
	if(n<0) return 0;
	else if(n==0) return 1;
	else {
		for(i=1;i<n+1;i++) {
			fact *=i;
		}
		return fact;
	}
}
// ************************************************************************************ //
//									Matrix manipulation									//
// ************************************************************************************ //

// Real Matrix-Matrix multiplication 2D (mul=A.mul also works)
extern void mmprod2(double fst[2][2], double sec[2][2], double mul[2][2])
{
	int c,d,k,n=2;
	double temp[n][n];
	
	for (c = 0; c < n; c++) {
		for (d = 0; d < n; d++) {
			temp[c][d] = 0;
			for (k = 0; k < n; k++) {
				temp[c][d] += fst[c][k] * sec[k][d];
			}
		}
	}
	for (c = 0; c < n; c++) {
		for (d = 0; d < n; d++) {
			mul[c][d]=temp[c][d];
		}
	}
}

// Real Matrix-Matrix multiplication 4D (mul=A.mul also works)
extern void mmprod4(double fst[4][4], double sec[4][4], double mul[4][4])
{
	int c,d,k,n=4;
	double temp[n][n];
	
	for (c = 0; c < n; c++) {
		for (d = 0; d < n; d++) {
			temp[c][d] = 0;
			for (k = 0; k < n; k++) {
				temp[c][d] += fst[c][k] * sec[k][d];
			}
		}
	}
	for (c = 0; c < n; c++) {
		for (d = 0; d < n; d++) {
			mul[c][d]=temp[c][d];
		}
	}
}

// Matrix-vector 3D product V=A.U (V=A.V works also)
extern void mvprod3(double v[], double a[][3], double u[])
{
	int i, j, n=3;
	double tempo_v[n];
	
	for(i = 0; i< n; i++) {
		tempo_v[i] = 0;
		for(j = 0; j< n; j++) {
			tempo_v[i] += a[i][j]*u[j]; 
		}
	}
	
	for(i = 0; i< n; i++) v[i] = tempo_v[i];
}

// Matrix-vector 5D product V=A.U (V=A.V works also)
extern void mvprod4(double v[], double a[][4], double u[])
{
	int i, j, n=4;
	double tempo_v[n];
	
	for(i = 0; i< n; i++) {
		tempo_v[i] = 0;
		for(j = 0; j< n; j++) {
			tempo_v[i] += a[i][j]*u[j]; 
		}
	}
	
	for(i = 0; i< n; i++) v[i] = tempo_v[i];
}

// Matrix-vector 5D product V=A.U (V=A.V works also)
extern void mvprod5(double v[], double a[][5], double u[])
{
	int i, j, n=5;
	double tempo_v[n];
	
	for(i = 0; i< n; i++) {
		tempo_v[i] = 0;
		for(j = 0; j< n; j++) {
			tempo_v[i] += a[i][j]*u[j]; 
		}
	}
	
	for(i = 0; i< n; i++) v[i] = tempo_v[i];
}

//This function multiplies the entered matrices mult[r1][c2] = a[r1][c1] x b[r2][c2] with 1 dimension matrix
extern void matrix_dot_1dim(double complex a[], double complex b[],double complex mult[],int r1,int c1,int r2,int c2)
{
	int i,j,k;
	double complex temp[r1*c2];
	
	if(c1!=r2) errorstop("cannot multiply matrices");
	for(i=0; i<r1; i++) {
		for(j=0; j<c2; j++) {
			temp[i+c2*j]=0;
			for(k=0; k<c1; k++) temp[i+c2*j] += a[k+c1*j]*b[i+c2*k];
		}
	}
	for(i=0; i<r1; i++) for(j=0; j<c2; j++) mult[i+c2*j] = temp[i+c2*j];
}

// determinant of a real matrix 4D
extern double matrix_det_4d(double mat[4][4], int n)
{
	int c, subi, i, j, subj;
	//int pi,pj;
	double d=0;
	double submat[4][4];  

	if (n == 1) return (mat[0][0]);
	else if (n == 2) return(mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1]);
	else {
		for(c = 0; c < n; c++) {
			subi = 0;
			for(i = 1; i < n; i++) {
				subj = 0;
				for(j = 0; j < n; j++) {
					if (j != c) {
						submat[subi][subj] = mat[i][j];
						subj++;
					}
				}
				subi++;
			}
			//printf("\nd=%le + %lf*\t",d, mat[0][c]);
			//for(pi=0;pi<n-1;pi++) {
			//	for(pj=0;pj<n-1;pj++) printf("%lf\t",submat[pi][pj]);
			//	printf("\n\t\t\t\t");
			//}
			d = d + (pow(-1 ,c) * mat[0][c] * matrix_det_4d(submat, n - 1));
			//printf("\n");
		}
	}
	return d;
}

// determinant of a real matrix 4D
extern double matrix_det_5d(double mat[5][5], int n)
{
	int c, subi, i, j, subj;
	//int pi,pj;
	double d=0;
	double submat[5][5];  

	if (n == 1) return (mat[0][0]);
	else if (n == 2) return(mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1]);
	else {
		for(c = 0; c < n; c++) {
			subi = 0;
			for(i = 1; i < n; i++) {
				subj = 0;
				for(j = 0; j < n; j++) {
					if (j != c) {
						submat[subi][subj] = mat[i][j];
						subj++;
					}
				}
				subi++;
			}
			//printf("\nd=%le + %lf*\t",d, mat[0][c]);
			//for(pi=0;pi<n-1;pi++) {
			//	for(pj=0;pj<n-1;pj++) printf("%lf\t",submat[pi][pj]);
			//	printf("\n\t\t\t\t");
			//}
			d = d + (pow(-1 ,c) * mat[0][c] * matrix_det_5d(submat, n - 1));
			//printf("\n");
		}
	}
	return d;
}

//transpose of a real matrix 4D
extern void matrix_transpose_4d(double m[4][4], double tr[4][4], int r)
{
	int i, j;
	for(i = 0; i < r; i++) for (j = 0; j < r; j++) tr[i][j] = m[j][i];
}

//cofactor of a real matrix 4D
extern void matrix_cofactor_4d(double num[4][4], double fac[4][4], int f)
{
	int p, q, m, n, i, j;
	double b[4][4];
	
	for(q = 0; q < f; q++) {
		for(p = 0; p < f; p++) {
			m = 0;
			n = 0;
			for(i = 0; i < f; i++) {
				for(j = 0; j < f; j++) {
					if(i != q && j != p) {
						b[m][n] = num[i][j];
						if(n < (f - 2)) n++;
						else {
							n = 0;
							m++;
						}
					}
				}
			}
			fac[q][p] = pow(-1, q+p)*matrix_det_4d(b,f-1);
		}
	}
}

// inverse real matrix num, 4D or less (numbered by f)
extern int matrix_inverse_4d(double num[4][4], double inverse[4][4])
{
	int i, j,f=4;
	double d, tr[4][4], fac[4][4];
	
	d = matrix_det_4d(num, f);
	if(d==0) {
		printf("matrix cannot be inversed\n");
		return FALSE;
	}
	matrix_cofactor_4d(num, fac, f);
	matrix_transpose_4d(fac, tr, f);
	for (i = 0; i < f; i++) for (j = 0; j < f; j++) inverse[i][j] = tr[i][j]/d;
	return TRUE;
}

// inverse complex matrix of dimension n
extern void matrix_complex_inverse(double complex *A, double complex *invA, int n)
{
	int LWORK=10*n;
	int *permutations;
	double complex *WORK, *tempA;
	tempA = (double complex*) malloc( n*n*sizeof(double complex) );
	permutations = (int*) malloc( 2*n*sizeof(int) );
	WORK = (double complex *)malloc(LWORK*sizeof(double complex));
	int INFO;
	
	zge_transpose(tempA,A,n);
	
	zgetrf_(&n, &n, tempA , &n, permutations , &INFO);
	
	if(INFO != 0) {
		printf("ComplexMatrixInverse: Error at zgetrf  \n"); 
		exit(0);
	}

	zgetri_(&n, tempA , &n, permutations , WORK, &LWORK, &INFO);

	if (INFO != 0) {
		printf("ComplexMatrixInverse: Error at zgetri  \n"); 
		exit(0);
	}

	zge_transpose(invA,tempA,n);

	free(WORK);
	free(tempA);
	free(permutations);
}

//invmat calculates the inverse of matrix. The algorithm used is the
//Gauss-Jordan algorithm described in Stoer, Numerische matematik, 1 Teil.
extern int invmat(double matrix[5][5], int nfree)
{
	double even,hv[5],mjk,rowmax;
	int evin,i,j,k,per[5],row;

	for(i = 0; i < nfree; i++) per[i] = i;	// set permutation array 
	for(j = 0; j < nfree; j++) {
		rowmax = fabs( matrix[j][j] ); // in j-th column, determine row with largest element.
		row = j;
		for(i = j + 1; i < nfree; i++) {
			if(fabs( matrix[i][j] ) > rowmax) {
				rowmax = fabs( matrix[i][j] );
				row = i;
			}
		}
		if(matrix[row][j] == 0.0) {
			printf("determinant is zero!");
			return FALSE;
		}
		if(row > j) {	// if largest element not on diagonal, then permutate rows.
			for(k = 0; k < nfree; k++) {
				even = matrix[j][k];
				matrix[j][k] = matrix[row][k];
				matrix[row][k] = even;
			}
			evin = per[j];	// keep track of permutation 
			per[j] = per[row];
			per[row] = evin;
		}
		even = 1.0 / matrix[j][j];	// modify column 
		for(i = 0; i < nfree; i++) matrix[i][j] *= even;
		matrix[j][j] = even;
		for(k = 0; k < j; k++) {
			mjk = matrix[j][k];
			for(i = 0; i < j; i++) matrix[i][k] -= matrix[i][j] * mjk;
			for(i = j + 1; i < nfree; i++) matrix[i][k] -= matrix[i][j] * mjk;
			matrix[j][k] = -even * mjk;
		}
		for(k = j + 1; k < nfree; k++) {
			mjk = matrix[j][k];
			for(i = 0; i < j; i++) matrix[i][k] -= matrix[i][j] * mjk;
			for(i = j + 1; i < nfree; i++) matrix[i][k] -= matrix[i][j] * mjk;
			matrix[j][k] = -even * mjk;
		}
	}
	for(i = 0; i < nfree; i++) {		// finally, repermute the columns
		for (k = 0; k < nfree; k++) hv[per[k]] = matrix[i][k];
		for (k = 0; k < nfree; k++) matrix[i][k] = hv[k];
	}
	return TRUE;
}

//inimat sets up the matrix to be inverted in inivec. 
//inimat returns the reduced chi-squared.
extern double inimat( double s[5][5], double rl[5], double *x, double *y, int n, double *p, double *e, int ip[5], int nfree )
{
	double chi = 0.0;	// return value 
	double conv_fac = 0.0174532925; //conversion factor 
	double cosi, cosp, sini, sinp;
	double cost, sint, u;
	int i, j, k;
	
	//initialize the matrix and the vector
	for(j = 0; j < nfree; j++) {
		rl[j] = 0.0;
		for(k = 0; k <= j; k++) s[j][k] = 0.0;
	}
	cosp = cos( conv_fac * p[4] );	// cosine of p.a.
	sinp = sin( conv_fac * p[4] );	// sine of p.a. 
	cosi = cos( conv_fac * p[1] );	// cosine of inclination 
	sini = sin( conv_fac * p[1] );	// sine of inclination 
	for (i = 0; i < n; i++) {
		cost=( -( x[i] - p[2] ) * sinp + ( y[i] - p[3] ) * cosp ) / p[0]; //Calculate the rotated x and y coordinates
		sint=(- ( x[i] - p[2] ) * cosp - ( y[i] - p[3] ) * sinp ) / p[0] / cosi;
		u = 1.0 - cost * cost - sint * sint;	// difference with model
		//Now calculate the partial derivatives
		e[0] = -2.0 * ( cost * cost + sint * sint ) / p[0];
		e[1] =  2.0 * conv_fac * sint * sint * sini / cosi;
		e[2] =  2.0 * ( cost * sinp + sint * cosp / cosi ) / p[0];
		e[3] = -2.0 * ( cost * cosp - sint * sinp / cosi ) / p[0];
		e[4] = -2.0 * conv_fac * sint * cost * sini * sini / cosi;
		chi = chi + u * u;	// add to reduced chi-squared
		for (j = 0; j < nfree; j++) { //Now we fill the matrix and the vector
			rl[j] += u * e[ip[j]];
			for(k = 0; k <= j; k++) s[j][k] += e[ip[j]] * e[ip[k]];
		}
	}
	return chi;				// return chi-squared
}

// transposed matrix in complex
extern void zge_transpose(double complex *Transposed, double complex *M ,int n)
{
	int i,j;
	for(i=0;i<n;i++) for(j=0;j<n;j++) Transposed[i+n*j] = M[i*n+j];
}

//  MatrixComplexEigensystem: computes the eigenvectors and eigenValues of input matrix A
//  The eigenvectors are stored in columns
extern void matrix_complex_eigensystem(double complex *eigenvectorsVR, double complex *eigenvaluesW, double complex *A, int N)
{
	int i;
	double complex *AT = (double complex*) malloc(N*N*sizeof(double complex));

	zge_transpose(AT, A , N);
	char JOBVL ='N';   // Compute Righv  t eigenvectors
	char JOBVR ='V';   // Do not compute Left eigenvectors
	double complex VL[1];
	int LDVL = 1; 
	int LDVR = N;
	int LWORK = 4*N; 
	double complex *WORK =  (double complex*)malloc(LWORK*sizeof(double complex));
	double complex *RWORK = (double complex*)malloc(2*N*sizeof(double complex));
	int INFO;

	zgeev(&JOBVL, &JOBVR, &N, AT, &N, eigenvaluesW, VL, &LDVL, eigenvectorsVR, &LDVR, WORK, &LWORK, RWORK, &INFO);
	zge_transpose(AT, eigenvectorsVR, N);

	for(i=0;i<N*N;i++) eigenvectorsVR[i]=AT[i];
	free(WORK);
	free(RWORK);
	free(AT);
}

//compute and print determinant of complex matrix
extern double comp_print_det_matrix_1dim(char *matrixname, double complex *m_1dim, int n, char *txtout)
{
	int i,j;
	double det;
	FILE *wfile;
	
	if(n<5) {
		double m_det[4][4];
		for(i=0;i<n;i++) for(j=0;j<n;j++) m_det[i][j] = creal(m_1dim[i*n+j]);
		det = matrix_det_4d(m_det, n);
	}
	else if(n==5) {
		double m_det[5][5];
		for(i=0;i<n;i++) for(j=0;j<n;j++) m_det[i][j] = creal(m_1dim[i*n+j]);
		det = matrix_det_5d(m_det, n);
	}
	else {
		errorstop("dimension of matrix >5 ?\n");
		return 0;
	}
	printf("det %s = %le\n", matrixname, det);
	if(txtout != NULL) {
		wfile = fopen(txtout, "a");
		fprintf(wfile, "det %s = %le\n\n", matrixname, det);
		fclose(wfile);
	}
	return det;
}

//compute eigenvectors and eigenvalues of matrix m_1dim
extern int compute_eigenvector_parzen(double complex *eigenvectorsVR, double complex *eigenvaluesW, double complex *m_1dim)
{
	int i,j,k=0,n=4, doyouprintf = NO;
	double norm, abs_val;
	double complex conj_mu, conj_change[n];//, normf;
	
	matrix_complex_eigensystem(eigenvectorsVR, eigenvaluesW, m_1dim, n); //compute eigenvectors and eigenvalues
	//if(doyouprintf==YES) {
	//	printf("original eigenvalues:\n");
	//	for(i=0;i<n;i++) printf("& \\mu_%i = %lf + %lf i\\\\ \n", i, creal(eigenvaluesW[i]), cimag(eigenvaluesW[i]));
	//	printf_mat_1dim_complex("original eigenvectorsVR:", eigenvectorsVR,n, "VR1",1, NULL);
	//	comp_print_det_matrix_1dim("eigenvectorsVR", eigenvectorsVR, n, NULL);
	//	//printf("original eigenvalues:\n");
	//	//for(i=0;i<n;i++) printf("& \\mu_%i = %lf + %lf i\\\\ \n", i, creal(eigenvaluesW[i]), cimag(eigenvaluesW[i]));
	//}
	// pair the conjugated eigenvectors
	for(i=1;i<n;i++) {
		if(cabs(conj(eigenvaluesW[0])-eigenvaluesW[i])<1.e-9) {
			k=i;
			break;
		}
	}
	if(k==0) {
		printf("problem in eigenvalues, not conjugates?\n");
		printf("original eigenvalues:\n");
		for(i=0;i<n;i++) printf("& \\mu_%i = %lf + %lf i\\\\ \n", i, creal(eigenvaluesW[i]), cimag(eigenvaluesW[i]));
		return FALSE;
	}
	conj_mu = eigenvaluesW[1];
	eigenvaluesW[1] = eigenvaluesW[k];
	eigenvaluesW[k] = conj_mu;
	for(i=0;i<n;i++) conj_change[i] = eigenvectorsVR[1+n*i];
	for(i=0;i<n;i++) eigenvectorsVR[1+n*i] = eigenvectorsVR[k+n*i];
	for(i=0;i<n;i++) eigenvectorsVR[k+n*i] = conj_change[i];
	
	//for(j=0;j<n;j+=2) { //force conjugate
	//	eigenvaluesW[j+1] = conj(eigenvaluesW[j]);
	//	for(i=0;i<n;i++) eigenvectorsVR[j+1+n*i] = conj(eigenvectorsVR[j+n*i]);
	//}
	
	for(j=0;j<n;j+=2) { //put the imaginary part>0 first
		//test unity absolute value
		abs_val = creal(eigenvaluesW[j])*creal(eigenvaluesW[j])+cimag(eigenvaluesW[j])*cimag(eigenvaluesW[j]);
		if(fabs(abs_val-1)>1.e-1) {
			printf("eigenvalue is not unity value: %le\n", abs_val);
			return FALSE;
		}
		if(fabs(cimag(eigenvaluesW[j]))>TINYLENGTH && cimag(eigenvaluesW[j])<0) {
			conj_mu = eigenvaluesW[j];
			eigenvaluesW[j] = eigenvaluesW[j+1];
			eigenvaluesW[j+1] = conj_mu;
			for(i=0;i<n;i++) conj_change[i] = eigenvectorsVR[j+n*i];
			for(i=0;i<n;i++) eigenvectorsVR[j+n*i] = eigenvectorsVR[j+1+n*i];
			for(i=0;i<n;i++) eigenvectorsVR[j+1+n*i] = conj_change[i];
		}
	}
	
	//for(j=0;j<2;j++) { // swap the order of the eigenvector pairs
	//	conj_mu = eigenvaluesW[j];
	//	eigenvaluesW[j] = eigenvaluesW[j+2];
	//	eigenvaluesW[j+2] = conj_mu;
	//	for(i=0;i<n;i++) conj_change[i] = eigenvectorsVR[j+n*i];
	//	for(i=0;i<n;i++) eigenvectorsVR[j+n*i] = eigenvectorsVR[j+2+n*i];
	//	for(i=0;i<n;i++) eigenvectorsVR[j+2+n*i] = conj_change[i];
	//}
	/*if(doyouprintf==YES) {
		printf_mat_1dim_complex("adjusted 1 eigenvectorsVR:", eigenvectorsVR,n, "VR",1, NULL);
		comp_print_det_matrix_1dim("eigenvectorsVR",eigenvectorsVR, n, NULL);
	}//*/
	
	for(j=0;j<n;j+=2) {
		norm = 0;
		for(i=0;i<n;i+=2) {
			//printf("%lf*%lf-%lf*%lf\n",creal(eigenvectorsVR[j+n*i]),cimag(eigenvectorsVR[j+n*(i+1)]), cimag(eigenvectorsVR[j+n*i]),creal(eigenvectorsVR[j+n*(i+1)]));
			norm += creal(eigenvectorsVR[j+n*i])*cimag(eigenvectorsVR[j+n*(i+1)]) - cimag(eigenvectorsVR[j+n*i])*creal(eigenvectorsVR[j+n*(i+1)]);
		}
		//printf("ori norm: %le\n", norm);
		if(norm<0) { //if norm<0, swap the order of the concerned eigenvector pair
			//printf("norm<0, swap of eigenvectors\n");
			norm = -norm;
			conj_mu = eigenvaluesW[j];
			eigenvaluesW[j] = eigenvaluesW[j+1];
			eigenvaluesW[j+1] = conj_mu;
			for(i=0;i<n;i++) conj_change[i] = eigenvectorsVR[j+n*i];
			for(i=0;i<n;i++) eigenvectorsVR[j+n*i] = eigenvectorsVR[j+1+n*i];
			for(i=0;i<n;i++) eigenvectorsVR[j+1+n*i] = conj_change[i];
		}
		norm = sqrt(norm);
		//if(doyouprintf==YES) printf("norm = %le\n", norm);
		
		if(fabs(norm)<1.e-9) { //if norm too small, no normalisation done, (if the first pair is ok, the normalisation will be done for it)
			printf("normalisation value<1.e-9\n");
			printf("adjusted 1 eigenvalues:\n");
			for(i=0;i<n;i++) printf("& \\mu_%i = %lf + %lf i\\\\ \n", i, creal(eigenvaluesW[i]), cimag(eigenvaluesW[i]));
			return FALSE;
		}
		//normf = (1+I)/csqrt(2.*norm);
		
		for(i=0;i<n;i++) {
			eigenvectorsVR[j+n*i] /= norm;// *= normf;
			eigenvectorsVR[j+1+n*i] /= norm;//*= conj(normf);
		}
	}
	if(doyouprintf==YES) {
		printf("eigenvectors:\n");
		printf_mat_1dim_complex("normalised eigenvectorsVR:", eigenvectorsVR,n, "VR",1, NULL);
		//comp_print_det_matrix_1dim("eigenvectorsVR", eigenvectorsVR, n, NULL);
		printf("adjusted eigenvalues (imaginary>0 first):\n");
		for(i=0;i<n;i++) printf("& \\mu_%i = %lf + %lf i\\\\ \n", i, creal(eigenvaluesW[i]), cimag(eigenvaluesW[i]));
	}
	return TRUE;
}

// decouple real matrix m following Parzen method
extern int decouple_matrix_parzen(double m[4][4], double decoup_m[4][4], double decoup_to_coup[4][4], double output[6], char *txtout)
{
	int i, j, doyoudecouple, n=4, doyouprintf = YES;
	double beta, alpha, phase, tune;
	double complex ratio, exp_phi, exp_mphi, eigenvaluesW[n], m_1dim[n*n], eigenvectorsVR[n*n], par_t_evector_inv[n*n];
	double complex decoup_to_coup_1dim[n*n], decoup_to_coup_inv_1dim[n*n], decoup_m_1dim[n*n];
	double complex conj_mu, conj_change[n];
	//double complex par_t_evector[n*n];
	//double complex test[n*n];
	FILE *wfile;
	
	for(i=0;i<n;i++) {
		for(j=0;j<n;j++) {
			decoup_m[i][j] = 0.;
			//par_t_evector[i*n+j] = 0.;
			par_t_evector_inv[i*n+j] = 0.;
			m_1dim[i*n+j] = m[i][j];//normalized : times (-2I)^-(1/2)
		}
	}
	if(doyouprintf==YES) {
		printf_mat_1dim_complex("M=R^-1*T*R=", m_1dim, n, "m_1dim",3, txtout);
		comp_print_det_matrix_1dim("m", m_1dim, n, txtout);
	}
	doyoudecouple = compute_eigenvector_parzen(eigenvectorsVR, eigenvaluesW, m_1dim);
	
	if(doyoudecouple == FALSE) {
		printf("warning: degenerated eigenvalues, no decoupling done\n");
		for(i=0;i<n;i++) for(j=0;j<n;j++) decoup_m[i][j] = m[i][j];
		return FALSE;
	}
	if(txtout != NULL) wfile = fopen(txtout, "a");
	for(j=0;j<n;j+=2) {
		//printf("eig[%i,%i]=%lf+i*%lf, eig[%i,%i]=%lf+i*%lf\n", j+1,j, creal(eigenvectorsVR[j+n*(j+1)]), cimag(eigenvectorsVR[j+n*(j+1)]), j,j, creal(eigenvectorsVR[j+n*j]), cimag(eigenvectorsVR[j+n*j]));
		ratio = eigenvectorsVR[j+n*(j+1)]/eigenvectorsVR[j+n*j];
		beta = 1/cimag(ratio);
		printf("j=%i, ratio: (%lf,%lf)\n", j, creal(ratio), cimag(ratio));
		printf("j=%i, beta: %lf [m]\n", j, beta);
		if(beta<0) {
			printf("j=%i, beta<0, swapping of eigenvectors and eigenvalues\n",j);
			//printf_mat_1dim_complex("before swapping", eigenvectorsVR, n, "e",1, txtout);
			conj_mu = eigenvaluesW[j];
			eigenvaluesW[j] = eigenvaluesW[j+1];
			eigenvaluesW[j+1] = conj_mu;
			for(i=0;i<n;i++) conj_change[i] = eigenvectorsVR[j+n*i];
			for(i=0;i<n;i++) eigenvectorsVR[j+n*i] = eigenvectorsVR[j+1+n*i];
			for(i=0;i<n;i++) eigenvectorsVR[j+1+n*i] = conj_change[i];
			//printf_mat_1dim_complex("eigenvectorsVR", eigenvectorsVR, n, "eig",1, txtout);
			ratio = eigenvectorsVR[j+n*(j+1)]/eigenvectorsVR[j+n*j];
			beta = 1/cimag(ratio);
			if(beta<0) {
				printf("beta still <0, no decoupling done\n");
				if(txtout != NULL) fclose(wfile);
				return FALSE;
			}
		}
		phase = carg(eigenvectorsVR[j+n*j]);
		exp_phi = cexp(I*phase);
		exp_mphi = cexp(-I*phase);
		alpha = -beta*creal(ratio);
		tune = atan_ratio(cimag(eigenvaluesW[j]), creal(eigenvaluesW[j]))/(2.*PI);
		if(j==0) {
			output[0] = tune;
			output[2] = beta;
			output[3] = alpha;
		}
		else if(j==2) {
			output[1] = tune;
			output[4] = beta;
			output[5] = alpha;
		}
		if(doyouprintf==YES) {
			//printf("VR[%i]=%le +i %le\n", j+n*(j+1), creal(eigenvectorsVR[j+n*(j+1)]), cimag(eigenvectorsVR[j+n*(j+1)]));
			//printf("VR[%i]=%le +i %le\n", j+n*j, creal(eigenvectorsVR[j+n*j]), cimag(eigenvectorsVR[j+n*j]));
			printf("j=%i, ratio: (%lf,%lf)\n", j, creal(ratio), cimag(ratio));
			printf("j=%i, beta: %lf [m]\n", j, beta);
			printf("j=%i, phase: %lf [deg]\n", j, phase*180.0/PI);
			printf("j=%i, tune:%.15e, 1-tune:%.15e\n", j, tune, 1.-fabs(tune));
			printf("j=%i, exp_phi: (%lf, %lf)\n", j,creal(exp_phi), cimag(exp_phi));
			printf("j=%i, exp_mphi: (%lf, %lf)\n", j,creal(exp_mphi), cimag(exp_mphi));
			printf("j=%i, alpha: %lf\n", j, alpha);
		}
		//if(beta<0) {
		//	printf("beta<0, no decoupling done\n");
		//	if(txtout != NULL) fclose(wfile);
		//	return FALSE;
		//}
		if(txtout != NULL) {
			//fprintf(wfile, "j=%i, ratio: (%lf,%lf)\n", j, creal(ratio), cimag(ratio));
			fprintf(wfile, "j=%i, beta: %lf [m]\n", j, beta);
			fprintf(wfile, "j=%i, alpha: %lf\n", j, alpha);
			fprintf(wfile, "j=%i, phase: %lf [deg]\n", j, phase*180.0/PI);
			fprintf(wfile, "j=%i, tune:%lf\n", j, tune);
			//fprintf(wfile, "j=%i, exp_phi: (%lf, %lf)\n", j,creal(exp_phi), cimag(exp_phi));
			//fprintf(wfile, "j=%i, exp_mphi: (%lf, %lf)\n", j,creal(exp_mphi), cimag(exp_mphi));
		}
		par_t_evector_inv[j+n*j] = (-alpha-I)*exp_mphi/(sqrt(beta));
		par_t_evector_inv[j+n*(j+1)] = -(-alpha+I)*exp_phi/(sqrt(beta));
		par_t_evector_inv[j+1+n*j] = -sqrt(beta)*exp_mphi;
		par_t_evector_inv[j+1+n*(j+1)] = sqrt(beta)*exp_phi;
	}
	if(txtout != NULL) fclose(wfile);
	printf_mat_1dim_complex("par_t_evector_inv", par_t_evector_inv,n, "P-inv",1, txtout);
	
	
	matrix_dot_1dim(eigenvectorsVR, par_t_evector_inv, decoup_to_coup_1dim, n,n,n,n); //compute R matrix in Parzen paper 
	printf_mat_1dim_complex("fefore final norm", decoup_to_coup_1dim,n, "r_temp",1, txtout);
	for(i=0;i<n*n;i++) decoup_to_coup_1dim[i] /= -2.*I; //finalize normalisation!!
	if(doyouprintf==YES) {
		printf_mat_1dim_complex("R=", decoup_to_coup_1dim,n, "R",3, txtout);
		//printf_mat_1dim_complex("R = eigenvectorsVR * par_t_evector_inv", decoup_to_coup_1dim,n, "R",2, txtout);
		comp_print_det_matrix_1dim("R", decoup_to_coup_1dim, n, txtout);
	}
	matrix_complex_inverse(decoup_to_coup_1dim, decoup_to_coup_inv_1dim, n);
	if(doyouprintf==YES) {
		printf_mat_1dim_complex("R^-1=", decoup_to_coup_inv_1dim,n, "R^-1",3, txtout);
	}
	
	matrix_dot_1dim(decoup_to_coup_inv_1dim, m_1dim, decoup_to_coup_inv_1dim, n,n,n,n);
	matrix_dot_1dim(decoup_to_coup_inv_1dim, decoup_to_coup_1dim, decoup_m_1dim, n,n,n,n);
	if(doyouprintf==YES) {
	//	printf_mat_1dim_complex("T=", decoup_m_1dim,n, "T", 3,txtout);
		//printf_mat_1dim_complex("decoupled matrix = R_inv * m * R", decoup_m_1dim,n, "T", 2, txtout);
	//	comp_print_det_matrix_1dim("T", decoup_m_1dim, n, txtout);
	}
	for(i=0;i<n;i++) { //1dim matrix to 2dim vectors
		for(j=0;j<n;j++) {
			decoup_m[i][j] = creal(decoup_m_1dim[i*n+j]);
			decoup_to_coup[i][j] = creal(decoup_to_coup_1dim[i*n+j]);
		}
	}
	return TRUE;
}

// ************************************************************************************ //
//									Random number generator								//
// ************************************************************************************ //
//"Minimal" random number generator of Park and Miller with Bays_Durham shuffle and added safeguards.
//Returns a uniform random deviate between 0.0 and 1.0. (exclusive of endpoint values))
//Call with a negative integer to initialize; thereafter, do not alter idum between successive deviates in a sequence.
//RNMX should approximate the largest floating value that is less that 1.

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define RNMX (1.0-1.2e-7)  //RNMX should approximate the largest floating value that is less that 1. in simple precision the mantisse is store with 23 digits, 1.2e-7 ~ 2^-23

float ran1(long *idum)
{
	int j;
	long k;
	static long iy = 0;
	static long iv[NTAB];
	float temp;
	
	if(*idum <= 0 || !iy) {				//Initialize
		if (-(*idum) < 1) *idum = 1;	//Be sure to prevent idum = 0
		else *idum = -(*idum);
		for(j=NTAB+7;j>=0;j--) {		//Load the shuffle table (after 8 warm-ups)
			k=(*idum)/IQ;
			*idum = IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;						//Start here when not initializing
	*idum=IA*(*idum-k*IQ)-IR*k;			//Compute idum=(IA*idum)%IM without over-flows by Schrage's method
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;							//Will be in the range 0..NTAB-1
	iy=iv[j];							//Output previously stored value and refill the shuffle table
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX; //Because usere don't expect endpoint values
	else return temp;	
}

//gasdev retruns a normally (i.e. Gaussian) distributed with zero mean and unit variance, using ran1(idum) as the source of uniform derivates.
float gasdev(long *idum)
{
	static int iset = 0;
	static float gset;
	float fac,rsq,v1,v2;
	
	if(*idum < 0) iset = 0;		//Reinitialize
	if(iset == 0) {					//We don't have an extra deviate handy so
		do{
			v1=2.0*ran1(idum)-1.0;	//pick two uniform numbers in the square
			v2=2.0*ran1(idum)-1.0;	//extending from -1 to +1 in each direction
			rsq=v1*v1 + v2*v2;		//see if tey are in the unit circle,
		} while (rsq >= 1.0 || rsq == 0.0); //and if they are not, try again.
		fac=sqrt(-2.0*log(rsq)/rsq);		//Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time.
		gset=v1*fac;
		iset=1;						//Set flag.
		return v2*fac;
	} else {						//We have an extra deviate handy,
		iset=0;						//so unset the flag,
		return gset;				//and return it.
	}
}

