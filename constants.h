/*
 *  constant.h
 *  ringdesign
 *
 *  Copyright 2009 Kyoto University. All rights reserved.
 *
 */

#ifndef CONSTANTS
#define CONSTANTS

#include <stdio.h>
#include <math.h>


/* ******************** universal constants - SI Units ******************** */
#define PI				3.1415926535898
#define TWOPI			6.2831853071796
#define UNITCHARGE		1.602177330e-19				/*[C]*/
#define PROTONM			1.6726231e-27				/*[kg]*/
#define ELECTRONM		9.1093897e-31				/*[kg]*/
#define MUONM			1.8835323e-28				/*[kg]*/
#define CLIGHT			2.99792458e8				/*[m/s]*/
#define CLIGHT2			8.987551787368176e16		/*[m*m/s/s]*/
#define EULERGAMMA		0.57721566490153

/* ******************** simulation's parameters ******************** */
/* arrays size ------------------------------------------------------------------------------------------ */
#define KEYWORDL		25				/* maximum characters to define keyword string (see struct Cell in common.h) */
#define SIZE_MPARA		7				/* maximum number of parameters to charaterize the main body of a component (see struct Cell in common.h) */
#define SIZE_EFB		9				/* maximum number of parameters to charaterize one effective field boundary (EFB) of a component (see struct Cell in common.h) */
#define SIZE_ALIERR		11				/* maximum number of parameters to charaterize the alignment errors (see struct Cell in common.h) */

/* Input files parameter(s) ______________________________________________________________________________ */
#define MAX_CHARINLINE	300				/* maximum number of characters in each line of input files */
#define FILENAME_L		100				/* maximum number of character to define a name of file, such as a map file name */

/* Very small and very large numbers ______________________________________________________________________________ */
#define TINYLENGTH		1.e-12			/* [m] > 0!!! a very small length (used to check if a particle is right on a boundary for instance)*/
#define TINYDIMLESS		1.e-12			/* [\] */
#define EPS_MAPPOSI		1.e-5			/* [m] maximum allowed discrepancy between the expected node coordinates and its coodrinates in the field map data file */
#define NOOVFLOW_ENGE	100				/* [\] magic number used to avoid overflowing the largest double in field calculations with Enge type fringe field (see get_field.c)*/

/* find closed orbit ______________________________________________________________________________ */
#define FCLO_TRIES		10000			/* [\int] maximum number of try to find the closed orbit within the specified precision (see track.c)*/
#define FCLO_CROSS		1				/* [\int] number of times the Lattice is cross in each "try" (see track.c)*/

/* Tune calculation ______________________________________________________________________________ */
#define TRY_PERIOD_DISP	100					/* number of tries to get periodic dispersion in get_periodic_disp */
//------------------------------------------------------------------------------------------------------------------------------//
/* Runge-Kutta integration parameters ----------------------------------------------------------------------------------------- */
#define NVAR			7				/* [\int] number of unknowns in the differential equations to integrate*/
#define MAXSTP			(int)(5.0e7)	/* [\int] maximum number of steps to be used to cross one Cell*/
#define STPSZ_MAXFAC	15				/* [\int] if you use adaptative step size, step size cannot be STPSZ_MAXFAC times the specifyed step size (i.e. struct Cell の stepsize)*/


/* ******************** particles status ******************** */
#define ALIVE			0				/* Particle is not lost */
#define LOST_H_MIN		11				/* Particle has been lost on a horizontal collimator */
#define LOST_H_MAX		12				/* Particle has been lost on a horizontal collimator */
#define LOST_FARADAY	13				/* Particle has been lost on a horizontal collimator */
#define LOST_V_MIN		21				/* Particle has been lost on a vertical collimator */
#define LOST_V_MAX		22				/* Particle has been lost on a vertical collimator */
#define LOST_BACK		3				/* Particle is going backward */
#define LOST_L			4				/* Particle is more than 50% out of the cell longitudinally */
#define LOST_ERR_RF		990				/* Error has been done computing rf kick (probably Ekin is very close to zero, wich makes this error possible)*/
#define LOST_ERR1		991				/* Particle has been marked as lost because an error append in getready_part */
#define LOST_ERR2		992				/* Particle has been marked as lost because an error append in update_part */

/* ******************** Booleans! ******************** */
#define TRUE			0				/* [\int] better to keep 0 as it is related to particle status alive in acceptanceu and acceptancev*/
#define FALSE			1				/* [\int] has to be an integer!! */
#define YES				0				/* [\int] has to be an integer!! */
#define NO				1				/* [\int] has to be an integer!! */
#define TR				2				/* [\int] TOP RIGHT in title for easyplot */
#define TL				3				/* [\int] TOP LEFT in title for easyplot */
#define BR				4				/* [\int] BOTTOM RIGHT in title for easyplot */
#define BL				5				/* [\int] BOTTOM LEFT in title for easyplot */

/* ******************** Different types of instruments ******************** */
//Note: N0, defined above, means "no instrument in this cell". You have to choose in this section integers different than the one used for NO.
#define PICKUP			2				/* [\int] has to be an integer!! */
#define CUP				3				/* [\int] has to be an integer!! */


/* ******************** Colors in terminal ******************** */

// erase screen, and place the cursor on the left
//#define clrscr() printf("\033[H\033[2J")
#define CLRSCR() printf("\r")
//#define CLRSCR() printf("\033[2K\r")
//#define clrscr() printf("\033[2K\r")

// colors selection
#define COLOR(param) printf("\033[%sm",param)
//   param has to be const char *, empty (identical to "0") or made by values separated with ; 
// 0  initialize         1  high intensity (font)
// 5  blinking             7 inversed video
// 30, 31, 32, 33, 34, 35, 36, 37 letters color
// 40, 41, 42, 43, 44, 45, 46, 47 background color
// black, red, green, yellow, blue, magenta, cyan and white 

#endif

