/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nrutil.h.  Do not confuse this file with the same-named
   file nrutil.h that is supplied in the 'misc' subdirectory.
   *That* file is the one from the book, and contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */
   
//static float maxarg1,maxarg2;
//#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

//static float minarg1,minarg2;
//#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ? (minarg1) : (minarg2))

/* sign of a * b */
//#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#ifndef NRUTIL 
#define NRUTIL

void nrerror(char error_text[]);
double *vector(long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
void free_vector(double *v, long nl);
void free_matrix(double **m, long nrl, long ncl);


#endif
