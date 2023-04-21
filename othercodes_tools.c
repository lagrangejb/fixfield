/*
 *  othercodes_tools.c
 *  ringdesign
 *
 *  Copyright 2009 Kyoto University. All rights reserved.
 *
 */


#include "othercodes_tools.h"


extern void plot_mice_inspec(char *txtfilename, char *txt_mom1, char *txt_mom2, char *txt_mom3, char *outname, double mom_ref, double mom_spread)
{
	char buf[MAX_CHARINLINE], nameeps[200];
	int pid, nbpart;
	double T,X,Y,Z,r,Px,Py,Pz,Pt,P,E, pav;
	FILE *rfile=NULL;
	FILE *wfile1, *wfile2, *wfile3;
	
	pav=0.;
	nbpart=0;
	rfile = fopen(txtfilename,"r");
	if(rfile==NULL) errorstop("cannot open rfile");
	wfile1 = fopen(txt_mom1, "w");
	wfile2 = fopen(txt_mom2, "w");
	wfile3 = fopen(txt_mom3, "w");
	newline(rfile);
	if(fgets(buf, MAX_CHARINLINE-1, rfile) == NULL) errorstop("strange in plot_mice_inspec");
	if(sscanf(buf, "%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %i", &T,&X,&Y,&Z,&r,&Px,&Py,&Pz,&Pt,&P,&E,&pid) != 12) errorstop("strange2 in plot_mice_inspec");
	
	while(!feof(rfile)) {
	    if(fgets(buf, MAX_CHARINLINE-1, rfile) == NULL) break;
	    buf[MAX_CHARINLINE-1]='\0';

	    if(sscanf(buf, "%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %i", &T,&X,&Y,&Z,&r,&Px,&Py,&Pz,&Pt,&P,&E,&pid) != 12) {
	    	printf("why?\n");
			continue;
	    }
		pav+=P;
		nbpart++;
		if(P<=mom_ref-mom_spread) fprintf(wfile1, "%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %i\n", T,X,Y,Z,r,Px,Py,Pz,Pt,P,E,pid);
		else if(P>mom_ref-mom_spread && P<mom_ref+mom_spread) fprintf(wfile2, "%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %i\n", T,X,Y,Z,r,Px,Py,Pz,Pt,P,E,pid);
		else if(P>=mom_ref+mom_spread) fprintf(wfile3, "%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %i\n", T,X,Y,Z,r,Px,Py,Pz,Pt,P,E,pid);

	}
	fclose(rfile);
	fclose(wfile1);
	fclose(wfile2);
	fclose(wfile3);
	pav = pav/((double)nbpart);
	printf("pav=%lf\n",pav);
	
	//x-px
	sprintf(nameeps, "%s_x_px.eps",outname);
	easyplot3(txt_mom1, txt_mom2, txt_mom3, "2", "6", "2", "6", "2", "6", "points pt 7 ps 0.4 lc 3","points pt 7 ps 0.4 lc 0", "points pt 7 ps 0.4 lc 1", "x [mm]", "Px [MeV/c]", NULL, NULL, nameeps, NULL);
	sprintf(nameeps, "%s_x_pz.eps",outname);
	easyplot3(txt_mom1, txt_mom2, txt_mom3, "2", "8", "2", "8", "2", "8", "points pt 7 ps 0.4 lc 3","points pt 7 ps 0.4 lc 0", "points pt 7 ps 0.4 lc 1", "x [mm]", "Pz [MeV/c]", NULL, NULL, nameeps, NULL);
	sprintf(nameeps, "%s_x_px_mom1.eps",outname);
	easyplot(txt_mom1, "2", "6", "points pt 7 ps 0.4 lc 3", "x [mm]", "Px [MeV/c]", NULL, NULL, nameeps, NULL);
	sprintf(nameeps, "%s_x_px_mom2.eps",outname);
	easyplot(txt_mom2, "2", "6", "points pt 7 ps 0.4 lc 0", "x [mm]", "Px [MeV/c]", NULL, NULL, nameeps, NULL);
	sprintf(nameeps, "%s_x_px_mom3.eps",outname);
	easyplot(txt_mom3, "2", "6", "points pt 7 ps 0.4 lc 1", "x [mm]", "Px [MeV/c]", NULL, NULL, nameeps, NULL);
	
	/*sprintf(nameeps, "%s_x_px.eps",outname);
	easyplot3(txt_mom1, txt_mom2, txt_mom3, "2", "6", "2", "6", "2", "6", "points pt 7 ps 0.4 lc 1","points pt 7 ps 0.4 lc 2", "points pt 7 ps 0.4 lc 3", "x [mm]", "Px [MeV/c]", "[-250:250]", "[-50:50]", nameeps, NULL);
	sprintf(nameeps, "%s_x_px_mom1.eps",outname);
	easyplot(txt_mom1, "2", "6", "points pt 7 ps 0.4 lc 1", "x [mm]", "Px [MeV/c]", "[-250:250]", "[-50:50]", nameeps, NULL);
	sprintf(nameeps, "%s_x_px_mom2.eps",outname);
	easyplot(txt_mom2, "2", "6", "points pt 7 ps 0.4 lc 2", "x [mm]", "Px [MeV/c]", "[-250:250]", "[-50:50]", nameeps, NULL);
	sprintf(nameeps, "%s_x_px_mom3.eps",outname);
	easyplot(txt_mom3, "2", "6", "points pt 7 ps 0.4 lc 3", "x [mm]", "Px [MeV/c]", "[-250:250]", "[-50:50]", nameeps, NULL);
	//*/
	
	/*//y-py
	sprintf(nameeps, "%s_y_py.eps",outname);
	easyplot3(txt_mom1, txt_mom2, txt_mom3, "3", "7", "3", "7", "3", "7", "points pt 7 ps 0.4 lc 1","points pt 7 ps 0.4 lc 2", "points pt 7 ps 0.4 lc 3", "y [mm]", "Py [MeV/c]", NULL, NULL, nameeps, NULL);
	//r-pt
	sprintf(nameeps, "%s_r_pt.eps",outname);
	easyplot3(txt_mom1, txt_mom2, txt_mom3, "5", "9", "5", "9", "5", "9", "points pt 7 ps 0.4 lc 1","points pt 7 ps 0.4 lc 2", "points pt 7 ps 0.4 lc 3", "r [mm]", "Pt [MeV/c]", NULL, NULL, nameeps, NULL);
	//x-y
	sprintf(nameeps, "%s_x_y.eps",outname);
	easyplot3(txt_mom1, txt_mom2, txt_mom3, "2", "3", "2", "3", "2", "3", "points pt 7 ps 0.4 lc 1","points pt 7 ps 0.4 lc 2", "points pt 7 ps 0.4 lc 3", "x [mm]", "y [mm]", NULL, NULL, nameeps, NULL);
	//px-py
	sprintf(nameeps, "%s_px_py.eps",outname);
	easyplot3(txt_mom1, txt_mom2, txt_mom3, "6", "7", "6", "7", "6", "7", "points pt 7 ps 0.4 lc 1","points pt 7 ps 0.4 lc 2", "points pt 7 ps 0.4 lc 3", "Px [MeV/c]", "Py [MeV/c]", NULL, NULL, nameeps, NULL);
	//*/
}

extern void delta_beta_mice(char *txtfile_center, char *txtfile_off, char *outputfile)
{
	//char buf_c[MAX_CHARINLINE], buf_o[MAX_CHARINLINE];
	int i,io;
	double X, Y, Z, Px, Py, Pz, Pt, Bz, Energy, Emittance, Beta, Alpha, EmittanceLong, BetaLong, AlphaLong, Emittance6D, Weight;
	double Xo, Yo, Zo, Pxo, Pyo, Pzo, Pto, Bzo, Energyo, Emittanceo, Betao, Alphao, EmittanceLongo, BetaLongo, AlphaLongo, Emittance6Do, Weighto;
	double deltabeta;
	FILE *wfile;
	FILE *rfile_c=NULL;
	FILE *rfile_o=NULL;
	
	wfile = fopen(outputfile, "w");
	rfile_c = fopen(txtfile_center, "r");
	if(rfile_c==NULL) errorstop("cannot open rfile_c");
	rfile_o = fopen(txtfile_off, "r");
	if(rfile_o==NULL) errorstop("cannot open rfile_o");
	newline(rfile_c);
	newline(rfile_o);
	
	while(!feof(rfile_c)) {
	   /* if(fgets(buf_c, MAX_CHARINLINE-1, rfile_c) == NULL) break;
	    buf_c[MAX_CHARINLINE-1]='\0';
	    if(fgets(buf_o, MAX_CHARINLINE-1, rfile_o) == NULL) break;
	    buf_o[MAX_CHARINLINE-1]='\0';
	    if(sscanf(buf_c, "%i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &i, &X, &Y, &Z, &Px, &Py, &Pz, &Pt, &Bz, &Energy, &Emittance, &Beta, &Alpha, &EmittanceLong, &BetaLong, &AlphaLong, &Emittance6D, &Weight) != 18) {
	    	printf("why_c?\n");
			continue;
	    }
	    if(sscanf(buf_o, "%i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &io, &Xo, &Yo, &Zo, &Pxo, &Pyo, &Pzo, &Pto, &Bzo, &Energyo, &Emittanceo, &Betao, &Alphao, &EmittanceLongo, &BetaLongo, &AlphaLongo, &Emittance6Do, &Weighto) != 18) {
	    	printf("why_o?\n");
			continue;
	    }//*/
		fscanf(rfile_c, "%i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &i, &X, &Y, &Z, &Px, &Py, &Pz, &Pt, &Bz, &Energy, &Emittance, &Beta, &Alpha, &EmittanceLong, &BetaLong, &AlphaLong, &Emittance6D, &Weight);
		printf("z=%lf\n",Z);
		fscanf(rfile_o, "%i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &io, &Xo, &Yo, &Zo, &Pxo, &Pyo, &Pzo, &Pto, &Bzo, &Energyo, &Emittanceo, &Betao, &Alphao, &EmittanceLongo, &BetaLongo, &AlphaLongo, &Emittance6Do, &Weighto);
		if(fabs(Z-Zo)>TINYLENGTH) errorstop("z don't match");
		else deltabeta = Betao-Beta;
		fprintf(wfile, "%lf	%le\n", Z, deltabeta);
	}
	
	fclose(wfile);
	fclose(rfile_c);
	fclose(rfile_o);
	
}

extern void plot_magnetic_field_mice(char *name1, char *name2)
{
	int i;
	double zpos1[161], zpos2[161], bx1[161], bx2[161], by1[161], by2[161], bz1[161], bz2[161], difx[161], dify[161], difz[161];
	FILE *rfile1 = NULL;
	FILE *rfile2 = NULL;
	FILE *wfile;
	rfile1 = fopen(name1, "r");
	rfile2 = fopen(name2, "r");
	if(rfile1 == NULL) errorstop("cannot open 1");
	if(rfile2 == NULL) errorstop("cannot open 2");
	wfile = fopen("difference.dat", "w");
	for(i=0;i<161;i++) {
		fscanf(rfile1, "%lf %lf %lf %lf\n", &(zpos1[i]), &(bx1[i]), &(by1[i]), &(bz1[i]));
		fscanf(rfile2, "%lf %lf %lf %lf\n", &(zpos2[i]), &(bx2[i]), &(by2[i]), &(bz2[i]));
		if(fabs(zpos1[i]-zpos2[i]) < 1.e-9) {
			difx[i] = bx1[i]-bx2[i];
			dify[i] = by1[i]-by2[i];
			difz[i] = bz1[i]-bz2[i];
			printf("i = %i\n", i);
			fprintf(wfile,"%le \t %le \t %le \t %le \t %le\n", zpos1[i], difx[i], dify[i], difz[i], difz[i]/bz1[i]);
		}
		else {
			printf("zpos1 = %le, zpos2 = %le\n", zpos1[i], zpos2[i]);
			errorstop("different z!!");
		}
		
	}
	
	fclose(rfile1);
	fclose(rfile2);
	fclose(wfile);
	easyplot("difference.dat", "($1*1.e-3)", "($4*1e3)", "lines lt 1 lw 4", "z [m]", "delta B_z [T]", NULL, NULL, "DifferenceMagneticField.eps", "grid");
	easyplot("difference.dat", "($1*1.e-3)", "($5)", "lines lt 1 lw 4", "z [m]", "delta B_z/B_z []", NULL, NULL, "prop_dif.eps", "grid");
	easyplot("difference.dat", "($1*1.e-3)", "($5)", "lines lt 1 lw 4", "z [m]", "delta B_z/B_z []", NULL, NULL, "prop_dif2.eps", "grid");
	easyplot2(name1, name2, "($1*1.e-3)", "($4*1000)", "($1*1.e-3)", "($4*1000)", "lines lt 1 lw 4 lc 3", "lines lt 1 lw 4 lc 1", NULL, NULL, "z [m]", "B_z [T]", NULL, NULL, "MagneticField.eps", "grid");
	
}

extern void plot_magnetic_field_sam(char *name1, char *name2)
{
	char shell_write[300], name[300];
	int i, n, nblines;
	FILE *shell_read;
	sprintf(shell_write,"wc -l %s", name1);
	shell_read = popen(shell_write,"r");
	fscanf(shell_read, "%i %s", &nblines, name);
	pclose(shell_read);
	double xpos1[nblines], xpos2[nblines], ypos1[nblines], ypos2[nblines],zpos1[nblines], zpos2[nblines], bx1[nblines], bx2[nblines], by1[nblines], by2[nblines], bz1[nblines], bz2[nblines], difx[nblines], dify[nblines], difz[nblines];
	FILE *rfile1 = NULL;
	FILE *rfile2 = NULL;
	FILE *wfile;
	FILE *wfile2;
	rfile1 = fopen(name1, "r");
	rfile2 = fopen(name2, "r");
	if(rfile1 == NULL) errorstop("cannot open 1");
	newline(rfile1);
	if(rfile2 == NULL) errorstop("cannot open 2");
	wfile = fopen("difference.dat", "w");
	wfile2 = fopen("sam_dif.dat", "w");
	
	for(i=0;i<nblines;i++) {
		fscanf(rfile1, "%lf  %lf  %lf  %lf  %lf  %lf\n", &(xpos1[i]), &(ypos1[i]), &(zpos1[i]), &(bx1[i]), &(by1[i]), &(bz1[i]));
		xpos1[i] *= 1.e-2;
		ypos1[i] *= 1.e-2;
		zpos1[i] *= 1.e-2;
		fscanf(rfile2, "%lf  %lf  %lf  %lf  %lf  %lf\n", &(xpos2[i]), &(ypos2[i]), &(zpos2[i]), &(bx2[i]), &(by2[i]), &(bz2[i]));
		xpos2[i] -=36.2;
		n=0;
		while(fabs(zpos1[i]-zpos2[i]) > 1.e-3 || fabs(ypos1[i]-ypos2[i]) > 1.e-3 || fabs(xpos1[i]-xpos2[i]) > 1.e-3) {
			fscanf(rfile1, "%lf  %lf  %lf  %lf  %lf  %lf\n", &(xpos1[i]), &(ypos1[i]), &(zpos1[i]), &(bx1[i]), &(by1[i]), &(bz1[i]));
			xpos1[i] *= 1.e-2;
			ypos1[i] *= 1.e-2;
			zpos1[i] *= 1.e-2;
			printf("i=%i, xdif = %le, ydif = %le, zdif = %le\n", i, fabs(xpos1[i]-xpos2[i]), fabs(ypos1[i]-ypos2[i]), fabs(zpos1[i]-zpos2[i]));
			n++;
			if(n>20) errorstop("problem in while");
		}	
		difx[i] = bx1[i]-bx2[i];
		dify[i] = by1[i]-by2[i];
		difz[i] = bz1[i]-bz2[i];
		fprintf(wfile,"%le \t %le \t %le \t %le \t %le \t %le\n", xpos2[i], ypos2[i], zpos2[i], difx[i], dify[i], difz[i]);
		if(fabs(ypos2[i]-2.5) < 1.e-4) fprintf(wfile2, "%le  %le\n", xpos2[i], difz[i]);
	}
	fclose(rfile1);
	fclose(rfile2);
	fclose(wfile);
	fclose(wfile2);
	easyplot("difference.dat", "($2)", "($6)", "lines lt 1 lw 4", "y [m]", "delta B_z [T]", NULL, NULL, "DifferenceMagneticField.eps", "grid");
	easyplot("sam_dif.dat", "($1)", "($2)", "lines lt 1 lw 4", "y [m]", "delta B_z [T]", NULL, NULL, "sam_dif.eps", "grid");
	//easyplot("difference.dat", "($1*1.e-3)", "($5)", "lines lt 1 lw 4", "z [m]", "delta B_z/B_z []", "[-4:-0.01]", NULL, "prop_dif.eps", "grid");
	//easyplot("difference.dat", "($1*1.e-3)", "($5)", "lines lt 1 lw 4", "z [m]", "delta B_z/B_z []", "[0.01:4]", NULL, "prop_dif2.eps", "grid");
	//easyplot2(name1, name2, "($1*1.e-3)", "($4*1000)", "($1*1.e-3)", "($4*1000)", "lines lt 1 lw 4 lc 3", "lines lt 1 lw 4 lc 1", NULL, NULL, "z [m]", "B_z [T]", NULL, NULL, "MagneticField.eps", "grid");
}

extern void change_input_zgoubi(char *txtfilename, char *output, int nblines)
{
	int n, n_str=0, n_half=0, flag_str=NO;
	char label[20];
	double s[50987], muy[50987], betay[50987], alphay[50987], gammay[50987], dispy[50987], disppy[50987], muz[50987], betaz[50987], alphaz[50987], gammaz[50987], dispz[50987], disppz[50987];
	FILE *rfile;
	FILE *wfile;
	
	rfile = fopen(txtfilename,"r");
	wfile = fopen(output, "w");
	for(n=0;n<nblines-1;n++) {
		newline(rfile);
		//printf("n=%i\n", n);
		//fscanf(rfile, "%lf", &s[n]);
		fscanf(rfile, "%lf	%s	%lf	%lf	%lf	%lf %lf	%lf %lf	%lf	%lf	%lf %lf	%lf", &s[n], label, &muy[n], &betay[n], &alphay[n], &gammay[n], &dispy[n], &disppy[n], &muz[n], &betaz[n], &alphaz[n], &gammaz[n], &dispz[n], &disppz[n]);
		//printf("%lf\n", s[n]);
		//if(strcmp(label, "arc_m") == 0) printf("arc_m\n");
		//else printf("no");
		if((strcmp(label, "s_qd1") == 0) && flag_str == NO) {
			n_str=n;
			printf("n_str=%i\n",n_str);
			flag_str = YES;
		}
		//if(flag_str == YES && (strcmp(label, "arc_m") == 0)) {
		//	n_half = n;
		//	printf("nhalf=%i\n",n_half);
		//	break;
		//}
	}
	n_half = floor((nblines+n_str)/2);
	
	fclose(rfile);
	for(n=n_half;n<nblines-1;n++) fprintf(wfile, "%lf	%lf	%lf	%lf	%lf %lf	%lf %lf	%lf	%lf	%lf %lf	%lf\n", s[n]-s[n_half], muy[n], betay[n], alphay[n], gammay[n], dispy[n], disppy[n], muz[n], betaz[n], alphaz[n], gammaz[n], dispz[n], disppz[n]);
	for(n=0;n<n_half;n++) fprintf(wfile, "%lf	%lf	%lf	%lf	%lf %lf	%lf %lf	%lf	%lf	%lf %lf	%lf\n", s[n]+s[nblines-2]-s[n_half], muy[n], betay[n], alphay[n], gammay[n], dispy[n], disppy[n], muz[n], betaz[n], alphaz[n], gammaz[n], dispz[n], disppz[n]);

	
	fclose(wfile);  
}

extern void change_input_zgoubi_bfield(char *txtfilename, char *output)
{
	int n, n_str=0;
	//int flag_str=NO;
	double s[50987], b[50987], btemp;
	FILE *rfile=NULL;
	FILE *wfile;
	
	rfile = fopen(txtfilename,"r");
	if(rfile==NULL) errorstop("cannot open rfile");
	wfile = fopen(output, "w");
	newline(rfile);
	
	for(n=0;n<50987-1;n++) {
		newline(rfile);
		//printf("n=%i\n", n);
		//fscanf(rfile, "%lf", &s[n]);
		fscanf(rfile, "%le	%le", &s[n], &btemp);
		if(n<14950) {
			b[n]=-btemp;
		}
		else b[n]=btemp;
		//printf("%lf\n", s[n]);
		//if(strcmp(label, "arc_m") == 0) printf("arc_m\n");
		//else printf("no");
		//if(s[n]>25469.9451103/2. && flag_str == NO) {
		//	n_str=n;
		//	printf("n_str=%i\n",n_str);
		//	flag_str = YES;
		//}
		//if(flag_str == YES && (strcmp(label, "arc_m") == 0)) {
		//	n_half = n;
		//	printf("nhalf=%i\n",n_half);
		//	break;
		//}
	}
	n_str=32968;
	
	fclose(rfile);
	for(n=n_str;n<50987-1;n++) fprintf(wfile, "%le	%le\n", s[n]-s[n_str], b[n]);
	for(n=0;n<n_str;n++) fprintf(wfile, "%le	%le\n", s[n]+s[50985]-s[n_str], b[n]);

	
	fclose(wfile);  
}

extern void change_input_zgoubi_acceptance(char *txtfilename, char *output)
{
	int n;
	double x[202], xp[202];
	FILE *rfile=NULL;
	FILE *wfile;
	
	rfile = fopen(txtfilename,"r");
	if(rfile==NULL) errorstop("cannot open rfile");
	wfile = fopen(output, "w");
	
	for(n=0;n<201;n++) fscanf(rfile, "%lf", &x[n]);
	newline(rfile);
	for(n=0;n<201;n++) fscanf(rfile, "%lf", &xp[n]);
	
	
	fclose(rfile);
	for(n=0;n<201;n++) fprintf(wfile, "%lf	%lf\n", x[n], xp[n]);
	
	fclose(wfile);  
}

extern void map_opal_circ(struct Cell *cell, char *mapfilename, double rstart, double rend, int nbr, double thstart, double thend, int nbth, double zstart, double zend, int nbz, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int ir,ith,iz;
	double rstep, thstep, zstep, r, th, z, x, y, br, bth, bx, by, bz; 
	double lscale = 100.; //cm
	double bscale = 10.; //kG
	double thscale = 180./PI; //deg
	FILE *wfile;
	
	wfile = fopen(mapfilename, "w");
	fprintf(wfile, "cylindrical mesh, cylindrical coordinates\n");
	fprintf(wfile, "units - coordinates cm and deg., field kG\n");
	fprintf(wfile, "loop order (outer to inner): r	z(vertical)	theta\n");
	fprintf(wfile, "list order: r	z(vertical)	theta	Br	Bz	Bth\n");
	fprintf(wfile, "r from %le to %le, %i nodes (%i steps)\n", rstart*lscale, rend*lscale, nbr,nbr-1);
	fprintf(wfile, "z from %le to %le, %i nodes (%i steps)\n", zstart*lscale, zend*lscale, nbz,nbz-1);
	fprintf(wfile, "theta from %le to %le, %i nodes (%i steps)\n\n", thstart*thscale, thend*thscale, nbth,nbth-1);
	
	if(nbz>1) zstep = (zend-zstart)/(nbz-1);
	else zstep = 0.;
	if(nbr>1) rstep = (rend-rstart)/(nbr-1);
	else rstep = 0.;
	if(nbth>1) thstep = (thend-thstart)/(nbth-1);
	else thstep = 0.;
	
	for(ir=0;ir<nbr;ir++) {
		r = rstart + ir*rstep;
		for(iz=0;iz<nbz;iz++) {
			z = zstart + iz*zstep;
			for(ith=0;ith<nbth;ith++) {
				th = thstart + ith*thstep;
				x = r*cos(th);
				y = r*sin(th);
				if(get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp)!=ALIVE) errorstop("bfield not alive in map_opal_circ");
				br = bx*cos(th) + by*sin(th);
				bth = by*cos(th) - bx*sin(th);
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", r*lscale,z*lscale,th*thscale,br*bscale,bz*bscale,bth*bscale);
			}
		}
	}
	
	fclose(wfile);
}

extern void map_opal_cart(struct Cell *cell, char *mapfilename, double xstart, double xend, int nbx, double ystart, double yend, int nby, double zstart, double zend, int nbz, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	int ix,iy,iz;
	double xstep, ystep, zstep, x, y, z, bx, by, bz; 
	double lscale = 100.; //cm
	double bscale = 10.; //kG
	FILE *wfile;
	
	wfile = fopen(mapfilename, "w");
	fprintf(wfile, "cartesian mesh, cartesian coordinates\n");
	fprintf(wfile, "units - coordinates cm, field kG\n");
	fprintf(wfile, "loop order (outer to inner): x	z(vertical)	y (longitudinal)\n");
	fprintf(wfile, "list order: x	z(vertical)	y (longitudinal)	Bx	Bz	By\n");
	fprintf(wfile, "x from %le to %le, %i nodes (%i steps)\n", xstart*lscale, xend*lscale, nbx, nbx-1);
	fprintf(wfile, "z from %le to %le, %i nodes (%i steps)\n", zstart*lscale, zend*lscale, nbz, nbz-1);
	fprintf(wfile, "y from %le to %le, %i nodes (%i steps)\n\n", ystart*lscale, yend*lscale, nby, nby-1);
	
	zstep = comp_step(zstart, zend, nbz);
	xstep = comp_step(xstart, xend, nbx);
	ystep = comp_step(ystart, yend, nby);
	
	for(ix=0;ix<nbx;ix++) {
		x = xstart + ix*xstep;
		for(iz=0;iz<nbz;iz++) {
			z = zstart + iz*zstep;
			for(iy=0;iy<nby;iy++) {
				y = ystart + iy*ystep;
				if(get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp)!=ALIVE) errorstop("bfield not alive in map_opal_cart");
				fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", x*lscale,z*lscale,y*lscale,bx*bscale,bz*bscale,by*bscale);
			}
		}
	}
	fclose(wfile);
}

////////////////////////////////??///////////

extern void adjust_b0_opal(double x, double y, double z, double bz_opal, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), double precision, int doyouf)
{
	int i,j, n=500;
	double sign,bx,by,bz;
	
	if(doyouf==YES) {
		sign = -1.;
		printf("adjust B0F\n");
	}
	else {
		sign = 1.;
		printf("adjust B0D\n");
	}
	
	for(i=0;i<n;i++) {
		if(get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp)!=ALIVE) errorstop("bfield not alive");
		printf("i=%i, bz = %le, bz_opal = %le, b01=%le, b02=%le\n", i, bz, bz_opal, cell->mpara[0][2], cell->mpara[cell->nbcomp-1][2]);
		if(fabs(bz-bz_opal)<precision) break;
		for(j=0;j<cell->nbcomp;j++) {
			if(sign*cell->mpara[j][2] > 0) cell->mpara[j][2] -= (bz-bz_opal)/2.;
		}
	}
	for(j=0;j<cell->nbcomp;j++) {
		if(sign*cell->mpara[j][2] > 0) {
			printf("B0 = %le\n", cell->mpara[j][2]);
			break;
		}
	}
}

extern void opal_magnet_position_angle(char *opal_file, char *opal_gnuplot)
{
	int i, nblines;
	double x, y, x_n, y_n;
	FILE *rfile=NULL;
	FILE *wfile;
	
	
	nblines = get_nb_lines_file(opal_file);
	rfile = fopen(opal_file,"r");
	if(rfile==NULL) errorstop("cannot open rfile");
	wfile = fopen(opal_gnuplot, "w");
	
	for(i=0;i<nblines+1;i++) {
		fscanf(rfile,"%le	%le	%le	%le", &x, &y, &x_n, &y_n);
		fprintf(wfile, "%le	%le\n", x, y);
		//fprintf(wfile, "%le	%le\n", x+x_n, y+y_n);
		fprintf(wfile, "%le	%le\n", x-y_n, y+x_n);
		fprintf(wfile, "\n");
	}
	fclose(rfile);
	fclose(wfile);
}

extern void compare_field_opal_vffa(char *file_dif_prefix, double z_1, struct Cell *cell1, void(*add_contribution_comp1)(double,double,double,double*,double*,double*,struct Cell*,int), 
double z_2, struct Cell *cell2, void(*add_contribution_comp2)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	char namex[300], namey[300], namez[300]; 
	int i,j, iz, ic;
	double x,y, newx, x_0, y_0;
	FILE *wfilex,*wfiley,*wfilez;
	sprintf(namex,  "%s_bx.dat", file_dif_prefix);
	sprintf(namey,  "%s_by.dat", file_dif_prefix);
	sprintf(namez,  "%s_bz.dat", file_dif_prefix);
	emptyfile(namex);
	emptyfile(namey);
	emptyfile(namez);
	
	iz = (int)((z_2 - cell2->map.mapdim[4])/cell2->map.stepsize[2]);
	printf("iz = %i, z2 = %le, zmap = %le\n", iz, z_2, cell2->map.node[0][0][iz].coord[2]);
	for(i=0;i<cell2->map.nnodes[0];i++) {
		for(j=0;j<cell2->map.nnodes[1];j++) {
			x = cell2->map.node[i][j][iz].coord[0]*cos(cell2->map.node[i][j][iz].coord[1]);
			y = cell2->map.node[i][j][iz].coord[0]*sin(cell2->map.node[i][j][iz].coord[1]);
			if(cell2->map.node[i][j][iz].coord[1] < 6*PI/180.) ic = 0;
			else if(cell2->map.node[i][j][iz].coord[1] < 18*PI/180.) ic = 1;
			else ic = 2;
			x_0 = cell1->mpara[ic][1]*cos(cell1->mpara[ic][0]);
			y_0 = cell1->mpara[ic][1]*sin(cell1->mpara[ic][0]);
			newx = (y-y_0)*sin(cell1->mpara[ic][0]+cell1->mpara[ic][5]) + (x-x_0)*cos(cell1->mpara[ic][0]+cell1->mpara[ic][5]);
			if(fabs(newx)<cell1->mpara[ic][6]) compare_field_point(file_dif_prefix, x, y, z_1, cell1, add_contribution_comp1, x, y, z_2, cell2, add_contribution_comp2,NO);
			else {
				wfilex = fopen(namex,"a");
				fprintf(wfilex,"%lf	%lf	0\n", y, x);
				fclose(wfilex);
				wfiley = fopen(namey,"a");
				fprintf(wfiley,"%lf	%lf	0\n", y, x);
				fclose(wfiley);
				wfilez = fopen(namez,"a");
				fprintf(wfilez,"%lf	%lf	0\n", y, x);
				fclose(wfilez);
			}
		}
		wfilex = fopen(namex,"a");
		fprintf(wfilex,"\n");
		fclose(wfilex);
		wfiley = fopen(namey,"a");
		fprintf(wfiley,"\n");
		fclose(wfiley);
		wfilez = fopen(namez,"a");
		fprintf(wfilez,"\n");
		fclose(wfilez);
		printf("%i/%i done\n", i, cell2->map.nnodes[0]);
	}
	
}

extern void compare_field_map_opal(char *file_opal, char *fileout, struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int))
{
	char namex[300],namey[300],namez[300],outx[300], outy[300], outz[300];
	int i,j,nbx,nby;
	double x,y,z, bx_opal,by_opal,bz_opal, bx,by,bz, difbx,difby,difbz, r,th, xstep,ystep,xmax,ymax;
	double difminx=+10.,difmaxx=-10.,difminy=+10.,difmaxy=-10.,difminz=+10.,difmaxz=-10.;
	FILE *rfile=NULL;
	FILE *wfilex,*wfiley,*wfilez;
	
	xstep = 0.005;
	nbx = (int) 6/xstep +1;
	ystep = 0.005;
	nby = (int) 6/ystep +1;
	
	printf("nbx = %i, nby = %i\n",nbx,nby);
	
	rfile = fopen(file_opal, "r");
	if(rfile==NULL) errorstop("cannot open opal file\n");
	for(i=0;i<8;i++) newline(rfile);
	sprintf(namex,"%s_x.dat", fileout);
	wfilex = fopen(namex,"w");
	sprintf(namey,"%s_y.dat", fileout);
	wfiley = fopen(namey,"w");
	sprintf(namez,"%s_z.dat", fileout);
	wfilez = fopen(namez,"w");
	for(i=0;i<nbx;i++) {
		for(j=0;j<nby;j++) {
			fscanf(rfile, "%le	%le	%le	%le	%le	%le",&x,&y,&z,&bx_opal,&by_opal,&bz_opal);
			//printf("%le	%le	%le	%le	%le	%le\n",x,y,z,bx_opal,by_opal,bz_opal);
			if(fabs(x-i*xstep)>EPS_MAPPOSI || fabs(y-j*ystep)>EPS_MAPPOSI) errorstop("in compare_field_map_opal, dimensions not ok\n");
			if(j==0 || j==nby-1) ;
			else {
				r = sqrt(x*x + y*y);
				th = atan_ratio(y, x);
				if(r>cell->collim.rmin && r<cell->collim.rmax && th>0 && th<cell->boun.thmax) {
					if(get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp)!=ALIVE) errorstop("bfield not alive");
					if(fabs(bx)>TINYDIMLESS) difbx = (bx_opal/10. - bx);///bx;
					else difbx=0;
					if(fabs(by)>TINYDIMLESS) difby = (by_opal/10. - by);///by;
					else difby=0;
					if(fabs(bz)>TINYDIMLESS) difbz = (bz_opal/10. - bz);///bz;
					else difbz=0;
					difminx = MIN(difbx,difminx);
					difmaxx = MAX(difbx,difmaxx);
					difminy = MIN(difby,difminy);
					difmaxy = MAX(difby,difmaxy);
					//difminz = MIN(difbz,difminz);
					//difmaxz = MAX(difbz,difmaxz);
					if(difminz>difbx) difminz = difbz;
					if(difmaxz<difbz) {
						difmaxz = difbz;
						xmax = x;
						ymax= y;
					}//*/
					fprintf(wfilex, "%le	%le	%le\n", x,y,difbx);
					fprintf(wfiley, "%le	%le	%le\n", x,y,difby);
					fprintf(wfilez, "%le	%le	%le\n", x,y,difbz);
				}
			}
		}
		fprintf(wfilex,"\n");
		fprintf(wfiley,"\n");
		fprintf(wfilez,"\n");
	}
	printf("difx(min,max) = (%le,%le)\n",difminx,difmaxx);
	printf("dify(min,max) = (%le,%le)\n",difminy,difmaxy);
	printf("difz(min,max) = (%le,%le)\n",difminz,difmaxz);
	printf("(x,y)max = (%le,%le)\n",xmax,ymax);
	fclose(rfile);
	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
	sprintf(outx,"%s_x.eps",fileout);
	sprintf(outy,"%s_y.eps",fileout);
	sprintf(outz,"%s_z.eps",fileout);
	easyplot3dmap(namex, "x [m]", "y [m]", "dif", NULL, "[-0.1:3.5]", outx, "size ratio -1",difminx,difmaxx);
	easyplot3dmap(namey, "x [m]", "y [m]", "dif", NULL, "[-0.1:3.5]", outy, "size ratio -1",difminy,difmaxy);
	easyplot3dmap(namez, "x [m]", "y [m]", "dif", NULL, "[-0.1:3.5]", outz, "size ratio -1",difminz,difmaxz);
}

extern void compare_field_map_opal2(char *file_opal, char *fileout, struct Cell *cell, double z, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), 
double xmin, double xmax, double xstep, double ymin, double ymax, double ystep, double zmin, double zmax, double zstep)
{
	char namex[300],namey[300],namez[300],outx[300], outy[300], outz[300];
	int i,j,k,nbx,nby, nbz, nblines, skip_lines=8, check_flag;
	double x,y,z_opal, bx_opal,by_opal,bz_opal, bx,by,bz, difbx,difby,difbz, r,th;
	double difminx=+1.e9,difmaxx=-1.e9,difminy=+1.e9,difmaxy=-1.e9,difminz=+1.e9,difmaxz=-1.e9;
	FILE *rfile=NULL;
	FILE *wfilex,*wfiley,*wfilez;
	
	nblines = get_nb_lines_file(file_opal);
	nblines -= skip_lines;
	if(fabs(xstep)>TINYLENGTH) nbx = (int) ((xmax - xmin)/xstep) +1;
	else nbx = 1;
	if(fabs(ystep)>TINYLENGTH) nby = (int) ((ymax - ymin)/ystep) +1;
	else nby = 1;
	if(fabs(zstep)>TINYLENGTH) nbz = (int) ((zmax - zmin)/zstep) +1;
	else nbz = 1;
	if(nblines-nbx*nby*nbz !=0) errorstop("strange nb of lines");
	
	rfile = fopen(file_opal, "r");
	if(rfile==NULL) errorstop("cannot open opal file\n");
	for(i=0;i<skip_lines;i++) newline(rfile);
	sprintf(namex,"%s_x.dat", fileout);
	wfilex = fopen(namex,"w");
	sprintf(namey,"%s_y.dat", fileout);
	wfiley = fopen(namey,"w");
	sprintf(namez,"%s_z.dat", fileout);
	wfilez = fopen(namez,"w");
	
	for(i=0;i<nbx;i++) {
		for(j=0;j<nby;j++) {
			for(k=0;k<nbz;k++) {
				fscanf(rfile, "%le	%le	%le	%le	%le	%le",&x,&y,&z_opal,&bx_opal,&by_opal,&bz_opal);
				if(cell->boun.thmax != 0) {
					r = sqrt(x*x+y*y);
					if(x!=0) th = atan2(y,x);
					else th = sign(y)*PI/2;
					if(r>cell->collim.rmin && r<cell->collim.rmax && th>0 && th<cell->boun.thmax && fabs(z-z_opal)<TINYLENGTH) check_flag = YES;
					else check_flag = NO;
				}
				else {
					if(x>cell->collim.rmin && x<cell->collim.rmax && y>0 && y<cell->boun.ymax && fabs(z-z_opal)<TINYLENGTH) check_flag = YES;
					else check_flag = NO;
				}
				if(check_flag == YES) {
					if(get_bfield(x, y, z, &bx, &by, &bz, cell, add_contribution_comp)!=ALIVE) errorstop("bfield not alive");
					if(fabs(bx)>1.e-8) difbx = (bx_opal/10. - bx);///bx;
					else difbx=0;
					if(fabs(by)>1.e-8) difby = (by_opal/10. - by);///by;
					else difby=0;
					if(fabs(bz)>1.e-8) difbz = (bz_opal/10. - bz);///bz;
					else difbz=0;
					//if(x>4.7 && x<4.72 && y>0.5 && y<0.6) printf("%lf, %lf, %lf: %le, %le, %le: %le, %le, %le\n",x,y,z,bx,by,bz, difbx,difby,difbz);
					difminx = MIN(difbx,difminx);
					difmaxx = MAX(difbx,difmaxx);
					difminy = MIN(difby,difminy);
					difmaxy = MAX(difby,difmaxy);
					difminz = MIN(difbz,difminz);
					difmaxz = MAX(difbz,difmaxz);
					fprintf(wfilex, "%le	%le	%le\n", y,x,difbx);
					fprintf(wfiley, "%le	%le	%le\n", y,x,difby);
					fprintf(wfilez, "%le	%le	%le\n", y,x,difbz);
				}
			}
		}
		fprintf(wfilex,"\n");
		fprintf(wfiley,"\n");
		fprintf(wfilez,"\n");
	}
	printf("difx(min,max) = (%le,%le)\n",difminx,difmaxx);
	printf("dify(min,max) = (%le,%le)\n",difminy,difmaxy);
	printf("difz(min,max) = (%le,%le)\n",difminz,difmaxz);
	//printf("(x,y)max = (%le,%le)\n",xmax,ymax);
	fclose(rfile);
	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
	sprintf(outx,"%s_x.eps",fileout);
	sprintf(outy,"%s_y.eps",fileout);
	sprintf(outz,"%s_z.eps",fileout);
	//easyplot3dmap(namex, "y [m]", "x [m]", "dif Bx", NULL, NULL, outx, "size ratio -1",difminx,difmaxx);
	//easyplot3dmap(namey, "y [m]", "x [m]", "dif By", NULL, NULL, outy, "size ratio -1",difminy,difmaxy);
	//easyplot3dmap(namez, "y [m]", "x [m]", "dif Bz", NULL, NULL, outz, "size ratio -1",difminz,difmaxz);
	//easyplot3dmap_plusfile("data/aspectview.dat", "2", "1", "lines lt 1 lc 0 lw 2", namex, "y (long) [m]", "x (hor) [m]", "dif Bx [T]", NULL, NULL, outx, "size ratio -1",-1.,1.);//difminx,difmaxx);
	//easyplot3dmap_plusfile("data/aspectview.dat", "2", "1", "lines lt 1 lc 0 lw 2", namey, "y (long) [m]", "x (hor) [m]", "dif By [T]", NULL, NULL, outy, "size ratio -1",-1.,1.);//difminy,difmaxy);
	//easyplot3dmap_plusfile("data/aspectview.dat", "2", "1", "lines lt 1 lc 0 lw 2", namez, "y (long) [m]", "x (hor) [m]", "dif Bz [T]", NULL, NULL, outz, "size ratio -1",-1.,1.);//difminz,difmaxz);
	
}

extern void create_trackout_from_opal_file(char *opal_ori, char *trackoutfile, double rmin, double rmax, double th0_deg, double thmax_deg)
{
	char idpart[10];
	int i, nblines, jump_head = 2;
	double x,xp,y,yp,z,zp,r,th;
	FILE *rfile=NULL;
	FILE *wfile;
	
	nblines = get_nb_lines_file(opal_ori);
	rfile = fopen(opal_ori,"r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(i=0;i<jump_head;i++) newline(rfile);
	wfile = fopen(trackoutfile, "w");
	
	for(i=0;i<nblines-jump_head;i++) {
		fscanf(rfile,"%s %le %le %le %le %le %le", idpart, &x,&xp,&y,&yp,&z,&zp);
		r = sqrt(x*x + y*y);
		th = atan_ratio(y, x);
		//printf("i=%i: %s %le %le %le %le %le %le\n", i, idpart, x,xp,y,yp,z,zp);
		//printf("i=%i, r=%le,th=%le\n", i, r, th);
		if(strcmp(idpart, "ID1") == 0) {
			if(r>rmin && r<rmax && th>th0_deg*PI/180. && th<thmax_deg*PI/180.) {
				fprintf(wfile,"%le	%le	%le	%le	%le\n", x,y,z,r,th);
			}
		}
	}
	fclose(rfile);
	fclose(wfile);
}

extern void adjust_map_opal(char *opal_origin_file, char *file_for_fixfield, int nb_skipped_lines, int nbx, int nby, int nbz)
{
	char buf[MAX_CHARINLINE];
	int i,nblines;
	double x, y, z, t, bx, by, bz, ex, ey, ez;
	FILE *rfile=NULL;
	FILE *wfile;
	
	nblines = get_nb_lines_file(opal_origin_file);
	printf("in file %s, # lines: %i\n", opal_origin_file, nblines);
	if(nblines-nb_skipped_lines-nbx*nby*nbz !=0) errorstop("strange nb of lines");
	nblines -= nb_skipped_lines;
	rfile = fopen(opal_origin_file, "r");
	if(rfile==NULL) errorstop("cannot open opal origin file");
	wfile = fopen(file_for_fixfield,"w");
	
	for(i=0;i<nb_skipped_lines;i++) {
		if(fgets(buf, MAX_CHARINLINE-1, rfile) == NULL) errorstop("strange in adjust_map_opal");
		if(i==0) fprintf(wfile, "%i	%i	%i\n", nbx,nby,nbz);
		else fprintf(wfile, "%s", buf);
	}
	for(i=0;i<nblines;i++) {
		fscanf(rfile, "%le	%le	%le	%le	%le	%le	%le	%le	%le	%le", &x, &y, &z, &t, &bx, &by, &bz, &ex, &ey, &ez);
		fprintf(wfile, "%le	%le	%le	%le	%le	%le\n", x, y, z, bx, by, bz);
	}
	fclose(rfile);
	fclose(wfile);
}

extern void plot_field_map_opal(char *mapname)
{
	char namex[300], namey[300], namez[300], outx[300], outy[300], outz[300];
	int i,n=103042;
	double x,y,z,bx,by,bz;
	double bxmin,bxmax,bymin,bymax,bzmin,bzmax;
	FILE *rfile=NULL;
	FILE *wfilex,*wfiley,*wfilez;
	rfile=fopen(mapname,"r");
	if(rfile==NULL) errorstop("cannot open mapfile");
	sprintf(namex,"data/map_opal_x.dat");
	sprintf(namey,"data/map_opal_y.dat");
	sprintf(namez,"data/map_opal_z.dat");
	sprintf(outx,"output/map_opal_x.eps");
	sprintf(outy,"output/map_opal_y.eps");
	sprintf(outz,"output/map_opal_z.eps");
	
	wfilex=fopen(namex, "w");
	wfiley=fopen(namey, "w");
	wfilez=fopen(namez, "w");
	bxmin = +1.e9;
	bxmax = -1.e9;
	bymin = +1.e9;
	bymax = -1.e9;
	bzmin = +1.e9;
	bzmax = -1.e9;
	
	for(i=0;i<n;i++) {
		fscanf(rfile,"%lf	%lf	%lf	%lf	%lf	%lf", &x,&y,&z,&bx,&by,&bz);
		fprintf(wfilex,"%le\t%le\t%le\n", y, x, bx);
		if(bx>bxmax) bxmax = bx;
		if(bx<bxmin) bxmin = bx;
		
		fprintf(wfiley,"%le\t%le\t%le\n", y, x, by);
		if(by>bymax) bymax = by;
		if(by<bymin) bymin = by;
		
		fprintf(wfilez,"%le\t%le\t%le\n", y, x, bz);
		if(bz>bzmax) bzmax = bz;
		if(bz<bzmin) bzmin = bz;
		if(fabs(y-5.33333)<TINYLENGTH) {
			fprintf(wfilex,"\n");
			fprintf(wfiley,"\n");
			fprintf(wfilez,"\n");
		}
	}
	printf("bxmin,max=%le,%le\n",bxmin,bxmax);
	printf("bymin,max=%le,%le\n",bymin,bymax);
	printf("bzmin,max=%le,%le\n",bzmin,bzmax);
	fclose(rfile);
	fclose(wfilex);
	fclose(wfiley);
	fclose(wfilez);
	easyplot3dmap(namex, "x [m]", "y [m]", "bx [T]", NULL, NULL, outx, "size ratio -1", bxmin,bxmax);
	easyplot3dmap(namey, "x [m]", "y [m]", "by [T]", NULL, NULL, outy, "size ratio -1", bymin,bymax);
	easyplot3dmap(namez, "x [m]", "y [m]", "bz [T]", NULL, NULL, outz, "size ratio -1", bzmin,bzmax);
}

extern void para_vffa_fodo_scode_to_fixfield(double cell_length, double theta_cell_deg, double magnet_length, double small_drift, double shift, double tilt_d_deg, double tilt_f_deg, double b0d, double b0f, double fringe_add, double m, double z0, int order_expansion)
{
	double r0d, thetad1_deg, thetad2_deg, phid1_deg, phid2_deg, b0d_ff, fringe_ff;
	double r0f, thetaf1_deg, thetaf2_deg, phif1_deg, phif2_deg, b0f_ff;
	double ac, oh, ok, d1d2, f1f2, d1oh, f1ok, quarter_theta_cell, tilt_d, tilt_f;//, oa;
	
	ac = cell_length/2.;
	quarter_theta_cell = theta_cell_deg*PI/180./4.;
	tilt_d = tilt_d_deg*PI/180.;
	tilt_f = tilt_f_deg*PI/180.;
	d1d2 = magnet_length+small_drift;
	f1f2 = d1d2;
	//oa = ac/(2*sin(quarter_theta_cell));
	oh = ac/(2*tan(quarter_theta_cell))-shift;
	d1oh = atan(d1d2/(2.*oh));
	ok = ac/(2*tan(quarter_theta_cell))+shift;
	f1ok = atan(f1f2/(2.*ok));
	
	r0d = 2*oh*oh/(sqrt(4*oh*oh+d1d2*d1d2));
	thetad1_deg = (quarter_theta_cell - d1oh)*180./PI;
	thetad2_deg = (quarter_theta_cell + d1oh)*180./PI;
	phid1_deg = (tilt_d + d1oh)*180./PI;
	phid2_deg = -phid1_deg;
	
	r0f = 2*ok*ok/(sqrt(4*ok*ok+f1f2*f1f2));
	thetaf1_deg = (quarter_theta_cell*3. - f1ok)*180./PI;
	thetaf2_deg = (quarter_theta_cell*3. + f1ok)*180./PI;
	phif1_deg = (-tilt_f + f1ok)*180./PI;
	phif2_deg = -phif1_deg;
	
	find_b0_add_to_mult(b0d, -magnet_length/2., magnet_length/2., fringe_add, &b0d_ff, 1.e-7);
	find_b0_add_to_mult(b0f, -magnet_length/2., magnet_length/2., fringe_add, &b0f_ff, 1.e-7);
	fringe_ff = fringe_add/2.;
	
	printf("(F1, F2), D1, D2, F1, F2, (D1, D2):\n"); //-29.26521	3.789131	 0.21103232		0	1.58	6.26521
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf1_deg-theta_cell_deg, r0f, b0f_ff, z0, m, phif1_deg); //F1
	printf("\t%lf\t%lf\t%i\t1\n", -magnet_length/2., fringe_ff, order_expansion); //F1
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet_length/2., fringe_ff, order_expansion); //F1
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf2_deg-theta_cell_deg, r0f, b0f_ff, z0, m, phif2_deg); //F2
	printf("\t%lf\t%lf\t%i\t1\n", -magnet_length/2., fringe_ff, order_expansion); //F2
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet_length/2., fringe_ff, order_expansion); //F2
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetad1_deg, r0d, b0d_ff, z0, m, phid1_deg); //D1
	printf("\t%lf\t%lf\t%i\t1\n", -magnet_length/2., fringe_ff, order_expansion); //D1
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet_length/2., fringe_ff, order_expansion); //D1
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetad2_deg, r0d, b0d_ff, z0, m, phid2_deg); //D2
	printf("\t%lf\t%lf\t%i\t1\n", -magnet_length/2., fringe_ff, order_expansion); //D2
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet_length/2., fringe_ff, order_expansion); //D2
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf1_deg, r0f, b0f_ff, z0, m, phif1_deg); //F1
	printf("\t%lf\t%lf\t%i\t1\n", -magnet_length/2., fringe_ff, order_expansion); //F1
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet_length/2., fringe_ff, order_expansion); //F1
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf2_deg, r0f, b0f_ff, z0, m, phif2_deg); //F2
	printf("\t%lf\t%lf\t%i\t1\n", -magnet_length/2., fringe_ff, order_expansion); //F2
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet_length/2., fringe_ff, order_expansion); //F2
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetad1_deg+theta_cell_deg, r0d, b0d_ff, z0, m, phid1_deg); //D1
	printf("\t%lf\t%lf\t%i\t1\n", -magnet_length/2., fringe_ff, order_expansion); //D1
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet_length/2., fringe_ff, order_expansion); //D1
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetad2_deg+theta_cell_deg, r0d, b0d_ff, z0, m, phid2_deg); //D2
	printf("\t%lf\t%lf\t%i\t1\n", -magnet_length/2., fringe_ff, order_expansion); //D2
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet_length/2., fringe_ff, order_expansion); //D2
	
}

extern void para_vffa_triplet_scode_to_fixfield(double cell_length, double theta_cell_deg, double magnet_length, double small_drift, double shift, double tilt_f_deg, double b0d, double b0f, double fringe_add, double m, double z0, int order_expansion)
{
	double r0d, thetad, thetad_deg, phid_deg, b0d_ff, fringe_ff;
	double r0f, thetaf1, thetaf1_deg, thetaf2_deg, phif1_deg, phif2_deg, b0f_ff;
	double ac, half_theta_cell, eb, oa, ob;//, tilt_f;
	
	ac = cell_length;
	half_theta_cell = theta_cell_deg*PI/180./2.;
	//tilt_f = tilt_f_deg*PI/180.;
	eb = magnet_length+small_drift;
	oa = ac/(2*sin(half_theta_cell));
	ob = ac/(2*tan(half_theta_cell));
	
	thetad = half_theta_cell;
	thetad_deg = thetad*180./PI;
	r0d = ob-shift;
	phid_deg = 0.;
	
	thetaf1 = thetad - atan(eb/(ob+shift));
	thetaf1_deg = thetaf1*180./PI;
	r0f = (ob+shift)/(cos(thetad-thetaf1));
	phif1_deg = tilt_f_deg + thetad_deg - thetaf1_deg;
	
	thetaf2_deg = theta_cell_deg - thetaf1_deg;
	phif2_deg = -phif1_deg;
	
	find_b0_add_to_mult(b0d, -magnet_length/2., magnet_length/2., fringe_add, &b0d_ff, 1.e-7);
	find_b0_add_to_mult(b0f, -magnet_length/2., magnet_length/2., fringe_add, &b0f_ff, 1.e-7);
	fringe_ff = fringe_add/2.;
	printf("OA = %le\n", oa);
	printf("F1, D, F2:\n"); //-29.26521	3.789131	 0.21103232		0	1.58	6.26521
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf1_deg, r0f, b0f_ff, z0, m, phif1_deg); //F1
	printf("\t%lf\t%lf\t%i\t1\n", -magnet_length/2., fringe_ff, order_expansion); //F1
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet_length/2., fringe_ff, order_expansion); //F1
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetad_deg, r0d, b0d_ff, z0, m, phid_deg); //D
	printf("\t%lf\t%lf\t%i\t1\n", -magnet_length/2., fringe_ff, order_expansion); //D
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet_length/2., fringe_ff, order_expansion); //D
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf2_deg, r0f, b0f_ff, z0, m, phif2_deg); //F2
	printf("\t%lf\t%lf\t%i\t1\n", -magnet_length/2., fringe_ff, order_expansion); //F2
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet_length/2., fringe_ff, order_expansion); //F2
	
}

//corrected r0d with shift (wrong sign)
extern void para_vffa_triplet_scode_to_fixfield2(double cell_length, double theta_cell_deg, double magnet13_length, double small_drift, double magnet2_length, double shift, double tilt_mag1_deg, double b0f, double b0d, double fringe_add, double m, double z0, int order_expansion)
{
	double r0d, thetad, thetad_deg, phid_deg, b0d_ff, fringe_ff;
	double r0f, thetaf1, thetaf1_deg, thetaf2_deg, phif1_deg, phif2_deg, b0f_ff;
	double ac, half_theta_cell, eb, oa, ob;//, tilt_f;
	
	ac = cell_length;
	half_theta_cell = theta_cell_deg*PI/180./2.;
	//tilt_f = tilt_mag1_deg*PI/180.;
	eb = magnet13_length/2. + magnet2_length/2. + small_drift;
	oa = ac/(2*sin(half_theta_cell));
	ob = ac/(2*tan(half_theta_cell));
	
	thetad = half_theta_cell;
	thetad_deg = thetad*180./PI;
	r0d = ob+shift;
	phid_deg = 0.;
	
	thetaf1 = thetad - atan(eb/(ob+shift));
	thetaf1_deg = thetaf1*180./PI;
	r0f = (ob+shift)/(cos(thetad-thetaf1));
	phif1_deg = tilt_mag1_deg + thetad_deg - thetaf1_deg;
	
	thetaf2_deg = theta_cell_deg - thetaf1_deg;
	phif2_deg = -phif1_deg;
	
	find_b0_add_to_mult(b0d, -magnet2_length/2., magnet2_length/2., fringe_add, &b0d_ff, 1.e-8);
	find_b0_add_to_mult(b0f, -magnet13_length/2., magnet13_length/2., fringe_add, &b0f_ff, 1.e-8);
	fringe_ff = fringe_add/2.;
	printf("OA = %le\n", oa);
	printf("F1, D, F2:\n"); //-29.26521	3.789131	 0.21103232		0	1.58	6.26521
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf1_deg, r0f, b0f_ff, z0, m, phif1_deg); //F1
	printf("\t%lf\t%lf\t%i\t1\n", -magnet13_length/2., fringe_ff, order_expansion); //F1
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet13_length/2., fringe_ff, order_expansion); //F1
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetad_deg, r0d, b0d_ff, z0, m, phid_deg); //D
	printf("\t%lf\t%lf\t%i\t1\n", -magnet2_length/2., fringe_ff, order_expansion); //D
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet2_length/2., fringe_ff, order_expansion); //D
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf2_deg, r0f, b0f_ff, z0, m, phif2_deg); //F2
	printf("\t%lf\t%lf\t%i\t1\n", -magnet13_length/2., fringe_ff, order_expansion); //F2
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet13_length/2., fringe_ff, order_expansion); //F2
	
}

extern void para_vffa_triplet_scode_to_fixfield3(double cell_length, double theta_cell_deg, double magnet13_length, double small_drift, double magnet2_length, double shiftf, double shiftd, double tilt_mag1_deg, double b0f, double b0d, double fringe_add, double m, double z0, int order_expansion)
{
	double r0d, thetad, thetad_deg, phid_deg, b0d_ff, fringe_ff;
	double r0f, thetaf1, thetaf1_deg, thetaf2_deg, phif1_deg, phif2_deg, b0f_ff;
	double ac, half_theta_cell, eb, oa, ob;//, tilt_f;
	
	ac = cell_length;
	half_theta_cell = theta_cell_deg*PI/180./2.;
	//tilt_f = tilt_mag1_deg*PI/180.;
	eb = magnet13_length/2. + magnet2_length/2. + small_drift;
	oa = ac/(2*sin(half_theta_cell));
	ob = ac/(2*tan(half_theta_cell));//=r_D
	
	thetad = half_theta_cell;
	thetad_deg = thetad*180./PI;
	r0d = ob+shiftd;
	phid_deg = 0.;
	
	thetaf1 = thetad - atan(eb/(ob+shiftf));
	thetaf1_deg = thetaf1*180./PI;
	r0f = (ob+shiftf)/(cos(thetad-thetaf1));
	phif1_deg = tilt_mag1_deg + thetad_deg - thetaf1_deg;
	
	thetaf2_deg = theta_cell_deg - thetaf1_deg;
	phif2_deg = -phif1_deg;
	
	//find_b0_add_to_mult(b0d, -magnet2_length/2., magnet2_length/2., fringe_add, &b0d_ff, 1.e-8);
	b0d_ff = find_b0_add_to_mult2(b0d, -magnet2_length/2., magnet2_length/2., fringe_add);
	//find_b0_add_to_mult(b0f, -magnet13_length/2., magnet13_length/2., fringe_add, &b0f_ff, 1.e-8);
	b0f_ff = find_b0_add_to_mult2(b0f, -magnet13_length/2., magnet13_length/2., fringe_add);
	fringe_ff = fringe_add/2.;
	printf("OA = %le\n", oa);
	printf("F2, F1, D, F2, F1:\n"); //-29.26521	3.789131	 0.21103232		0	1.58	6.26521
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf2_deg-theta_cell_deg, r0f, b0f_ff, z0, m, phif2_deg); //F2
	printf("\t%lf\t%lf\t%i\t1\n", -magnet13_length/2., fringe_ff, order_expansion); //F2
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet13_length/2., fringe_ff, order_expansion); //F2
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf1_deg, r0f, b0f_ff, z0, m, phif1_deg); //F1
	printf("\t%lf\t%lf\t%i\t1\n", -magnet13_length/2., fringe_ff, order_expansion); //F1
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet13_length/2., fringe_ff, order_expansion); //F1
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetad_deg, r0d, b0d_ff, z0, m, phid_deg); //D
	printf("\t%lf\t%lf\t%i\t1\n", -magnet2_length/2., fringe_ff, order_expansion); //D
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet2_length/2., fringe_ff, order_expansion); //D
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf2_deg, r0f, b0f_ff, z0, m, phif2_deg); //F2
	printf("\t%lf\t%lf\t%i\t1\n", -magnet13_length/2., fringe_ff, order_expansion); //F2
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet13_length/2., fringe_ff, order_expansion); //F2
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf1_deg+theta_cell_deg, r0f, b0f_ff, z0, m, phif1_deg); //F1
	printf("\t%lf\t%lf\t%i\t1\n", -magnet13_length/2., fringe_ff, order_expansion); //F1
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet13_length/2., fringe_ff, order_expansion); //F1
	
}

// report 200203_vffa_fixfield
extern void para_vffa_triplet_scode_to_fixfield_add(double cell_length, double theta_cell_deg, double magnet13_length, double small_drift, double magnet2_length, double shiftf, double shiftd, double tilt_mag1_deg, double b0f, double b0d, double fringe_add, double m, double z0, int order_expansion)
{
	double r0d, thetad, thetad_deg, phid_deg, fringe_ff;
	double r0f, thetaf1, thetaf1_deg, thetaf2_deg, phif1_deg, phif2_deg;
	double ac, half_theta_cell, eb, oa, ob;//, tilt_f;
	
	ac = cell_length;
	half_theta_cell = theta_cell_deg*PI/180./2.;
	//tilt_f = tilt_mag1_deg*PI/180.;
	eb = magnet13_length/2. + magnet2_length/2. + small_drift;
	oa = ac/(2*sin(half_theta_cell));
	ob = ac/(2*tan(half_theta_cell));
	
	thetad = half_theta_cell;
	thetad_deg = thetad*180./PI;
	r0d = ob+shiftd;
	phid_deg = 0.;
	
	thetaf1 = thetad - atan(eb/(ob+shiftf));
	thetaf1_deg = thetaf1*180./PI;
	r0f = (ob+shiftf)/(cos(thetad-thetaf1));
	phif1_deg = tilt_mag1_deg + thetad_deg - thetaf1_deg;
	
	thetaf2_deg = theta_cell_deg - thetaf1_deg;
	phif2_deg = -phif1_deg;
	
	fringe_ff = fringe_add/2.;
	printf("OA = %le\n", oa);
	//printf("F2, F1, D, F2, F1:\n"); //-29.26521	3.789131	 0.21103232		0	1.58	6.26521
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf2_deg-theta_cell_deg, r0f, b0f, z0, m, phif2_deg); //F2
	printf("\t%lf\t%lf\t%i\t1\n", -magnet13_length/2., fringe_ff, order_expansion); //F2
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet13_length/2., fringe_ff, order_expansion); //F2
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf1_deg, r0f, b0f, z0, m, phif1_deg); //F1
	printf("\t%lf\t%lf\t%i\t1\n", -magnet13_length/2., fringe_ff, order_expansion); //F1
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet13_length/2., fringe_ff, order_expansion); //F1
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetad_deg, r0d, b0d, z0, m, phid_deg); //D
	printf("\t%lf\t%lf\t%i\t1\n", -magnet2_length/2., fringe_ff, order_expansion); //D
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet2_length/2., fringe_ff, order_expansion); //D
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf2_deg, r0f, b0f, z0, m, phif2_deg); //F2
	printf("\t%lf\t%lf\t%i\t1\n", -magnet13_length/2., fringe_ff, order_expansion); //F2
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet13_length/2., fringe_ff, order_expansion); //F2
	
	printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", thetaf1_deg+theta_cell_deg, r0f, b0f, z0, m, phif1_deg); //F1
	printf("\t%lf\t%lf\t%i\t1\n", -magnet13_length/2., fringe_ff, order_expansion); //F1
	printf("\t%lf\t%lf\t%i\t1\n\n", magnet13_length/2., fringe_ff, order_expansion); //F1
	
}

// get the b0 for multiplication of the fringe field from parameters of addition of the fringe field, lambda_add=2*lambda_FF
extern int find_b0_add_to_mult(double b0_add, double efbe, double efbs, double lambda_add, double *b0_mult, double prec)
{
	int n;
	double x, dif;
	x = 0.5*(efbe+efbs);
	*b0_mult = b0_add;
	for(n=0;n<100;n++) {
		dif = *b0_mult/((1+exp(2.*(x-efbs)/lambda_add))*(1+exp(2.*(efbe-x)/lambda_add)))-b0_add/2.*(tanh((x-efbe)/lambda_add) + tanh((efbs-x)/lambda_add));
		if(fabs(dif)<prec) {
			printf("\n\tb0 found: %.8f, b0_add = %.8f\n", *b0_mult, b0_add);
			return TRUE;
		}
		*b0_mult -= dif;
		//printf("n=%i, dif = %le, b0_mult=%le\n", n, dif, *b0_mult);
	}
	printf("\n\tb0 not found: dif=%le, b0_add = %le\n", dif, b0_add);
	return FALSE;
}

//b0 for multiplication of the fringe field from parameters of addition of the fringe field, lambda_add=2*lambda_FF
extern double find_b0_add_to_mult2(double b0_add, double efbe, double efbs, double lambda_add)
{
	double b0_mult, lambda_mult = lambda_add/2.;
	
	b0_mult = b0_add*(1+exp((efbe-efbs)/(2*lambda_mult)))*(1+exp((efbe-efbs)/(2*lambda_mult)))*tanh((efbs-efbe)/(2*lambda_add));
	printf("\n\tb0_mult: %.8f, b0_add = %.8f\n", b0_mult, b0_add);
	return b0_mult;
}

extern void add_theta_opal_file(struct Cell *cell, void(*add_contribution_comp)(double,double,double,double*,double*,double*,struct Cell*,int), char *opal_ori, char *newfile)
{
	char idpart[10];
	int i, nblines, jump_head = 2;
	double x,xp,y,yp,z,zp,r,th,bx,by,bz,ex,ey,ez;
	double bx_home,by_home,bz_home;
	FILE *rfile=NULL;
	FILE *wfile;
	FILE *wfile2;
	
	nblines = get_nb_lines_file(opal_ori);
	rfile = fopen(opal_ori,"r");
	if(rfile==NULL) errorstop("cannot open rfile");
	for(i=0;i<jump_head;i++) newline(rfile);
	wfile = fopen(newfile, "w");
	wfile2 = fopen("data/opal_erit_coor_homemade_field.dat", "w");
	
	for(i=0;i<nblines-jump_head;i++) {
		fscanf(rfile,"%s %le %le %le %le %le %le %le %le %le %le %le %le", idpart, &x,&xp,&y,&yp,&z,&zp,&bx,&by,&bz,&ex,&ey,&ez);
		r = sqrt(x*x + y*y);
		th = atan_ratio(y, x);
		printf("i=%i, r=%le,th=%le\n", i, r, th);
		if(r>cell->collim.rmin && r<cell->collim.rmax && th<cell->boun.thmax) {
			fprintf(wfile,"%le	%le	%le	%le	%le	%le	%le	%le\n", x,y,z,r,th,bx,by,bz);
			if(get_bfield(x, y, z, &bx_home, &by_home, &bz_home, cell, add_contribution_comp)!=ALIVE) errorstop("bfield not alive");
			fprintf(wfile2,"%le	%le	%le	%le	%le	%le	%le	%le\n", x,y,z,r,th,bx_home,by_home,bz_home);
		}
		else break;
	}
	fclose(rfile);
	fclose(wfile);
	fclose(wfile2);
}

