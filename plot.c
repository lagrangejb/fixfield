/*
 *  plot.c
 *  ringdesign
 *
 *  Copyright 2009 Kyoto University. All rights reserved.
 *
 */

#include "plot.h"

extern void easyplot(char *txtfilename, char *xcolumn, char *ycolumn, char *with, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,30.'\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set lmargin 15\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,30.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -5,0 font 'helvetica,30.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color\n");
		fprintf(gp," set output '%s' \n", psfilename);
		fprintf(gp,"plot '%s' using %s:%s with %s\n", txtfilename, xcolumn, ycolumn, with);
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using %s:%s with %s\n", txtfilename, xcolumn, ycolumn, with);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot2(char *txtfilename1, char *txtfilename2, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *with1, char *with2, char *title1, char *title2, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	gp=popen("gnuplot -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,30.'\n");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,30.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -5,0 font 'helvetica,30.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	if(title1==NULL && title2==NULL) fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color\n");
		fprintf(gp," set output '%s' \n", psfilename);
		if(title1==NULL) fprintf(gp,"plot '%s' using %s:%s with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, with1);
		else fprintf(gp,"plot '%s' using %s:%s with %s title '%s', ", txtfilename1, xcolumn1, ycolumn1, with1, title1);
		if(title2==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle\n", txtfilename2, xcolumn2, ycolumn2, with2);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s'\n", txtfilename2, xcolumn2, ycolumn2, with2, title2);
		//if(title1==NULL && title2==NULL) fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2);
		//else if(title1==NULL) {
		//	fprintf(gp,"plot '%s' using %s:%s with %s notitle, '%s' using %s:%s with %s title '%s'\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2,title2);
		//}
		//else if(title2==NULL) {
		//	fprintf(gp,"plot '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s notitle\n", txtfilename1, xcolumn1, ycolumn1, with1, title1, txtfilename2, xcolumn2, ycolumn2, with2);
		//}
		//else {	
		//	fprintf(gp,"plot '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s'\n", txtfilename1, xcolumn1, ycolumn1, with1, title1, txtfilename2, xcolumn2, ycolumn2, with2,title2);
		//}
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot3(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *with1, char *with2, char *with3, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	gp=popen("gnuplot -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color\n");
		fprintf(gp," set output '%s' \n", psfilename);
		fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot3p(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *with1, char *with2, char *with3, char *title1, char *title2, char *title3, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//fprintf(gp, "set colorsequence classic\n");
	gp=popen("gnuplot -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");	
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color\n");
		fprintf(gp," set output '%s' \n", psfilename);
		//fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
		//fprintf(gp,"plot '%s' using %s:%s with %s notitle, '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s'\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, title2, txtfilename3, xcolumn3, ycolumn3, with3, title3);
		if(title1==NULL) fprintf(gp,"plot '%s' using %s:%s with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, with1);
		else fprintf(gp,"plot '%s' using %s:%s with %s title '%s', ", txtfilename1, xcolumn1, ycolumn1, with1, title1);
		if(title2==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename2, xcolumn2, ycolumn2, with2);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename2, xcolumn2, ycolumn2, with2, title2);
		if(title3==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle\n", txtfilename3, xcolumn3, ycolumn3, with3);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s'\n", txtfilename3, xcolumn3, ycolumn3, with3, title3);
		//fprintf(gp,"plot '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s'\n", txtfilename1, xcolumn1, ycolumn1, with1, title1, txtfilename2, xcolumn2, ycolumn2, with2, title2, txtfilename3, xcolumn3, ycolumn3, with3, title3);
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
	//fprintf(gp,"plot '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s'\n", txtfilename1, xcolumn1, ycolumn1, with1, title1, txtfilename2, xcolumn2, ycolumn2, with2, title2, txtfilename3, xcolumn3, ycolumn3, with3, title3);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot4(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *xcolumn4, char *ycolumn4, char *with1, char *with2, char *with3, char *with4, char *title1, char *title2, char *title3, char *title4, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	gp=popen("gnuplot -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	if(title1 == NULL && title2 == NULL && title3 == NULL && title4 == NULL) fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color \n");
		fprintf(gp," set output '%s' \n", psfilename);
		if(title1==NULL) fprintf(gp,"plot '%s' using %s:%s with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, with1);
		else fprintf(gp,"plot '%s' using %s:%s with %s title '%s', ", txtfilename1, xcolumn1, ycolumn1, with1, title1);
		if(title2==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename2, xcolumn2, ycolumn2, with2);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename2, xcolumn2, ycolumn2, with2, title2);
		if(title3==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename3, xcolumn3, ycolumn3, with3);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename3, xcolumn3, ycolumn3, with3, title3);
		if(title4==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle\n", txtfilename4, xcolumn4, ycolumn4, with4);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s'\n", txtfilename4, xcolumn4, ycolumn4, with4, title4);
		
		//if(title1 == NULL && title2 == NULL && title3 == NULL && title4 == NULL) fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3, txtfilename4, xcolumn4, ycolumn4, with4);
		//else if(title1 == NULL) fprintf(gp,"plot '%s' using %s:%s with %s notitle, '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s'\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, title2, txtfilename3, xcolumn3, ycolumn3, with3, title3, txtfilename4, xcolumn4, ycolumn4, with4, title4);
		//else fprintf(gp,"plot '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s'\n", txtfilename1, xcolumn1, ycolumn1, with1, title1, txtfilename2, xcolumn2, ycolumn2, with2, title2, txtfilename3, xcolumn3, ycolumn3, with3, title3, txtfilename4, xcolumn4, ycolumn4, with4, title4);
		//fprintf(gp," set output\n");
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3, txtfilename4, xcolumn4, ycolumn4, with4);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot5(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *txtfilename5, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *xcolumn4, char *ycolumn4, char *xcolumn5, char *ycolumn5, char *with1, char *with2, char *with3, char *with4, char *with5, char *title1, char *title2, char *title3, char *title4, char *title5, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	gp=popen("gnuplot -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	if(title1 == NULL && title2 == NULL && title3 == NULL && title4 == NULL && title5 == NULL)fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color \n");
		fprintf(gp," set output '%s' \n", psfilename);
		if(title1==NULL) fprintf(gp,"plot '%s' using %s:%s with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, with1);
		else fprintf(gp,"plot '%s' using %s:%s with %s title '%s', ", txtfilename1, xcolumn1, ycolumn1, with1, title1);
		if(title2==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename2, xcolumn2, ycolumn2, with2);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename2, xcolumn2, ycolumn2, with2, title2);
		if(title3==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename3, xcolumn3, ycolumn3, with3);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename3, xcolumn3, ycolumn3, with3, title3);
		if(title4==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename4, xcolumn4, ycolumn4, with4);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename4, xcolumn4, ycolumn4, with4, title4);
		if(title5==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle\n", txtfilename5, xcolumn5, ycolumn5, with5);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s'\n", txtfilename5, xcolumn5, ycolumn5, with5, title5);
		//if(title1 == NULL && title2 == NULL && title3 == NULL && title4 == NULL && title5 == NULL) fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3, txtfilename4, xcolumn4, ycolumn4, with4, txtfilename5, xcolumn5, ycolumn5, with5);
		//else fprintf(gp,"plot '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s'\n", txtfilename1, xcolumn1, ycolumn1, with1, title1, txtfilename2, xcolumn2, ycolumn2, with2, title2, txtfilename3, xcolumn3, ycolumn3, with3, title3, txtfilename4, xcolumn4, ycolumn4, with4, title4, txtfilename5, xcolumn5, ycolumn5, with5, title5);
		//fprintf(gp," set output\n");
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s , '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3, txtfilename4, xcolumn4, ycolumn4, with4, txtfilename5, xcolumn5, ycolumn5, with5);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot6(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *txtfilename5, char *txtfilename6, 
char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *xcolumn4, char *ycolumn4, char *xcolumn5, char *ycolumn5, char *xcolumn6, char *ycolumn6, 
char *with1, char *with2, char *with3, char *with4, char *with5, char *with6, 
char *title1, char *title2, char *title3, char *title4, char *title5, char *title6, 
char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	gp=popen("gnuplot -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	if(title1 == NULL && title2 == NULL && title3 == NULL && title4 == NULL && title5 == NULL && title6 == NULL) fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color \n");
		fprintf(gp," set output '%s' \n", psfilename);
		if(title1==NULL) fprintf(gp,"plot '%s' using %s:%s with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, with1);
		else fprintf(gp,"plot '%s' using %s:%s with %s title '%s', ", txtfilename1, xcolumn1, ycolumn1, with1, title1);
		if(title2==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename2, xcolumn2, ycolumn2, with2);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename2, xcolumn2, ycolumn2, with2, title2);
		if(title3==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename3, xcolumn3, ycolumn3, with3);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename3, xcolumn3, ycolumn3, with3, title3);
		if(title4==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename4, xcolumn4, ycolumn4, with4);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename4, xcolumn4, ycolumn4, with4, title4);
		if(title5==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename5, xcolumn5, ycolumn5, with5);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename5, xcolumn5, ycolumn5, with5, title5);
		if(title6==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle\n", txtfilename6, xcolumn6, ycolumn6, with6);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s'\n", txtfilename6, xcolumn6, ycolumn6, with6, title6);
		//if(title1 == NULL && title2 == NULL && title3 == NULL && title4 == NULL && title5 == NULL  && title6 == NULL) fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3, txtfilename4, xcolumn4, ycolumn4, with4, txtfilename5, xcolumn5, ycolumn5, with5, txtfilename6, xcolumn6, ycolumn6, with6);
		//else fprintf(gp,"plot '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s'\n", txtfilename1, xcolumn1, ycolumn1, with1, title1, txtfilename2, xcolumn2, ycolumn2, with2, title2, txtfilename3, xcolumn3, ycolumn3, with3, title3, txtfilename4, xcolumn4, ycolumn4, with4, title4, txtfilename5, xcolumn5, ycolumn5, with5, title5, txtfilename6, xcolumn6, ycolumn6, with6, title6);
		//fprintf(gp," set output\n");
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3, txtfilename4, xcolumn4, ycolumn4, with4, txtfilename5, xcolumn5, ycolumn5, with5, txtfilename6, xcolumn6, ycolumn6, with6);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot8(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *txtfilename5, char *txtfilename6,  char *txtfilename7, char *txtfilename8,
char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *xcolumn4, char *ycolumn4, char *xcolumn5, char *ycolumn5, char *xcolumn6, char *ycolumn6, char *xcolumn7, char *ycolumn7, char *xcolumn8, char *ycolumn8, 
char *with1, char *with2, char *with3, char *with4, char *with5, char *with6,  char *with7, char *with8,
char *title1, char *title2, char *title3, char *title4, char *title5, char *title6, char *title7, char *title8, 
char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	gp=popen("gnuplot -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	if(title1 == NULL && title2 == NULL && title3 == NULL && title4 == NULL && title5 == NULL && title6 == NULL) fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color \n");
		fprintf(gp," set output '%s' \n", psfilename);
		if(title1==NULL) fprintf(gp,"plot '%s' using %s:%s with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, with1);
		else fprintf(gp,"plot '%s' using %s:%s with %s title '%s', ", txtfilename1, xcolumn1, ycolumn1, with1, title1);
		if(title2==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename2, xcolumn2, ycolumn2, with2);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename2, xcolumn2, ycolumn2, with2, title2);
		if(title3==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename3, xcolumn3, ycolumn3, with3);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename3, xcolumn3, ycolumn3, with3, title3);
		if(title4==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename4, xcolumn4, ycolumn4, with4);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename4, xcolumn4, ycolumn4, with4, title4);
		if(title5==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename5, xcolumn5, ycolumn5, with5);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename5, xcolumn5, ycolumn5, with5, title5);
		if(title6==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename6, xcolumn6, ycolumn6, with6);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename6, xcolumn6, ycolumn6, with6, title6);
		if(title7==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle, ", txtfilename7, xcolumn7, ycolumn7, with7);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s', ", txtfilename7, xcolumn7, ycolumn7, with7, title7);
		if(title8==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle\n", txtfilename8, xcolumn8, ycolumn8, with8);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s'\n", txtfilename8, xcolumn8, ycolumn8, with8, title8);
		//if(title1 == NULL && title2 == NULL && title3 == NULL && title4 == NULL && title5 == NULL  && title6 == NULL) fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3, txtfilename4, xcolumn4, ycolumn4, with4, txtfilename5, xcolumn5, ycolumn5, with5, txtfilename6, xcolumn6, ycolumn6, with6);
		//else fprintf(gp,"plot '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s', '%s' using %s:%s with %s title '%s'\n", txtfilename1, xcolumn1, ycolumn1, with1, title1, txtfilename2, xcolumn2, ycolumn2, with2, title2, txtfilename3, xcolumn3, ycolumn3, with3, title3, txtfilename4, xcolumn4, ycolumn4, with4, title4, txtfilename5, xcolumn5, ycolumn5, with5, title5, txtfilename6, xcolumn6, ycolumn6, with6, title6);
		//fprintf(gp," set output\n");
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s, '%s' using %s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3, txtfilename4, xcolumn4, ycolumn4, with4, txtfilename5, xcolumn5, ycolumn5, with5, txtfilename6, xcolumn6, ycolumn6, with6, txtfilename7, xcolumn7, ycolumn7, with7, txtfilename8, xcolumn8, ycolumn8, with8);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}
//old 3dmap
/*extern void easyplot3dmap(char *txtfilename, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,20.'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 10\n");
	fprintf(gp,"set pm3d map\n");
	fprintf(gp,"set palette rgbformulae 22,13,-31\n");
	//fprintf(gp,"set palette defined (-35 'blue', 0 'white', 15 'red')\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript color\n");
		fprintf(gp," set output '%s' \n", psfilename);
		fprintf(gp,"splot '%s'\n", txtfilename);
		fprintf(gp," set output\n");
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"splot '%s'\n", txtfilename);
	pclose(gp);
}//*/

//extern void easyplot3dmap(char *txtfilename, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption, int dens_max)
/*extern void easyplot3dmap(char *txtfilename, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_min, double dens_max)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,20.'\n");
	//fprintf(gp,"set bmargin 5\n");
	//fprintf(gp,"set rmargin 10\n");
	fprintf(gp,"set pm3d map\n");
	//fprintf(gp,"set palette defined (0 'white', 6 'yellow', 12 'orange', 18 'red')\n");//,dens_max/3., dens_max/2.,dens_max);
	//fprintf(gp,"set palette rgbformulae 33,13,10\n");
	//fprintf(gp,"set palette rgbformulae 22,13,-31\n");
	//fprintf(gp,"set palette defined (%lf 'blue', 0 'white', %lf 'red')\n", dens_min, dens_max);
	//if(dens_max>100) fprintf(gp,"set palette defined (0 'white', 5 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n", dens_max/5, dens_max*2/5, dens_max*3/5, dens_max*4/5, dens_max);
	//else if(dens_max>10) fprintf(gp,"set palette defined (0 'white', %i 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n", (dens_max/25+1), dens_max/5, dens_max*2/5, dens_max*3/5, dens_max*4/5, dens_max);
	//else if(dens_max>5) fprintf(gp,"set palette defined (0 'white', 1 'blue', %lf 'cyan', %lf 'green', %lf 'yellow', %lf 'orange', %i 'red')\n", (double) (dens_max/5.0), (double) (dens_max*2.0/5.0), (double) (dens_max*3.0/5.0), (double) (dens_max*4.0/5.0), dens_max);
	//else fprintf(gp,"set palette defined (0 'white', 1 'blue', 2 'green', 3 'yellow', 4 'orange', 5 'red')\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	//if(dens_max>5) fprintf(gp, " set cbrange [0:%i]\n", dens_max);
	//else fprintf(gp, " set cbrange [0:5]\n");
	//fprintf(gp, " set cbrange [%lf:%lf]\n", dens_min, dens_max);
	fprintf(gp, " set cbrange [%lf:%lf]\n", dens_min, dens_max);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		//fprintf(gp," set terminal postscript color\n");
		fprintf(gp," set terminal postscript eps enhanced color \n");
		fprintf(gp," set output '%s' \n", psfilename);
		fprintf(gp,"splot '%s'\n", txtfilename);
		fprintf(gp," set output\n");
		//fprintf(gp," set terminal x11\n");
	}
	//fprintf(gp,"splot '%s'\n", txtfilename);
	pclose(gp);
}//*/

extern void easyplot3dmap(char *txtfilename, char *xlabel, char *ylabel, char *cblabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_min, double dens_max)
{
	double scale;
	FILE *gp;
	scale = dens_max-dens_min;
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,20.'\n");
	//fprintf(gp,"set bmargin 5\n");
	//fprintf(gp,"set rmargin 10\n");
	fprintf(gp,"set pm3d map\n");
	//fprintf(gp,"set palette defined (0 'white', 6 'yellow', 12 'orange', 18 'red')\n");//,dens_max/3., dens_max/2.,dens_max);
	//fprintf(gp,"set palette rgbformulae 33,13,10\n");
	//fprintf(gp,"set palette rgbformulae 22,13,-31\n");
	if(dens_max>0) {
		if(dens_min<0 && fabs(dens_min)>scale/100) {
			if(fabs(dens_min)<scale/5. || dens_max<scale/5.) {
				fprintf(gp,"set palette defined (%le 'blue', 0 'white', %le 'red')\n", -scale, scale);
				fprintf(gp, " set cbrange [%le:%le]\n", -scale, scale);
				//printf("palette defined (%le 'blue', 0 'white', %le 'red')\n", -scale, scale);
			}
			else {
				fprintf(gp,"set palette defined (%le 'blue', 0 'white', %le 'red')\n", dens_min, dens_max);
				fprintf(gp, " set cbrange [%le:%le]\n", dens_min, dens_max);
				//printf("palette defined (%le 'blue', 0 'white', %le 'red')\n", dens_min, dens_max);
			}
		}
		else {
			fprintf(gp,"set palette defined (0 'white', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %le 'red')\n", scale/5., scale*2/5., scale*3/5.,scale*4/5, scale);
			fprintf(gp, " set cbrange [0:%le]\n", dens_max);
			//printf("palette defined (0 'white', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %le 'red')\n", scale/5., scale*2/5., scale*3/5.,scale*4/5, scale);
		}
	}
	else {
		fprintf(gp,"set palette defined (%le 'red', %le 'orange', %le 'yellow', %le 'green', %le 'cyan', 0 'white')\n", scale, scale*4/5., scale*3/5.,scale*2/5, scale/5);
		fprintf(gp, " set cbrange [%le:0]\n", dens_min);
		//printf("palette defined (%le 'red', %le 'orange', %le 'yellow', %le 'green', %le 'cyan', 0 'white')\n", scale, scale*4/5., scale*3/5.,scale*2/5, scale/5);
	}
	
	//fprintf(gp,"set palette defined (0 'white', %le 'blue', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %le 'red')\n", scale/6., scale*2/6., scale*3/6.,scale*4/6, scale*5/6, scale);
	//fprintf(gp,"set palette defined (0 'white', 5 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n"
	//if(dens_max>100) fprintf(gp,"set palette defined (0 'white', 5 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n", dens_max/5, dens_max*2/5, dens_max*3/5, dens_max*4/5, dens_max);
	//else if(dens_max>10) fprintf(gp,"set palette defined (0 'white', %i 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n", (dens_max/25+1), dens_max/5, dens_max*2/5, dens_max*3/5, dens_max*4/5, dens_max);
	//else if(dens_max>5) fprintf(gp,"set palette defined (0 'white', 1 'blue', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %i 'red')\n", (double) (dens_max/5.0), (double) (dens_max*2.0/5.0), (double) (dens_max*3.0/5.0), (double) (dens_max*4.0/5.0), dens_max);
	//else fprintf(gp,"set palette defined (0 'white', 1 'blue', 2 'green', 3 'yellow', 4 'orange', 5 'red')\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set cblabel '%s' offset 3,0 font 'helvetica,25.'\n",cblabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	//if(dens_max>5) fprintf(gp, " set cbrange [0:%i]\n", dens_max);
	//else fprintf(gp, " set cbrange [0:5]\n");
	//fprintf(gp, " set cbrange [%le:%le]\n", dens_min, dens_max);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		//fprintf(gp," set terminal postscript color\n");
		fprintf(gp," set terminal postscript eps enhanced color \n");
		fprintf(gp," set output '%s' \n", psfilename);
		fprintf(gp,"splot '%s'\n", txtfilename);
		fprintf(gp," set output\n");
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"splot '%s'\n", txtfilename);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot3dmap_gen(char *txtfilename, char *xcolumn, char *ycolumn, char *ccolumn, char *xlabel, char *ylabel, char *cblabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_min, double dens_max)
{
	double scale;
	FILE *gp;
	scale = dens_max-dens_min;
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,20.'\n");
	//fprintf(gp,"set bmargin 5\n");
	//fprintf(gp,"set rmargin 10\n");
	fprintf(gp,"set pm3d map\n");
	//fprintf(gp,"set palette defined (0 'white', 6 'yellow', 12 'orange', 18 'red')\n");//,dens_max/3., dens_max/2.,dens_max);
	//fprintf(gp,"set palette rgbformulae 33,13,10\n");
	//fprintf(gp,"set palette rgbformulae 22,13,-31\n");
	if(dens_max>0) {
		if(dens_min<0 && fabs(dens_min)>scale/100) {
			if(fabs(dens_min)<scale/5. || dens_max<scale/5.) {
				fprintf(gp,"set palette defined (%le 'blue', 0 'white', %le 'red')\n", -scale, scale);
				fprintf(gp, " set cbrange [%le:%le]\n", -scale, scale);
				//printf("palette defined (%le 'blue', 0 'white', %le 'red')\n", -scale, scale);
			}
			else {
				fprintf(gp,"set palette defined (%le 'blue', 0 'white', %le 'red')\n", dens_min, dens_max);
				fprintf(gp, " set cbrange [%le:%le]\n", dens_min, dens_max);
				//printf("palette defined (%le 'blue', 0 'white', %le 'red')\n", dens_min, dens_max);
			}
		}
		else {
			fprintf(gp,"set palette defined (0 'white', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %le 'red')\n", scale/5., scale*2/5., scale*3/5.,scale*4/5, scale);
			fprintf(gp, " set cbrange [0:%le]\n", dens_max);
			//printf("palette defined (0 'white', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %le 'red')\n", scale/5., scale*2/5., scale*3/5.,scale*4/5, scale);
		}
	}
	else {
		fprintf(gp,"set palette defined (%le 'red', %le 'orange', %le 'yellow', %le 'green', %le 'cyan', 0 'white')\n", scale, scale*4/5., scale*3/5.,scale*2/5, scale/5);
		fprintf(gp, " set cbrange [%le:0]\n", dens_min);
		//printf("palette defined (%le 'red', %le 'orange', %le 'yellow', %le 'green', %le 'cyan', 0 'white')\n", scale, scale*4/5., scale*3/5.,scale*2/5, scale/5);
	}
	
	//fprintf(gp,"set palette defined (0 'white', %le 'blue', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %le 'red')\n", scale/6., scale*2/6., scale*3/6.,scale*4/6, scale*5/6, scale);
	//fprintf(gp,"set palette defined (0 'white', 5 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n"
	//if(dens_max>100) fprintf(gp,"set palette defined (0 'white', 5 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n", dens_max/5, dens_max*2/5, dens_max*3/5, dens_max*4/5, dens_max);
	//else if(dens_max>10) fprintf(gp,"set palette defined (0 'white', %i 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n", (dens_max/25+1), dens_max/5, dens_max*2/5, dens_max*3/5, dens_max*4/5, dens_max);
	//else if(dens_max>5) fprintf(gp,"set palette defined (0 'white', 1 'blue', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %i 'red')\n", (double) (dens_max/5.0), (double) (dens_max*2.0/5.0), (double) (dens_max*3.0/5.0), (double) (dens_max*4.0/5.0), dens_max);
	//else fprintf(gp,"set palette defined (0 'white', 1 'blue', 2 'green', 3 'yellow', 4 'orange', 5 'red')\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set cblabel '%s' offset 3,0 font 'helvetica,25.'\n",cblabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	//if(dens_max>5) fprintf(gp, " set cbrange [0:%i]\n", dens_max);
	//else fprintf(gp, " set cbrange [0:5]\n");
	//fprintf(gp, " set cbrange [%le:%le]\n", dens_min, dens_max);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		//fprintf(gp," set terminal postscript color\n");
		fprintf(gp," set terminal postscript eps enhanced color \n");
		fprintf(gp," set output '%s' \n", psfilename);
		fprintf(gp,"splot '%s' using %s:%s:%s\n", txtfilename, xcolumn, ycolumn, ccolumn);
		fprintf(gp," set output\n");
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"splot '%s' using %s:%s:%s\n", txtfilename, xcolumn, ycolumn, ccolumn);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot3dmap_plusfile(char *txtfile_2d, char *xcolumn_2d, char *ycolumn_2d, char *with_2d, char *txtfilename, char *xlabel, char *ylabel, char *cblabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_min, double dens_max)
{
	double scale;
	FILE *gp;
	scale = dens_max-dens_min;
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,20.'\n");
	//fprintf(gp,"set bmargin 5\n");
	//fprintf(gp,"set rmargin 10\n");
	fprintf(gp,"set pm3d map\n");
	//fprintf(gp,"set palette defined (0 'white', 6 'yellow', 12 'orange', 18 'red')\n");//,dens_max/3., dens_max/2.,dens_max);
	//fprintf(gp,"set palette rgbformulae 33,13,10\n");
	//fprintf(gp,"set palette rgbformulae 22,13,-31\n");
	printf("dens_min=%le, dens_max=%le, scale=%le\n", dens_min, dens_max, scale);
	if(dens_max>0) {
		if(dens_min<0 && fabs(dens_min)>scale/100) {
			if(fabs(dens_min)<scale/5. || dens_max<scale/5.) {
				fprintf(gp,"set palette defined (%le 'blue', 0 'white', %le 'red')\n", -scale, scale);
				fprintf(gp, " set cbrange [%le:%le]\n", -scale, scale);
				//printf("palette defined (%le 'blue', 0 'white', %le 'red')\n", -scale, scale);
			}
			else {
				fprintf(gp,"set palette defined (%le 'blue', 0 'white', %le 'red')\n", dens_min, dens_max);
				fprintf(gp, " set cbrange [%le:%le]\n", dens_min, dens_max);
				//printf("palette defined (%le 'blue', 0 'white', %le 'red')\n", dens_min, dens_max);
			}
		}
		else {
			fprintf(gp,"set palette defined (0 'white', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %le 'red')\n", scale/5., scale*2/5., scale*3/5.,scale*4/5, scale);
			fprintf(gp, " set cbrange [0:%le]\n", dens_max);
			//printf("palette defined (0 'white', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %le 'red')\n", scale/5., scale*2/5., scale*3/5.,scale*4/5, scale);
		}
	}
	else {
		fprintf(gp,"set palette defined (%le 'red', %le 'orange', %le 'yellow', %le 'green', %le 'cyan', 0 'white')\n", scale, scale*4/5., scale*3/5.,scale*2/5, scale/5);
		fprintf(gp, " set cbrange [%le:0]\n", dens_min);
		//printf("palette defined (%le 'red', %le 'orange', %le 'yellow', %le 'green', %le 'cyan', 0 'white')\n", scale, scale*4/5., scale*3/5.,scale*2/5, scale/5);
	}
	
	//fprintf(gp,"set palette defined (0 'white', %le 'blue', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %le 'red')\n", scale/6., scale*2/6., scale*3/6.,scale*4/6, scale*5/6, scale);
	//fprintf(gp,"set palette defined (0 'white', 5 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n"
	//if(dens_max>100) fprintf(gp,"set palette defined (0 'white', 5 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n", dens_max/5, dens_max*2/5, dens_max*3/5, dens_max*4/5, dens_max);
	//else if(dens_max>10) fprintf(gp,"set palette defined (0 'white', %i 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n", (dens_max/25+1), dens_max/5, dens_max*2/5, dens_max*3/5, dens_max*4/5, dens_max);
	//else if(dens_max>5) fprintf(gp,"set palette defined (0 'white', 1 'blue', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %i 'red')\n", (double) (dens_max/5.0), (double) (dens_max*2.0/5.0), (double) (dens_max*3.0/5.0), (double) (dens_max*4.0/5.0), dens_max);
	//else fprintf(gp,"set palette defined (0 'white', 1 'blue', 2 'green', 3 'yellow', 4 'orange', 5 'red')\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set cblabel '%s' offset 3,0 font 'helvetica,25.'\n",cblabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	//if(dens_max>5) fprintf(gp, " set cbrange [0:%i]\n", dens_max);
	//else fprintf(gp, " set cbrange [0:5]\n");
	//fprintf(gp, " set cbrange [%le:%le]\n", dens_min, dens_max);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		//fprintf(gp," set terminal postscript color\n");
		fprintf(gp," set terminal postscript eps enhanced color \n");
		fprintf(gp," set output '%s' \n", psfilename);
		fprintf(gp,"set parametric\n");
		fprintf(gp,"splot '%s' using 1:2:3 with pm3d, '%s' using %s:%s:(0.0) with %s\n", txtfilename, txtfile_2d, xcolumn_2d, ycolumn_2d, with_2d);
		
		//fprintf(gp, "unset object\n");
		//fprintf(gp,"plot '%s' using %s:%s with %s\n", txtfilename1, xcolumn, ycolumn, with);
		fprintf(gp," set output\n");
		//fprintf(gp," set terminal x11\n");
	}
	fprintf(gp," set terminal x11\n");
	fprintf(gp,"splot '%s' using 1:2:3 with pm3d, '%s' using %s:%s:(0.0) with %s\n", txtfilename, txtfile_2d, xcolumn_2d, ycolumn_2d, with_2d);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot3dmap_plus2files(char *txtfile_2d, char *xcolumn_2d, char *ycolumn_2d, char *with_2d, char *txtfile_2d2, char *xcolumn_2d2, char *ycolumn_2d2, char *with_2d2, char *txtfilename, char *xlabel, char *ylabel, char *cblabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_min, double dens_max)
{
	double scale;
	FILE *gp;
	scale = dens_max-dens_min;
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,20.'\n");
	//fprintf(gp,"set bmargin 5\n");
	//fprintf(gp,"set rmargin 10\n");
	fprintf(gp,"set pm3d map\n");
	//fprintf(gp,"set palette defined (0 'white', 6 'yellow', 12 'orange', 18 'red')\n");//,dens_max/3., dens_max/2.,dens_max);
	//fprintf(gp,"set palette rgbformulae 33,13,10\n");
	//fprintf(gp,"set palette rgbformulae 22,13,-31\n");
	printf("dens_min=%le, dens_max=%le, scale=%le\n", dens_min, dens_max, scale);
	if(dens_max>0) {
		if(dens_min<0 && fabs(dens_min)>scale/100) {
			if(fabs(dens_min)<scale/5. || dens_max<scale/5.) {
				fprintf(gp,"set palette defined (%le 'blue', 0 'white', %le 'red')\n", -scale, scale);
				fprintf(gp, " set cbrange [%le:%le]\n", -scale, scale);
				//printf("palette defined (%le 'blue', 0 'white', %le 'red')\n", -scale, scale);
			}
			else {
				fprintf(gp,"set palette defined (%le 'blue', 0 'white', %le 'red')\n", dens_min, dens_max);
				fprintf(gp, " set cbrange [%le:%le]\n", dens_min, dens_max);
				//printf("palette defined (%le 'blue', 0 'white', %le 'red')\n", dens_min, dens_max);
			}
		}
		else {
			fprintf(gp,"set palette defined (0 'white', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %le 'red')\n", scale/5., scale*2/5., scale*3/5.,scale*4/5, scale);
			fprintf(gp, " set cbrange [0:%le]\n", dens_max);
			//printf("palette defined (0 'white', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %le 'red')\n", scale/5., scale*2/5., scale*3/5.,scale*4/5, scale);
		}
	}
	else {
		fprintf(gp,"set palette defined (%le 'red', %le 'orange', %le 'yellow', %le 'green', %le 'cyan', 0 'white')\n", scale, scale*4/5., scale*3/5.,scale*2/5, scale/5);
		fprintf(gp, " set cbrange [%le:0]\n", dens_min);
		//printf("palette defined (%le 'red', %le 'orange', %le 'yellow', %le 'green', %le 'cyan', 0 'white')\n", scale, scale*4/5., scale*3/5.,scale*2/5, scale/5);
	}
	
	//fprintf(gp,"set palette defined (0 'white', %le 'blue', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %le 'red')\n", scale/6., scale*2/6., scale*3/6.,scale*4/6, scale*5/6, scale);
	//fprintf(gp,"set palette defined (0 'white', 5 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n"
	//if(dens_max>100) fprintf(gp,"set palette defined (0 'white', 5 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n", dens_max/5, dens_max*2/5, dens_max*3/5, dens_max*4/5, dens_max);
	//else if(dens_max>10) fprintf(gp,"set palette defined (0 'white', %i 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n", (dens_max/25+1), dens_max/5, dens_max*2/5, dens_max*3/5, dens_max*4/5, dens_max);
	//else if(dens_max>5) fprintf(gp,"set palette defined (0 'white', 1 'blue', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %i 'red')\n", (double) (dens_max/5.0), (double) (dens_max*2.0/5.0), (double) (dens_max*3.0/5.0), (double) (dens_max*4.0/5.0), dens_max);
	//else fprintf(gp,"set palette defined (0 'white', 1 'blue', 2 'green', 3 'yellow', 4 'orange', 5 'red')\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set cblabel '%s' offset 3,0 font 'helvetica,25.'\n",cblabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	//if(dens_max>5) fprintf(gp, " set cbrange [0:%i]\n", dens_max);
	//else fprintf(gp, " set cbrange [0:5]\n");
	//fprintf(gp, " set cbrange [%le:%le]\n", dens_min, dens_max);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		//fprintf(gp," set terminal postscript color\n");
		fprintf(gp," set terminal postscript eps enhanced color \n");
		fprintf(gp," set output '%s' \n", psfilename);
		fprintf(gp,"set parametric\n");
		fprintf(gp,"splot '%s' using 1:2:3 with pm3d, '%s' using %s:%s:(0.0) with %s, '%s' using %s:%s:(0.0) with %s\n", txtfilename, txtfile_2d, xcolumn_2d, ycolumn_2d, with_2d, txtfile_2d2, xcolumn_2d2, ycolumn_2d2, with_2d2);
		
		//fprintf(gp, "unset object\n");
		//fprintf(gp,"plot '%s' using %s:%s with %s\n", txtfilename1, xcolumn, ycolumn, with);
		fprintf(gp," set output\n");
		//fprintf(gp," set terminal x11\n");
	}
	//fprintf(gp,"splot '%s'\n", txtfilename);
	fprintf(gp," set terminal x11\n");
	fprintf(gp,"splot '%s' using 1:2:3 with pm3d, '%s' using %s:%s:(0.0) with %s, '%s' using %s:%s:(0.0) with %s\n", txtfilename, txtfile_2d, xcolumn_2d, ycolumn_2d, with_2d, txtfile_2d2, xcolumn_2d2, ycolumn_2d2, with_2d2);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot3dmap_2d(char *filemap, char *file2d, char *para2d, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_min, double dens_max)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,20.'\n");
	//fprintf(gp,"set bmargin 5\n");
	//fprintf(gp,"set rmargin 10\n");
	fprintf(gp,"set pm3d map\n");
	//fprintf(gp,"set palette defined (0 'white', 6 'yellow', 12 'orange', 18 'red')\n");//,dens_max/3., dens_max/2.,dens_max);
	//fprintf(gp,"set palette rgbformulae 33,13,10\n");
	//fprintf(gp,"set palette rgbformulae 22,13,-31\n");
	fprintf(gp,"set palette defined (%le 'blue', 0 'white', %le 'red')\n", dens_min, dens_max);
	//if(dens_max>100) fprintf(gp,"set palette defined (0 'white', 5 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n", dens_max/5, dens_max*2/5, dens_max*3/5, dens_max*4/5, dens_max);
	//else if(dens_max>10) fprintf(gp,"set palette defined (0 'white', %i 'blue', %i 'cyan', %i 'green', %i 'yellow', %i 'orange', %i 'red')\n", (dens_max/25+1), dens_max/5, dens_max*2/5, dens_max*3/5, dens_max*4/5, dens_max);
	//else if(dens_max>5) fprintf(gp,"set palette defined (0 'white', 1 'blue', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %i 'red')\n", (double) (dens_max/5.0), (double) (dens_max*2.0/5.0), (double) (dens_max*3.0/5.0), (double) (dens_max*4.0/5.0), dens_max);
	//else fprintf(gp,"set palette defined (0 'white', 1 'blue', 2 'green', 3 'yellow', 4 'orange', 5 'red')\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	//if(dens_max>5) fprintf(gp, " set cbrange [0:%i]\n", dens_max);
	//else fprintf(gp, " set cbrange [0:5]\n");
	//fprintf(gp, " set cbrange [%le:%le]\n", dens_min, dens_max);
	fprintf(gp, " set cbrange [%le:%le]\n", dens_min, dens_max);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		//fprintf(gp," set terminal postscript color\n");
		fprintf(gp," set terminal postscript eps enhanced color \n");
		fprintf(gp," set output '%s' \n", psfilename);
		fprintf(gp,"splot '%s'\n", filemap);
		fprintf(gp, "unset object\n");
		fprintf(gp,"plot '%s' %s\n", file2d, para2d);
		fprintf(gp," set output\n");
		//fprintf(gp," set terminal x11\n");
	}
	//fprintf(gp,"splot '%s'\n", filemap);
	//fprintf(gp,"plot '%s'\n", file2d);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

/*extern void easyplot3d(char *txtfilename, char *with, char *xlabel, char *ylabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'Times-Roman,20.'\n");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 10\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(zrange != NULL) fprintf(gp," set zrange %s \n", zrange);
	
	fprintf(gp, "set dgrid3d 50,50\n");
	fprintf(gp, "set hidden3d\n");
	fprintf(gp, "set view 60, 60\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set terminal postscript eps enhanced color \n");
	fprintf(gp," set output '%s' \n", psfilename);
	fprintf(gp,"splot '%s' with %s\n", txtfilename, with);
	fprintf(gp," set output\n");
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}//*/

extern void easyplot3d(char *txtfilename, char *xcolumn, char *ycolumn, char *zcolumn, char *with, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set zlabel '%s' offset -3,0 font 'helvetica,25.'\n", zlabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(zrange != NULL) fprintf(gp," set zrange %s \n", zrange);
	
	//fprintf(gp, "set dgrid3d 50,50\n");
	//fprintf(gp, "set hidden3d\n");
	//fprintf(gp, "set view 60, 60\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set terminal postscript eps enhanced color \n");
	fprintf(gp," set output '%s' \n", psfilename);
	fprintf(gp,"splot '%s' using %s:%s:%s with %s\n", txtfilename, xcolumn, ycolumn, zcolumn, with);
	//fprintf(gp," set output\n");
	fprintf(gp," set terminal x11\n");
	fprintf(gp,"splot '%s' using %s:%s:%s with %s\n", txtfilename, xcolumn, ycolumn, zcolumn, with);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot2_3d(char *txtfilename1, char *txtfilename2, char *xcolumn1, char *ycolumn1, char *zcolumn1, char *xcolumn2, char *ycolumn2, char *zcolumn2, char *with1, char *with2, char *title1, char *title2, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set zlabel '%s' offset -3,0 font 'helvetica,25.'\n", zlabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(zrange != NULL) fprintf(gp," set zrange %s \n", zrange);
	if(title1==NULL && title2==NULL) fprintf(gp,"unset key\n");
	
	//fprintf(gp, "set dgrid3d 50,50\n");
	//fprintf(gp, "set hidden3d\n");
	//fprintf(gp, "set view 60, 60\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set terminal postscript eps enhanced color \n");
	fprintf(gp," set output '%s' \n", psfilename);
	if(title1==NULL) fprintf(gp,"splot '%s' using %s:%s:%s with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1);
	else fprintf(gp,"splot '%s' using %s:%s:%s with %s title '%s', ", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, title1);
	if(title2==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle\n", txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s'\n", txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2, title2);
	//fprintf(gp," set output\n");
	fprintf(gp," set terminal x11\n");
	fprintf(gp,"splot '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot3_3d(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *xcolumn1, char *ycolumn1, char *zcolumn1, char *xcolumn2, char *ycolumn2, char *zcolumn2, char *xcolumn3, char *ycolumn3, char *zcolumn3, 
	char *with1, char *with2, char *with3, char *title1, char *title2, char *title3, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set zlabel '%s' offset -3,0 font 'helvetica,25.'\n", zlabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(zrange != NULL) fprintf(gp," set zrange %s \n", zrange);
	if(title1==NULL && title2==NULL && title3==NULL) fprintf(gp,"unset key\n");
	
	//fprintf(gp, "set dgrid3d 50,50\n");
	//fprintf(gp, "set hidden3d\n");
	//fprintf(gp, "set view 60, 60\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set terminal postscript eps enhanced color \n");
	fprintf(gp," set output '%s' \n", psfilename);
	if(title1==NULL) fprintf(gp,"splot '%s' using %s:%s:%s with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1);
	else fprintf(gp,"splot '%s' using %s:%s:%s with %s title '%s', ", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, title1);
	if(title2==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2, title2);
	if(title3==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle\n", txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s'\n", txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3, title3);
	//fprintf(gp," set output\n");
	fprintf(gp," set terminal x11\n");
	fprintf(gp,"splot '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s\n", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2, txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot4_3d(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *xcolumn1, char *ycolumn1, char *zcolumn1, char *xcolumn2, char *ycolumn2, char *zcolumn2, char *xcolumn3, char *ycolumn3, char *zcolumn3, char *xcolumn4, char *ycolumn4, char *zcolumn4, 
	char *with1, char *with2, char *with3, char *with4, char *title1, char *title2, char *title3, char *title4, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set zlabel '%s' offset -3,0 font 'helvetica,25.'\n", zlabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(zrange != NULL) fprintf(gp," set zrange %s \n", zrange);
	if(title1==NULL && title2==NULL && title3==NULL && title4==NULL) fprintf(gp,"unset key\n");
	
	
	//fprintf(gp, "set dgrid3d 50,50\n");
	//fprintf(gp, "set hidden3d\n");
	//fprintf(gp, "set view 60, 60\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set terminal postscript eps enhanced color \n");
	fprintf(gp," set output '%s' \n", psfilename);
	if(title1==NULL) fprintf(gp,"splot '%s' using %s:%s:%s with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1);
	else fprintf(gp,"splot '%s' using %s:%s:%s with %s title '%s', ", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, title1);
	if(title2==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2, title2);
	if(title3==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3, title3);
	if(title4==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle\n", txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s'\n", txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4, title4);
	//fprintf(gp," set output\n");
	fprintf(gp," set terminal x11\n");
	fprintf(gp,"splot '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s\n", 
	txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2, txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3, txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot5_3d(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *txtfilename5, char *xcolumn1, char *ycolumn1, char *zcolumn1, char *xcolumn2, char *ycolumn2, char *zcolumn2, char *xcolumn3, char *ycolumn3, char *zcolumn3, char *xcolumn4, char *ycolumn4, char *zcolumn4, char *xcolumn5, char *ycolumn5, char *zcolumn5, 
	char *with1, char *with2, char *with3, char *with4, char *with5, char *title1, char *title2, char *title3, char *title4, char *title5, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set zlabel '%s' offset -3,0 font 'helvetica,25.'\n", zlabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(zrange != NULL) fprintf(gp," set zrange %s \n", zrange);
	if(title1==NULL && title2==NULL && title3==NULL && title4==NULL && title5==NULL) fprintf(gp,"unset key\n");
	
	//fprintf(gp, "set dgrid3d 50,50\n");
	//fprintf(gp, "set hidden3d\n");
	//fprintf(gp, "set view 60, 60\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set terminal postscript eps enhanced color \n");
	fprintf(gp," set output '%s' \n", psfilename);
	if(title1==NULL) fprintf(gp,"splot '%s' using %s:%s:%s with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1);
	else fprintf(gp,"splot '%s' using %s:%s:%s with %s title '%s', ", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, title1);
	if(title2==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2, title2);
	if(title3==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3, title3);
	if(title4==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4, title4);
	if(title5==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle\n", txtfilename5, xcolumn5, ycolumn5, zcolumn5, with5);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s'\n", txtfilename5, xcolumn5, ycolumn5, zcolumn5, with5, title5);
	//fprintf(gp," set output\n");
	fprintf(gp," set terminal x11\n");
	fprintf(gp,"splot '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s\n", 
	txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2, txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3, txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4, txtfilename5, xcolumn5, ycolumn5, zcolumn5, with5);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot6_3d(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *txtfilename5, char *txtfilename6, char *xcolumn1, char *ycolumn1, char *zcolumn1, char *xcolumn2, char *ycolumn2, char *zcolumn2, char *xcolumn3, char *ycolumn3, char *zcolumn3, char *xcolumn4, char *ycolumn4, char *zcolumn4, char *xcolumn5, char *ycolumn5, char *zcolumn5, char *xcolumn6, char *ycolumn6, char *zcolumn6, 
	char *with1, char *with2, char *with3, char *with4, char *with5, char *with6, char *title1, char *title2, char *title3, char *title4, char *title5, char *title6, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set lmargin 15\n");
	fprintf(gp,"set bmargin 25\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xtics offset 0,-1\n");
	fprintf(gp," set xlabel '%s' offset 0,-2 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset 5,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set zlabel '%s' offset -6,0 font 'helvetica,25.'\n", zlabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(zrange != NULL) fprintf(gp," set zrange %s \n", zrange);
	if(title1==NULL && title2==NULL && title3==NULL && title4==NULL && title5==NULL && title6==NULL) fprintf(gp,"unset key\n");
	
	//fprintf(gp, "set dgrid3d 50,50\n");
	//fprintf(gp, "set hidden3d\n");
	//fprintf(gp, "set view 60, 60\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set terminal postscript eps enhanced color \n");
	fprintf(gp," set output '%s' \n", psfilename);
	if(title1==NULL) fprintf(gp,"splot '%s' using %s:%s:%s with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1);
	else fprintf(gp,"splot '%s' using %s:%s:%s with %s title '%s', ", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, title1);
	if(title2==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2, title2);
	if(title3==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3, title3);
	if(title4==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4, title4);
	if(title5==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename5, xcolumn5, ycolumn5, zcolumn5, with5);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename5, xcolumn5, ycolumn5, zcolumn5, with5, title5);
	if(title6==NULL) fprintf(gp,"'%s' using %s:%s:%s with %s notitle\n", txtfilename6, xcolumn6, ycolumn6, zcolumn6, with6);
	else fprintf(gp,"'%s' using %s:%s:%s with %s title '%s'\n", txtfilename6, xcolumn6, ycolumn6, zcolumn6, with6, title6);
	//fprintf(gp," set output\n");
	fprintf(gp," set terminal x11\n");
	fprintf(gp,"splot '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s\n", 
	txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2, txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3, txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4, txtfilename5, xcolumn5, ycolumn5, zcolumn5, with5,  txtfilename6, xcolumn6, ycolumn6, zcolumn6, with6);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot7_3d(double x0, double z0, char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *txtfilename5, char *txtfilename6, char *txtfilename7, char *xcolumn1, char *ycolumn1, char *zcolumn1, char *xcolumn2, char *ycolumn2, char *zcolumn2, char *xcolumn3, char *ycolumn3, char *zcolumn3, char *xcolumn4, char *ycolumn4, char *zcolumn4, char *xcolumn5, char *ycolumn5, char *zcolumn5, char *xcolumn6, char *ycolumn6, char *zcolumn6, char *xcolumn7, char *ycolumn7, char *zcolumn7, 
	char *with1, char *with2, char *with3, char *with4, char *with5, char *with6, char *with7, char *title1, char *title2, char *title3, char *title4, char *title5, char *title6, char *title7, char *xlabel, char *ylabel, char *zlabel, char *xrange, char *yrange, char *zrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	gp=popen("gnuplot -persist","w");
	//fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp," set tics font 'helvetica,20.'\n");
	//fprintf(gp," set tics font 'helvetica,40.'\n");
	fprintf(gp,"set lmargin 25\n");
	fprintf(gp,"set bmargin 25\n");
	fprintf(gp,"set rmargin 25\n");
	//fprintf(gp," set xlabel '%s' offset 0,-3 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,20.' \n", xlabel);
	//fprintf(gp," set ylabel '%s' offset -1,-2 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set ylabel '%s' offset 0,-1 font 'helvetica,20.'\n", ylabel);
	//fprintf(gp," set zlabel '%s' offset -7,0 font 'helvetica,25.'\n", zlabel);
	fprintf(gp," set zlabel '%s' offset -4,0 font 'helvetica,20.'\n", zlabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(zrange != NULL) fprintf(gp," set zrange %s \n", zrange);
	if(title1==NULL && title2==NULL && title3==NULL && title4==NULL && title5==NULL && title6==NULL  && title7==NULL) fprintf(gp,"unset key\n");
	else fprintf(gp, "set key font 'helvetica,20.'\n");
	//fprintf(gp, "set dgrid3d 50,50\n");
	//fprintf(gp, "set hidden3d\n");
	//fprintf(gp, "set view 60, 60\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set terminal postscript eps enhanced color \n");
	fprintf(gp," set output '%s' \n", psfilename);
	if(title1==NULL) {
		fprintf(gp,"splot '%s' using %s:%s:%s with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1);
		fprintf(gp,"'%s' using %s:%s:(%le) with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, z0, with1);
		fprintf(gp,"'%s' using (%le):%s:%s with %s notitle, ", txtfilename1, x0, ycolumn1, zcolumn1, with1);
	}
	else {
		fprintf(gp,"splot '%s' using %s:%s:%s with %s title '%s', ", txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, title1);
		fprintf(gp,"'%s' using %s:%s:(%le) with %s notitle, ", txtfilename1, xcolumn1, ycolumn1, z0, with1);
		fprintf(gp,"'%s' using (%le):%s:%s with %s notitle, ", txtfilename1, x0, ycolumn1, zcolumn1, with1);
	}
	if(title2==NULL) {
		fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2);
		fprintf(gp,"'%s' using %s:%s:(%le) with %s notitle, ", txtfilename2, xcolumn2, ycolumn2, z0, with2);
		fprintf(gp,"'%s' using (%le):%s:%s with %s notitle, ", txtfilename2, x0, ycolumn2, zcolumn2, with2);
	}
	else {
		fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2, title2);
		//fprintf(gp,"'%s' using %s:%s:(%le) with %s notitle, ", txtfilename2, xcolumn2, ycolumn2, z0, with2);
		fprintf(gp,"'%s' using (%le):%s:%s with %s notitle, ", txtfilename2, x0, ycolumn2, zcolumn2, with2);
	}
	if(title3==NULL) {
		fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3);
		fprintf(gp,"'%s' using %s:%s:(%le) with %s notitle, ", txtfilename3, xcolumn3, ycolumn3, z0, with3);
		fprintf(gp,"'%s' using (%le):%s:%s with %s notitle, ", txtfilename3, x0, ycolumn3, zcolumn3, with3);
	}
	else {
		fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3, title3);
		//fprintf(gp,"'%s' using %s:%s:(%le) with %s notitle, ", txtfilename3, xcolumn3, ycolumn3, z0, with3);
		fprintf(gp,"'%s' using (%le):%s:%s with %s notitle, ", txtfilename3, x0, ycolumn3, zcolumn3, with3);
	}
	if(title4==NULL) {
		fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4);
		fprintf(gp,"'%s' using %s:%s:(%le) with %s notitle, ", txtfilename4, xcolumn4, ycolumn4, z0, with4);
		fprintf(gp,"'%s' using (%le):%s:%s with %s notitle, ", txtfilename4, x0, ycolumn4, zcolumn4, with4);
	}
	else {
		fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4, title4);
		//fprintf(gp,"'%s' using %s:%s:(%le) with %s notitle, ", txtfilename4, xcolumn4, ycolumn4, z0, with4);
		fprintf(gp,"'%s' using (%le):%s:%s with %s notitle, ", txtfilename4, x0, ycolumn4, zcolumn4, with4);
	}
	if(title5==NULL) {
		fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename5, xcolumn5, ycolumn5, zcolumn5, with5);
		fprintf(gp,"'%s' using %s:%s:(%le) with %s notitle, ", txtfilename5, xcolumn5, ycolumn5, z0, with5);
		fprintf(gp,"'%s' using (%le):%s:%s with %s notitle, ", txtfilename5, x0, ycolumn5, zcolumn5, with5);
	}
	else {
		fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename5, xcolumn5, ycolumn5, zcolumn5, with5, title5);
		//fprintf(gp,"'%s' using %s:%s:(%le) with %s notitle, ", txtfilename5, xcolumn5, ycolumn5, z0, with5);
		fprintf(gp,"'%s' using (%le):%s:%s with %s notitle, ", txtfilename5, x0, ycolumn5, zcolumn5, with5);
	}
	if(title6==NULL) {
		fprintf(gp,"'%s' using %s:%s:%s with %s notitle, ", txtfilename6, xcolumn6, ycolumn6, zcolumn6, with6); //remove for pete's picture
		fprintf(gp,"'%s' using %s:%s:(%le) with %s notitle, ", txtfilename6, xcolumn6, ycolumn6, z0, with6);
		//fprintf(gp,"'%s' using (%le):%s:%s with %s notitle, ", txtfilename6, x0, ycolumn6, zcolumn6, with6);
	}
	else {
		fprintf(gp,"'%s' using %s:%s:%s with %s title '%s', ", txtfilename6, xcolumn6, ycolumn6, zcolumn6, with6, title6);  //remove for pete's picture
		fprintf(gp,"'%s' using %s:%s:(%le) with %s title '%s', ", txtfilename6, xcolumn6, ycolumn6, z0, with6, title6);
		//fprintf(gp,"'%s' using (%le):%s:%s with %s notitle, ", txtfilename6, x0, ycolumn6, zcolumn6, with6);
	}
	if(title7==NULL) {
		fprintf(gp,"'%s' using %s:%s:%s with %s notitle,", txtfilename7, xcolumn7, ycolumn7, zcolumn7, with7);  //remove for pete's picture
		fprintf(gp,"'%s' using %s:%s:(%le) with %s notitle\n", txtfilename7, xcolumn7, ycolumn7, z0, with7);
		//fprintf(gp,"'%s' using (%le):%s:%s with %s notitle\n", txtfilename7, x0, ycolumn7, zcolumn7, with7);
	}
	else {
	//	fprintf(gp,"'%s' using %s:%s:%s with %s title '%s',", txtfilename7, xcolumn7, ycolumn7, zcolumn7, with7, title7);  //remove for pete's picture
		fprintf(gp,"'%s' using %s:%s:(%le) with %s title '%s'\n", txtfilename7, xcolumn7, ycolumn7, z0, with7, title7);
		//fprintf(gp,"'%s' using (%le):%s:%s with %s notitle\n", txtfilename7, x0, ycolumn7, zcolumn7, with7);
	}
	//fprintf(gp," set output\n");
	fprintf(gp," set terminal x11\n");
	//fprintf(gp,"splot '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s, '%s' using %s:%s:%s with %s\n", 
	//txtfilename1, xcolumn1, ycolumn1, zcolumn1, with1, txtfilename2, xcolumn2, ycolumn2, zcolumn2, with2, txtfilename3, xcolumn3, ycolumn3, zcolumn3, with3, txtfilename4, xcolumn4, ycolumn4, zcolumn4, with4, txtfilename5, xcolumn5, ycolumn5, zcolumn5, with5,  txtfilename6, xcolumn6, ycolumn6, zcolumn6, with6);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot2_left_right(char *txtfilename1, char *txtfilename2, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *with1, char *with2, char *xlabel, char *ylabel, char *y2label, char *xrange, char *yrange, char *y2range, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 15\n");
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(y2range != NULL) fprintf(gp," set y2range %s \n", y2range);
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set y2label \"%s\" offset 3,0 font 'helvetica,25.'\n", y2label);
	fprintf(gp,"set ytics nomirror\n set mxtics 5\nset mytics 5\nset my2tics 5\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color\n");
		fprintf(gp," set output '%s' \n", psfilename);
		//fprintf(gp,"plot '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
		fprintf(gp,"plot '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2);
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2);
	//fprintf(gp,"plot '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
	//fprintf(gp,"plot 'data/efb_betafunc.dat' using ($1):($2*2) with lines lt 1 lc 0 lw 1 axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot3_left_right(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *with1, char *with2, char *with3, char *title1, char *title2, char *title3, char *xlabel, char *ylabel, char *y2label, char *xrange, char *yrange, char *y2range, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 15\n");
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(y2range != NULL) fprintf(gp," set y2range %s \n", y2range);
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set y2label \"%s\" offset 3,0 font 'helvetica,25.'\n", y2label);
	fprintf(gp,"set ytics nomirror\n set mxtics 5\nset mytics 5\nset my2tics 5\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp," set tics font 'helvetica,25.'\n");
	//fprintf(gp,"unset key\n");
	fprintf(gp,"set key left bottom\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color\n");
		fprintf(gp," set output '%s' \n", psfilename);
		//fprintf(gp,"plot '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
		if(title1==NULL) fprintf(gp,"plot '%s' using %s:%s with %s notitle axis x1y1, ", txtfilename1, xcolumn1, ycolumn1, with1);
		else fprintf(gp,"plot '%s' using %s:%s with %s title '%s' axis x1y1, ", txtfilename1, xcolumn1, ycolumn1, with1, title1);
		if(title2==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle axis x1y1, ", txtfilename2, xcolumn2, ycolumn2, with2);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s' axis x1y1, ", txtfilename2, xcolumn2, ycolumn2, with2, title2);
		if(title3==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle axis x1y2, ", txtfilename3, xcolumn3, ycolumn3, with3);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s' axis x1y2\n", txtfilename3, xcolumn3, ycolumn3, with3, title3);
		//fprintf(gp,"plot '%s' using %s:%s with %s title '%s' axis x1y1, '%s' using %s:%s with %s title '%s' axis x1y1, '%s' using %s:%s with %s title '%s' axis x1y2, '%s' using %s:%s with %s title '%s' axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, title1, txtfilename2, xcolumn2, ycolumn2, with2, title2, txtfilename3, xcolumn3, ycolumn3, with3, title3, txtfilename4, xcolumn4, ycolumn4, with4, title4);
		fprintf(gp," set terminal x11\n");
	}
	//fprintf(gp,"plot '%s' using %s:%s with %s title '%s' axis x1y1, '%s' using %s:%s with %s title '%s' axis x1y1, '%s' using %s:%s with %s title '%s' axis x1y2, '%s' using %s:%s with %s title '%s' axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, title1, txtfilename2, xcolumn2, ycolumn2, with2, title2, txtfilename3, xcolumn3, ycolumn3, with3, title3, txtfilename4, xcolumn4, ycolumn4, with4, title4);
	fprintf(gp,"plot '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", 
	txtfilename1, xcolumn1, ycolumn1, with1, 
	txtfilename2, xcolumn2, ycolumn2, with2, 
	txtfilename3, xcolumn3, ycolumn3, with3);
	//fprintf(gp,"plot 'data/efb_betafunc.dat' using ($1):($2*2) with lines lt 1 lc 0 lw 1 axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void easyplot4_left_right(char *txtfilename1, char *txtfilename2, char *txtfilename3, char *txtfilename4, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *xcolumn4, char *ycolumn4, char *with1, char *with2, char *with3, char *with4, char *title1, char *title2, char *title3, char *title4, char *xlabel, char *ylabel, char *y2label, char *xrange, char *yrange, char *y2range, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 15\n");
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(y2range != NULL) fprintf(gp," set y2range %s \n", y2range);
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set y2label \"%s\" offset 3,0 font 'helvetica,25.'\n", y2label);
	fprintf(gp,"set ytics nomirror\n set mxtics 5\nset mytics 5\nset my2tics 5\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp," set tics font 'helvetica,25.'\n");
	//fprintf(gp,"unset key\n");
	fprintf(gp,"set key left bottom\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color\n");
		fprintf(gp," set output '%s' \n", psfilename);
		//fprintf(gp,"plot '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
		if(title1==NULL) fprintf(gp,"plot '%s' using %s:%s with %s notitle axis x1y1, ", txtfilename1, xcolumn1, ycolumn1, with1);
		else fprintf(gp,"plot '%s' using %s:%s with %s title '%s' axis x1y1, ", txtfilename1, xcolumn1, ycolumn1, with1, title1);
		if(title2==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle axis x1y1, ", txtfilename2, xcolumn2, ycolumn2, with2);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s' axis x1y1, ", txtfilename2, xcolumn2, ycolumn2, with2, title2);
		if(title3==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle axis x1y2, ", txtfilename3, xcolumn3, ycolumn3, with3);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s' axis x1y2, ", txtfilename3, xcolumn3, ycolumn3, with3, title3);
		if(title4==NULL) fprintf(gp,"'%s' using %s:%s with %s notitle axis x1y2\n", txtfilename4, xcolumn4, ycolumn4, with4);
		else fprintf(gp,"'%s' using %s:%s with %s title '%s' axis x1y2\n", txtfilename4, xcolumn4, ycolumn4, with4, title4);
		
		//fprintf(gp,"plot '%s' using %s:%s with %s title '%s' axis x1y1, '%s' using %s:%s with %s title '%s' axis x1y1, '%s' using %s:%s with %s title '%s' axis x1y2, '%s' using %s:%s with %s title '%s' axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, title1, txtfilename2, xcolumn2, ycolumn2, with2, title2, txtfilename3, xcolumn3, ycolumn3, with3, title3, txtfilename4, xcolumn4, ycolumn4, with4, title4);
		fprintf(gp," set terminal x11\n");
	}
	//fprintf(gp,"plot '%s' using %s:%s with %s title '%s' axis x1y1, '%s' using %s:%s with %s title '%s' axis x1y1, '%s' using %s:%s with %s title '%s' axis x1y2, '%s' using %s:%s with %s title '%s' axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, title1, txtfilename2, xcolumn2, ycolumn2, with2, title2, txtfilename3, xcolumn3, ycolumn3, with3, title3, txtfilename4, xcolumn4, ycolumn4, with4, title4);
	fprintf(gp,"plot '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename2, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3, txtfilename4, xcolumn4, ycolumn4, with4);
	//fprintf(gp,"plot 'data/efb_betafunc.dat' using ($1):($2*2) with lines lt 1 lc 0 lw 1 axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

	// ************************************************************************************ //
	//							    "plot beta & disp functions								//
	// ************************************************************************************ //

extern void easyplot_beta_disp(char *txtfilename1, char *txtfilename3, char *txt_filename_efb, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *with1, char *with2, char *with3, char *xlabel, char *ylabel, char *y2label, char *xrange, char *yrange, char *y2range, char *psfilename, char *setoption, double efb_scale, int title)
{
	double bmin;
	FILE *gp;
	
	get_betamin_file(txtfilename1, &bmin);
	bmin = floor(bmin);
	printf("bmin = %lf\n", bmin);
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp,"set lmargin 15\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 15\n");
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(y2range != NULL) fprintf(gp," set y2range %s \n", y2range);
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set y2label \"%s\" offset 3,0 font 'helvetica,25.'\n", y2label);
	fprintf(gp,"set ytics nomirror\n set mxtics 5\nset mytics 5\nset my2tics 5\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp," set tics offset 0.3,0 font 'helvetica,25.'\n"); 
	//fprintf(gp,"unset key\n");
	if(title == NO) fprintf(gp,"unset key\n");
	else if(title == TL) fprintf(gp, "set key top left\n");
	else if(title == BL) fprintf(gp, "set key bottom left\n");
	else if(title == BR) fprintf(gp, "set key bottom right\n");
	else if(title == TR) fprintf(gp, "set key top right\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color\n");
		fprintf(gp," set output '%s' \n", psfilename);
		//fprintf(gp,"plot '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
		if(title==NO)fprintf(gp,"plot '%s' using 1:($2*%lf) with lines lt 1 lc 0 lw 1 axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txt_filename_efb, efb_scale, txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
		//else fprintf(gp,"plot '%s' using 1:($2*%lf) with lines lt 1 lc 0 lw 1 axis x1y1 title '', '%s' using %s:%s with %s axis x1y1 title ' {/Symbol b}_x', '%s' using %s:%s with %s axis x1y1 title ' {/Symbol b}_z', '%s' using %s:%s with %s axis x1y2 title ' {/Symbol h}'\n", txt_filename_efb, efb_scale, txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
		else fprintf(gp,"plot '%s' using 1:($2*%lf+%lf) with lines lt 1 lc 0 lw 1 axis x1y1 title '', '%s' using %s:%s with %s axis x1y1 title ' {/Symbol b}_x', '%s' using %s:%s with %s axis x1y1 title ' {/Symbol b}_z', '%s' using %s:%s with %s axis x1y2 title ' B_z'\n", txt_filename_efb, efb_scale, bmin, txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using 1:($2*%lf+%lf) with lines lt 1 lc 0 lw 1 axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txt_filename_efb, efb_scale, bmin, txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
	//fprintf(gp,"plot '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
	//fprintf(gp,"plot 'data/efb_betafunc.dat' using ($1):($2*2) with lines lt 1 lc 0 lw 1 axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}
/*
extern void easyplot_beta_bz(char *txtfilename1, char *txtfilename3, char *txt_filename_efb, char *xcolumn1, char *ycolumn1, char *xcolumn2, char *ycolumn2, char *xcolumn3, char *ycolumn3, char *with1, char *with2, char *with3, char *xlabel, char *ylabel, char *y2label, char *xrange, char *yrange, char *y2range, char *psfilename, char *setoption, double efb_scale, int title)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp,"set lmargin 15\n");
	//fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 12\n");
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(y2range != NULL) fprintf(gp," set y2range %s \n", y2range);
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set y2label \"%s\" offset 3,0 font 'helvetica,25.'\n", y2label);
	fprintf(gp,"set ytics nomirror\n set mxtics 5\nset mytics 5\nset my2tics 5\n");
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp," set tics font 'helvetica,25.'\n"); 
	//fprintf(gp,"unset key\n");
	if(title == NO) fprintf(gp,"unset key\n");
	else if(title == TL) fprintf(gp, "set key top left\n");
	else if(title == BL) fprintf(gp, "set key bottom left\n");
	else if(title == BR) fprintf(gp, "set key bottom right\n");
	else if(title == TR) fprintf(gp, "set key top right\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color\n");
		fprintf(gp," set output '%s' \n", psfilename);
		//fprintf(gp,"plot '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
		if(title==NO)fprintf(gp,"plot '%s' using 1:($2*%lf) with lines lt 1 lc 0 lw 1 axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txt_filename_efb, efb_scale, txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
		else fprintf(gp,"plot '%s' using 1:($2*%lf) with lines lt 1 lc 0 lw 1 axis x1y1 title '', '%s' using %s:%s with %s axis x1y1 title ' {/Symbol b}_x', '%s' using %s:%s with %s axis x1y1 title ' {/Symbol b}_z', '%s' using %s:%s with %s axis x1y2 title 'Bz'\n", txt_filename_efb, efb_scale, txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using 1:($2*%lf) with lines lt 1 lc 0 lw 1 axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txt_filename_efb, efb_scale, txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
	//fprintf(gp,"plot '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
	//fprintf(gp,"plot 'data/efb_betafunc.dat' using ($1):($2*2) with lines lt 1 lc 0 lw 1 axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y1, '%s' using %s:%s with %s axis x1y2\n", txtfilename1, xcolumn1, ycolumn1, with1, txtfilename1, xcolumn2, ycolumn2, with2, txtfilename3, xcolumn3, ycolumn3, with3);
	pclose(gp);
}//*/
//old, not working with spi and several cells
/*extern void write_efb_txtfile(char *filename, struct Particle *reference, struct Lattice *latt)
{
	int i,j;
	double s, s0, efb_ymax, x, y, spi_shift, r_comp, r_cell, th_comp;
	struct Particle test_part, test_part2;
	struct Lattice tempo_latt;
	FILE *wfile;
	wfile = fopen(filename, "w");
	
	copy_latt(&tempo_latt, latt);
	test_part = *reference;
	efb_ymax = 1;
	
	s0 = 0.;
	for(i = 0; i < latt->nbcell; i++) {
		if(test_cell_map(&(latt->cell[i])) == YES ||
			strcmp(latt->cell[i].keyword, "drift") == 0 ||
		strcmp(latt->cell[i].keyword, "rf-thingap") == 0 ||
		strcmp(latt->cell[i].keyword, "collimator") == 0) continue;
		tempo_latt.cell[i].instrutype = CUP;
		for(j = 0; j < latt->cell[i].nbcomp; j++) {
			spi_shift = 0.0;
			r_cell = 0.0;
			if(latt->cell[i].boun.thmax!=0) {
				//track particle to the pole to get the radius:
				test_part2 = test_part; 
				if(strcmp(latt->cell[i].keyword, "ffag-spi-lin") == 0 || 
				strcmp(latt->cell[i].keyword, "ffag-spi-he") == 0 ||
				strcmp(latt->cell[i].keyword, "ffag-spi-enge") == 0) {
					th_comp = latt->cell[i].mpara[j][0] + (latt->cell[i].efben[j][0]+latt->cell[i].efbex[j][0])/2. + tan(latt->cell[i].mpara[j][4])*log(test_part.x/latt->cell[i].mpara[j][1]);
					tempo_latt.cell[i].instru.thmax = th_comp;
					part_cross_cell(&test_part2, &(tempo_latt.cell[i]),NULL);
					r_comp = sqrt(pow(test_part2.x,2) + pow(test_part2.y,2));
					spi_shift = tan(latt->cell[i].mpara[j][4])*log(r_comp/latt->cell[i].mpara[j][1]);
					printf("spi_shift = %le, x = %le\n",spi_shift,r_comp);
				}
				else if(strcmp(latt->cell[i].keyword, "ffag-tilt-lin") == 0 || 
				strcmp(latt->cell[i].keyword, "ffag-tilt-he") == 0 ||
				strcmp(latt->cell[i].keyword, "ffag-tilt-enge") == 0) {
					th_comp = latt->cell[i].mpara[j][0] + (latt->cell[i].efben[j][0]+latt->cell[i].efbex[j][0])/2. + latt->cell[i].mpara[j][5] - asin(latt->cell[i].mpara[j][1]*sin(latt->cell[i].mpara[j][5])/test_part.x);
					tempo_latt.cell[i].instru.thmax = th_comp;
					part_cross_cell(&test_part2, &(tempo_latt.cell[i]),NULL);
					r_comp = sqrt(pow(test_part2.x,2) + pow(test_part2.y,2));
					spi_shift = latt->cell[i].mpara[j][5] - asin(latt->cell[i].mpara[j][1]*sin(latt->cell[i].mpara[j][5])/r_comp);
				}
				else {
					th_comp = latt->cell[i].mpara[j][0] + (latt->cell[i].efben[j][0]+latt->cell[i].efbex[j][0])/2.;
					tempo_latt.cell[i].instru.thmax = th_comp;
					part_cross_cell(&test_part2, &(tempo_latt.cell[i]),NULL);
					r_comp = sqrt(pow(test_part2.x,2) + pow(test_part2.y,2));
				}
				
				r_cell += r_comp; //for computing s (below)
				
				x = s0 + (r_comp+latt->cell[i].deltar)*(latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] + spi_shift);
				y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				
				y = efb_ymax;
				fprintf(wfile, "%le  %le\n",x,y);
				
				x = s0 + (test_part.x+latt->cell[i].deltar)*(latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0] + spi_shift);
				fprintf(wfile, "%le  %le\n",x,y);
				
				y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				
				x = s0 + (r_comp+latt->cell[i].deltar)*(latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] + spi_shift);
				fprintf(wfile, "%le  %le\n",x,y);
				
				fprintf(wfile, "\n");
			}
			else if(latt->cell[i].boun.ymax!=0) {
				if(strcmp(latt->cell[i].keyword, "ffag-sti-lin") == 0 || 
				strcmp(latt->cell[i].keyword, "ffag-sti-he") == 0 ||
				strcmp(latt->cell[i].keyword, "ffag-sti-enge") == 0) {
					spi_shift = tan(latt->cell[i].mpara[j][4])*(test_part.x-latt->cell[i].mpara[j][1]);
				}
				x = s0 + (latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] + spi_shift);
				y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				
				y = efb_ymax;
				fprintf(wfile, "%le  %le\n",x,y);
				
				x = s0 + (latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0] + spi_shift);
				fprintf(wfile, "%le  %le\n",x,y);
				
				y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				
				x = s0 + (latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] + spi_shift);
				fprintf(wfile, "%le  %le\n",x,y);
				
				fprintf(wfile, "\n");
			}
		}
		if(latt->cell[i].boun.thmax!=0) {
			r_cell = r_cell/latt->cell[i].nbcomp;
			tempo_latt.cell[i].instru.thmax = latt->cell[i].boun.thmax;
			s = (r_cell+tempo_latt.cell[i].deltar)*latt->cell[i].boun.thmax;
		}
		else if(latt->cell[i].boun.ymax!=0) s = latt->cell[i].boun.ymax;
		else s=0;
		test_part = *reference;
		part_cross_latt(&test_part, &tempo_latt,NULL); //(only to calculate s from sqrt(pow(test_part.x,2) + pow(test_part.y,2)))
		tempo_latt.cell[i].instrutype = NO;
		s0 += s;
	}
	fclose(wfile);
	free(tempo_latt.cell);
}//*/

//works only with spiral, not tilt
/*extern void write_efb_txtfile(char *filename, struct Particle *reference, struct Lattice *latt)
{
	int i,j,m, flag_en, flag_ex;
	double efb_ymax, x, y, th_comp_en, th_comp_ex, th_comp_old, s_en, s_ex, spi_angle;
	struct Particle test_part, test_part2;
	struct Lattice tempo_latt;
	FILE *wfile;
	wfile = fopen(filename, "w");
	
	copy_latt(&tempo_latt, latt);
	test_part = *reference;
	efb_ymax = 1;
	
	for(i = 0; i < latt->nbcell; i++) {
		if(test_cell_map(&(latt->cell[i])) == YES ||
			strcmp(latt->cell[i].keyword, "drift") == 0 ||
		strcmp(latt->cell[i].keyword, "rf-thingap") == 0 ||
		strcmp(latt->cell[i].keyword, "collimator") == 0) continue;
		tempo_latt.cell[i].instrutype = CUP;
		for(j = 0; j < latt->cell[i].nbcomp; j++) {
			if(strcmp(latt->cell[i].keyword, "ffag-spi-lin") == 0 || 
			strcmp(latt->cell[i].keyword, "ffag-spi-he") == 0 ||
			strcmp(latt->cell[i].keyword, "ffag-spi-enge") == 0 ||
			strcmp(latt->cell[i].keyword, "ffag-sti-lin") == 0 || 
			strcmp(latt->cell[i].keyword, "ffag-sti-he") == 0 ||
			strcmp(latt->cell[i].keyword, "ffag-sti-enge") == 0) spi_angle = latt->cell[i].mpara[j][4];
			else if(strcmp(latt->cell[i].keyword, "ffag-tiltpol-lin") == 0 ||
				strcmp(latt->cell[i].keyword, "ffag-tiltpol-lin") == 0) spi_angle = latt->cell[i].mpara[j][2];
			else if(strcmp(latt->cell[i].keyword, "ffag-tiltpol-lin") == 0 ||
			strcmp(latt->cell[i].keyword, "ffag-tiltpol-he") == 0 ||
			strcmp(latt->cell[i].keyword, "ffag-tiltpol-enge") == 0) spi_angle = latt->cell[i].mpara[j][2];
			else spi_angle = 0.;
			if(latt->cell[i].boun.thmax!=0) {
				//track particle to the pole to get s:
				
				th_comp_en = latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] + tan(spi_angle)*log(test_part.x/latt->cell[i].mpara[j][1]);
				for(m=0;m<10;m++) {
					th_comp_old = th_comp_en;
					test_part2 = test_part;
					tempo_latt.cell[i].instru.thmax = th_comp_en;
					part_cross_cell(&test_part2, &(tempo_latt.cell[i]),NULL);
					th_comp_en = latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] + tan(spi_angle)*log(test_part2.x/latt->cell[i].mpara[j][1]);
					printf("m=%i, th_comp_old=%le, th_comp_en=%le\n", m, th_comp_old*180./PI, th_comp_en*180./PI);
					if(fabs(th_comp_en-th_comp_old)<1.e-4) break;
					if(m==9) errorstop("in write_efb_txtfile_spi, th_comp_en not found");
				}
				if(th_comp_en>0 && th_comp_en<latt->cell[i].boun.thmax) {
					flag_en = 1;
					s_en = test_part2.s;
				}
				else {
					flag_en = 0;
					s_en = test_part.s;
				}
				//printf("flag_en = %i, s_en = %le\n", flag_en,s_en);
				
				th_comp_ex = latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0] + tan(spi_angle)*log(test_part.x/latt->cell[i].mpara[j][1]);
				for(m=0;m<10;m++) {
					th_comp_old = th_comp_ex;
					test_part2 = test_part; 
					tempo_latt.cell[i].instru.thmax = th_comp_ex;
					part_cross_cell(&test_part2, &(tempo_latt.cell[i]),NULL);
					th_comp_ex = latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0] + tan(spi_angle)*log(test_part2.x/latt->cell[i].mpara[j][1]);
					//printf("m=%i, th_comp_old=%le, th_comp_ex=%le\n", m, th_comp_old*180./PI, th_comp_ex*180./PI);
					if(fabs(th_comp_ex-th_comp_old)<1.e-4) break;
					if(m==9) errorstop("in write_efb_txtfile_spi, th_comp_ex not found");
				}
				s_ex = test_part2.s;
				if(th_comp_ex>0 && th_comp_ex<latt->cell[i].boun.thmax) flag_ex = 1;
				else flag_ex = 0;
				//printf("flag_ex = %i, s_ex = %le\n", flag_ex,s_ex);
			}
			else if(latt->cell[i].boun.ymax!=0) {
				th_comp_en = latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] + tan(spi_angle)*(test_part.x-latt->cell[i].mpara[j][1]);
				for(m=0;m<5;m++) {
					th_comp_old = th_comp_en;
					test_part2 = test_part;
					tempo_latt.cell[i].instru.ymax = th_comp_en;
					part_cross_cell(&test_part2, &(tempo_latt.cell[i]),NULL);
					th_comp_en = latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] + tan(spi_angle)*(test_part.x-latt->cell[i].mpara[j][1]);
					//printf("m=%i, y_comp_old=%le, y_comp_en=%le\n", m, th_comp_old, th_comp_en);
					if(fabs(th_comp_en-th_comp_old)<1.e-4) break;
					if(m==4) errorstop("in write_efb_txtfile, y_comp_en not found");
				}
				if(th_comp_en>0 && th_comp_en<latt->cell[i].boun.ymax) {
					flag_en = 1;
					s_en = test_part2.s;
				}
				else {
					flag_en = 0;
					s_en = test_part.s;
				}
				//printf("flag_en = %i, s_en = %le\n", flag_en,s_en);
				th_comp_ex = latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0] + tan(spi_angle)*(test_part.x-latt->cell[i].mpara[j][1]);
				for(m=0;m<5;m++) {
					th_comp_old = th_comp_ex;
					test_part2 = test_part; 
					tempo_latt.cell[i].instru.thmax = th_comp_ex;
					part_cross_cell(&test_part2, &(tempo_latt.cell[i]),NULL);
					th_comp_ex = latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0] + tan(spi_angle)*(test_part.x-latt->cell[i].mpara[j][1]);
					//printf("m=%i, y_comp_old=%le, y_comp_ex=%le\n", m, th_comp_old, th_comp_ex);
					if(fabs(th_comp_ex-th_comp_old)<1.e-4) break;
					if(m==4) errorstop("in write_efb_txtfile, y_comp_ex not found");
				}
				s_ex = test_part2.s;
				if(th_comp_ex>0 && th_comp_ex<latt->cell[i].boun.ymax) flag_ex = 1;
				else flag_ex = 0;
				//printf("flag_ex = %i, s_ex = %le\n", flag_ex,s_ex);
			}
			if(flag_en+flag_ex>0) {
				x = s_en;
				y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				y = efb_ymax;
				fprintf(wfile, "%le  %le\n",x,y);
				x = s_ex;
				fprintf(wfile, "%le  %le\n",x,y);
				y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				x = s_en;
				fprintf(wfile, "%le  %le\n",x,y);
				fprintf(wfile, "\n");
			}
		}
		if(latt->cell[i].boun.thmax!=0) tempo_latt.cell[i].instru.thmax = latt->cell[i].boun.thmax;
		else if(latt->cell[i].boun.ymax!=0) tempo_latt.cell[i].instru.ymax = latt->cell[i].boun.ymax;
		else errorstop("strange in write_efb_txtfile");
		test_part = *reference;
		part_cross_latt(&test_part, &tempo_latt, NULL); //calculate s0
		tempo_latt.cell[i].instrutype = NO;
		//printf("s0 = %le\n",test_part.s);
	}
	fclose(wfile);
	free(tempo_latt.cell);
}//*/

extern void write_efb_txtfile(char *filename, struct Particle *reference, struct Lattice *latt)
{
	int i,j,m, flag_en, flag_ex, flag_spi=0, mmax=20;
	double efb_ymax, x, y, th_comp_en, th_comp_ex, th_comp_old, s_en, s_ex, spi_angle, newx, x_cell, y_cell, x_0, y_0;
	struct Particle test_part, test_part2;
	struct Lattice tempo_latt;
	FILE *wfile;
	wfile = fopen(filename, "w");
	
	copy_latt(&tempo_latt, latt);
	test_part = *reference;
	test_part.s = 0;
	efb_ymax = 1;
	
	for(i = 0; i < latt->nbcell; i++) {
		if(test_cell_map(&(latt->cell[i])) == YES ||
			strcmp(latt->cell[i].keyword, "drift") == 0 || strcmp(latt->cell[i].keyword, "rf-thingap") == 0 || strcmp(latt->cell[i].keyword, "collimator") == 0) continue;
		tempo_latt.cell[i].instrutype = CUP;
		for(j = 0; j < latt->cell[i].nbcomp; j++) { 
			//printf("comp %i\n", j);
			// flag_spi: 1=spiral, 3=straight tilted, 4=tilted sector, 5=tilted sector taylor, 6=rectangle vffa
			if(strcmp(latt->cell[i].keyword, "ffag-spi-lin") == 0 || 
			strcmp(latt->cell[i].keyword, "ffag-spi-he") == 0 ||
			strcmp(latt->cell[i].keyword, "ffag-spi-enge") == 0 ||
			strcmp(latt->cell[i].keyword, "ffag-spi-fullenge") == 0) flag_spi = 1; 
			else if(strcmp(latt->cell[i].keyword, "ffag-sti-lin") == 0 || 
			strcmp(latt->cell[i].keyword, "ffag-sti-he") == 0 ||
			strcmp(latt->cell[i].keyword, "ffag-sti-enge") == 0) flag_spi = 3;
			else if(strcmp(latt->cell[i].keyword, "ffag-tilt-lin") == 0 ||
			strcmp(latt->cell[i].keyword, "ffag-tilt-he") == 0 ||
			strcmp(latt->cell[i].keyword, "ffag-tilt-enge") == 0) flag_spi = 4;
			else if(strcmp(latt->cell[i].keyword, "ffag-tiltpol-lin") == 0 ||
			strcmp(latt->cell[i].keyword, "ffag-tiltpol-he") == 0 ||
			strcmp(latt->cell[i].keyword, "ffag-tiltpol-enge") == 0) flag_spi = 5;
			else if(strcmp(latt->cell[i].keyword, "vffa-rect-lin") == 0 || 
			strcmp(latt->cell[i].keyword, "vffa-rect-he") == 0 ||
			strcmp(latt->cell[i].keyword, "vffa-rect-enge") == 0 ||
			strcmp(latt->cell[i].keyword, "vffa-rect-enge-add") == 0 ||
			strcmp(latt->cell[i].keyword, "vffa-rect-atan") == 0) flag_spi = 6;
			if(flag_spi==5) spi_angle = latt->cell[i].mpara[j][2];
			else if(flag_spi==6) {
				x_0 = latt->cell[i].mpara[j][1]*cos(latt->cell[i].mpara[j][0]);
				y_0 = latt->cell[i].mpara[j][1]*sin(latt->cell[i].mpara[j][0]);
				spi_angle = latt->cell[i].mpara[j][0] + latt->cell[i].mpara[j][5];
			}
			else if(flag_spi>0) spi_angle = latt->cell[i].mpara[j][4];
			else spi_angle = 0.;
			if(latt->cell[i].boun.thmax!=0) {
				//track particle to the pole to get s:
				//if(flag_spi < 2) th_comp_en = latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] + tan(spi_angle)*log(sqrt(test_part.x*test_part.x+test_part.y*test_part.y)/latt->cell[i].mpara[j][1]); //spiral angle
				if(flag_spi < 2) th_comp_en = latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] + tan(spi_angle)*log(test_part.x/latt->cell[i].mpara[j][1]); //spiral angle
				else if(flag_spi == 6) th_comp_en = latt->cell[i].mpara[j][0] + atan_ratio(latt->cell[i].efben[j][0]*cos(latt->cell[i].mpara[j][5]),latt->cell[i].mpara[j][1]);
				else th_comp_en = latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] - asin(latt->cell[i].mpara[j][1]*sin(spi_angle)/test_part.x) + spi_angle; //tilted sector
				for(m=0;m<mmax;m++) {
					th_comp_old = th_comp_en;
					test_part2 = test_part;
					tempo_latt.cell[i].instru.thmax = th_comp_en;
					part_cross_cell(&test_part2, &(tempo_latt.cell[i]),NULL);
					
					if(flag_spi < 2) th_comp_en = latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] + tan(spi_angle)*log(test_part2.x/latt->cell[i].mpara[j][1]);
					else if(flag_spi == 6) {
						newx = (test_part2.y-y_0)*sin(spi_angle) + (test_part2.x-x_0)*cos(spi_angle);
						y_cell = y_0 + latt->cell[i].efben[j][0]*cos(spi_angle) + newx*sin(spi_angle);
						x_cell = x_0 - latt->cell[i].efben[j][0]*sin(spi_angle) + newx*cos(spi_angle);
						th_comp_en = atan_ratio(y_cell, x_cell);
					}
					else th_comp_en = latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] - asin(latt->cell[i].mpara[j][1]*sin(spi_angle)/test_part2.x) + spi_angle;
					//printf("m=%i, th_comp_old=%le, th_comp_en=%le\n", m, th_comp_old*180./PI, th_comp_en*180./PI);
					if(fabs(th_comp_en-th_comp_old)<1.e-4) break;
					if(m==mmax-1) errorstop("in write_efb_txtfile, th_comp_en not found");
				}
				if(th_comp_en>0 && th_comp_en<latt->cell[i].boun.thmax) {
					flag_en = 1;
					s_en = test_part2.s;
				}
				else {
					flag_en = 0;
					s_en = test_part.s;
				}
				//printf("flag_en = %i, s_en = %le\n", flag_en,s_en);
				
				if(flag_spi < 2) th_comp_ex = latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0] + tan(spi_angle)*log(test_part.x/latt->cell[i].mpara[j][1]);
				else if(flag_spi == 6) th_comp_en = latt->cell[i].mpara[j][0] + atan_ratio(latt->cell[i].efbex[j][0]*cos(latt->cell[i].mpara[j][5]),latt->cell[i].mpara[j][1]);
				else th_comp_ex = latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0] - asin(latt->cell[i].mpara[j][1]*sin(spi_angle)/test_part.x) + spi_angle;
				for(m=0;m<mmax;m++) {
					th_comp_old = th_comp_ex;
					test_part2 = test_part; 
					tempo_latt.cell[i].instru.thmax = th_comp_ex;
					part_cross_cell(&test_part2, &(tempo_latt.cell[i]),NULL);
					if(flag_spi < 2) th_comp_ex = latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0] + tan(spi_angle)*log(test_part2.x/latt->cell[i].mpara[j][1]);
					else if(flag_spi == 6) {
						newx = (test_part2.y-y_0)*sin(spi_angle) + (test_part2.x-x_0)*cos(spi_angle);
						y_cell = y_0 + latt->cell[i].efbex[j][0]*cos(spi_angle) + newx*sin(spi_angle);
						x_cell = x_0 - latt->cell[i].efbex[j][0]*sin(spi_angle) + newx*cos(spi_angle);
						th_comp_ex = atan_ratio(y_cell, x_cell);
					}
					else th_comp_ex = latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0] - asin(latt->cell[i].mpara[j][1]*sin(spi_angle)/test_part2.x) + spi_angle;
					//printf("m=%i, th_comp_old=%le, th_comp_ex=%le\n", m, th_comp_old*180./PI, th_comp_ex*180./PI);
					if(fabs(th_comp_ex-th_comp_old)<1.e-4) break;
					if(m==mmax-1) errorstop("in write_efb_txtfile, th_comp_ex not found");
				}
				s_ex = test_part2.s;
				if(th_comp_ex>0 && th_comp_ex<latt->cell[i].boun.thmax) flag_ex = 1;
				else flag_ex = 0;
				//printf("flag_ex = %i, s_ex = %le\n", flag_ex,s_ex);
			}
			else if(latt->cell[i].boun.ymax!=0) {
				th_comp_en = latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] + tan(spi_angle)*(test_part.x-latt->cell[i].mpara[j][1]);
				for(m=0;m<mmax;m++) {
					th_comp_old = th_comp_en;
					test_part2 = test_part;
					tempo_latt.cell[i].instru.ymax = th_comp_en;
					part_cross_cell(&test_part2, &(tempo_latt.cell[i]),NULL);
					th_comp_en = latt->cell[i].mpara[j][0] + latt->cell[i].efben[j][0] + tan(spi_angle)*(test_part.x-latt->cell[i].mpara[j][1]);
					//printf("m=%i, y_comp_old=%le, y_comp_en=%le\n", m, th_comp_old, th_comp_en);
					if(fabs(th_comp_en-th_comp_old)<1.e-4) break;
					if(m==mmax-1) errorstop("in write_efb_txtfile, y_comp_en not found");
				}
				if(th_comp_en>0 && th_comp_en<latt->cell[i].boun.ymax) {
					flag_en = 1;
					s_en = test_part2.s;
				}
				else {
					flag_en = 0;
					s_en = test_part.s;
				}
				//printf("flag_en = %i, s_en = %le\n", flag_en,s_en);
				th_comp_ex = latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0] + tan(spi_angle)*(test_part.x-latt->cell[i].mpara[j][1]);
				for(m=0;m<mmax;m++) {
					th_comp_old = th_comp_ex;
					test_part2 = test_part; 
					tempo_latt.cell[i].instru.ymax = th_comp_ex;
					part_cross_cell(&test_part2, &(tempo_latt.cell[i]),NULL);
					th_comp_ex = latt->cell[i].mpara[j][0] + latt->cell[i].efbex[j][0] + tan(spi_angle)*(test_part.x-latt->cell[i].mpara[j][1]);
					//printf("m=%i, y_comp_old=%le, y_comp_ex=%le\n", m, th_comp_old, th_comp_ex);
					if(fabs(th_comp_ex-th_comp_old)<1.e-4) break;
					if(m==mmax-1) errorstop("in write_efb_txtfile, y_comp_ex not found");
				}
				s_ex = test_part2.s;
				if(th_comp_ex>0 && th_comp_ex<latt->cell[i].boun.ymax) flag_ex = 1;
				else flag_ex = 0;
				//printf("flag_ex = %i, s_ex = %le\n", flag_ex,s_ex);
			}
			if(flag_en+flag_ex>0) {
				x = s_en;
				y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				y = efb_ymax;
				fprintf(wfile, "%le  %le\n",x,y);
				x = s_ex;
				fprintf(wfile, "%le  %le\n",x,y);
				y = 0;
				fprintf(wfile, "%le  %le\n",x,y);
				x = s_en;
				fprintf(wfile, "%le  %le\n",x,y);
				fprintf(wfile, "\n");
			}
		}
		if(latt->cell[i].boun.thmax!=0) tempo_latt.cell[i].instru.thmax = latt->cell[i].boun.thmax;
		else if(latt->cell[i].boun.ymax!=0) tempo_latt.cell[i].instru.ymax = latt->cell[i].boun.ymax;
		else errorstop("strange in write_efb_txtfile");
		test_part = *reference;
		part_cross_latt(&test_part, &tempo_latt, NULL); //calculate s0
		tempo_latt.cell[i].instrutype = NO;
		//printf("s0 = %le\n",test_part.s);
	}
	fclose(wfile);
	free(tempo_latt.cell);
}

// Plot Twiss parameter gamma from beta file (needs to generate first).
extern void plot_twiss_gamma(char *filename, char *with1, char *with2, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,30.'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 12\n");
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,30.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,30.'\n", ylabel);
	fprintf(gp,"set mxtics 5\nset mytics 5\n");
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color \n");
		fprintf(gp," set output '%s' \n", psfilename);
		fprintf(gp,"plot '%s' using 1:(1/$2*(1+$3**2)) with %s, '%s' using 1:(1/$4*(1+$5**2)) with %s\n", filename, with1, filename, with2);
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using 1:(1/$2*(1+$3**2)) with %s, '%s' using 1:(1/$4*(1+$5**2)) with %s\n", filename, with1, filename, with2);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

// Plot Twiss parameter gamma from beta file (needs to generate first).
extern void plot_twiss_gamma_w_efb(char *filename, char *file_efb, double efb_scale, char *with1, char *with2, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,30.'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 12\n");
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,30.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,30.'\n", ylabel);
	fprintf(gp,"set mxtics 5\nset mytics 5\n");
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		fprintf(gp," set terminal postscript eps enhanced color \n");
		fprintf(gp," set output '%s' \n", psfilename);
		fprintf(gp,"plot '%s' using 1:(1/$2*(1+$3**2)) with %s, '%s' using 1:(1/$4*(1+$5**2)) with %s, '%s' using 1:($2*%lf) with lines lt 1 lc 0 lw 1\n", filename, with1, filename, with2, file_efb, efb_scale);
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using 1:(1/$2*(1+$3**2)) with %s, '%s' using 1:(1/$4*(1+$5**2)) with %s, '%s' using 1:($2*%lf) with lines lt 1 lc 0 lw 1\n", filename, with1, filename, with2, file_efb, efb_scale);
	//fprintf(gp,"plot '%s' using 1:(1/$2*(1+$3**2)) with %s, '%s' using 1:(1/$4*(1+$5**2)) with %s\n", filename, with1, filename, with2);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

	// ************************************************************************************ //
	//							     "plot histogram" functions								//
	// ************************************************************************************ //

extern void plot_one_histogram(char *txtfilename, double xmin, double xmax, double nb_bin, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	if(psfilename != NULL) fprintf(gp," set output '%s' \n", psfilename);
	fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	//fprintf(gp,"width = %lf\n", width);
	//fprintf(gp,"hist(x,width)=width*floor(x/width)+width/2.0;");
	fprintf(gp,"min = %lf\n", xmin);
	fprintf(gp,"max = %lf\n", xmax);
	fprintf(gp,"n = %lf\n", nb_bin);
	fprintf(gp,"width = (max-min)/n\n");
	fprintf(gp,"hist(x,width)=width*(floor((x-min)/width)+0.5) + min;");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	//fprintf(gp,"set offset graph 0.05,0.05,0.05,0.0\n");
	fprintf(gp,"set style fill solid 1.0\n");
	fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red'\n", txtfilename);
	//fprintf(gp,"plot '%s' u (hist(($1*1.e-9),width)):(1.0) smooth freq with boxes lc rgb'red'\n", txtfilename);
	fprintf(gp," set output\n");
	fprintf(gp," set terminal x11\n");
	fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red'\n", txtfilename);
	//fprintf(gp,"plot '%s' u (hist(($1*1.e-9),width)):(1.0) smooth freq with boxes lc rgb'red'\n", txtfilename);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

/*extern void plot_two_histogram(char *txtfilename, char *txtfilename2, double xmin, double xmax, double nb_bin, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
 	FILE *gp;
 
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set terminal x11\n");
	if(psfilename != NULL) fprintf(gp," set output '%s' \n", psfilename);
	fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set multiplot\n");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp,"min = %lf\n", xmin);
	fprintf(gp,"max = %lf\n", xmax);
	fprintf(gp,"n = %lf\n", nb_bin);
	fprintf(gp,"width = (max-min)/n\n");
	fprintf(gp,"hist(x,width)=width*(floor((x-min)/width)+0.5) + min;");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	//fprintf(gp,"set offset graph 0.05,0.05,0.05,0.0\n");
	fprintf(gp,"set style fill solid 1.0 \n");
	fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red'\n", txtfilename, txtfilename2);
	fprintf(gp," set output\n");
	//fprintf(gp," set multiplot\n");
	//fprintf(gp," set terminal x11\n");
	//fprintf(gp,"set offset graph 0.05,0.05,0.05,0.0\n");
	//fprintf(gp,"set style fill solid 1.0 \n");
	//fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red'\n", txtfilename, txtfilename2);
	pclose(gp);
}//*/

extern void plot_two_histogram(char *txtfilename, char *title1, char *txtfilename2, char *title2, double xmin, double xmax, double nb_bin, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
 	FILE *gp;
 
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set terminal x11\n");
	if(psfilename != NULL) fprintf(gp," set output '%s' \n", psfilename);
	fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set multiplot\n");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp,"min = %lf\n", xmin);
	fprintf(gp,"max = %lf\n", xmax);
	fprintf(gp,"n = %lf\n", nb_bin);
	fprintf(gp,"width = (max-min)/n\n");
	fprintf(gp,"hist(x,width)=width*(floor((x-min)/width)+0.5) + min;");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	if(title1 == NULL && title2 == NULL) fprintf(gp,"unset key\n");
	//fprintf(gp,"set offset graph 0.05,0.05,0.05,0.0\n");
	fprintf(gp,"set style fill solid 1.0 \n");
	//fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red'\n", txtfilename, txtfilename2);
	if(title1 == NULL && title2 == NULL) fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red'\n", txtfilename, txtfilename2);
	else {
		printf("with title.. %s and %s\n", title1, title2);
		fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue' title '%s', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red' title '%s'\n", txtfilename, title1, txtfilename2, title2);
	}
	//fprintf(gp," set output\n");
	//fprintf(gp," set multiplot\n");
	//fprintf(gp," set terminal x11\n");
	//fprintf(gp,"set offset graph 0.05,0.05,0.05,0.0\n");
	//fprintf(gp,"set style fill solid 1.0 \n");
	//fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red'\n", txtfilename, txtfilename2);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void plot_three_histogram(char *txtfilename, char *title1, char *txtfilename2, char *title2, char *txtfilename3, char *title3, double xmin, double xmax, double nb_bin, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
 	FILE *gp;
 
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set terminal x11\n");
	if(psfilename != NULL) fprintf(gp," set output '%s' \n", psfilename);
	fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set multiplot\n");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp,"min = %lf\n", xmin);
	fprintf(gp,"max = %lf\n", xmax);
	fprintf(gp,"n = %lf\n", nb_bin);
	fprintf(gp,"width = (max-min)/n\n");
	fprintf(gp,"hist(x,width)=width*(floor((x-min)/width)+0.5) + min;");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	if(title1 == NULL && title2 == NULL && title3 == NULL)fprintf(gp,"unset key\n");
	fprintf(gp,"set offset graph 0.05,0.05,0.05,0.0\n");
	fprintf(gp,"set style fill solid 1.0 \n");
	if(title1 == NULL && title2 == NULL && title3 == NULL) fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'green'\n", txtfilename, txtfilename2, txtfilename3);
	else fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue' title '%s', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red' title '%s', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'forest-green' title '%s'\n", txtfilename, title1, txtfilename2, title2, txtfilename3, title3);
	fprintf(gp," set output\n");
	//fprintf(gp," set multiplot\n");
	//fprintf(gp," set terminal x11\n");
	//fprintf(gp,"set offset graph 0.05,0.05,0.05,0.0\n");
	//fprintf(gp,"set style fill solid 1.0 \n");
	//fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red'\n", txtfilename, txtfilename2);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void plot_four_histogram(char *txtfilename, char *title1, char *txtfilename2, char *title2, char *txtfilename3, char *title3, char *txtfilename4, char *title4, double xmin, double xmax, double nb_bin, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
 	FILE *gp;
 
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set terminal x11\n");
	if(psfilename != NULL) fprintf(gp," set output '%s' \n", psfilename);
	fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set multiplot\n");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp,"min = %lf\n", xmin);
	fprintf(gp,"max = %lf\n", xmax);
	fprintf(gp,"n = %lf\n", nb_bin);
	fprintf(gp,"width = (max-min)/n\n");
	fprintf(gp,"hist(x,width)=width*(floor((x-min)/width)+0.5) + min;");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	if(title1 == NULL && title2 == NULL && title3 == NULL && title4 == NULL)fprintf(gp,"unset key\n");
	fprintf(gp,"set offset graph 0.05,0.05,0.05,0.0\n");
	fprintf(gp,"set style fill solid 1.0 \n");
	if(title1 == NULL && title2 == NULL && title3 == NULL && title4 == NULL) fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'forest-green', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'orange'\n", txtfilename, txtfilename2, txtfilename3, txtfilename4);
	else fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue' title '%s', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red' title '%s', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'forest-green' title '%s', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'orange' title '%s'\n", txtfilename, title1, txtfilename2, title2, txtfilename3, title3, txtfilename4, title4);
	fprintf(gp," set output\n");
	//fprintf(gp," set multiplot\n");
	//fprintf(gp," set terminal x11\n");
	//fprintf(gp,"set offset graph 0.05,0.05,0.05,0.0\n");
	//fprintf(gp,"set style fill solid 1.0 \n");
	//fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red'\n", txtfilename, txtfilename2);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void plot_five_histogram(char *txtfilename, char *title1, char *txtfilename2, char *title2, char *txtfilename3, char *title3, char *txtfilename4, char *title4, char *txtfilename5, char *title5, double xmin, double xmax, double nb_bin, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
 	FILE *gp;
 
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set terminal x11\n");
	if(psfilename != NULL) fprintf(gp," set output '%s' \n", psfilename);
	fprintf(gp," set terminal postscript color\n");
	fprintf(gp," set multiplot\n");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp,"min = %lf\n", xmin);
	fprintf(gp,"max = %lf\n", xmax);
	fprintf(gp,"n = %lf\n", nb_bin);
	fprintf(gp,"width = (max-min)/n\n");
	fprintf(gp,"hist(x,width)=width*(floor((x-min)/width)+0.5) + min;");
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	if(title1 == NULL && title2 == NULL && title3 == NULL && title4 == NULL && title5 == NULL)fprintf(gp,"unset key\n");
	fprintf(gp,"set offset graph 0.05,0.05,0.05,0.0\n");
	fprintf(gp,"set style fill solid 1.0 \n");
	if(title1 == NULL && title2 == NULL && title3 == NULL && title4 == NULL && title5 == NULL) fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'magenta', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'forest-green', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'orange'\n", txtfilename, txtfilename2, txtfilename3, txtfilename4, txtfilename5);
	else fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue' title '%s', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'magenta' title '%s', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red' title '%s', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'forest-green' title '%s', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'orange' title '%s'\n", txtfilename, title1, txtfilename2, title2, txtfilename3, title3, txtfilename4, title4, txtfilename5, title5);
	fprintf(gp," set output\n");
	//fprintf(gp," set multiplot\n");
	//fprintf(gp," set terminal x11\n");
	//fprintf(gp,"set offset graph 0.05,0.05,0.05,0.0\n");
	//fprintf(gp,"set style fill solid 1.0 \n");
	//fprintf(gp,"plot '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'blue', '%s' u (hist($1,width)):(1.0) smooth freq with boxes lc rgb'red'\n", txtfilename, txtfilename2);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

	// ************************************************************************************ //
	//							    "plot trajectory" functions								//
	// ************************************************************************************ //

extern void plot_traj(struct Lattice *latt, char *trackout, char *aspectfile, char *xcolumn_out, char *ycolumn_out, char *xcolumn_asp, char *ycolumn_asp, char *with_out, char *with_asp, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	aspectview(aspectfile, latt);
	easyplot2(trackout, aspectfile, xcolumn_out, ycolumn_out, xcolumn_asp, ycolumn_asp, with_out, with_asp, NULL, NULL, xlabel, ylabel, xrange,yrange, psfilename, setoption);
	//aspect_fringe_view("data/aspect_en.dat", "data/aspect_ex.dat", latt);
	//easyplot4(trackout, aspectfile, "data/aspect_en.dat", "data/aspect_ex.dat", xcolumn_out, ycolumn_out, xcolumn_asp, ycolumn_asp, "2", "1", "2", "1", with_out, with_asp, "lines dt 2 lc 2 lw 2", "lines dt 2 lc 3 lw 2", NULL, NULL, NULL, NULL, xlabel, ylabel, xrange,yrange, psfilename, setoption);
}

extern void aspectview(char *filename, struct Lattice *latt)
{
	int i, j;
	struct Lattice tempo_latt;
	
	//copy a Lattice, with cell[] stored at a diffent memory address, 
	//!!!and DO NOT FORGET to free the allocated memory before exiting aspectview!!!
	copy_latt(&tempo_latt, latt);
	emptyfile(filename);
	for(i = 0; i < tempo_latt.periodicity; i++) {
		
		for(j = 0; j < tempo_latt.nbcell ;j++) {
			if (strcmp(tempo_latt.cell[j].keyword, "ffag-r-he") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "ffag-r-lin") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "ffag-r-enge") == 0) aspect_radial(filename, &(tempo_latt.cell[j]));
			else if (strcmp(tempo_latt.cell[j].keyword, "ffag-tilt-he") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "ffag-tilt-lin") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "ffag-tilt-enge") == 0 ||
				strcmp(tempo_latt.cell[j].keyword, "ffag-tiltpol-he") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "ffag-tiltpol-lin") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "ffag-tiltpol-enge") == 0) aspect_tilt(filename, &(tempo_latt.cell[j]));
			else if (strcmp(tempo_latt.cell[j].keyword, "ffag-spi-he") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "ffag-spi-lin") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "ffag-spi-enge") == 0 ||
				strcmp(tempo_latt.cell[j].keyword, "ffag-spi-fullenge") == 0) aspect_spiral(filename, &(tempo_latt.cell[j]));
			else if (strcmp(tempo_latt.cell[j].keyword, "ffag-s-lin") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "ffag-s-he") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "ffag-s-enge") == 0) aspect_straight(filename, &(tempo_latt.cell[j]));
			else if (strcmp(tempo_latt.cell[j].keyword, "ffag-sti-lin") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "ffag-sti-he") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "ffag-sti-enge") == 0) aspect_straighttilt(filename, &(tempo_latt.cell[j]));
			else if (strcmp(tempo_latt.cell[j].keyword, "realbend") == 0) aspect_realbend(filename, &(tempo_latt.cell[j]));
			else if (strcmp(tempo_latt.cell[j].keyword, "purebend-he") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "purebend-lin") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "purebend-enge") == 0) aspect_purebend(filename, &(tempo_latt.cell[j]));
			else if (test_cell_map(&(tempo_latt.cell[j])) == YES) aspect_boun(filename, &(tempo_latt.cell[j]));
			else if (strcmp(tempo_latt.cell[j].keyword, "quad-lin") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "quad-he") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "quad-enge") == 0) aspect_straight(filename, &(tempo_latt.cell[j]));
			else if(strcmp(latt->cell[i].keyword, "vffa-rect-enge") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-enge-add") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-enge-bx") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-enge-separ") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-str-enge") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-str-lin") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-lin") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-str-he") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-he") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-atan-str") == 0 ||
				strcmp(latt->cell[i].keyword, "vffa-rect-atan") == 0) aspect_rect_vffa(filename, &(tempo_latt.cell[j]));
			else if(strcmp(tempo_latt.cell[j].keyword, "vffa-sect-lin") == 0 || 
				strcmp(tempo_latt.cell[j].keyword, "vffa-sect-enge") == 0) aspect_sect_vffa(filename, &(tempo_latt.cell[j]));
		}
		//move framework
		latt_move_exit_to_entrance(&tempo_latt);
	}
	//free allocated memory
	free(tempo_latt.cell);
}

static void change_coord_vffa_sect(double *x, double *y, double r, double th, int comp_nb, struct Cell *cell)
{
	int sign_rho;
	double rc, thc, rho, xm, ym;
	
	thc = cell->mpara[comp_nb][0]; //[rad]
	rc = cell->mpara[comp_nb][1];
	rho = cell->mpara[comp_nb][5];
	sign_rho = sign(cell->mpara[comp_nb][5]);
	
	xm = r*cos(th); //rho, ffbe
	ym = r*sin(th);
	*x = (xm*sign_rho+rc-rho)*cos(thc) - ym*sin(thc);
	*y = (xm*sign_rho+rc-rho)*sin(thc) + ym*cos(thc);
	fwk_pt_loctoglob(x, y, &(cell->framework));
}

static void aspect_sect_vffa(char *filename, struct Cell *cell)
{
	int i, j, n;
	double rho, conv_lim, x, y, ffbe, ffbs, rmin, rmax, th;
	FILE *wfile;
	n = 20;
	
	wfile = fopen(filename, "a");
	for(j = 0; j < cell->nbcomp; j++) {
		if(cell->mpara[j][0]>0-TINYLENGTH && cell->mpara[j][0]<cell->boun.thmax+TINYLENGTH) {
			rho = fabs(cell->mpara[j][5]);
			conv_lim = cell->mpara[j][6]; //convergence limit [m]
			if(conv_lim<TINYLENGTH) errorstop("convergence limit <=0, have you initialized it?");
			ffbe = cell->efben[j][0]; //[rad]
			ffbs = cell->efbex[j][0]; //[rad]
			rmin = MAX(0,rho-conv_lim);
			rmax = rho+conv_lim;
			
			change_coord_vffa_sect(&x, &y, rho, ffbe, j, cell);
			fprintf(wfile, "%le  %le\n",x,y);
			for(i=0;i<n+1;i++) { //rmax, ffbe->ffbs
				th = ffbe + i*(ffbs-ffbe)/n;
				change_coord_vffa_sect(&x, &y, rmax, th, j, cell);
				fprintf(wfile, "%le  %le\n",x,y);
			}
			fprintf(wfile, "%le  %le\n",x,y);
			for(i=0;i<n+1;i++) { //rho, ffbs->ffbe
				th = ffbs + i*(ffbe-ffbs)/n;
				change_coord_vffa_sect(&x, &y, rho, th, j, cell);
				fprintf(wfile, "%le  %le\n",x,y);
			}
			change_coord_vffa_sect(&x, &y, rmin, ffbe, j, cell);
			fprintf(wfile, "%le  %le\n",x,y);
			if(rmin>TINYLENGTH) { //rmin, ffbe->ffbs
				for(i=1;i<n+1;i++) {
					th = ffbe + i*(ffbs-ffbe)/n;
					change_coord_vffa_sect(&x, &y, rmin, th, j, cell);
					fprintf(wfile, "%le  %le\n",x,y);
				}
			}
			change_coord_vffa_sect(&x, &y, rho, ffbs, j, cell);
			fprintf(wfile, "%le  %le\n",x,y);
			fprintf(wfile, "\n\n");
		}
	}
	fclose(wfile);
}

/*static void aspect_sect_vffa(char *filename, struct Cell *cell)
{
	int i, j, n, sign_rho;
	double thc, rc, rho, conv_lim, xm, ym, x, y, ffbe, ffbs, rmin, rmax, th;
	FILE *wfile;
	n = 20;
	
	wfile = fopen(filename, "a");
	for(j = 0; j < cell->nbcomp; j++) {
		if(cell->mpara[j][0]>0-TINYLENGTH && cell->mpara[j][0]<cell->boun.thmax+TINYLENGTH) {
			thc = cell->mpara[j][0]; //[rad]
			rc = cell->mpara[j][1];
			sign_rho = sign(cell->mpara[j][5]);
			rho = cell->mpara[j][5];
			conv_lim = cell->mpara[j][6]; //convergence limit [m]
			if(conv_lim<TINYLENGTH) errorstop("convergence limit <=0, have you initialized it?");
			ffbe = cell->efben[j][0]; //[rad]
			ffbs = cell->efbex[j][0]; //[rad]
			rmin = MAX(0,fabs(rho)-conv_lim);
			rmax = fabs(rho)+conv_lim;
			
			xm = fabs(rho)*cos(ffbe); //rho, ffbe
			ym = fabs(rho)*sin(ffbe);
			x = (xm*sign_rho+rc-rho)*cos(thc) - ym*sin(thc);
			y = (xm*sign_rho+rc-rho)*sin(thc) + ym*cos(thc);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			for(i=0;i<n+1;i++) { //rmax, ffbe->ffbs
				th = ffbe + i*(ffbs-ffbe)/n;
				xm = rmax*cos(th);
				ym = rmax*sin(th);
				x = (xm*sign_rho+rc-rho)*cos(thc) - ym*sin(thc);
				y = (xm*sign_rho+rc-rho)*sin(thc) + ym*cos(thc);
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le\n",x,y);
			}
			xm = fabs(rho)*cos(ffbs); //rho, ffbs
			ym = fabs(rho)*sin(ffbs);
			x = (xm*sign_rho+rc-rho)*cos(thc) - ym*sin(thc);
			y = (xm*sign_rho+rc-rho)*sin(thc) + ym*cos(thc);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			for(i=1;i<n+1;i++) { //rho, ffbs->ffbe
				th = ffbs + i*(ffbe-ffbs)/n;
				xm = fabs(rho)*cos(th);
				ym = fabs(rho)*sin(th);
				x = (xm*sign_rho+rc-rho)*cos(thc) - ym*sin(thc);
				y = (xm*sign_rho+rc-rho)*sin(thc) + ym*cos(thc);
				fwk_pt_loctoglob(&x, &y, &(cell->framework));
				fprintf(wfile, "%le  %le\n",x,y);
			}
			xm = rmin*cos(ffbe); //rmin, ffbe
			ym = rmin*sin(ffbe);
			x = (xm*sign_rho+rc-rho)*cos(thc) - ym*sin(thc);
			y = (xm*sign_rho+rc-rho)*sin(thc) + ym*cos(thc);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			if(rmin>TINYLENGTH) { //rmin, ffbe->ffbs
				for(i=1;i<n+1;i++) {
					th = ffbe + i*(ffbs-ffbe)/n;
					xm = rmin*cos(th);
					ym = rmin*sin(th);
					x = (xm*sign_rho+rc-rho)*cos(thc) - ym*sin(thc);
					y = (xm*sign_rho+rc-rho)*sin(thc) + ym*cos(thc);
					fwk_pt_loctoglob(&x, &y, &(cell->framework));
					fprintf(wfile, "%le  %le\n",x,y);
				}
			}
			xm = fabs(rho)*cos(ffbs); //rho, ffbs
			ym = fabs(rho)*sin(ffbs);
			x = (xm*sign_rho+rc-rho)*cos(thc) - ym*sin(thc);
			y = (xm*sign_rho+rc-rho)*sin(thc) + ym*cos(thc);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			fprintf(wfile, "\n");
		}
	}
	fclose(wfile);
}
//*/


static void aspect_rect_vffa(char *filename, struct Cell *cell)
{
	int j, doyouplotcomp=YES, doyoulimitcomp=NO;
	double x, y, x_0, y_0, ffbe, ffbs, lambda_e, angle;//, lambda_s;
	FILE *wfile;
	
	
	wfile = fopen(filename, "a");
	for(j = 0; j < cell->nbcomp; j++) {
		if(doyoulimitcomp==YES) {
			if(cell->boun.thmax!=0) {
				if(cell->mpara[j][0]>0-TINYLENGTH && cell->mpara[j][0]<cell->boun.thmax+TINYLENGTH) doyouplotcomp = YES;
				else doyouplotcomp = NO;
			}
			else {
				if(cell->mpara[j][0] > 0-TINYLENGTH && cell->mpara[j][0] < cell->boun.ymax+TINYLENGTH) doyouplotcomp = YES;
				else doyouplotcomp = NO;
			}
		}
		//printf("comp %i, doyouplot=%i, YES=%i\n", j, doyouplotcomp, YES);
		if(doyouplotcomp==YES) {
			if(cell->boun.thmax!=0) {
				x_0 = cell->mpara[j][1]*cos(cell->mpara[j][0]);
				y_0 = cell->mpara[j][1]*sin(cell->mpara[j][0]);
			}
			else {
				x_0 = cell->mpara[j][1];
				y_0 = cell->mpara[j][0];
			}
			ffbe = cell->efben[j][0];
			ffbs = cell->efbex[j][0];
			//lambda_e = cell->efben[j][1];
			//lambda_s = cell->efbex[j][1];
			lambda_e = cell->mpara[j][6];
			angle = cell->mpara[j][5];
			if(cell->boun.thmax!=0) angle +=cell->mpara[j][0];
			
			x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle);
			y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = x_0 - ffbe*sin(angle) + lambda_e*cos(angle);
			y = y_0 + ffbe*cos(angle) + lambda_e*sin(angle);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			
			x = x_0 - ffbe*sin(angle);
			y = y_0 + ffbe*cos(angle);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			
			x = x_0 - ffbs*sin(angle);
			y = y_0 + ffbs*cos(angle);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			
			x = x_0 - ffbe*sin(angle);
			y = y_0 + ffbe*cos(angle);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			
			x = x_0 - ffbe*sin(angle) - lambda_e*cos(angle);
			y = y_0 + ffbe*cos(angle) - lambda_e*sin(angle);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = x_0 - ffbs*sin(angle) - lambda_e*cos(angle);
			y = y_0 + ffbs*cos(angle) - lambda_e*sin(angle);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = x_0 - ffbs*sin(angle) + lambda_e*cos(angle);
			y = y_0 + ffbs*cos(angle) + lambda_e*sin(angle);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);//*/
			fprintf(wfile, "\n\n");
		}
	}
	fclose(wfile);
}

static void aspect_boun(char *filename, struct Cell *cell)
{
	int i, n;
	double x, y, rmin, rmax, th, thmax, thstep;
	FILE *wfile;
	
	n = 50;
	
	rmin = cell->collim.rmin;
	rmax = cell->collim.rmax;
	
	//test = (rmax - 1/n*(rmax-rmin));
	//printf("test = %lf\n", test);
	wfile = fopen(filename, "a");
	
	x = rmin;
	y = 0.;
	fwk_pt_loctoglob(&x, &y, &(cell->framework));
	fprintf(wfile, "%lf  %lf\n",x,y);
	
	if(cell->boun.thmax != 0) {
		thmax = cell->boun.thmax;
		thstep = thmax/n;
		//fprintf(wfile, "0  0\n");
		for(i = 0; i < n+1; i++) {
			th = i*thstep;
			x = rmax*cos(th);
			y = rmax*sin(th);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		}
		//fprintf(wfile, "0  0\n\n");
		x = rmin*cos(thmax);
		y = rmin*sin(thmax);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%lf  %lf\n",x,y);
		for(i = n; i > -1; i--) {
			th = i*thstep;
			x = rmin*cos(th);
			y = rmin*sin(th);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		}
	}
	else {
		x = rmax;
		y = 0.;
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %lf\n",x,y);
		x = rmax;
		y = cell->boun.ymax;
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		x = rmin;
		y = cell->boun.ymax;
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		x = rmin;
		y = 0.;
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%lf  %lf\n",x,y);
	}
	fprintf(wfile, "\n\n");
	
	if(cell->instrutype != NO) {
		printf("instrutype != NO\n");
	}
	fclose(wfile);
}

static void aspect_radial(char *filename, struct Cell *cell)
{
	int j;
	double x, y, rmin, rmax;
	FILE *wfile;
	
	rmin = cell->collim.rmin;
	rmax = cell->collim.rmax;
	
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		if(cell->mpara[j][0]>0-TINYLENGTH && cell->mpara[j][0]<cell->boun.thmax+TINYLENGTH) {
			x = rmin*cos(cell->mpara[j][0] + cell->efben[j][0]);
			y = rmin*sin(cell->mpara[j][0] + cell->efben[j][0]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = rmax*cos(cell->mpara[j][0] + cell->efben[j][0]);
			y = rmax*sin(cell->mpara[j][0] + cell->efben[j][0]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = rmax*cos(cell->mpara[j][0] + cell->efbex[j][0]);
			y = rmax*sin(cell->mpara[j][0] + cell->efbex[j][0]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = rmin*cos(cell->mpara[j][0] + cell->efbex[j][0]);
			y = rmin*sin(cell->mpara[j][0] + cell->efbex[j][0]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = rmin*cos(cell->mpara[j][0] + cell->efben[j][0]);
			y = rmin*sin(cell->mpara[j][0] + cell->efben[j][0]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			fprintf(wfile, "\n\n");
		}
	}
	
	if(cell->instrutype != NO) {
		x = (rmin-0.2)*cos(cell->instru.thmax);
		y = (rmin-0.2)*sin(cell->instru.thmax);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = (rmax+0.2)*cos(cell->instru.thmax);
		y = (rmax+0.2)*sin(cell->instru.thmax);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		fprintf(wfile, "\n\n");
	}
	
	fclose(wfile);
}

static void aspect_straight(char *filename, struct Cell *cell)
{
	int j;
	double x, y, xmin, xmax;
	FILE *wfile;
	
	xmin = cell->collim.rmin;
	xmax = cell->collim.rmax;
	
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		if(cell->mpara[j][0]>0 && cell->mpara[j][0]<cell->boun.ymax) {
			x = xmin;
			y = cell->mpara[j][0] + cell->efben[j][0];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = xmax;
			y = cell->mpara[j][0] + cell->efben[j][0];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = xmax;
			y = cell->mpara[j][0] + cell->efbex[j][0];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = xmin;
			y = cell->mpara[j][0] + cell->efbex[j][0];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = xmin;
			y = cell->mpara[j][0] + cell->efben[j][0];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			fprintf(wfile, "\n\n");
		}
	}
	
	fclose(wfile);
}

static void aspect_realbend(char *filename, struct Cell *cell)
{
	double x, y, x0, y0, ffbe, ffbs, xm, xp;
	FILE *wfile;
	
	x0 = cell->mpara[0][0]*cos(cell->boun.thmax/2.);
	y0 = cell->mpara[0][0]*sin(cell->boun.thmax/2.);
	ffbe = y0 + cell->efben[0][0];
	ffbs = y0 + cell->efbex[0][0];
	xm = x0 - cell->mpara[0][1];
	xp = x0 + cell->mpara[0][1];
	
	wfile = fopen(filename, "a");
	
	x = x0+(xm-x0)*cos(-cell->boun.thmax/2.)+(ffbe-y0)*sin(-cell->boun.thmax/2.);
	y = y0-(xm-x0)*sin(-cell->boun.thmax/2.)+(ffbe-y0)*cos(-cell->boun.thmax/2.);
	fwk_pt_loctoglob(&x, &y, &(cell->framework));
	fprintf(wfile, "%le  %le\n",x,y);
		
	x = x0+(xm-x0)*cos(-cell->boun.thmax/2.)+(ffbs-y0)*sin(-cell->boun.thmax/2.);
	y = y0-(xm-x0)*sin(-cell->boun.thmax/2.)+(ffbs-y0)*cos(-cell->boun.thmax/2.);
	fwk_pt_loctoglob(&x, &y, &(cell->framework));
	fprintf(wfile, "%le  %le\n",x,y);
		
	x = x0+(xp-x0)*cos(-cell->boun.thmax/2.)+(ffbs-y0)*sin(-cell->boun.thmax/2.);
	y = y0-(xp-x0)*sin(-cell->boun.thmax/2.)+(ffbs-y0)*cos(-cell->boun.thmax/2.);
	fwk_pt_loctoglob(&x, &y, &(cell->framework));
	fprintf(wfile, "%le  %le\n",x,y);
		
	x = x0+(xp-x0)*cos(-cell->boun.thmax/2.)+(ffbe-y0)*sin(-cell->boun.thmax/2.);
	y = y0-(xp-x0)*sin(-cell->boun.thmax/2.)+(ffbe-y0)*cos(-cell->boun.thmax/2.);
	fwk_pt_loctoglob(&x, &y, &(cell->framework));
	fprintf(wfile, "%le  %le\n",x,y);
		
	x = x0+(xm-x0)*cos(-cell->boun.thmax/2.)+(ffbe-y0)*sin(-cell->boun.thmax/2.);
	y = y0-(xm-x0)*sin(-cell->boun.thmax/2.)+(ffbe-y0)*cos(-cell->boun.thmax/2.);
	fwk_pt_loctoglob(&x, &y, &(cell->framework));
	fprintf(wfile, "%le  %le\n",x,y);
	fprintf(wfile, "\n\n");
	
	/*
	x = 0;
	y = 0;
	fwk_pt_loctoglob(&x, &y, &(cell->framework));
	fprintf(wfile, "%le  %le\n",x,y);
	
	x = cell->collim.rmax;
	y = 0;
	fwk_pt_loctoglob(&x, &y, &(cell->framework));
	fprintf(wfile, "%le  %le\n",x,y);
	
	x = cell->collim.rmax*cos(cell->boun.thmax);
	y = cell->collim.rmax*sin(cell->boun.thmax);
	fwk_pt_loctoglob(&x, &y, &(cell->framework));
	fprintf(wfile, "%le  %le\n",x,y);
	
	x = 0;
	y = 0;
	fwk_pt_loctoglob(&x, &y, &(cell->framework));
	fprintf(wfile, "%le  %le\n",x,y);
	fprintf(wfile, "\n");
	//*/
	fclose(wfile);
}

static void aspect_purebend(char *filename, struct Cell *cell)
{
	int j;
	double x, y, rmin, rmax, rc;
	FILE *wfile;
	
	rmin = cell->collim.rmin;
	rmax = cell->collim.rmax;
	rc = cell->mpara[0][3];
	
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		x = rc*cos(cell->mpara[j][0]) + (rmin-rc)*cos(cell->mpara[j][0] + cell->efben[j][0]);
		y = rc*sin(cell->mpara[j][0]) + (rmin-rc)*sin(cell->mpara[j][0] + cell->efben[j][0]);
		//x = rc*cos(cell->mpara[j][0]);
		//y = rc*sin(cell->mpara[j][0]);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rc*cos(cell->mpara[j][0]) + (rmax-rc)*cos(cell->mpara[j][0] + cell->efben[j][0]);
		y = rc*sin(cell->mpara[j][0]) + (rmax-rc)*sin(cell->mpara[j][0] + cell->efben[j][0]);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rc*cos(cell->mpara[j][0]) + (rmax-rc)*cos(cell->mpara[j][0] + cell->efbex[j][0]);
		y = rc*sin(cell->mpara[j][0]) + (rmax-rc)*sin(cell->mpara[j][0] + cell->efbex[j][0]);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rc*cos(cell->mpara[j][0]) + (rmin-rc)*cos(cell->mpara[j][0] + cell->efbex[j][0]);
		y = rc*sin(cell->mpara[j][0]) + (rmin-rc)*sin(cell->mpara[j][0] + cell->efbex[j][0]);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rc*cos(cell->mpara[j][0]) + (rmin-rc)*cos(cell->mpara[j][0] + cell->efben[j][0]);
		y = rc*sin(cell->mpara[j][0]) + (rmin-rc)*sin(cell->mpara[j][0] + cell->efben[j][0]);
		//x = rc*cos(cell->mpara[j][0]);
		//y = rc*sin(cell->mpara[j][0]);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		fprintf(wfile, "\n\n");
		/*
		x = 0.;
		y = 0.;
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmax;
		y = 0.;
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmax*cos(cell->boun.thmax);
		y = rmax*sin(cell->boun.thmax);;
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = 0.;
		y = 0.;
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		fprintf(wfile, "\n");//*/
	}
	
	if(cell->instrutype != NO) {
		errorstop("don't know what's going on in aspect_purebend");
		x = (rmin-0.2)*cos(cell->instru.thmax);
		y = (rmin-0.2)*sin(cell->instru.thmax);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = (rmax+0.2)*cos(cell->instru.thmax);
		y = (rmax+0.2)*sin(cell->instru.thmax);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		fprintf(wfile, "\n\n");
	}
	
	fclose(wfile);
}

static void aspect_spiral(char *filename, struct Cell *cell)
{
	int j, i, n;
	double x, y, rmin, rmax;
	//double test;
	FILE *wfile;
	
	n = 20;
	
	rmin = cell->collim.rmin;
	rmax = cell->collim.rmax;
	//test = (rmax - 1/n*(rmax-rmin));
	//printf("test = %lf\n", test);
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		for(i = 0; i < n+1; i++) {
			x = (rmin + i*(rmax-rmin)/n)*cos(cell->mpara[j][0] + cell->efben[j][0] + tan(cell->mpara[j][4])*log((rmin + i*(rmax-rmin)/n)/cell->mpara[j][1]));
			y = (rmin + i*(rmax-rmin)/n)*sin(cell->mpara[j][0] + cell->efben[j][0] + tan(cell->mpara[j][4])*log((rmin + i*(rmax-rmin)/n)/cell->mpara[j][1]));
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		}
		for(i = 0; i < n+1; i++) {
			x = (rmax - i*(rmax-rmin)/n)*cos(cell->mpara[j][0] + cell->efbex[j][0] + tan(cell->mpara[j][4])*log((rmax - i*(rmax-rmin)/n)/cell->mpara[j][1]));
			y = (rmax - i*(rmax-rmin)/n)*sin(cell->mpara[j][0] + cell->efbex[j][0] + tan(cell->mpara[j][4])*log((rmax - i*(rmax-rmin)/n)/cell->mpara[j][1]));
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		}
		
		x = rmin*cos(cell->mpara[j][0] + cell->efben[j][0] + tan(cell->mpara[j][4])*log((rmin)/cell->mpara[j][1]));
		y = rmin*sin(cell->mpara[j][0] + cell->efben[j][0] + tan(cell->mpara[j][4])*log((rmin)/cell->mpara[j][1]));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		fprintf(wfile, "\n\n");
	}
	
	if(cell->instrutype != NO) {
		printf("instrutype != NO\n");
	}
	
	fclose(wfile);
}

static void aspect_tilt(char *filename, struct Cell *cell)
{
	int j;
	double x, y, rmin, rmax, angle;
	//double test;
	FILE *wfile;
	
	rmin = cell->collim.rmin;
	rmax = cell->collim.rmax;
	
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		if (strcmp(cell->keyword, "ffag-tilt-he") == 0 || 
		strcmp(cell->keyword, "ffag-tilt-lin") == 0 || 
		strcmp(cell->keyword, "ffag-tilt-enge") == 0) angle = cell->mpara[j][5];
		if(strcmp(cell->keyword, "add_field_comp_FFAGtiltpol_lin") == 0 ||
		strcmp(cell->keyword, "add_field_comp_FFAGtiltpol_he") == 0 ||
		strcmp(cell->keyword, "add_field_comp_FFAGtiltpol_enge") == 0)  angle = cell->mpara[j][2];
		else angle = cell->mpara[j][2];
		//test = (cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin);
		//printf("(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin) = %lf\n", test);
		x = rmin*cos(cell->mpara[j][0] + cell->efben[j][0] + angle - asin(cell->mpara[j][1]*sin(angle)/rmin));
		y = rmin*sin(cell->mpara[j][0] + cell->efben[j][0] + angle - asin(cell->mpara[j][1]*sin(angle)/rmin));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmax*cos(cell->mpara[j][0] + cell->efben[j][0] + angle - asin(cell->mpara[j][1]*sin(angle)/rmax));
		y = rmax*sin(cell->mpara[j][0] + cell->efben[j][0] + angle - asin(cell->mpara[j][1]*sin(angle)/rmax));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmax*cos(cell->mpara[j][0] + cell->efbex[j][0] + angle - asin(cell->mpara[j][1]*sin(angle)/rmax));
		y = rmax*sin(cell->mpara[j][0] + cell->efbex[j][0] + angle - asin(cell->mpara[j][1]*sin(angle)/rmax));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmin*cos(cell->mpara[j][0] + cell->efbex[j][0] + angle - asin(cell->mpara[j][1]*sin(angle)/rmin));
		y = rmin*sin(cell->mpara[j][0] + cell->efbex[j][0] + angle - asin(cell->mpara[j][1]*sin(angle)/rmin));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmin*cos(cell->mpara[j][0] + cell->efben[j][0] + angle - asin(cell->mpara[j][1]*sin(angle)/rmin));
		y = rmin*sin(cell->mpara[j][0] + cell->efben[j][0] + angle - asin(cell->mpara[j][1]*sin(angle)/rmin));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		fprintf(wfile, "\n\n");
	}
	
	if(cell->instrutype != NO) {
		printf("instrutype != NO\n");
	}
	
	fclose(wfile);
}

static void aspect_straighttilt(char *filename, struct Cell *cell)
{
	int j;
	double x, y, xmin, xmax;
	FILE *wfile;
	
	xmin = cell->collim.rmin;
	xmax = cell->collim.rmax;
	
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		x = xmin;
		y = cell->mpara[j][0] + cell->efben[j][0] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmin);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = xmax;
		y = cell->mpara[j][0] + cell->efben[j][0] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmax);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = xmax;
		y = cell->mpara[j][0] + cell->efbex[j][0] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmax);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = xmin;
		y = cell->mpara[j][0] + cell->efbex[j][0] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmin);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = xmin;
		y = cell->mpara[j][0] + cell->efben[j][0] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmin);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		fprintf(wfile, "\n\n");
	}
	
	fclose(wfile);
}

extern void aspect_fringe_view(char *filename_en, char *filename_ex, struct Lattice *latt)
{
	int i, j;
	struct Lattice tempo_latt;
	
	//copy a Lattice, with cell[] stored at a diffent memory address, 
	//!!!and DO NOT FORGET to free the allocated memory befor exiting aspectview!!!
	copy_latt(&tempo_latt, latt);
	 
	emptyfile(filename_en);
	emptyfile(filename_ex);
	for(i = 0; i < tempo_latt.periodicity; i++) {
		
		for(j = 0; j < tempo_latt.nbcell ;j++) {
			if (strcmp(tempo_latt.cell[j].keyword, "ffag-r-lin") == 0) {// || 
			//strcmp(tempo_latt.cell[j].keyword, "ffag-r-enge") == 0) {
				aspect_fringe_en_radial(filename_en, &(tempo_latt.cell[j]));
				aspect_fringe_ex_radial(filename_ex, &(tempo_latt.cell[j]));
			}
			else if (strcmp(tempo_latt.cell[j].keyword, "ffag-tilt-lin") == 0) {// || 
			//strcmp(tempo_latt.cell[j].keyword, "ffag-tilt-enge") == 0) {
				aspect_fringe_en_tilt(filename_en, &(tempo_latt.cell[j]));
				aspect_fringe_ex_tilt(filename_ex, &(tempo_latt.cell[j]));
			}
			else if (strcmp(tempo_latt.cell[j].keyword, "ffag-spi-lin") == 0) {// || 
			//strcmp(tempo_latt.cell[j].keyword, "ffag-spi-enge") == 0) {
				aspect_fringe_en_spiral(filename_en, &(tempo_latt.cell[j]));
				aspect_fringe_ex_spiral(filename_ex, &(tempo_latt.cell[j]));
			}
			else if (strcmp(tempo_latt.cell[j].keyword, "ffag-s-lin") == 0) {// || 
			//strcmp(tempo_latt.cell[j].keyword, "ffag-s-enge") == 0) {
				aspect_fringe_en_straight(filename_en, &(tempo_latt.cell[j]));
				aspect_fringe_ex_straight(filename_ex, &(tempo_latt.cell[j]));
			}
			else if (strcmp(tempo_latt.cell[j].keyword, "ffag-sti-lin") == 0){// || 
			//strcmp(tempo_latt.cell[j].keyword, "ffag-sti-enge") == 0) {
				aspect_fringe_en_straighttilt(filename_en, &(tempo_latt.cell[j]));
				aspect_fringe_ex_straighttilt(filename_ex, &(tempo_latt.cell[j]));
			}
			//else if (strcmp(tempo_latt.cell[j].keyword, "realbend") == 0) aspect_realbend(filename, &(tempo_latt.cell[j]));
			//else if (strcmp(tempo_latt.cell[j].keyword, "purebend-lin") == 0 || 
			//	strcmp(tempo_latt.cell[j].keyword, "purebend-enge") == 0) aspect_purebend(filename, &(tempo_latt.cell[j]));
			else if (strcmp(tempo_latt.cell[j].keyword, "quad-lin") == 0) {// || 
			//strcmp(tempo_latt.cell[j].keyword, "quad-enge") == 0) {
				aspect_fringe_en_straight(filename_en, &(tempo_latt.cell[j]));
				aspect_fringe_ex_straight(filename_ex, &(tempo_latt.cell[j]));
			}
		}
		//move framework
		latt_move_exit_to_entrance(&tempo_latt);
	}
	
	//free allocated memory
	free(tempo_latt.cell);
}

static void aspect_fringe_en_radial(char *filename, struct Cell *cell)
{
	int j;
	double x, y, rmin, rmax;
	FILE *wfile;
	
	rmin = cell->collim.rmin;
	rmax = cell->collim.rmax;
	
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		if(cell->mpara[j][0]>0-TINYLENGTH && cell->mpara[j][0]<cell->boun.thmax+TINYLENGTH) {
			//entrance
			x = rmin*cos(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1]);
			y = rmin*sin(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = rmax*cos(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1]);
			y = rmax*sin(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			
			x = rmin*cos(cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1]);
			y = rmin*sin(cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = rmax*cos(cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1]);
			y = rmax*sin(cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			
			x = rmin*cos(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1]);
			y = rmin*sin(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			fprintf(wfile, "\n\n");
		}
	}
	
	
	fclose(wfile);
}

static void aspect_fringe_ex_radial(char *filename, struct Cell *cell)
{
	int j;
	double x, y, rmin, rmax;
	FILE *wfile;
	
	rmin = cell->collim.rmin;
	rmax = cell->collim.rmax;
	
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		if(cell->mpara[j][0]>0-TINYLENGTH && cell->mpara[j][0]<cell->boun.thmax+TINYLENGTH) {
			//exit
			x = rmax*cos(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1]);
			y = rmax*sin(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = rmin*cos(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1]);
			y = rmin*sin(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			
			x = rmax*cos(cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1]);
			y = rmax*sin(cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = rmin*cos(cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1]);
			y = rmin*sin(cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = rmin*cos(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1]);
			y = rmin*sin(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1]);
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			fprintf(wfile, "\n\n");
		}
	}
	
	
	fclose(wfile);
}

static void aspect_fringe_en_straight(char *filename, struct Cell *cell)
{
	int j;
	double x, y, xmin, xmax;
	FILE *wfile;
	
	xmin = cell->collim.rmin;
	xmax = cell->collim.rmax;
	
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		if(cell->mpara[j][0]>0 && cell->mpara[j][0]<cell->boun.ymax) {
			x = xmin;
			y = cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = xmax;
			y = cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = xmax;
			y = cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = xmin;
			y = cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = xmin;
			y = cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			fprintf(wfile, "\n\n");
		}
	}
	fclose(wfile);
}

static void aspect_fringe_ex_straight(char *filename, struct Cell *cell)
{
	int j;
	double x, y, xmin, xmax;
	FILE *wfile;
	
	xmin = cell->collim.rmin;
	xmax = cell->collim.rmax;
	
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		if(cell->mpara[j][0]>0 && cell->mpara[j][0]<cell->boun.ymax) {
			x = xmin;
			y = cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = xmax;
			y = cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = xmax;
			y = cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = xmin;
			y = cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		
			x = xmin;
			y = cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1];
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
			fprintf(wfile, "\n\n");
		}
	}
	fclose(wfile);
}

static void aspect_fringe_en_spiral(char *filename, struct Cell *cell)
{
	int j, i, n;
	double x, y, rmin, rmax;
	//double test;
	FILE *wfile;
	
	n = 20;
	
	rmin = cell->collim.rmin;
	rmax = cell->collim.rmax;
	//test = (rmax - 1/n*(rmax-rmin));
	//printf("test = %lf\n", test);
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		for(i = 0; i < n+1; i++) {
			x = (rmin + i*(rmax-rmin)/n)*cos(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1] + tan(cell->mpara[j][4])*log((rmin + i*(rmax-rmin)/n)/cell->mpara[j][1]));
			y = (rmin + i*(rmax-rmin)/n)*sin(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1] + tan(cell->mpara[j][4])*log((rmin + i*(rmax-rmin)/n)/cell->mpara[j][1]));
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		}
		for(i = 0; i < n+1; i++) {
			x = (rmax - i*(rmax-rmin)/n)*cos(cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1] + tan(cell->mpara[j][4])*log((rmax - i*(rmax-rmin)/n)/cell->mpara[j][1]));
			y = (rmax - i*(rmax-rmin)/n)*sin(cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1] + tan(cell->mpara[j][4])*log((rmax - i*(rmax-rmin)/n)/cell->mpara[j][1]));
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		}
		
		x = rmin*cos(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1] + tan(cell->mpara[j][4])*log((rmin)/cell->mpara[j][1]));
		y = rmin*sin(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1] + tan(cell->mpara[j][4])*log((rmin)/cell->mpara[j][1]));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		fprintf(wfile, "\n\n");
	}
	
	if(cell->instrutype != NO) {
		printf("instrutype != NO\n");
	}
	
	fclose(wfile);
}

static void aspect_fringe_ex_spiral(char *filename, struct Cell *cell)
{
	int j, i, n;
	double x, y, rmin, rmax;
	//double test;
	FILE *wfile;
	
	n = 20;
	
	rmin = cell->collim.rmin;
	rmax = cell->collim.rmax;
	//test = (rmax - 1/n*(rmax-rmin));
	//printf("test = %lf\n", test);
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		for(i = 0; i < n+1; i++) {
			x = (rmin + i*(rmax-rmin)/n)*cos(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1] + tan(cell->mpara[j][4])*log((rmin + i*(rmax-rmin)/n)/cell->mpara[j][1]));
			y = (rmin + i*(rmax-rmin)/n)*sin(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1] + tan(cell->mpara[j][4])*log((rmin + i*(rmax-rmin)/n)/cell->mpara[j][1]));
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		}
		for(i = 0; i < n+1; i++) {
			x = (rmax - i*(rmax-rmin)/n)*cos(cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1] + tan(cell->mpara[j][4])*log((rmax - i*(rmax-rmin)/n)/cell->mpara[j][1]));
			y = (rmax - i*(rmax-rmin)/n)*sin(cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1] + tan(cell->mpara[j][4])*log((rmax - i*(rmax-rmin)/n)/cell->mpara[j][1]));
			fwk_pt_loctoglob(&x, &y, &(cell->framework));
			fprintf(wfile, "%le  %le\n",x,y);
		}
		
		x = rmin*cos(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1] + tan(cell->mpara[j][4])*log((rmin)/cell->mpara[j][1]));
		y = rmin*sin(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1] + tan(cell->mpara[j][4])*log((rmin)/cell->mpara[j][1]));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		fprintf(wfile, "\n\n");
	}
	
	if(cell->instrutype != NO) {
		printf("instrutype != NO\n");
	}
	
	fclose(wfile);
}

static void aspect_fringe_en_tilt(char *filename, struct Cell *cell)
{
	int j;
	double x, y, rmin, rmax;
	//double test;
	FILE *wfile;
	
	rmin = cell->collim.rmin;
	rmax = cell->collim.rmax;
	
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		//test = (cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin);
		//printf("(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin) = %lf\n", test);
		x = rmin*cos(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin));
		y = rmin*sin(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmax*cos(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmax));
		y = rmax*sin(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmax));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmax*cos(cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmax));
		y = rmax*sin(cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmax));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmin*cos(cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin));
		y = rmin*sin(cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmin*cos(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin));
		y = rmin*sin(cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		fprintf(wfile, "\n\n");
	}
	
	if(cell->instrutype != NO) {
		printf("instrutype != NO\n");
	}
	
	fclose(wfile);
}

static void aspect_fringe_ex_tilt(char *filename, struct Cell *cell)
{
	int j;
	double x, y, rmin, rmax;
	//double test;
	FILE *wfile;
	
	rmin = cell->collim.rmin;
	rmax = cell->collim.rmax;
	
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		//test = (cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin);
		//printf("(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin) = %lf\n", test);
		x = rmin*cos(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin));
		y = rmin*sin(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmax*cos(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmax));
		y = rmax*sin(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmax));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmax*cos(cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmax));
		y = rmax*sin(cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmax));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmin*cos(cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin));
		y = rmin*sin(cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = rmin*cos(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin));
		y = rmin*sin(cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1] + cell->mpara[j][5] - asin(cell->mpara[j][1]*sin(cell->mpara[j][5])/rmin));
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		fprintf(wfile, "\n\n");
	}
	
	if(cell->instrutype != NO) {
		printf("instrutype != NO\n");
	}
	
	fclose(wfile);
}

static void aspect_fringe_en_straighttilt(char *filename, struct Cell *cell)
{
	int j;
	double x, y, xmin, xmax;
	FILE *wfile;
	
	xmin = cell->collim.rmin;
	xmax = cell->collim.rmax;
	
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		x = xmin;
		y = cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmin);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = xmax;
		y = cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmax);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = xmax;
		y = cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmax);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = xmin;
		y = cell->mpara[j][0] + cell->efben[j][0] + cell->efben[j][1] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmin);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = xmin;
		y = cell->mpara[j][0] + cell->efben[j][0] - cell->efben[j][1] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmin);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		fprintf(wfile, "\n\n");
	}
	
	fclose(wfile);
}

static void aspect_fringe_ex_straighttilt(char *filename, struct Cell *cell)
{
	int j;
	double x, y, xmin, xmax;
	FILE *wfile;
	
	xmin = cell->collim.rmin;
	xmax = cell->collim.rmax;
	
	wfile = fopen(filename, "a");
	
	for(j = 0; j < cell->nbcomp; j++) {
		x = xmin;
		y = cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmin);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = xmax;
		y = cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmax);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = xmax;
		y = cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmax);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = xmin;
		y = cell->mpara[j][0] + cell->efbex[j][0] + cell->efbex[j][1] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmin);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		
		x = xmin;
		y = cell->mpara[j][0] + cell->efbex[j][0] - cell->efbex[j][1] - tan(cell->mpara[j][4])*(cell->mpara[j][1]-xmin);
		fwk_pt_loctoglob(&x, &y, &(cell->framework));
		fprintf(wfile, "%le  %le\n",x,y);
		fprintf(wfile, "\n\n");
	}
	
	fclose(wfile);
}

	// ************************************************************************************ //
	//									plot simple ellipse									//
	// ************************************************************************************ //

extern void gnuplot_fit(char *txtfilename1, char *xcolumn1, char *ycolumn1, char *a, char *b, char *c, char *with1, 
				 char *with2, char *xlabel, char *ylabel, char *title1, char *title2, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	char fitname[400], name[500];
	FILE *gp;
	
	sprintf(fitname ,"%s-fit.log",psfilename);
	printf("couc\n");
	remove(fitname);
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,30.'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,30.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -5,0 font 'helvetica,30.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		//fprintf(gp," set terminal postscript color\n");
		fprintf(gp," set terminal postscript eps enhanced color\n");
		fprintf(gp," set output '%s.eps' \n", psfilename);
		fprintf(gp,"set fit logfile '%s-fit.log'\n", psfilename);
		
		//fprintf(gp," f(x) = a*exp(-(x-b)*(x-b)/c/c)\n");
		//fprintf(gp," f(x) = a/2*(tanh((b-x)/c)+tanh((b+x)/c))\n");
		fprintf(gp," f(x) = a/2/c/1.48*(1/(cosh((b+x)/c)*cosh((b+x)/c))-1/(cosh((b-x)/c)*cosh((b-x)/c)))\n");
		
		//fprintf(gp," f(x) = a/pi*(atan((b-x)/c)+atan((b+x)/c))\n");
		//fprintf(gp," f(x) = a/pi/1.48/c*(1/(1+((b+x)/c)*((b+x)/c))-1/(1+((b-x)/c)*((b-x)/c)))\n");
		
		
		
		fprintf(gp,"a = %s\n", a);
		fprintf(gp,"b = %s\n", b);
		fprintf(gp,"c = %s\n", c);
		fprintf(gp,"fit f(x) '%s' using %s:%s via a, b, c\n", txtfilename1, xcolumn1, ycolumn1);
		//fprintf(gp,"fit f(x) '%s' using %s:%s via c\n", txtfilename1, xcolumn1, ycolumn1);
		//fprintf(gp,"fit f(x) '%s' using %s:%s via a, c\n", txtfilename1, xcolumn1, ycolumn1);
		if(title1 != NULL) fprintf(gp,"plot '%s' using %s:%s with %s title '%s', ", txtfilename1, xcolumn1, ycolumn1, with1, title1);
		else fprintf(gp,"plot '%s' using %s:%s with %s, ", txtfilename1, xcolumn1, ycolumn1, with1);
		if(title2 != NULL) fprintf(gp,"f(x) with %s title '%s'\n", with2, title2);
		else fprintf(gp,"f(x) with %s\n", with2);

		//fprintf(gp," set output\n");
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using %s:%s with %s, f(x) with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, with2);
	pclose(gp);
	//sprintf(name, "output/%s.eps", psfilename);
	sprintf(name, "%s.eps", psfilename);
	if(psfilename != NULL) fix_boundingbox(name, YES);
}

extern void gnuplot_fit_general(char *txtfilename, char *xcolumn1, char *ycolumn1, char *fit_function_of_x, char *variables_name, char *init_variables, char *with1, 
				 char *with2, char *xlabel, char *ylabel, char *title1, char *title2, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	char fitname[400], name[500];
	FILE *gp;
	
	sprintf(fitname ,"%s-fit.log", psfilename);
	remove(fitname);
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,30.'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,30.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -5,0 font 'helvetica,30.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		//fprintf(gp," set terminal postscript color\n");
		fprintf(gp," set terminal postscript eps enhanced color\n");
		fprintf(gp," set output '%s.eps' \n", psfilename);
		fprintf(gp,"set fit logfile '%s'\n", fitname);
		fprintf(gp,"f(x) = %s\n", fit_function_of_x);
		fprintf(gp, "%s\n", init_variables);
		fprintf(gp,"fit f(x) '%s' using %s:%s via %s\n", txtfilename, xcolumn1, ycolumn1, variables_name);
		//fprintf(gp,"fit f(x) '%s' using %s:%s via c\n", txtfilename1, xcolumn1, ycolumn1);
		//fprintf(gp,"fit f(x) '%s' using %s:%s via a, c\n", txtfilename1, xcolumn1, ycolumn1);
		if(title1 != NULL) fprintf(gp,"plot '%s' using %s:%s with %s title '%s', ", txtfilename, xcolumn1, ycolumn1, with1, title1);
		else fprintf(gp,"plot '%s' using %s:%s with %s notitle, ", txtfilename, xcolumn1, ycolumn1, with1);
		if(title2 != NULL) fprintf(gp,"f(x) with %s title '%s'\n", with2, title2);
		else fprintf(gp,"f(x) with %s notitle\n", with2);

		//fprintf(gp," set output\n");
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using %s:%s with %s, f(x) with %s\n", txtfilename, xcolumn1, ycolumn1, with1, with2);
	pclose(gp);
	//sprintf(name, "output/%s.eps", psfilename);
	sprintf(name, "%s.eps", psfilename);
	if(psfilename != NULL) fix_boundingbox(name, YES);
}

/*extern void gnuplot_fit2(char *txtfilename1, char *xcolumn1, char *ycolumn1, char *a, char *b, char *c, char *with1, 
				 char *with2, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *setoption)
{
	char fitname[400], name[500];
	FILE *gp;
	
	sprintf(fitname ,"%s-fit.log",psfilename);
	printf("couc\n");
	remove(fitname);
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,30.'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,30.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -5,0 font 'helvetica,30.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	//fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		//fprintf(gp," set terminal postscript color\n");
		fprintf(gp," set terminal postscript eps enhanced color\n");
		fprintf(gp," set output '%s.eps' \n", psfilename);
		fprintf(gp,"set fit logfile '%s-fit.log'\n", psfilename);
		
		//fprintf(gp," f(x) = a*exp(-(x-b)*(x-b)/c/c)\n");
		//fprintf(gp," f(x) = a/2*(tanh((b-x)/c)+tanh((b+x)/c))\n");
		//fprintf(gp," f(x) = a/2/c/1.48*(1/(cosh((b+x)/c)*cosh((b+x)/c))-1/(cosh((b-x)/c)*cosh((b-x)/c)))\n");
		
		fprintf(gp," f(x) = a/pi*(atan((b-x)/c)+atan((b+x)/c))\n");
		
		
		
		fprintf(gp,"a = %s\n", a);
		fprintf(gp,"b = %s\n", b);
		fprintf(gp,"c = %s\n", c);
		fprintf(gp,"fit f(x) '%s' using %s:%s via a, b, c\n", txtfilename1, xcolumn1, ycolumn1);
		//fprintf(gp,"fit f(x) '%s' using %s:%s via c\n", txtfilename1, xcolumn1, ycolumn1);
		//fprintf(gp,"fit f(x) '%s' using %s:%s via a, c\n", txtfilename1, xcolumn1, ycolumn1);
		fprintf(gp,"plot '%s' using %s:%s with %s, f(x) with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, with2);
		//fprintf(gp," set output\n");
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot '%s' using %s:%s with %s, f(x) with %s\n", txtfilename1, xcolumn1, ycolumn1, with1, with2);
	pclose(gp);
	//sprintf(name, "output/%s.eps", psfilename);
	sprintf(name, "%s.eps", psfilename);
	if(psfilename != NULL) fix_boundingbox(name, YES);
}//*/

extern void tune_diag(char *textfilename, int supersym, double qxmin, double qxmax, double qzmin, double qzmax, int maxorder_s, int maxorder_ns, int doyounorm, int doyouskew, char *output_name)
{
	int n, m, p, incr, order, i, npts = 4;
	double pmax, xtics, ytics;
	FILE *fill2norm,*fill4norm,*fill6norm,*fill8norm, *fill10norm,
	*fill1norm,*fill3norm,*fill5norm,*fill7norm, *fill9norm,
	*fill2skew,*fill4skew,*fill6skew,*fill8skew, *fill10skew,
	*fill1skew,*fill3skew,*fill5skew,*fill7skew, *fill9skew;
	
	
	fill1norm = fopen("nonstruct_dip_norm","w");
	fill1skew = fopen("nonstruct_dip_skew","w");
	fill2norm = fopen("struct_dip_norm","w");
	fill2skew = fopen("struct_dip_skew","w");
	fill3norm = fopen("nonstruct_quad_norm","w");
	fill3skew = fopen("nonstruct_quad_skew","w");
	fill4norm = fopen("struct_quad_norm","w");
	fill4skew = fopen("struct_quad_skew","w");
	fill5norm = fopen("nonstruct_sext_norm","w");
	fill5skew = fopen("nonstruct_sext_skew","w");
	fill6norm = fopen("struct_sext_norm","w");
	fill6skew = fopen("struct_sext_skew","w");
	fill7norm = fopen("nonstruct_oct_norm","w");
	fill7skew = fopen("nonstruct_oct_skew","w");
	fill8norm = fopen("struct_oct_norm","w");
	fill8skew = fopen("struct_oct_skew","w");
	fill9norm = fopen("nonstruct_dec_norm","w");
	fill9skew = fopen("nonstruct_dec_skew","w");
	fill10norm = fopen("struct_dec_norm","w");
	fill10skew = fopen("struct_dec_skew","w");
	
	fprintf(fill1norm, "0 0 \n");
	fprintf(fill1skew, "0 0 \n");
	fprintf(fill2norm, "0 0 \n");
	fprintf(fill2skew, "0 0 \n");
	fprintf(fill3norm, "0 0 \n");
	fprintf(fill3skew, "0 0 \n");
	fprintf(fill4norm, "0 0 \n");
	fprintf(fill4skew, "0 0 \n");
	fprintf(fill5norm, "0 0 \n");
	fprintf(fill5skew, "0 0 \n");
	fprintf(fill6norm, "0 0 \n");
	fprintf(fill6skew, "0 0 \n");
	fprintf(fill7norm, "0 0 \n");
	fprintf(fill7skew, "0 0 \n");
	fprintf(fill8norm, "0 0 \n");
	fprintf(fill8skew, "0 0 \n");
	fprintf(fill9norm, "0 0 \n");
	fprintf(fill9skew, "0 0 \n");
	fprintf(fill10norm, "0 0 \n");
	fprintf(fill10skew, "0 0 \n");
	
	pmax = qxmax*maxorder_s;
	if(qzmax > qxmax) pmax = qzmax*maxorder_s;
	pmax = supersym*(1. + floor(pmax/(supersym)));	
	incr = supersym;
	for(n = 0; n <= maxorder_s; n++) {
		for(m = -1*maxorder_s; m <= maxorder_s; m++) {
			for(p = -pmax; p <= pmax; p += incr) {
				order = n + abs(m);
				if(order <= maxorder_s) {
					if(m == 0) {
						if(order == 1) {
							for(i = 0; i < npts; i++) fprintf(fill2norm, "%lf %lf %i %i\n", 1.*p/n, i*(qzmax-qzmin)/(npts - 1) + qzmin, m, n);
							fprintf(fill2norm, "\n");
						}
						if(order == 2) {
							for(i = 0; i < npts; i++) fprintf(fill4norm, "%lf %lf %i %i\n", 1.*p/n, i*(qzmax-qzmin)/(npts - 1) + qzmin, m, n);
							fprintf(fill4norm, "\n");
						}
						if(order == 3) {
							for(i = 0; i < npts; i++) fprintf(fill6norm, "%lf %lf %i %i\n", 1.*p/n, i*(qzmax-qzmin)/(npts - 1) + qzmin, m, n);
							fprintf(fill6norm, "\n");
						}
						if(order == 4) {
							for(i = 0; i < npts; i++) fprintf(fill8norm, "%lf %lf %i %i\n", 1.*p/n, i*(qzmax-qzmin)/(npts - 1) + qzmin, m, n);
							fprintf(fill8norm, "\n");
						}
						if(order == 5) {
							for(i = 0; i < npts; i++) fprintf(fill10norm, "%lf %lf %i %i\n", 1.*p/n, i*(qzmax-qzmin)/(npts - 1) + qzmin, m, n);
							fprintf(fill10norm, "\n");
						}
					}
					else {
						if( (m % 2) == 0) {
							if(order == 1) {
								for(i = 0; i < npts; i++) fprintf(fill2norm, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill2norm, "\n");
							}
							if(order == 2) {
								for(i = 0; i < npts; i++) fprintf(fill4norm, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill4norm, "\n");
							}
							if(order == 3) {
								for(i = 0; i < npts; i++) fprintf(fill6norm, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill6norm, "\n");
							}
							if(order == 4) {
								for(i = 0; i < npts; i++) fprintf(fill8norm, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill8norm, "\n");
							}
							if(order == 5) {
								for(i = 0; i < npts; i++) fprintf(fill10norm, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill10norm, "\n");
							}
						}
						else {
							if(order == 1) {
								for(i = 0; i < npts; i++) fprintf(fill2skew, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill2skew, "\n");
							}
							if(order == 2) {
								for(i = 0; i < npts; i++) fprintf(fill4skew, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill4skew, "\n");
							}
							if(order == 3) {
								for(i = 0; i < npts; i++) fprintf(fill6skew, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill6skew, "\n");
							}
							if(order == 4) {
								for(i = 0; i < npts; i++) fprintf(fill8skew, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill8skew, "\n");
							}
							if(order == 5) {
								for(i = 0; i < npts; i++) fprintf(fill10skew, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill10skew, "\n");
							}
						}
					}
				}
			}
		}
	}
	
	pmax = qxmax*maxorder_ns;
	if(qzmax > qxmax) pmax = qzmax*maxorder_ns;
	incr = 1;	
	for(n = 0; n <= maxorder_ns; n++) {
		for(m = -1*maxorder_ns; m <= maxorder_ns; m++) {
			for(p = -pmax; p <= pmax; p += incr) {
				order = n + abs(m);
				if(order <= maxorder_ns) {
					if(m == 0) {
						if(order == 1) {
							for(i = 0; i < npts; i++) fprintf(fill1norm, "%lf %lf %i %i\n", 1.*p/n, i*(qzmax-qzmin)/(npts - 1) + qzmin, m, n);
							fprintf(fill1norm, "\n");
						}
						if(order == 2) {
							for(i = 0; i < npts; i++) fprintf(fill3norm, "%lf %lf %i %i\n", 1.*p/n, i*(qzmax-qzmin)/(npts - 1) + qzmin, m, n);
							fprintf(fill3norm, "\n");
						}
						if(order == 3) {
							for(i = 0; i < npts; i++) fprintf(fill5norm, "%lf %lf %i %i\n", 1.*p/n, i*(qzmax-qzmin)/(npts - 1) + qzmin, m, n);
							fprintf(fill5norm, "\n");
						}
						if(order == 4) {
							for(i = 0; i < npts; i++) fprintf(fill7norm, "%lf %lf %i %i\n", 1.*p/n, i*(qzmax-qzmin)/(npts - 1) + qzmin, m, n);
							fprintf(fill7norm, "\n");
						}
						if(order == 5) {
							for(i = 0; i < npts; i++) fprintf(fill9norm, "%lf %lf %i %i\n", 1.*p/n, i*(qzmax-qzmin)/(npts - 1) + qzmin, m, n);
							fprintf(fill9norm, "\n");
						}
					}
					
					else {
						if( (m % 2) == 0) {
							if(order == 1) {
								for(i = 0; i < npts; i++) fprintf(fill1norm, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill1norm, "\n");
							}
							if(order == 2) {
								for(i = 0; i < npts; i++) fprintf(fill3norm, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill3norm, "\n");
							}
							if(order == 3) {
								for(i = 0; i < npts; i++) fprintf(fill5norm, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill5norm, "\n");
							}
							if(order == 4) {
								for(i = 0; i < npts; i++) fprintf(fill7norm, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill7norm, "\n");
							}
							if(order == 5) {
								for(i = 0; i < npts; i++) fprintf(fill9norm, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill9norm, "\n");
							}
						}
						else {
							if(order == 1) {
								for(i = 0; i < npts; i++) fprintf(fill1skew, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill1skew, "\n");
							}
							if(order == 2) {
								for(i = 0; i < npts; i++) fprintf(fill3skew, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill3skew, "\n");
							}
							if(order == 3) {
								for(i = 0; i < npts; i++) fprintf(fill5skew, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill5skew, "\n");
							}
							if(order == 4) {
								for(i = 0; i < npts; i++) fprintf(fill7skew, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill7skew, "\n");
							}
							if(order == 5) {
								for(i = 0; i < npts; i++) fprintf(fill9skew, "%lf %lf %i %i\n", qxmin + 1.*i*(qxmax-qxmin)/(npts - 1), -1.*n/m*(qxmin + 1.*i*(qxmax-qxmin)/(npts - 1)) + 1.*p/m, m, n);
								fprintf(fill9skew, "\n");
							}
						}
					}
				}
			}
		}
	}
	
	fclose(fill1norm);
	fclose(fill2norm);
	fclose(fill3norm);
	fclose(fill4norm);
	fclose(fill5norm);
	fclose(fill6norm);
	fclose(fill7norm);
	fclose(fill8norm);
	fclose(fill9norm);
	fclose(fill10norm);
	fclose(fill1skew);
	fclose(fill2skew);
	fclose(fill3skew);
	fclose(fill4skew);
	fclose(fill5skew);
	fclose(fill6skew);
	fclose(fill7skew);
	fclose(fill8skew);
	fclose(fill9skew);
	fclose(fill10skew);
	
	xtics = (qxmax - qxmin)/4.;
	ytics = (qzmax - qzmin)/4.;
	//Gnuplot
	FILE *gp;
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//fprintf(gp, "set colorsequence classic\n");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp, "set term x11\n");
	fprintf(gp," set xrange [%f:%f]\n",qxmin,qxmax);
	fprintf(gp," set yrange [%f:%f]\n",qzmin,qzmax);
	//fprintf(gp," set tics font 'Arial,15.'\n");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp," set xlabel 'Qx' offset 0,-1 font 'helvetica,25.' \n");
	fprintf(gp," set ylabel 'Qz' offset -3,0 font 'helvetica,25.'\n");
	//fprintf(gp,"unset tics\n");
	fprintf(gp,"unset key\n");
	fprintf(gp," set multiplot\n");
	fprintf(gp," set sample 500\n");
	fprintf(gp," set size ratio -1\n");
	if(doyounorm == YES) {
		//fprintf(gp,"plot 'nonstruct_dec_norm' using 1:2 with lines lt 2 lw 1 lc 5, 'struct_dec_norm' using 1:2 with lines lt 1 lw 8 lc 5\n");
		//fprintf(gp,"plot 'nonstruct_oct_norm' using 1:2 with lines lt 2 lw 1 lc 4, 'struct_oct_norm' using 1:2 with lines lt 1 lw 8 lc 4\n");
		//fprintf(gp,"plot 'nonstruct_sext_norm' using 1:2 with lines lt 2 lw 1 lc 3, 'struct_sext_norm' using 1:2 with lines lt 1 lw 8 lc 3\n");
		//fprintf(gp,"plot 'nonstruct_quad_norm' using 1:2 with lines lt 2 lw 1 lc 2, 'struct_quad_norm' using 1:2 with lines lt 1 lw 8 lc 2\n");
		//fprintf(gp,"plot 'nonstruct_dip_norm' using 1:2 with lines lt 2 lw 1 lc 1, 'struct_dip_norm' using 1:2 with lines lt 1 lw 8 lc 1\n");
		fprintf(gp,"plot 'nonstruct_dec_norm' using 1:2 with lines dt 1 lw 1 lc 3, 'struct_dec_norm' using 1:2 with lines dt 1 lw 8 lc 3\n");
		fprintf(gp,"plot 'nonstruct_oct_norm' using 1:2 with lines dt 1 lw 1 lc 1, 'struct_oct_norm' using 1:2 with lines dt 1 lw 8 lc 1\n");
		fprintf(gp,"plot 'nonstruct_sext_norm' using 1:2 with lines dt 1 lw 1 lc 6, 'struct_sext_norm' using 1:2 with lines dt 1 lw 8 lc 6\n");
		fprintf(gp,"plot 'nonstruct_quad_norm' using 1:2 with lines dt 1 lw 1 lc 2, 'struct_quad_norm' using 1:2 with lines dt 1 lw 8 lc 2\n");
		fprintf(gp,"plot 'nonstruct_dip_norm' using 1:2 with lines dt 1 lw 1 lc 7, 'struct_dip_norm' using 1:2 with lines dt 1 lw 8 lc 7\n");
	}
	if(doyouskew == YES) {
		//fprintf(gp,"plot 'nonstruct_dec_skew' using 1:2 with lines lt 2 lw 1 lc 5, 'struct_dec_skew' using 1:2 with lines lt 1 lw 4 lc 5\n");
		//fprintf(gp,"plot 'nonstruct_oct_skew' using 1:2 with lines lt 2 lw 1 lc 4, 'struct_oct_skew' using 1:2 with lines lt 1 lw 4 lc 4\n");
		//fprintf(gp,"plot 'nonstruct_sext_skew' using 1:2 with lines lt 2 lw 1 lc 3, 'struct_sext_skew' using 1:2 with lines lt 1 lw 4 lc 3\n");
		//fprintf(gp,"plot 'nonstruct_quad_skew' using 1:2 with lines lt 2 lw 1 lc 2, 'struct_quad_skew' using 1:2 with lines lt 1 lw 4 lc 2\n");
		//fprintf(gp,"plot 'nonstruct_dip_skew' using 1:2 with lines lt 2 lw 1 lc 1, 'struct_dip_skew' using 1:2 with lines lt 1 lw 4 lc 1\n");
		fprintf(gp,"plot 'nonstruct_dec_skew' using 1:2 with lines dt 2 lw 1 lc 3, 'struct_dec_skew' using 1:2 with lines dt 2 lw 4 lc 3\n");
		fprintf(gp,"plot 'nonstruct_oct_skew' using 1:2 with lines dt 2 lw 1 lc 1, 'struct_oct_skew' using 1:2 with lines dt 2 lw 4 lc 1\n");
		fprintf(gp,"plot 'nonstruct_sext_skew' using 1:2 with lines dt 2 lw 1 lc 6, 'struct_sext_skew' using 1:2 with lines dt 2 lw 4 lc 6\n");
		fprintf(gp,"plot 'nonstruct_quad_skew' using 1:2 with lines dt 2 lw 1 lc 2, 'struct_quad_skew' using 1:2 with lines dt 2 lw 4 lc 2\n");
		fprintf(gp,"plot 'nonstruct_dip_skew' using 1:2 with lines dt 2 lw 1 lc 7, 'struct_dip_skew' using 1:2 with lines dt 2 lw 4 lc 7\n");
	}
	//fprintf(gp,"plot 'whitesquare.dat' using 1:2  with points pt 4 ps 17.3 lc 0 lw 3\n");
	//fprintf(gp,"plot 'whitesquare.dat' using 1:2  with points pt 5 ps 17 lc rgb '#FFFFFF'\n");
	//fprintf(gp,"plot 'tunepoints_scan.dat' using 1:2:(sqrt($3)/10)  with points pt 7 pointsize variable\n");
	//fprintf(gp,"plot 'tunepoints_scan.dat' using 1:2:(sqrt($3)/10)  with points pt 6 pointsize variable lc 0\n");
	//fprintf(gp,"plot 'tunepoints.dat' using 1:2  with points lt 2 pt 7 ps 2 linecolor 0\n");
	//if(textfilename!=NULL) fprintf(gp,"plot '%s' using 2:3:($4/10000)  with points pt 2\n", textfilename);
	//if(textfilename!=NULL) fprintf(gp,"plot '%s' using 4:5:($4/10000)  with points pt 2\n", textfilename);	
	if(textfilename!= NULL) fprintf(gp,"plot '%s' using 2:3  with linespoints lt 1 pt 7 ps 1 lw 4 linecolor 0\n",textfilename);
	if(textfilename!= NULL) fprintf(gp,"plot '%s' using 4:5  with linespoints lt 1 pt 7 ps 1 lw 4 linecolor 0\n",textfilename);
	
	//fprintf(gp,"plot '%s' using 6:7  with points pt 5 lc 0\n", "data/vffa_err_study/1sttry_shift0.000010_twist0.000000_results.dat");	
	//fprintf(gp,"plot '%s' using 6:7  with points pt 2 lc 7\n", "data/vffa_err_study/1sttry_shift0.000000_twist0.000000_results.dat");	
	//fprintf(gp,"plot '%s' using 2:3:($4/10000)  with points pt 2\n", textfilename);
	//fprintf(gp,"plot '%s' using 4:5:($4/10000)  with points pt 2\n", textfilename);
	
	pclose(gp);
//	
	//EPS save
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//fprintf(gp, "set colorsequence classic\n");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set xrange [%f:%f]\n",qxmin,qxmax);
	fprintf(gp," set yrange [%f:%f]\n",qzmin,qzmax);
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set tics font 'helvetica,25.'\n");
	fprintf(gp,"set xtics %f offset 0,0\n set mxtics 5\n", xtics);
	fprintf(gp,"set ytics %f offset 0,0\n set mytics 5\n", ytics);
	//fprintf(gp," set xlabel 'Q_u' offset 0,-1 font 'helvetica,25.' \n");
	//fprintf(gp," set ylabel 'Q_v' offset -4,0 font 'helvetica,25.'\n");
	fprintf(gp," set xlabel 'Q_x' offset 0,-1 font 'helvetica,25.' \n");
	fprintf(gp," set ylabel 'Q_z' offset -4,0 font 'helvetica,25.'\n");
	//fprintf(gp,"unset tics\n");
	fprintf(gp,"unset key\n");
	fprintf(gp," set terminal postscript eps enhanced color\n");
	fprintf(gp," set output '%s' \n", output_name);
	fprintf(gp," set multiplot\n");
	fprintf(gp," set sample 500\n");
	fprintf(gp," set size ratio -1\n");
	if(doyouskew == YES) {
		fprintf(gp,"plot 'nonstruct_dec_skew'  using 1:2 with lines dt 2 lw 1 lc 3, 'struct_dec_skew'  using 1:2 with lines dt 2 lw 4 lc 3\n");
		fprintf(gp,"plot 'nonstruct_oct_skew'  using 1:2 with lines dt 2 lw 1 lc 1, 'struct_oct_skew'  using 1:2 with lines dt 2 lw 4 lc 1\n");
		fprintf(gp,"plot 'nonstruct_sext_skew' using 1:2 with lines dt 2 lw 1 lc 6, 'struct_sext_skew' using 1:2 with lines dt 2 lw 4 lc 6\n");
		fprintf(gp,"plot 'nonstruct_quad_skew' using 1:2 with lines dt 2 lw 1 lc 2, 'struct_quad_skew' using 1:2 with lines dt 2 lw 4 lc 2\n");
		fprintf(gp,"plot 'nonstruct_dip_skew'  using 1:2 with lines dt 2 lw 1 lc 7, 'struct_dip_skew'  using 1:2 with lines dt 2 lw 4 lc 7\n");
		//fprintf(gp,"plot 'nonstruct_dec_skew' using 1:2 with lines lt 2 lw 0.5 lc 5, 'struct_dec_skew' using 1:2 with lines lt 1 lw 8 lc 5\n");
		//fprintf(gp,"plot 'nonstruct_oct_skew' using 1:2 with lines lt 2 lw 0.5 lc 4, 'struct_oct_skew' using 1:2 with lines lt 1 lw 8 lc 4\n");
		//fprintf(gp,"plot 'nonstruct_sext_skew' using 1:2 with lines lt 2 lw 0.5 lc 3, 'struct_sext_skew' using 1:2 with lines lt 1 lw 8 lc 3\n");
		//fprintf(gp,"plot 'nonstruct_quad_skew' using 1:2 with lines lt 2 lw 0.5 lc 2, 'struct_quad_skew' using 1:2 with lines lt 1 lw 8 lc 2\n");
		//fprintf(gp,"plot 'nonstruct_dip_skew' using 1:2 with lines lt 2 lw 0.5 lc 1, 'struct_dip_skew' using 1:2 with lines lt 1 lw 8 lc 1\n");
		//fprintf(gp,"plot 'nonstruct_dec_skew'  using 1:2 with lines lt 5 dt 2 lw 0.5, 'struct_dec_skew'  using 1:2 with lines lt 5 dt 2 lw 8\n");
		//fprintf(gp,"plot 'nonstruct_oct_skew'  using 1:2 with lines lt 4 dt 2 lw 0.5, 'struct_oct_skew'  using 1:2 with lines lt 4 dt 2 lw 8\n");
		//fprintf(gp,"plot 'nonstruct_sext_skew' using 1:2 with lines lt 3 dt 2 lw 0.5, 'struct_sext_skew' using 1:2 with lines lt 3 dt 2 lw 8\n");
		//fprintf(gp,"plot 'nonstruct_quad_skew' using 1:2 with lines lt 2 dt 2 lw 0.5, 'struct_quad_skew' using 1:2 with lines lt 2 dt 2 lw 8\n");
		//fprintf(gp,"plot 'nonstruct_dip_skew'  using 1:2 with lines lt 1 dt 2 lw 0.5, 'struct_dip_skew'  using 1:2 with lines lt 1 dt 2 lw 8\n");
	}
	if(doyounorm == YES) {
		fprintf(gp,"plot 'nonstruct_dec_norm'  using 1:2 with lines dt 1 lw 1 lc 3, 'struct_dec_norm'  using 1:2 with lines dt 1 lw 8 lc 3\n");
		fprintf(gp,"plot 'nonstruct_oct_norm'  using 1:2 with lines dt 1 lw 1 lc 1, 'struct_oct_norm'  using 1:2 with lines dt 1 lw 8 lc 1\n");
		fprintf(gp,"plot 'nonstruct_sext_norm' using 1:2 with lines dt 1 lw 1 lc 6, 'struct_sext_norm' using 1:2 with lines dt 1 lw 8 lc 6\n");
		fprintf(gp,"plot 'nonstruct_quad_norm' using 1:2 with lines dt 1 lw 1 lc 2, 'struct_quad_norm' using 1:2 with lines dt 1 lw 8 lc 2\n");
		fprintf(gp,"plot 'nonstruct_dip_norm'  using 1:2 with lines dt 1 lw 1 lc 7, 'struct_dip_norm'  using 1:2 with lines dt 1 lw 8 lc 7\n");
		//fprintf(gp,"plot 'nonstruct_dec_norm' using 1:2 with lines lt 2 lw 1 lc 5\n");
		//fprintf(gp,"plot 'struct_dec_norm' using 1:2 with lines lt 1 lw 8 lc 5\n");
		//fprintf(gp,"plot 'nonstruct_oct_norm' using 1:2 with lines lt 2 lw 1 lc 4\n");
		//fprintf(gp,"plot 'struct_oct_norm' using 1:2 with lines lt 1 lw 8 lc 4\n");
		//fprintf(gp,"plot 'nonstruct_sext_norm' using 1:2 with lines lt 2 lw 1 lc 3\n");
		//fprintf(gp,"plot 'struct_sext_norm' using 1:2 with lines lt 1 lw 8 lc 3\n");
		//fprintf(gp,"plot 'nonstruct_quad_norm' using 1:2 with lines lt 2 lw 1 lc 2\n");
		//fprintf(gp,"plot 'struct_quad_norm' using 1:2 with lines lt 1 lw 8 lc 2\n");
		//fprintf(gp,"plot 'nonstruct_dip_norm' using 1:2 with lines lt 2 lw 1 lc 1\n");
		//fprintf(gp,"plot 'struct_dip_norm' using 1:2 with lines lt 1 lw 8 lc 1\n");
		//fprintf(gp,"plot 'nonstruct_dec_norm'  using 1:2 with lines lt 5 dt 1 lw 1, 'struct_dec_norm'  using 1:2 with lines lt 5 dt 1 lw 8\n");
		//fprintf(gp,"plot 'nonstruct_oct_norm'  using 1:2 with lines lt 4 dt 1 lw 1, 'struct_oct_norm'  using 1:2 with lines lt 4 dt 1 lw 8\n");
		//fprintf(gp,"plot 'nonstruct_sext_norm' using 1:2 with lines lt 3 dt 1 lw 1, 'struct_sext_norm' using 1:2 with lines lt 3 dt 1 lw 8\n");
		//fprintf(gp,"plot 'nonstruct_quad_norm' using 1:2 with lines lt 2 dt 1 lw 1, 'struct_quad_norm' using 1:2 with lines lt 2 dt 1 lw 8\n");
		//fprintf(gp,"plot 'nonstruct_dip_norm'  using 1:2 with lines lt 1 dt 1 lw 1, 'struct_dip_norm'  using 1:2 with lines lt 1 dt 1 lw 8\n");
		//fprintf(gp,"set xtics 0.2\n");
		//fprintf(gp,"set ytics 0.2\n");
	}
	
	//if(textfilename!= NULL) fprintf(gp,"plot 'tunepoints.dat' using 2:3  with linespoints lt 1 pt 7 ps 1 lw 4 linecolor 0\n");
	//if(textfilename!= NULL) fprintf(gp,"plot 'tunepoints.dat' using 4:5  with linespoints lt 1 pt 7 ps 1 lw 4 linecolor 0\n");
	if(textfilename!= NULL) fprintf(gp,"plot '%s' using 2:3  with linespoints lt 1 pt 7 ps 1 lw 4 linecolor 0\n",textfilename);
	if(textfilename!= NULL) fprintf(gp,"plot '%s' using 4:5  with linespoints lt 1 pt 7 ps 1 lw 4 linecolor 0\n",textfilename);
	//fprintf(gp,"plot 'tunevar+.dat' using 1:2  with linespoints lt 1 pt 7 ps 1 lw 4 linecolor 0\n");
	//if(textfilename!= NULL) fprintf(gp,"plot '%s' using 4:5  with linespoints dt 2 lw 2 pt 2 ps 1 linecolor 0\n",textfilename);
	//if(textfilename!= NULL) fprintf(gp,"plot '%s' using 2:3  with linespoints dt 2 lw 2 pt 2 ps 1 linecolor 0\n",textfilename);
	pclose(gp);
	
	if(output_name != NULL) fix_boundingbox(output_name, YES);
	remove("nonstruct_dec_skew");
	remove("nonstruct_dec_norm");
	remove("nonstruct_oct_skew");
	remove("nonstruct_oct_norm");
	remove("nonstruct_sext_skew");
	remove("nonstruct_sext_norm");
	remove("nonstruct_quad_skew");
	remove("nonstruct_quad_norm");
	remove("nonstruct_dip_skew");
	remove("nonstruct_dip_norm");
	remove("struct_dec_skew");
	remove("struct_dec_norm");
	remove("struct_oct_skew");
	remove("struct_oct_norm");
	remove("struct_sext_skew");
	remove("struct_sext_norm");
	remove("struct_quad_skew");
	remove("struct_quad_norm");
	remove("struct_dip_skew");
	remove("struct_dip_norm");
}

extern void plotfunction(char *function, char *xlabel, char *ylabel, char *xrange, char *yrange, char *psfilename, char *title, char *setoption)
{
	FILE *gp;
	
	//gp=popen("/Applications/gnuplot.app/gnuplot.command -persist","w");
	//gp=popen("/usr/local/bin/gnuplot -persist","w");
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,30.'\n");
	fprintf(gp,"set key font 'helvetica,20'\n");
	fprintf(gp,"set bmargin 5\n");
	fprintf(gp,"set rmargin 5\n");
	fprintf(gp," set xlabel \"%s\" offset 0,-1 font 'helvetica,30.' \n", xlabel);
	fprintf(gp," set ylabel \"%s\" offset -5,0 font 'helvetica,30.'\n", ylabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	if(title == NULL) fprintf(gp,"unset key\n");
	if(psfilename != NULL) {
		//fprintf(gp," set terminal postscript color\n");
		fprintf(gp," set terminal postscript eps enhanced color\n");
		fprintf(gp," set output '%s' \n", psfilename);
		if(title == NULL) fprintf(gp,"plot %s\n", function);
		else fprintf(gp,"plot %s title '%s'\n", function, title);
		//fprintf(gp,"plot (-x-(0.33*0.006-x*x)**(0.5))/0.33\n");
		//fprintf(gp," set output\n");
		fprintf(gp," set terminal x11\n");
	}
	fprintf(gp,"plot %s\n", function);
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
}

extern void gnuplot_3d_scatered_points(char *txtfilename, char *xcolumn, char *ycolumn, char *cbcolumn, char *xlabel, char *ylabel, char *cblabel, char *xrange, char *yrange, char *psfilename, char *setoption, double dens_max)
{
	FILE *gp;
	
	gp=popen("gnuplot -persist","w");
	fprintf(gp," set tics font 'helvetica,20.'\n");
	fprintf(gp,"set pm3d map\n");
	
	fprintf(gp,"set palette defined (0 'white', %le 'cyan', %le 'green', %le 'yellow', %le 'orange', %le 'red')\n", dens_max/5., dens_max*2/5., dens_max*3/5.,dens_max*4/5, dens_max);
	fprintf(gp, " set cbrange [0:%le]\n", dens_max);
	fprintf(gp," set xlabel '%s' offset 0,-1 font 'helvetica,25.' \n", xlabel);
	fprintf(gp," set ylabel '%s' offset -3,0 font 'helvetica,25.'\n", ylabel);
	fprintf(gp," set cblabel '%s' offset 3,0 font 'helvetica,25.'\n",cblabel);
	if(xrange != NULL) fprintf(gp," set xrange %s \n", xrange);
	if(yrange != NULL) fprintf(gp," set yrange %s \n", yrange);
	if(setoption != NULL) fprintf(gp," set %s \n", setoption);
	fprintf(gp,"unset key\n");
	fprintf(gp," set terminal postscript eps enhanced color \n");
	fprintf(gp," set output '%s' \n", psfilename);
	fprintf(gp,"splot '%s' u %s:%s:%s with points pt 7 ps 3 linecolor palette\n", txtfilename, xcolumn, ycolumn, cbcolumn);
	fprintf(gp," set output\n");
	fprintf(gp," set terminal x11\n");
	fprintf(gp,"splot '%s' u %s:%s:%s with points pt 7 ps 3 linecolor palette\n", txtfilename, xcolumn, ycolumn, cbcolumn);
	
	pclose(gp);
	if(psfilename != NULL) fix_boundingbox(psfilename, YES);
	
}
