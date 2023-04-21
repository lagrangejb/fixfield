/*
 *  init.h
 *  ringdesign
 *
 *  Copyright 2009 Kyoto University. All rights reserved.
 *
 */

#ifndef INIT 
#define INIT

#include "constants.h"
#include "common.h"

static void genecell_collimator(FILE *rfile, int i, struct Lattice *latt);
static void genecell_float_faradaycup(FILE *rfile, int i, struct Lattice *latt);
static void genecell_drift(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_thindipole(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_realdipole(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_puredipole(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_ffagr(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_ffagstr(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_ffagspiral(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_cavity(FILE *rfile, int i, int periodicity, struct Lattice *latt);
static void load_map(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void load_map_legacy(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_ffagtilt(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_ffagstrtilt(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_quad(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_sext(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_oct(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_vffa_rect(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_vffa_rect_str(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void genecell_vffa_sect(FILE *rfile, struct Framework *fwk, int i, struct Lattice *latt);
static void duplicate_cell_fmod(int nb_refcell, struct Framework *fwk, int nb_idcells, struct Lattice *latt);



#endif
