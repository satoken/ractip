/*
 * $Id$
 * 
 * Copyright (C) 2010 Kengo Sato
 *
 * This file is part of RactIP.
 *
 * RactIP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * RactIP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with RactIP.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#undef HAVE_CONFIG_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>

#define MIN2(A, B)      ((A) < (B) ? (A) : (B))
#define MAX2(A, B)      ((A) > (B) ? (A) : (B))
void  *space(unsigned size) /*@ensures MaxSet(result) == (size-1);@*/;
void  *xrealloc(/*@null@*/ /*@only@*/ /*@out@*/ /*@returned@*/ void *p, unsigned size) /*@modifies *p @*/ /*@ensures MaxSet(result) == (size-1) @*/;

#include <ViennaRNA/energy_par.h>
#include <ViennaRNA/fold_vars.h>
#include <ViennaRNA/fold.h>
#include <ViennaRNA/pair_mat.h>
#include <ViennaRNA/params.h>
#ifdef HAVE_VIENNA20
#include <ViennaRNA/loop_energies.h>
#endif
#include "pf_duplex.h"

#define PUBLIC
#define PRIVATE static

PRIVATE void  encode_seqs(const char *s1, const char *s2);
PRIVATE short *encode_seq(const char *seq);

PRIVATE paramT *P = NULL;

PRIVATE double **fw;      /* energy array, given that i-j pair */
PRIVATE double **bk;      /* energy array, given that i-j pair */
PUBLIC double **pr_duplex;      /* energy array, given that i-j pair */
PUBLIC double **pr_duplex2;
PRIVATE short  *S1, *SS1, *S2, *SS2;
PRIVATE int   n1,n2;    /* sequence lengths */
#ifndef HAVE_VIENNA20
extern  int  LoopEnergy(int n1, int n2, int type, int type_2,
                        int si1, int sj1, int sp1, int sq1);
#endif
PRIVATE double pf_duplex_fw();
PRIVATE double pf_duplex_bk();

inline
double
LogAdd(double x, double y)
{
  if (x<=log(0.0)) return y;
  if (y<=log(0.0)) return x;
  return x>y ? log1p(exp(y-x))+x : log1p(exp(x-y))+y;
}

double
pf_duplex(const char *s1, const char *s2)
{
  int i, j;
  double Esum;
  double kT = (temperature+K0)*GASCONST;
  n1 = (int) strlen(s1);
  n2 = (int) strlen(s2);
  
  if ((!P) || (fabs(P->temperature - temperature)>1e-6)) {
    update_fold_params();  P = scale_parameters();
    make_pair_matrix();
  }

  encode_seqs(s1, s2);

  Esum = pf_duplex_fw();
  pf_duplex_bk();

  pr_duplex = (double **) space(sizeof(double *) * (n1+1));
  for (i=1; i<=n1; i++) pr_duplex[i] = (double *) space(sizeof(double) * (n2+1));
  pr_duplex2 = (double **) space(sizeof(double *) * (n1+1));
  for (i=1; i<=n1; i++) pr_duplex2[i] = (double *) space(sizeof(double) * (n2+1));

  for (i=1; i<=n1; i++) {
    for (j=n2; j>0; j--) {
      int type, type2, E;
      type = pair[S1[i]][S2[j]];
      pr_duplex[i][j] = exp(fw[i][j]+bk[i][j]-Esum);
      if (i>1 && j<n2) {
        type2 = pair[S1[i-1]][S2[j+1]];
        if (type2) {
#ifdef HAVE_VIENNA20
          E = E_IntLoop(0, 0, type2, rtype[type],
                        SS1[i], SS2[j], SS1[i-1], SS2[j+1], P);
#else
	  E = LoopEnergy(0, 0, type2, rtype[type],
                         SS1[i], SS2[j], SS1[i-1], SS2[j+1]);
#endif
          pr_duplex2[i][j] = exp(fw[i-1][j+1]+bk[i][j]-E*10./kT-Esum);
        }        
      }
    }
  }

  for (i=1; i<=n1; i++) {
    free(fw[i]);
    free(bk[i]);
  }
  free(fw);
  free(bk);
  
  return Esum;
}

void
free_pf_duplex()
{
  int i;
  for (i=1; i<=n1; i++) free(pr_duplex[i]);
  free(pr_duplex);
  for (i=1; i<=n1; i++) free(pr_duplex2[i]);
  free(pr_duplex2);
}

PRIVATE
double
pf_duplex_fw()
{
  int i, j;
  double Esum=log(0.0);
  double kT = (temperature+K0)*GASCONST;
  
  fw = (double **) space(sizeof(double *) * (n1+1));
  for (i=1; i<=n1; i++) fw[i] = (double *) space(sizeof(double) * (n2+1));
  
  for (i=1; i<=n1; i++) {
    for (j=n2; j>0; j--) {
      int type, type2, E, k,l;
      type = pair[S1[i]][S2[j]];
      fw[i][j] = log(0.0);
      if (!type) continue;
      E = P->DuplexInit;
      if (i>1)  E += P->dangle5[type][SS1[i-1]];
      if (j<n2) E += P->dangle3[type][SS2[j+1]];
      if (type>2) E += P->TerminalAU;
      fw[i][j] = -E*10./kT;
      for (k=i-1; k>0 && k>i-MAXLOOP-2; k--) {
	for (l=j+1; l<=n2; l++) {
	  if (i-k+l-j-2>MAXLOOP) break;
	  type2 = pair[S1[k]][S2[l]];
	  if (!type2) continue;
#ifdef HAVE_VIENNA20
	  E = E_IntLoop(i-k-1, l-j-1, type2, rtype[type],
                        SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1], P);
#else
	  E = LoopEnergy(i-k-1, l-j-1, type2, rtype[type],
			    SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1]);
#endif
	  fw[i][j] = LogAdd(fw[i][j], fw[k][l]-E*10./kT);
	}
      }
      E = 0.;
      if (i<n1) E += P->dangle3[rtype[type]][SS1[i+1]];
      if (j>1)  E += P->dangle5[rtype[type]][SS2[j-1]];
      if (type>2) E += P->TerminalAU;
      Esum = LogAdd(Esum, fw[i][j]-E*10./kT);
    }
  }

  return Esum;
}

PRIVATE
double
pf_duplex_bk()
{
  int i, j;
  double Esum=log(0.0);
  double kT = (temperature+K0)*GASCONST;
  
  bk = (double **) space(sizeof(double *) * (n1+1));
  for (i=1; i<=n1; i++) {
    bk[i] = (double *) space(sizeof(double) * (n2+1));
    for (j=1; j<=n2; j++) bk[i][j] = log(0.0);
  }
  
  for (i=n1; i>0; i--) {
    for (j=1; j<=n2; j++) {
      int type, type2, E, k,l;
      type = pair[S1[i]][S2[j]];
      if (!type) continue;
      E = 0.;
      if (i<n1) E += P->dangle3[rtype[type]][SS1[i+1]];
      if (j>1)  E += P->dangle5[rtype[type]][SS2[j-1]];
      if (type>2) E += P->TerminalAU;
      bk[i][j] = LogAdd(bk[i][j], -E*10./kT);

      for (k=i-1; k>0 && k>i-MAXLOOP-2; k--) {
	for (l=j+1; l<=n2; l++) {
	  if (i-k+l-j-2>MAXLOOP) break;
	  type2 = pair[S1[k]][S2[l]];
	  if (!type2) continue;
#ifdef HAVE_VIENNA20
	  E = E_IntLoop(i-k-1, l-j-1, type2, rtype[type],
                        SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1], P);
#else
	  E = LoopEnergy(i-k-1, l-j-1, type2, rtype[type],
			    SS1[k+1], SS2[l-1], SS1[i-1], SS2[j+1]);
#endif
          bk[k][l] = LogAdd(bk[k][l], bk[i][j]-E*10./kT);
	}
      }

      E = P->DuplexInit;
      if (i>1)  E += P->dangle5[type][SS1[i-1]];
      if (j<n2) E += P->dangle3[type][SS2[j+1]];
      if (type>2) E += P->TerminalAU;
      Esum = LogAdd(Esum, bk[i][j]-E*10./kT);
    }
  }

  return Esum;
}

PRIVATE short * encode_seq(const char *sequence) {
  unsigned int i,l;
  short *S;
  l = strlen(sequence);
  S = (short *) space(sizeof(short)*(l+2));
  S[0] = (short) l;

  /* make numerical encoding of sequence */
  for (i=1; i<=l; i++)
    S[i]= (short) encode_char(toupper(sequence[i-1]));

  /* for circular folding add first base at position n+1 */
  S[l+1] = S[1];

  return S;
}

PRIVATE void encode_seqs(const char *s1, const char *s2) {
  unsigned int i,l;

  l = strlen(s1);
  S1 = encode_seq(s1);
  SS1= (short *) space(sizeof(short)*(l+1));
  /* SS1 exists only for the special X K and I bases and energy_set!=0 */
  
  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    SS1[i] = alias[S1[i]];   /* for mismatches of nostandard bases */
  }

  l = strlen(s2);
  S2 = encode_seq(s2);
  SS2= (short *) space(sizeof(short)*(l+1));
  /* SS2 exists only for the special X K and I bases and energy_set!=0 */
  
  for (i=1; i<=l; i++) { /* make numerical encoding of sequence */
    SS2[i] = alias[S2[i]];   /* for mismatches of nostandard bases */
  }
}
