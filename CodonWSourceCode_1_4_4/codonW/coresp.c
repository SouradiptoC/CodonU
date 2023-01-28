/**************************************************************************/
/* CodonW codon usage analysis package                                    */
/* Copyright (C) 2005            John F. Peden                            */
/* This program is free software; you can redistribute                    */
/* it and/or modify it under the terms of the GNU General Public License  */
/* as published by the Free Software Foundation; version 2 of the         */
/* License,                                                               */
/*                                                                        */
/* This program is distributed in the hope that it will be useful, but    */
/* WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the           */
/* GNU General Public License for more details.                           */
/* You should have received a copy of the GNU General Public License along*/
/* with this program; if not, write to the Free Software Foundation, Inc.,*/
/* 675 Mass Ave, Cambridge, MA 02139, USA.                                */
/*                                                                        */
/*                                                                        */
/* The author can be contacted by email (jfp#hanson-codonw@yahoo.com Anti-*/
/* Spam please change the # in my email to an _)                          */
/*                                                                        */
/* For the latest version and information see                             */
/* http://codonw.sourceforge.net 					  */
/**************************************************************************/
/* This file contains source code for                                     */
/* the core functions involved in correspondence                          */
/* analysis, this code was originally written                             */
/* by Jean Thioulouse                                                     */
/* ADE software: multivariate analysis and graphical                      */
/* display of environmental data                                          */
/* IN Guariso,G and Rizzoli, A (eds),                                     */
/* Software per l'Ambiente. Patron editor, Bolonia, pp.57-62.             */
/*                                                                        */
/* and is used with kind permission                                       */
/*                                                                        */
/* It has however been extensively modified to integrate it               */
/* as seamlessly as practical into CodonW and as such can no              */
/* longer be considered as a stand alone package                          */   
/*                                                                        */
/* Originally written as a general Multivariate analysis (MVA)            */
/* package, it is now hardwired specifically for codon or amino           */
/* acid usage analysis                                                    */
/*                                                                        */
/* All unnecessary functions have been removed                            */
/* Originally each data file had an associated resource file              */
/* which described required parameters                                    */
/* The need for these files has been removed                              */
/*                                                                        */
/**************************************************************************/
/* Functions                                                              */
/* textbin      converts codon usage to binary data file                  */
/*                                                                        */
/**************************************************************************/


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "codonW.h"


/*************** textbin          *****************************************/
/* examines the struct pcoa to see which codons/amino acids are to be inc */
/* in the analysis. It then writes this data to a binary file             */
/* it also counts the amino acid and codon usage of each gene             */
/**************************************************************************/

void textbin(char *fileread, char *fileout) 
{
  double     *vlec;
  int         v2;
  int i,j,x;
  
  pcoa->colm=0;
  
  if ( pm->coa=='a' ) {
    for ( i=1; i<22;i++) 
      if ( pcoa->amino[i] ) pcoa->colm++;  /* number of colms in analysis */
  }else {
    for ( i=1; i<65;i++) 
      if ( pcoa->codons[i]) pcoa->colm++;  /* number of colms in analysis */    
  }
  
  vecalloc(&vlec, pcoa->colm);             /* allocate an array           */
  
  /* open output files                                                    */
  if ( (pm->fcoa_in = open_file( "", fileread, "r", FALSE)) == NULL )  {
    fprintf(pm->my_err,"(txt2bin)");
    my_exit(1,"txt2bin");
  }
  if ( (pm->fcoa_out = open_file( "",fileout, "wb", FALSE)) == NULL )  {
    fprintf(pm->my_err,"(txt2bin)");
    my_exit(6,"fileout");
  }     
  
  for (i=1;i<=pcoa->rows;i++) {            /* pcoa-rows is the No of genes */
    fscanf(pm->fcoa_in,"%s",pm->junk);

    /* read the data from coa_raw into the array vlec                      */
    switch (pm->coa){ 
    case 'a': 
      for (j=1,x=1;j<21;j++) {
        fscanf(pm->fcoa_in,"%i",&v2);
        if ( pcoa->amino[j] )
         vlec[x++] = (double) v2;
        }                                       
        fscanf(pm->fcoa_in,"%i\n",&v2);
        if ( pcoa->amino[j] ) 
          vlec[pcoa->colm]  = (double) v2;
        if ( x != pcoa->colm ) my_exit (99,"Fatal Error in txt2bin");
        break;
    case 'c':
      for (j=1,x=1;j<64;j++) { 
         fscanf(pm->fcoa_in,"%i",&v2);
        if( pcoa->codons[j] )
            vlec[x++] = (double) v2;
        }                       
        fscanf(pm->fcoa_in,"%i\n",&v2);
        if(pcoa->codons[j] ) 
        vlec[pcoa->colm] = (double) v2;
        if ( x != pcoa->colm ) my_exit (99,"Fatal Error in txt2bin");
      break;
    case 'r':
      clean_up ( ncod , naa );   
      for (j=1,x=1;j<64;j++) { 
        fscanf(pm->fcoa_in,"%i",&v2);
        naa[pcu->ca[j]]+=v2;                         /* count amino acids */
        ncod[j]         =v2;                         /* count codons      */
      }                       
      fscanf(pm->fcoa_in,"%i\n",&v2);                /* read last codon   */
      naa[pcu->ca[j]]+=v2;
      ncod[j]         =v2;

      for (j=1,x=0;j<=64;j++) { 
      if(pcoa->codons[j] ) {
     ++x;
         vlec[x] = (double) ((naa[pcu->ca[j]])? 
                  (float) ncod[j]/naa[pcu->ca[j]]*(float)( *(ds+j) ):
            0.00);
     } 
    }
      break;
    
#ifdef DEBUG
    default:
      fprintf(pm->my_err,"error in textbin %c unknown \n",pm->coa );
      break;
#endif
    }                                                          /* end if */
    writevec(vlec, pm->fcoa_out);    
  }
                           /* close files and release memory and return  */
  fileclose(&pm->fcoa_in);
  fileclose(&pm->fcoa_out);  
  free  (vlec); 
}

/*************** colmout          *****************************************/
/* The user has already decided how many axis to be recorded to file      */
/* this value is stored in pcoa->axis. After the analysis is complete the */
/* output data is stored in several binary formatted file. In this case   */
/* nfice and nfics points at the file names.                              */  
/* For each axis that has been requested to be recorded, the position     */
/* of each column (either amino or codon ) is read from the binary file   */
/* and converted into an easily read text file, which is pointed          */
/* at by nfics and the summary file pointed at by summary.                */
/**************************************************************************/
void colmout(char *nfice, char *nfics,AMINO_STRUCT *ppaa, FILE *summary)
{
    double  *vlec;
    int         col, lig=0;
    FILE        *fice=NULL, *fics=NULL;
    float       v2;
    int x,i,j;
    char sp=pm->seperator;

    lig=pcoa->colm;
 
  col=pcoa->axis;                               /* number of axis        */

    vecalloc(&vlec, col);

if( (fice=open_file("",nfice,"rb",FALSE))==NULL) my_exit(6,"nfice2");
if( (fics=open_file("",nfics, "w",FALSE))==NULL) my_exit(1,"nfics2");

fprintf(summary,"\n\nThe position of each %s by axis \n"
	"also see %s for seperate output\n", 
	(pm->coa=='a')? "amino acid":"codon",nfics);

fprintf(fics   , "%s","label");
fprintf(summary, "%-20.20s","label");


for (j=1;j<=col;j++) {
  fprintf(fics   , "%c%s%d",sp,"Axis",j);
  fprintf(summary, "%c%9s%d",sp, "Axis",j);
}  
fprintf(fics   , "\n");
fprintf(summary, "\n");


i=0;
x=1;
  while( x<=lig ) {
                                /* only write out for the columns analysed */
   if( pm->coa == 'a' ) { 
      
       while  ( !pcoa->amino[++i] );                /* skip amino if false */

      fprintf(fics   , "%s%c",ppaa->aa3[i],sp );
      fprintf(summary, "%-20.20s%c",ppaa->aa3[i],sp );
      x++;
    }else{
      
      while  ( !pcoa->codons[++i] );                /* skip codon if false */   
      
      fprintf(fics    , "%s%c",ppaa->cod[i],sp);
      fprintf(summary , "%-20.20s%c",ppaa->cod[i],sp);
      x++;
    }  
        readvec(vlec, fice);
        for (j=1;j<col;j++) {
            v2 = (float) vlec[j];
            fprintf(fics   , "%f%c", v2,sp); 
            fprintf(summary, "%10.5f%c", v2,sp);
        }
        v2 = (float) vlec[col];
        fprintf(fics   , "%f\n", v2); 
        fprintf(summary, "%10.5f\n", v2);
    }
    fileclose(&fics);
    fileclose(&fice);
    free(vlec);
}
/*************** rowout           *****************************************/
/* The position of each gene on each of the principle axis as given by    */
/* pcoa->axis is converted from a binary text file to an ASCII file as    */
/* well as the summary file                                               */
/**************************************************************************/
void rowout(char *nfice, char *nfics, char *ncout, FILE *summary)
{
  double    *vlec, *ax1;
  int           col, lig,*sortax1;
  FILE      *fice=NULL, *fics=NULL, *fnam=NULL;
  float     v2;
  int i,j;
  char sp=pm->seperator;
  
  lig=pcoa->rows;
  col=pcoa->axis;

  vecalloc(&vlec, col);
  vecalloc(&ax1 , lig);
  if( (sortax1= (int *) calloc(lig+1,sizeof(int)))==NULL) 
      my_exit(3,"sortax1"); 
 
  if( (fice=open_file("",nfice,"rb",FALSE))==NULL) my_exit(6,"nfice3");
  if( (fics=open_file("",nfics, "w",FALSE))==NULL) my_exit(1,"nfics3");
  if( (fnam=open_file("",ncout, "r",FALSE))==NULL) my_exit(6,"ncout3");

  fprintf(summary,"\n\nThe position of each gene by axis \n" 
      "(see also %s)\n",nfics);

  fprintf(fics   , "%s%c","label",sp);
  fprintf(summary, "%-20.20s%c","label",sp);

  for (j=1;j<=col;j++) {
    fprintf(fics   , "%s%d%c","Axis",j,sp);
    fprintf(summary, "%9s%d%c", "Axis",j,sp);
  }  
    fprintf(fics   , "\n");
    fprintf(summary, "\n");

  for (i=1;i<=lig;i++) {

    fgets(pm->junk,BUFSIZ,fnam);
    pm->junk[35]='\0';
    for ( j=35 ; j>=0; j--) 
      if ( isspace( (int) pm->junk[j]) ) pm->junk[j]='\0';
 
   fprintf(fics   , "%s%c",pm->junk,sp);
    fprintf(summary, "%-20.20s%c",pm->junk,sp);
    
    readvec(vlec, fice);
    for (j=1;j<col;j++) {
      
      if (j==1)  ax1[i]=vlec[j];                        /* first factors */
      
      v2 = (float) vlec[j];
      fprintf(fics   , "%f%c", v2,sp);
      fprintf(summary, "%10.5f%c", v2,sp);
    }
    v2 = (float) vlec[col];
    fprintf(fics    , "%f\n", v2);
    fprintf(summary , "%10.5f\n", v2);
  } 

 if ( pm->coa != 'a' )   {  
  sorted_by_axis1  ( ax1, sortax1, lig);
  gen_cusort_fop ( sortax1, lig, fnam, summary );
 }
  fileclose(&fics);
  fileclose(&fice);
  fileclose(&fnam);
  free(ax1);
  free(sortax1);
  free(vlec);
}

/************** vecalloc          *****************************************/
/* Allocate memory for a vector of size n and assign that memory to the   */
/* pointer to a pointer vac                                               */
/**************************************************************************/
void vecalloc (double **vec, int n)
{
    if ( (*vec = (double *) calloc(n+1, sizeof(double))) != NULL) {
        **vec = n;
        return;
    } else 
        my_exit(3,"vecalloc");
}

/************** writevec          *****************************************/
/* Write out the value of the vector v1 to a binary file fic              */
/**************************************************************************/
void writevec(double *v1, FILE *fic)
{
    float   v2;
    int     i, c1;

    c1 = (int) v1[0];                       /* Num of vectors             */

    for (i=1;i<=c1;i++) {
        v2 = (float) v1[i];
        if ( fwrite((const char *)&v2, 4, 1, fic) != 1)
        my_exit(4,"writevec");
    }
}

/************** PrepAFC           *****************************************/
/* Calculated Distance matrix for values in contingency table             */
/* Values are first scaled by n (where n is the total usage of a row or   */
/* column                                                                 */
/**************************************************************************/

void PrepAFC(char *nfic)
{
    char bid[17];
    int    i, j;
    double  **w;
    double  *poili, *poico;
    double  a1, a2, x1, n;

/*-------------------------------------------------------------------------*/

    vecalloc(&poili, pcoa->rows);
    vecalloc(&poico, pcoa->colm);
    taballoc(&w, pcoa->rows, pcoa->colm);

    lecmat(w, nfic);

    n = 0;
    for (i=1;i<=pcoa->rows;i++) {
        a1 = 0.0;
        a2 = 0.0;
        for (j=1;j<=pcoa->colm;j++) {
            x1 = w[i][j];
            a1 = a1 + x1;
            poico[j] = poico[j] + x1;
        }
        n = n + a1;
        poili[i] = a1;
    }    
/* scale the vectors, and matrix                                           */
    scalvec(poili, 1.0/n);
    scalvec(poico, 1.0/n);
    scalmat(w, 1.0/n);
    strcpy(bid,"cbfcpl"); 
    ecrvec(poili, bid);
    strcpy(bid,"cbfcpc");
    ecrvec(poico, bid);

/*-------------------------------------------------------------------------*/

    for (i=1;i<=pcoa->rows;i++) {
        a1 = poili[i];
        if (a1 != 0.0) {
            for (j=1;j<=pcoa->colm;j++) {
                a2 = poico[j];
                if (a2 != 0) w[i][j] = w[i][j] / a1 / a2 - 1;
            }
        }
    }
    strcpy(bid,"cbfcta");
    ecrmat(w, bid);

/*-------------------------------------------------------------------------*/
    freetab(w);
    freevec(poili);
    freevec(poico);
    pcoa->inertia = (float) inertot ();
}

/************** inertot         ********************************************/
/* Calculate total data inertia                                            */
/***************************************************************************/

double inertot ( void )
{
    int     i, j; 
    double      **tab;
    double  *pl, *pc;
    double  a1, s1, inertia;
    taballoc (&tab, pcoa->rows, pcoa->colm);
    vecalloc (&pc, pcoa->colm);
    vecalloc (&pl, pcoa->rows);

    lecmat (tab,"cbfcta");
    lecvec(pl, "cbfcpl");
    lecvec(pc, "cbfcpc");
    inertia = 0;
    for (i=1;i<=pcoa->rows;i++) {
        a1 = pl[i];
        for (j=1;j<=pcoa->colm;j++) {
            s1 = tab[i][j];
            inertia = inertia + s1 * s1 * a1 * pc[j];
        }
    }   
    freetab(tab);
    freevec(pl);
    freevec(pc);
    
    return inertia;
}

/************** lecmat            *****************************************/
/* Opens binary file nfic, reads the values it contains and records them  */
/* in the matrix pointed to by tab                                        */
/**************************************************************************/
void lecmat (double **tab, char *nfic)
{
    int     i, j, l1, c1;
    float   v2;
    FILE    *fic=NULL;
    
    l1 = (int) tab[0][0];
    c1 = (int) tab[1][0];

    if( (fic=open_file("",nfic,"rb",FALSE))==NULL) my_exit(1,"lecmat");


    for (i=1;i<=l1;i++) {
        for (j=1;j<=c1;j++) {
            if ( fread((char *)&v2, 4, 1, fic) != 1) {
            fprintf(pm->my_err,"Error: can't read matrix (lecmat)");
            my_exit(5,"lecmat");
            }
            tab[i][j] = v2;
        }
    }  
    fileclose(&fic);
}

/************** freetab           *****************************************/
/* Releases memory dynamically allocated to a table tab(x,y)              */
/**************************************************************************/
void freetab (double **tab)
{
    int     i, n;
    n = (int) *(*(tab));                /* number of rows in table        */
    for (i=0;i<=n;i++) {
            free((char *) *(tab+i) );
    }
    free((char *) tab);
}

/************** freevec           *****************************************/
/* Releases memory dynamically allocated to a vector                      */
/**************************************************************************/
void freevec (double *vec)
{   
    free((char *) vec); 
}  

/************** taballoc          *****************************************/
/* Dynamically allocates memory to the table tab(l1,c1)                   */
/**************************************************************************/
void taballoc (double ***tab, int l1, int c1)
{
    int     i;
    
    if ( (*tab = (double **) calloc(l1+1, sizeof(double *))) != NULL) {
        for (i=0;i<=l1;i++) {
            if ( (*(*tab+i)=(double *) calloc(c1+1, sizeof(double))) == NULL ) {
            fprintf(pm->my_err,"(taballoc)");
            my_exit(3,"taballoc");
            }           
        }
    } else {fprintf(pm->my_err,"(taballoc)");
            my_exit(3,"taballoc2");
            }

    **(*tab) = l1;
    **(*tab+1) = c1;
}

/**************   lecvec          *****************************************/
/* Reads vectors from filename *nfic and assigns them to a vector         */
/**************************************************************************/
void lecvec (double *v1, char *nfic)
{
    float   v2;
    int     i, c1;
    FILE    *fic=NULL;
    
    if( (fic=open_file("",nfic,"rb",FALSE))==NULL) my_exit(6,"lecvec");

    c1 = (int) v1[0];
    for (i=1;i<=c1;i++) {
        if ( fread((char *)&v2, 4, 1, fic) != 1){
        fprintf(pm->my_err,"(lecvec)");
        my_exit(5,"lecvec");
        }
        v1[i] = v2;
    }    
    fileclose(&fic);
}

/************** ecrmat           ******************************************/
/* Writes the table pointed to by **tab to the binary filename *nfic      */
/**************************************************************************/
void ecrmat (double **tab, char *nfic)
{
    int     i, j, l1, c1;
    float   v2;
    FILE    *fic=NULL;
    
    l1 = (int)tab[0][0];
    c1 = (int)tab[1][0];

    if( (fic=open_file("",nfic,"wb",FALSE))==NULL) my_exit(1,"ecrmat");

    for (i=1;i<=l1;i++) {
        for (j=1;j<=c1;j++) {
            v2 = (float)tab[i][j];
            if ( fwrite((const char *)&v2, 4, 1, fic) != 1)  {
            fprintf(pm->my_err,"(ecrmat)");
            my_exit(4,"ecrmat");         
            }
        }
    }
    
    fileclose(&fic);
}
/************** ecrvec           ******************************************/
/* Writes the pointer pointed to by *v1 to the binary file *nfic          */
/**************************************************************************/
void ecrvec (double *v1, char *nfic)
{
    float   v2;
    int     i, c1;
    FILE    *fic=NULL;

    c1 = (int)v1[0];

    if( (fic=open_file("",nfic,"wb",FALSE))==NULL) my_exit(1,"ecrvec");
    
    
    for (i=1;i<=c1;i++) {
        v2 = (float)v1[i];
        if ( fwrite((const char *)&v2, 4, 1, fic) != 1){
            fprintf(pm->my_err,"(ecrvec)");
            my_exit(4,"ecrvec");         
            }
    }
    
    fileclose(&fic);
}

/************** scalmat          ******************************************/
/* Scale the matrix pointed to by **tab by r                              */
/**************************************************************************/
void scalmat (double **tab, double r)
{
    int l1, c1, i, j;

    l1 = (int) tab[0][0];
    c1 = (int) tab[1][0];
    for (i=1;i<=l1;i++) {
        for (j=1;j<=c1;j++) {
            tab[i][j] = tab[i][j] * r;
        }
    }
}

/************** scalvec          ******************************************/
/* Scale the vector pointed to by *v1 by r                                */
/**************************************************************************/
void scalvec (double *v1, double r)
{
    int i, c1;
    
    c1 = (int) v1[0];
    
    for (i=1;i<=c1;i++) {
        v1[i] = v1[i] * r;
    }
}

/************** DiagoRC          ******************************************/
/* This function generates/calculates the correspondence analysis factors */
/**************************************************************************/
void DiagoRC ( FILE *summary)
{
  int           lcmin, rang, f1, i, j, k;
  double           **w, **ctab, **auxi, **vp1, **vp2;
  double        *poili, *poico, *l;
  double         s, s1, a1, inertotal;

    
  lcmin = pcoa->colm;
  if (pcoa->rows < pcoa->colm) lcmin = pcoa->rows;
  taballoc(&w, pcoa->rows, pcoa->colm);
  taballoc(&ctab, lcmin, lcmin);
  taballoc(&auxi, lcmin, 2);
  vecalloc(&poili, pcoa->rows);
  vecalloc(&poico, pcoa->colm);
  vecalloc(&l, lcmin);

  lecvec(poili, "cbfcpl");
  sqrvec(poili);
  lecvec(poico, "cbfcpc");
  sqrvec(poico);
  lecmat(w, "cbfcta");
  
  inertotal=0;
  for (i=1;i<=pcoa->rows;i++) {
    a1 = poili[i];
    for (j=1;j<=pcoa->colm;j++) {
      s1 = w[i][j] * a1 * poico[j];
      w[i][j] = s1;
      s1 = s1 * s1;
      inertotal = inertotal + s1;
    }
  }
  
  fprintf(summary,"The total inertia of the data was %f\n",inertotal); 
  fprintf(summary, "\nExplanation of the variation by axis "
                    "(see also eigen.coa)\n");


/*  prodmatAAtB and prodmatAtAB calc product of the scaled distance matrix */
/*  DiagoComp diagnolises the product matrix ctab                          */
/*  editvalpro output the eigen values                                     */

dot(1,10);
    if (pcoa->rows < pcoa->colm) {

        prodmatAAtB(w, ctab);                            
        DiagoComp(pcoa->rows, ctab, l, &rang);
        f1=pcoa->axis;
        editvalpro(summary, l, pcoa->rows, inertotal);
        for (j=1;j<=pcoa->rows;j++) {
            auxi[j][1] = l[j];
            auxi[j][2] = l[j]/inertotal;
        }
        sqrvec(l);
    } else {
        prodmatAtAB(w, ctab);
        DiagoComp(pcoa->colm, ctab, l, &rang);
        f1=pcoa->axis;
        editvalpro(summary, l, pcoa->colm, inertotal);
        for (j=1;j<=pcoa->colm;j++) {
            auxi[j][1] = l[j];
            auxi[j][2] = l[j]/inertotal;
        }
        sqrvec(l);
    }   

    if (f1==0) {
        if (lcmin == 1) f1 = 1;
        else f1 = 2;
    }

   /* output the relative inertia values                                   */
    ecrmat(auxi, "cbfcvp");

   /* Calculate the factorial coordinates                                  */

    if (pcoa->rows < pcoa->colm) {
        taballoc(&vp2, pcoa->colm, f1);
        for (j=1;j<=pcoa->colm;j++) {
            for (k=1;k<=f1;k++) {
                s = 0;
                for (i=1;i<=pcoa->rows;i++) {
                    s = s + w[i][j] * ctab[i][k];
                }
            vp2[j][k] = s;
            }       
        }
        for (i=1;i<=pcoa->colm;i++) {
            if (poico[i] != 0) {
                for (j=1;j<=f1;j++) {
                    vp2[i][j] = vp2[i][j] / poico[i];
                }
            }
        }
        for (i=1;i<=pcoa->rows;i++) {
            if (poili[i] != 0) {
                for (j=1;j<=pcoa->rows;j++) {
                    ctab[i][j] = ctab[i][j] * l[j] / poili[i];
                }
            }
        }
        ecrmatred(ctab, f1, "cbfcli");
        ecrmatred(vp2, f1,  "cbfcco");
        freetab(vp2);
    } else {
        taballoc(&vp1, pcoa->colm, f1);
        taballoc(&vp2, pcoa->rows, f1);
        for (i=1;i<=pcoa->colm;i++) {
            for (j=1;j<=f1;j++) {
                vp1[i][j] = ctab[i][j];
            }
        }
        prodmatABC(w, vp1, vp2);
        for (i=1;i<=pcoa->rows;i++) {
            if (poili[i] != 0.0) {
                for (j=1;j<=f1;j++) {
                    vp2[i][j] = vp2[i][j] / poili[i];
                }
            }
        }
        for (i=1;i<=pcoa->colm;i++) {
            if (poico[i] != 0) {
                for (j=1;j<=rang;j++) {
                    ctab[i][j] = ctab[i][j] * l[j] / poico[i];
                }
            }
        }
        ecrmat(vp2, "cbfcli");
        ecrmatred(ctab, f1, "cbfcco");
        freetab(vp1);
        freetab(vp2);
    }

    goto fin;

/* free memory                                                            */

fin:
    freetab(w);
    freetab(ctab);
    freetab(auxi);
    freevec(poili);
    freevec(poico);
    freevec(l);
    
}                                                       /* End of DiagoRC */
 
/************** sqrvec           ******************************************/
/* This function calculates the square root of a vector                   */
/**************************************************************************/
void sqrvec (double *v1)
{
    int i, c1;
    double v2;
    
    c1 = (int) v1[0];
    
    for (i=1;i<=c1;i++) {
        v2 = v1[i];
        if (v2 < 0.0) {
        fprintf(pm->my_err,"Error: Square root of negative number (sqrvec)");
        my_exit(99,"sqrvec");        
        }
        v2 = sqrt(v2);
        v1[i] = v2;
    }
}

/************** prodmatAAtB         ***************************************/
/* Calculate the product of matrix a*a and return it as matrix b          */
/**************************************************************************/
void prodmatAAtB (double **a, double **b)
{
    int j, k, i, lig, col;
    double s;
    
    lig = (int) a[0][0];
    col = (int) a[1][0];

    for (j=1;j<=lig;j++) { 
    dot ( 1 , 10 );
        for (k=j;k<=lig;k++) {
            s = 0;
            for (i=1;i<=col;i++) {
                s = s + a[j][i] * a[k][i];
            }
        b[j][k] = s;
        b[k][j] = s;
        }       
    }
}

/************** prodmatABC          ***************************************/
/* Calculate the product of matrix a*b and return it as matrix c          */
/**************************************************************************/
void prodmatABC (double **a, double **b, double **c)
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = (int) a[0][0];
    col = (int) a[1][0];
    
    col2 = (int) b[1][0];

    for (i=1;i<=lig;i++) {
    dot(1,10);
        for (k=1;k<=col2;k++) {
            s = 0;
            for (j=1;j<=col;j++) {
                s = s + a[i][j] * b[j][k];
            }
        c[i][k] = s;
        }       
    }
}

/************** prodmatAtAB         ***************************************/
/* Calculate the product of matrix a*A and return it as matrix b          */
/**************************************************************************/
void prodmatAtAB (double **a, double **b)
{
    int j, k, i, lig, col;
    double s;
    
    lig = (int) a[0][0];
    col = (int) a[1][0];

    for (j=1;j<=col;j++) { 
     dot(1,100);
        for (k=j;k<=col;k++) {
            s = 0;
            for (i=1;i<=lig;i++) {
                s = s + a[i][k] * a[i][j];
            }
        b[j][k] = s;
        b[k][j] = s;
        }       
    }
}

/**************  editvalpro         ***************************************/
/* Calculate eigenvalues, relative inertia and Sum of inertia for each    */
/* factor and record this to eigen.coa and summary.coa                    */
/**************************************************************************/
void editvalpro (FILE *ficlist, double *vp, int n, double s)
{
  double        sc1, sc2;
  int           i, n1;
  float         v2, v3, v4;
  FILE *eigen=NULL;  
  char sp;

  sp=pm->seperator;

  if ( (eigen=open_file("","eigen.coa","w",FALSE))==NULL ) 
          my_exit(1,"editvalpro");


  sc1 = 0.0;
  for (i=1;i<=n;i++) {
    if (vp[i] < 0.0) {
      v2 = (float) vp[i];
      fprintf(ficlist, "Eigenvalue number %d is negative : %+.4E\n", i, v2);
      vp[i] = 0.0;
    }
  }
  n1 = (n > 40) ? 40 : n;
  fprintf(ficlist, "Num. Eigenval.   R.Iner.  R.Sum    "
			"|Num. Eigenval.   R.Iner.  R.Sum  |");
  fprintf(ficlist, "\n");
  for (i=1;i<=n1;i=i+2) {
    sc1 = sc1 + vp[i];
    if (i < n1) {
      sc2 = sc1 + vp[i+1];
      v2 = (float) vp[i];
      v3 = (float)vp[i]/(float)s;
      v4 = (float)sc1/(float)s;
      fprintf(ficlist, "%.2d   %+.4E %+.4f %+.4f ", i, v2, v3, v4);
      fprintf(eigen ,"%.2d%c%.4E%c%.4f%c%.4f\n",i,sp,v2,sp,v3,sp,v4);
      v2 = (float)vp[i+1];
      v3 = (float)vp[i+1]/(float)s;
      v4 = (float)sc2/(float)s;
      fprintf(ficlist, "  |%.2d   %+.4E %+.4f %+.4f |", i+1, v2, v3, v4);
      fprintf(eigen ,"%.2d%c%.4E%c%.4f%c%.4f\n",i+1,sp,v2,sp,v3,sp,v4);
    } else {
      v2 = (float)vp[i];
      v3 = (float)vp[i]/(float)s;
      v4 = (float)sc1/(float)s;
      fprintf(ficlist, "%.2d   %+.4E %+.4f %+.4f ", i, v2, v3, v4);
      fprintf(eigen ,"%.2d%c%.4E%c%.4f%c%.4f\n",i,sp,v2,sp,v3,sp,v4);
    }
    sc1 = sc2;
    fprintf(ficlist, "\n");
  }
  fprintf(ficlist, "\n");
fileclose(&eigen);
}

/**************  ecrmatred        *****************************************/
/* Output c1 columns of matrix tab to filename *nfic                      */
/**************************************************************************/
void ecrmatred (double **tab, int c1, char *nfic)
{
    int     i, j, l1;
    float   v2;
    FILE    *fic=NULL;
    
    l1 = (int) tab[0][0];

    if( (fic=open_file("",nfic,"wb",FALSE))==NULL) my_exit(1,"ecrmatred");

    for (i=1;i<=l1;i++) {
        for (j=1;j<=c1;j++) {
            v2 = (float) tab[i][j];
            if ( fwrite((const char *)&v2, 4, 1, fic) != 1){
        fprintf(pm->my_err,"(ecrmatred)");
        my_exit(4,"ecrmatred");     
        }
        }
    }
    
    fileclose(&fic);
}

/**************  readvec            ***************************************/
/* read vector v1 from filehandle fic                                     */
/**************************************************************************/
void readvec (double *v1, FILE *fic)
{
    float   v2;
    int     i, c1;

    c1 = (int) v1[0];

    for (i=1;i<=c1;i++) {
        if ( fread((char *)&v2, 4, 1, fic) != 1) {
        fprintf(pm->my_err,"(readvec)");
        my_exit(5,"readvec");     
        }
        v1[i] = v2;
    }
}

/**************  DiagoComp         ***************************************/
/* Diagnolisation of matrix w                                            */
/* T. FOUCART Analyse factorielle de tableaux multiples,                 */
/* Masson, Paris 1984,185p., p. 62. D'aprhs VPROP et TRIDI,              */
/* de LEBART et coll.                                                    */
/* Lots of nasty goto jumps ... ported from Fortran                      */
/*************************************************************************/
void DiagoComp (int n0, double **w, double *d, int *rang)
{
    double          *s;
    double          a, b, c, x, xp, q, bp, ab, ep, h, t, u , v;
    double          dble;
    int             ni, i, i2, j, k, jk, ijk, ij, l, ix, m, m1, isnou;
    
    vecalloc(&s, n0);
    a = 0.000000001;
    ni = 100;
    if (n0 == 1) {
        d[1] = w[1][1];
        w[1][1] = 1.0;
        *rang = 1;
        freevec (s);
        return;
    }
    
    for (i2=2;i2<=n0;i2++) {
       
        b=0.0;
        c=0.0;
        i=n0-i2+2;
        k=i-1;
        if (k < 2) goto Et1;
        for (l=1;l<=k;l++) {
            c = c + fabs((double) w[i][l]);
        }
        if (c != 0.0) goto Et2;
        
Et1:    s[i] = w[i][k];
        goto Etc;
        
Et2:    for (l=1;l<=k;l++) {
            x = w[i][l] / c;
            w[i][l] = x;
            b = b + x * x;
        }
        xp = w[i][k];
        ix = 1;
        if (xp < 0.0) ix = -1;
        
/*      q = -sqrt(b) * ix; */
        dble = b;
        dble = -sqrt(dble);
        q = dble * ix;

        s[i] = c * q;
        b = b - xp * q;
        w[i][k] = xp - q;
        xp = 0;
        for (m=1;m<=k;m++) {
            w[m][i] = w[i][m] / b / c;
            q = 0;
            for (l=1;l<=m;l++) {
                q = q + w[m][l] * w[i][l];
            }
            m1 = m + 1;
            if (k < m1) goto Et3;
            for (l=m1;l<=k;l++) {
                q = q + w[l][m] * w[i][l];
            }
            
Et3:        s[m] = q / b;
            xp = xp + s[m] * w[i][m];
        }
        bp = xp * 0.5 / b;
        for (m=1;m<=k;m++) {
            xp = w[i][m];
            q = s[m] - bp * xp;
            s[m] = q;
            for (l=1;l<=m;l++) {
                w[m][l] = w[m][l] - xp * s[l] - q * w[i][l];
            }
        }
        for (l=1;l<=k;l++) {
            w[i][l] = c * w[i][l];
        }
        
Etc:    d[i] = b;
    } /* for (i2=2;i2<n0;i2++) */
    
    s[1] = 0.0;
    d[1] = 0.0;
    
    for (i=1;i<=n0;i++) {
     dot(1,100);
        k = i - 1;
        if (d[i] == 0.0) goto Et4;
        for (m=1;m<=k;m++) {
            q = 0.0;
            for (l=1;l<=k;l++) {
                q = q + w[i][l] * w[l][m];
            }
            for (l=1;l<=k;l++) {
                w[l][m] = w[l][m] - q * w[l][i];
            }
        }
        
Et4:    d[i] = w[i][i];
        w[i][i] = 1.0;
        if (k < 1) goto Et5;
        for (m=1;m<=k;m++) {
            w[i][m] = 0.0;
            w[m][i] = 0.0;
        }

Et5:;
    }
    
    for (i=2;i<=n0;i++) {
        s[i-1] = s[i];
    }
    s[n0] = 0.0;
    for (k=1;k<=n0;k++) {
        m = 0;

Et6:    for (j=k;j<=n0;j++) {
     dot(1,100);
            if (j == n0) goto Et7;
            ab = fabs((double) s[j]);
            ep = a * (fabs((double) d[j]) + fabs((double) d[j+1]));
            if (ab < ep) goto Et7;
        }
    
Et7:    isnou = 1;
        h = d[k];
        if (j == k) goto Eta;
        if (m < ni) goto Etd;
        
        fprintf(pm->my_err,"Error: can't compute matrix eigenvalues");
        my_exit(99,"corresp");
        
Etd:    m = m + 1;
        q = (d[k+1]-h) * 0.5 / s[k];
        
/*      t = sqrt(q * q + 1.0); */
        dble = q * q + 1.0;
        dble = sqrt(dble);
        t = dble;
        
        if (q < 0.0) isnou = -1;
        q = d[j] - h + s[k] / (q + t * isnou);
        u = 1.0;
        v = 1.0;
        h = 0.0;
        jk = j-k;
        for (ijk=1;ijk<=jk;ijk++) {
    dot(1,100);
            i = j - ijk;
            xp = u * s[i];
            b = v * s[i];
            if (fabs((double) xp) < fabs((double) q)) goto Et8;
            u = xp / q;
            
/*          t = sqrt(u * u + 1); */
            dble = u * u + 1.0;
            dble = sqrt(dble);
            t = dble;
            
            s[i+1] = q * t;
            v = 1 / t;
            u = u * v;
            goto Et9;

Et8:        v = q / xp;

/*          t = sqrt(1 + v * v); */
            dble = 1.0 + v * v;
            dble = sqrt(dble);
            t = dble;
            
            s[i+1] = t * xp;
            u = 1 / t;
            v = v * u;

Et9:
            q = d[i+1] - h;
            t = (d[i] - q) * u + 2.0 * v * b;
            h = u * t;
            d[i+1] = q + h;
            q = v * t - b;
            for (l=1;l<=n0;l++) {
                xp = w[l][i+1];
                w[l][i+1] = u * w[l][i] + v * xp;
                w[l][i] = v * w[l][i] - u * xp;
            }
        }
        d[k] = d[k] - h;
        s[k] = q;
        s[j] = 0.0;
        goto Et6;

Eta:;
    } /* for (k=1;k<=n0;k++) */
    
    for (ij=2;ij<=n0;ij++) {
     dot(1,300);
        i = ij - 1;
        l = i;
        h = d[i];
        for (m=ij;m<=n0;m++) {
            if (d[m] >= h) {
                l = m;
                h = d[m];
            }
        }
        if (l == i) {
            goto Etb;
        } else {
            d[l] = d[i];
            d[i] = h;
        }
        for (m=1;m<=n0;m++) {
            h = w[m][i];
            w[m][i] = w[m][l];
            w[m][l] = h;
        }

Etb:;
    } /* for (ij=2;ij<=n0;ij++) */

    *rang = 0;
    for (i=1;i<=n0;i++) {
        if (d[i] / d[1] < 0.00001) d[i] = 0.0;
        if (d[i] != 0.0) *rang = *rang + 1;
    }
    freevec(s);
} /* DiagoComp */

/************** inertialig         ***************************************/
/* Called when advanced correspondence analysis option has been selected */
/* This analyses and reports the absolute and relative contributions of  */
/* each gene to the inertia of the principal factors (by default the     */
/* first 4 axis)                                                         */
/*************************************************************************/
void inertialig( char *inertia_out, char *ncout, FILE *summary)
{
    int     i, j, k, f1, l1,c1,lcmin;
    double  **cooli, **w;
    double      *vtab, *conli, *poili, *poico;
    double      l0, inertotal, a1, a2, m2, m3, s1;
    double      temp1=0,temp2=0;
    FILE *inert_out=NULL,*fnam=NULL;
    
   l1       =pcoa->rows;
   c1       =pcoa->colm;
   f1       =pcoa->axis;
   inertotal  =pcoa->inertia;
   
    if( (inert_out=open_file( "",inertia_out,"w",FALSE))==NULL) 
	      my_exit(1,"inertia out");
    
    lcmin = c1; if (l1<lcmin) lcmin=l1;
    taballoc (&w, l1,c1);
    vecalloc(&poili, l1);
    vecalloc(&poico, c1);
    taballoc(&cooli, l1, f1);
    vecalloc(&conli, l1);
    vecalloc(&vtab, lcmin);

    lecvec(poili, "cbfcpl");
    sqrvec(poili);
    lecvec(poico, "cbfcpc");
    sqrvec(poico);
    lecmat(cooli, "cbfcli");
    selectcol("cbfcvp", vtab, 2);
    lecmat(w, "cbfcta");

    fprintf(summary, "\n\nNumber of rows: %d, columns: %d\n", l1, c1);
    fprintf(summary, "Total inertia: %8.6G - Number of axes: %d\n\n", 
	                  inertotal, f1);
    fprintf(summary, "Contributions of each gene to the recorded factors "
                     "A.K.A axes\n");

/* calculate the contribution                                              */

    for (i=1;i<=l1;i++) {
        a1 = poili[i];
        for (j=1;j<=c1;j++) {
            s1 = w[i][j] * a1 * poico[j];
            s1 = s1 * s1;
            conli[i] = conli[i] + s1;
        }
    }

/* scale the vectors by 1/inertia total                                    */

    scalvec(conli, 1.0/inertotal);    

    
  if( (fnam=open_file("",ncout, "r",FALSE))==NULL) my_exit(6,"inertialgn");

    fprintf(summary, "Row inertia\n");
    fprintf(summary, "All contributions are in 1/10000\n\n");
    fprintf(summary, "----------Absolute contributions----------\n");
    fprintf(summary, "Short_Gene_Name|Num  |");
    for (k=1;k<=f1;k++) {
        fprintf(summary, "Fac%2d|", k);
        
    }
    fprintf(summary  , "\n");
    fprintf(inert_out, "\n");
    for (i=1;i<=l1;i++) {

      fgets(pm->junk,BUFSIZ,fnam);
      pm->junk[35]='\0';
      for ( j=35 ; j>=0; j--) if ( isspace((int)pm->junk[j]) ) 
	       pm->junk[j]='\0';
      
      fprintf(inert_out   ,"%-.15s%c",pm->junk,pm->seperator);
      fprintf(summary, "%-15.15s",pm->junk);
      
      fprintf(summary  ,"|%5d|", i);
      fprintf(inert_out,"%d%c", i,pm->seperator);
      
      l0 = poili[i]*poili[i]/inertotal;

      for (j=1;j<=f1;j++) {
	temp1=(cooli[i][j] * cooli[i][j]);         /* bug fix for Think C      */
	temp2=(l0 / vtab[j]);                      /* need to split calculation*/
	a1 = temp1 * temp2;
	fprintf(summary, "%5d|", (int) (a1 * 10000));
	fprintf(inert_out,"%d%c",(int) (a1 * 10000),pm->seperator);
      }
        fprintf(summary, "\n");
        fprintf(inert_out,"\n");
    }
    fprintf(summary, "\n\nRelative contributions\nThis is the variation \n"
	    "in the %s usage of each gene that is \n"
	    "explained by each axis/factor\n"
	    "see also %s \n",
	    (pm->coa=='a')?"amino acid":"codon",inertia_out);

    fclose(fnam);
    if( (fnam=open_file("",ncout, "r",FALSE))==NULL) 
        my_exit(6,"inertialgn");

    fprintf(summary, "----------Relative contributions----------\n");
    fprintf(summary, "Short_gene_name|Num  |");
    for (k=1;k<=f1;k++) {
        fprintf(summary, "Fac%2d|", k);
    }
    fprintf(summary, "|Remains| Weight | Cont.|");
    fprintf(summary, "\n");
    fprintf(inert_out,"\n");
    
    for (i=1;i<=l1;i++) {

      fgets(pm->junk,BUFSIZ,fnam);
      pm->junk[35]='\0';
      for ( j=35 ; j>=0; j--) if ( isspace( (int) pm->junk[j]) ) 
          pm->junk[j]='\0';
      
      fprintf(inert_out   , "%-.15s%c",pm->junk,pm->seperator);
      fprintf(summary, "%-15.15s",pm->junk);


        fprintf(summary, "|%5d|", i);  
        fprintf(inert_out,"%d%c", i,pm->seperator);
        a2 = 0.;
        m3 = poili[i]*poili[i]/inertotal;
        m2 = conli[i];
        if (m2 == 0.) m2 = 1.;
        for (j=1;j<=f1;j++) {
            a1 = cooli[i][j] * cooli[i][j] * m3 / m2;
            a2 = a2 + a1;
            fprintf(summary, "%5d|", (int) (a1 * 10000)); 
            fprintf(inert_out,"%d%c",(int) (a1 * 10000),pm->seperator);
        }
        fprintf(summary, "|%5d  ", (int) ((1-a2) * 10000));
        fprintf(summary, "|%5d   |%5d |\n", (int) (inertotal * m3 * 10000), 
            (int) (m2 * 10000));
        fprintf(inert_out, "\n");
    }
    fprintf(summary  , "\n");
    fprintf(inert_out, "\n");
    
                                                        /* free memory    */  
    freetab(w);
    freevec(poili);
    freevec(poico);
    freetab(cooli);
    freevec(conli);
    freevec(vtab);
    fileclose(&inert_out);
    fileclose(&fnam);

}                                                       /* End of Inertia */

/************** inertiacol         ****************************************/
/* Called when advanced correspondence analysis option has been selected  */
/* This analyses and reports the absolute and relative contributions of   */
/* each codon or amino acid to the inertia of the principal factors (by   */
/* default the first 4 axis)                                              */
/**************************************************************************/
void inertiacol(char *inertia_out, FILE *summary )
{
    int             x,i, j, k, f1, l1,c1, lcmin;
    double  **cooco, **w;
    double      *vtab, *conco, *poili, *poico;
    double      l0, inertotal, a1, a2, m2, m3, s1;
    FILE *inert_out=NULL;
    
    if( (inert_out=open_file( "",inertia_out,"a",FALSE))==NULL) 
		my_exit(1,"inertia out2");

    l1      =pcoa->rows;
    c1      =pcoa->colm;
    f1      =pcoa->axis;
    inertotal =pcoa->inertia;

    lcmin = c1; if (l1<lcmin) lcmin=l1;

    taballoc (&w, l1,c1);
    vecalloc(&poili, l1);
    vecalloc(&poico, c1);
    taballoc(&cooco, c1, f1);
    vecalloc(&conco, c1);
    vecalloc(&vtab, lcmin);

    lecvec(poili, "cbfcpl");
    sqrvec(poili);
    lecvec(poico, "cbfcpc");
    sqrvec(poico);
    lecmat(cooco, "cbfcco");
    selectcol("cbfcvp", vtab, 2);
    lecmat(w, "cbfcta");

    fprintf(summary, "\n\nColumn inertia\nNumber of genes: %d, columns: "
	                 "%d\n\n", l1, c1);
    fprintf(summary, "This is the fraction of the total inertia that is\n"
	    "explained for each %s by each of the recorded\n"
            "factors or axes\n\n\n",(pm->coa=='a')? "amino acids":"codons");


    for (i=1;i<=l1;i++) {
        a1 = poili[i];
        for (j=1;j<=c1;j++) {
            s1 = w[i][j] * a1 * poico[j];
            s1 = s1 * s1;
            conco[j] = conco[j] + s1;
        }
    }

    /* scale the vectors by 1/inertia total                                     */
    scalvec(conco, 1.0/inertotal);
    
    fprintf(summary, "\n\nColumn inertia\n");
    fprintf(summary, "All contributions are in 1/10000\n\n");
    fprintf(summary, "----------Absolute contributions----------\n");
    fprintf(summary, "Key|Num  |");
    for (k=1;k<=f1;k++) {
        fprintf(summary, "Fac%2d|", k);
    }
    fprintf(summary, "\n");
    for (x=0,i=1;i<=c1;i++) {

      if (pm->coa == 'a' ){
        
          while(pcoa->amino[++x] == FALSE);
       
        fprintf(summary, "%s", paa->aa3[x]); 
        fprintf(inert_out,"%s%c",paa->aa3[x],pm->seperator);      
      }else{ 
	
          while(pcoa->codons[++x] == FALSE);
        
          fprintf(summary, "%s", paa->cod[x]); 
        fprintf(inert_out,"%s%c",paa->cod[x],pm->seperator);      	
      }
      
      fprintf(summary, "|%5d|", i); 
      fprintf(inert_out,"%d%c",i,pm->seperator);
      
      l0 = poico[i]*poico[i]/inertotal;
      for (j=1;j<=f1;j++) {
	a1 = cooco[i][j] * cooco[i][j] * l0 / vtab[j];
	fprintf(summary,  "%5d|", (int) (a1 * 10000));
	fprintf(inert_out,"%i%c", (int) (a1 * 10000),pm->seperator );
      }
        fprintf(summary,  "\n");
        fprintf(inert_out,"\n"); 
    }
    fprintf(summary,  "\n"); 
    fprintf(inert_out,"\n"); 
    fprintf(summary, "----------Relative contributions----------\n");
    fprintf(summary, "Key|Num  |");
    for (k=1;k<=f1;k++) {
        fprintf(summary, "Fac%2d|", k);
    }
    fprintf(summary, "|Remains| Weight | Cont.|");
    fprintf(summary, "\n");
    for (x=0,i=1;i<=c1;i++) {


      if (pm->coa == 'a' ){
        
          while(pcoa->amino[++x] == FALSE);
        
        fprintf(summary, "%s", paa->aa3[x]); 
        fprintf(inert_out,"%s%c",paa->aa3[x],pm->seperator);      
	}else{ 
	
          while(pcoa->codons[++x] == FALSE);
        
        fprintf(summary, "%s", paa->cod[x]); 
        fprintf(inert_out,"%s%c",paa->cod[x],pm->seperator);      	
	}

        fprintf(summary, "|%5d|", i); 
        fprintf(inert_out,"%d%c",i,pm->seperator);
        a2 = 0.;
        m3 = poico[i]*poico[i]/inertotal;
        m2 = conco[i];
        if (m2 == 0.) m2 = 1.;
        for (j=1;j<=f1;j++) {
            a1 = cooco[i][j] * cooco[i][j] * m3 / m2;
            a2 = a2 + a1;
            fprintf(summary, "%5d|", (int) (a1 * 10000)); 
            fprintf(inert_out,"%d%c",(int) (a1 * 10000),pm->seperator);
        }
        fprintf(summary, "|%5d  ", (int) ((1-a2) * 10000));
        fprintf(summary, "|%5d   |%5d |\n", 
			(int) (inertotal * m3 * 10000), (int) (m2 * 10000));
        fprintf(inert_out,"\n");
    }
    fprintf(summary, "\n");
    
    freetab(w);
    freetab(cooco);
    freevec(poili);
    freevec(poico);
    freevec(conco);
    freevec(vtab);    
} 									/* End of Inertia */       

/**************  selectcol         ***************************************/
/* extract a column from the file *nfic, column has the dimension of the */
/* number of genes. If these disagree it will about. Col is the number of*/
/* the column to extract.                                                */
/*************************************************************************/
void selectcol (char *nfic , double *col, int numcol)
{
    FILE    *fic=NULL;
    int i, c1,l1;
    double *vlec;   
    
    c1=2;
    l1=( pcoa->rows < pcoa->colm)? pcoa->rows:pcoa->colm;
    
    vecalloc(&vlec, c1);


    if (numcol>c1) {
        fprintf (pm->my_err,"fatal input-output error numcol>c1 (selectcol");
        my_exit(99,"corresp");
    }
    
if( (fic=open_file( "",nfic,"rb",FALSE))==NULL) my_exit(6,"nfic4");
    for (i=1;i<=l1;i++) {
        readvec(vlec, fic);
        col[i] = vlec[numcol];
    }
    
    fileclose(&fic);
    freevec(vlec);
}
 
/**************  suprow            ***************************************/
/* This sub adds supplementary genes after the correspondence analysis   */
/* has completed for an initial set of genes. The supplementary genes are*/
/* read in and processed up to the point of the generation of factors    */
/* at which point the factors for the initial analysis are used to calc  */
/* the position of the supplementary genes on the originally identified  */
/* axis                                                                  */
/*************************************************************************/
void suprow (int num_seq, char *nficvp, char *nfictasup, char *nficlisup, 
char*option , FILE *summary)
{
    int         l1,c1,l2,c2,i,j,k;
    double      **compos, **tabsup;
    double      *vp, *poico;
    double      *moy, *var;
    double      a1, a2;
    FILE        *ficlisup=NULL;
    FILE        *fnam=NULL;
    
    l2=num_seq;
    c2=pcoa->colm;
    l1=pcoa->rows;
    c1=pcoa->colm;

   if( (fnam=open_file("",option, "r",FALSE))==NULL) 
       my_exit(6,"sup row corresp");
      
                
    taballoc(&tabsup, l2, c2);
    lecmat(tabsup, nfictasup);

    taballoc(&compos, c1, pcoa->axis);
    lecmat(compos, "cbfcco");
    vecalloc(&moy, c1);
    vecalloc(&var, c1);

    vecalloc(&vp, pcoa->axis);
    lecvalpro(vp, nficvp);

    vecalloc(&poico, c1);
    lecvec(poico, "cbfcpc");
    
    for (j=1;j<=pcoa->axis;j++) {
        vp[j] = sqrt((double)vp[j]);
        a1 = vp[j];
        for (i=1;i<=c1;i++) {
            compos[i][j] = compos[i][j] / a1;
        }
    }
    for (i=1;i<=c1;i++) {
        a1 = poico[i];
        for (j=1;j<=pcoa->axis;j++) {
            compos[i][j] = compos[i][j] * a1;
        }
    }
     
 /* Transform genes with the initial factor                               */

    for (i=1;i<=l2;i++) {
            a1 = 0.0;
            for (j=1;j<=c1;j++) {
                a1 = a1 + tabsup[i][j];
            }
            if (a1 != 0.) {
                for (j=1;j<=c1;j++) {
                    a2 = tabsup[i][j] / a1;
                    if (poico[j]!=0) {tabsup[i][j] = a2 / poico[j];}
                }
            }
        }
    
  /* Position the suppli. genes on the original factors                    */
  
  if( (ficlisup = open_file("",nficlisup,"a",FALSE))==NULL ) 
         my_exit(1,"nficlisup");

   fprintf(summary,"\n\nThe position of each additional gene by axis " 
                "(see also %s )\n",option);

    fprintf(summary, "Additional genes added after COA: \n");
    fprintf(summary, "Number of genes: %d, columns: %d\n\n", l1, c1);

      
    for (i=1;i<=l2;i++) { 
      fgets(pm->junk,BUFSIZ,fnam);
      pm->junk[35]='\0';
      for ( j=35 ; j>=0; j--) 
          if ( isspace((int)pm->junk[j]) ) pm->junk[j]='\0';
       fprintf(ficlisup, "%s%c",pm->junk,pm->seperator);
       fprintf(summary , "%s%c",pm->junk,pm->seperator);
    
        for (k=1;k<=pcoa->axis;k++) {
            a1 = 0.;
            for (j=1;j<=c1;j++) {
                a1 = a1 + tabsup[i][j] * compos[j][k];
            }
            fprintf(ficlisup,"%f%c",(float)a1,pm->seperator);
            fprintf(summary ,"%10.5f%c",(float)a1,pm->seperator);
        }
    fprintf(ficlisup,"\n");
    fprintf(summary ,"\n");
    }
    fclose(ficlisup);

    freetab (tabsup);
    freetab (compos);
    freevec (vp);
    freevec(poico);
    freevec(moy);
    freevec(var);
    fileclose(&fnam);
}

/**************  lecvalpro         ***************************************/
/* Read a vector from a binary formatted file                            */
/*************************************************************************/
void lecvalpro (double *v1, char *nfic)
{
    float   v2;
    int     i, c1;
    FILE    *fic=NULL;
    
    if ( (fic=open_file("",nfic,"rb",FALSE))==NULL) my_exit(6,"lecvalpro");

    c1 = (int) v1[0];
    for (i=1;i<=c1;i++) {
        if ( fread((char *)&v2, 4, 1, fic) != 1) {
        fprintf(pm->my_err,"(lecvalpro)");    
        my_exit(5,"lecvalpro");
        }
        v1[i] = v2;
        if ( fread((char *)&v2, 4, 1, fic) != 1)  {
        fprintf(pm->my_err,"(lecvalpro)");
        my_exit(5,"lecvalpro2");
        }
    }   
    fileclose(&fic);
}

